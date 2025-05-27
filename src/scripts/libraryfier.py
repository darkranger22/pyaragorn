import argparse
import copy
import os
import sys
import itertools
import pathlib
import tempfile
import subprocess

import pycparser
import pycparser.c_ast
from pycparser.c_ast import NodeVisitor
from pycparser.c_generator import CGenerator

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=pathlib.Path, required=True, help="path to the source file")
parser.add_argument("-o", "--output", type=pathlib.Path, required=True, help="path to the output folder")
parser.add_argument("--cpp", default="cpp", help="name of the C preprocessor executable")
args = parser.parse_args()

# record contents, includes & defines
contents = ["typedef void FILE;\n"]  # needed for pycparser to be happy about 
includes = ["#include <stdint.h>\n"]
defines = []

# load file, extract defines and includes
with args.input.open() as f:
    for line in map(str.lstrip, f):
        if line.startswith(("#define", "#ifndef", "#endif")):
            defines.append(line)
        elif line.startswith("#include"):
            includes.append(line)
        else:
            contents.append(line)

# preprocess contents and load AST
with tempfile.TemporaryDirectory() as folder:
    source = pathlib.Path(folder).joinpath("aragorn.c")
    with source.open("w") as tmp:
        tmp.writelines(contents)
    flag = "/EP" if os.name == "nt" else "-E"
    proc = subprocess.run([args.cpp, flag, os.fspath(source)], stdout=subprocess.PIPE)
    proc.check_returncode()
    parser = pycparser.CParser()
    ast = parser.parse("\n".join(proc.stdout.decode().splitlines()))

# record typedefs, constant arrays, functions
typedefs = []
functions = []
arrays = []
data = []

# process AST, extract functions, constant arrays, typedefs
for i, (_, node) in enumerate(itertools.islice(ast.children(), 1, None)):
    # extract typedef
    if isinstance(node, pycparser.c_ast.Typedef):
        typedefs.append(node)
    # extract constant array
    elif isinstance(node, pycparser.c_ast.Decl) and isinstance(node.type, pycparser.c_ast.ArrayDecl):
        arrays.append(node)
    # extract function definition
    elif isinstance(node, pycparser.c_ast.FuncDef):
        functions.append(node)

# --- Transform AST ------------------------------------------------------------

# add "genes" field to the `csw` struct
node = next(ty for ty in typedefs if ty.name == "csw")
node.type.type.decls.insert(1, 
    pycparser.c_ast.Decl(
        name='genes', 
        quals=[], 
        align=[], 
        storage=[], 
        funcspec=[], 
        init=None,
        bitsize=None, 
        type=pycparser.c_ast.TypeDecl(
            declname='genes',
            quals=[],
            align=None,
            type=pycparser.c_ast.IdentifierType(names=['gene*'])
        )
    )
)

# replace all uses of ts
class FunctionPatcher(NodeVisitor):

    @property
    def _replacement(self):
        return pycparser.c_ast.StructRef(
            name=pycparser.c_ast.ID(name='sw'),
            type='->',
            field=pycparser.c_ast.ID(name='genes')
        )

    def visit_Assignment(self, node):
        if isinstance(node.lvalue, pycparser.c_ast.ID) and node.lvalue.name == "ts":
            node.lvalue = self._replacement
        super().generic_visit(node)

    def visit_ArrayRef(self, node):
        if isinstance(node.name, pycparser.c_ast.ID) and node.name.name == "ts":
            node.name = self._replacement
        super().generic_visit(node)

    def visit_FuncCall(self, node):
        if node.name.name == "init_gene":
            node.args.exprs.insert(0, self._replacement)
        elif node.name.name == "nearest_tmrna_gene":
            node.args.exprs.append(pycparser.c_ast.ID(name="sw"))
        elif node.args is not None:
            for i, expr in enumerate(node.args.exprs):
                if isinstance(expr, pycparser.c_ast.ID) and expr.name == "ts":
                    node.args.exprs[i] = self._replacement

        super().generic_visit(node)

    def visit_BinaryOp(self, node):
        if isinstance(node.left, pycparser.c_ast.ID) and node.left.name == "ts":
            node.left = pycparser.c_ast.StructRef(
                name=pycparser.c_ast.ID(name='sw'),
                type='->',
                field=pycparser.c_ast.ID(name='genes')
            )
        super().generic_visit(node)

    def visit_Cast(self, node):
        if isinstance(node.expr, pycparser.c_ast.ID) and node.expr.name == "ts":
            node.expr = self._replacement
        super().generic_visit(node)

for node in functions:
    if node.decl.name != "init_gene":
        FunctionPatcher().visit(node)

# add sw argument to `nearest_tmrna_gene`
node = next(f for f in functions if f.decl.name == "nearest_tmrna_gene")
node.decl.type.args.params.append(
    pycparser.c_ast.Decl(
        name='sw', 
        quals=[], 
        align=[], 
        storage=[], 
        funcspec=[], 
        init=None,
        bitsize=None, 
        type=pycparser.c_ast.TypeDecl(
            declname='sw',
            quals=[],
            align=None,
            type=pycparser.c_ast.IdentifierType(names=['csw*'])
        )
    )
)

# add "ts" argument to `init_gene`
node = next(f for f in functions if f.decl.name == "init_gene")
node.decl.type.args.params.insert(0,
    pycparser.c_ast.Decl(
        name='ts', 
        quals=[], 
        align=[], 
        storage=[], 
        funcspec=[], 
        init=None,
        bitsize=None, 
        type=pycparser.c_ast.TypeDecl(
            declname='ts',
            quals=[],
            align=None,
            type=pycparser.c_ast.IdentifierType(names=['gene*'])
        )
    )
)

# reduce size of aamap array (int to int8_t)
aamap = next(a for a in arrays if a.name == "aamap")
aamap.type.type.type.type = pycparser.c_ast.ID(name="int8_t")

# extract the map array from move_forward
move_forward = next(f for f in functions if f.decl.name == "move_forward")
map_array = move_forward.body.block_items.pop(2)
map_array.storage = []
arrays.append(map_array)

# remove main function
main = next(f for f in functions if f.decl.name == "main")
functions.remove(main)

# remove helpmenu array and function
helpmenu = next(a for a in arrays if a.name == "helpmenu")
arrays.remove(helpmenu)
aragorn_help_menu = next(f for f in functions if f.decl.name == "aragorn_help_menu")
functions.remove(aragorn_help_menu)


# --- Emit code ----------------------------------------------------------------

args.output.mkdir(exist_ok=True)
gen = CGenerator()

# write common header
with args.output.joinpath("aragorn.h").open("w") as dst:
    dst.write("#ifndef ARAGORN_H\n")
    dst.write("#define ARAGORN_H\n")
    # write includes and defines from the original file
    dst.writelines(includes)
    dst.writelines(defines)
    # write typedefs
    for node in typedefs:
        dst.write(gen.visit(node) + ";\n")
    # write arrays
    for node in arrays:
        node = copy.deepcopy(node)
        node.init = None
        node.storage = ["extern"]
        dst.write(gen.visit(node) + ";\n")
    # write function signatures
    for node in functions:
        node = copy.deepcopy(node.decl)
        node.storage = ["extern"]
        dst.write(gen.visit(node) + ";\n")
    dst.write("#endif\n")


# write common source
with args.output.joinpath("aragorn.c").open("w") as dst:
    dst.write(f'#include "aragorn.h"\n')
    # write function definitions
    for node in functions:
        dst.write(gen.visit(node)+"\n")

# with args.output.joinpath("main.c").open("w") as dst:
#     dst.write('#include "aragorn.h"\n')
#     dst.write(gen.visit(main)+"\n")

# write data file
with args.output.joinpath("data.c").open("w") as dst:
    dst.write('#include "aragorn.h"\n')
    # write arrays
    for node in arrays:
        dst.write(gen.visit(node) + ";\n")
