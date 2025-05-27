from . import lib
from .lib import (
    __version__,
    Gene,
    TRNAGene,
    TMRNAGene,
    RNAFinder,
    ARAGORN_VERSION,
    TRANSLATION_TABLES,
)

__doc__ = lib.__doc__
__all__ = [
    "Gene",
    "TRNAGene",
    "TMRNAGene",
    "RNAFinder",
    "ARAGORN_VERSION",
    "TRANSLATION_TABLES",
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPL-3.0-or-later"

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pyaragorn.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )
