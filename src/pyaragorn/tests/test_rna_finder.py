import itertools
import re
import unittest

from .. import RNAFinder
from . import data


_TRNA_RX = re.compile(r"^(\d+)\s+tRNA-([A-Za-z]{3})\s+(c?)\[(\d+),(\d+)\]\s+([\d.]+)\s+(\d+)\s+\(([a-z]{2,4})\)")


class TestMeta(unittest.TestCase):
    
    def test_default(self):
        record = data.load_record("CP001621.fna.gz")
        lines = data.load_text("CP001621.default.txt").splitlines()
        
        finder = RNAFinder(translation_table=11)
        genes = finder.find_rna(str(record.seq))

        for gene, expected in itertools.zip_longest(genes, itertools.batched(lines[2:], 3)):
            self.assertIsNotNone(gene)
            self.assertIsNotNone(expected)
            result, seq, ss = expected
            if gene.type == "tRNA":
                matched = _TRNA_RX.match(result)
                _, aa, complement, begin, end, energy, offset, anticodon = matched.groups()
                self.assertEqual(gene.amino_acid, aa)
                self.assertEqual(gene.begin, int(begin))
                self.assertEqual(gene.end, int(end))
                self.assertEqual(gene.anticodon_offset, int(offset))
                self.assertEqual(gene.anticodon, anticodon)
                self.assertEqual(gene.strand, -1 if complement == "c" else +1)
                self.assertAlmostEqual(gene.energy, float(energy), places=1)
                self.assertEqual(gene.sequence().lower(), seq)


