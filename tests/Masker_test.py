##########################################################################
"""unit testing module for the Tree.py class."""

import cgat.Masker as Masker
import unittest


class SegCheck(unittest.TestCase):

    mMasker = Masker.MaskerSeg()

    def testEmpty(self):
        """test empty input."""
        self.assertEqual(self.mMasker(""), "")

    def testProtein(self):
        """test protein input."""
        self.assertEqual(self.mMasker(
            "ACDEFGHIKLWWWWWWWWWWWWWWwwwwwwwwwwwacdefghikl"),
            "ACDEFGHIKLXXXXXXXXXXXXXXxxxxxxxxxxxacdefghikl")

    def testCoding(self):
        """test coding sequence input."""
        self.assertEqual(self.mMasker(
            "GCCTGCGACGAGTTCGGCCACATCAAGCT"
            "GTGGTGGTGGTGGTGGTGGTGGTGGTGGT"
            "GGTGGTGGTGGTGGTGGTGGTGGTGGTGG"
            "tggtggtggtggtggtgggcctgcgacga"
            "gttcggccacatcaagctg"),
            "GCCTGCGACGAGTTCGGCCACATCAAGCT"
            "GNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "nnnnnnnnnnnnnnnnnngcctgcgacga"
            "gttcggccacatcaagctg")


class DustMaskerCheck(unittest.TestCase):
    mMasker = Masker.MaskerDustMasker()

if __name__ == "__main__":
    unittest.main()
