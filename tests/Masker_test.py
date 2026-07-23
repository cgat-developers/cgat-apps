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

    def test_get_command(self):
        """test argv list is built without shell interpolation."""
        cmd = self.mMasker.getCommand("/tmp/test.fa")
        self.assertEqual(cmd[0], "dustmasker")
        self.assertIn("/tmp/test.fa", cmd)

if __name__ == "__main__":
    unittest.main()
