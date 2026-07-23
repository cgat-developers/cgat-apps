##########################################################################
"""unit testing module for the Tree.py class."""

import shutil

import cgat.Masker as Masker
import unittest

HAS_SEGMASKER = shutil.which("segmasker") is not None


class SegCheck(unittest.TestCase):

    mMasker = Masker.MaskerSeg()

    def testEmpty(self):
        """test empty input."""
        self.assertEqual(self.mMasker(""), "")

    @unittest.skipUnless(HAS_SEGMASKER, "requires segmasker")
    def testProtein(self):
        """test protein input."""
        self.assertEqual(self.mMasker(
            "ACDEFGHIKLWWWWWWWWWWWWWWwwwwwwwwwwwacdefghikl"),
            "ACDEFGHIKLXXXXXXXXXXXXXXxxxxxxxxxxxacdefghikl")

    @unittest.skipUnless(HAS_SEGMASKER, "requires segmasker")
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
