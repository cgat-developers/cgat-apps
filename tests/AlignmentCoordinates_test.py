"""Unit tests for cgat.AlignmentCoordinates."""

import unittest

from cgat.AlignmentCoordinates import (
    Alignment,
    addDiagonal2Alignment,
    formatEmissions,
    getAlignmentIdentity,
    getAlignmentOverlap,
    getAlignmentShortestDistance,
    makeAlignmentBlocks,
    makeAlignmentVector,
)


class AlignmentBuildCheck(unittest.TestCase):

    def test_add_diagonal(self):
        aln = Alignment()
        aln.addDiagonal(10, 15, 5)
        self.assertEqual(aln.getNumAligned(), 5)
        self.assertEqual(aln.mapRowToCol(10), 15)
        self.assertEqual(aln.mapRowToCol(14), 19)
        self.assertEqual(aln.mapRowToCol(15), -1)
        self.assertEqual(aln.getRowFrom(), 10)
        self.assertEqual(aln.getRowTo(), 15)  # exclusive end
        self.assertEqual(aln.getColFrom(), 15)
        self.assertEqual(aln.getColTo(), 20)  # exclusive end

    def test_add_blocks(self):
        aln = Alignment()
        aln.addBlocks([0, 20], [100, 130], [10, 5])
        self.assertEqual(aln.getNumAligned(), 15)
        self.assertEqual(aln.mapRowToCol(0), 100)
        self.assertEqual(aln.mapRowToCol(9), 109)
        self.assertEqual(aln.mapRowToCol(20), 130)
        self.assertEqual(aln.mapRowToCol(24), 134)
        self.assertEqual(aln.mapRowToCol(10), -1)

    def test_add_pair_explicit(self):
        aln = makeAlignmentVector()
        aln.addPairExplicit(1, 2, 0.0)
        aln.addPair(3, 4)
        self.assertEqual(aln.mapRowToCol(1), 2)
        self.assertEqual(aln.mapRowToCol(3), 4)

    def test_add_diagonal2alignment_helper(self):
        aln = makeAlignmentBlocks()
        addDiagonal2Alignment(aln, 0, 3, 10)
        self.assertEqual(list(aln.getBlocks()), [(0, 10, 3)])


class AlignmentAmbiguityCheck(unittest.TestCase):

    def test_conflicting_cols_raise(self):
        aln = Alignment()
        aln.addDiagonal(0, 5, 0)
        aln.addDiagonal(2, 4, 10)
        with self.assertRaises(RuntimeError):
            aln.mapRowToCol(0)

    def test_idempotent_add_ok(self):
        aln = Alignment()
        aln.addDiagonal(0, 5, 0)
        aln.addDiagonal(0, 5, 0)
        self.assertEqual(aln.mapRowToCol(2), 2)
        self.assertEqual(aln.getNumAligned(), 5)


class AlignmentBlocksCheck(unittest.TestCase):

    def test_get_blocks_merges_diagonal(self):
        aln = Alignment()
        for r in range(0, 5):
            aln.addPair(r, r + 10)
        for r in range(10, 13):
            aln.addPair(r, r + 20)
        self.assertEqual(aln.getBlocks(), [(0, 10, 5), (10, 30, 3)])

    def test_copy_and_get_new(self):
        src = Alignment()
        src.addDiagonal(5, 8, 1)
        dst = src.getNew()
        self.assertEqual(dst.getNumAligned(), 0)
        dst.copy(src)
        self.assertEqual(dst.getNumAligned(), 3)
        self.assertEqual(dst.mapRowToCol(6), 7)


class AlignmentCompareCheck(unittest.TestCase):

    def test_identity_and_overlap(self):
        a = Alignment()
        b = Alignment()
        a.addDiagonal(0, 10, 0)
        b.addDiagonal(5, 15, 0)
        # rows 5..9 mapped in both to the same col
        self.assertEqual(getAlignmentOverlap(a, b), 5)
        self.assertEqual(getAlignmentIdentity(a, b), 5)

        c = Alignment()
        c.addDiagonal(5, 15, 1)
        self.assertEqual(getAlignmentOverlap(a, c), 5)
        self.assertEqual(getAlignmentIdentity(a, c), 0)

    def test_shortest_distance(self):
        a = Alignment()
        b = Alignment()
        a.addDiagonal(0, 10, 0)
        b.addDiagonal(15, 20, 0)
        self.assertEqual(getAlignmentShortestDistance(a, b), 5)
        self.assertEqual(getAlignmentShortestDistance(a, a), 0)

    def test_format_emissions(self):
        aln = Alignment()
        aln.addDiagonal(1, 6, 0)
        text = formatEmissions(aln)
        self.assertIn("+5+5", text)


if __name__ == "__main__":
    unittest.main()
