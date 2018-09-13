"""unit testing module for the Tree.py class."""

import cgat.Intervals as Intervals
import unittest


class TruncateCheck(unittest.TestCase):

    def testEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.truncate([], []), [])

    def testHalfEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.truncate([], [(0, 5)]), [])
        self.assertEqual(Intervals.truncate([(0, 5)], []), [(0, 5)])

    def testSingle(self):
        """test empty input."""
        self.assertEqual(Intervals.truncate([(0, 5)], [(0, 5)]), [])
        self.assertEqual(Intervals.truncate([(0, 5)], [(0, 3)]), [(3, 5)])
        self.assertEqual(Intervals.truncate([(0, 3)], [(0, 5)]), [])
        self.assertEqual(Intervals.truncate([(0, 5)], [(3, 5)]), [(0, 3)])
        self.assertEqual(Intervals.truncate([(3, 5)], [(0, 5)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(5, 10)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(5, 20)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(0, 10)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(0, 10)]), [])
        self.assertEqual(Intervals.truncate([(5, 10)], [(0, 20)]), [])

    def testMultiple(self):
        """test empty input."""
        self.assertEqual(
            Intervals.truncate([(0, 5), (10, 15)], [(0, 5)]), [(10, 15)])
        self.assertEqual(
            Intervals.truncate([(0, 5), (10, 15)], [(0, 10)]), [(10, 15)])
        self.assertEqual(Intervals.truncate([(0, 5), (10, 15)], [(0, 15)]), [])
        self.assertEqual(Intervals.truncate([(0, 5), (5, 10)], [(0, 10)]), [])
        self.assertEqual(
            Intervals.truncate([(0, 5), (5, 10)], []), [(0, 5), (5, 10)])

    def testNoOverlap(self):
        """test empty input."""
        self.assertEqual(
            Intervals.truncate([(0, 5), (10, 15)], [(5, 10)]), [(0, 5), (10, 15)])
        self.assertEqual(
            Intervals.truncate([(5, 10)], [(0, 5), (10, 15)]), [(5, 10)])
        self.assertEqual(
            Intervals.truncate([(0, 5), (5, 10)], [(10, 15)]), [(0, 5), (5, 10)])


class IntersectCheck(unittest.TestCase):

    def testEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.intersect([], []), [])

    def testHalfEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.intersect([(0, 5)], []), [])
        self.assertEqual(Intervals.intersect([], [(0, 5)]), [])

    def testSingle(self):
        """test empty input."""
        self.assertEqual(Intervals.intersect([(0, 5)], [(0, 5)]), [(0, 5)])
        self.assertEqual(Intervals.intersect([(0, 5)], [(0, 3)]), [(0, 3)])
        self.assertEqual(Intervals.intersect([(0, 3)], [(0, 5)]), [(0, 3)])
        self.assertEqual(Intervals.intersect([(0, 5)], [(3, 5)]), [(3, 5)])
        self.assertEqual(Intervals.intersect([(3, 5)], [(0, 5)]), [(3, 5)])
        self.assertEqual(Intervals.intersect([(5, 10)], [(5, 20)]), [(5, 10)])
        self.assertEqual(Intervals.intersect([(5, 10)], [(0, 20)]), [(5, 10)])

    def testMultiple(self):
        """test empty input."""
        self.assertEqual(
            Intervals.intersect([(0, 5), (10, 15)], [(0, 5)]), [(0, 5)])
        self.assertEqual(
            Intervals.intersect([(0, 5), (10, 15)], [(0, 10)]), [(0, 5)])
        self.assertEqual(
            Intervals.intersect([(0, 5), (10, 15)], [(0, 15)]), [(0, 5), (10, 15)])
        self.assertEqual(
            Intervals.intersect([(0, 5), (5, 10)], [(0, 10)]), [(0, 5), (5, 10)])

    def testNoOverlap(self):
        """test empty input."""
        self.assertEqual(
            Intervals.intersect([(0, 5), (10, 15)], [(5, 10)]), [])
        self.assertEqual(
            Intervals.intersect([(5, 10)], [(0, 5), (10, 15)]), [])


class FromArrayCheck(unittest.TestCase):

    def testEmpty(self):
        """test empty input."""
        self.assertEqual(Intervals.fromArray([]), [])

    def testArray1(self):
        """test simple array."""
        a = [1, 1, 1, 0, 0, 0, 1, 1, 1]
        self.assertEqual(Intervals.fromArray(a), [(0, 3), (6, 9)])
        self.assertEqual(Intervals.fromArray([not x for x in a]), [(3, 6)])

    def testArray2(self):
        """test longer array."""
        a = [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]
        self.assertEqual(Intervals.fromArray(a), [(0, 3), (6, 9), (12, 15)])
        self.assertEqual(
            Intervals.fromArray([not x for x in a]), [(3, 6), (9, 12)])


if __name__ == "__main__":
    unittest.main()
