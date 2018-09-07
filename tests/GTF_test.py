import unittest
import os
import cgatcore.iotools as IOTools
import cgat.GTF as GTF


class TestIteration(unittest.TestCase):

    filename = os.path.join(os.path.dirname(__file__), "data", "hg19.small.gtf.gz")

    def test_number_of_intervals_is_correct(self):

        with IOTools.open_file(self.filename) as inf:
            records = list(GTF.iterator(inf))

        self.assertEqual(len(records),
                         100)


if __name__ == "__main__":
    unittest.main()
