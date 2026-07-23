import unittest
import io
import os
import cgatcore.iotools as iotools
import cgat.GTF as GTF
import cgat.Histogram as Histogram


class TestIteration(unittest.TestCase):

    filename = os.path.join(os.path.dirname(__file__), "data", "hg19.small.gtf.gz")

    def test_number_of_intervals_is_correct(self):

        with iotools.open_file(self.filename) as inf:
            records = list(GTF.iterator(inf))

        self.assertEqual(len(records),
                         100)

    def test_iterator_overlapping_genes_does_not_crash(self):
        """Regression test for tuple append bug in iterator_overlapping_genes."""
        gtf_data = io.StringIO(
            "chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id \"A\"; transcript_id \"A1\";\n"
            "chr1\tHAVANA\texon\t150\t250\t.\t+\t.\tgene_id \"B\"; transcript_id \"B1\";\n"
        )
        records = list(GTF.iterator(gtf_data))
        groups = list(GTF.iterator_overlapping_genes(iter(records)))
        self.assertIsInstance(groups, list)


class TestHistogram(unittest.TestCase):

    def test_histogram_mode_sorts_by_count(self):
        """Regression test for Py2 cmp/sort usage in histogram()."""
        result = Histogram.histogram([1, 1, 1, 2, 2, 3], mode=True)
        self.assertEqual(result, [(1, 3), (2, 2), (3, 1)])


if __name__ == "__main__":
    unittest.main()
