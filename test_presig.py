import unittest
from presig.presig import *

class PSTests(unittest.TestCase):


# classify_SBS_feature(ref_allele, alt_allele, ref_context_fiveprime, ref_context_threeprime):
    def test_sbs_t1(self):
        sbs = classify_SBS_feature("A", "C", "AAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAA")
        self.assertEqual(sbs, "A[T>C]A")
    def test_sbs_t2(self):
        sbs = classify_SBS_feature("C", "G", "AAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAA")
        self.assertEqual(sbs, "A[C>G]A")
    def test_sbs_t3(self):
        sbs = classify_SBS_feature("C", "G", "T", "C")
        self.assertEqual(sbs, "T[C>G]C")

#   @unittest.skip("not finished yet")

    def test_rpt_t1(self):
        #C -       CAGACTCTTGACCTTAGGCAGTTTG       CCACCTCAGCCTCCTAAAGTGCTGG
        rpt, rpt_len = detect_repeat("C", "-", 1, 0, "CAGACTCTTGACCTTAGGCAGTTTG", "CCACCTCAGCCTCCTAAAGTGCTGG")
        self.assertEqual(rpt_len, 2)
        self.assertTrue(rpt)
    def test_rpt_t2(self):
        #C -       CAGACTCTTGACCTTAGGCAGTTTG       CCACCTCAGCCTCCTAAAGTGCTGG
        rpt, rpt_len = detect_repeat("GC", "G-", 2, 1, "CAGACTCTTGACCTTAGGCAGTTTG", "CCACCTCAGCCTCCTAAAGTGCTGG")
        self.assertEqual(rpt_len, 2)
        self.assertTrue(rpt)

    def test_mh_t1(self):
        #AAAT    -       GGCTCAAGGACTGGGTGAGCCACTC       AAAGACCAAGAGTTTTTCTAGGCCC
        mh, mhlen = detect_microhomology("AAAT", "-",4,1, "GGCTCAAGGACTGGGTGAGCCACTC", "AAAGACCAAGAGTTTTTCTAGGCCC")
        self.assertEqual(mhlen, 3)
        self.assertTrue(mh)
        mh, mhlen = detect_microhomology("AAAT", "-",4,1, "GGCTCAAGGACTGGGTGAGCCCAAT", "ACCGACCAAGAGTTTTTCTAGGCCC")
        rpt, rpt_len = detect_repeat("AAAT", "-",4,1, "GGCTCAAGGACTGGGTGAGCCCAAT", "ACCGACCAAGAGTTTTTCTAGGCCC")
        self.assertTrue(mh, "5-prime MH is properly detected")
        self.assertEqual(mhlen, 3, "5-prime MH is the correct length")
        self.assertFalse(rpt)
        self.assertEqual(rpt_len, 0)

    def test_mh_2(self):
        rpt, rpt_len = detect_repeat("T", "-", 1, 0, "AGCCACCATGCCCGGCATATTTTTG", "TTCAATAAAATACAGCCAAGGTTTC")
        self.assertTrue(rpt)
        self.assertEqual(rpt_len, 2)

        rpt, rpt_len = detect_repeat("T", "-", 1, 0, "AGCCACCATGCCCGGCATATTTTTT", "TTCAATAAAATACAGCCAAGGTTTC")
        self.assertTrue(rpt)
        self.assertEqual(rpt_len, 5, "5-prime single base repeat motifs are properly detected")
    
    def test_ins_1(self):
        #"1:Ins:T:5       111283679       111283680       INS     -       T       CTCACTCAAAGCTTCTTGGGTTTGT       TTTTTCCCAGGCCATGGTCATTCCT"
        rpt, rpt_len = detect_repeat("-", "T", 0, 1, "CTCACTCAAAGCTTCTTGGGTTTGT", "TTTTTCCCAGGCCATGGTCATTCCT")
        self.assertTrue(rpt)
        self.assertEqual(rpt_len, 5)

    def test_ins_2(self):
        rpt, rpt_len = detect_repeat("-", "TT", 0, 2, "CTCACTCAAAGCTTCTTGGGTTTGT", "TTTTTCCCAGGCCATGGTCATTCCT")
        self.assertTrue(rpt)
        self.assertEqual(rpt_len, 2)

    def test_ins_3(self):
        rpt, rpt_len = detect_repeat("-", "TC", 0, 2, "CTCACTCAAAGCTTCTTGGGTTTGT", "TCTCTCCCAGGCCATGGTCATTCCT")
        self.assertTrue(rpt)
        self.assertEqual(rpt_len, 3, "Motifs of length two are detected as repeats.")

    def test_ins_4(self):
        rpt, rpt_len = detect_repeat("-", "TCA", 0, 2, "CTCACTCAAAGCTTCTTGGGTTTCA", "TCATCACCAGGCCATGGTCATTCCT")
        self.assertTrue(rpt)
        self.assertEqual(rpt_len, 3, "5-prime motif sequences are properly detected for insertions.")
    
    def test_ins_mh_1(self):
        mh, mh_len = detect_microhomology("-", "TCGA", 0, 4, "CTCACTCAAAGCTTCTTGGGTTCCGA", "TCATCACCAGGCCATGGTCATTCCT")
        self.assertTrue(mh)
        self.assertEqual(mh_len, 3, "Microhomology is properly detected at insertions.")
    
    def test_ins_mh_2(self):
        mh, mh_len = detect_microhomology("A-", "ATGTTGTTT", 1, 9, "AGTTAATCCAGATAGAAGCACATTT", "GTTTGTTTTGTTTGTTTGTTTGTTT")
        self.assertTrue(mh)
        self.assertEqual(mh_len,3)

    def test_calculate_end_1(self):
        end = calculate_end_position(192092, "A", "T")
        self.assertEqual(end, 192092)
    
    def test_calculate_end_2(self):
        end = calculate_end_position(10000, "A", "AAA")
        self.assertEqual(end, 10000)
    
    def test_calculate_end_3(self):
        end = calculate_end_position(100, "AC", "GG")
        self.assertEqual(end, 101)

    def test_calculate_end_4(self):
        end = calculate_end_position(100, "AA", "A")
        self.assertEqual(end, 101)
        
if __name__ == '__main__':
    from pycotap import TAPTestRunner

    suite = unittest.TestLoader().loadTestsFromTestCase(PSTests)
    TAPTestRunner().run(suite)
