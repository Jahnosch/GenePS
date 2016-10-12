#!/usr/bin/env python3
import unittest
import os
from run_command import tempdir
from run_GenePS import ResultsObject, ExonerateError, ScoreError, judge_score, coverage_filter, score_prediction, \
    make_prediction
from find_regions import run_tblastn, parse_blastdb, make_blast_db, BlastObject, HspListObject
from exonerate_parser import run_exonerate, remove_non_letter_signs, grap_values, aacode_3to1, ExonerateObject


script_path = os.path.dirname(os.path.realpath(__file__))
folder_path = os.path.join(script_path, "test_data")
test_data = os.path.join(folder_path, "group1.makeGenePS")
test_db = os.path.join(folder_path, "c_elegans.PRJNA13758.WS254.genomic.fa")


class TestResultsObject(unittest.TestCase):

    group_result = ResultsObject(test_data)
    group_result.read_gene_ps_consensus_file()

    def test_path(self):
        self.assertTrue(self.group_result.path == test_data)

    def test_group(self):
        self.assertTrue(self.group_result.group_name == "group1")

    def test_group_size(self):
        self.assertEqual(self.group_result.group_size, 2)

    def test_seq_length(self):
        self.assertEqual(self.group_result.seq_length["20_seqs_eef_test"], 850)

    def test_consensus(self):
        self.assertEqual(len(self.group_result.consensus["20_seqs_eef_test"]), 850)

    def test_phmm(self):
        self.assertIn(self.group_result.phmm["20_seqs_eef_test"],
                      "/home/jgravemeyer/Dropbox/MSc_project/data/test_out/20_seqs_eef_test.hmmGenePS")

    def test_score_list(self):
        self.assertListEqual(self.group_result.score_list["20_seqs_eef_test"], [1981, 1808, 1893, 1981, 1981, 1889])

    def test_consensus_to_fa_correct_path(self):
        group_cons = self.group_result.consensus_to_fa("20_seqs_eef_test", folder_path + "/test")
        self.assertTrue(os.path.exists(group_cons))
        os.remove(group_cons)

    def test_consensus_to_fa_correct_header_content(self):
        header = self.group_result.consensus.keys()
        group_cons = self.group_result.consensus_to_fa(header, folder_path + "/test")
        file_list = [line.strip() for line in open(group_cons)]
        self.assertIn(">20_seqs_eef_test", file_list)
        os.remove(group_cons)

    def test_consensus_to_fa_correct_size_content(self):
        header = self.group_result.consensus.keys()
        group_cons = self.group_result.consensus_to_fa(header, folder_path + "/test")
        file_list = [line.strip() for line in open(group_cons)]
        self.assertEqual(len(file_list), 4)
        os.remove(group_cons)

    def test_consensus_to_fa_correct_seq_content(self):
        header = self.group_result.consensus.keys()
        group_cons = self.group_result.consensus_to_fa(header, folder_path + "/test")
        file_list = [line.strip() for line in open(group_cons)]
        self.assertEqual(len(file_list[3]) + len(file_list[1]), 1012+850)
        os.remove(group_cons)


class TestBlast(unittest.TestCase):

    test_blast = os.path.join(folder_path, "test_tblastn.out")
    single_hit_blast = os.path.join(folder_path, "single_hit_test_blast.out")

    def test_make_blast_db(self):
        db_dir = make_blast_db(test_db, folder_path)
        self.assertIn(db_dir, test_db)

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_make_blast_db_command(self):
        with tempdir() as tmp:
            db_dir = make_blast_db(test_db, tmp)
            self.assertIn(db_dir, tmp + "/c_elegans.PRJNA13758.WS254.genomic.fa")

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_parse_blast_db_command(self):
        blast_parse_out = parse_blastdb(test_db, "I", 9158949, 9171695)
        self.assertIn(">I;9158949-9171695\nAGGAGCA", blast_parse_out)

    def test_parse_blast_db_without_command(self):
        test_file = os.path.join(folder_path, "blast_region.fasta")
        command = ["less", test_file]
        blast_parse_out = parse_blastdb(test_db, "I", 9158949, 9171695, command)
        self.assertIn(">I;9158949-9171695\nAGGAGCA", blast_parse_out)

    def test_run_tblastn_all_columns_of_row_present(self):
        test_cmd = ["less", self.test_blast]
        blast_obj = run_tblastn("no_db", "no_query", test_cmd)
        self.assertEqual(len(blast_obj.blast_out["20_seqs_eef_test"]["I"][0]), 8)

    def test_run_tblastn_all_rows_of_contig_present(self):
        test_cmd = ["less", self.test_blast]
        blast_obj = run_tblastn("no_db", "no_query", test_cmd)
        self.assertEqual(len(blast_obj.blast_out["20_seqs_eef_test"]["I"]), 6)

    def test_merging_dist_flanking_dist(self):
        test_cmd = ["less", self.test_blast]
        blast_obj = run_tblastn("no_db", "no_query", test_cmd)
        self.assertEqual(blast_obj.flanking_distance, 5000)
        self.assertEqual(blast_obj.merging_distance, 10000)

    def test_infer_regions_chunck_coverage(self):
        test_cmd = ["less", self.test_blast]
        blast_obj = run_tblastn("no_db", "no_query", test_cmd)
        blast_obj.infer_regions()
        self.assertEqual(blast_obj.inferred_regions["20_seqs_eef_test"]["I"][1].chunk_cov, 93)

    def test_infer_regions_query_coverage(self):
        test_cmd = ["less", self.test_blast]
        blast_obj = run_tblastn("no_db", "no_query", test_cmd)
        blast_obj.infer_regions()
        self.assertEqual(blast_obj.inferred_regions["20_seqs_eef_test"]["I"][1].query_cov, 63)

    def test_infer_regions_query_s_start(self):
        test_cmd = ["less", self.test_blast]
        blast_obj = run_tblastn("no_db", "no_query", test_cmd)
        blast_obj.infer_regions()
        self.assertEqual(blast_obj.inferred_regions["20_seqs_eef_test"]["I"][0].s_start, 9158949)

    def test_infer_regions_single_hit_blast(self):
        test_cmd = ["less", self.single_hit_blast]
        blast_obj = run_tblastn("no_db", "no_query", test_cmd)
        blast_obj.infer_regions()
        self.assertEqual(blast_obj.inferred_regions["20_seqs_eef_test"]["I"][0].s_start, 9164821 - 5000)


class TestHspListObject(unittest.TestCase):

    hsp_list = [{'s_end': 9166695, 's_start': 9164821,
                 'contig': 'I', 'strand': '+', 'evalue': 0.0, 'q_end': 850, 'q_start': 244, 'q_len': 850},
                {'s_end': 9164804, 's_start': 9164208,
                 'contig': 'I', 'strand': '+', 'evalue': 2.18e-99, 'q_end': 253, 'q_start': 73, 'q_len': 850},
                {'s_end': 9164164, 's_start': 9163949,
                 'contig': 'I', 'strand': '+', 'evalue': 2.18e-99, 'q_end': 73, 'q_start': 2, 'q_len': 850},
                {'s_end': 9865045, 's_start': 9864005,
                 'contig': 'I', 'strand': '+', 'evalue': 2.09e-37, 'q_end': 416, 'q_start': 55, 'q_len': 850},
                {'s_end': 9863915, 's_start': 9863820,
                 'contig': 'I', 'strand': '+', 'evalue': 2.09e-37, 'q_end': 41, 'q_start': 10, 'q_len': 850},
                {'s_end': 9866277, 's_start': 9865864,
                 'contig': 'I', 'strand': '+', 'evalue': 3.83e-09, 'q_end': 579, 'q_start': 439, 'q_len': 850}]

    def test_merge_regions(self):
        hits = HspListObject(self.hsp_list, 10000)
        hits.sort_hsp_list()
        all_merged_regions = hits.merge_to_region()
        self.assertListEqual(all_merged_regions, [[0, 1, 2], [3, 4, 5]])

    def test_merge_regions_2_single_regions(self):
        hsp_list_single = [{'s_end': 9166695, 's_start': 9164821,
                 'contig': 'I', 'strand': '+', 'evalue': 0.0, 'q_end': 850, 'q_start': 244, 'q_len': 850},
                           {'s_end': 9866277, 's_start': 9865864,
                 'contig': 'I', 'strand': '+', 'evalue': 3.83e-09, 'q_end': 579, 'q_start': 439, 'q_len': 850}]
        hits = HspListObject(hsp_list_single, 10000)
        hits.sort_hsp_list()
        all_merged_regions = hits.merge_to_region()
        self.assertListEqual(all_merged_regions, [[0], [1]])

    def test_merge_regions_wrong_strand(self):
        hsp_list_single = [{'s_end': 9166695, 's_start': 9164821,
                 'contig': 'I', 'strand': '+', 'evalue': 0.0, 'q_end': 850, 'q_start': 244, 'q_len': 850},
                {'s_end': 9164804, 's_start': 9164208,
                 'contig': 'I', 'strand': '-', 'evalue': 2.18e-99, 'q_end': 253, 'q_start': 73, 'q_len': 850}]
        hits = HspListObject(hsp_list_single, 10000)
        hits.sort_hsp_list()
        all_merged_regions = hits.merge_to_region()
        self.assertNotEqual(all_merged_regions, [[0, 1]])

    def test_merge_regions_merging_distance(self):
        hsp_list_single = [{'s_end': 9166695, 's_start': 9164821,
                 'contig': 'I', 'strand': '+', 'evalue': 0.0, 'q_end': 850, 'q_start': 244, 'q_len': 850},
                {'s_end': 9164804, 's_start': 9164208,
                 'contig': 'I', 'strand': '-', 'evalue': 2.18e-99, 'q_end': 253, 'q_start': 73, 'q_len': 850}]
        hits = HspListObject(hsp_list_single, 1)
        hits.sort_hsp_list()
        all_merged_regions = hits.merge_to_region()
        self.assertNotEqual(all_merged_regions, [[0, 1]])

    def test_attributes(self):
        hits = HspListObject(self.hsp_list, 10000)
        hits.sort_hsp_list()
        hits.merge_to_region()
        self.assertEqual(hits.q_end[2], 850)
        self.assertEqual(hits.q_start[0], 2)
        self.assertEqual(hits.s_end[2], 9166695)
        self.assertEqual(hits.s_start[0], 9163949)
        self.assertTrue(hits.strand[0] == "+")
        self.assertEqual(hits.q_len[0], 850)

    def test_compute_cov(self):
        hits = HspListObject(self.hsp_list, 10000)
        hits.sort_hsp_list()
        hits.merge_to_region()
        q_start_pos = hits.q_start[0:3]
        q_end_pos = hits.q_end[0:3]
        chunck_cov, query_cov = hits.compute_coverage(q_start_pos, q_end_pos, hits.q_len[0])
        self.assertEqual(query_cov, 100)
        self.assertEqual(chunck_cov, 100)

class TestExonerateParser(unittest.TestCase):

    exonerate_file = folder_path + "/test_exonerate.out"
    exo_obj = ExonerateObject(exonerate_file)

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_run_exonerate(self):
        region_file = folder_path + "/blast_region.fasta"
        query_file = folder_path + "/test_consensus_eef.fa"
        with tempdir() as tmp:
            exo_ob = run_exonerate("test_exo", tmp, region_file, query_file)
            self.assertTrue(os.path.exists(os.path.join(tmp, "test_exo")))
            self.assertEqual(len(exo_ob.header), 1)

    def test_grap_values(self):
        self.assertTrue(type(grap_values(self.exo_obj.target_prot)) == list)

    def test_exonerate_processor_attributes(self):
        self.assertTrue(grap_values(self.exo_obj.header)[0]["query"] == "20_seqs_eef_test")
        self.assertTrue(self.exo_obj.path == self.exonerate_file)
        self.assertIn("MVNFTVDEIRA", grap_values(self.exo_obj.target_prot)[0])
        self.assertIn("ATGgtagGTCAACTT", grap_values(self.exo_obj.target_dna)[0])
        self.assertIn("MetValAsnPheThrValAspGluIle", grap_values(self.exo_obj.query_prot)[0])
        self.assertIn(['gene', '4706', '7747', '+'], grap_values(self.exo_obj.gff)[0])

    def test_aacode_3to1(self):
        seq = aacode_3to1("MetValAsnPheThrValAspGluIle")
        self.assertIn("MVNFTVDEI", seq)

    def test_remove_non_letter_signs_string_with_signs(self):
        seq = "MV*NFTVDEI--"
        new_seq = remove_non_letter_signs(seq)
        self.assertFalse(seq == new_seq)

    def test_remove_non_letter_signs_normal_seq(self):
        seq = "MVNFTVDEI"
        new_seq = remove_non_letter_signs(seq)
        self.assertTrue(seq == new_seq)


class TestGlobalGenePS(unittest.TestCase):

    exonerate_file = folder_path + "/test_exonerate.out"
    exo_obj = ExonerateObject(exonerate_file)
    hmm = folder_path + "/test_hmm_eef.hmmGenePS"
    test_blast = os.path.join(folder_path, "test_tblastn.out")
    single_hit_blast = os.path.join(folder_path, "single_hit_test_blast.out")
    test_consensus = folder_path + "/test_consensus_eef.fa"

    def test_judge_score_in_range(self):
        score_list = [1, 2, 3, 4, 5]
        test_score = 0
        result_score = judge_score(score_list, test_score)
        self.assertTrue(result_score["adj_aver"] > result_score["conf_int"][0])

    def test_judge_score_out_range(self):
        score_list = [1, 2, 3, 4, 5]
        test_score = 10
        self.assertRaises(ScoreError, lambda: judge_score(score_list, test_score))

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_score_prediction(self):
        score = score_prediction(self.exo_obj, self.hmm)
        self.assertEqual(score, 1899.9)

    @unittest.skipIf("travis" in os.getcwd(), "does not run on travis ci")
    def test_make_prediction_corret_output(self):
        test_cmd = ["less", self.test_blast]
        blast_obj = run_tblastn("no_db", "no_query", test_cmd)
        blast_obj.infer_regions()
        region = blast_obj.inferred_regions["20_seqs_eef_test"]["I"][0]
        with tempdir() as tmp:
            exo_obj = make_prediction("test_precition", self.test_consensus, tmp, region, test_db)
            self.assertTrue(grap_values(exo_obj.header)[0]["query"] == "20_seqs_eef_test")

if __name__ == '__main__':
    unittest.main()
