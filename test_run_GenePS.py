#!/usr/bin/env python3
import unittest
import os
from unittest.mock import patch
from run_command import tempdir
import run_GenePS
from find_regions import read_blast_output, HspListObject
from exonerate_parser import remove_non_letter_signs, clear_hashed_bases, aacode_3to1, ExonerateObject

script_path = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(script_path, "test_data")


class TestBlastObject(unittest.TestCase):

    blast_file = os.path.join(test_data, "run_geneps/remanei_OG1685.blast")
    blast = read_blast_output(blast_file, os.path.join(test_data, "databases/c_brenneri.PRJNA20035.WS249.genomic.fa"))
    blast.infer_regions()

    def test_inferred_regions_single_hsp_in_hsp_list(self):
        single_hsp_region = self.blast.inferred_regions["Cbre_Contig419"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region.s_start, 60160)
        self.assertEqual(single_hsp_region.s_end, 70339)
        self.assertEqual(single_hsp_region.strand, "+")
        self.assertEqual(single_hsp_region.chunk_cov, 100)
        self.assertEqual(single_hsp_region.query_cov, 11)
        self.assertEqual(single_hsp_region.q_len, 466)

    def test_inferred_regions_merging_reverse_compliment(self):
        single_hsp_region = self.blast.inferred_regions["Cbre_Contig6"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region.s_start, 1479061)
        self.assertEqual(single_hsp_region.s_end, 1492025)
        self.assertEqual(single_hsp_region.strand, "-")
        self.assertEqual(single_hsp_region.chunk_cov, 82)
        self.assertEqual(single_hsp_region.query_cov, 82)
        self.assertEqual(single_hsp_region.q_len, 466)

    def test_inferred_regions_merging_plus_strand(self):
        single_hsp_region = self.blast.inferred_regions["Cbre_Contig231"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region.s_start, 10349)
        self.assertEqual(single_hsp_region.s_end, 24214)
        self.assertEqual(single_hsp_region.strand, "+")
        self.assertEqual(single_hsp_region.chunk_cov, 82)
        self.assertEqual(single_hsp_region.query_cov, 82)
        self.assertEqual(single_hsp_region.q_len, 466)

    def test_unmerged_blast_output(self):
        single_hsp_region = self.blast.blast_out["Cbre_Contig419"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertEqual(single_hsp_region['s_start'], 65160)
        self.assertEqual(single_hsp_region['s_end'], 65339)

    def test_unmerged_blast_reverse_strand(self):
        blast = read_blast_output(self.blast_file, "fake_db")
        single_hsp_region = blast.blast_out["Cbre_Contig6"]["OrthologousGroups_I3.5.OGoverlapp.txt"][0]
        self.assertLess(single_hsp_region['s_end'], single_hsp_region['s_start'])

    def test_compare_merged_unmerged_number_regions(self):
        blast = read_blast_output(self.blast_file, "fake_db")
        number_unmerged_regions = len(blast.blast_out["Cbre_Contig231"]["OrthologousGroups_I3.5.OGoverlapp.txt"])
        self.assertEqual(number_unmerged_regions, 4)
        number_merged_regions = len(self.blast.inferred_regions["Cbre_Contig231"]["OrthologousGroups_I3.5.OGoverlapp.txt"])
        self.assertEqual(number_merged_regions, 1)

    def test_compute_coverage(self):
        test_obj = HspListObject("fake", "fake")
        chunk, query = test_obj.compute_coverage([35, 40, 80], [50, 60, 100], 100)
        self.assertEqual(chunk, 69)
        self.assertEqual(query, 45)
        chunk_single, query_single = test_obj.compute_coverage([35], [50], 100)
        self.assertEqual(chunk_single, 100)
        self.assertEqual(query_single, 15)


class TestExonerateObject(unittest.TestCase):

    # with hashes in sequence
    exonerate_obj = ExonerateObject(os.path.join(test_data, "run_geneps/elegans_eef_true.exonerate"))
    target_hash = "I;9859005-9873867"
    target_range_hash = ('4791', '9857')
    target_prot_hash = exonerate_obj.target_prot[target_hash][target_range_hash]
    target_dna_hash = exonerate_obj.target_dna[target_hash][target_range_hash]

    # no hashes in sequence
    target2 = 'I;9159307-9171689'
    target_range2 = ("819", "10142")
    target_prot2 = exonerate_obj.target_prot[target2][target_range2]
    target_dna2 = exonerate_obj.target_dna[target2][target_range2]

    def test_target_target_range(self):
        self.assertIn(self.target_hash, self.exonerate_obj.query_prot)
        self.assertIn(self.target_range_hash, list(self.exonerate_obj.query_prot[self.target_hash].keys()))

    def test_sequences_dividability(self):
        self.assertEqual(len(self.target_dna_hash) / 3, len(self.target_prot_hash))
        self.assertEqual(len(self.target_dna2) / 3, len(self.target_prot2))

    def test_protein_length(self):
        self.assertEqual(827, len(self.target_prot_hash))
        self.assertEqual(961, len(self.target_prot2))

    def test_dna_length(self):
        self.assertEqual(2481, len(self.target_dna_hash))
        self.assertEqual(2883, len(self.target_dna2))

    def test_gff_before_rewriting(self):
        first_element = 'I;9859005-9873867\texonerate:protein2genome:bestfit\tgene\t4792\t9857\t485\t+\t.\tgene_id 0 ; sequence eef_3.5 ; gene_orientation +'
        self.assertEqual(self.exonerate_obj.gff[self.target_hash][self.target_range_hash][0], first_element)

    def test_gff_cds_phase(self):
        self.assertEqual(self.exonerate_obj.gff[self.target_hash][self.target_range_hash][11].split("\t")[-2], '1')

    def test_gff_mRNA(self):
        self.assertEqual(self.exonerate_obj.gff[self.target_hash][self.target_range_hash][1].split("\t")[2], 'mRNA')

    def test_gff_no_splice_no_similarity(self):
        self.assertNotIn("similarity", self.exonerate_obj.gff[self.target_hash][self.target_range_hash][-1].split("\t")[2])
        self.assertIn("exon", self.exonerate_obj.gff[self.target_hash][self.target_range_hash][-1].split("\t")[2])
        self.assertNotIn("splice", self.exonerate_obj.gff[self.target_hash][self.target_range_hash][3].split("\t")[2])

    def test_aa3_to_1_coding_dict(self):
        test_string = 'Val***XaaUnkSer'
        out_string = 'VXXXS'
        failed_string = "ValUnk*"
        self.assertEqual(out_string, aacode_3to1(test_string))
        self.assertRaises(Exception, aacode_3to1, failed_string)

    def test_clear_hashed_bases(self):
        testprot = "        sHisValGluLeuLeu##LysSer------------------SerGlyIleSerLeuLeuCy"
        testdna = " 5030 : ACACGTGGAATTACTATGAAATCC------------------AGTGGAATCAGTTTGCTCTG : 5073"
        out_dna = " 5030 : ACACGTGGAATTACTAAAATCC------------------AGTGGAATCAGTTTGCTCTG : 5073"
        self.assertEqual(out_dna, clear_hashed_bases(testprot, testdna))

    def test_remove_non_letter_signs(self):
        test_string = "AdfA%y*8l#+-w"
        out_string = "AdfAy*lw"
        self.assertEqual(remove_non_letter_signs(test_string), out_string)


class TestDataProviderObject(unittest.TestCase):

    data_base = run_GenePS.DataProviderObject(test_data + "/run_geneps")

    def test_all_attributes(self):
        self.assertEqual(self.data_base.cluster_scope, 5)
        self.assertEqual(self.data_base.group_names[0], "next_best_blast_eef")
        self.assertEqual(self.data_base.group_to_single_copy_ortholog["next_best_blast_eef"], "True")
        self.assertEqual(self.data_base.group_to_group_size["next_best_blast_eef"], 5)
        self.assertEqual(self.data_base.group_by_cluster_to_consensus_length["next_best_blast_eef"]["eef_3.5"], 838)
        self.assertTrue(self.data_base.group_by_cluster_to_hmm["next_best_blast_eef"]["eef_3.5"].split(".")[-1] == "hmmGenePS")
        self.assertEqual(self.data_base.group_by_cluster_to_score_cutoff["next_best_blast_eef"]["eef_3.5"], 1.1826699418473992)
        self.assertEqual(self.data_base.group_by_cluster_to_len_confidence["next_best_blast_eef"]["eef_3.5"], ['303.97390927100224', '913.1505146921315'])


class TestOverseerPredictAllRegions(unittest.TestCase):
    elegans_db = os.path.join(test_data, "databases/c_elegans.PRJNA13758.WS254.genomic.fa")
    overseer = run_GenePS.Overseer("c_elegans", elegans_db, "prediction_location")
    blast = read_blast_output(os.path.join(test_data, "run_geneps/elegans_blast_eef_nextBest.txt"), elegans_db)
    blast.infer_regions()
    overseer.group_to_blast_obj["next_best_blast_eef"] = blast
    overseer.merged_regions = 21
    run_GenePS.data_base = run_GenePS.DataProviderObject(test_data + "/run_geneps")
    with tempdir() as tmp:
        valid_predictions = overseer.predict_on_all_regions(tmp)

    def test_predict_all_regions_number_valid_predictions(self):
        self.assertEqual(4, self.valid_predictions)
        self.assertEqual(17, self.overseer.filter_count)

    def test_predict_all_regions_correct_cluster(self):
        cluster_list = ["OrthologousGroups_I3.5.OG0000365.txt", "OrthologousGroups_I3.5.OG0002137.txt", "OrthologousGroups_I3.5.OG0000365.txt", "eef_3.5"]
        cluster_test_list = []
        for cluster in self.overseer.group_by_cluster_by_contig_to_valid_prediction["next_best_blast_eef"]:
            cluster_test_list.append(cluster)
        self.assertSetEqual(set(cluster_test_list), set(cluster_list))

    def test_predict_all_regions_minus_strand_prediction(self):
        check_list = ['III', 1191666, 1185388, '-']
        for pred in self.overseer.group_by_cluster_by_contig_to_valid_prediction["next_best_blast_eef"]["OrthologousGroups_I3.5.OG0000365.txt"]["III"]:
            self.assertListEqual(check_list, [pred.contig, pred.gene_start, pred.gene_end, pred.strand])

    def test_predict_all_regions_plus_strand_prediction(self):
        check_list = ['II', 13119742, 13121820, '+']
        for pred in self.overseer.group_by_cluster_by_contig_to_valid_prediction["next_best_blast_eef"]["OrthologousGroups_I3.5.OG0000365.txt"]["II"]:
            self.assertListEqual(check_list, [pred.contig, pred.gene_start, pred.gene_end, pred.strand])

    def test_predict_all_regions_no_exonerate_fails(self):
        fail_count = 0
        for group, cluster in self.overseer.group_by_cluster_no_prediction_region.items():
            for pred_list in self.overseer.group_by_cluster_no_prediction_region[group][cluster]:
                fail_count += len(pred_list)
        self.assertEqual(fail_count, 0)


class TestOverlapp(unittest.TestCase):

    def setUp(self, size=5):
        self.overseer = run_GenePS.Overseer("test", "fake", "wayne")
        names = ["zero", "one", "two", "three", "four"]
        for score in range(0, size):
            pred_obj = run_GenePS.PredictionObject("test", "fake", "wayne", "still_test", names[score], "fake2")
            pred_obj.score = score
            pred_obj.contig = "I"
            self.overseer.group_by_contig_to_passed_prediction_list["group"]["I"].append(pred_obj)

    def overlapp_side_effect(self, pred):
        overlapping_pred = ["zero", "two", "four"]
        overlapp_flag = False
        if self.cluster != pred.cluster:
            if pred.cluster in overlapping_pred:
                if self.cluster in overlapping_pred:
                    overlapp_flag = True
        return overlapp_flag

    @patch('run_GenePS.PredictionObject.check_for_overlapp_if_same_contig', return_value=True)
    def test_all_overlapp(self, overlapp):
        count_overlapp = self.overseer.contigwise_overlapping_control("group", "I")
        for cluster in self.overseer.group_by_cluster_by_contig_to_valid_prediction["group"]:
            for pred_obj in self.overseer.group_by_cluster_by_contig_to_valid_prediction["group"][cluster]["I"]:
                highest_value = pred_obj.score
                self.assertEqual(highest_value, 4)
        self.assertEqual(count_overlapp, 5)

    @patch('run_GenePS.PredictionObject.check_for_overlapp_if_same_contig', return_value=False)
    @patch("builtins.print", autospec=True, return_value=None)
    def test_no_overlapp(self, printfunction, overlapp):
        count_overlapp = self.overseer.contigwise_overlapping_control("group", "I")
        count_valid = 0
        for cluster in self.overseer.group_by_cluster_by_contig_to_valid_prediction["group"]:
            count_valid += len(self.overseer.group_by_cluster_by_contig_to_valid_prediction["group"][cluster]["I"])
        self.assertEqual(count_valid, 5)
        self.assertEqual(count_overlapp, 0)

    @patch('run_GenePS.PredictionObject.check_for_overlapp_if_same_contig', autospec=True, side_effect=overlapp_side_effect)
    @patch("builtins.print", autospec=True, return_value=None)
    def test_three_overlapp(self, print_function, overlapp):
        count_valid = 0
        count_overlapp = self.overseer.contigwise_overlapping_control("group", "I")
        self.assertEqual(2, count_overlapp)     # 2 because although 3 are overlapping, just two regions get filtered
        for cluster in self.overseer.group_by_cluster_by_contig_to_valid_prediction["group"]:
            count_valid += len(self.overseer.group_by_cluster_by_contig_to_valid_prediction["group"][cluster]["I"])
        self.assertEqual(3, count_valid)

    @patch('run_GenePS.PredictionObject.check_for_overlapp_if_same_contig', autospec=True, side_effect=overlapp_side_effect)
    @patch("builtins.print", autospec=True, return_value=None)
    def test_three_overlapp_correct_pick(self,  print_function, overlapp):
        self.overseer.contigwise_overlapping_control("group", "I")
        winner_names = ["four", "three", "one"]
        for cluster in self.overseer.group_by_cluster_by_contig_to_valid_prediction["group"]:
            for obj in self.overseer.group_by_cluster_by_contig_to_valid_prediction["group"][cluster]["I"]:
                self.assertIn(obj.cluster, winner_names)

        #Region(contig='III', s_start=8424615, s_end=8434740, strand='+', chunk_cov=100, query_cov=5, q_len=968)

if __name__ == '__main__':
    unittest.main()
