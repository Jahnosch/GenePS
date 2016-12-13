#!/usr/bin/env python3
import unittest
import os
from unittest.mock import patch
import make_Datasets
import tempfile as tmp
import compute_msa
from run_command import tempdir
from compute_msa import msa_operations, MsaObject
from unittest.mock import mock_open

script_path = os.path.dirname(os.path.realpath(__file__))
test_data = os.path.join(script_path, "test_data/compile_script")
file_name = "eef_3.5_no_cea.fa"


class TestHashFasta(unittest.TestCase):

    fa_hash = make_Datasets.hash_fasta(os.path.join(test_data, file_name))

    @patch("make_Datasets.hash_fasta", return_value={">a": 1, ">b": 2})
    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_check_and_hash_fasta(self, logger, fasta_hash):
        fasta_hash = make_Datasets.check_and_hash_fasta(os.path.join(test_data, file_name), "filename")
        self.assertEqual(fasta_hash, None)

    @patch("make_Datasets.hash_fasta", return_value={})
    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_check_and_hash_fasta(self, logger, fasta_hash):
        fasta_hash = make_Datasets.check_and_hash_fasta(os.path.join(test_data, file_name), "filename")
        self.assertEqual(fasta_hash, None)

    @patch("make_Datasets.hash_fasta", return_value={">a": 1, ">b": 2, ">c": 3})
    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_check_and_hash_fasta(self, logger, fasta_hash):
        fasta_hash = make_Datasets.check_and_hash_fasta(os.path.join(test_data, file_name), "filename")
        self.assertEqual(fasta_hash, {">a": 1, ">b": 2, ">c": 3})

    def test_hash_fasta_equal_number_header_seq(self):
        number_seq = len(self.fa_hash.keys())
        number_header = len(self.fa_hash.values())
        self.assertEqual(number_header, number_seq)

    def test_hash_fasta_number_of_entries(self):
        number_entries = len(self.fa_hash.keys())
        self.assertEqual(number_entries, 190)

    def test_hash_fasta_correct_seq_length(self):
        single_seq = len(self.fa_hash[">AMELL.GB42352-PA"][0])
        self.assertEqual(single_seq, 942)


class TestCleanFastaHash(unittest.TestCase):

    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_remove_two_keys(self, logger):
        clean_hash = make_Datasets.clean_fasta_hash({">a": 1, ">b": 2, ">c": 3}, [">c"])
        self.assertEqual(clean_hash, {">c": 3})

    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_key_not_in_hash(self, logger):
        clean_hash = make_Datasets.clean_fasta_hash({">a": 1, ">b": 2, ">c": 3}, [">t", ">p"])
        self.assertEqual(clean_hash, {})

    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_empty_dictionary(self, logger):
        clean_hash = make_Datasets.clean_fasta_hash({}, [">t", ">p"])
        self.assertEqual(clean_hash, {})


class TestWriteHashToFasta(unittest.TestCase):

    def test_fasta_dictionary(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"})
            self.assertSetEqual(set(open(file_path).readlines()), {'>2\n', 'ADFAT\n', '>1\n', 'ATGSAD\n'})

    def test_empty_hash(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {})
            self.assertEqual(file_path, None)

    def test_rigth_optional_formatting(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"}, line_style="{}\t{}\t")
            self.assertSetEqual(set(open(file_path).readlines()[0].split("\t")), {'>2', 'ADFAT', '>1', 'ATGSAD', ''})

    def test_formatting_without_dictionary_value(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"}, line_style="{}")
            self.assertSetEqual(set(open(file_path).readlines()[0].split(">")), {'', '2', '1'})

    def test_formatting_None_if_more_then_two_brackets(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1" : "ATGSAD", ">2": "ADFAT"}, line_style="{}\n{}\n{}\n")
            self.assertEqual(file_path, None)

    def test_value_is_a_list(self):
        with tmp.NamedTemporaryFile() as hash_tmp:
            file_path = make_Datasets.write_hash_to_fasta(hash_tmp.name, {">1": ["ATGSAD"], ">2": "ADFAT"})
            self.assertSetEqual(set(open(file_path).readlines()), {'>2\n', 'ADFAT\n', '>1\n', 'ATGSAD\n'})


class TestIntersectionPoint(unittest.TestCase):

    def test_no_point_smaller_tp_mean(self):
        interpoint = make_Datasets.find_intersection(0.585692592593, 0.397771428571, 0.324498504734, 0.21228085651)
        self.assertEqual(interpoint, 0.729)

    def test_extremly_high_point_in_list(self):
        interpoint = make_Datasets.find_intersection(1.10492361111, 0.449, 0.232010418823, 0.327981237737)
        self.assertEqual(interpoint, 0.874)

    def test_best_point_pick_between_too_large_and_perfect(self):
        interpoint = make_Datasets.find_intersection(1.31034831461, 0.620403669725, 0.261504825592, 0.268836948615)
        self.assertEqual(interpoint, 1.038)

    def test_identical_distributions(self):
        interpoint = make_Datasets.find_intersection(500, 500, 300, 300)
        self.assertEqual(interpoint, None)


class TestLengthBinnedFasta(unittest.TestCase):

    def test_estimate_bin_number_never_below_amount_data(self):
        self.assertEqual(make_Datasets.estimate_bin_number(4, 3, 2, forced_bins=3), 2)

    def test_estimate_bin_force_over_guess(self):
        self.assertEqual(make_Datasets.estimate_bin_number(5, 3, 10, forced_bins=4), 4)

    def test_estimate_bin_emit_guess_if_possible(self):
        self.assertEqual(make_Datasets.estimate_bin_number(5, 3, 10), 5)

    def test_estimate_bin_minimum_before_guess(self):
        self.assertEqual(make_Datasets.estimate_bin_number(5, 6, 10), 6)

    def test_bin_sequence_length_into_three(self):
        bin_dict = make_Datasets.bin_sequence_lengths({"0": 0, "4": 4, "12": 12, "16a": 16, "16b": 16, "18": 18, "24": 24, "26": 26, "28": 28}, 3)
        self.assertEqual(bin_dict, {1: '4', 2: '18', 3: '28'})

    def test_huge_window_size_reduces_bin_number(self):
        self.assertEqual(make_Datasets.bin_sequence_lengths({"0": 0, "4": 4, "12": 12, "16a": 16, "16b": 16, "18": 18, "24": 24, "26": 26, "100": 100}, 3), {1: '26', 3: '100'})

    def test_write_length_binned_fasta_correct_length_hash(self):
        with tmp.NamedTemporaryFile() as tempf:
            length_dict = make_Datasets.write_length_binned_fasta({">first" : ["ADFADF"], ">second" : ["DFAF"]}, "test", tempf.name)
            self.assertEqual(length_dict, {'>second': 4, '>first': 6})

    def test_write_length_binned_fasta_raise_error_if_seq_not_list(self):
        with tmp.NamedTemporaryFile() as tempf:
            self.assertRaises(AttributeError, make_Datasets.write_length_binned_fasta, {">first" : "ADFADF", ">second" : "DFAF"}, "test", tempf.name)

    def test_write_length_binned_fasta_converts_from_bin_to_header_into_bin_to_sequence_correctly(self):
        with tmp.NamedTemporaryFile() as tempf:
            length_dict = make_Datasets.write_length_binned_fasta({">first" : ["ADFADF"], ">second" : ["DFAF"]}, "test", tempf.name)
            self.assertSetEqual(set(open(tempf.name).readlines()), {'>1_test\n', 'DFAF\n', '>2_test\n', 'ADFADF\n'})


class TestGetpHmmScores(unittest.TestCase):

    hmm_search_score_file = os.path.join(test_data, "domtblout.txt")

    def test_hmmsearch_score_dict_length(self):
        score_dict = make_Datasets.parse_hmmer_domain_table(self.hmm_search_score_file)
        self.assertEqual(len(score_dict), 190)

    def test_hmmsearch_one_domain_score(self):
        score_dict = make_Datasets.parse_hmmer_domain_table(self.hmm_search_score_file)
        self.assertEqual(1583, score_dict[">LOLOA1.EN70_1811"])

    def test_hmmsearch_two_domain_score(self):
        score_dict = make_Datasets.parse_hmmer_domain_table(self.hmm_search_score_file)
        self.assertEqual(667, score_dict[">SSCAP.L892_g29260.t1"])


class TestScoreObject(unittest.TestCase):

    def setUp(self, size=5):
        self.fasta_dict = {">ACRAS.cds.Contig10403m.5077": ["QTFRRIVENINVIIATYGDDDGPMGPIMVDPALGNVGFGSGLHGWAFTLKQFAEMYASKFGVQVDKLMKNLWGDRFFNMKTKKRSTSQEDGAVRGFTQFVLDPIFKVF*"],
                  ">ALUMB.ALUE_0000951001-mRNA-1": ["XVCVQTETVLRQAIAERIKPVLFMNKMDRALLELQLGQEELYQTFQRIVENTNVIIATYGDDDGPMGQIMVDPAIGNVGFGSGLHGWAFTLKQFAEMYSEKFGVQ"],
                  ">ACRAS.cds.Contig3658m.1561": ["MEMLYEGPHDDEVAVAIKNCDPNGPLMMYVSKMVPTSDKGRFYAFGRVFSGKVATGMKARIQGPNYVPGKKEDLYEKTIQRTILMMGRYVEPIEDIPSGNIAGLVGVDQYLVKGGTITTFKDAHNLRVMKFSVSPVVRVAVEPKNAGDLPKLVEGLKRLSKSDPMV"]}
        self.score_obj = make_Datasets.ScoreObject(self.fasta_dict, os.path.join(test_data, "eef.hmm"))

    def test_iterative_scoring_hand_checked_score(self):
        score_dict = self.score_obj.iterative_score_computation()
        self.assertEqual(score_dict[">ALUMB.ALUE_0000951001-mRNA-1"], 17)
        self.assertNotIn(">ACRAS.cds.Contig3658m.1561", score_dict)

    def test_bulk_scoring_hand_checked_score(self):
        score_dict = self.score_obj.bulk_score_computation()
        self.assertEqual(score_dict[">ALUMB.ALUE_0000951001-mRNA-1"], 190)

    def test_length_range(self):
        length_dict = {">ACRAS.cds.Contig10403m.5077": 10, ">ALUMB.ALUE_0000951001-mRNA-1": 20, ">ACRAS.cds.Contig3658m.1561": 30}
        self.assertListEqual(list(make_Datasets.calculate_length_range(length_dict)), [10.0, 30.0])

    def test_score_distribution_return_minimum(self):
        self.score_obj.score_dict = {"4": 4, "12": 12, "16": 16, "18": 18, "24": 24, "26": 26}
        self.assertEqual(4, self.score_obj.calculate_score_distribution_parameters())

    def test_score_distribution_return_intersection(self):
        self.score_obj.score_dict = {"4": 4, "12": 12, "16": 16, "18": 18, "24": 24, "26": 26}
        self.assertEqual(13.35, self.score_obj.calculate_score_distribution_parameters(true_negative_scores=[2, 3, 9, 10, 5, 11]))

    def test_score_distribution_too_less_values(self):
        self.score_obj.score_dict = {"4": 4, "12": 12, "16": 16, "18": 18, "24": 24, "26": 26}
        self.assertEqual(4, self.score_obj.calculate_score_distribution_parameters(true_negative_scores=[2, 3, 9]))


class TestMsaObject(unittest.TestCase):

    msa_path = open(os.path.join(test_data, "eef.aln"))
    msa_list = msa_path.readlines()
    msa_path.close()

    @patch("compute_msa.run_cmd", return_value=msa_list)
    def test_msa_operations(self, run_cmd_mock):
        msa_list = compute_msa.msa_operations("blaaa")
        self.assertEqual(len(msa_list)/2, 190)


class Overseer(unittest.TestCase):

    @patch("make_Datasets.output_dir", return_value="")
    @patch("builtins.print", autospec=True, return_value=None)
    def test_initalize_with_single_file(self, printfunction, outdir):
        overseer_obj = make_Datasets.Overseer(os.path.join(test_data, file_name))
        overseer_obj.initialize_input_data()
        self.assertEqual(overseer_obj.group_by_file_to_filepath["compile_script"]["eef_3.5_no_cea"], "/home/jgravemeyer/Dropbox/MSc_project/src/GenePS/test_data/compile_script/eef_3.5_no_cea.fa")

    @patch("make_Datasets.output_dir", return_value="")
    @patch("builtins.print", autospec=True, return_value=None)
    @patch("make_Datasets.logger_Filtered", return_value=None)
    def test_initalize_with_directory(self, logger, printfunction, outdir):
        overseer_obj = make_Datasets.Overseer(test_data)
        overseer_obj.initialize_input_data()
        self.assertEqual(overseer_obj.group_by_file_to_filepath["compile_script"]["eef_3.5_no_cea"], "/home/jgravemeyer/Dropbox/MSc_project/src/GenePS/test_data/compile_script/eef_3.5_no_cea.fa")
        self.assertNotIn("eef.hmm", overseer_obj.group_by_file_to_filepath["compile_script"])

    def test_filter_file_list_remove_two(self):
        overseer_obj = make_Datasets.Overseer("test")
        overseer_obj.group_to_file_list["testgroup"] = ["1", "2", "3", "4"]
        overseer_obj.valid_input_scope = 4
        valid_files = overseer_obj.remove_filtered_files({"testgroup" : ["1", "2"]})
        self.assertListEqual(overseer_obj.group_to_file_list["testgroup"], ['3', '4'])
        self.assertEqual(4-2, valid_files)

    def test_filter_file_list_nothing_to_remove(self):
        overseer_obj = make_Datasets.Overseer("test")
        overseer_obj.group_to_file_list["testgroup"] = ["1", "2", "3", "4"]
        overseer_obj.valid_input_scope = 4
        valid_files = overseer_obj.remove_filtered_files({})
        self.assertListEqual(overseer_obj.group_to_file_list["testgroup"], ["1", "2", "3", "4"])
        self.assertEqual(4, valid_files)
'''
    @patch("make_Datasets.output_dir", return_value="")
    @patch("make_Datasets.generate_hmm", return_value="Test_HMM")
    @patch("make_Datasets.write_length_binned_fasta", return_value={"a": 1})
    @patch("make_Datasets.print_progress", return_value="")
    def test_output_hmm_fasta_correct_file_writing_after_filtering(self, progress_print, binned_fasta, hmmbuild, outputdir):
        overseer_obj = make_Datasets.Overseer(os.path.join(test_data, file_name))
        overseer_obj.initialize_input_data()
        print(overseer_obj.group_by_file_to_filepath)
        print(overseer_obj.group_to_file_list)
        print(len(overseer_obj.group_by_file_to_cluster_hash))
        with tempdir() as msa_dir:
            overseer_obj.output_HMM_and_fasta(msa_dir)
            print(overseer_obj.group_by_file_to_msa_obj['compile_script']['eef_3.5_no_cea'].size_history)
        # check if same path holds later on a fasta file and now aln
            # aln afer realining same size as history as well as the length of the purged fasta hash
'''
# 190, 129
'''
def remove_filtered_files(self, removal_group_to_file_list):
    for group, file_list in removal_group_to_file_list.items():
        for file_name in file_list:
            self.group_to_file_list[group].remove(file_name)
            self.valid_input_scope -= 1
    return self.valid_input_scope
'''

if __name__ == '__main__':
    unittest.main()



