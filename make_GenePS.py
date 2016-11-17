#!/usr/bin/env python3
###############
# MAKE GenePS !
###############

"""
Usage: make_GenePS.py                         -i <DIR> -o <DIR> [-t <FILE>] [--keep] [--sco]

    Options:
        -h, --help                            show this screen.

        General
        -i, --input <DIR>                     either single input file or directory with files or sub-folders
        -o, --output <DIR>                    directory where to save the output file(s)
        -t, --translation_files <FILE>        tab separated file in style of: dir_path blast <line> file_path speciesID <line> file_path protein.fa <line> file_path sequenceID
        --keep                                command to safe intermediate files
        --sco                                 only if all input clusters are considered to be single copy orthologs
"""

import os
import sys
import tempfile as tmp
from collections import defaultdict
from run_command import run_cmd, tempdir, check_programs
from compute_msa import generate_msa, MsaObject
import_errors = []

try:
    from docopt import docopt
except ImportError:
    import_errors.append("[ERROR] : Module \'Docopt\' was not found. Please install \'Docopt\' using \'pip install docopt\'")
try:
    from numpy import roots, log, mean, std
except ImportError:
    import_errors.append("[ERROR] : Module \'numpy\' was not found. Please install \'numpy\' using \'pip install numpy\'")
if import_errors:
    sys.exit("\n".join(import_errors))

########################################################################################################################
# Global Functions: generic functions
########################################################################################################################


def hash_fasta(fasta_file):
    fasta = {}
    active_sequence_name = ""
    with open(fasta_file) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line.split(" ")[0]
                if active_sequence_name not in fasta:
                    fasta[active_sequence_name] = []
                continue
            sequence = line
            fasta[active_sequence_name].append(sequence)
    file_one.close()
    return fasta


def check_and_hash_fasta(fasta_file, file_name):
    fasta_dict = hash_fasta(fasta_file)
    if len(fasta_dict) <= 2:
        print("\t[!] INPUT ERROR: '{}' - has less than 2 sequences\n".format(file_name))
        return None
    else:
        return fasta_dict


def clean_fasta_hash(old_hash, to_keep_list):
    clean_dict = {}
    for header in to_keep_list:
        clean_dict[header] = old_hash[header]
    return clean_dict


def write_hash_to_fasta(hash_fa, location):
    string_list = []
    with open(location, "w") as hash_f:
        for header in hash_fa:
            string_list.append(header)
            string_list.append("".join(hash_fa[header]))
        hash_f.write("\n".join(string_list))
    return location


def get_outdir(out_directory, add_dir=""):
    if type(out_directory) != str:
        print("\t[!] {} is NOT a directory!\n".format(out_directory))
        print("\t[!] Please specify an output directory\n")
        sys.exit()
    elif os.path.isfile(out_directory):
        print("\t[!] {} is a File\n".format(out_directory))
        print("\t[!] Please specify an output directory\n")
        sys.exit()
    elif not os.path.exists(os.path.join(out_directory, add_dir)):
        os.mkdir(os.path.join(out_directory, add_dir))
        return os.path.abspath(os.path.join(out_directory, add_dir))
    else:
        return os.path.abspath(os.path.join(out_directory, add_dir))


def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    format_str = "{0:." + str(decimals) + "f}"
    percents = format_str.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.write('\r')
    if iteration == total:
        sys.stdout.write("\n")
    sys.stdout.flush()


def write_to_tempfile(tmp_name, string):
    new_file = open(tmp_name, "w")
    new_file.write(string)
    new_file.seek(0)
    new_file.close()


def find_intersection(m_p, m_n, std_p, std_n):
    a = 1/(2 * std_p ** 2) - 1 / (2 * std_n ** 2)
    b = m_n / (std_n ** 2) - m_p / (std_p ** 2)
    c = m_p ** 2 / (2 * std_p ** 2) - m_n ** 2 / (2 * std_n ** 2) - log(std_n / (std_p * 2))
    intersection = roots([a, b, c])
    sorted_no_negative = sorted([x for x in intersection if not x < 0], reverse=True)
    for point in sorted_no_negative:
        if point < m_p:
            return point


########################################################################################################################
# Global Functions: HMM related functions
########################################################################################################################


def generate_hmm(hmm_path, msa_path):
    command = ["hmmbuild", hmm_path, msa_path]
    run_cmd(command=command, wait=True)
    return hmm_path


def get_phmm_score(hmm_file, query_file, cluster_dict=None):
    command = ["hmmsearch", "--noali", hmm_file, query_file]
    read_counter = 0
    name_to_score = {}
    for line in run_cmd(command=command, wait=False):
        if "E-value" in line or read_counter == 1:
            read_counter += 1
        elif read_counter == 2:
            line = line.strip("\n").split()
            if len(line) > 0:
                try:
                    sequence_header = ">" + line[8].strip()
                    if cluster_dict:
                        length = cluster_dict[sequence_header]
                        name_to_score[sequence_header] = (round(float(line[1])) / length)
                    else:
                        name_to_score[sequence_header] = (round(float(line[1])))
                except ValueError:
                    print("VALUE ERROR")
                    print(line)
                    return name_to_score
            else:
                return name_to_score


def get_consensus(hmm_file):
    command = "hmmemit -c " + hmm_file
    cons_list = []
    for line in run_cmd(command=command, wait=False):
        if not line.startswith(">"):
            cons_list.append(line.strip("\n"))
    return "".join(cons_list)


########################################################################################################################
# Global Functions: True Negative Score computation
########################################################################################################################

def parse_true_negative_arg(tn_arg):
    arg_dict = {}
    should_be = ["blast", "sequenceID", "speciesID", "protein"]
    error_list = []
    with open(tn_arg) as file_dict:
        for line in file_dict:
            line = line.strip("\n").split("\t")
            arg_dict[line[1]] = line[0].strip()
    for argument in should_be:
        if argument not in arg_dict:
            error_list.append("[!] Could not find '{}' - please check format of directory file".format(argument))
        else:
            if not os.path.exists(arg_dict[argument]):
                error_list.append("[!] specified {} - directory does not exist!".format(argument))
    if error_list:
        print("\n".join(error_list))
        sys.exit()
    else:
        return arg_dict


def get_blast_files(input_dir):
    """returns the directory of the folder storing the blast files and a set of all blast file names"""
    dir_list = []
    if not os.path.isdir(input_dir):
        print(input_dir, " : has to be directory")
        sys.exit()
    else:
        for subdir, dirs, files in os.walk(input_dir):
            for single_file in files:
                if not single_file.startswith("Blast"):
                    print("\t[!] not only blast files found in: ", input_dir)
                    sys.exit()
                dir_list.append(single_file)
                if not dir_list:
                    print("\t[!]Blast files could either not be read or found in: ", input_dir)
                    sys.exit()
            return os.path.abspath(subdir), set(dir_list)


def hash_blast_files(species_id_set):
    for species_id in species_id_set:
        b_file_name = "Blast{}_{}.txt".format(species_id, species_id)
        b_file_path = os.path.join(blast_path, b_file_name)
        blast_specID_protID_hitList[species_id] = blast_prot_id_hits(b_file_path)
    return len(blast_specID_protID_hitList)


def hash_sequence_translation_file(sequence_id_name_file):
    """ hashes translation file and returns hash-list in style of: [species_protein_name] = species_protein_id
     as well as an reversed hash: [species_protein_id] = species_protein_name """
    translation_dict = {}
    translation_dict_rev = {}
    with open(sequence_id_name_file) as tf:
        for line in tf:
            line = line.strip("\n").split(":")
            prot_sp_name = line[1].strip()
            prot_sp_id = line[0]
            translation_dict[prot_sp_name] = prot_sp_id
            translation_dict_rev[prot_sp_id] = prot_sp_name
    return translation_dict, translation_dict_rev


def hash_species_translation_file(species_id_file):
    id_list = []
    with open(species_id_file) as tf:
        for line in tf:
            line = line.strip("\n").split(":")
            id_list.append(line[0])
    return set(id_list)


def blast_prot_id_hits(blast_file):
    """ hashes within species blast and returns hash in style of: [protein_id] = list(hits) """
    protid_hits_hash = defaultdict(list)
    with open(blast_file) as blast_f:
        for line in blast_f:
            line = line.split("\t")
            spec_id, prot_id = line[0].split("_")   # line[0] = self, line[1] = other
            protid_hits_hash[prot_id].append(line[1].split("_")[1])
    return protid_hits_hash


def hash_spec_id_to_prot_id_list(header_list):
    """ hashes cluster file and returns hash-list in style of: [species_id] = list(protein_id)
    needs idpair to name pair as translation dict"""
    cluster_dict = defaultdict(list)
    for header in header_list:
        if header.startswith(">"):
            header = header.strip("\n")
            protein_sp_name = header.strip(">")
            try:
                protein_sp_id = idPair_2_namePair[protein_sp_name]
            except KeyError:
                return None
            sp_id, protein_id = protein_sp_id.split("_")
            cluster_dict[sp_id].append(protein_id)
    return cluster_dict


def next_best_blast_hit(forbidden, self_hits):
    """forbidden = list of all protein ids of one species in the same cluster
    self_hits = list of all the hits from blasting a protein within a species"""
    for hit in self_hits:
        if hit not in forbidden:
            return hit
    return None


def get_twin_hash(cluster_specID_to_protList, speciesID_to_forbidden_proteins):
    """generates a True negative fasta hash of a cluster_file by searching for each proteins next best blast hit"""
    fasta_hash = {}
    for species_id, prot_list in cluster_specID_to_protList.items():
        for prot_id in prot_list:
            # two proteins of the same species are allowed to have the same next best
            protID_selfHits = blast_specID_protID_hitList[species_id][prot_id]
            next_best = next_best_blast_hit(speciesID_to_forbidden_proteins[species_id], protID_selfHits)
            if next_best is not None:
                idx_pair = "{}_{}".format(species_id, next_best)
                name_pair = ">" + namePair_2_idPair[idx_pair]
                try:
                    sequence = all_protein_fasta_dict[name_pair][0]
                except KeyError:
                    log_not_in_protein_file.append(name_pair)
                    continue
                fasta_hash[name_pair] = [sequence]
    return fasta_hash


########################################################################################################################
# Class ScoreObject - Positive Scores
########################################################################################################################

class ScoreObject:
    def __init__(self, fasta_dict, hmm_path):
        self.fasta_hash = fasta_dict
        self.hmm_path = hmm_path
        self.score_dict = {}
        self.score_distribution_parameters = None
        self.length_dict = {}
        self.length_distribution_parameters = None

    def query_for_fasta(self, query):
        query = query.split()[0]
        return query + "\n" + "".join(self.fasta_hash[query])

    def seq_length_hash(self):
        for header in self.fasta_hash:
            self.length_dict[header] = len(self.fasta_hash[header][0])

    def calculate_length_distribution_parameters(self):
        std_len = std(list(self.length_dict.values()))
        average_len = mean(list(self.length_dict.values()))
        self.length_distribution_parameters = (average_len - std_len, average_len + std_len)
        return self.length_distribution_parameters

    def calculate_score_distribution_parameters(self, true_negative_scores=None):
        tp_std = std(list(self.score_dict.values()))
        tp_average = mean(list(self.score_dict.values()))
        if true_negative_scores:
            if len(true_negative_scores) >= 4:
                tn_std = std(true_negative_scores)
                tn_average = mean(true_negative_scores)
                self.score_distribution_parameters = find_intersection(tp_average, tn_average, tp_std, tn_std)
            else:
                if tp_std >= (tp_average/3):
                    self.score_distribution_parameters = tp_average - tp_std
                else:
                    self.score_distribution_parameters = tp_average - (2 * tp_std)
        else:
            self.score_distribution_parameters = (tp_average, tp_std)
        return self.score_distribution_parameters

    def generate_msa_string(self, rest_prot):
        seq_list = []
        for header in rest_prot:
            header = header.split()[0]
            seq_list.append(header)
            seq_list.append("".join(self.fasta_hash[header]))
        with tmp.NamedTemporaryFile() as r_tmp:
            write_to_tempfile(r_tmp.name, "\n".join(seq_list))
            list_msa = generate_msa(r_tmp.name)
        return "\n".join(list_msa)

    def iterative_score_computation(self, length_normalized=False):
        if length_normalized:
            self.seq_length_hash()
        for idx in range(0, len(self.fasta_hash.keys())):
            rest_prot = list(self.fasta_hash.keys())
            query = rest_prot.pop(idx)
            with tmp.NamedTemporaryFile() as q_tmp:
                write_to_tempfile(q_tmp.name, self.query_for_fasta(query))
                msa_string = self.generate_msa_string(rest_prot)
                with tmp.NamedTemporaryFile() as msa_tmp:
                    write_to_tempfile(msa_tmp.name, msa_string)
                    with tmp.NamedTemporaryFile() as hmm_tmp:
                        generate_hmm(hmm_tmp.name, msa_tmp.name)
                        try:
                            if length_normalized is True:
                                score_dict = get_phmm_score(hmm_tmp.name, q_tmp.name, self.length_dict)
                            else:
                                score_dict = get_phmm_score(hmm_tmp.name, q_tmp.name)
                        except IndexError:
                            continue
            self.score_dict.update(score_dict)
        return self.score_dict

    # not needed while re-aligned
    def compute_full_phmm(self, location):
        msa_string = self.generate_msa_string(self.fasta_hash.keys())
        with tmp.NamedTemporaryFile() as msa_tmp:
            write_to_tempfile(msa_tmp.name, msa_string)
            generate_hmm(location, msa_tmp.name)
        return location

    def bulk_score_computation(self, length_normalized=False):
        if length_normalized:
            self.seq_length_hash()
        with tmp.NamedTemporaryFile() as q_tmp:
            seq_list = []
            for header in self.fasta_hash.keys():
                seq_list.append(header)
                seq_list.append(self.fasta_hash[header][0])
            write_to_tempfile(q_tmp.name, "\n".join(seq_list))
            if length_normalized:
                self.score_dict = get_phmm_score(self.hmm_path, q_tmp.name, self.length_dict)
            else:
                self.score_dict = get_phmm_score(self.hmm_path, q_tmp.name)
        return self.score_dict


########################################################################################################################
# Class: Overseer - stores / controls the data flow and orchestrates the generation of other Objects
########################################################################################################################

class Overseer:

    def __init__(self, input_dir):
        self.input_dir = input_dir
        self.group_by_file_to_filepath = defaultdict(dict)
        self.group_to_file_list = defaultdict(list)
        self.group_to_result_path = {}
        self.group_by_file_to_cluster_hash = defaultdict(dict)
        self.group_by_file_to_fasta_path = defaultdict(dict)
        self.group_by_file_to_msa_obj = defaultdict(dict)
        self.group_by_file_to_hmm = defaultdict(dict)
        self.group_by_file_to_score_obj = defaultdict(dict)
        self.removed_group_to_file_list = defaultdict(list)
        self.input_scope = 0
        self.filtered_input_scope = 0

        # True negative cluster attributes
        self.species_id_set = None
        self.group_by_file_to_twin_hash = defaultdict(dict)
        self.group_by_file_to_twin_score_obj = defaultdict(dict)
        self.group_by_file_to_twin_hmm = defaultdict(dict)
        self.untranslatable_files = []

    ####################################################################################################################
    # Overseer - Functions to specify input/output directories and hash cluster files
    ####################################################################################################################

    def walk_and_hash_input(self):
        if os.path.isfile(self.input_dir):
            group_name = self.input_dir.split("/")[-2]
            file_name = ".".join(self.input_dir.split("/")[-1].strip().split(".")[0:-1])
            single_file = os.path.abspath(self.input_dir)
            self.group_to_result_path[group_name] = os.path.join(output_dir, group_name + ".GenePS")
            fasta_hash = check_and_hash_fasta(single_file, file_name)
            if fasta_hash is not None:
                self.group_by_file_to_cluster_hash[group_name][file_name] = fasta_hash
                self.group_by_file_to_filepath[group_name][file_name] = single_file
                self.group_to_file_list[group_name].append(file_name)
                self.filtered_input_scope = 1
            else:
                self.removed_group_to_file_list[group_name].append(file_name)
                log_filtered_files.append(file_name)
            self.input_scope = 1
            return self.filtered_input_scope
        else:
            for subdir, dirs, files in os.walk(self.input_dir):
                group_name = subdir.split("/")[-1]
                self.group_to_result_path[group_name] = os.path.join(output_dir, group_name + ".GenePS")
                for single_file in files:
                    file_path = os.path.join(os.path.abspath(subdir), single_file)
                    file_name = ".".join(single_file.strip().split(".")[0:-1])
                    fasta_hash = check_and_hash_fasta(file_path, file_name)
                    if fasta_hash is not None:
                        self.group_by_file_to_cluster_hash[group_name][file_name] = fasta_hash
                        self.group_by_file_to_filepath[group_name][file_name] = file_path
                        self.group_to_file_list[group_name].append(file_name)
                        self.filtered_input_scope += 1
                    else:
                        self.removed_group_to_file_list[group_name].append(file_name)
                        print("\t[!] Not enough entries\n".format(file_name))
                        log_filtered_files.append(file_name)
                    self.input_scope += 1
            return self.filtered_input_scope

    def initialize_input_data(self):
        if not os.path.exists(self.input_dir):
            print("\t[!] FATAL ERROR: '{}' is neither File nor Directory\n".format(self.input_dir))
            sys.exit()
        valid_clusters = self.walk_and_hash_input()
        print("\n\t{} group(s) - containing {} file(s)\n\t{}\n".format(str(len(self.group_by_file_to_filepath)), str(self.input_scope), "-"*33))
        if valid_clusters == 0:
            print("\t[!] FATAL ERROR: No valid input files given\n")
            sys.exit()
        else:
            if self.input_scope - self.filtered_input_scope > 0:
                print("\t[!] {} cluster(s) too small\n".format(self.input_scope - self.filtered_input_scope))
            return self.filtered_input_scope

    ####################################################################################################################
    # Overseer - MSA/HMM generation and normal true negative score computation
    ####################################################################################################################

    def remove_filtered_from_file_list(self):
        for group, file_list in self.removed_group_to_file_list.items():
            for file_name in file_list:
                self.group_to_file_list[group].remove(file_name)
                self.filtered_input_scope -= 1
        return self.filtered_input_scope

    def compute_msa_and_hmm(self, directory):
        count = 1
        for group, file_list in self.group_to_file_list.items():
            for file_name in file_list:
                msa_list = generate_msa(self.group_by_file_to_filepath[group][file_name])
                msa_obj = MsaObject(msa_list, file_name, directory)
                msa_obj.msa_to_fasta()
                msa_obj.trim_remove()
                if msa_obj.check_msa_size_and_length() is True:
                    self.group_by_file_to_cluster_hash[group][file_name] = clean_fasta_hash(self.group_by_file_to_cluster_hash[group][file_name], msa_obj.all_header())
                    self.group_by_file_to_fasta_path[group][file_name] = write_hash_to_fasta(self.group_by_file_to_cluster_hash[group][file_name], os.path.join(output_dir, file_name + ".fa_GenePS"))
                    if msa_obj.size_history[0] != msa_obj.size_history[-1]:
                        msa_obj.re_align(self.group_by_file_to_fasta_path[group][file_name])
                    self.group_by_file_to_msa_obj[group][file_name] = msa_obj
                    self.group_by_file_to_hmm[group][file_name] = generate_hmm(os.path.join(output_dir, file_name + ".hmmGenePS"), msa_obj.path)
                else:
                    self.removed_group_to_file_list[group].append(file_name)
                    log_filtered_files.append(file_name)
                print_progress(count, self.filtered_input_scope, prefix='\tGenerating Hidden Markov Models:\t', suffix='Complete', bar_length=30)
                count += 1
        return self.remove_filtered_from_file_list()

    def choose_scoring_method(self, cluster_hash, group, file_name, normalize=False):
        hmm = self.group_by_file_to_hmm[group][file_name]
        filtered_proteins = self.group_by_file_to_msa_obj[group][file_name].all_header()
        scoring_obj = ScoreObject(cluster_hash, hmm)
        if len(filtered_proteins) < 20:
            score_hash = scoring_obj.iterative_score_computation(length_normalized=normalize)
        else:
            score_hash = scoring_obj.bulk_score_computation(length_normalized=normalize)
        return scoring_obj

    def compute_all_hmm_scores(self, length_normalized=False):
        count = 1
        print("\n")
        for group, file_list in self.group_to_file_list.items():
            for file_name in file_list:
                fasta_hash = self.group_by_file_to_cluster_hash[group][file_name]
                score_obj = self.choose_scoring_method(fasta_hash, group, file_name, normalize=length_normalized)
                if keep:
                    with open(os.path.join(keep_dir, "{}_{}_scores.txt".format(group, file_name)), "w") as score_f:
                        for protein in score_obj.score_dict:
                            score_f.write("{}\t{}\n".format(protein, score_obj.score_dict[protein]))
                self.group_by_file_to_score_obj[group][file_name] = score_obj
                print_progress(count, self.filtered_input_scope, prefix='\tComputing HMM Score Distributions:\t', suffix='Complete', bar_length=30)
                count += 1
        return self.group_by_file_to_score_obj

    ####################################################################################################################
    # Overseer - Functions for generating True negative scores
    ####################################################################################################################

    def make_cluster_specific_true_negativ_hash(self, group, file_name):
        unfiltered_header = self.group_by_file_to_cluster_hash[group][file_name].keys()
        filtered_header = self.group_by_file_to_msa_obj[group][file_name].all_header()
        speciesID_to_forbidden_proteins = hash_spec_id_to_prot_id_list(unfiltered_header)
        if speciesID_to_forbidden_proteins is None:
            print("\t[!] WARNING - omitting True negative computation for '{}'\n".format(file_name))
            self.untranslatable_files.append(file_name)
            log_files_omitted_for_tn.append(file_name)
            return None
        cluster_specID_to_protList = hash_spec_id_to_prot_id_list(filtered_header)
        fasta_hash = get_twin_hash(cluster_specID_to_protList, speciesID_to_forbidden_proteins)
        self.group_by_file_to_twin_hash[group][file_name] = fasta_hash
        return fasta_hash

    def compute_true_negative_hmm_scores(self, length_normalized=False):
        count = 1
        print("\n")
        for group, file_list in self.group_to_file_list.items():
            for file_name in file_list:
                hmm = self.group_by_file_to_hmm[group][file_name]
                fasta_hash = self.make_cluster_specific_true_negativ_hash(group, file_name)
                if fasta_hash is None:
                    continue
                scoring_obj = ScoreObject(fasta_hash, hmm)
                score_hash = scoring_obj.bulk_score_computation(length_normalized=length_normalized)
                if keep:
                    with open(os.path.join(keep_dir, "{}_{}_TrueNegativeScores.txt".format(group, file_name)), "w") as score_f:
                        for protein in score_hash:
                            score_f.write("{}\t{}\n".format(protein, score_hash[protein]))
                self.group_by_file_to_twin_score_obj[group][file_name] = scoring_obj
                self.group_by_file_to_twin_hmm[group][file_name] = scoring_obj.compute_full_phmm(os.path.join(output_dir, file_name + ".TN.hmmGenePS"))
                print_progress(count, self.filtered_input_scope, prefix='\tTrue Negative Score Distributions:\t', suffix='Complete', bar_length=30)
                count += 1
        return self.group_by_file_to_twin_score_obj


########################################################################################################################
# main
########################################################################################################################

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    infile = args['--input']
    output_dir = get_outdir(args['--output'])
    keep = args['--keep']
    true_negative_file = args['--translation_files']
    single_copy_ortholog = args['--sco']
    check_programs("hmmsearch", "hmmemit", "hmmbuild", "mafft", "trimal")

    print("\n{}\n# MAKE GenePS #\n{}\n\nPreparing Files...\n".format("#"*15, "#"*15))
    # if true negative translation files provided
    blast_specID_protID_hitList = defaultdict(lambda: defaultdict(list))
    idPair_2_namePair, namePair_2_idPair = {}, {}
    all_protein_fasta_dict = {}
    if true_negative_file:
        tn_args = parse_true_negative_arg(true_negative_file)
        blast_path, blast_file_set = get_blast_files(tn_args["blast"])
        idPair_2_namePair, namePair_2_idPair = hash_sequence_translation_file(tn_args["sequenceID"])
        all_protein_fasta_dict = hash_fasta(tn_args["protein"])
        species_ids = hash_species_translation_file(tn_args["speciesID"])
        number_blast_files = hash_blast_files(species_ids)
        if not number_blast_files - len(species_ids) == 0:
            print("\t[!] FATAL ERROR: Not all Blast files could be hashed\n")
            sys.exit()
    if keep:
        keep_dir = get_outdir(output_dir, add_dir="intermediate_files")

    # lists for log files
    log_filtered_files = ["\n\n# Files with either to less proteins or which got filtered based on their MSAs:\n\n"]
    log_not_in_protein_file = ["\n\n# Sequences unable to extract from protein file:\n\n"]
    log_untranslatable_species = ["\n\n# Species names unable to translate:\n\n"]
    log_files_omitted_for_tn = ["\n\n# Files containing untranslatable Species: \n\n"]

    # initialize the overseer and hash and check files
    overseer_obj = Overseer(infile)
    filtered_data_scope = overseer_obj.initialize_input_data()

    with tempdir() as temp_dir:
        # msa and hmm for all cluster
        filtered_data_scope = overseer_obj.compute_msa_and_hmm(temp_dir)
        if not filtered_data_scope > 0:
            print("\t[!] FATAL ERROR: NO Multiple Sequence Alignments computable\n")
            sys.exit()
        # score distributions for all cluster
        all_score_hashes = overseer_obj.compute_all_hmm_scores(length_normalized=True)
        # compute true negative scores
        if true_negative_file:
            all_tn_scores_hash = overseer_obj.compute_true_negative_hmm_scores(length_normalized=True)
        # write results
        print("\n")
        read_count = 1
        for name_group, all_files in overseer_obj.group_to_file_list.items():
            with open(overseer_obj.group_to_result_path[name_group], "w") as results_file:
                results_file.write("group: {}\ngroup_size: {}\nsingle_copy_ortholog: {}\n".format(name_group, str(len(all_files)), single_copy_ortholog))
                for cluster_name in all_files:
                    consensus = get_consensus(overseer_obj.group_by_file_to_hmm[name_group][cluster_name])
                    len_mean_std = ["-", "-"]
                    if true_negative_file:
                        tn_scores = list(overseer_obj.group_by_file_to_twin_score_obj[name_group][cluster_name].score_dict.values())
                        if not tn_scores:
                            tn_scores = [0]
                        tp_score_cut_off = overseer_obj.group_by_file_to_score_obj[name_group][cluster_name].calculate_score_distribution_parameters(true_negative_scores=tn_scores)
                        if single_copy_ortholog:
                            len_mean_std = overseer_obj.group_by_file_to_score_obj[name_group][cluster_name].calculate_length_distribution_parameters()
                        distributions_str = "#score_cut_off: {}\n#length_mean_std: {},{}".format(tp_score_cut_off, len_mean_std[0], len_mean_std[1])
                    else:
                        tp_score_mean_std = overseer_obj.group_by_file_to_score_obj[name_group][cluster_name].calculate_score_distribution_parameters()
                        if single_copy_ortholog:
                            len_mean_std = overseer_obj.group_by_file_to_score_obj[name_group][cluster_name].calculate_length_distribution_parameters()
                        distributions_str = "#score_mean_std: {},{}\n#length_mean_std: {},{}".format(tp_score_mean_std[0], tp_score_mean_std[1], len_mean_std[0], len_mean_std[1])
                    results_file.write("#name: {}\n{}\n{}\n".format(cluster_name, distributions_str, consensus))
                    print_progress(read_count, overseer_obj.filtered_input_scope, prefix='\tWriting Results to Files:\t\t', suffix='Complete', bar_length=30)
                    read_count += 1
                results_file.close()

        # writing a log file
        with open(output_dir + "/LOG.txt", "w") as log:
            log.write("\n# {} groups and {} files\n".format(len(overseer_obj.group_to_file_list), str(overseer_obj.input_scope)))
            log.write("\n# {} filtered files\n".format(overseer_obj.filtered_input_scope))
            log.write("\n".join(log_filtered_files))
            log.write("\n".join(log_untranslatable_species))
            log.write("\n".join(log_files_omitted_for_tn))
            log.write("\n".join(log_not_in_protein_file))
    print("\nDONE!\n")
