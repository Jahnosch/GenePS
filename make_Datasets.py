#!/usr/bin/env python3
###############
# MAKE GenePS !
###############

"""
Usage: make_Datasets.py                         -i <DIR> -o <DIR> [-t <FILE>] [--keep]

    Options:
        -h, --help                            show this screen.

        General
        -i, --input <DIR>                     either single input file or directory with files or sub-folders
        -o, --output <DIR>                    directory where to save the output file(s)
        -t, --translation_files <FILE>        tab separated file in style of: dir_path blast <line> file_path speciesID <line> file_path protein.fa <line> file_path sequenceID
        --keep                                command to safe intermediate files
"""

import os
import sys
import tempfile as tmp
import logging
from operator import itemgetter
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
    """takes sting to fasta file and returns a dictionary in style of: >header:[ADFASDF]"""
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
    """wrapper for hash_fasta function to check if hash has sufficient entries"""
    fasta_dict = hash_fasta(fasta_file)
    if len(fasta_dict) <= 2:
        logger_Filtered.warning("INPUT ERROR: '{}' - has less than 2 sequences\n".format(file_name))
        return None
    else:
        return fasta_dict


def clean_fasta_hash(old_hash, key_list):
    """returns a new hash which takes all entries from the old hash if they are present in the key_list"""
    clean_dict = {}
    for header in key_list:
        try:
            clean_dict[header] = old_hash[header]
        except KeyError:
            continue
    return clean_dict


def write_hash_to_fasta(file_path, dictionary, line_style="{}\n{}\n"):
    """writes a hash to a file in the given line_style. Line style default writes key and value to separate lines.
    If value is a list it joins elements without spaces."""
    if not dictionary:
        return None
    try:
        with open(file_path, "w") as hash_file:
            for header, value in dictionary.items():
                if type(value) is list:
                    value = "".join(value)
                hash_file.write(line_style.format(header, value))
        return file_path
    except IndexError or KeyError or AttributeError or TypeError:
        return None


def estimate_bin_number(k_guess, k_min, amount_data, forced_bins=None):
    """takes a guess for k = bin number, a specified 'minimum bin number' and the total amount of data points.
    Returns forced bin or guess if not smaller then min or amount of data points."""
    if k_min > amount_data:
        return amount_data
    elif forced_bins:
        return forced_bins
    elif k_guess < k_min:
        return k_min
    else:
        return k_guess


def bin_sequence_lengths(length_dict, bins=None):
    """takes a dictionary with fasta header as key and sequence lengths as values. k = number of bins, default max/25
    since the smallest protein is a bit above 20.Returns a collapsed dict using bins as keys, pointing to the header
    of the longest sequence representing that bin."""
    highest_member_dict = {}
    min_L, max_L = min(length_dict.values()), max(length_dict.values())
    k = estimate_bin_number(round(max_L / 25), 3, len(length_dict), forced_bins=bins)
    b_width = round((max_L - min_L) / k)
    for header, L in sorted(length_dict.items(), key=itemgetter(1)):
        for bin_x in range(1, k+1):
            if bin_x == 1 and L <= min_L + b_width:
                highest_member_dict[bin_x] = header
            elif bin_x == k and L > min_L + (k * b_width):
                highest_member_dict[bin_x] = header
            elif min_L + ((bin_x-1) * b_width) < L <= min_L + (bin_x * b_width):
                highest_member_dict[bin_x] = header
    return highest_member_dict


def write_length_binned_fasta(fasta_dict, cluster_name, location):
    """takes a dictionary in in style of {'fasta_header':sequence} and outputs
    a fasta file with header in style of {'>cluster_name_bin': sequence} at the specified 'location'."""
    length_hash = {header: len(seq[0]) for header, seq in fasta_dict.items() if type(seq) is list}
    if not length_hash:
        raise AttributeError("! No length hash")
    bin_to_header = bin_sequence_lengths(length_hash)
    bin_to_sequence = {bin_x: fasta_dict[header] for bin_x, header in bin_to_header.items()}
    write_hash_to_fasta(location, bin_to_sequence, ">{}_" + cluster_name + "\n" + "{}\n")
    return length_hash


def get_outdir(out_directory, add_dir=""):
    """generates output directory in case it does not exist."""
    if type(out_directory) != str:
        print("\t[!] {} is NOT a directory! Please specify an output directory\n".format(out_directory))
        sys.exit()
    elif os.path.isfile(out_directory):
        print("\t[!] {} is a File! Please specify an output directory\n".format(out_directory))
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
    """takes tmp_name, opens it as a real file and writes string into it"""
    new_file = open(tmp_name, "w")
    new_file.write(string)
    new_file.seek(0)
    new_file.close()


def find_intersection(m_p, m_n, std_p, std_n):
    """find the intersection between to distributions assuming they are normally distributed. If > 1 insersection points
    it returns the biggest inter.-point which is smaller then the mean of the true positive distribution"""
    a = 1/(2 * std_p ** 2) - 1 / (2 * std_n ** 2)
    b = m_n / (std_n ** 2) - m_p / (std_p ** 2)
    c = m_p ** 2 / (2 * std_p ** 2) - m_n ** 2 / (2 * std_n ** 2) - log(std_n / (std_p * 2))
    intersection = roots([a, b, c])
    sorted_no_negative = sorted([round(x, 3) for x in intersection if not x < 0], reverse=True)
    for point in sorted_no_negative:
        if point < m_p:
            return point
    try:
        return max(sorted_no_negative)
    except ValueError:
        return None


########################################################################################################################
# Global Functions: HMM related functions
########################################################################################################################


def generate_hmm(hmm_path, msa_path):
    """builds hidden markov model from MSA"""
    command = ["hmmbuild", hmm_path, msa_path]
    run_cmd(command=command, wait=True)
    return hmm_path


def get_phmm_score(hmm_file, query_file, header_to_length=None):
    """aligns query-fasta-file against HMM and returns a score hash in style of >header:score.
    If dictionary >header:len(seq) is provided, it directly length-normalizes the score for each protein"""
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
                    if header_to_length:
                        length = header_to_length[sequence_header]
                        name_to_score[sequence_header] = (round(float(line[1])) / length)
                    else:
                        name_to_score[sequence_header] = (round(float(line[1])))
                except ValueError or IndexError:
                    print("VALUE ERROR or IndexError\n", line)
                    return name_to_score
            else:
                return name_to_score


def get_consensus(hmm_file):
    """infers a header-less consensus sequence from HMM"""
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
    """input control for optional translation file (holding 4 other files for generating true negative cluster files)"""
    arg_dict = {}
    needed_arguments = ["blast", "sequenceID", "speciesID", "protein"]
    error_list = []
    with open(tn_arg) as file_dict:
        for line in file_dict:
            line = line.strip("\n").split("\t")
            arg_dict[line[1]] = line[0].strip()
    for argument in needed_arguments:
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
    """returns the directory of the blast files folder and a set of all blast file names"""
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


def hash_single_blast_file(blast_file):
    """ hashes single within-species-blast and returns hash in style of: [protein_id] = list(hits). Hits are
     sorted because they are sorted by E-value in the blast file."""
    protid_hits_hash = defaultdict(list)
    with open(blast_file) as blast_f:
        for line in blast_f:
            line = line.split("\t")
            spec_id, prot_id = line[0].split("_")   # line[0] = self, line[1] = other
            protid_hits_hash[prot_id].append(line[1].split("_")[1])
    return protid_hits_hash


def hash_all_blast_files(species_id_set):
    """opens all blast files for which species_ids could be found. Calls hash_single_blast_file and fills the global blast dictionary in style of
    [species_id][protein_id]=list of all blast hits for this protein id"""
    for species_id in species_id_set:
        b_file_name = "Blast{}_{}.txt".format(species_id, species_id)
        b_file_path = os.path.join(blast_path, b_file_name)
        blast_specID_protID_hitList[species_id] = hash_single_blast_file(b_file_path)
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
    """takes file with all species id and returns them as a set"""
    id_list = []
    with open(species_id_file) as tf:
        for line in tf:
            line = line.strip("\n").split(":")
            id_list.append(line[0])
    return set(id_list)


def translate_cluster_hash(header_list):
    """ hashes cluster file and returns hash-list in style of: [species_id] = list(protein_id)
    needs idpair_to_name pair as translation dict"""
    cluster_dict = defaultdict(list)
    for header in header_list:
        header = header.strip("\n")
        protein_sp_name = header.strip(">")
        try:
            protein_sp_id = idPair_2_namePair[protein_sp_name]
            sp_id, protein_id = protein_sp_id.split("_")
            cluster_dict[sp_id].append(protein_id)
        except KeyError:
            return None
    return cluster_dict


def next_best_blast_hit(unfiltered_protein_ids, self_hits):
    """unfiltered_protein_ids -> list of all protein ids of one species belonging to one cluster - CAN NOT BE TN
    self_hits -> [protein]= list of corresponding blast hits within one species, sorted by evalue."""
    for hit in self_hits:
        if hit not in unfiltered_protein_ids:
            return hit
    return None


def get_twin_hash(translated_filtered_cluster_header, translated_unfiltered_cluster_header):
    """generates a True negative fasta hash of a cluster_file by searching for each proteins next best blast hit.
    Takes as input files in style of [species_id]=list(protein_ids). One for the filtered data one for unfiltered"""
    fasta_hash = {}
    for species_id, prot_list in translated_filtered_cluster_header.items():
        for prot_id in prot_list:
            # two proteins of the same species are allowed to have the same next best
            protID_selfHits = blast_specID_protID_hitList[species_id][prot_id]
            next_best = next_best_blast_hit(translated_unfiltered_cluster_header[species_id], protID_selfHits)
            if next_best is not None:
                idx_pair = "{}_{}".format(species_id, next_best)
                name_pair = ">" + namePair_2_idPair[idx_pair]
                try:
                    sequence = all_protein_fasta_dict[name_pair][0]
                except KeyError:
                    logger_TN_Warning.info("could not translate {}".format(name_pair))
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
        self.length_dict = {}
        self.score_distribution_parameters = None
        self.length_distribution_parameters = None

    def query_for_fasta(self, query):
        """returns a fasta string for one single given header"""
        query = query.split()[0]
        return query + "\n" + "".join(self.fasta_hash[query])

    def seq_length_hash(self):
        for header in self.fasta_hash:
            self.length_dict[header] = len(self.fasta_hash[header][0])

    def calculate_length_distribution_parameters(self):
        """returns tuple (mean-std, mean+std) of length"""
        std_len = std(list(self.length_dict.values()), ddof=1)
        average_len = mean(list(self.length_dict.values()))
        self.length_distribution_parameters = (average_len - std_len, average_len + std_len)
        return self.length_distribution_parameters

    def calculate_score_distribution_parameters(self, true_negative_scores=None):
        """calculates score parameters either for TP distribution or TN and TP. If both, cutoff will be intersection point
        between two distributions."""
        if not self.score_dict:
            raise AttributeError("Empty Score List, exciting ... ")
        score_list = list(self.score_dict.values())
        tp_std = std(score_list, ddof=1)
        tp_average = mean(score_list)
        if true_negative_scores and len(true_negative_scores) >= 4:
                tn_std = std(true_negative_scores, ddof=1)
                tn_average = mean(true_negative_scores)
                inter_point = find_intersection(tp_average, tn_average, tp_std, tn_std)
                if inter_point:
                    self.score_distribution_parameters = inter_point
                else:
                    logger_Filtered.warning("No Intersection Point")
                    self.score_distribution_parameters = min(score_list)
        else:
            self.score_distribution_parameters = min(score_list)
        return self.score_distribution_parameters

    def generate_msa_string(self, rest_prot):
        """returns a MSA as string in fasta format, derived from the a list of headers"""
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
        if not self.length_dict:
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

    def compute_full_phmm(self, location):
        """realignes the fasta file stored in this object but to a different location"""
        msa_string = self.generate_msa_string(self.fasta_hash.keys())
        with tmp.NamedTemporaryFile() as msa_tmp:
            write_to_tempfile(msa_tmp.name, msa_string)
            generate_hmm(location, msa_tmp.name)
        return location

    def bulk_score_computation(self, length_normalized=False):
        if not self.length_dict:
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
        self.group_by_file_to_length_hash = defaultdict(dict)
        self.group_by_file_to_msa_obj = defaultdict(dict)
        self.group_by_file_to_hmm = defaultdict(dict)
        self.group_by_file_to_score_obj = defaultdict(dict)
        self.input_scope = 0
        self.valid_input_scope = 0

        # True negative cluster attributes
        self.group_by_file_to_unfiltered_header_set = defaultdict(dict)
        self.group_by_file_to_twin_hash = defaultdict(dict)
        self.group_by_file_to_twin_score_obj = defaultdict(dict)
        self.group_by_file_to_twin_hmm = defaultdict(dict)

    ####################################################################################################################
    # Overseer - Functions to specify input/output directories and hash cluster files
    ####################################################################################################################
    def feed_in_to_overseer(self, group_name, file_name, file_path):
        fasta_hash = check_and_hash_fasta(file_path, file_name)
        if fasta_hash is not None:
            self.group_by_file_to_cluster_hash[group_name][file_name] = fasta_hash
            self.group_by_file_to_unfiltered_header_set[group_name][file_name] = set(fasta_hash.keys())
            self.group_by_file_to_filepath[group_name][file_name] = file_path
            self.group_to_file_list[group_name].append(file_name)
            self.valid_input_scope += 1
        else:
            logger_Filtered.warning("Not enough entries in {}\n".format(file_name))
        return 1

    def hash_single_file_input(self):
        group_name = os.path.split(self.input_dir)[-2]
        file_name = ".".join(self.input_dir.split("/")[-1].strip().split(".")[0:-1])
        single_file = os.path.abspath(self.input_dir)
        self.group_to_result_path[group_name] = os.path.join(output_dir, group_name + ".GenePS")
        self.input_scope = self.feed_in_to_overseer(group_name, file_name, single_file)
        return self.valid_input_scope

    def hash_directory_input(self):
        for subdir, dirs, files in os.walk(self.input_dir):
            group_name = os.path.split(subdir.strip("/"))[-1]
            self.group_to_result_path[group_name] = os.path.join(output_dir, group_name + ".GenePS")
            for single_file in files:
                file_path = os.path.join(os.path.abspath(subdir), single_file)
                file_name = ".".join(single_file.strip().split(".")[0:-1])
                self.input_scope += self.feed_in_to_overseer(group_name, file_name, file_path)
        return self.valid_input_scope

    def initialize_input_data(self):
        if not os.path.exists(self.input_dir):
            print("\t[!] FATAL ERROR: '{}' is neither File nor Directory\n".format(self.input_dir))
            sys.exit()
        elif os.path.isfile(self.input_dir):
            valid_clusters = self.hash_single_file_input()
        else:
            valid_clusters = self.hash_directory_input()
        print("\n\t{} group(s) - containing {} file(s)\n\t{}\n".format(str(len(self.group_by_file_to_filepath)), str(self.input_scope), "-"*33))
        if valid_clusters == 0:
            print("\t[!] FATAL ERROR: No valid input files given\n")
            sys.exit()
        elif self.input_scope - self.valid_input_scope > 0:
            logging.info("[!] {} cluster(s) too small".format(self.input_scope - self.valid_input_scope))
            return self.valid_input_scope

    ####################################################################################################################
    # Overseer - MSA/HMM generation and normal true negative score computation
    ####################################################################################################################

    def remove_filtered_files(self, removal_group_to_file_list):
        for group, file_list in removal_group_to_file_list.items():
            for file_name in file_list:
                self.group_to_file_list[group].remove(file_name)
                self.valid_input_scope -= 1
        return self.valid_input_scope

    def output_HMM_and_fasta(self, directory):
        count = 1
        removed_group_to_file_list = defaultdict(list)
        for group, file_list in self.group_to_file_list.items():
            for file_name in file_list:
                msa_list = generate_msa(self.group_by_file_to_filepath[group][file_name])
                msa_obj = MsaObject(msa_list, file_name, directory)
                msa_obj.msa_to_fasta()
                msa_obj.trim_remove()
                if msa_obj.check_msa_size_and_length() is True:
                    if msa_obj.size_history[0] != msa_obj.size_history[-1]:
                        self.group_by_file_to_cluster_hash[group][file_name] = clean_fasta_hash(self.group_by_file_to_cluster_hash[group][file_name], msa_obj.all_header())
                        same_msa_path = write_hash_to_fasta(msa_obj.file_path, self.group_by_file_to_cluster_hash[group][file_name], ">{}\n{}\n")
                        msa_obj.re_align(same_msa_path)
                    self.group_by_file_to_msa_obj[group][file_name] = msa_obj
                    self.group_by_file_to_length_hash[group][file_name] = write_length_binned_fasta(self.group_by_file_to_cluster_hash[group][file_name], file_name, os.path.join(output_dir, file_name + ".fasta"))
                    self.group_by_file_to_hmm[group][file_name] = generate_hmm(os.path.join(output_dir, file_name + ".hmm"), msa_obj.file_path)
                else:
                    removed_group_to_file_list[group].append(file_name)
                    logger_Filtered.warning("Filtered due to MSA benchmarking {}".format(file_name))
                print_progress(count, self.valid_input_scope, prefix='\tGenerating Hidden Markov Models:\t', suffix='Complete', bar_length=30)
                count += 1
        return self.remove_filtered_files(removed_group_to_file_list)

    def compute_all_hmm_scores(self, length_normalized=False):
        count = 1
        print("\n")
        for group, file_list in self.group_to_file_list.items():
            for file_name in file_list:
                fasta_hash = self.group_by_file_to_cluster_hash[group][file_name]
                self.group_by_file_to_score_obj[group][file_name] = ScoreObject(fasta_hash, self.group_by_file_to_hmm[group][file_name])
                self.group_by_file_to_score_obj[group][file_name].length_dict = self.group_by_file_to_length_hash[group][file_name]
                if len(fasta_hash) < 20:
                    score_hash = self.group_by_file_to_score_obj[group][file_name].iterative_score_computation(length_normalized=length_normalized)
                else:
                    score_hash = self.group_by_file_to_score_obj[group][file_name].bulk_score_computation(length_normalized=length_normalized)
                print_progress(count, self.valid_input_scope, prefix='\tComputing HMM Score Distributions:\t', suffix='Complete', bar_length=30)
                count += 1
                if keep:
                    keep_file = write_hash_to_fasta(os.path.join(keep_dir, "{}_{}_scores.txt".format(group, file_name)), score_hash, "{}\t{}\n")
        return self.group_by_file_to_score_obj

    ####################################################################################################################
    # Overseer - Functions for generating True negative scores
    ####################################################################################################################

    def make_cluster_specific_TN_hash(self, group, file_name):
        """searches for next best blast hit for each protein in the filtered fasta file. The next best hits must not
        be part of the unfiltered fasta file to avoid having a TP in the TN set!"""
        translated_unfiltered_cluster_header = translate_cluster_hash(self.group_by_file_to_unfiltered_header_set[group][file_name])
        if translated_unfiltered_cluster_header is None:
            logger_TN_Warning.warning("Omitting True negative computation for {}".format(file_name))
            return None
        translated_filtered_cluster_header = translate_cluster_hash(self.group_by_file_to_cluster_hash[group][file_name].keys())
        fasta_hash = get_twin_hash(translated_filtered_cluster_header, translated_unfiltered_cluster_header)
        self.group_by_file_to_twin_hash[group][file_name] = fasta_hash
        return fasta_hash

    def compute_true_negative_hmm_scores(self, length_normalized=False):
        count = 1
        print("\n")
        for group, file_list in self.group_to_file_list.items():
            for file_name in file_list:
                hmm = self.group_by_file_to_hmm[group][file_name]
                fasta_hash = self.make_cluster_specific_TN_hash(group, file_name)
                if fasta_hash is not None:
                    scoring_obj = ScoreObject(fasta_hash, hmm)
                    score_hash = scoring_obj.bulk_score_computation(length_normalized=length_normalized)
                    self.group_by_file_to_twin_score_obj[group][file_name] = scoring_obj
                    self.group_by_file_to_twin_hmm[group][file_name] = scoring_obj.compute_full_phmm(os.path.join(output_dir, file_name + ".TN_hmm"))
                    print_progress(count, self.valid_input_scope, prefix='\tTrue Negative Score Distributions:\t', suffix='Complete', bar_length=30)
                    if keep:
                        keep_file = write_hash_to_fasta(os.path.join(keep_dir, "{}_{}_TrueNegativeScores.txt".format(group, file_name)), score_hash, "{}\t{}\n")
                else:
                    pass
                count += 1
        return self.group_by_file_to_twin_score_obj


##################
# global variables
##################

console = logging.StreamHandler()
console.setLevel(logging.INFO)
logger_Filtered = logging.getLogger("Filtered")
logger_TN_Warning = logging.getLogger("TranslationFiles")

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
        number_blast_files = hash_all_blast_files(species_ids)
        if not number_blast_files - len(species_ids) == 0:
            print("\t[!] FATAL ERROR: Not all Blast files could be hashed\n")
            sys.exit()
    if keep:
        keep_dir = get_outdir(output_dir, add_dir="intermediate_files")
    log_path = os.path.join(output_dir, "LOG.txt")
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', filename=log_path, filemode='w')

    # Start
    overseer_obj = Overseer(infile)
    filtered_data_scope = overseer_obj.initialize_input_data()
    logging.info("# {} groups and {} files\n".format(len(overseer_obj.group_to_file_list), str(overseer_obj.input_scope)))
    with tempdir() as temp_dir:
        filtered_data_scope = overseer_obj.output_HMM_and_fasta(temp_dir)
        if not filtered_data_scope > 0:
            print("\t[!] FATAL ERROR: NO Multiple Sequence Alignments computable\n")
            sys.exit()
        # score distributions for all cluster
        all_score_hashes = overseer_obj.compute_all_hmm_scores(length_normalized=True)
        if true_negative_file:
            all_tn_scores_hash = overseer_obj.compute_true_negative_hmm_scores(length_normalized=True)
        print("\n")
        read_count = 1
        for name_group, all_files in overseer_obj.group_to_file_list.items():
            with open(overseer_obj.group_to_result_path[name_group], "w") as results_file:
                results_file.write("group: {}\ngroup_size: {}\n".format(name_group, str(len(all_files))))
                for cluster_name in all_files:
                    consensus = get_consensus(overseer_obj.group_by_file_to_hmm[name_group][cluster_name])
                    length_range = overseer_obj.group_by_file_to_score_obj[name_group][cluster_name].calculate_length_distribution_parameters()
                    if true_negative_file:
                        try:
                            tn_scores = list(overseer_obj.group_by_file_to_twin_score_obj[name_group][cluster_name].score_dict.values())
                        except AttributeError:
                            tn_scores = [0]
                        score_cut_off = overseer_obj.group_by_file_to_score_obj[name_group][cluster_name].calculate_score_distribution_parameters(true_negative_scores=tn_scores)
                    else:
                        score_cut_off = overseer_obj.group_by_file_to_score_obj[name_group][cluster_name].calculate_score_distribution_parameters()
                    results_file.write("#name: {}\n#score_cut_off: {}\n#length_range: {},{}\n{}\n".format(cluster_name, score_cut_off, length_range[0], length_range[1], consensus))
                    print_progress(read_count, overseer_obj.valid_input_scope, prefix='\tWriting Results to Files:\t\t', suffix='Complete', bar_length=30)
                    read_count += 1
                results_file.close()
    logging.info("# {} filtered files\n".format(overseer_obj.valid_input_scope))
    print("\nDONE!\n")
