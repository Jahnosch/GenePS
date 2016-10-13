#!/usr/bin/env python3
###############
# MAKE GenePS !
###############

"""
Usage: make_GenePS.py                         -i <DIR> -o <DIR> [-b <DIR>] [-t <FILE>] [-p <FILE>]

    Options:
        -h, --help                            show this screen.

        General
        -i, --input <DIR>                     either single input file or directory with files or sub-folders
        -o, --output <DIR>                    directory where to save the output file(s)
        -b, --blast <DIR>                     directory containing within-species-BLAST text files (one per species)
        -t, --translation_file <FILE>         file to translate species/protein names to ids (69_25600: MJAVA.g25693.t1)
        -p, --protein_fasta <FILE>            fasta-file providing all protein sequences with name as header

"""

import os
import sys
import tempfile as tmp
import operator
from collections import defaultdict

from run_command import run_cmd, tempdir, check_programs
from compute_msa import generate_msa, MsaObject

import_errors = []
try:
    from docopt import docopt
except ImportError:
    import_errors.append("[ERROR] : Module \'Docopt\' was not found. Please install \'Docopt\' using \'pip install docopt\'")
try:
    import matplotlib.pyplot as plt
except ImportError:
    import_errors.append("[ERROR] : Module \'matplotlib.pyplot\' was not found. Please install \'matplotlib\' using your package manager")
try:
    import seaborn as sns
    sns.set(color_codes=True)
except ImportError:
    import_errors.append("[ERROR] : Module \'seaborn\' was not found. Please install \'seaborn\' using \'pip install seaborn\'")
if import_errors:
    sys.exit("\n".join(import_errors))

########################################################################################################################
# General functions
########################################################################################################################


class MsaSizeError(Exception):
    pass


class MsaLengthError(Exception):
    pass


class InputFileError(Exception):
    pass


def walk_through_dir(input_dir):
    dir_hash = defaultdict(list)
    if os.path.isfile(input_dir):
        dir_hash[os.path.abspath(os.path.dirname(input_dir))].append(os.path.split(input_dir)[-1])
        return dir_hash
    else:
        for subdir, dirs, files in os.walk(input_dir):
            for single_file in files:
                single_file = single_file
                single_folder = os.path.abspath(subdir)
                dir_hash[single_folder].append(single_file)
        return dir_hash


def check_optional_arguments():
    if args['--blast'] and args['--translation_file'] and args['--protein_fasta']:
        return True
    else:
        return False


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


def check_for_sufficient_taxa(fasta_dict):
    if len(fasta_dict) <= 2:
        raise InputFileError
    else:
        pass


def check_size_len_after_trimal(trimmed_msa):
    if trimmed_msa.size[-1] < 2 \
            or ((trimmed_msa.size[0] - trimmed_msa.size[-1]) / trimmed_msa.size[0]) > 50:
        raise MsaSizeError

    elif trimmed_msa.lengths[-1] < 20:
        raise MsaLengthError
    else:
        pass


def get_outdir(output_dir, add_dir=""):
    if os.path.isfile(output_dir):
        raise SystemExit
    elif not os.path.exists(os.path.join(output_dir, add_dir)):
        os.mkdir(os.path.join(output_dir, add_dir))
        return os.path.abspath(os.path.join(output_dir, add_dir))    # absolute paths
    else:
        return os.path.abspath(os.path.join(output_dir, add_dir))


def generate_hmm (hmm_path, msa_path):
    command = ["hmmbuild", hmm_path, msa_path]
    run_cmd(command=command, wait=True)
    return hmm_path


def get_phmm_score(hmm_file, query_file, test_this_function=None):
    if test_this_function is None:
        command = ["hmmsearch", "--noali", hmm_file, query_file]
    else:
        command = test_this_function
    read_count = 0
    score_list = []
    for line in run_cmd(command=command, wait=False):
        if "E-value" in line or read_count == 1:
            read_count += 1
        elif read_count == 2:
            line = line.strip("\n").split()     # line[8] is protein name
            if len(line) > 0:
                try:
                    score_list.append(round(float(line[1])))
                except ValueError:
                    print("VALUE ERROR")
                    print(line)
                    return score_list
            else:
                return score_list


def get_consensus(hmm_file):
    command = "hmmemit -c " + hmm_file
    cons_list = []
    for line in run_cmd(command=command, wait=False):
        if not line.startswith(">"):
            cons_list.append(line.strip("\n"))
    return "".join(cons_list)


# future use
def progress(iteration, steps, max_value):
    if int(iteration) == int(max_value):
        sys.stdout.write('\r')
        print ("[PROGRESS]\t: %d%%" % (100))
    elif int(iteration) % int(steps+1) == 0:
        sys.stdout.write('\r')
        print ("[PROGRESS]\t: %d%%" % (float(int(iteration)/int(max_value))*100))
        sys.stdout.flush()
    else:
        pass


def write_to_tempfile (tmp_name, string):
    new_file = open(tmp_name, "w")
    new_file.write(string)
    new_file.seek(0)
    new_file.close()


########################################################################################################################
# Functions: True Negative Score computation
########################################################################################################################

def get_blast_files(input_dir):
    """returns the directory of the folder storing the blast files and a set of all blast file names"""
    dir_list = []
    if not os.path.isdir(input_dir):
        print (input_dir, " : has to be directory")
        sys.exit()
    else:
        for subdir, dirs, files in os.walk(input_dir):
            for single_file in files:
                if not single_file.startswith("Blast"):
                    print("not only blast files found in: ", input_dir)
                    sys.exit()
                dir_list.append(single_file)
            return subdir, set(dir_list)


def hash_translator(name_id_file):
    """ hashes translation file and returns hash-list in style of: [species_protein_name] = species_protein_id
     as well as an reversed hash: [species_protein_id] = species_protein_name """
    translation_dict = {}
    translation_dict_rev = {}
    with open(name_id_file) as tf:
        for line in tf:
            line = line.strip("\n").split(":")
            prot_sp_name = line[1].strip()
            prot_sp_id = line[0]
            translation_dict[prot_sp_name] = prot_sp_id
            translation_dict_rev[prot_sp_id] = prot_sp_name
    return translation_dict, translation_dict_rev


def make_cluster_name_to_id_hash(in_file, translation_dict):
    """ hashes cluster file and returns hash-list in style of: [species_id] = list(protein_id) """
    cluster_dict = defaultdict(list)
    with open(in_file) as cf:
        for line in cf:
            if line.startswith(">"):
                line = line.strip("\n")
                prot_sp_name = line.strip(">")
                prot_sp_id = translation_dict[prot_sp_name]
                sp_id, prot_id = prot_sp_id.split("_")
                cluster_dict[sp_id].append(prot_id)
    return cluster_dict


def blast_prot_id_hits(blast_file):
    """ hashes within species blast and returns hash in style of: [protein_id] = list(hits) """
    protid_hits_hash = defaultdict(list)
    with open(blast_file) as blast_f:
        for line in blast_f:
            line = line.split("\t")
            spec_id, prot_id = line[0].split("_")   # line[0] = self, line[1] = other
            protid_hits_hash[prot_id].append(line[1].split("_")[1])
    return protid_hits_hash


def twin_hmm_score(header, seq, hmm):
    """Hmm score against the "real" cluster pHMM"""
    fasta_str = header + "\n" + seq
    with tmp.NamedTemporaryFile() as temp:
        write_to_tempfile(temp.name, fasta_str)
        try:
            tn_score = get_phmm_score(hmm, temp.name)
        except IndexError:
            return None
    return tn_score


def get_next_best(forbidden, self_hits):
    """forbidden = list of all protein ids of one species in the same cluster
    self_hits = list of all the hits from blasting a protein within a species"""
    for hit in self_hits:
        if hit not in forbidden:
            return hit
    return None


def get_twin_file_scores(cluster_file, cluster_hmm_path, forbidden_proteins_file):
    """give the next best blast hits for protein in the cluster-file. Needs a global dict(dict(list))"""

    # hash all files
    all_proteins = hash_fasta(all_cluster_all_protein)
    idPair_2_namePair, namePair_2_idPair = hash_translator(translation_file)
    infile_specID_protList = make_cluster_name_to_id_hash(cluster_file, idPair_2_namePair)
    forbidden_specID_protList = make_cluster_name_to_id_hash(forbidden_proteins_file, idPair_2_namePair)
    blast_path, blast_file_set = get_blast_files(blast_dir)

    # generate fasta_string
    results_fasta = []
    results_score = []
    for species_id, prot_list in infile_specID_protList.items():
        b_file_name = "Blast{}_{}.txt".format(species_id, species_id)
        if b_file_name in blast_file_set:
            if not blast_specID_protID_hitList[species_id]:
                b_file_path = os.path.join(blast_path, b_file_name)
                blast_specID_protID_hitList[species_id] = blast_prot_id_hits(b_file_path)
            for prot_id in prot_list:
                if prot_id in blast_specID_protID_hitList[species_id]:
                    # two proteins of the same species are allowed to have the same next best
                    protID_selfHits = blast_specID_protID_hitList[species_id][prot_id]
                    next_best = get_next_best(forbidden_specID_protList[species_id], protID_selfHits)
                    if next_best is not None:
                        idx_pair = "{}_{}".format(species_id, next_best)
                        name_pair = ">" + namePair_2_idPair[idx_pair]
                        try:
                            sequence = all_proteins[name_pair][0]
                        except KeyError:
                            continue
                        results_fasta.append(name_pair)
                        results_fasta.append(sequence)
                        tn_score = twin_hmm_score(name_pair, sequence, cluster_hmm_path)
                        if tn_score is not None:
                            results_score.append(tn_score)
    return "\n".join(results_fasta), results_score


def call_true_negative_twin(cluster_headers, cluster_hmm, forbidden_names, verbose = None):
    """calls the main get_twin_file_scores function but checks for errors. Forbidden names = set of proteins
    before trimming of spurious alignments"""
    TN_score_list = None
    seq_fasta = None
    if optional_arguments:
        header_fasta = "\n".join(cluster_headers)
        forbidden_fasta = "\n".join(forbidden_names)
        with tmp.NamedTemporaryFile() as header_f:
            write_to_tempfile(header_f.name, header_fasta)
            with tmp.NamedTemporaryFile() as forbidden_f:
                write_to_tempfile(forbidden_f.name, forbidden_fasta)
                seq_fasta, TN_score_list = get_twin_file_scores(header_f.name, cluster_hmm, forbidden_f.name)
    return TN_score_list, seq_fasta


########################################################################################################################
# Class ScoreObject - Positive Scores
########################################################################################################################

class ScoreObject:
    def __init__(self, fasta_dict, taxa, name, directory):
        self.fasta_hash = fasta_dict
        self.left_taxa = taxa
        self.dir = directory
        self.name = name + ".hmmGenePS"
        self.hmm_path = None
        self.score_list = []

    def query_for_fasta(self, query):
        query = query.split()[0]
        return query + "\n" + "".join(self.fasta_hash[query])

    def generate_msa_string(self, rest_prot):
        seq_list = []
        for taxa in rest_prot:
            taxa = taxa.split()[0]
            seq_list.append(taxa)
            seq_list.append("".join(self.fasta_hash[taxa]))
        with tmp.NamedTemporaryFile() as r_tmp:
            write_to_tempfile(r_tmp.name, "\n".join(seq_list))
            list_msa = generate_msa(r_tmp.name)
        return "\n".join(list_msa)

    def iterative_score_computation(self):
        score_dict = {}
        for idx in range(0, len(self.left_taxa)):
            rest_prot = self.left_taxa[:]
            query = rest_prot.pop(idx)
            with tmp.NamedTemporaryFile() as q_tmp:
                write_to_tempfile(q_tmp.name, self.query_for_fasta(query))
                msa_string = self.generate_msa_string(rest_prot)
                with tmp.NamedTemporaryFile() as msa_tmp:
                    write_to_tempfile(msa_tmp.name, msa_string)
                    with tmp.NamedTemporaryFile() as hmm_tmp:
                        generate_hmm(hmm_tmp.name, msa_tmp.name)
                        try:
                            score = get_phmm_score(hmm_tmp.name, q_tmp.name)[0]
                        except IndexError:
                            continue
        #    score_dict[query] = int(score)
            self.score_list.append(int(score))
       # sorted_dict = sorted(score_dict.items(), key=operator.itemgetter(1))
       # unique_species = set([])
       # for x in sorted_dict:
           # species = x[0].split(".")[0]
           # if species not in unique_species:
            #    unique_species.add(species)
                #print(x[0], ",", x[1])
            #else:
                #print(x[0], ",", x[1], ",\t X")
        return self.score_list

    def compute_full_phmm(self):
        msa_string = self.generate_msa_string(self.left_taxa)
        with tmp.NamedTemporaryFile() as msa_tmp:
            write_to_tempfile(msa_tmp.name, msa_string)
            self.hmm_path = os.path.join(self.dir, self.name)
            generate_hmm(self.hmm_path, msa_tmp.name)
        return self.hmm_path

    def bulk_score_computation(self):
        with tmp.NamedTemporaryFile() as q_tmp:
            seq_list = []
            for header in self.left_taxa:
                seq_list.append(header)
                seq_list.append(self.fasta_hash[header][0])
            write_to_tempfile(q_tmp.name, "\n".join(seq_list))
            self.compute_full_phmm()
            score_list = get_phmm_score(self.hmm_path, q_tmp.name)
        return score_list, self.hmm_path


########################################################################################################################
# Function: Invokes Functions for computation of single cluster
########################################################################################################################

def main_process(tmp_dir, out_dir, cluster_file, folder, verbose = None):
    """summary of all function needed to process one single file"""

    # preparing input file
    file_name = cluster_file.strip().split(".")[0]
    file_path = os.path.join(folder, cluster_file)
    print("\tanalyzing:\t {}\n".format(file_name))
    try:
        fasta_hash = hash_fasta(file_path)
        check_for_sufficient_taxa(fasta_hash)
    except InputFileError:
        print("\t[!] {} does not contain enough entries\n".format(file_name))
        return None

    # first MSA for pHMM consensus
    print("\tComputing Multiple Sequence Alignment\n")
    msa_list = generate_msa(file_path)
    msa_file = MsaObject(msa_list, file_name, tmp_dir)
    original_protein_set = msa_file.all_header()
    msa_file.msa_to_fasta()
    msa_file.trim_remove()
    msa_file.re_align_to_fasta(fasta_hash)
    try:
        check_size_len_after_trimal(msa_file)
    except MsaSizeError:
        print("\t[!] {} : NO MSA computable - "
              "only {} taxa remained MSA after trimming\n".format(msa_file.name, msa_file.size))
        return None
    except MsaLengthError:
        print("\t[!] {} : NO MSA computable - "
              "MSA length too short after trimming\n".format(msa_file.name))
        return None

    # pHMM consensus
    print("\tGenerating Hidden Markov Model\n")
    cons_hmm = os.path.join(tmp_dir, file_name + ".chmm")
    generate_hmm(cons_hmm, msa_file.path)
    consensus_seq = get_consensus(cons_hmm)

    # compute scores
    print("\tComputing HMM-Score Distribution\n")
    left_proteins = msa_file.all_header()
    scoring_obj = ScoreObject(fasta_hash, left_proteins, file_name, out_dir)
    if len(left_proteins) < 20:
        score_list = scoring_obj.iterative_score_computation()
        phmm_path = scoring_obj.compute_full_phmm()
    else:
        score_list, phmm_path = scoring_obj.bulk_score_computation()
    print("\tHMM Scores: ", score_list, "\n")

    # score True negative "twin file"
    TrueN_scores, TrueN_fasta = call_true_negative_twin(left_proteins, phmm_path, original_protein_set)

    return file_name, phmm_path, score_list, consensus_seq, TrueN_scores


########################################################################################################################
# main
########################################################################################################################

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    infile = args['--input']
    out_directory = args['--output']
    optional_arguments = check_optional_arguments()
    if optional_arguments:
        blast_dir = os.path.abspath(args['--blast'])
        translation_file = os.path.abspath(args['--translation_file'])
        all_cluster_all_protein = os.path.abspath(args['--protein_fasta'])
        blast_specID_protID_hitList = defaultdict(lambda: defaultdict(list))
        gather_next_best_scores = []
        gather_intra_cluster_scores = []

    check_programs("hmmsearch", "hmmemit", "hmmbuild", "mafft", "trimal")

    # other function, looks ugly
    dir_tree = walk_through_dir(infile)
    try:
        output_dir = get_outdir(out_directory)
    except SystemExit:
        print(out_directory, " is NOT a directory!")
        print("Please specify an output directory")
        sys.exit()
    with tempdir() as temp_dir:
        counter = 1
        number_groups = len(dir_tree)
        print("#" * 27)
        print("# {} - group(s) given as input".format(str(number_groups)))
        print("#" * 27 + "\n")
        # folder = group; file = cluster
        for group_folder, file_list in dir_tree.items():
            number_files = len(file_list)
            folder_name = group_folder.split("/")[-1]
            print("[{}] starting with {} -> containing {} files\n".format(str(counter), folder_name, str(number_files)))
            # single results file per group/folder
            results_path = os.path.join(output_dir, folder_name + ".makeGenePS")
            results_file = open(results_path, "w")
            results_file.write("#group: {}\n#group_size: {}\n".format(folder_name, str(number_files)))
            # run main function processing each file
            for each_file in file_list:
                file__name, hmm_path, scores, consensus, TN_scores = main_process(temp_dir, output_dir, each_file, group_folder)
                results_file.write(">name: {}\n>phmm_dir: {}\n>score_list: {}\n{}\n".format(file__name, hmm_path, scores, consensus))
                if optional_arguments:
                    sns.distplot(scores, hist=False, rug=True, color="r", label="Scores within cluster")
                    sns.distplot(TN_scores, hist=False, rug=True, color="b", label="Scores next best BLAST hits")
                    plt.title("HMM-Score Distributions", size=18, weight="bold")
                    plt.xlabel("Score", size=14)
                    plt.ylabel("Density", size=14)
                    plt.show()
            counter += 1
            results_file.close()

