#!/usr/bin/env python3
###############
# MAKE GenePS !
###############

"""
Usage: make_GenePS.py                         -i <DIR> -o <DIR>

    Options:
        -h, --help                            show this screen.

        General
        -i, --input <DIR>                     either single input file or directory with files or sub-folders
        -o, --output <DIR>                    directory where to save the output file(s)

"""

import os
import sys
import tempfile as tmp
from collections import defaultdict
from docopt import docopt

from run_command import run_cmd, tempdir, check_programs
from compute_msa import generate_msa, MsaObject


##################
# Global functions
##################


class MsaSizeError(Exception):
    pass


class MsaLengthError(Exception):
    pass


class InputFileError(Exception):
    pass


def walk_through_input(input_dir):
    dir_hash = defaultdict(list)
    if not os.path.isdir(input_dir):
        dir_hash[os.path.dirname(input_dir)].append(os.path.split(input_dir)[-1])
        return dir_hash
    else:
        for subdir, dirs, files in os.walk(input_dir):
            for single_file in files:
                single_file = single_file
                single_folder = subdir
                dir_hash[single_folder].append(single_file)
        return dir_hash


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
    if len(fasta_dict) <= 4:
        raise InputFileError
    else:
        pass


def check_size_len_after_trimal(trimmed_msa):
    if trimmed_msa.size[-1] < 4 \
            or ((trimmed_msa.size[0] - trimmed_msa.size[-1]) / trimmed_msa.size[0]) > 50:
        raise MsaSizeError

    elif trimmed_msa.lengths[-1] < 20:
        raise MsaLengthError
    else:
        pass


def get_outdir(output_dir, add_dir=""):
    if not os.path.isdir(output_dir):
        raise SystemExit
    elif not os.path.exists(os.path.join(output_dir, add_dir)):
        os.mkdir(os.path.join(output_dir, add_dir))
        return os.path.join(output_dir, add_dir)
    else:
        return os.path.join(output_dir, add_dir)


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
    for line in run_cmd(command=command, wait=False):
        if "E-value" in line or read_count == 1:
            read_count += 1
        elif read_count == 2:
            line = line.strip("\n").split()
            return float(line[1])


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


###############
# Scoring class
###############

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

    def compute_scores(self):
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
                        score = get_phmm_score(hmm_tmp.name, q_tmp.name)
            self.score_list.append(int(score))
        return self.score_list

    def compute_full_phmm(self):
        msa_string = self.generate_msa_string(self.left_taxa)
        with tmp.NamedTemporaryFile() as msa_tmp:
            write_to_tempfile(msa_tmp.name, msa_string)
            self.hmm_path = os.path.join(self.dir, self.name)
            generate_hmm(self.hmm_path, msa_tmp.name)
        return self.hmm_path


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    infile = args['--input']
    out_directory = args['--output']

    check_programs("hmmsearch", "hmmemit", "hmmbuild", "mafft", "trimal")

    try:
        out_dir = get_outdir(out_directory)
    except SystemExit:
        print(out_directory, " is NOT a directory!")
        print("Please specify an output directory")
        sys.exit()

    # tmp_dir automatically cleaned up
    with tempdir() as tmp_dir:
        counter = 1
        dir_tree = walk_through_input(infile)
        number_groups = len(dir_tree)
        print("#" * 27)
        print("# {} - groups given as input".format(str(number_groups)))
        print("#" * 27 + "\n")

        # folder = group; file = cluster
        for folder, file_list in dir_tree.items():
            number_files = len(file_list)
            folder_name = folder.split("/")[-1]
            print("[{}] starting with {} -> containing {} files\n"
                  .format(str(counter), folder_name, str(number_files)))

            # single results file per group
            results_path = os.path.join(out_dir, folder_name + ".makeGenePS")
            results_file = open(results_path, "w")
            results_file.write("#group: {}\n#group_size: {}\n".format(folder_name, str(number_files)))

            # run script for all files
            for file in file_list:
                file_name = file.strip().split(".")[0]
                file_path = os.path.join(folder, file)
                print("\tanalyzing:\t {}\n".format(file_name))
                try:
                    fasta_hash = hash_fasta(file_path)
                    check_for_sufficient_taxa(fasta_hash)
                except InputFileError:
                    print("\t[!] {} does not contain enough entries\n".format(file_name))
                    continue

                # first MSA for pHMM consensus
                msa_list = generate_msa(file_path)
                msa_file = MsaObject(msa_list, file_name, tmp_dir)
                msa_file.msa_to_fasta()
                msa_file.trim_remove()
                msa_file.trim_length()
                try:
                    check_size_len_after_trimal(msa_file)
                except MsaSizeError:
                    print("\t[!] {} : NO MSA computable - "
                          "only {} taxa remained MSA after trimming\n".format(msa_file.name, msa_file.size))
                    continue
                except MsaLengthError:
                    print("\t[!] {} : NO MSA computable - "
                          "MSA length too short after trimming\n".format(msa_file.name))
                    continue

                # pHMM consensus
                cons_hmm = os.path.join(tmp_dir, file_name + ".chmm")
                generate_hmm(cons_hmm, msa_file.path)
                consensus_seq = get_consensus(cons_hmm)

                # compute scores
                left_taxa = msa_file.all_header()
                scoring_obj = ScoreObject(fasta_hash, left_taxa, file_name, out_dir)
                score_list = scoring_obj.compute_scores()
                phmm_path = scoring_obj.compute_full_phmm()
                print("\tscores: ", score_list, "\n")

                # write to results
                results_file.write(">name: {}\n>phmm_dir: {}\n>score_list: {}\n{}\n".
                                   format(file_name, phmm_path, score_list, consensus_seq))
            counter += 1
            results_file.close()



# to add:
    # verbose and "keep" options (exchanging tmpdir)
