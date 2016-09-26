#!/usr/bin/env python3
from run_command import run_cmd
import os


def msa_operations(command):
    read_flag = 0
    seq, msa_list = [], []
    for line in run_cmd(command=command, wait=False):
        if read_flag == 0:
            if line.startswith(">"):
                read_flag = 1
        if read_flag == 1:
            if line.startswith(">"):
                if len(msa_list) > 0:
                    msa_list.append("".join(seq))
                msa_list.append(line.rstrip("\n"))
                seq = []
            else:
                seq.append(line.rstrip("\n"))
    msa_list.append("".join(seq))
    return msa_list


def generate_msa(fasta):
    command = "mafft " + fasta
    return msa_operations(command)
    # vll einführen das es automatisch ein boject wird


class MsaObject:
    def __init__(self, msa_list, cluster_name, output_dir):
        self.msa_list = msa_list
        self.cluster_name = cluster_name
        self.name = cluster_name + ".msaGenePS"
        self.path = os.path.join(output_dir, self.name)

        self.size = [len(self.msa_list)/2]
        self.lengths = [len(self.msa_list[1])]
        self.cmd_trim_remove = "trimal -in " + self.path + " -resoverlap 0.6 -seqoverlap 80"
        self.cmd_trim_length = "trimal -in " + self.path + " -automated1"

    def msa_to_fasta(self):
        with open(self.path, "w") as m:
            m.write("\n".join(self.msa_list) + "\n")

    def trim_remove(self):
        self.msa_list = msa_operations(self.cmd_trim_remove)
        self.size.append(len(self.msa_list)/2)
        if self.size[-1] < 4 or \
                        ((self.size[0] - self.size[-1]) / self.size[0]) > 50:
            return print(self.name, " : NO MSA computable")
        self.msa_to_fasta()

    def trim_length(self):
        self.msa_list = msa_operations(self.cmd_trim_length)
        self.lengths.append(len(self.msa_list[1]))
        if self.lengths[-1] < 20:
            # print statement to a log file and on screen
            # return a key word which jumps to the next cluster
            return print(self.name, " : NO MSA computable")
        self.msa_to_fasta()

    def all_header(self):
        return [self.msa_list[x] for x in range(0, len(self.msa_list), 2)]

    def all_aln(self):
        return [self.msa_list[x] for x in range(1, len(self.msa_list), 2)]



