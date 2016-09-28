#!/usr/bin/env python3
import re
import os
from collections import defaultdict
from run_command import run_cmd


aa3_to_1_coding_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 'Ter':'*',
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M','Xaa':'X'}


def aacode_3to1(seq):
    '''Turn a three letter protein into a one letter protein.
The 3 letter code can be upper, lower, or any mix of cases
The seq input length should be a factor of 3 or else results
in an error
'''

    if len(seq) % 3 == 0:
        single_seq = []
        for i in range(0, len(seq), 3):
            single_seq.append(aa3_to_1_coding_dict.get(seq[i:i+3]))
        return "".join(single_seq)
    else:
        print("ERROR: Sequence was not a factor of 3 in length!")


def remove_non_letter_signs(string):
    '''removes all non-alphabetic letters from string'''
    regex = re.compile('[^a-zA-Z]')
    return regex.sub("", string)


def grap_values(attribute):
    out_list = []
    for query, trange in attribute.items():
        for i, last_value in trange.items():
            out_list.append(last_value)
    return out_list


class ExonerateObject:

    def __init__(self, exonerate_file):
        self.path = os.path.realpath(exonerate_file)
        self.header = defaultdict(dict)
        self.query_prot = defaultdict(lambda: defaultdict(str)) #three leafs dict
        self.target_dna = defaultdict(lambda: defaultdict(str))
        self.target_prot = defaultdict(lambda: defaultdict(str))
        self.gff = defaultdict(lambda: defaultdict(list))
        self.exonerate_processor(exonerate_file)

    def exonerate_processor(self, exonerate_file):
        ''' protein against 1 or more genome regions (targets),
        file with aln and/or gff output '''
        del_intron = lambda seq: seq.replace("TargetIntron","")                  # no lambda
        prep_line = lambda line: line.strip("\n").split(": ")[1].strip(" ")
        next_line = lambda: prep_line(next(ex))
        next_block = lambda size: [next_line() for x in range(0,size)]

        read_flag = 0
        with open(exonerate_file, "r") as ex:
            for line in ex:
                if line.strip():
                    if read_flag == 0:
                        if line.startswith("C4 Alignment"):
                            read_flag = 1
                    elif read_flag == 1:
                        query, target, model, score= next_block(4)
                        qrange, trange = [tuple(x.split(" -> ")) for x in next_block(2)]
                        self.header[target][trange] = {"model": model, "score": score,
                                                       "qrange": qrange, "query": query}
                        read_flag = 2
                    elif read_flag == 2:
                        if not line.startswith("#"):
                            self.query_prot[target][trange] += (del_intron(remove_non_letter_signs(line)))
                            next(ex)
                            self.target_prot[target][trange] += (remove_non_letter_signs(next(ex)))
                            self.target_dna[target][trange] += (remove_non_letter_signs(next(ex)))
                        elif "GFF" in line:
                            read_flag = 3
                        else:
                            read_flag = 0
                    elif read_flag == 3:
                        if not line.startswith("#"):
                            line = line.strip("\n").split("\t")[2:7]
                            line.pop(3)
                            self.gff[target][trange].append(line)
                        elif "END OF GFF DUMP" in line:
                            read_flag = 0
                    else:
                        pass
        for target, trange in self.target_prot.items():
            for trange, seq in trange.items():
                self.target_prot[target][trange] = aacode_3to1(seq)


def run_exonerate(name, directory, region, query):
    cmd = "exonerate -m p2g:b --softmaskquery no -E yes -Q protein -T dna -n 1 " \
      "--softmasktarget yes --showvulgar no --minintron 20 --maxintron 50000 " \
      "--showalignment yes --showtargetgff yes -q {} -t {}".format(query, region)
    with open(os.path.join(directory, name), "w") as ex:
        for line in run_cmd(command=cmd, wait=False):
            ex.write(line)
    return ExonerateObject(os.path.join(directory, name))


if __name__ == "__main__":
    test = run_exonerate("test_exonerate.out", "/home/jgravemeyer/Dropbox/MSc_project/data",
                  "/home/jgravemeyer/Dropbox/MSc_project/src/GenePS/test_data/blast_region.fasta",
                  "/home/jgravemeyer/Dropbox/MSc_project/src/GenePS/test_data/test_consensus_eef.fa")

    print(grap_values(test.header)[0]["query"])
    print(grap_values(test.target_prot))
    print(grap_values(test.target_dna))
    print(grap_values(test.query_prot))
    print(grap_values(test.gff))
