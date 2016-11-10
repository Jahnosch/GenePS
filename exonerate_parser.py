#!/usr/bin/env python3
import re
import os
from collections import defaultdict
from run_command import run_cmd
# if is going to be a parser, deal with softmasking

aa3_to_1_coding_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', '***': 'X',
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M','Xaa': 'X', 'Unk': 'X'}


def aacode_3to1(seq):
    """Turn a three letter protein into a one letter protein.
The 3 letter code can be upper, lower, or any mix of cases
The seq input length should be a factor of 3 or else results
in an error
"""
    if len(seq) % 3 == 0:
        single_seq = []
        for i in range(0, len(seq), 3):
            single_seq.append(aa3_to_1_coding_dict.get(seq[i:i+3]))
        return "".join(single_seq)
    else:
        raise False



def kill_char(string, n):
    if n > len(string):
        raise IndexError
    begin = string[:n]
    end = string[n+1:]
    return begin + end


def find_hashes(protein_string):
    for idx, sign in enumerate(protein_string):
        if sign == "#":
            yield idx


def clear_hashed_bases(protein_string, dna_string):
    """removes nucleotides which are opposite to #'s in the protein sequence (non-coding-nucleotides)"""
    if "#" in protein_string:
        for idx in sorted(find_hashes(protein_string), reverse=True):
            dna_string = kill_char(dna_string, idx)
    return dna_string


def remove_non_letter_signs(string):
    """removes all non-alphabetic letters from string"""
    regex = re.compile('[^a-zA-Z*]')
    return regex.sub("", string)


def remove_lower(text_string):
    return re.sub('[a-z]', '', text_string)


def grap_values(attribute):
    out_list = []
    for query, trange in attribute.items():
        for i, last_value in trange.items():
            out_list.append(last_value)
    return out_list


def write_exonerate_gff(gff_list):
    last_phase = 0
    new_gff_list = []
    mrna_copy = None
    for line in gff_list:
        line = line.split("\t")
        if "cds" in line[2]:
            cds_start = int(line[3])
            cds_end = int(line[4])
            current_phase = ((cds_end - cds_start + 1) + last_phase) % 3
            line[7] = str(last_phase)
            last_phase = current_phase
        elif "gene" in line[2]:
            mrna_copy = line[:]
            mrna_copy[2] = "mRNA"
        elif "similarity" in line[2] or "splice" in line[2]:
            continue
        new_gff_list.append("\t".join(line))
    new_gff_list.insert(1, "\t".join(mrna_copy))
    return new_gff_list


class ExonerateObject:

    def __init__(self, exonerate_file):
        self.path = os.path.realpath(exonerate_file)
        self.header = defaultdict(dict)
        self.query_prot = defaultdict(lambda: defaultdict(list)) #three leafs dict
        self.target_dna = defaultdict(lambda: defaultdict(list))
        self.target_prot = defaultdict(lambda: defaultdict(list))
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
                            self.query_prot[target][trange].append(del_intron(remove_non_letter_signs(line)))
                            next(ex)
                            target_prot = next(ex)
                            target_dna = clear_hashed_bases(target_prot, next(ex))
                            self.target_prot[target][trange].append(remove_non_letter_signs(target_prot))
                            self.target_dna[target][trange].append(remove_lower(remove_non_letter_signs(target_dna)))
                        else:
                            self.query_prot[target][trange] = "".join(self.query_prot[target][trange])
                            self.target_dna[target][trange] = "".join(self.target_dna[target][trange])
                            self.target_prot[target][trange] = aacode_3to1("".join(self.target_prot[target][trange]))
                            if "GFF" in line:
                                read_flag = 3
                            else:
                                read_flag = 0
                    elif read_flag == 3:
                        if not line.startswith("#"):
                            line = line.strip("\n")
                            self.gff[target][trange].append(line)
                        elif "END OF GFF DUMP" in line:
                            self.gff[target][trange] = write_exonerate_gff(self.gff[target][trange])
                            read_flag = 0
                    else:
                        pass
            ex.seek(0)


def make_exonerate_command(model, query_file, region_file):
    cmd = "exonerate {} --softmaskquery no -Q protein -T dna -n 1 " \
      "--softmasktarget no --showvulgar no --minintron 20 --maxintron 50000 " \
      "--showalignment yes --showtargetgff yes -q {} -t {}".format(model, query_file, region_file)
    return cmd


def safe_exonerate_output(output_path, command):
    line_count = 0
    with open(output_path, "w") as ex:
        for line in run_cmd(command=command, wait=False):
            ex.write(line)
            line_count += 1
    if line_count < 20:
        return None
    else:
        return ExonerateObject(output_path)


def run_exonerate(name, directory, region, query):
    cmd = make_exonerate_command("-m p2g:b -E yes", query, region)
    out_file = os.path.join(directory, name)
    exonerate_obj = safe_exonerate_output(out_file, cmd)
    if safe_exonerate_output(out_file, cmd) is None:
        cmd = make_exonerate_command("-m p2g -E no", query, region)
        exonerate_obj = safe_exonerate_output(out_file, cmd)
    return exonerate_obj


if __name__ == "__main__":
#    test = run_exonerate("test_exonerate.out", "/home/jgravemeyer/Dropbox/MSc_project/data",
#                  "/home/jgravemeyer/Desktop/blast_region.fasta",
#                  "/home/jgravemeyer/Desktop/test_consensus_eef.fa")

    test = run_exonerate("test_exonerate.out", "/home/jgravemeyer/Dropbox/MSc_project/data",
              "/home/jgravemeyer/Desktop/eef_region_elegans.fa",
              "/home/jgravemeyer/Desktop/eef_consensus.fa")


    for target in test.query_prot:
        print(list(test.query_prot[target].keys())[0])
        for trange in test.query_prot[target]:
            print(">" + test.header[target][trange]["query"])
            print(test.target_prot[target][trange])
            print("\n".join(test.gff[target][trange]))
'''
    print(grap_values(test.header)[0].values())
    print(grap_values(test.target_prot))
    print(grap_values(test.target_dna))
    print(grap_values(test.query_prot))
    print(grap_values(test.gff))
'''
