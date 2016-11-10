#!/usr/bin/env python3
# plotting true negatives vs true positives

"""
Usage: plotting_tn_tp.py                      -i <DIR> -o <DIR> [--membership]

    Options:
        -h, --help                            show this screen.

        General
        -i, --input_dir <DIR>                 directory with TruePositive and corresponding True Negative Score files
        -o, --out_dir <DIR>                   output directory
        --membership                          command to safe intermediate files
"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from docopt import docopt

from run_command import run_cmd

import_errors = []
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


def grap_output (protein_id, orthologous_file):
    command = ["grep", """{}""".format(protein_id), orthologous_file]
    for line in run_cmd(command=command, wait=False):
        group = line.strip().split(":")[0]
        return group


def hash_score_files(score_file):
    score_hash = {}
    with open(score_file) as sf:
        for line in sf:
            if not line.startswith("#"):
                line = line.strip("\n").split("\t")
                score_hash[line[0]] = round(float(line[1]),3)
    return score_hash, list(score_hash.values())


def extract_scores(score_dict, normalize=False):
    score_list = []
    if not normalize:
        for protein in score_dict:
            score = round(float(score_dict[protein][0]))
            score_list.append(score)
        return score_list
    else:
        for protein in score_dict:
            score = float(score_dict[protein][0]) / float(score_dict[protein][1])
            score_list.append(score)
        return score_list



def analise_membership(score_hash):
    group_hash = {}
    for protein in score_hash:
        protein = protein.strip(">")
        o_group = grap_output(protein, ortho_group_file)
        if o_group not in group_hash:
            group_hash[o_group] = 1
        else:
            group_hash[o_group] += 1
    return group_hash


def find_intersection(m_p, m_n, std_p, std_n):
    a = 1/(2 * std_p ** 2) - 1 / (2 * std_n ** 2)
    b = m_n / (std_n ** 2) - m_p / (std_p ** 2)
    c = m_p ** 2 / (2 * std_p ** 2) - m_n ** 2 / (2 * std_n ** 2) - np.log(std_n / (2*std_p))
    intersection = np.roots([a, b, c])
    sorted_no_negative = sorted([x for x in intersection if not x < 0], reverse=True)
    for point in sorted_no_negative:
        if point < m_p:
            return point


def tp_vs_tn_plot(tp_list, tn_list):
    m_tp, std_tp = np.mean(tp_list), np.std(tp_list)
    m_tn, std_tn = np.mean(tn_list), np.std(tn_list)
    intersection = find_intersection(m_tp, m_tn, std_tp, std_tn)
    if len(tn_list) >= 4:
        single_inter = intersection
    else:
        if std_tp >= (m_tp/3):
            single_inter = m_tp - std_tp
        else:
            single_inter = m_tp - (2*std_tp)
    print(m_tp, std_tp)
    print(single_inter)
    sns.distplot(tp_list, hist=False, rug=True, color="r", label="Scores within cluster")
    sns.distplot(tn_list, hist=False, rug=True, color="b", label="Scores next best BLAST hits")
    plt.axvline(x=single_inter, linewidth=2, color="k")
    plt.title("HMM-Score Distributions", size=18, weight="bold")
    plt.xlabel("Score", size=14)
    plt.ylabel("Density", size=14)
    plt.savefig(out_name)
    plt.close()


def plot_summary(tp_list, tn_list):
    sns.distplot(tp_list, hist=False, color="r")
    sns.distplot(tn_list, hist=False, color="b")
    plt.title("All HMM-Score Distributions", size=18, weight="bold")
    plt.xlabel("Score", size=14)
    plt.ylabel("Density", size=14)
    plt.savefig(summary_plot)


def plot_tn_group_belonging(true_negative_hash):
    tn_memberships = analise_membership(true_negative_hash)

    name_list = []
    count_list = []
    for name in tn_memberships:
        name_list.append(name)
        count_list.append(tn_memberships[name])

    ax = plt.subplot(111)

    width = 0.3
    bins = [x - (width/2) for x in range(1, len(name_list)+1)]

    ax.bar(bins, count_list, width=width)
    ax.set_xticks([x + (width/2) for x in range(1, len(name_list)+1)])
    ax.set_xticklabels(name_list, rotation=30, size=10, rotation_mode="anchor", ha="right", weight="bold")
    plt.title("Next Best Blast Hits assigned to their Orthologous Groups", size=12, weight="bold")
    plt.ylabel("Group Size", size=14)
    plt.savefig(out_name)
    plt.close()


def walk_and_hash_input(input_dir):
    tn_pathes = {}
    tp_pathes = {}
    for subdir, dirs, files in os.walk(input_dir):
        for single_file in files:
            file_path = os.path.join(os.path.abspath(subdir), single_file)
            single_file = ".".join(single_file.strip().split(".")[0:-1]).split("_")
            file_belonging = single_file[-1]
            file_name = "_".join(single_file[:-1])
            if file_belonging == "scores":
                tp_pathes[file_name] = file_path
            elif file_belonging == "TrueNegativeScores":
                tn_pathes[file_name] = file_path
            else:
                print("path Error")
                pass
    if len(tn_pathes) != len(tp_pathes):
        print("not as much TN files as TP files")
    return tn_pathes, tp_pathes


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    input_dir = args['--input_dir']
    output_dir = args['--out_dir']
    membership = args['--membership']
    ortho_group_file = "/home/jgravemeyer/Dropbox/MSc_project/data/orthologousGroups/OrthologousGroups_I3.5.txt"

    tn_files, tp_files = walk_and_hash_input(input_dir)

    summary_plot = os.path.join(output_dir, "summary_TNvsTP.pdf")
    for name in tn_files:
        print(name)
        out_name = os.path.join(output_dir, name + "_TNvsTP.pdf")
        tp_file = tp_files[name]
        tn_file = tn_files[name]
        tp_hash, tp_scores = hash_score_files(tp_file)
        tn_hash, tn_scores = hash_score_files(tn_file)

        if len(tn_scores) > 1 and len(tp_scores) > 1:
            tp_vs_tn_plot(tp_scores, tn_scores)
            #if len(tn_scores) > 15:
            #    plot_summary(tp_scores, tn_scores)

            if membership:
                out_name = os.path.join(output_dir, name + "_membership.pdf")
                plot_tn_group_belonging(tn_hash)
        else:
            print("{} not enough values to plot".format(name))


