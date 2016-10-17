# plotting true negatives vs true positives

import sys
from run_command import run_cmd
import matplotlib.pyplot as plt

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
                score_hash[line[0]] = (line[1], line[2])
    return score_hash


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
    for protein in tn_hash:
        o_group = grap_output(protein, ortho_group_file)
        if o_group not in group_hash:
            group_hash[o_group] = 1
        else:
            group_hash[o_group] += 1
    return group_hash




tp_file = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/TPostive_eef_proteins_names_scores_length_3.5.txt"
tn_file = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/TNegative_eef_names_scores_length_3.5.txt"
ortho_group_file = "/home/jgravemeyer/Dropbox/MSc_project/data/orthologousGroups/OrthologousGroups_I3.5.txt"

tp_hash = hash_score_files(tp_file)
tn_hash = hash_score_files(tn_file)

'''
tp_scores = extract_scores(tp_hash, normalize = True)
tn_scores = extract_scores(tn_hash, normalize = True)

sns.distplot(tp_scores, hist=False, rug=True, color="r", label="Scores within cluster")
sns.distplot(tn_scores, hist=False, rug=True, color="b", label="Scores next best BLAST hits")
plt.title("HMM-Score Distributions", size=18, weight="bold")
plt.xlabel("Score", size=14)
plt.ylabel("Density", size=14)
plt.show()
'''

tn_memberships = analise_membership(tn_hash)

name_list = []
count_list = []
for name in tn_memberships:
    name_list.append(name)
    count_list.append(tn_memberships[name])

#ax = plt.bar(range(1, len(name_list)+1), count_list)
#ax.set_xticklabels(name_list)

ax = plt.subplot(111)

width = 0.3
bins = [x - (width/2) for x in range(1, len(name_list)+1)]

ax.bar(bins, count_list, width=width)
ax.set_xticks([x + (width/2) for x in range(1, len(name_list)+1)])
ax.set_xticklabels(name_list, rotation=0, size=12, rotation_mode="anchor", ha="right", weight="bold")
plt.title("EEF-2 Next Best Blast Hits assigned to their Orthologous Groups", size=14, weight="bold")
plt.ylabel("Relative Group Size", size=14)
#plt.xlabel("Number of Protein Sequences", size=14)

plt.show()
