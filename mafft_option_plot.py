import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from plotting_tn_tp import extract_scores, hash_score_files

einsi_old = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/mafft_option_test/einsi_oldgenafpair_test/intermediate_files/next_best_blast_eef_OrthologousGroups_I3.5.OG0001685.txt_scores.txt"
einsi = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/mafft_option_test/einsi_test/intermediate_files/next_best_blast_eef_OrthologousGroups_I3.5.OG0001685.txt_scores.txt"
mafft_default = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/mafft_option_test/mafft_default_test/intermediate_files/next_best_blast_eef_OrthologousGroups_I3.5.OG0001685.txt_scores.txt"
ffnsi = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/mafft_option_test/fftnsi_option_test/intermediate_files/next_best_blast_eef_OrthologousGroups_I3.5.OG0001685.txt_scores.txt"

einsi_scores = hash_score_files(einsi)[0]
einsi_old_scores = hash_score_files(einsi_old)[0]
mafft_scores = hash_score_files(mafft_default)[0]
ffnsi_scores = hash_score_files(ffnsi)[0]

print(np.std(list(einsi_old_scores.values())))
print(np.std(list(einsi_scores.values())))
print(np.std(list(mafft_scores.values())))
print(np.std(list(ffnsi_scores.values())))

sns.distplot(list(einsi_scores.values()), hist=False, rug=False, color="r", label="einsi")
sns.distplot(list(mafft_scores.values()), hist=False, rug=False, color="b", label="mafft_default")
sns.distplot(list(ffnsi_scores.values()), hist=False, rug=False, color="k", label="ffnsi")
sns.distplot(list(einsi_old_scores.values()), hist=False, rug=False, color="g", label="einsi_old")
plt.title("HMM-Score Distributions", size=18, weight="bold")
plt.xlabel("Score", size=14)
plt.ylabel("Density", size=14)
plt.show()
plt.close()

