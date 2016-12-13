#!/usr/bin/env python3
import compute_msa
import make_Datasets
import tempfile as tmp
import os
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#true_header = set(make_Datasets.hash_fasta(eef_inf5_clean).keys())
eef_inf5_polluted = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf5/trimming_test/eef_polluted_20eftu2.mafft"

eef_inf5_clean = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf5/trimming_test/eef_clean_inf5.0.mafft"
eef5_start_header = set(make_Datasets.hash_fasta(eef_inf5_clean).keys())
eef_4_0_path = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf5/trimming_test/eef_inf4.0.mafft"
eef_3_5_path = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf5/trimming_test/eef_inf3.5.mafft"
eef_2_0_path = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf5/trimming_test/eef_inf2.0.mafft"

W04_5_0_path = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf5/trimming_test/OrthologousGroups_I5.0.OG0003142_W04E12.mafft"
W04_5_start_header = set(make_Datasets.hash_fasta(W04_5_0_path).keys())
W04_4_0_path = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf5/trimming_test/OrthologousGroups_I4.0.OG0003238_W04E12.mafft"
W04_3_5_path = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf5/trimming_test/OrthologousGroups_I3.5.OG0003256_W04E12.mafft"
W04_2_0_path = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf5/trimming_test/OrthologousGroups_I2.0.OG0002785_W04E12.mafft"

'''
mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
z = range(50, 70, 10)
x = range(50, 70, 10)
y = range(50, 70, 10)
ax.plot(x, y, z, label='bla')
ax.legend()

plt.show()
'''

resoverlap = range(20, 60, 5)#range(10, 100, 10)
seqoverlap = range(50, 100, 5)

x_values, y_values = [], []
trim_eef_clean, trim_remaining_true_polluted = [],[]
W04E12_inf5, W04E12_overlap, W04E12_inf2_left = [], [], []
eef_inf5, eef_overlap, eef_inf2_left, eef_polluted_overlap = [], [], [], []
for resover in resoverlap:
    resover = (resover / 100)
    for seqover in seqoverlap:
        x_values.append(resover)
        y_values.append(seqover)
        print(resover, seqover)
        eef_5 = compute_msa.msa_operations("trimal -in {} -resoverlap {} -seqoverlap {}".format(eef_inf5_clean, resover, seqover))
        eef_5 = [eef_5[x] for x in range(0, len(eef_5), 2)]
        eef_2 = compute_msa.msa_operations("trimal -in {} -resoverlap {} -seqoverlap {}".format(eef_2_0_path, resover, seqover))
        eef_2 = [eef_2[x] for x in range(0, len(eef_2), 2)]
        eef_polluted = compute_msa.msa_operations("trimal -in {} -resoverlap {} -seqoverlap {}".format(eef_inf5_polluted, resover, seqover))
        eef_polluted = [eef_polluted[x] for x in range(0, len(eef_polluted), 2)]
        eef_polluted_overlap.append(len(set(eef5_start_header) & set(eef_polluted))/len(eef_polluted))
        eef_inf5.append(len(eef_5)/237)
        eef_overlap.append(len(set(eef5_start_header) & set(eef_2))/len(eef_2))
        eef_inf2_left.append(len(eef_2)/410)
        print("eef", len(eef_5)/237, len(set(eef5_start_header) & set(eef_2))/len(eef_2))

        W04_5 = compute_msa.msa_operations("trimal -in {} -resoverlap {} -seqoverlap {}".format(W04_5_0_path, resover, seqover))
        W04_5 = [W04_5[x] for x in range(0, len(W04_5), 2)]
        W04_2 = compute_msa.msa_operations("trimal -in {} -resoverlap {} -seqoverlap {}".format(W04_2_0_path, resover, seqover))
        W04_2 = [W04_2[x] for x in range(0, len(W04_2), 2)]
        W04E12_inf5.append(len(W04_5)/142)
        W04E12_overlap.append(len(set(W04_5_start_header) & set(W04_2))/len(W04_2))
        W04E12_inf2_left.append(len(W04_2)/158)
        print("W04E12", len(W04_5)/142, len(set(W04_5_start_header) & set(W04_2))/len(W04_2))


x = range(0, len(x_values), 1)
plt.plot(x, W04E12_inf5, label="W04E12 remaining in inf5", color="r")
plt.plot(x, W04E12_overlap, label="W04E12 overlap inf5 and inf2", color="r", linestyle="dashed")
plt.plot(x, W04E12_inf2_left, label="W04E12 inf2 left", color="r", linestyle="dotted")
plt.plot(x, eef_inf5, label="EFalpha remaining in inf5", color="b")
plt.plot(x, eef_overlap, label="EFalpha overlap inf5 and inf2", color="b", linestyle="dashed")
plt.plot(x, eef_inf2_left, label="EFalpha inf2 left", color="b", linestyle="dotted")
plt.plot(x, eef_polluted_overlap, label="EFalpha plus wrong sequences", color="g", linestyle="dotted")
plt.show()
