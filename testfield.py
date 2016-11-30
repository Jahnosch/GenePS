#!/usr/bin/env python3
from make_Datasets import hash_fasta, write_hash_to_fasta
from operator import itemgetter
trimed_fasta = "/home/jgravemeyer/Dropbox/MSc_project/data/testing_GenePS/inf3.5/lethal_clusters_out/glb1_3.5.fa_GenePS"
fasta_hash = hash_fasta(trimed_fasta)


def estimate_bin_number(k_guess, k_min, amount_data, forced_bins=None):
    """takes a guess for k = bin number, a specified 'minimum bin number' and the total amount of data points.
    The function checks if more bins then actual data points are demanded and that k_guess is higher then k_min.
    Binning can also be forced but not below #data points."""
    if k_min > amount_data:
        return amount_data
    elif forced_bins:
        return forced_bins
    elif k_guess < k_min:
        return k_min
    else:
        return k_guess

#os.path.join(output_dir, file_name + ".fa_GenePS")
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
    length_hash = {header: len(seq[0]) for header, seq in fasta_hash.items()}
    bin_to_header = bin_sequence_lengths(length_hash)
    bin_to_sequence = {bin_x: fasta_dict[header] for bin_x, header in bin_to_header.items()}
    write_hash_to_fasta(location, bin_to_sequence, ">{}_" + cluster_name + "\n" + "{}\n")
    return length_hash


print(bin_sequence_lengths({"0":0, "4": 4, "12": 12, "16a": 16, "16b": 16, "18": 18, "24": 24, "26": 26, "28": 28}, 3))
#binned_dict_to_fasta(member_dict, fasta_hash, "glb", "adf")


#print(write_length_binned_fasta(fasta_hash, "glb", "/home/jgravemeyer/Desktop/glb.fasta"))
#print(estimate_bin_number(4, 3, 9))


