# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import numpy as np
import pandas as pd
import hdbscan
from Bio import SeqIO


def hdbscan_clustering(input_tSNE_file, output_dir, output_prefix, excluding_contig_file=None, min_cluster_size=10):

    # read in tSNE file
    tSNE_df =  pd.read_csv(input_tSNE_file, sep='\t', header=0, index_col=0)

    # exclude contigs
    if excluding_contig_file:
        excluding_set = None
        try:
            contigs = pd.read_csv(excluding_contig_file, header=None)
            contigs.columns = ['Contig_ID']
            excluding_set = set(contigs.Contig_ID)
            idx_rows = []
            for idx in tSNE_df.index:
                if idx in excluding_set:
                    idx_rows.append(idx)
            tSNE_df.drop(idx_rows, axis=0, inplace=True)
        except Exception as e:
            print("WARNING: excluding contigs failed: {e}".format(e=e))

    # hdbscan clustering
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
    #clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, prediction_data=True)
    #soft_membership_vec = hdbscan.all_points_membership_vectors(clusterer=clusterer)
    cluster_labels = clusterer.fit_predict(tSNE_df)
    cluster_df = pd.DataFrame(data=cluster_labels.transpose(),
                              columns=['cluster'],
                              index=tSNE_df.index)
    # write output
    output_hdbscan_file = os.path.join(output_dir, output_prefix+"_hdbscan.tsv")
    cluster_df.to_csv(output_hdbscan_file, sep="\t", header=True, index=True)


def generate_bins(hdbscan_cluster_file, input_contig_file, output_dir, bin_folder):

    # output folder
    output_folder = os.path.join(output_dir, bin_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # parse contig file
    contig2seq = {} # {contig_id : seq_rec}
    for rec in SeqIO.parse(input_contig_file, "fasta"):
        Contig_ID = rec.name
        if Contig_ID not in contig2seq:
            contig2seq[Contig_ID] = rec

    # parse hdbscan cluster, get sequences of each cluster
    cluster2seq = {} # {cluster:[seq_rec]}
    cluster_df =  pd.read_csv(hdbscan_cluster_file, sep='\t', header=0, index_col=0)
    for Contig_ID in cluster_df.index:
        Cluster = cluster_df.loc[Contig_ID, 'cluster']
        seq = contig2seq.get(Contig_ID, None)
        if not seq:
            raise Exception("{ID} was not found in input contig!".format(ID=Contig_ID))
        if Cluster not in cluster2seq:
            cluster2seq[Cluster] = [seq]
        else:
            cluster2seq[Cluster].append(seq)

    # write bins
    for cluster_id, seq_list in cluster2seq.items():
        if cluster_id == -1:
            bin_file = os.path.join(output_folder, "unbinned.fa")
        else:
            bin_file = os.path.join(output_folder, "bin_" + str(cluster_id) + ".fa")
        with open(bin_file, "w") as oh:
            for seq in seq_list:
                oh.write(">"+seq.name+"\n")
                oh.write(str(seq.seq)+"\n")


def main():

    # main parser
    parser = argparse.ArgumentParser(description="cluster contigs using hdbscan with t-SNE coordinates")
    parser.add_argument("input_tSNE_file", help="input coordinate file of t-SNE")
    parser.add_argument("input_contig_file", help="input metagenome contig file")
    parser.add_argument("-e", "--excluding_contig_file", help="file contains a list of contig names which will be excluded in the following binning step")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--output_dir", help="output directory, default=./", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output file")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    if not args.prefix:
        basename = os.path.basename(args.input_tSNE_file)
        args.prefix = os.path.splitext(basename)[0]
    output_hdbscan_file = os.path.join(args.output_dir, args.prefix+"_hdbscan.tsv")
    if os.path.exists(output_hdbscan_file):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # hdbscan
    hdbscan_clustering(args.input_tSNE_file, args.output_dir, args.prefix, excluding_contig_file=args.excluding_contig_file, min_cluster_size=10)

    # generate bins
    bin_folder= args.prefix + "_bins"
    generate_bins(output_hdbscan_file, args.input_contig_file, args.output_dir, bin_folder)


if __name__ == "__main__":
    main()
