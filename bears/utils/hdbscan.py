# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging
import subprocess
import numpy as np
import pandas as pd
import hdbscan
from Bio import SeqIO
from collections.abc import Iterable


_logger = logging.getLogger("bears")


def hdbscan_clustering(input_tSNE_file, output_dir, prefix, excluding_contig_file_or_list=None, min_cluster_size=10):

    # read in tSNE file
    tSNE_df =  pd.read_csv(input_tSNE_file, sep='\t', header=0, index_col=0)

    # exclude contigs
    excluding_set = {}
    if isinstance(excluding_contig_file_or_list, Iterable):
        excluding_set = set(excluding_contig_file_or_list)
    elif sinstance(excluding_contig_file_or_list, str) and os.path.isfile(excluding_contig_file_or_list):
        try:
            contigs = pd.read_csv(excluding_contig_file_or_list, header=None)
            contigs.columns = ['Contig_ID']
            excluding_set = set(contigs.Contig_ID)
        except Exception as e:
            _logger.error("read excluding_contig_file_or_list failed: {0}".format(e))
    else:
        _logger.error("the excluding_contig_file_or_list has to be a file or a list!")

    if excluding_set:
        idx_rows = []
        for idx in tSNE_df.index:
            if idx in excluding_set:
                idx_rows.append(idx)
        try:
            tSNE_df.drop(idx_rows, axis=0, inplace=True)
        except Exception as e:
            _logger.error("excluding contigs failed: {e}".format(e=e))

    # hdbscan clustering
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
    #clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, prediction_data=True)
    #soft_membership_vec = hdbscan.all_points_membership_vectors(clusterer=clusterer)
    cluster_labels = clusterer.fit_predict(tSNE_df)
    cluster_df = pd.DataFrame(data=cluster_labels.transpose(),
                              columns=['cluster'],
                              index=tSNE_df.index)
    # write output
    output_hdbscan_file = os.path.join(output_dir, prefix+"_hdbscan.tsv")
    cluster_df.to_csv(output_hdbscan_file, sep="\t", header=True, index=True)

    return output_hdbscan_file


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

