# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from .binner import Binner
from .binner import IterBinner
#from .binner import generate_bins
#from .binner import scaffolds_to_bins
from sklearn.cluster import DBSCAN
from collections.abc import Iterable
#from blendit.utils.common import scaffolds_to_bins
import warnings
warnings.filterwarnings(action='ignore', category=DeprecationWarning, module='sklearn')  # message="divide by zero encountered in divide")


_logger = logging.getLogger("BlendIt")


class BinnerDBSCAN(Binner):

    name = 'dbscan'

    def bin(self, eps=0.5, min_samples=3, metric='euclidean', metric_params=None, algorithm='auto', leaf_size=30,
            p=None, n_jobs=10, *args, **kwargs):
        """ dbscan clustering
        """
        clusterer = DBSCAN(eps=eps, min_samples=min_samples, metric=metric, metric_params=metric_params,
                           algorithm=algorithm, leaf_size=leaf_size, p=p, n_jobs=n_jobs)
        cluster_labels = clusterer.fit_predict(self.embedding_df)
        cluster_df = pd.DataFrame(data=cluster_labels.transpose(),
                                  columns=['cluster'],
                                  index=self.embedding_df.index)
        # write output
        output_cluster_file = os.path.join(self.output_dir, self.cluster_file)
        cluster_df.to_csv(output_cluster_file, sep="\t", header=True, index=True)

        return cluster_df


def dbscan_clustering(assembly, embedding_df, length_df, output_dir, prefix, min_length=2000, embedding_method='tsne',
                      binner_param_dict={'n_jobs': 10}):

    # initialize a binner
    dbscan = BinnerDBSCAN(assembly, length_df, embedding_df, output_dir, prefix, min_length=min_length,
                 embedding_method=embedding_method)

    # generate bin seqs and scaffold2bin file
    bin_folder, scaffold2bin = dbscan(binner_param_dict)

    return bin_folder, scaffold2bin


def iterative_dbscan_clustering_old(embeddings, assembly, contig_length_file, output_dir, prefix,
                              min_length_x=2000, min_length_y=10000, length_step=1000):

    binner_param_dict = {'eps': 0.5, 'min_samples': 3, 'metric': 'euclidean',
                         'metric_params': None, 'algorithm': 'auto',
                         'leaf_size': 30, 'p': None, 'n_jobs': 10}

    # iterative binning
    scaffold2bin_files = []
    bin_folders = []
    for embedding_file in embeddings:
        dimred_method = embedding_file.split('.')[-1]
        embedding_df = pd.read_csv(embedding_file, sep='\t', header=0, index_col=0)
        _logger.info("clustering using dbscan based on {0} embedding file ...".format(dimred_method))
        length_df = pd.read_csv(contig_length_file, sep="\t", index_col=0, header=0)
        for min_length in range(min_length_x, min_length_y+1, length_step):
            _logger.info("clustering using dbscan with contigs >= {l} ...".format(l=min_length))
            bin_folder, scaffold2bin = dbscan_clustering(assembly, embedding_df, length_df, output_dir, prefix,
                                                         min_length=min_length, embedding_method=dimred_method,
                                                         binner_param_dict=binner_param_dict)
            scaffold2bin_files.append(scaffold2bin)
            bin_folders.append(bin_folder)

    return scaffold2bin_files, bin_folders


def iterative_dbscan_clustering(embeddings, assembly, contig_length_file, output_dir, prefix,
                                min_length_x=2000, min_length_y=10000, length_step=1000):

    binner_param_dict = {'eps': 0.5, 'min_samples': 3, 'metric': 'euclidean',
                         'metric_params': None, 'algorithm': 'auto',
                         'leaf_size': 30, 'p': None, 'n_jobs': 10}

    # initialize IterBinner instance 
    iterbinner = IterBinner(binner_cls=BinnerDBSCAN, embeddings=embeddings, assembly=assembly, 
                            contig_length_file=contig_length_file, output_dir=output_dir, 
                            prefix=prefix, min_length_x=min_length_x, min_length_y=min_length_y, 
                            length_step=length_step, binner_param_dict=binner_param_dict)

    # run 
    scaffold2bin_files, bin_folders = iterbinner()

    return scaffold2bin_files, bin_folders





