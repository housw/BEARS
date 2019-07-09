# -*- coding: utf-8 -*-

import os
import sys
import pymer
import argparse
import logging
import numpy as np
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp
from sklearn import preprocessing
from sklearn import decomposition
from MulticoreTSNE import MulticoreTSNE as TSNE
import umap

import warnings
warnings.filterwarnings(action="ignore", category=DeprecationWarning, module='sklearn')  # message="divide by zero encountered in divide")


_logger = logging.getLogger("bears")


def compute_PCA_tSNE_coordinates(input_feature_table, output_dir, prefix, threads=20, n_components=100):
    """ tSNE dimension reduction using tSNE, return tSNE coordinates
    """
    output_file = os.path.join(output_dir, prefix + ".tSNE")

    df = pd.read_csv(input_feature_table, sep="\t", index_col=0, header=0)

    # run PCA to take the first 100 dimensions
    df[df == np.inf] = np.nan
    df.fillna(df.mean(), inplace=True)
    if n_components > len(df.columns):
        n_components = len(df.columns) - 1
    pca = decomposition.PCA(n_components=n_components)
    pca.fit(df)
    pca_arr = pca.transform(df)
    pca_df = pd.DataFrame(data=pca_arr, columns=['PCA_{i}'.format(i=n) for n in range(n_components)], index=df.index)

    # run t-SNE
    arr = np.array(pca_df)
    tSNE_coordinates = TSNE(n_jobs=threads).fit_transform(arr)
    tSNE_df = pd.DataFrame(data=tSNE_coordinates, index=df.index, columns=['Feature_tSNE_X', 'Feature_tSNE_Y'])
    tSNE_df.to_csv(output_file, sep="\t", header=True, index=True, float_format='%.6f')

    return output_file


def compute_tSNE_coordinates(freq_file, output_dir, prefix, threads=20, force=False, columns=['tSNE_X', 'tSNE_Y']):
    """ tSNE dimension reduction using tSNE, return tSNE coordinates
    """

    output_tsne_file = os.path.join(output_dir, prefix+"_tSNE.tsv")

    # run t-SNE
    df = pd.read_csv(freq_file, sep="\t", index_col=0, header=0)
    arr = np.array(df)
    tSNE_coordinates = TSNE(n_jobs=threads).fit_transform(arr)
    tSNE_df = pd.DataFrame(data=tSNE_coordinates, index=df.index, columns=columns)
    tSNE_df.to_csv(output_tsne_file, sep="\t", header=True, index=True, float_format='%.6f')

    return output_tsne_file


def compute_feature_UMAP_coordinates(input_feature_table, output_file, n_components=100):
    df = pd.read_csv(input_feature_table, sep="\t", index_col=0, header=0)
    df[df == np.inf] = np.nan
    df.fillna(df.mean(), inplace=True)
    if n_components > len(df.columns):
        n_components = len(df.columns) - 1
    pca = decomposition.PCA(n_components=n_components)
    pca.fit(df)
    pca_arr = pca.transform(df)
    pca_df = pd.DataFrame(data=pca_arr, columns=['PCA_{i}'.format(i=n) for n in range(n_components)], index=df.index)

    # run umap
    arr = np.array(pca_df)
    reducer = umap.UMAP(random_state=42)
    umap_coordinates = reducer.fit_transform(arr)
    umap_df = pd.DataFrame(data=umap_coordinates, index=df.index, columns=['Feature_UMAP_X', 'Feature_UMAP_Y'])
    umap_df.to_csv(output_file, sep="\t", header=True, index=True, float_format='%.6f')
