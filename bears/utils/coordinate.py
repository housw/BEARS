# -*- coding: utf-8 -*-

import os
import sys
import pymer
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp
from sklearn import preprocessing
from sklearn import decomposition
from MulticoreTSNE import MulticoreTSNE as TSNE
import umap


def combine_feature_tables(feature_file_list, output_file):
    all_feature_dfs = []
    for feature_file in feature_file_list:
        curr_df = pd.read_csv(feature_file, sep='\t', header=0, index_col=0)
        curr_df.sort_index(inplace=True)
        all_feature_dfs.append(curr_df)

    combined_df = pd.concat(all_feature_dfs, axis=1)
    combined_df.index.name = "Contig_ID"
    combined_df.to_csv(output_file, sep="\t", header=True, index=True)


def compute_feature_tSNE_coordinates(input_feature_table, output_file, threads=20, n_components=100):
    """ tSNE dimension reduction using tSNE, return tSNE coordinates
    """
    df = pd.read_csv(input_feature_table, sep="\t", index_col=0, header=0)

    # run PCA ?
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


def compute_tSNE_coordinates(codon_freq_file, output_dir, prefix, threads=20, force=False, columns=['tSNE_X', 'tSNE_Y']):
    """ tSNE dimension reduction using tSNE, return tSNE coordinates
    """

    output_codon_tsne_file = os.path.join(output_dir, prefix+"_tSNE.tsv")

    # run t-SNE
    df = pd.read_csv(codon_freq_file, sep="\t", index_col=0, header=0)
    arr = np.array(df)
    tSNE_coordinates = TSNE(n_jobs=threads).fit_transform(arr)
    tSNE_df = pd.DataFrame(data=tSNE_coordinates, index=df.index, columns=columns)
    tSNE_df.to_csv(output_codon_tsne_file, sep="\t", header=True, index=True, float_format='%.6f')

    return output_codon_tsne_file



def compute_kmer_tSNE_coordinates(input_kmer_freq, output_file, threads=20):
    """ tSNE dimension reduction using tSNE, return tSNE coordinates
    """

    # run t-SNE
    df = pd.read_csv(input_kmer_freq, sep="\t", index_col=0, header=0)
    arr = np.array(df)
    tSNE_coordinates = TSNE(n_jobs=threads).fit_transform(arr)
    tSNE_df = pd.DataFrame(data=tSNE_coordinates, index=df.index, columns=['KmerFreq_tSNE_X', 'KmerFreq_tSNE_Y'])
    tSNE_df.to_csv(output_file, sep="\t", header=True, index=True, float_format='%.6f')




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


def main():

    # main parser
    parser = argparse.ArgumentParser(description="combine normalized contig feature tables for metagenomic binning")
    parser.add_argument("input_feature_file", nargs='+', help="input normalized feature files")
    parser.add_argument("-p", "--prefix", default='combined_feature_table', help="output prefix [combined_feature_table]")
    parser.add_argument("-o", "--output_dir", help="output directory, default=./", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output file")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    combined_feature_file = os.path.join(args.output_dir, args.prefix+".tsv")
    combined_tSNE_file = os.path.join(args.output_dir, args.prefix+"_tSNE.tsv")
    combined_umap_file = os.path.join(args.output_dir, args.prefix+"_UMAP.tsv")
    if os.path.exists(combined_feature_file):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # combine feature files
    combine_feature_tables(args.input_feature_file, combined_feature_file)

    # compute t-SNE coordinates
    compute_feature_tSNE_coordinates(combined_feature_file, combined_tSNE_file, threads=20)

    # compute UMAP coordinates
    #compute_feature_UMAP_coordinates(combined_feature_file, combined_umap_file)


if __name__ == "__main__":
    main()
