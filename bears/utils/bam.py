# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn import preprocessing
from MulticoreTSNE import MulticoreTSNE as TSNE


def get_depth_per_bam_file(bam_file, output_depth, min_read_len=30, min_MQ=0, min_BQ=0):

    # command for bamcov
    cmd = ['bamcov', '--output', output_depth,
                     '--min-read-len', str(min_read_len),
                     '--min-MQ', str(min_MQ),
                     '--min-BQ', str(min_BQ),
                     bam_file]
    print("[bamcov] {c}".format(c=" ".join(cmd)))
    try:
        ret_code = subprocess.check_call(cmd, shell=False)
    except Exception as e:
        print("It seems bamcov call failed, the error message is: {e}".format(e=e))
        sys.exit(1)


def calculate_contig_depth_from_bam_files(input_bam_file_list, output_dir, output_prefix):

    depth_files = []
    sample_names = []

    # run bamcov for each bam file
    for bam_file in input_bam_file_list:
        basename = os.path.basename(bam_file)
        filestem = os.path.splitext(basename)[0]
        depth_file = os.path.join(output_dir, filestem+"_bamcov.tsv")
        get_depth_per_bam_file(bam_file, depth_file)
        sample_names.append(filestem)
        depth_files.append(depth_file)

    # write contig length file
    output_length_file = os.path.join(output_dir, output_prefix+"_length.tsv")
    first_df = pd.read_csv(depth_files[0], sep='\t', header=0, index_col=0)
    length_df = first_df[['endpos']]
    length_df.rename(columns={'endpos':'Length'}, inplace=True)
    length_df.index.name = "Contig_ID"
    length_df.to_csv(output_length_file, sep='\t', header=True, index=True)

    # merge depth profiles and write depth file
    output_depth_file = os.path.join(output_dir, output_prefix+"_depth.tsv")
    all_depth_dfs = []
    shape = None
    for i, depth_file in enumerate(depth_files):
        curr_depth_df = pd.read_csv(depth_file, sep='\t', header=0, index_col=0).meandepth.rename(sample_names[i])
        if not shape:
            shape = curr_depth_df.shape
        else:
            assert curr_depth_df.shape == shape, "input depth files has different dimensions!"
        all_depth_dfs.append(curr_depth_df)
    depth_df = pd.concat(all_depth_dfs, axis=1)
    depth_df.index.name="Contig_ID"
    depth_df.to_csv(output_depth_file, sep="\t", header=True, index=True)


def normalization_contig_depth(input_length_file, input_depth_file, output_dir, output_prefix, scale_func=np.log10, read_length=250):
    """ 1) add one read to each contig as prior
        2) scale the data with numpy function
        3) normalize the data using Normalizer
    """
    length_df = pd.read_csv(input_length_file, sep="\t", index_col=0, header=0)
    length_df.sort_index(inplace=True)
    depth_df = pd.read_csv(input_depth_file, sep="\t", index_col=0, header=0)
    depth_df.sort_index(inplace=True)

    # add 1 read prior to depth_df, log-scale the data
    prior_df = read_length / length_df.Length
    for col in depth_df.columns:
        depth_df[col] += prior_df
    log_depth_df = depth_df.apply(func=scale_func)

    # scale each value using sklearn.preprocessing.Normalizer()
    scaler = preprocessing.Normalizer()
    scaled_log_df = scaler.fit_transform(log_depth_df)
    scaled_log_df = pd.DataFrame(scaled_log_df, columns=log_depth_df.columns, index=log_depth_df.index)
    output_norm_depth_file = os.path.join(output_dir, output_prefix+".tsv")
    scaled_log_df.to_csv(output_norm_depth_file, sep="\t", header=True, index=True)


def compute_depth_tSNE_coordinates(input_depth_file, output_dir, output_prefix, threads=20):
    """ tSNE dimension reduction using tSNE, return tSNE coordinates
    """

    # run t-SNE
    df = pd.read_csv(input_depth_file, sep="\t", index_col=0, header=0)
    arr = np.array(df)
    tSNE_coordinates = TSNE(n_jobs=threads).fit_transform(arr)
    tSNE_df = pd.DataFrame(data=tSNE_coordinates, index=df.index, columns=['Depth_tSNE_X', 'Depth_tSNE_Y'])
    output_file = os.path.join(output_dir, output_prefix+".tsv")
    tSNE_df.to_csv(output_file, sep="\t", header=True, index=True)


def main():

    # main parser
    parser = argparse.ArgumentParser(description="get depth of input contigs")
    parser.add_argument("input_bam_files", nargs='+', help="input coordinate sorted bam files")
    parser.add_argument("-p", "--prefix", default='contig_bamcov', help="output prefix [contig_depth]")
    parser.add_argument("-o", "--output_dir", help="output directory, default=./", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output file")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        sys.stderr.write("\nERROR: Not enough parameters were provided, please refer to the usage.\n")
        sys.stderr.write(parser.format_help())
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    output_depth_file = os.path.join(args.output_dir, args.prefix+"_depth.tsv")
    output_length_file = os.path.join(args.output_dir, args.prefix+"_length.tsv")
    output_norm_file = os.path.join(args.output_dir, args.prefix+"_depth_norm.tsv")
    output_norm_tSNE_file = os.path.join(args.output_dir, args.prefix+"_depth_norm_tSNE.tsv")
    if os.path.exists(output_depth_file):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # calculate depth
    calculate_contig_depth_from_bam_files(args.input_bam_files, args.output_dir, args.prefix)

    # normalization
    normalization_contig_depth(output_length_file, output_depth_file, args.output_dir, args.prefix+"_depth_norm", read_length=250)

    # compute tSNE
    compute_depth_tSNE_coordinates(output_norm_file, args.output_dir, args.prefix+"_depth_norm_tSNE", threads=20)

if __name__ == "__main__":
    main()
