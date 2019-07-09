# -*- coding: utf-8 -*-


import os
import pymer
import logging
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp
from bears.utils.common import folder_exists
from bears.utils.common import create_directory
from bears.utils.common import CommandException


_logger = logging.getLogger("bears")


def get_kmer_count_per_contig(contig, ksize=5):

    contig_name = contig.name
    contig_seq = str(contig.seq)

    # initialize ExactKmerCounter
    kc = pymer.ExactKmerCounter(ksize)
    kc.consume(contig_seq)
    kmer_dict = kc.to_dict(sparse=False)

    return (contig_name, kmer_dict)


def get_kmer_counts_for_contigs(input_contig_file, output_dir, prefix, k=5, cpus=10, force=False):

    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)
    output_kmer_freq = os.path.join(output_dir, prefix + ".kmerfreq")

    if os.path.exists(output_kmer_freq):
        if not force:
            err_msg = "{0} exists, use --force if you want to re-generate kmer frequency".format(output_kmer_freq)
            _logger.error(err_msg)
            raise CommandException(err_msg)
        else:
            warn_msg = "re-generate kmer frequency, {0} will be over-writen".format(output_kmer_freq)
            _logger.warn(warn_msg)


    # do kmer count with multiple cores using mp
    contig_iter = SeqIO.parse(input_contig_file, "fasta")
    pool = mp.Pool(processes=cpus)
    result_iter = pool.imap(get_kmer_count_per_contig, contig_iter)

    # write
    first_contig_name, first_kmer_dict = next(result_iter)
    header = sorted(first_kmer_dict.keys())
    first_count = [str(first_kmer_dict[kmer]) for kmer in header]

    with open(output_kmer_freq, 'w') as oh:
        oh.write("Contig_ID" +"\t"+ "\t".join(header)+"\n")
        oh.write(first_contig_name +"\t"+ "\t".join(first_count)+"\n")
        # ask for forgiveness than permission
        while True:
            try:
                contig_name, kmer_dict = next(result_iter)
                count = [str(kmer_dict[kmer]) for kmer in header]
                oh.write(contig_name +"\t"+ "\t".join(count)+"\n")
            except IndexError:
                break
            except StopIteration:
                break
    return output_kmer_freq


def get_kmer_frequencies_for_contigs(input_kmer_count_file, output_dir, output_prefix):
    """ divide each kmer count by the sum of kmer counts in each row
    """

    kmercount = pd.read_csv(input_kmer_count_file, sep='\t', header=0, index_col=0)
    kmerfreq = kmercount.div(kmercount.sum(axis=1), axis=0)
    output_kmerfreq_file = os.path.join(output_dir, output_prefix+".tsv")
    kmerfreq.to_csv(output_kmerfreq_file, sep="\t", header=True, index=True, float_format='%.6f')
