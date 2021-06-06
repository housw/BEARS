# -*- coding: utf-8 -*-


import os
import logging
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp
from functools import partial
from itertools import product
from blendit.utils.common import folder_exists
from blendit.utils.common import create_directory
from blendit.utils.common import normalizer
from blendit.utils.common import emit_file_exist_warning
from blendit.utils.kmercounter import kmer_counter


_logger = logging.getLogger("BlendIt")


def get_kmer_counts_for_contigs(input_contig_file, output_file, k=5, cpus=10, canonical=True):

    # do kmer counting with multiple cores using mp
    name_iter = (rec.name for rec in SeqIO.parse(input_contig_file, "fasta"))
    seq_iter = (rec.seq for rec in SeqIO.parse(input_contig_file, "fasta"))
    pool = mp.Pool(processes=cpus)
    func = partial(kmer_counter, ksize=int(k), canonical=canonical)
    result_iter = pool.imap(func, seq_iter)

    # write first record 
    first_contig_name = next(name_iter)
    first_kmer_dict = next(result_iter)
    header = sorted(first_kmer_dict.keys())
    first_count = [str(first_kmer_dict[kmer]) for kmer in header]

    with open(output_file, 'w') as oh:
        oh.write("Contig_ID" +"\t"+ "\t".join(header)+"\n")
        oh.write(first_contig_name +"\t"+ "\t".join(first_count)+"\n")
        # ask for forgiveness than permission
        while True:
            try:
                contig_name = next(name_iter)
                kmer_dict = next(result_iter)
                count = [str(kmer_dict[kmer]) for kmer in header]
                oh.write(contig_name +"\t"+ "\t".join(count)+"\n")
            except IndexError:
                break
            except StopIteration:
                break
    return output_file


def get_kmer_frequencies_for_contigs(input_contig_file, output_dir, prefix, k=5, cpus=10, canonical=True, 
                                    force=False, scale_func='cbrt'):
    """ divide each kmer count by the sum of kmer counts in each row
    """

    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)
    output_kmer_count = os.path.join(output_dir, prefix + "_kmercount.tsv")
    output_kmer_freq = os.path.join(output_dir, prefix + "_kmerfreq.tsv")
    output_kmer_freq_norm = os.path.join(output_dir, prefix + "_kmerfreq_norm.tsv")

    _logger.info("calculating kmer count ...")
    try:
        emit_file_exist_warning(filename=output_kmer_count, force=force)
    except Exception as e:
        _logger.info(e)
        get_kmer_counts_for_contigs(input_contig_file, output_kmer_count, k=k, cpus=cpus, canonical=canonical)
    _logger.info("calculating kmer frequency ...")
    try:
        output_kmer_freq = emit_file_exist_warning(filename=output_kmer_freq, force=force)
    except Exception as e:
        _logger.info(e)
        kmercount = pd.read_csv(output_kmer_count, sep='\t', header=0, index_col=0)
        kmerfreq = kmercount.div(kmercount.sum(axis=1), axis=0)
        kmerfreq.to_csv(output_kmer_freq, sep="\t", header=True, index=True, float_format='%.6f')

    _logger.info("normalizing kmer frequency ...")
    try:
        emit_file_exist_warning(filename=output_kmer_freq_norm, force=force)
    except Exception as e:
        _logger.info(e)
        normalizer(input_freq_file=output_kmer_freq, output_file=output_kmer_freq_norm, scale_func=scale_func)

    return output_kmer_freq_norm
