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


_logger = logging.getLogger("BlendIt")

base_fwd = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_fwd, base_rev)

def initialize_kmer_dict(ksize=4, canonical=True):
    kmer_dict = {}
    for lett in product("ACGT", repeat=ksize):
        kmer_fwd = "".join(lett)
        if canonical:
            kmer_rev = kmer_fwd.translate(comp_tab)[::-1]
            if kmer_fwd < kmer_rev: 
                kmer = kmer_fwd
            else:
                kmer = kmer_rev
            kmer_dict[kmer] = 0
        else:
            kmer_dict[kmer_fwd] = 0
    return kmer_dict


def get_kmer_count_per_contig(contig, ksize=5, canonical=True):
    """ count kmer for input contig
    """

    contig_name = contig.name
    contig_seq = str(contig.seq).upper()
    kmer_dict = initialize_kmer_dict(ksize, canonical)
    for i in range(len(contig_seq) - ksize + 1):
        kmer_fwd = contig_seq[i:(i+ksize)]
        if 'N' in kmer_fwd: 
            continue
        kmer_rev = kmer_fwd.translate(comp_tab)[::-1]
        if canonical:
            if kmer_fwd < kmer_rev: 
                kmer = kmer_fwd
            else:
                kmer = kmer_rev
            kmer_dict[kmer] += 1
        else:
            kmer_dict[kmer_fwd] += 1
            kmer_dict[kmer_rev] += 1

    return (contig_name, kmer_dict)


def get_kmer_counts_for_contigs(input_contig_file, output_file, k=5, cpus=10, canonical=True):

    # do kmer counting with multiple cores using mp
    contig_iter = SeqIO.parse(input_contig_file, "fasta")
    pool = mp.Pool(processes=cpus)
    func = partial(get_kmer_count_per_contig, ksize=int(k), canonical=canonical)
    result_iter = pool.imap(func, contig_iter)

    # write
    first_contig_name, first_kmer_dict = next(result_iter)
    header = sorted(first_kmer_dict.keys())
    first_count = [str(first_kmer_dict[kmer]) for kmer in header]

    with open(output_file, 'w') as oh:
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
