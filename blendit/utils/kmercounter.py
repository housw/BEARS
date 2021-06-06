# -*- coding: utf-8 -*-


import os
import logging
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections.abc import Iterable
import multiprocessing as mp
from functools import partial
from itertools import product


_logger = logging.getLogger("BlendIt")

base_fwd = "ACGT"
base_rev = "TGCA"
ascii_map = str.maketrans(base_fwd, base_rev)

def initialize_kmer_dict(ksize=4, canonical=True):
    _kmer_dict = {}
    for lett in product("ACGT", repeat=ksize):
        kmer_fwd = "".join(lett)
        if canonical:
            kmer_rev = kmer_fwd.translate(ascii_map)[::-1]
            if kmer_fwd < kmer_rev:
                kmer = kmer_fwd
            else:
                kmer = kmer_rev
            _kmer_dict[kmer] = 0
        else:
            _kmer_dict[kmer_fwd] = 0
    return _kmer_dict


def _count(seq, kmer_dict, ksize=5, canonical=True):
    """ count kmer for input seq
    """
    # convert to upper case
    if isinstance(seq, SeqRecord):
        seq = str(seq.seq).upper()
    else:
        seq = str(seq).upper()
    # count 
    for i in range(len(seq) - ksize + 1):
        kmer_fwd = seq[i:(i+ksize)]
        if 'N' in kmer_fwd:
            continue
        kmer_rev = kmer_fwd.translate(ascii_map)[::-1]
        if canonical:
            if kmer_fwd < kmer_rev:
                kmer = kmer_fwd
            else:
                kmer = kmer_rev
            kmer_dict[kmer] += 1
        else:
            kmer_dict[kmer_fwd] += 1
            kmer_dict[kmer_rev] += 1
    return kmer_dict


def kmer_counter(seq_or_list, ksize=5, canonical=True):
    """ count kmer for input seq_or_list, can ba one seq, 
        or a list of seq records. In both cases, return 
        kmer counts as one record
    """
    # initialize kmer dict
    kmer_dict = initialize_kmer_dict(ksize, canonical)
    # count
    if isinstance(seq_or_list, (list, tuple, Iterable)) and not isinstance(seq_or_list, (str, SeqRecord)):
        for seq in seq_or_list:
            kmer_dict = _count(seq, kmer_dict, ksize, canonical)
    else:
        kmer_dict = _count(seq_or_list, kmer_dict, ksize, canonical)
    return kmer_dict

