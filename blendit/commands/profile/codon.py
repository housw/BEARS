# -*- coding: utf-8 -*-

import os
import logging
from Bio import SeqIO
from collections import OrderedDict
import multiprocessing as mp
from functools import partial
from freqgen import k_mer_frequencies, codon_frequencies, genetic_codes
from blendit.utils.common import folder_exists
from blendit.utils.common import create_directory
from blendit.utils.common import CommandException
from blendit.utils import run_prodigal
from blendit.utils.common import normalizer
from blendit.utils.common import emit_file_exist_warning


_logger = logging.getLogger("BlendIt")


def get_codon_frequency_per_contig(contig2seqs, cpus=10, mode="absolute", genetic_code=11):

    # do codon counting with multiple cores using mp
    seqs_iter = contig2seqs.values()
    pool = mp.Pool(processes=cpus)
    func = partial(codon_frequencies, mode=mode, genetic_code=genetic_code)
    result_iter = pool.imap(func, seqs_iter)

    # get frequencies
    contig2frequencies = {}
    contigs = iter(contig2seqs.keys())
    while True:
        try:
            contig = next(contigs)
            frequency = next(result_iter)
            contig2frequencies[contig] = frequency
        except IndexError:
            break
        except StopIteration:
            break

    return contig2frequencies


def get_genes_per_contig(input_prodigal_nucl_file):
    """ return a dict which contains list of genes for each contig
    """

    contig2seqs = OrderedDict() # {contig1:[seq1, seq2], contig2:[seq1, seq2]}

    for gene in SeqIO.parse(input_prodigal_nucl_file, "fasta"):
        header = gene.name
        contig = "_".join(header.split("_")[0:-1])
        seq = gene.seq
        # don't use genes if there are 'N's
        #if "N" in seq:
        #    print("[WARNING]: {header} contains 'N' in the seq, ignored!".format(header=header))
        #    continue
        if len(seq)  % 3 != 0:
            _logger.warn("the length of {header} is not divisible by 3, ignored!".format(header=header))
            continue
        if contig not in contig2seqs:
            contig2seqs[contig] = [seq]
        else:
            contig2seqs[contig].append(seq)

    return contig2seqs


def get_codon_frequencies_for_contigs(input_contig_file, output_dir, prefix, genetic_code=11, force=False,
                                      scale_func='cbrt', cpus=10):

    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)
    # prodigal files
    prodigal_file = os.path.join(output_dir, prefix + ".gbk")
    prodigal_gene = os.path.join(output_dir, prefix + ".gene")
    prodigal_prot = os.path.join(output_dir, prefix + ".prot")
    output_codon_freq = os.path.join(output_dir, prefix + "_codonfreq.tsv")
    output_codon_freq_norm = os.path.join(output_dir, prefix + "_codonfreq_norm.tsv")

    _logger.info("predicting coding genes using prodigal ...")
    try:
        emit_file_exist_warning(filename=prodigal_prot, force=force)
    except Exception as e:
        _logger.info(e)
        prodigal_file, prodigal_gene, prodigal_prot = \
            run_prodigal(input_contig_file, prefix=prefix, output_dir=output_dir, output_fmt='gbk', flags=['m'], force=force)

    _logger.info("calculating codon frequency ...")
    try:
        emit_file_exist_warning(filename=output_codon_freq, force=force)
    except Exception as e:
        _logger.info(e)
        contig2genes = get_genes_per_contig(prodigal_gene)
        contig2frequencies = get_codon_frequency_per_contig(contig2seqs=contig2genes, cpus=cpus,
                                                            mode="absolute", genetic_code=genetic_code)
        with open(output_codon_freq, "w") as oh:
            oh.write("Contig_ID"+"\t"+"\t".join(genetic_codes[genetic_code])+"\n")
            for contig, frequency_dict in contig2frequencies.items():
                line = [contig]
                for codon in genetic_codes[genetic_code]:
                    _frequency = frequency_dict.get(codon, 0)
                    frequency = "{:.6f}".format(_frequency)
                    line.append(frequency)
                oh.write("\t".join(line)+"\n")

    _logger.info("normalizing codon frequency ...")
    try:
        emit_file_exist_warning(filename=output_codon_freq_norm, force=force)
    except Exception as e:
        _logger.info(e)
        output_codon_freq_norm = normalizer(input_freq_file=output_codon_freq, output_file=output_codon_freq_norm,
                                            scale_func=scale_func)

    return output_codon_freq_norm
