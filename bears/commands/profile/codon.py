# -*- coding: utf-8 -*-

import os
import logging
from Bio import SeqIO
from freqgen import k_mer_frequencies, codon_frequencies, genetic_codes
from bears.utils.common import folder_exists
from bears.utils.common import create_directory
from bears.utils.common import CommandException


_logger = logging.getLogger("bears")


def get_codon_frequency_per_contig(contig2seqs):

    contig2frequencies = {} # {contig1:{'ATG':001, ...}, contig2:{'ATG':002, ... }}

    for contig, seqs in contig2seqs.items():
        frequency = codon_frequencies(seq=seqs, mode="absolute", genetic_code=11)
        contig2frequencies[contig] = frequency

    return contig2frequencies


def get_genes_per_contig(input_prodigal_nucl_file):
    """ return a dict which contains list of genes for each contig
    """

    contig2seqs = {} # {contig1:[seq1, seq2], contig2:[seq1, seq2]}

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


def get_codon_frequencies_for_contigs(input_prodigal_nucl_file, output_dir, prefix,
                                          genetic_code=11, force=False):

    # output file handling
    if not folder_exists(output_dir):
        create_directory(output_dir)
    output_codon_freq = os.path.join(output_dir, prefix + ".codonfreq")
    if os.path.exists(output_codon_freq):
        if not force:
            err_msg = "{0} exists, use --force if you want to re-generate codon frequency".format(output_codon_freq)
            _logger.error(err_msg)
            raise CommandException(err_msg)
        else:
            warn_msg = "re-generate codon frequency, {0} will be over-writen".format(output_codon_freq)
            _logger.warn(warn_msg)

    # get seqs per contig
    contig2genes = get_genes_per_contig(input_prodigal_nucl_file)

    # get frequencies
    contig2frequencies = get_codon_frequency_per_contig(contig2genes)

    # write to output
    with open(output_codon_freq, "w") as oh:
        oh.write("Contig_ID"+"\t"+"\t".join(genetic_codes[genetic_code])+"\n")
        for contig, frequency_dict in contig2frequencies.items():
            line = [contig]
            for codon in genetic_codes[genetic_code]:
                frequency = frequency_dict.get(codon, 0)
                line.append(str(frequency))
            oh.write("\t".join(line)+"\n")

    return output_codon_freq
