# -*- coding: utf-8 -*-

import os
import logging
import pandas as pd
from Bio import SeqIO
from abc import ABC, abstractmethod
from blendit.utils.common import emit_file_exist_warning
import warnings
warnings.filterwarnings(action='ignore', category=DeprecationWarning)


_logger = logging.getLogger("BlendIt")


class Binner(ABC):
    """Abstract Binner class.
    Note:
        Do not initialize this abstract class,
        instead, please inherit this abstract class and implement the `bin` and `save` abstract methods.
    """

    name = 'binner'

    def __init__(self, assembly, length_df, embedding_df, output_dir, prefix, min_length=2000,
                 embedding_method='tsne', *args, **kwargs):
            """ To initialize a Binner object
            """
            self.assembly = assembly
            self.length_df = length_df
            self.embedding_df = self._filter_embedding_df(length_df, embedding_df, min_length)
            self.min_length = min_length
            self.output_dir = output_dir
            self.embedding_method = embedding_method
            self.prefix = prefix
            # output handling
            self.cluster_file = prefix + "_{0}_min_{1}_{2}.tsv".format(embedding_method, min_length, self.name)
            self.bin_folder = prefix + "_{0}_min_{1}".format(embedding_method, min_length)
            self.scaffold2bin = prefix + "_{0}_min_{1}_scaffold2bin.tsv".format(embedding_method, min_length)

    def _filter_embedding_df(self, length_df, embedding_df, min_length):
        """ filter embedding dataframe by length
        :return: embedding dataframe
        """
        _logger.info("filtering embedding dataframe to keep contigs >= {l} ...".format(l=min_length))
        drop_length_df = length_df[length_df.Length < min_length]
        excluding_set = set(drop_length_df.index)
        idx_rows = []
        for idx in embedding_df.index:
            if idx in excluding_set:
                idx_rows.append(idx)
        embedding_df = embedding_df.drop(idx_rows, axis=0, inplace=False)
        return embedding_df

    @abstractmethod
    def bin(self, binner_param_dict, **kwargs):
        """bin should take binner_param_dict as input, and should return a cluster dataframe"""
        pass

    def _write_bin_seqs(self, cluster_df, suffix="fa"):

        # output folder
        output_folder = os.path.join(self.output_dir, self.bin_folder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # parse contig file
        contig2seq = {}  # {contig_id : seq_rec}
        for rec in SeqIO.parse(self.assembly, "fasta"):
            Contig_ID = rec.name
            if Contig_ID not in contig2seq:
                contig2seq[Contig_ID] = rec

        # parse cluster file, get sequences of each cluster
        cluster2seq = {}  # {cluster:[seq_rec]}
        #cluster_df = pd.read_csv(cluster_file, sep='\t', header=0, index_col=0)
        for Contig_ID in cluster_df.index:
            Cluster = cluster_df.loc[Contig_ID, 'cluster']
            seq = contig2seq.get(Contig_ID, None)
            if not seq:
                raise Exception("{ID} was not found in input contig!".format(ID=Contig_ID))
            if Cluster not in cluster2seq:
                cluster2seq[Cluster] = [seq]
            else:
                cluster2seq[Cluster].append(seq)

        # write bins
        for cluster_id, seq_list in cluster2seq.items():
            if cluster_id == -1:
                bin_file = os.path.join(output_folder, "unbinned.fa")
            else:
                bin_file = os.path.join(output_folder, "bin_" + str(cluster_id) + ".{0}".format(suffix))
            with open(bin_file, "w") as oh:
                for seq in seq_list:
                    oh.write(">" + seq.name + "\n")
                    oh.write(str(seq.seq) + "\n")
        return output_folder

    def _write_scaffold2bin_file(self, suffix="fa"):
        scaffold2bin_file = os.path.join(self.output_dir, self.scaffold2bin)
        with open(scaffold2bin_file, "w") as oh:
            bin_folder = os.path.join(self.output_dir, self.bin_folder)
            for f in os.listdir(bin_folder):
                if f.endswith(suffix):
                    bin_nr = f.strip("." + suffix)
                    with open(os.path.join(bin_folder, f), "r") as ih:
                        for line in ih:
                            if line.startswith(">"):
                                scaffold = line[1:].split()[0].strip()
                                oh.write(scaffold + "\t" + bin_nr + "\n")
        return scaffold2bin_file

    def __call__(self, binner_param_dict):

        # binning
        _logger.info("clustering using {0} for contigs >= {1} ...".format(self.name, self.min_length))
        cluster_file = self.bin(**binner_param_dict)

        # write bin sequences
        _logger.info("generating bins for contigs >= {0} ...".format(self.min_length))
        bin_folder = self._write_bin_seqs(cluster_file)

        # write scaffold2bin files
        _logger.info("generating scaffold2bin file for contigs >= {0} ...".format(self.min_length))
        scaffold2bin = self._write_scaffold2bin_file()

        return bin_folder, scaffold2bin


class IterBinner(object):

    def __init__(self, binner_cls, embeddings, assembly, contig_length_file, output_dir, prefix, force=False,
                 min_length_x=2000, min_length_y=10000, length_step=1000, binner_param_dict={}):
        self.binner_cls = binner_cls
        self.binner_name = binner_cls.name
        self.embeddings = embeddings
        self.assembly = assembly
        self.contig_length_file = contig_length_file
        self.output_dir = output_dir
        self.prefix = prefix
        self.force = force
        self.min_length_x = min_length_x
        self.min_length_y = min_length_y
        self.length_step = length_step
        self.binner_param_dict = binner_param_dict

    def bin(self, embedding_df, length_df, min_length=2000, embedding_method='tsne', binner_param_dict={'n_jobs': 10}):

        # initialize a binner
        binner = self.binner_cls(assembly=self.assembly, length_df=length_df, embedding_df=embedding_df,
                                output_dir=self.output_dir, prefix=self.prefix, min_length=min_length,
                                embedding_method=embedding_method)
        # generate bin seqs and scaffold2bin file
        bin_folder, scaffold2bin = binner(binner_param_dict)

        return bin_folder, scaffold2bin

    def __call__(self, cpus=10):

        # iterative binning
        scaffold2bin_files = []
        bin_folders = []
        for embedding_file in self.embeddings:
            dimred_method = embedding_file.split('.')[-1]
            embedding_df = pd.read_csv(embedding_file, sep='\t', header=0, index_col=0)
            _logger.info("clustering using {0} based on {1} embedding file ...".format(self.binner_name, dimred_method))
            length_df = pd.read_csv(self.contig_length_file, sep="\t", index_col=0, header=0)
            for min_length in range(self.min_length_x, self.min_length_y+1, self.length_step):
                _logger.info("clustering using {0} with contigs >= {1} ...".format(self.binner_name, min_length))
                bin_folder = os.path.join(self.output_dir, self.prefix + "_{0}_min_{1}".format(dimred_method, min_length))
                scaffold2bin = bin_folder + "_scaffold2bin.tsv"
                try:
                    scaffold2bin = emit_file_exist_warning(filename=scaffold2bin, force=self.force)
                except Exception as e:
                    _logger.info(e)
                    bin_folder, scaffold2bin = self.bin(embedding_df, length_df, min_length=min_length,
                                                        embedding_method=dimred_method,
                                                        binner_param_dict=self.binner_param_dict)
                scaffold2bin_files.append(scaffold2bin)
                bin_folders.append(bin_folder)

        return scaffold2bin_files, bin_folders
