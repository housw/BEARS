# -*- coding: utf-8 -*-


import os
import logging
import pandas as pd
import hdbscan
from .binner import Binner
from .binner import IterBinner
import warnings
warnings.filterwarnings(action='ignore', category=DeprecationWarning)

_logger = logging.getLogger("BlendIt")


class BinnerHDBSCAN(Binner):

    name = 'hdbscan'

    def bin(self, algorithm='best', alpha=1.0, approx_min_span_tree=True,
            gen_min_span_tree=False, leaf_size=30, metric='euclidean',
            min_cluster_size=5, min_samples=3, p=None, core_dist_n_jobs=4,
            cluster_selection_method='eom', allow_single_cluster=False,
            match_reference_implementation=False, prediction_data=False,**kwargs):
        """ hdbscan clustering
        """

        # hdbscan clustering
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, metric=metric,
                                    alpha=alpha, p=p, algorithm=algorithm, leaf_size=leaf_size,
                                    approx_min_span_tree=approx_min_span_tree, core_dist_n_jobs=core_dist_n_jobs,
                                    cluster_selection_method=cluster_selection_method,
                                    match_reference_implementation=match_reference_implementation,
                                    allow_single_cluster=allow_single_cluster, prediction_data=prediction_data)

        cluster_labels = clusterer.fit_predict(self.embedding_df)
        cluster_df = pd.DataFrame(data=cluster_labels.transpose(),
                                  columns=['cluster'],
                                  index=self.embedding_df.index)
        # write output
        output_cluster_file = os.path.join(self.output_dir, self.prefix + "_hdbscan.tsv")
        cluster_df.to_csv(output_cluster_file, sep="\t", header=True, index=True)

        return cluster_df


def iterative_hdbscan_clustering(embeddings, assembly, contig_length_file, output_dir, prefix, force=False,
                                min_length_x=2000, min_length_y=10000, length_step=1000, binner_params={}):

    binner_param_dict = {'algorithm':'best', 'alpha':1.0, 'approx_min_span_tree':True,
                         'gen_min_span_tree':False, 'leaf_size':30, 'metric':'euclidean',
                         'min_cluster_size':10, 'min_samples':3, 'p':None, 'core_dist_n_jobs':4,
                         'cluster_selection_method':'eom', 'allow_single_cluster':False,
                         'match_reference_implementation':False, 'prediction_data':False}
    binner_param_dict.update(binner_params)

    # initialize IterBinner instance
    iterbinner = IterBinner(binner_cls=BinnerHDBSCAN, embeddings=embeddings, assembly=assembly,
                            contig_length_file=contig_length_file, output_dir=output_dir,
                            prefix=prefix, force=force, min_length_x=min_length_x, min_length_y=min_length_y,
                            length_step=length_step, binner_param_dict=binner_param_dict)

    # run iterbinner
    scaffold2bin_files, bin_folders = iterbinner()

    return scaffold2bin_files, bin_folders
