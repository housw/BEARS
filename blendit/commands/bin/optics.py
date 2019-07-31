# -*- coding: utf-8 -*-

import os
import logging
import numpy as np
import pandas as pd
from .binner import Binner
from .binner import IterBinner
from sklearn.cluster import OPTICS
import warnings
warnings.filterwarnings(action='ignore', category=DeprecationWarning)


_logger = logging.getLogger("BlendIt")


class BinnerOPTICS(Binner):

    name = 'optics'

    def bin(self, min_samples=5, max_eps=np.inf, metric='euclidean', p=2, metric_params=None, cluster_method='xi',
            eps=None, xi=0.05, predecessor_correction=True, min_cluster_size=None, algorithm='auto',
            leaf_size=30, n_jobs=10, **kwargs):
        """ optics clustering
        """

        # optics clustering
        clusterer = OPTICS(min_samples=min_samples, max_eps=max_eps, metric=metric, p=p, metric_params=metric_params,
                           cluster_method=cluster_method, eps=eps, xi=xi, predecessor_correction=predecessor_correction,
                           min_cluster_size=min_cluster_size, algorithm=algorithm, leaf_size=leaf_size, n_jobs=n_jobs)
        cluster_labels = clusterer.fit_predict(self.embedding_df)
        cluster_df = pd.DataFrame(data=cluster_labels.transpose(),
                                  columns=['cluster'],
                                  index=self.embedding_df.index)
        # write output
        output_cluster_file = os.path.join(self.output_dir, self.prefix + "_optics.tsv")
        cluster_df.to_csv(output_cluster_file, sep="\t", header=True, index=True)

        return cluster_df


def iterative_optics_clustering(embeddings, assembly, contig_length_file, output_dir, prefix, force=False,
                                min_length_x=2000, min_length_y=10000, length_step=1000, binner_params={}):

    binner_param_dict = {'min_samples': 5, 'max_eps': np.inf, 'metric': 'euclidean', 'p': 2,
                         'metric_params': None, 'cluster_method': 'xi', 'eps': None, 'xi': 0.05,
                         'predecessor_correction': True, 'min_cluster_size': None, 'algorithm': 'auto',
                         'leaf_size': 30, 'n_jobs': 10}
    binner_param_dict.update(binner_params)

    # initialize IterBinner instance
    iterbinner = IterBinner(binner_cls=BinnerOPTICS, embeddings=embeddings, assembly=assembly,
                            contig_length_file=contig_length_file, output_dir=output_dir,
                            prefix=prefix, force=force, min_length_x=min_length_x, min_length_y=min_length_y,
                            length_step=length_step, binner_param_dict=binner_param_dict)

    # run iterbinner
    scaffold2bin_files, bin_folders = iterbinner()

    return scaffold2bin_files, bin_folders
