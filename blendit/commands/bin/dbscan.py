# -*- coding: utf-8 -*-

import os
import logging
import pandas as pd
from .binner import Binner
from .binner import IterBinner
from sklearn.cluster import DBSCAN
import warnings
warnings.filterwarnings(action='ignore', category=DeprecationWarning)


_logger = logging.getLogger("BlendIt")


class BinnerDBSCAN(Binner):

    name = 'dbscan'

    def bin(self, eps=0.5, min_samples=3, metric='euclidean', metric_params=None, algorithm='auto', leaf_size=30,
            p=None, n_jobs=10, *args, **kwargs):
        """ dbscan clustering
        """
        clusterer = DBSCAN(eps=eps, min_samples=min_samples, metric=metric, metric_params=metric_params,
                           algorithm=algorithm, leaf_size=leaf_size, p=p, n_jobs=n_jobs)
        cluster_labels = clusterer.fit_predict(self.embedding_df)
        cluster_df = pd.DataFrame(data=cluster_labels.transpose(),
                                  columns=['cluster'],
                                  index=self.embedding_df.index)
        # write output
        output_cluster_file = os.path.join(self.output_dir, self.cluster_file)
        cluster_df.to_csv(output_cluster_file, sep="\t", header=True, index=True)

        return cluster_df


def iterative_dbscan_clustering(embeddings, assembly, contig_length_file, output_dir, prefix, force=False,
                                min_length_x=2000, min_length_y=10000, length_step=1000, binner_params={}):

    binner_param_dict = {'eps': 0.5, 'min_samples': 3, 'metric': 'euclidean',
                         'metric_params': None, 'algorithm': 'auto',
                         'leaf_size': 30, 'p': None, 'n_jobs': 10}
    binner_param_dict.update(binner_params)

    # initialize IterBinner instance
    iterbinner = IterBinner(binner_cls=BinnerDBSCAN, embeddings=embeddings, assembly=assembly,
                            contig_length_file=contig_length_file, output_dir=output_dir,
                            prefix=prefix, force=force, min_length_x=min_length_x, min_length_y=min_length_y,
                            length_step=length_step, binner_param_dict=binner_param_dict)

    # run iterbinner
    scaffold2bin_files, bin_folders = iterbinner()

    return scaffold2bin_files, bin_folders
