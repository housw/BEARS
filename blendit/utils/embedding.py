# -*- coding: utf-8 -*-

import os
import sys
import logging
import warnings
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn import decomposition
from MulticoreTSNE import MulticoreTSNE as multiTSNE
import umap
from .common import emit_file_exist_warning
from numba import errors
warnings.simplefilter('ignore', category=errors.NumbaWarning)
warnings.simplefilter('ignore', category=errors.NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=errors.NumbaPendingDeprecationWarning)
warnings.simplefilter('ignore', category=errors.NumbaPerformanceWarning)
warnings.filterwarnings(action='ignore', category=DeprecationWarning, module='sklearn')  # message="divide by zero encountered in divide")

_logger = logging.getLogger("BlendIt")


def compute_PCA_components(df, n_components=100):
    """
    :param df: input pandas dataframe
    :param n_components: number of PCA components to keep
    :return: pandas dataframe after PCA
    """

    _logger.info("tranforming with PCA ...")
    if n_components > len(df.columns):
        n_components = len(df.columns) - 1
    pca = decomposition.PCA(n_components=n_components, random_state=1)
    pca.fit(df)
    pca_arr = pca.transform(df)
    pca_df = pd.DataFrame(data=pca_arr, columns=['PCA_{i}'.format(i=n) for n in range(n_components)], index=df.index)

    return pca_df


def compute_tSNE_embeddings(df, threads=20, n_components=3):
    """ tSNE dimension reduction using tSNE, return tSNE embeddings
    """

    # MulticoreTSNE
    _logger.info("embedding with tSNE ...")
    arr = np.array(df)
    tSNE_embeddings = multiTSNE(n_jobs=threads, n_components=n_components, random_state=1).fit_transform(arr)
    tSNE_df = pd.DataFrame(data=tSNE_embeddings, index=df.index,
                           columns=['tSNE_{i}'.format(i=n) for n in range(n_components)])
    '''
    # openTSNE
    from openTSNE import TSNE, TSNEEmbedding, affinity, initialization
    init = initialization.pca(X=arr, n_components=n_components, random_state=1)
    affinities_multiscale_mixture = affinity.Multiscale(data=arr, perplexities=[50, 500], method="approx",
                                                        metric="cosine", n_jobs=threads, random_state=1)
    embedding = TSNEEmbedding(embedding=init, affinities=affinities_multiscale_mixture, random_state=1,
                              negative_gradient_method="fft", n_jobs=threads)
    embedding.optimize(n_iter=250, exaggeration=12, momentum=0.5, inplace=True)
    embedding.optimize(n_iter=750, exaggeration=1, momentum=0.8, inplace=True)
    embedding_multiscale = embedding.view(np.ndarray)
    tSNE_df = pd.DataFrame(data=embedding_multiscale, index=df.index,
                           columns=['tSNE_{i}'.format(i=n) for n in range(n_components)])
    '''

    return tSNE_df


def compute_UMAP_embeddings(df, n_components=3, random_state=42):
    """
    :param df: input pandas dataframe
    :param n_components: number of components to keep after UMAP dimension reduction
    :return: pandas dataframe after UMAP
    """

    _logger.info("embedding with UMAP ...")
    arr = np.array(df)
    reducer = umap.UMAP(random_state=random_state, n_components=n_components)
    umap_embeddings = reducer.fit_transform(arr)
    umap_df = pd.DataFrame(data=umap_embeddings, index=df.index, columns=['UMAP_{i}'.format(i=n) for n in range(n_components)])

    return umap_df


def compute_PCA_tSNE_embeddings(input_feature_table, output_dir, prefix, threads=20, n_components=10, pca_components=100):
    """ dimension reduction using tSNE, return tSNE embeddings
    """
    output_file = os.path.join(output_dir, prefix + ".tsne")

    df = pd.read_csv(input_feature_table, sep="\t", index_col=0, header=0)
    df[df == np.inf] = np.nan
    df.fillna(df.mean(), inplace=True)

    # run PCA to take the first 100 dimensions
    pca_df = compute_PCA_components(df, n_components=pca_components)

    # run t-SNE
    tSNE_df = compute_tSNE_embeddings(pca_df, threads=threads, n_components=n_components)
    tSNE_df.to_csv(output_file, sep="\t", header=True, index=True, float_format='%.6f')

    return output_file


def compute_PCA_UMAP_embeddings(input_feature_table, output_dir, prefix, n_components=10, pca_components=100):
    """ dimension reduction using UMAP, return UMAP embeddings
    """

    output_file = os.path.join(output_dir, prefix + ".umap")

    df = pd.read_csv(input_feature_table, sep="\t", index_col=0, header=0)
    df[df == np.inf] = np.nan
    df.fillna(df.mean(), inplace=True)

    # run PCA to take the first 100 dimensions
    pca_df = compute_PCA_components(df, n_components=pca_components)

    # run UMAP
    umap_df = compute_UMAP_embeddings(pca_df, n_components=n_components)
    umap_df.to_csv(output_file, sep="\t", header=True, index=True, float_format='%.6f')

    return output_file


def compute_embeddings(merged_profile, output_dir, prefix, threads, n_components, pca_components, dimred=['tsne', 'umap', 'both'], force=False):
    """ a wrap function to compute tsne, umap or both  embeddings
    """
    embeddings = []
    if dimred in ('tsne', 'both'):
        # run t-SNE
        _logger.info("computing t-SNE embeddings after PCA ...")
        tsne_file =  os.path.join(output_dir, prefix+"_merged_{}d.tsne".format(n_components))
        try:
            tsne_file = emit_file_exist_warning(filename=tsne_file, force=force)
        except Exception as e:
            _logger.info(e)
            tsne_file = compute_PCA_tSNE_embeddings(merged_profile, output_dir, prefix+"_merged_{}d".format(n_components), threads=threads, n_components=n_components, pca_components=pca_components)
        embeddings.append(tsne_file)
    if dimred in ('umap', 'both'):
        # run UMAP
        _logger.info("computing UMAP embeddings after PCA ...")
        umap_file = os.path.join(output_dir, prefix + "_merged_{}d.umap".format(n_components))
        try:
            umap_file = emit_file_exist_warning(filename=umap_file, force=force)
        except Exception as e:
            _logger.info(e)
            umap_file = compute_PCA_UMAP_embeddings(merged_profile, output_dir, prefix+"_merged_{}d".format(n_components), n_components=n_components, pca_components=pca_components)
        embeddings.append(umap_file)

    return embeddings
