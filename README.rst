
**WARNING**: This software is still under active development, please wait until the first release.



=======
BlendIt
=======


**BlendIt**: Binning metagenomic contigs via length-dependent Iterative clustering and integration

* Free software: GNU General Public License v3



Installation
------------

:install from github:

First, follow the instructions to install `DAS_Tool <https://github.com/cmks/DAS_Tool>`_, then install `BlendIt`
from github:

::

    $ git clone https://github.com/housw/BlendIt.git && \
      cd BlendIt && pip install -r requirements_dev.txt && python setup.py install

:install using docker:

::

    $ git clone https://github.com/housw/BlendIt.git && \
      cd BlendIt/docker && \
      docker build -t 'blendit:0.1.0' .

:install using singularity:

::

    $ git clone https://github.com/housw/BlendIt.git && \
      cd BlendIt/singularity && \
      sudo singularity build --notest BlendIt.sif Singularity.BlendIt.def


Usage
-----

:commands:

::

    $ blendit --help

    Usage: blendit [OPTIONS] COMMAND [ARGS]...

    Blendit: Binning metagenomic contigs via length-dependent iterative clustering and integration

    The program implements the `profile` and `bin` steps of  metagenomic data
    analysis. You can invoke each step using the  following subcommands:

    -> blendit profile COMMAND [ARGS] ...

        -> blendit profile codon:  run codon usage profiling
        -> blendit profile cov:    run read coverage profiling
        -> blendit profile kmer:   run k-mer frequency profiling

    -> blendit bin COMMAND [ARGS] ...
        -> blendit bin hdbscan:    run hdbscan binning
        -> blendit bin dbscan:     run dbscan binning
        -> blendit bin optics:     run optics binning

    -> blendit pipe COMMAND [ARGS] ...
        -> blendit pipe ph:        run all feature profilings and hdbscan clustering
        -> blendit pipe pd:        run all feature profilings and dbscan clustering

    Options:
        --help  Show this message and exit.

    Commands:
        bin
        pipe
        profile


:profile subcommands:

- ``blendit profile kmer``

::

    $ blendit profile kmer --help

    Usage: blendit profile kmer [OPTIONS] ASSEMBLY

        k-mer frequency profiling

    Options:

        -k, --kmer_size [4|5]           k-mer size  [default: 5]
        -c, --cpus INTEGER              number of cores to use for kmer counting [default: 20]
        -p, --prefix TEXT               output prefix  [default: assembly]
        -o, --output_dir TEXT           output directory  [default: ./blendit_results]
        -f, --force                     force to overwrite the output file
        -l, --loglevel [critical|error|warning|info|debug] [default: debug]
        --version                       Show the version and exit.
        --help                          Show this message and exit.

- ``blendit profile codon``

::

    $ blendit profile codon --help

    Usage: blendit profile codon [OPTIONS] ASSEMBLY

        codon usage profiling

    Options:
        -p, --prefix TEXT               output prefix  [default: assembly]
        -o, --output_dir TEXT           output directory  [default:./blendit_results]
        -f, --force                     force to overwrite the output file
        -l, --loglevel [critical|error|warning|info|debug] [default: debug]
        --version                       Show the version and exit.
        --help                          Show this message and exit.

- ``blendit profile cov``

::

    $ blendit profile cov --help

    Usage: blendit profile cov [OPTIONS] [BAM_FILES]...

        read coverage profiling

    Options:
        -l, --read_length INTEGER       read length for log-scaled transformation [default: 250]
        -p, --prefix TEXT               output prefix  [default: assembly]
        -o, --output_dir TEXT           output directory  [default:./blendit_results]
        -f, --force                     force to overwrite the output file
        -l, --loglevel [critical|error|warning|info|debug] [default: debug]
        --version                       Show the version and exit.
        --help                          Show this message and exit.


:bin subcommands:

- ``blendit bin hdbscan``

::

    $ blendit bin hdbscan --help

    Usage: blendit bin hdbscan [OPTIONS] KMERFREQ_FILE CODONFREQ_FILE DEPTH_FILE CONTIG_LENGTH_FILE ASSEMBLY

    hdbscan binning

    Options:
        -x, --min_length_x INTEGER      minimum contig length threshold x  [default: 2000]
        -y, --min_length_y INTEGER      minimum contig length threshold y  [default: 10000]
        -s, --length_step INTEGER       minimum contig length increasement step [default: 1000]
        -t, --threads INTEGER           maximum number of threads to use when available  [default: 20]
        -d, --dimred [tsne|umap|both]   dimension reduction methods, can be 'tsne', 'umap' or 'both'  [default: both]
        --dimensions INTEGER            number of dimensions to keep for embedding [default: 3]
        --components INTEGER            maximum PCA components to keep  [default:100]
        -p, --prefix TEXT               output prefix  [default: assembly]
        -o, --output_dir TEXT           output directory  [default:./blendit_results]
        -f, --force                     force to overwrite the output file
        -l, --loglevel [critical|error|warning|info|debug] [default: debug]
        --version                       Show the version and exit.
        --help                          Show this message and exit.

- ``blendit bin dbscan``

::

    $ blendit bin dbscan --help

    Usage: blendit bin dbscan [OPTIONS] KMERFREQ_FILE CODONFREQ_FILE DEPTH_FILE CONTIG_LENGTH_FILE ASSEMBLY

    dbscan binning

    Options:
        -x, --min_length_x INTEGER      minimum contig length threshold x  [default: 2000]
        -y, --min_length_y INTEGER      minimum contig length threshold y  [default: 10000]
        -s, --length_step INTEGER       minimum contig length increasement step [default: 1000]
        -t, --threads INTEGER           maximum number of threads to use when available  [default: 20]
        -d, --dimred [tsne|umap|both]   dimension reduction methods, can be 'tsne', 'umap' or 'both'  [default: both]
        --dimensions INTEGER            number of dimensions to keep for embedding [default: 3]
        --components INTEGER            maximum PCA components to keep  [default:100]
        -p, --prefix TEXT               output prefix  [default: assembly]
        -o, --output_dir TEXT           output directory  [default:./blendit_results]
        -f, --force                     force to overwrite the output file
        -l, --loglevel [critical|error|warning|info|debug] [default: debug]
        --version                       Show the version and exit.
        --help                          Show this message and exit.


:pipe subcommands:

- ``blendit pipe ph``

::


    $ blendit pipe ph --help

    Usage: blendit pipe ph [OPTIONS] ASSEMBLY [BAM_FILES]...

    run feature profiling and hdbscan clustering pipeline

    Options:
        -k, --kmer_size [4|5]           k-mer size  [default: 5]
        --kmerfreq_scale_func [none|sqrt|cbrt|log10] k-mer freq scale function  [default: cbrt]
        --codonfreq_scale_func [none|sqrt|cbrt|log10] codon freq scale function  [default: cbrt]
        --cov_scale_func [none|sqrt|cbrt|log10] coverage scale function  [default: log10]
        -x, --min_length_x INTEGER      minimum contig length threshold x  [default: 2000]
        -y, --min_length_y INTEGER      minimum contig length threshold y  [default: 10000]
        -s, --length_step INTEGER       minimum contig length increasement step [default: 1000]
        -t, --threads INTEGER           maximum number of threads to use when available  [default: 20]
        -d, --dimred [tsne|umap|both]   dimension reduction methods, can be 'tsne', 'umap' or 'both'  [default: both]
        --dimensions INTEGER            number of dimensions to keep for embedding [default: 3]
        --components INTEGER            maximum PCA components to keep  [default: 100]
        -l, --read_length INTEGER       read length for log-scaled transformation [default: 250]
        -p, --prefix TEXT               output prefix  [default: assembly]
        -o, --output_dir TEXT           output directory  [default: ./blendit_results]
        -f, --force                     force to overwrite the output file
        -l, --loglevel [critical|error|warning|info|debug] [default: debug]
        -t, --threads INTEGER           maximum number of threads/cpus to use when available  [default: 20]
        --version                       Show the version and exit.
        --help                          Show this message and exit.


Example
-------




TODO
----

:post subcommands:

- ``blendit post phylo``


:viz subcommands:

- ``blendit viz scatter``

- ``blendit viz tree``

