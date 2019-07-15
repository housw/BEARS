=======
BlendIt
=======


.. image:: https://img.shields.io/pypi/v/blendit.svg
        :target: https://pypi.python.org/pypi/blendit

.. image:: https://img.shields.io/travis/housw/blendit.svg
        :target: https://travis-ci.org/housw/blendit

.. image:: https://readthedocs.org/projects/blendit/badge/?version=latest
        :target: https://blendit.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


**BlendIt**: Binning metagenomic contigs via length-dependent Iterative clustering and integration

* Free software: GNU General Public License v3
* Documentation: https://blendit.readthedocs.io.



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

    Options:
        --help  Show this message and exit.

    Commands:
        bin
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


Example
-------




TODO
----

:post subcommands:

- ``blendit post phylo``


:viz subcommands:

- ``blendit viz scatter``

- ``blendit viz tree``


:pipe subcommands:

- ``blendit pipe tuh``

- ``blendit pipe th``

- ``blendit pipe ud``
