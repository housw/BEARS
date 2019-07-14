FROM shengwei/das_tool:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com


# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends cmake python python-dev python-pip python-wheel python-setuptools \
    libcurl3-dev prodigal mummer wget libhdf5-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# install miniconda and python 3.6
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && rm ~/miniconda.sh && \
    /opt/conda/bin/conda install conda-build python=3.6 && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
ENV PATH=/opt/conda/bin:$PATH





# ---- #
# dRep #
# ---- #




# ------- #
# BlendIt #
# ------- #

# bamcov
WORKDIR /
RUN git clone --recurse-submodules https://github.com/fbreitwieser/bamcov && \
    cd bamcov && make && make test && cd /
ENV PATH=/bamcov:$PATH

# python packages
RUN conda install -c bioconda -y cython numpy pandas scipy scikit-learn biopython && \
    conda clean -ya

# blendit
ENV BlendIt_DIR /BlendIt
WORKDIR /
RUN git clone https://github.com/housw/BlendIt.git && \
    cd ${BlendIt_DIR} && /opt/conda/bin/python setup.py install


# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
