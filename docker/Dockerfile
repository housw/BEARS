FROM shengwei/das_tool:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN add-apt-repository ppa:jonathonf/python-3.6 && \
    apt-get -qq update && \
    apt-get install --no-install-recommends -y ca-certificates libssl-dev libffi-dev python3.6 python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install miniconda
#RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#/bin/bash ~/miniconda.sh -b -p /opt/conda && \
#rm ~/miniconda.sh && \
#/opt/conda/bin/conda clean -tipsy && \
#ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
#ENV PATH=/opt/conda/bin:$PATH





# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

