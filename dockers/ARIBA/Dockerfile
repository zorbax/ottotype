FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

MAINTAINER Otoniel Maya <otto94@gmail.com>

ARG BOWTIE2_VERSION=2.2.9
ARG SPADES_VERSION=3.13.1
ARG ARIBA_TAG=master

RUN mkdir /data

RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse ssb && \
    echo "ssb\nssb\n" | passwd ssb && chown ssb:ssb /data

RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y \
    build-essential \
    cd-hit \
    curl \
    git \
    libbz2-dev \
    liblzma-dev \
    mummer \
    python3-dev \
    python3-setuptools \
    python3-pip \
    python3-tk \
    python3-matplotlib \
    unzip \
    wget \
    zlib1g-dev

RUN wget -q http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
  && unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
  && rm -f bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip

RUN wget -q https://github.com/ablab/spades/releases/download/v${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz \
  && tar -zxf SPAdes-${SPADES_VERSION}-Linux.tar.gz \
  && rm -f SPAdes-${SPADES_VERSION}-Linux.tar.gz

# Need MPLBACKEND="agg" to make matplotlib work without X11, otherwise get the error
# _tkinter.TclError: no display name and no $DISPLAY environment variable
ENV ARIBA_BOWTIE2=$PWD/bowtie2-${BOWTIE2_VERSION}/bowtie2 ARIBA_CDHIT=cdhit-est MPLBACKEND="agg"
ENV PATH=$PATH:$PWD/SPAdes-${SPADES_VERSION}-Linux/bin

RUN cd /usr/local/bin && ln -s /usr/bin/python3 python && cd

RUN git clone https://github.com/sanger-pathogens/ariba.git \
  && cd ariba \
  && git checkout ${ARIBA_TAG} \
  && rm -rf .git \
  && python3 setup.py test \
  && python3 setup.py install

CMD ariba

VOLUME /data

USER ssb

WORKDIR /data
