################## BASE IMAGE ######################

FROM ubuntu:23.04

################## METADATA ######################
LABEL base_image="ubuntu:23.04"
LABEL version="1.0.0"
LABEL software="RNA-Seq Alignment Pipeline"

################## MAINTAINER ######################

LABEL maintainer="Orfeas Gkourlias <o.gkourlias@umcg.nl>, Joost Bakker <j.bakker@umcg.nl>"


################## INSTALLATION ######################
ADD . /tmp/init
WORKDIR /tmp/init

ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV PIP_BREAK_SYSTEM_PACKAGES=1
ENV LANG=C.UTF-8
ENV TZ=Europe
ENV DEBIAN_FRONTEND=noninteractive

# Getting base apps & languages.
RUN apt-get update -y \
    && apt-get upgrade -y \
    && apt-get install -y \
        # Getters & VSC.
        wget \
        curl \
        git \
        # Languages.
        # Python
        python3 \
        python3-dev \
        python3-pip \
        cython3 \
        python-is-python3 \
        libpcap-dev \ 
        libpq-dev \
        # Java
        default-jre \
        # R
        r-base \
        r-cran-nloptr \
        # C++
        g++ \
        # Other compilers, builders & compressors.
        build-essential \
        make \
        gfortran \
        cmake \
        libssl-dev \
        libblas-dev \
        libgsl-dev \
        liblapack-dev \
        zlib1g-dev \
        # Tools
        samtools \
        fastqc \
        # Other
        ca-certificates

# Language packages
# R packages & requirements.
#RUN R --no-echo -e 'install.packages(c("RcppParallel", "StanHeaders", "rstan"), dependencies = TRUE)'

# Python packages & requirements.
RUN pip install numpy multiqc

# Non-apt Tools

# Nextflow
RUN wget -qO- https://get.nextflow.io > /usr/bin/nextflow \
    && chmod +x /usr/bin/nextflow \
    && nextflow

# RMats install
RUN mkdir /rmats_build \
    && cd /rmats_build \
    && git clone https://github.com/Xinglab/rmats-turbo.git \
    && cd rmats-turbo \
    && ./build_rmats \
    && mv rmats.py bamtools /usr/bin \
    && cd /tmp/init/

# STAR2.7.4a install
RUN wget https://github.com/alexdobin/STAR/archive/2.6.1c.tar.gz \
    && tar -xzf 2.6.1c.tar.gz \
    && cd STAR-2.6.1c/source \
    && make STAR \
    && cd ../bin/Linux_x86_64 \
    && mv STAR STARlong /usr/bin \
    && cd /tmp/init

# REgtools install
RUN git clone https://github.com/griffithlab/regtools \
    && cd regtools/ \
    && mkdir build \
    && cd build/ \
    && cmake .. \
    && make \
    && mv regtools /usr/bin \
    && cd /tmp/init

# SRAtools
RUN wget -O /tmp/init/sratools.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.7/sratoolkit.3.0.7-ubuntu64.tar.gz \
    && tar -xzvf /tmp/init/sratools.tar.gz -C /usr/bin/ --no-same-owner

# Picard
RUN wget -O /usr/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar

################## CLEANUP ######################
RUN cd /usr/bin

# Apt cleanup.
RUN apt-get clean \
    && apt-get autoremove -y

# Build files cleanup
RUN cd /tmp
# RUN rm -rf /tmp/init