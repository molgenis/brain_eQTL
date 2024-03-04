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
        wget=1.21.3-1ubuntu1 \
        git=1:2.39.2-1ubuntu1.1 \
        # Languages.
        # Python
        python3=3.11.2-1 \
        python3-dev=3.11.2-1 \
        python3-pip=23.0.1+dfsg-1ubuntu0.1 \
        cython3=0.29.32-2ubuntu2 \
        python-is-python3=3.11.1-3 \
        libpcap-dev=1.10.3-1 \ 
        libpq-dev=15.4-0ubuntu0.23.04.1 \
        # Java
        default-jre=2:1.17-74 \
        # R
        r-base=4.2.2.20221110-2build1 \
        r-cran-nloptr=2.0.3-1 \
        # C++
        g++=4:12.2.0-3ubuntu1 \
        # Other compilers, builders & compressors.
        build-essential=12.9ubuntu3 \
        make=4.3-4.1build1 \
        gfortran=4:12.2.0-3ubuntu1 \
        cmake=3.25.1-1ubuntu1 \
        libssl-dev=3.0.8-1ubuntu1.2 \
        libblas-dev=3.11.0-2 \
        libgsl-dev=2.7.1+dfsg-3 \
        liblapack-dev=3.11.0-2 \
        zlib1g-dev=1:1.2.13.dfsg-1ubuntu4 \
        # Tools
        samtools=1.16.1-1 \
        fastqc=0.11.9+dfsg-6 \
        # Other
        ca-certificates=20230311ubuntu0.23.04.1

# Python packages & requirements.
RUN pip install numpy==1.26.0 multiqc==1.16

# Non-apt Tools
# RMats 4.1.2
RUN mkdir /tmp/init/rmats_build \
    && cd /tmp/init/rmats_build \
    && wget -O rmats-turbo.tar.gz https://github.com/Xinglab/rmats-turbo/releases/download/v4.1.2/rmats_turbo_v4_1_2.tar.gz \
    && tar -xzf rmats-turbo.tar.gz \
    && cd rmats_turbo_v4_1_2 \
    && ./build_rmats \
    && cd .. \
    && mv rmats_turbo_v4_1_2/ /usr/bin \
    && cd /tmp/init/

# STAR2.6.1c
RUN wget https://github.com/alexdobin/STAR/archive/2.6.1c.tar.gz \
    && tar -xzf 2.6.1c.tar.gz \
    && cd STAR-2.6.1c/source \
    && make STAR \
    && cd ../bin/Linux_x86_64 \
    && mv STAR STARlong /usr/bin \
    && cd /tmp/init

# REgtools 1.0.0
RUN wget -O regtools.tar.gz https://github.com/griffithlab/regtools/archive/refs/tags/1.0.0.tar.gz \
    && tar -xzf regtools.tar.gz \
    && cd regtools-1.0.0/ \
    && mkdir build \
    && cd build/ \
    && cmake .. \
    && make \
    && mv regtools /usr/bin \
    && cd /tmp/init

# SRAtools3.0.7
RUN wget -O /tmp/init/sratools.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.7/sratoolkit.3.0.7-ubuntu64.tar.gz \
    && tar -xzvf /tmp/init/sratools.tar.gz -C /usr/bin/ --no-same-owner

# Picard3.1.0
RUN wget -O /usr/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar

################## CLEANUP ######################
RUN cd /usr/bin

# Apt cleanup.
RUN apt-get clean \
    && apt-get autoremove -y

# Build files cleanup
RUN cd /tmp
RUN rm -rf /tmp/init