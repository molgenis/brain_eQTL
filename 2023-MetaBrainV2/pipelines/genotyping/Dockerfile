################## BASE IMAGE ######################

FROM ubuntu:23.04

################## METADATA ######################
LABEL base_image="ubuntu:23.04"
LABEL version="1.0.0"
LABEL software="RNA-Seq Genotyping Pipeline"

################## MAINTAINER ######################

LABEL maintainer="Orfeas Gkourlias <o.gkourlias@umcg.nl>"


################## INSTALLATION ######################
ADD . /tmp/init
WORKDIR /tmp/init

ENV SHELL=/bin/bash
ENV LC_ALL=C
ENV PIP_BREAK_SYSTEM_PACKAGES=1
ENV LANG=C.UTF-8
ENV TZ=Europe
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH "/usr/bin/gatk-4.4.0.0:$PATH"
ENV CLASSPATH "/usr/bin/beagle.jar:$CLASSPATH"

# Getting base apps & languages.
RUN apt-get update -y \
    && apt-get upgrade -y \
    && apt-get install -y \
        # Getters & VSC.
        wget=1.21.3-1ubuntu1 \
        git=1:2.39.2-1ubuntu1.1 \
        # Languages.
        # Java
        default-jre=2:1.17-74 \
        # C++
        g++=4:12.2.0-3ubuntu1 \
        # Python
        python-is-python3 \
        # Other compilers, builders & compressors.
        make=4.3-4.1build1 \
        bzip2=1.0.8-5build1 \
        unzip \
        # Tools
        fastqc=0.11.9+dfsg-6 \
        # Other
        liblzma-dev \
        libbz2-dev \
        libclang-dev \
        libcurl4-gnutls-dev \
        ca-certificates=20230311ubuntu0.23.04.1 \
        zlib1g-dev=1:1.2.13.dfsg-1ubuntu4

# Non-apt Tools

# VCFTools 0.1.13 Install
RUN cd /tmp/init \
    && wget -O bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
    && tar -xf bcftools.tar.bz2 \
    && cd bcftools-1.19/ \
    && make \
    && make install \
    && cd /tmp/init

# Gatk4.4.0.0 install
RUN wget -O gatk.zip https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip \
    && unzip gatk.zip -d /usr/bin 
    ## Final path: /usr/bin/gatk-4.4.0.0/gatk

################## CLEANUP ######################
# Apt cleanup.
RUN apt-get clean \
    && apt-get autoremove -y

# Build files cleanup
RUN cd /
RUN rm -rf /tmp/init