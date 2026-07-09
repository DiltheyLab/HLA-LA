# work from latest LTS ubuntu release
FROM ubuntu:22.04

# set the environment variables
ENV hla_la_version=1.0.4 \
    samtools_version=1.21 \
    bwa_version=0.7.19 \
    picard_version=2.25.7 \
    bamtools_version=2.5.2

# install required system packages
RUN apt-get update -y && apt-get install --no-install-recommends -y \
    build-essential \
    curl \
    unzip \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libnss-sss \
    libbz2-dev \
    liblzma-dev \
    vim \
    less \
    libcurl4-openssl-dev \
    wget \
    libz-dev \
    openjdk-11-jre \
    libboost-all-dev \
    cmake \
    libjsoncpp-dev \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

# install samtools
WORKDIR /usr/local/bin
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
 && tar -xjf samtools-${samtools_version}.tar.bz2 \
 && cd samtools-${samtools_version} \
 && ./configure \
 && make \
 && make install

# install picard
RUN mkdir -p /usr/local/bin/picard && \
    wget -O /usr/local/bin/picard/picard.jar https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar && \
    chmod 0644 /usr/local/bin/picard/picard.jar && \
    echo '#!/bin/bash\nexec java -jar /usr/local/bin/picard/picard.jar "$@"' > /usr/local/bin/run-picard && \
    chmod +x /usr/local/bin/run-picard

# install bwa
RUN curl -SL https://github.com/lh3/bwa/archive/v${bwa_version}.zip -o v${bwa_version}.zip \
 && unzip v${bwa_version}.zip \
 && cd bwa-${bwa_version} \
 && make \
 && ln -s /usr/local/bin/bwa-${bwa_version}/bwa /usr/local/bin/bwa

# install bamtools
RUN wget https://github.com/pezmaster31/bamtools/archive/v${bamtools_version}.zip \
 && unzip v${bamtools_version}.zip \
 && mkdir -p bamtools-${bamtools_version}/build \
 && cd bamtools-${bamtools_version}/build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local/bin/bamtools-${bamtools_version} .. \
 && make \
 && make install


RUN mkdir -p /usr/local/bin/HLA-LA/bin \
    /usr/local/bin/HLA-LA/src \
    /usr/local/bin/HLA-LA/obj \
    /usr/local/bin/HLA-LA/temp \
    /usr/local/bin/HLA-LA/working \
    /usr/local/bin/HLA-LA/graphs
COPY . /usr/local/bin/HLA-LA/src/
RUN sed -i 's@\$(BAMTOOLS_PATH)/lib64@\$(BAMTOOLS_PATH)/lib@' /usr/local/bin/HLA-LA/src/makefile \
 && make -C /usr/local/bin/HLA-LA/src all BOOST_PATH=/usr/include/boost BAMTOOLS_PATH=/usr/local/bin/bamtools-${bamtools_version}

# modify paths.ini for hla-la
RUN sed -i 's@samtools_bin=@samtools_bin=/usr/local/bin/samtools@' /usr/local/bin/HLA-LA/src/paths.ini \
 && sed -i 's@bwa_bin=@bwa_bin=/usr/local/bin/bwa@' /usr/local/bin/HLA-LA/src/paths.ini \
 && sed -i 's@picard_sam2fastq_bin=.*@picard_sam2fastq_bin=/usr/local/bin/picard/picard.jar@' /usr/local/bin/HLA-LA/src/paths.ini

# set PATH
ENV PATH="/usr/local/bin/HLA-LA/bin:/usr/local/bin/HLA-LA/src:$PATH"

RUN groupadd app && useradd -m app -g app
RUN chown -R app:app /opt /usr/local/bin

USER app

WORKDIR /home/app
