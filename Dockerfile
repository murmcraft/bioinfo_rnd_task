# Using GATK image as a base
# docker pull broadinstitute/gatk
FROM broadinstitute/gatk

RUN apt-get update && \
    apt-get install -y vim git

# HTSlib, SAMtools and BCFtools 
RUN mkdir -p /tools 
WORKDIR /tools
RUN git clone git://github.com/samtools/htslib.git && \
    git clone git://github.com/samtools/bcftools.git && \
    cd bcftools && \
    make && \
    make install 

RUN export BCFTOOLS_PLUGINS=/tools/bcftools/plugins

RUN git clone git://github.com/samtools/samtools.git && \
    cd samtools && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install

ENV PATH=${PATH}:/tools

# Enter home directory to start working
WORKDIR /home
