# Using GATK image as a base
# docker pull broadinstitute/gatk
FROM broadinstitute/gatk

RUN apt-get update && \
    apt-get install -y vim git pandoc libcurl4-openssl-dev libxml2-dev

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

# Install additional R packages
RUN echo "install.packages('rmarkdown', repos='http://cran.us.r-project.org')" | R --no-save
RUN echo "install.packages('kableExtra', repos='http://cran.us.r-project.org')" | R --no-save


# Enter home directory to start working
WORKDIR /home
