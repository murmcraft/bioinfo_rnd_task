# Using GATK image as a base
# docker pull broadinstitute/gatk
FROM broadinstitute/gatk

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y vim git pandoc && \
    apt-get install -y libcurl4-openssl-dev libxml2-dev libssl-dev libcairo2-dev

# HTSlib, SAMtools and BCFtools 
RUN mkdir -p /tools 
WORKDIR /tools
RUN git clone git://github.com/samtools/htslib.git && \
    git clone git://github.com/samtools/bcftools.git && \
    cd bcftools && \
    make && \
    make install 

ENV BCFTOOLS_PLUGINS=/tools/bcftools/plugins

RUN git clone git://github.com/samtools/samtools.git && \
    cd samtools && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install

ENV PATH=${PATH}:/tools

# Install additional R packages
RUN echo "install.packages('ggplot2', repos='http://cran.us.r-project.org')" | R --no-save
RUN echo "install.packages('data.table', repos='http://cran.us.r-project.org')" | R --no-save
RUN echo "install.packages('rmarkdown', repos='http://cran.us.r-project.org')" | R --no-save
RUN echo "install.packages('kableExtra', repos='http://cran.us.r-project.org')" | R --no-save

# Copy workflow scripts
RUN mkdir -p /scripts

COPY scripts/* /scripts/
RUN chmod +x /scripts/*

ENV PATH=${PATH}:/scripts

# Enter home directory to start working
WORKDIR /home
