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
RUN R -e "install.packages('data.table', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('ggplot2', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('gridExtra', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('rmarkdown', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('kableExtra', repos='http://cran.us.r-project.org')"

# Copy workflow scripts
RUN mkdir -p /scripts

COPY scripts/* /scripts/
RUN chmod +x /scripts/*

ENV PATH=${PATH}:/scripts

# Add test data
ADD testdata.tar.gz /testdata/

# Enter home directory
WORKDIR /home
