FROM ubuntu:focal-20220531
LABEL Author "Kyler Anderson <kyler_anderson@dfci.harvard.edu>"

# Install basic linux libraries and dev packages
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils build-essential \
      gzip less curl wget file gnupg \
      apt-transport-https ca-certificates \
      zlib1g-dev

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | \
      tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
      gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get -qqy update && \
    apt-get -qy install --fix-missing google-cloud-sdk && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image

ENV DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_PRIORITY=critical
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      libncurses5-dev \
      libbz2-dev \
      liblzma-dev \
      libcurl4-openssl-dev \
      libssl-dev

# Install the suite
ARG SUITE_VERSION=1.18
RUN cd /opt/ && \
    wget https://github.com/samtools/htslib/releases/download/${SUITE_VERSION}/htslib-${SUITE_VERSION}.tar.bz2 && \
    tar xjf htslib-${SUITE_VERSION}.tar.bz2 && \
    cd htslib-${SUITE_VERSION} && \
    ./configure --enable-libcurl && \
    make && \
    make install && \
    cd -

RUN cd /opt/ && \
    wget https://github.com/samtools/samtools/releases/download/${SUITE_VERSION}/samtools-${SUITE_VERSION}.tar.bz2 && \
    tar xjf samtools-${SUITE_VERSION}.tar.bz2 && \
    cd samtools-${SUITE_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd -

RUN cd /opt/ && \
    wget https://github.com/samtools/bcftools/releases/download/${SUITE_VERSION}/bcftools-${SUITE_VERSION}.tar.bz2 && \
    tar xjf bcftools-${SUITE_VERSION}.tar.bz2 && \
    cd bcftools-${SUITE_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd -

# Entrypoint
CMD ["/bin/bash"]