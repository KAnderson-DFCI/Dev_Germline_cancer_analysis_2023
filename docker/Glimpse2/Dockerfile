FROM kylera/samtools-suite:latest

LABEL org.opencontainers.image.created="2024-10-30"
LABEL org.opencontainers.image.version="2.0.1"
LABEL org.opencontainers.image.title="kylera-glimpse"


WORKDIR /opt/

# Install required packages
RUN apt-get update && apt-get install -y build-essential libbz2-dev libcurl4-openssl-dev autoconf libssl-dev wget zlib1g-dev liblzma-dev libdeflate-dev git

# Download and build boost program_options and iostreams
RUN wget https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz && \
    tar -xf boost_1_86_0.tar.gz && \
    rm boost_1_86_0.tar.gz && \
    cd boost_1_86_0/ && \
    ./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=../boost && \
    ./b2 install && \
    cd .. && \
    cp boost/lib/libboost_iostreams.a boost/lib/libboost_program_options.a boost/lib/libboost_serialization.a /usr/local/lib/ && \
    cp -r boost/include/boost/ /usr/include/ && \
    rm -r boost_1_86_0 boost


RUN git clone https://github.com/odelaneau/GLIMPSE.git GLIMPSE && \
    cd GLIMPSE && \
    make clean && \
    make COMPILATION_ENV=docker && \
    mv chunk/bin/GLIMPSE2_chunk \
      split_reference/bin/GLIMPSE2_split_reference \
      phase/bin/GLIMPSE2_phase \
      ligate/bin/GLIMPSE2_ligate \
      concordance/bin/GLIMPSE2_concordance \
        /bin && \
    chmod +x /bin/GLIMPSE2* && \
    cd .. && \
    rm -rf GLIMPSE

WORKDIR /