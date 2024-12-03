FROM kylera/samtools-suite:gcloud as base

RUN apt-get update && \
    apt-get install -y \
      cmake libhts-dev libtabixpp-dev libtabixpp0 \
      pybind11-dev git python python-dev perl perl6

RUN cd /opt/ && \ 
    wget -q https://ziglang.org/download/0.10.1/zig-linux-x86_64-0.10.1.tar.xz && \
    tar xf zig-linux-x86_64*.tar.xz && \
    cd -

ENV PATH=$PATH:/opt/zig-linux-x86_64-0.10.1

RUN cd /opt/ && \
    git clone --recursive https://github.com/vcflib/vcflib.git && \
    cd vcflib && \
    mkdir -p build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DWFA_GITMODULE=ON .. && \
    cd .. && \
    cmake --build build/ --verbose -- -j 2 && \
    cd -

ENV PATH=$PATH:/opt/vcflib/build