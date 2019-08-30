FROM ubuntu:16.04
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn

ADD * /simulator/
RUN chmod -R 777 /simulator/
WORKDIR /simulator/

USER root

RUN apt-get update && apt-get install -yq --no-install-recommends \
    apt-utils \
    build-essential \
    gcc-multilib \
    gfortran \
    libgoogle-perftools-dev \
    libblas-dev \
    libboost-dev \
    libgsl-dev \
    git \
    gsl-bin \
    nano \
    python \
    python-pip \
    parallel \
    libreadline-dev \
    libssh2-1-dev \
    libssl-dev


RUN git config --global http.sslVerify false
RUN git clone https://github.com/corbett-lab/tRNA.git && cd tRNA && make
