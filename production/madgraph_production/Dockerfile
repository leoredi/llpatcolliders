FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Basic tools + Python + Fortran/OpenMP
RUN apt-get update && apt-get install -y \
    wget curl git ca-certificates \
    python3 python3-pip \
    build-essential gfortran g++ make \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# --- Install MadGraph5_aMC@NLO ---
WORKDIR /opt

# Use the local MG5_aMC_v3.6.6.tar.gz that is in the build context
COPY MG5_aMC_v3.6.6.tar.gz /opt/

RUN tar -xzf /opt/MG5_aMC_v3.6.6.tar.gz && \
    rm /opt/MG5_aMC_v3.6.6.tar.gz

ENV MG5_DIR=/opt/MG5_aMC_v3_6_6
ENV PATH="$MG5_DIR/bin:${PATH}"
ENV MG5_PATH=/opt/MG5_aMC_v3_6_6/bin/mg5_aMC

# Python deps needed by MadGraph + your pipeline
RUN pip3 install --no-cache-dir six numpy pandas pylhe

# Install Pythia8 + LHAPDF through MG (first time is slow)
RUN echo "install pythia8\ninstall lhapdf6\nquit" > /tmp/install_cards.mg5 && \
    $MG5_DIR/bin/mg5_aMC /tmp/install_cards.mg5 && \
    rm /tmp/install_cards.mg5

# Mount point for your repo
WORKDIR /work
