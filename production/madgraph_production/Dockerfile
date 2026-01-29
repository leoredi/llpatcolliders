FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Basic tools + Python + Fortran/OpenMP
RUN apt-get update && apt-get install -y \
    wget curl git ca-certificates rsync \
    python3 python3-pip python3-dev \
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

# Python deps needed by MadGraph + decay pipeline
RUN pip3 install --no-cache-dir \
    six numpy pandas pylhe \
    sympy mpmath numba particle scipy 'scikit-hep==0.4.0'

# Install LHAPDF through MG (Pythia8 is built separately below)
RUN echo "install lhapdf6\nquit" > /tmp/install_cards.mg5 && \
    $MG5_DIR/bin/mg5_aMC /tmp/install_cards.mg5 && \
    rm /tmp/install_cards.mg5

# Build Pythia8 with Python bindings (standalone, avoids MG5 hepmc install issues on ARM)
RUN ln -sf /usr/bin/python3-config /usr/bin/python-config && \
    wget -q https://pythia.org/download/pythia83/pythia8315.tgz -O /tmp/pythia8315.tgz && \
    tar -xzf /tmp/pythia8315.tgz -C /opt && \
    rm /tmp/pythia8315.tgz && \
    cd /opt/pythia8315 && \
    ./configure --prefix=/opt/pythia8 --with-python && \
    make -j2 && \
    make install

ENV PYTHIA8DATA="/opt/pythia8/share/Pythia8/xmldoc"
ENV LD_LIBRARY_PATH="/opt/pythia8/lib"
ENV PYTHONPATH="/opt/pythia8/lib/python3.10/site-packages"

# Copy the HNL model (SM_HeavyN_CKM_AllMasses_LO) into MadGraph models directory
COPY mg5/models/SM_HeavyN_CKM_AllMasses_LO $MG5_DIR/models/SM_HeavyN_CKM_AllMasses_LO

# Mount point for your repo
WORKDIR /work
