FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /tmp

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    python3 \
    python3-pip \
    python3-dev \
    gfortran \
    libxml2-dev \
    libssl-dev \
    libpcre3-dev \
    xlibmesa-glu-dev \
    libglew-dev \
    libftgl-dev \
    libfftw3-dev \
    libcfitsio-dev \
    libgraphviz-dev \
    libavahi-compat-libdnssd-dev \
    libldap2-dev \
    libkrb5-dev \
    libgsl-dev \
    libqt5webkit5-dev \
    qtwebengine5-dev \
    libxpm-dev \
    libxft-dev \
    libxext-dev \
    libx11-dev \
    libtbb-dev \
    libvdt-dev \
    libgif-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff-dev \
    vim \
    nano \
    less \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir \
    numpy \
    pandas \
    matplotlib \
    scipy \
    trimesh \
    shapely \
    tqdm \
    jupyter \
    notebook

# Build and install Pythia 8.312
ENV PYTHIA8=/opt/pythia8
ENV PYTHIA8DATA=/opt/pythia8/share/Pythia8/xmldoc

RUN mkdir -p /opt && \
    cd /opt && \
    wget https://pythia.org/download/pythia83/pythia8312.tgz && \
    tar xzf pythia8312.tgz && \
    mv pythia8312 pythia8 && \
    rm pythia8312.tgz && \
    cd pythia8 && \
    ./configure --prefix=/opt/pythia8 --enable-shared && \
    make -j$(nproc) && \
    make install

# Clone and build ROOT with Pythia8 support
ENV ROOTSYS=/opt/root
ENV PATH="${ROOTSYS}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${ROOTSYS}/lib:${PYTHIA8}/lib:${LD_LIBRARY_PATH}"
ENV PYTHONPATH="${ROOTSYS}/lib:${PYTHONPATH}"

RUN cd /opt && \
    git clone --branch latest-stable --depth 1 https://github.com/root-project/root.git root-source && \
    mkdir root-build && \
    cd root-build && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/root \
          -Dpythia8=ON \
          -DPYTHIA8_DIR=/opt/pythia8 \
          -Dpython3=ON \
          -Dbuiltin_glew=ON \
          -Dbuiltin_ftgl=ON \
          ../root-source && \
    cmake --build . --target install -j$(nproc) && \
    cd /opt && \
    rm -rf root-source root-build

# Set up environment to automatically source ROOT
RUN echo "source /opt/root/bin/thisroot.sh" >> /etc/bash.bashrc && \
    echo "export PYTHIA8=/opt/pythia8" >> /etc/bash.bashrc && \
    echo "export PYTHIA8DATA=/opt/pythia8/share/Pythia8/xmldoc" >> /etc/bash.bashrc

# Set up a non-root user (devcontainer best practice)
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Switch to non-root user
USER $USERNAME
WORKDIR /workspace

# Ensure ROOT is sourced for the user as well
RUN echo "source /opt/root/bin/thisroot.sh" >> ~/.bashrc && \
    echo "export PYTHIA8=/opt/pythia8" >> ~/.bashrc && \
    echo "export PYTHIA8DATA=/opt/pythia8/share/Pythia8/xmldoc" >> ~/.bashrc

CMD ["/bin/bash"]
