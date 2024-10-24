# Use Ubuntu as a base image
#FROM ubuntu:20.04
FROM --platform=linux/amd64 ubuntu:20.04

# Set the environment variables to non-interactive to avoid interaction during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install basic system dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    libpthread-stubs0-dev \
    libomp-dev \
    autoconf \
    libtool \
    pkg-config \
    python3 \
    python3-pip \
    ninja-build \
    vim \
    && apt-get clean

# Manually install CMake version 3.26 or higher
# Download and install CMake
RUN wget https://cmake.org/files/v3.30/cmake-3.30.0-linux-x86_64.sh \
    && chmod +x cmake-3.30.0-linux-x86_64.sh \
    && ./cmake-3.30.0-linux-x86_64.sh --skip-license --prefix=/usr/local \
    && rm cmake-3.30.0-linux-x86_64.sh

# Set up a working directory
WORKDIR /usr/src/mycppproject

# Copy the entire project to the working directory
COPY . .

# Run CMake to configure the project, which will trigger ExternalProject_Add
RUN mkdir -p build && \
    cd build && \
    cmake .. && \
    cmake --build .

# Set the path to the executables in the bin directory
ENV PATH="/usr/src/mycppproject/build/bin:$PATH"

# Define the default command to run when the container starts
CMD ["/bin/bash"]
