# Copyright 2020 Institut National de Recherche en Informatique et en Automatique (https://www.inria.fr/)
#                Poznan Supercomputing and Networking Center (https://www.psnc.pl/)

ARG release=22.04

FROM ubuntu:$release AS base

# OpenMPI workarounds
# Fix OMPI Vader: https://github.com/open-mpi/ompi/issues/4948
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
# Allow oversubscribing
ENV OMPI_MCA_rmaps_base_oversubscribe=1
# Disable busy waiting
ENV OMPI_MCA_mpi_yield_when_idle=1
# Make tzdata configure non-interactive
ENV TZ=Europe/Paris

ENV HOME=/home/docker

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
	echo $TZ >/etc/timezone && \
	apt-get update 

# Do not run testcases as root as mpirun does not like this...
ARG userid=1000
RUN adduser \
 	--uid=$userid \
 	--shell=/bin/bash \
 	--gecos='Docker,,,,' \
 	--disabled-password \
 	docker && \
	echo "docker:docker" | chpasswd

WORKDIR $HOME
COPY --chown=docker ./msolve-app ./msolve-app


FROM base AS installer
RUN apt-get --yes install \
	git \
	ca-certificates \
	cmake \
	build-essential \
	vim 

RUN apt-get --yes install \
# spack & melissa dependencies
	coreutils \
	curl \
	environment-modules \
	gfortran \
	gpg \
	lsb-release \
	python3-distutils \
	unzip \
	zip \
	libopenmpi-dev \
	python3 \
	python3-dev \
	python3-venv \
	libzmq3-dev
 

FROM installer as melissa-builder
# install spack and melissa
USER docker
RUN git clone -c feature.manyFiles=true https://github.com/spack/spack.git && \
	. ./spack/share/spack/setup-env.sh && \
	spack install melissa-api && \
	spack install py-melissa-core
 

FROM melissa-builder AS app
USER docker
WORKDIR $HOME/msolve-app
RUN	. $HOME/spack/share/spack/setup-env.sh && \
	spack load melissa-api && \
	spack load py-melissa-core && \
	cd $HOME/msolve-app/c-wrapper/ && \
	mkdir build && cd build && \
	cmake .. && make 

ENTRYPOINT ["/home/docker/msolve-app/setup-env.sh"]

# Command to keep the container running
#CMD ["tail", "-f", "/dev/null"]
