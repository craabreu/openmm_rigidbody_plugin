# Installation prefix:
#PREFIX ?= /usr/local/openmm
PREFIX ?= /home/charlles/Software/anaconda3/pkgs/openmm-7.2.1-py36_1
NPROC ?= 4

all:
	mkdir -p build/
	cd build && cmake -DOPENMM_DIR=${PREFIX} ..
	cd build && make -j ${NPROC}

.PHONY: install clean PythonInstall

test:
	cd build && make -j ${NPROC} test

install:
	cd build && make -j ${NPROC} install

PythonInstall:
	cd build && make -j ${NPROC} PythonInstall

clean:
	rm -rf build/*
