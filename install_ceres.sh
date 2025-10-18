#!/bin/bash

# simple script to install ceres on linux
# ripped from here: http://ceres-solver.org/installation.html#linux

echo "getting root..."
sudo sleep 1

# CMake
sudo apt-get install -y cmake
# google-glog + gflags
sudo apt-get install -y libgoogle-glog-dev libgflags-dev
# Use ATLAS for BLAS & LAPACK
sudo apt-get install -y libatlas-base-dev
# Eigen3
sudo apt-get install -y libeigen3-dev
# SuiteSparse (optional)
sudo apt-get install -y libsuitesparse-dev

wget http://ceres-solver.org/ceres-solver-2.2.0.tar.gz
tar zxf ceres-solver-2.2.0.tar.gz
rm -r -f ceres-bin
mkdir ceres-bin
cd ceres-bin
cmake ../ceres-solver-2.2.0
make -j3
make test
sudo make install
cd ..
rm -r -f ceres-solver-2.2.0.tar.gz
