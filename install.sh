#####################################################
# QMTIS QMTIMS QMTIS QMTIS QMTIS QMTIMS QMTIS QMTIS #
#####################################################
# Version 1.2
# Date: 2016-06-11
#####################################################
# Developed by: Pengfei Ji 
#####################################################

#!/bin/sh

# Make configuration for the QM, MD-TTM sub-programs
# QM
cd src/QM/
echo "The QM part will be installed in:" $(pwd)
./configure --prefix=$(pwd)
cd ../../
# MD-TTM
cd src/MD-TTM/src/
echo "The MD-TTM part will be installed in:" $(pwd)
./configure --prefix=$(pwd)
cd ../../../

# Executing the Makefiles for QM and MD-TTM
# QM
cd src/QM
make
cd ../../
# MD-TTM
cd src/MD-TTM/src/
make IMDSYS='Intel_CPU-gcc' imd-mpi-nvt-nve-ttm_Ce_Geph_Ke_Vf_AbI-eam-laser-stress
cd ../../../

# Create Symbol Links to the compiled files and BIN directory
cd src/QM/bin/
ln  ../../bin/* ./bin/*
cd ../../
cd src/MD-TTM/src/
ln  imd-mpi-nvt-nve-ttm_Ce_Geph_Ke_Vf_AbI-eam-laser-stress ../../../bin/imd-mpi-nvt-nve-ttm_Ce_Geph_Ke_Vf_AbI-eam-laser-stress
cd ../../../

#####################################################
