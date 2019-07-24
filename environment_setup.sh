#!/bin/bash

export GENIESUPPORTPATH=/home/afroditi/lamp/GENIESupport
export PYTHIA6=$GENIESUPPORTPATH/pythia6/v6_424/lib
export LD_LIBRARY_PATH=$GENIESUPPORTPATH/pythia6/v6_424/lib:$LD_LIBRARY_PATH
export GSLLIB=$GENIESUPPORTPATH/gsl/lib
export GSLINC=$GENIESUPPORTPATH/gsl/include
export LD_LIBRARY_PATH=$GENIESUPPORTPATH/gsl/lib:$LD_LIBRARY_PATH
export ROOTSYS=$GENIESUPPORTPATH/root
export PATH=$GENIESUPPORTPATH/root/bin:$PATH
export LD_LIBRARY_PATH=$GENIESUPPORTPATH/root/lib:$LD_LIBRARY_PATH
export LOG4CPP_INC=$GENIESUPPORTPATH/log4cpp/include
export LOG4CPP_LIB=$GENIESUPPORTPATH/log4cpp/lib
export LD_LIBRARY_PATH=$GENIESUPPORTPATH/log4cpp/lib:$LD_LIBRARY_PATH
export LHAPATH=$GENIESUPPORTPATH/lhapdf
export LHAPDF_INC=$GENIESUPPORTPATH/lhapdf/include
export LHAPDF_LIB=$GENIESUPPORTPATH/lhapdf/lib
export LD_LIBRARY_PATH=$GENIESUPPORTPATH/lhapdf/lib:$LD_LIBRARY_PATH

<<<<<<< cc1ddeaa571fd03d79c8e48ae0cebd21bef9401a
export GENIE=/home/afroditi/devel_srcrecoil
=======
export GENIE=/home/afroditi/devel_fermimover
>>>>>>> reverting to previous version of FermiMover
export PATH=$GENIE/bin:$PATH
export LD_LIBRARY_PATH=$GENIE/lib:$LD_LIBRARY_PATH
export XSECSPLINEDIR=$GENIE/data

# for the reweighting
export GENIE_REWEIGHT=/home/afroditi/reweight
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GENIE_REWEIGHT/lib
export PATH=$PATH:$GENIE_REWEIGHT/bin
