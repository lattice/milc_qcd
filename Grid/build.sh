#! /bin/bash

ARCH=$1
PK_CC=$2
PK_CXX=$3
GIT_BRANCH=feature/staggered-comms-compute

if [ -z ${PK_CXX} ]
then
  echo "Usage $0 <scalar|avx512|avx2> <PK_CC> <PK_CXX>"
  exit 1
fi

case ${ARCH} in
    scalar|avx512|avx2)
      ;;
    *)
      echo "Unsupported ARCH"
      echo "Usage $0 <scalar|avx512|avx2> <PK_CC> <PK_CXX>"
      exit 1
esac

MAKE="make -j4"

TOPDIR=`pwd`
SRCDIR=${TOPDIR}/Grid
BUILDDIR=${TOPDIR}/build-${ARCH}
INSTALLDIR=${TOPDIR}/install-${ARCH}

MAKE="make -j4"

if [ ! -d ${SRCDIR} ]
then
  echo "Fetching ${GIT_BRANCH} branch of Grid package from github"
  git clone https://github.com/paboyle/Grid -b ${GIT_BRANCH}
fi

# Fetch Eigen package, set up Make.inc files and create Grid configure
pushd ${SRCDIR}
./bootstrap.sh
popd

# Configure only if not already configured                                          
mkdir -p ${BUILDDIR}
pushd ${BUILDDIR}
if [ ! -f Makefile ]
then
  echo "Configuring Grid for ${ARCH} in ${BUILDDIR}"

  case ${ARCH} in

    scalar)

       module swap craype-mic-knl craype-haswell 

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --enable-precision=double \
            --enable-simd=GEN \
            --enable-comms=none \
	    --with-lime=${HOME}/scidac/install/qio-cori-omp-knl-icc \
            --with-openssl=/global/common/cori/software/openssl/1.1.0a/hsw \
            CXX="${PK_CXX}" \
            CXXFLAGS="-std=c++11" \

# 	    --with-hdf5=/opt/cray/pe/hdf5/1.10.0/INTEL/15.0 \

       status=$?
             ;;

    avx2)

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --enable-precision=double \
            --enable-simd=GEN \
            --enable-comms=mpi \
	    --with-lime=${HOME}/scidac/install/qio-cori-omp-knl-icc \
            --with-openssl=/global/common/cori/software/openssl/1.1.0a/hsw \
            CXX="${PK_CXX}" CC="${PK_CC}" \
            CXXFLAGS="-std=c++11 -xCORE-AVX2" \

#	    --with-hdf5=/opt/cray/pe/hdf5/1.10.0/INTEL/15.0 \
#            --enable-mkl=yes \

       status=$?
             ;;
    avx512)

       INCMKL="-I/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include"
       LIBMKL="-L/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin"

       ${SRCDIR}/configure \
            --prefix=${INSTALLDIR} \
            --enable-precision=double \
            --enable-simd=KNL \
            --enable-comms=mpi \
            --host=x86_64-unknown-linux-gnu \
	    --with-lime=${HOME}/scidac/install/qio-cori-omp-knl-icc \
            --with-openssl=/global/common/cori/software/openssl/1.1.0a/hsw \
            CXX="${PK_CXX}" CC="${PK_CC}" \
            CXXFLAGS="-std=c++11 -xMIC-AVX512" \

	    # --with-hdf5=/opt/cray/pe/hdf5/1.10.0.3/INTEL/16.0 \

       status=$?
       echo "Configure exit status $status"
             ;;
    *)
    echo "Unsupported ARCH ${ARCH}"
          exit 1;
  esac

  if [ $status -ne 0 ]
  then
      echo "Quitting because of configure errors"
  else
    echo "Building in ${BUILDDIR}"
    ${MAKE} -j4

    echo "Installing in ${INSTALLDIR}"
    ${MAKE} install
  fi

fi     
popd




