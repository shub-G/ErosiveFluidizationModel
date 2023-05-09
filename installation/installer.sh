#!/bin/bash
# usage: installer <directory>
#   directory: new directory in current directory where software is installed

# bash setup
set -x
set -e

# first we create a new directory, step into it and remember where we are
mkdir -p $1
cd $1
ROOT=$(pwd)

# make directory where to put external software
if [ ! "$EXTERNAL_HOME" ]; then
  EXTERNAL_HOME=$ROOT/external
  mkdir -p $EXTERNAL_HOME
fi

# extract default settings for the command names
if [ ! "$F77" ]; then
  F77=gfortran-9
fi
if [ ! "$CC" ]; then
CC=gcc-9
fi
if [ ! "$MPICC" ]; then
MPICC=mpicc
fi
if [ ! "$MPICXX" ]; then
MPICXX=mpicxx
fi
if [ ! "$CXX" ]; then
CXX=g++-9
fi
if [ ! "$CXXFLAGS" ]; then
CXXFLAGS="-O3 -DNDEBUG"
fi
CFLAGS="$CXXFLAGS"
if [ ! "$MAKE_FLAGS" ]; then
MAKE_FLAGS="-j"
fi

# check out all the required dune modules
checkout_dune_core ()
{
    # dune core modules
    for i in common geometry grid istl localfunctions grid-howto; do
	git clone https://gitlab.dune-project.org/core/dune-$i.git
	cd dune-$i
	git checkout master
	cd ..
    done

    # uggrid
    git clone https://gitlab.dune-project.org/staging/dune-uggrid.git
    cd dune-uggrid
    git checkout master
    cd ..

    # alugrid
    git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git
    cd dune-alugrid
    git checkout master
    cd ..

    # install whitespace hook
    ./dune-common/bin/dunecontrol vcsetup
}

checkout_dune_pdelab ()
{

    for i in typetree functions; do
    	git clone https://gitlab.dune-project.org/staging/dune-$i.git
    	cd dune-$i
    	git checkout master
    	cd ..
    done	

    for i in pdelab pdelab-tutorials; do
	git clone https://gitlab.dune-project.org/pdelab/dune-$i.git
	cd dune-$i
	git checkout master
	cd ..
    done

    # install whitespace hook
    ./dune-common/bin/dunecontrol vcsetup
}

# generate options file for Dune build system
generate_optsfile_release ()
{
    CC=$(which $CC)
    CXX=$(which $CXX)
    F77=$(which $F77)
echo "CMAKE_FLAGS=\"
-DCMAKE_C_COMPILER='$CC'
-DCMAKE_CXX_COMPILER='$CXX'
-DCMAKE_Fortran_COMPILER='$F77'
-DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -g0 -Wno-deprecated-declarations -funroll-loops'
-DCMAKE_BUILD_TYPE=Release
-DDUNE_SYMLINK_TO_SOURCE_TREE=1
\"" > release.opts
}

generate_optsfile_debug ()
{
    CC=$(which $CC)
    CXX=$(which $CXX)
    F77=$(which $F77)
echo "CMAKE_FLAGS=\"
-DCMAKE_C_COMPILER='$CC'
-DCMAKE_CXX_COMPILER='$CXX'
-DCMAKE_Fortran_COMPILER='$F77'
-DCMAKE_CXX_FLAGS_DEBUG='-g0 -Wno-deprecated-declarations'
-DCMAKE_BUILD_TYPE=Debug
-DDUNE_SYMLINK_TO_SOURCE_TREE=1
\"" > debug.opts
}

# create build script
generate_buildscript_release ()
{
echo '#!/bin/bash
ROOT=$(pwd)
BUILDDIR=release-build
OPTSFILE=release.opts
MODULES="common geometry uggrid grid alugrid istl localfunctions typetree functions pdelab grid-howto pdelab-tutorials"
for i in $MODULES; do
    echo build $i
    ./dune-common/bin/dunecontrol --builddir=$BUILDDIR  --opts=$OPTSFILE --only=dune-$i all
done
' > buildmodules_release.sh
chmod 755 buildmodules_release.sh
}

generate_buildscript_debug ()
{
echo '#!/bin/bash
ROOT=$(pwd)
BUILDDIR=debug-build
OPTSFILE=debug.opts
MODULES="common geometry uggrid grid alugrid istl localfunctions typetree functions pdelab grid-howto pdelab-tutorials"
for i in $MODULES; do
    echo build $i
    ./dune-common/bin/dunecontrol --builddir=$BUILDDIR  --opts=$OPTSFILE --only=dune-$i all
done
' > buildmodules_debug.sh
chmod 755 buildmodules_debug.sh
}

# create build script
generate_cleanscript_release ()
{
echo '#!/bin/bash
ROOT=$(pwd)
BUILDDIR=release-build
MODULES="common geometry uggrid grid alugrid istl localfunctions typetree functions pdelab grid-howto pdelab-tutorials"
for i in $MODULES; do
    echo clean $i
    cd dune-$i
    rm -rf $BUILDDIR
    cd ..
done
' > cleanmodules_release.sh
chmod 755 cleanmodules_release.sh
}

generate_cleanscript_debug ()
{
echo '#!/bin/bash
ROOT=$(pwd)
BUILDDIR=debug-build
MODULES="common geometry uggrid grid alugrid istl localfunctions typetree functions pdelab grid-howto pdelab-tutorials"
for i in $MODULES; do
    echo clean $i
    cd dune-$i
    rm -rf $BUILDDIR
    cd ..
done
' > cleanmodules_debug.sh
chmod 755 cleanmodules_debug.sh
}

# now run your selection of commands
checkout_dune_core
checkout_dune_pdelab
generate_optsfile_release
generate_buildscript_release
generate_cleanscript_release
generate_optsfile_debug
generate_buildscript_debug
generate_cleanscript_debug
