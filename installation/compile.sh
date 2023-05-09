#!/bin/bash

ROOT=$(pwd)
BUILDDIR=release-build
PROJECT=$1

cd $PROJECT/$BUILDDIR/src
make 
cd ../../..
