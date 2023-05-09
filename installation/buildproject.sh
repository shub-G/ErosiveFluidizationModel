#!/bin/bash

ROOT=$(pwd)
BUILDDIR=release-build
OPTSFILE=$ROOT/source/dune/release.opts
PROJECT=$1

./source/dune/dune-common/bin/dunecontrol --builddir=$BUILDDIR  --opts=$OPTSFILE --only=$PROJECT all

