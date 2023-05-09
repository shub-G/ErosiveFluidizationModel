#!/bin/bash

ROOT=$(pwd)
BUILDDIR=release-build
OPTSFILE=$ROOT/source/dune/release.opts
PROJECT=$1

./source/dune/dune-common/bin/duneproject $PROJECT

