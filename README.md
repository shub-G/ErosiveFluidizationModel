# ErosiveFluidizationModel
This repository contains source code for the model for simulation of pipes, chimneys and pockmarks in subsurface due to flow of over-pressured gas. The model is based on the mechanism of 'erosive fluidization'. This code was used in the simulations included in the manuscript titled "Mathematical modelling of the influence of erosive fluidization on the
morphology of fluid flow and escape structures", published in Mathematical Geosciences.

## Instructions for installation
* System requirements: 
  * Ububtu LTS 18.04 or higher
  * packages: git, auto-conf, cmake (>=3.13), clang, g++, gcc, gfortran, superlu (dev), lapack, libblas, libboost, metis (dev), parmetis (dev), gmsh, paraview, openmpi 
* Execute following commands in terminal (in given order):
  * mkdir dune_2_8
  * cd mkdir dune_2_8
  * mkdir source
  * cd mkdir source
  * cp /downloads/installer.sh .
  * chmod 755 installer.sh
  * ./installer.sh dune
  * cd dune
  * ./buildmobules.sh
  * cd ../..
  * chmod 755 newproject.sh
  * ./newproject.sh ErosiveFluidizationModel
    * On prompt for "2) Which modules should this module depend on?", enter: dune-common dune-geometry dune-uggrid dune-grid dune-localfunctions dune-istl dune-typetree dune-functions dune-alugrid dune-pdelab
    * On prompt for "3) Project/Module version?", enter 1
    * On promt for "4) Maintainer's email address?", enter your email address.
  * chmod 755 buildproject.sh
  * ./buildproject.sh ErosiveFluidizationModel
  * cd ErosiveFluidizationModel/src/
  * rm -rf ErosiveFluidizationModel.cc
  * rm -rf CMakeLists.txt
  * cp \_all_source_files_in_repo\_ .
  * cd ../..
  * chmod 755 compile.sh
  * ./compile.sh ErosiveFluidizationModel

## To run the simulations:
* Execute following commands in terminal (in given order):
  * cd \_HOME\_/dune_2_8/ErosiveFluidizationModel/release-build/src
  * ./main \_input-file\_  
    * In this case: sample

## Files included in this repo:
* installation files:
  * installation/installer.sh
  * installation/newproject.sh
  * installation/buildproject.sh
  * installation/compile.sh
* source files 
  * duneincludes.hh
  * problem09/main.cc
    * Implement the correct MODEL_PATH, INPUT_PATH, and OUTPUT_PATH before compiling.
  * problem09/driver_coupled.hh
  * problem09/driver_decoupled.hh
  * problem09/driver_2pflow.hh
  * parameters/pockmarks_study/parameter.hh
  * inputs/pockmarks_study/sample.ini
  * inputs/pockmarks_study/sample_input_file_generator.sh
  * mesh/E1_mesh.msh
  * operators/lop_2pflow.hh
  * operators/lop_coupled.hh
  * operators/lop_erosion.hh
  * operators/lop_l2projection.hh
* outputs
