/*
 * main.cc
 *
 *  Created on: Apr 6, 2021
 *      Author: sgupta
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
#include<stdlib.h>
#include<time.h>
#include<exception>
#include<chrono>
#include<array>

/*********************************************************************/
#define DIMENSION 2

//#define TWOPHASEFLOW
#define DECOUPLED

//#define USE_UG
/*********************************************************************/
#include"../duneincludes.hh"
#if defined(TWOPHASEFLOW)
#include"parameters/test_2pflow/parameters.hh"
#include"operators/lop_2pflow_debug.hh"
#include"driver_2pflow.hh"
//#elif defined(DECOUPLED)
#else
struct PrimaryVariables
{
  enum Type{ Pw=0, Sw=1, porosity=0, Cf=1, Pw_l2=0, Pn_l2=1 };
};
#include"parameters/pockmarks_study/parameters.hh"
#include"operators/lop_2pflow.hh"
#include"operators/lop_l2projection.hh"
#include"operators/lop_erosion.hh"
#include"driver_decoupled.hh"
//#else
//{
//  enum Type{ Pw=0, Sw=1, porosity=2, Cf=3, Pw_q1=4, Pn_q1=5 };
//};
//#include"parameters/pockmarks_study/parameters.hh"
//#include"operators/lop_coupled.hh"
//#include"driver_coupled.hh"
//struct PrimaryVariables
#endif
/*********************************************************************/

int main(int argc, char** argv)
{
	try{
	    // Maybe initialize MPI
	    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
	    if(helper.rank()==0){
		#ifdef TWOPHASEFLOW
	    std::cout << "Hello World! This is PROBLEM09 of project LandscapeEvolution." << std::endl;
	    std::cout << "This problem considers fully coupled fully implicit FV implementation "
	    			 "of unsaturated fluid flow unsing DELAUNEY TRIANGULATED MESH"
	    		  << std::endl;
		#else
	    std::cout << "Hello World! This is PROBLEM09 of project LandscapeEvolution." << std::endl;
	    std::cout << "This problem considers decoupled FV implementation "
	    			 "of unsaturated fluid flow and soil erosion."
	    		  << std::endl;
		#endif
	    }
	    if(Dune::MPIHelper::isFake){
	      std::cout<< "This is a sequential program." << std::endl;
	    }
	    else {
	    	if(helper.rank()==0){
	    		std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()<<" processes!"<<std::endl;
	    	}
	    }

		/**************************************************************************************************/
		// INPUTS
	    if (argc!=2)
	    {
	    	if(helper.rank()==0){
	    		std::cout << "usage: ./problem09 <input-file> " << std::endl;
	    		std::cout << "hint: sample (default) ..." << std::endl;
	    	}
	        return 1;
	    }

        /**************************************************************************************************/
	    // DUNE MODEL PATH
	    const std::string MODEL_PATH = "/home/sgupta/dune_2_8/LandscapeEvolution/src/";
	    // PROBLEM NAME
	    const std::string PROBLEM_NAME = "problem09/";
	    // INPUT PATH NAME
	    std::string INPUT_PATH  = MODEL_PATH + PROBLEM_NAME + "inputs/";
	    // OUTPUT PATH NAME
	    std::string OUTPUT_PATH = MODEL_PATH + "outputs/" + PROBLEM_NAME;
		#ifdef TWOPHASEFLOW
		INPUT_PATH  += "test_2pflow/";
		OUTPUT_PATH += "test_2pflow/";
		#else
		INPUT_PATH  += "pockmarks_study/";
		OUTPUT_PATH += "pockmarks_study/";
		#endif
	    // INI-FILE FOR USER-DEFINED INPUTS
	    char INI_FILE[100];
	    sscanf(argv[1],"%99s", INI_FILE);
	    std::string input_file = INPUT_PATH;
	    input_file += INI_FILE;
	    input_file += ".ini";
        if(helper.rank()==0){
	    std::cout<< "input file: " << input_file << std::endl ;
        }

        /**************************************************************************************************/
        // PARAMETER TREE
	    Dune::ParameterTree ptree;
	    Dune::ParameterTreeParser ptreeparser;
	    ptreeparser.readINITree(input_file,ptree);
	    ptreeparser.readOptions(argc,argv,ptree);

		/**************************************************************************************************/
		// MESH
#if DIMENSION==2
		const int dim = 2;
#elif DIMENSION==3
		const int dim = 3;
#else
		std::cout<< "Invalid DIMENSION:" << DIMENSION
				 << " in LINE:" << __LINE__
				 << " of FILE:" << __FILE__
				 << std::endl;
		exit(0);
#endif

		/**************************************************************************************************
		 * READ MESH
		 * AND GENERATE GRID VIEW
		 * WITH UG MESH
		 **************************************************************************************************/
#ifdef USE_UG
		using GridType = Dune::UGGrid<dim>;
		GridType grid_type;
		#if DIMENSION==2
		const std::string grid_name = ptree.get("domain.ug.name",(std::string)"default_1x1");
		#elif DIMENSION==3
		const std::string grid_name = ptree.get("domain.ug.name",(std::string)"default_1x1x1");
#else
		std::cout<< "Invalid DIMENSION:" << DIMENSION
				 << " in LINE:" << __LINE__
				 << " of FILE:" << __FILE__
				 << std::endl;
		exit(0);
#endif
		auto grid_file = MODEL_PATH + PROBLEM_NAME + "mesh/";
		#ifdef TWOPHASEFLOW
		grid_file  += "test_2pflow/";
		#else
		grid_file  += "pockmarks_study/";
		#endif
		grid_file += grid_name;
		grid_file += ".msh";
		Dune::GmshReader<GridType> gmshreader;
		std::shared_ptr<GridType> grid(gmshreader.read(grid_file,true,false));

		using GV = GridType::LeafGridView;
		GV gv = grid->leafGridView();
        grid->loadBalance();

#else
		Dune::FieldVector<double,dim> L(1.0);	// L represents the right top node of the rectangular/cuboidal domain
		std::array<int,dim> N;
		double xc = ptree.get("characteristic_values.length",(double)1.0);//in m
		double X = ptree.get("domain.length",(double)1.0);//in m
		X *= 1./xc; //ndim
		double Z = ptree.get("domain.height",(double)1.0);//in m
		Z *= 1./xc; //ndim
		const int X_cells = ptree.get("domain.yasp.nX",(int)100.0) ;
		const int Z_cells = ptree.get("domain.yasp.nZ",(int)100.0) ;
		L[0] = X;
		L[dim-1] = Z;
		N[0] = X_cells ;
		N[dim-1] = Z_cells ;
		#if DIMENSION==3
		L[1] = X;
		N[1] = X_cells ;
		#endif
		std::bitset<dim> periodic(false);
		int overlap=1;
		Dune::YaspGrid<dim> grid(L,N,periodic,overlap,helper.getCommunicator());
		using GV = Dune::YaspGrid<dim>::LeafGridView;
		const GV& gv=grid.leafGridView();
        grid.loadBalance();
#endif

        /**************************************************************************************************/
    	// MATERIAL PROPERTIES, NUMERICAL AND TEST PARAMETERS, CONSTITUTIVE RELATIONSHIPS
		#ifdef TWOPHASEFLOW
    	using Parameters = TwoPhaseFlowTestParameters<GV,Dune::ParameterTree>;
    	Parameters parameter(gv,ptree);
		#else
    	using Parameters = PockmarksStudyParameters<GV,Dune::ParameterTree>;
    	Parameters parameter(gv,ptree);
		#endif

		/**************************************************************************************************/
		// DRIVER
		driver( gv,
				ptree,
				parameter,
        		OUTPUT_PATH,
				helper);
		/**************************************************************************************************/
	}
	catch (Dune::Exception &e){
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...){
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
