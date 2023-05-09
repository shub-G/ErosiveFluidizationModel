/*
 * driver_coupled.hh
 *
 *  Created on: Apr 7, 2021
 *      Author: sgupta
 */

#ifndef PROBLEM09_DRIVER_COUPLED_HH_
#define PROBLEM09_DRIVER_COUPLED_HH_

template<typename GV,typename PTree,typename Parameters>
void driver( const GV& gv, 				// GridView
		 	 const PTree& ptree, 		// Input Parameters-Tree
			 Parameters& parameter,
			 std::string output_path,
			 Dune::MPIHelper& helper){

	std::string pathName = output_path;
	auto pathExt = ptree.get("output.path_name",(std::string)"test0");
	pathName += pathExt;
	pathName += "/";
	auto fileName = ptree.get("output.file_name",(std::string)"test_0");
	fileName += "_FC";

	using Coord = typename GV::Grid::ctype;
	const int dim = GV::dimension;

	double t_ast  = ptree.get("characteristic_value.time"	,(double)1.0);//s

	double time = 0.0;//ndim
	double dt = ptree.get("time.dt_initial",(double)0.001);//read in s
	dt *= 1./t_ast; //ndim
	// total simulation-time
	double t_END  = ptree.get("time.time_end",(double)1000.);//read in s
	t_END *= 1./t_ast; //ndim
	// output time interval
	double t_OP   = ptree.get("output.time_interval",(double)10.);//read in s
	t_OP *= 1./t_ast; //ndim
	//adaptive time control
	bool is_adaptive_time_control = ptree.get("adaptive_time_control.flag",(bool)true);
	double dt_min = ptree.get("adaptive_time_control.dt_min",(double)1.e-6);//read in s
	dt_min *= 1./t_ast; //ndim
	double dt_max = ptree.get("adaptive_time_control.dt_max",(double)1.);//read in s
	dt_max *= 1./t_ast; //ndim
	int maxAllowableIterations = ptree.get("adaptive_time_control.max_newton_steps",(int)6);
	int minAllowableIterations = ptree.get("adaptive_time_control.min_newton_steps",(int)4);

	double dtstart = dt; //ndim
	double time_op = time; //ndim
	double clock_time_elapsed = 0.; //s

	/****************************************************/
	// GRID FUNCTION SPACES
	/****************************************************/
	using CON0 = Dune::PDELab::P0ParallelConstraints;
	using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;

	// FEM SPACE P0-fem
//#ifdef USE_UG
//	auto gt = Dune::GeometryTypes::simplex(dim);
//#else
	auto gt = Dune::GeometryTypes::cube(dim);
//#endif
	typedef Dune::PDELab::P0LocalFiniteElementMap<Coord,double,dim> FEMP0;
	FEMP0 femp0(gt);
	using GFSP0 = Dune::PDELab::GridFunctionSpace<GV, FEMP0, CON0, VBE0>;
	GFSP0 gfsp0(gv, femp0);
	GFSP0 gfs_pw( gv, femp0);
	GFSP0 gfs_sw( gv, femp0);
	GFSP0 gfs_por(gv, femp0);
	GFSP0 gfs_cf( gv, femp0);

	const int degree = 1; //polynomial degree
	typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,double,degree> FEMQ1;
	FEMQ1 femq1(gv);
	typedef Dune::PDELab::GridFunctionSpace<GV,FEMQ1,CON0,VBE0> GFSQ1;
	GFSQ1 gfs_pwq1(gv,femq1);
	GFSQ1 gfs_swq1(gv,femq1);

	using VBE =Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
	using OrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
	using GFS = Dune::PDELab::CompositeGridFunctionSpace< VBE, OrderingTag,
														  GFSP0, GFSP0, GFSP0, GFSP0,
														  GFSQ1, GFSQ1 > ;
	GFS gfs(gfs_pw, gfs_sw, gfs_por, gfs_cf, gfs_pwq1, gfs_swq1);
	typedef typename GFS::template ConstraintsContainer<double>::Type C;
	C c;
	c.clear();
	gfs.update();
	std::cout << "degrees of freedom (two phase flow model): " << gfs.globalSize() << std::endl;

	/****************************************************/
	// VECTOR CONTAINERS
	/****************************************************/
	using U = Dune::PDELab::Backend::Vector<GFS,double>;
	U u_old(gfs);
	U u_new(gfs);

	/****************************************************/
	// INITIAL CONDITIONS
	/****************************************************/
	using Initial = InitialConditions<GV,PTree,Parameters>;
	Initial initial(gv,ptree,parameter);
	auto pw_ic_local = [&](const auto& e, const auto& x){return initial.Pw(e,x);};
	auto pw_ic = Dune::PDELab::makeGridFunctionFromCallable(gv,pw_ic_local);
	auto sw_ic_local = [&](const auto& e, const auto& x){return initial.Sw(e,x);};
	auto sw_ic = Dune::PDELab::makeGridFunctionFromCallable(gv,sw_ic_local);
	auto por_ic_local= [&](const auto& e, const auto& x){return initial.porosity(e,x);};
	auto por_ic = Dune::PDELab::makeGridFunctionFromCallable(gv,por_ic_local);
	auto cf_ic_local = [&](const auto& e, const auto& x){return initial.Cf(e,x);};
	auto cf_ic = Dune::PDELab::makeGridFunctionFromCallable(gv,cf_ic_local);
	using  InitialConditions = Dune::PDELab::CompositeGridFunction< decltype(pw_ic),
																	decltype(sw_ic),
																	decltype(por_ic),
																	decltype(cf_ic),
																	decltype(pw_ic),
																	decltype(sw_ic)>;
	InitialConditions ic( pw_ic,sw_ic,por_ic,cf_ic,pw_ic,sw_ic );

	// 	Initialize the solution at t=0 (uold) with the given initial values
	Dune::PDELab::interpolate( ic, gfs, u_old );
	u_new = u_old;

	/****************************************************/
	// BOUNDARY CONDITIONS
	/****************************************************/
	using WBC	= WettingPhaseBoundaryCondition< GV,PTree,Parameters>;
	WBC  w_bc( gv,ptree,parameter);
	using NBC	= NonwettingPhaseBoundaryCondition< GV,PTree,Parameters>;
	NBC  n_bc( gv,ptree,parameter);
	using FSBC	= FluidizedSoilBoundaryCondition< GV,PTree,Parameters>;
	FSBC fs_bc(gv,ptree,parameter);

	/****************************************************/
	// LOCAL OPERATORS
	/****************************************************/
	unsigned int intorder = ptree.get("CG_parameters.intorder",(int)4);
	using LOP = LocalOperatorCoupled<GV,Parameters,GFS,U,WBC,NBC,FSBC>;
	LOP lop( gv, parameter,
			 gfs, &u_old,
			 w_bc, n_bc, fs_bc,
			 &time, &dt,
			 intorder );

	/****************************************************/
	// GRID OPERATOR
	/****************************************************/
	using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
	MBE mbe(20); // Maximal number of nonzeroes per row

	using GOLOP = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,double,double,double,C,C>;
	GOLOP go_lop(gfs,c,gfs,c,lop,mbe);
	typename GOLOP::Traits::Jacobian jac(go_lop);
	if(helper.rank()==0){
		std::cout << jac.patternStatistics() << std::endl;
	}

	/****************************************************/
	// LINEAR SOLVER
	/****************************************************/
	using LS = Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GOLOP>;
	LS ls(gfs,1000,1);
	Dune::Amg::Parameters ls_params = ls.parameters();
	ls_params.setCoarsenTarget(100000);// max DoF at coarsest level
	ls.setParameters(ls_params);

	/****************************************************/
	// NON-LINEAR SOLVER
	/****************************************************/
	using PDESOLVER = Dune::PDELab::Newton<GOLOP,LS,U>;
	PDESOLVER pdesolver(go_lop,ls);
	pdesolver.setLineSearchStrategy(PDESOLVER::Strategy::noLineSearch);
	pdesolver.setReassembleThreshold(0.0);
	pdesolver.setVerbosityLevel(ptree.get("newton.verbosity",(int)2));
	pdesolver.setReduction(ptree.get("newton.reduction",(double)1e-6));
	pdesolver.setMinLinearReduction(ptree.get("newton.min_lin_reduction",(double)1e-9));
	pdesolver.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
	pdesolver.setForceIteration(ptree.get("newton.force_iteration",(bool)false));
	pdesolver.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-6));

	/****************************************************/
	// VTK
	/****************************************************/
	int subsampling = 1;
	using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
	VTKWRITER vtkwriter(gv,Dune::refinementIntervals(subsampling));
	using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
	VTKSEQUENCEWRITER vtkSequenceWriter( std::make_shared<VTKWRITER>(vtkwriter),fileName,pathName,"");

	using SUBGFS_Pw = typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Pw>> >;
	SUBGFS_Pw subgfs_pw(gfs);
	using DGF_Pw = Dune::PDELab::DiscreteGridFunction< SUBGFS_Pw, U >;
	DGF_Pw dgf_pw( subgfs_pw, u_new );
	vtkSequenceWriter.addCellData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_Pw> >(dgf_pw,"pw"));

	using SUBGFS_Sw = typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Sw>> >;
	SUBGFS_Sw subgfs_sw(gfs);
	using DGF_Sw = Dune::PDELab::DiscreteGridFunction< SUBGFS_Sw, U >;
	DGF_Sw dgf_sw( subgfs_sw, u_new );
	vtkSequenceWriter.addCellData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_Sw> >(dgf_sw,"sw"));

	using SUBGFS_Pw_q1 = typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Pw_q1>> >;
	SUBGFS_Pw_q1 subgfs_pw_q1(gfs);
	using DGF_Pw_q1 = Dune::PDELab::DiscreteGridFunction< SUBGFS_Pw_q1, U >;
	DGF_Pw_q1 dgf_pw_q1( subgfs_pw_q1, u_new );
	vtkSequenceWriter.addVertexData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_Pw_q1> >(dgf_pw_q1,"pw_q1"));

	using SUBGFS_Pn_q1 = typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Pn_q1>> >;
	SUBGFS_Pn_q1 subgfs_pn_q1(gfs);
	using DGF_Pn_q1 = Dune::PDELab::DiscreteGridFunction< SUBGFS_Pn_q1, U >;
	DGF_Pn_q1 dgf_pn_q1( subgfs_pn_q1, u_new );
	vtkSequenceWriter.addVertexData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_Pn_q1> >(dgf_pn_q1,"pn_q1"));

	using SUBGFS_por = typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::porosity>> >;
	SUBGFS_por subgfs_por(gfs);
	using DGF_por = Dune::PDELab::DiscreteGridFunction< SUBGFS_por, U >;
	DGF_por dgf_por( subgfs_por, u_new );
	vtkSequenceWriter.addCellData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_por> >(dgf_por,"porosity"));

	using SUBGFS_Cf = typename Dune::PDELab::GridFunctionSubSpace< GFS, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Cf>> >;
	SUBGFS_Cf subgfs_cf(gfs);
	using DGF_Cf = Dune::PDELab::DiscreteGridFunction< SUBGFS_Cf, U >;
	DGF_Cf dgf_cf( subgfs_cf, u_new );
	vtkSequenceWriter.addCellData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_Cf> >(dgf_cf,"Cf"));

	vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

	/****************************************************/
	// TIME-LOOP
	/****************************************************/
	int opcount = 1;
	double timecount = time;
	double dtLast = dtstart;
	int dtFlag = 0;
	bool exceptionCaught = false;
	int newton_iterations = 0;

	while( time < t_END - 1e-8 ){

		if( exceptionCaught==false ){
			dt = std::max(dt,dt_min);
		}

		if(helper.rank()==0){
		std::cout<< "_____________________________________________________" <<std::endl;
		std::cout<< " current opcount = " << opcount - 1 << std::endl;
		}

		clock_t start = clock();
		try{
			pdesolver.apply( u_new );

			newton_iterations = pdesolver.result().iterations;

			exceptionCaught = false;

		}catch ( Dune::Exception &e ) {
			exceptionCaught = true;
			if( dt > 1e-15 ){

				if(helper.rank()==0){
				std::cout << "Catched Error, Dune reported error: " << e << std::endl;
				}
				u_new = u_old;
				newton_iterations = 0;
				dt *= 0.5;

				continue;
			}
			else
			{
				if(helper.rank()==0){
					std::cout << "ABORTING, due to DUNE error: " << e << std::endl;
				}
				exit(0);
			}
		}
		clock_t end = clock();
		double clock_time_this_step = (double) (end-start) / CLOCKS_PER_SEC;
		clock_time_elapsed += clock_time_this_step;

		if(helper.rank()==0){
		std::cout<<"DONE"<<std::endl;
		std::cout<<"_____________________________________________________"<<std::endl;
		}

		/***********************************************
		 * OUTPUT
		 ***********************************************/
		/* At t_OP
		 *
		 */
		if( ( time+dt > t_OP*opcount - dt_min ) and ( time+dt < t_OP*opcount + dt_min )  )
		{
			// WRITE OUTPUT
			vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

			if(helper.rank()==0){
			std::cout<< " ******************************************************************* " << std::endl;
			std::cout<< " OUTPUT WRITTEN " << opcount << std::endl;
			std::cout<< " ******************************************************************* " << std::endl;
			std::cout<< std::flush;
			}

			timecount = time;
			opcount = opcount+1;
		}

		/***********************************************/
		// PREPARE FOR NEXT TIME INTEGRATION
		/***********************************************/
		//1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
		u_old = u_new;
		//2. ADVANCE TIME:
		time += dt;

		if(helper.rank()==0){
		std::cout<<" "<< std::endl;
		std::cout<< " time = " << time ;
		std::cout<< std::flush;
		}

		if( is_adaptive_time_control ){
			if(newton_iterations>maxAllowableIterations){
				dt=std::max(dt*0.75,dt_min);
			}
			else if(newton_iterations<=minAllowableIterations){
				dt=std::min(dt*1.25,dt_max);
			}
		}
		else{
			dt = dtstart;
		}

		if(helper.rank()==0){
		std::cout << " , time+dt = " << (time + dt)*t_ast
				  << " , opTime = "  << t_OP*t_ast * opcount ;
		std::cout<< std::flush;
		}

		if( time + dt  > t_OP * opcount){
			dtLast = dt;
			dt 	 = t_OP * opcount - time ;

			if(helper.rank()==0){
			std::cout<< " , because timeNext > opNext , dt set to : " << dt*t_ast << std::endl;
			std::cout<< std::flush;
			}

			dtFlag = 0;
		}
		dtFlag += 1;

		if(helper.rank()==0){
		std::cout<< " , dt  : " << dt*t_ast << std::endl;
		std::cout<<" "<< std::endl;
		std::cout << " READY FOR NEXT ITERATION. " << std::endl;
		std::cout<< std::flush;
		}

	}

}
#endif /* PROBLEM09_DRIVER_COUPLED_HH_ */
