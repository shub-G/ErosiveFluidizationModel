/*
 * driver.hh
 *
 *  Created on: Apr 6, 2021
 *      Author: sgupta
 */

#ifndef PROBLEM09_DRIVER_DECOUPLED_HH_
#define PROBLEM09_DRIVER_DECOUPLED_HH_

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
	using VBEFV =Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
	using OrderingTagFV = Dune::PDELab::EntityBlockedOrderingTag;

	// TWO PHASE FLOW MODEL
	using GFS_2P = Dune::PDELab::PowerGridFunctionSpace< GFSP0,2/*num eqns*/,VBEFV,OrderingTagFV >;
	GFS_2P gfs_2p(gfsp0);
	typedef typename GFS_2P::template ConstraintsContainer<double>::Type C_2P;
	C_2P c_2p;
	c_2p.clear();
	gfs_2p.update();
	std::cout << "degrees of freedom (two phase flow model): " << gfs_2p.globalSize() << std::endl;

	// INTERNAL EROSION MODEL
	using GFS_Er = Dune::PDELab::PowerGridFunctionSpace< GFSP0,2/*num eqns*/,VBEFV,OrderingTagFV >;
	GFS_Er gfs_er(gfsp0);
	typedef typename GFS_Er::template ConstraintsContainer<double>::Type C_Er;
	C_Er c_er;
	c_er.clear();
	gfs_er.update();
	std::cout << "degrees of freedom (internal erosion model): " << gfs_er.globalSize() << std::endl;

	// L2-PROJECTION OF P0-PRESSURES ONTO Q1-SPACE
	const int degree = 1; //polynomial degree
	typedef Dune::PDELab::QkLocalFiniteElementMap<GV,Coord,double,degree> FEMQ1;
	FEMQ1 femq1(gv);
	typedef Dune::PDELab::GridFunctionSpace<GV,FEMQ1,CON0,VBE0> GFSQ1;
	GFSQ1 gfsq1(gv,femq1);
	using VBEL2 = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
	using OrderingTagL2 = Dune::PDELab::EntityBlockedOrderingTag;
	typedef Dune::PDELab::PowerGridFunctionSpace< GFSQ1,2/*Pw,Pn*/,VBEL2,OrderingTagL2> GFS_L2;
	GFS_L2 gfs_l2(gfsq1);
	typedef typename GFS_L2::template ConstraintsContainer<double>::Type C_L2;
	C_L2 c_l2;
	c_l2.clear();

	/****************************************************/
	// VECTOR CONTAINERS
	/****************************************************/
	using U2P = Dune::PDELab::Backend::Vector<GFS_2P,double>;
	U2P u2p_old(gfs_2p);
	U2P u2p_new(gfs_2p);
	U2P u2p_dot(gfs_2p);
	using UEr = Dune::PDELab::Backend::Vector<GFS_Er,double>;
	UEr uer_old(gfs_er);
	UEr uer_new(gfs_er);
	UEr uer_dot(gfs_er);
	using UL2 = Dune::PDELab::Backend::Vector<GFS_L2,double>;
	UL2 ul2(gfs_l2);

	/****************************************************/
	// INITIAL CONDITIONS
	/****************************************************/
	using Initial = InitialConditions<GV,PTree,Parameters>;
	Initial initial(gv,ptree,parameter);
	auto pw_ic_local = [&](const auto& e, const auto& x){return initial.Pw(e,x);};
	auto pw_ic = Dune::PDELab::makeGridFunctionFromCallable(gv,pw_ic_local);
	auto sw_ic_local = [&](const auto& e, const auto& x){return initial.Sw(e,x);};
	auto sw_ic = Dune::PDELab::makeGridFunctionFromCallable(gv,sw_ic_local);
	using  InitialConditions2P = Dune::PDELab::CompositeGridFunction<decltype(pw_ic),decltype(sw_ic)>;
	InitialConditions2P ic_2p( pw_ic,sw_ic );
	auto cf_ic_local = [&](const auto& e, const auto& x){return initial.Cf(e,x);};
	auto cf_ic = Dune::PDELab::makeGridFunctionFromCallable(gv,cf_ic_local);
	auto por_ic_local= [&](const auto& e, const auto& x){return initial.porosity(e,x);};
	auto por_ic = Dune::PDELab::makeGridFunctionFromCallable(gv,por_ic_local);
	using  InitialConditionsEr = Dune::PDELab::CompositeGridFunction<decltype(por_ic),decltype(cf_ic)>;
	InitialConditionsEr ic_er( por_ic,cf_ic );

	// 	Initialize the solution at t=0 (uold) with the given initial values
	Dune::PDELab::interpolate( ic_2p, gfs_2p, u2p_old );
	u2p_new = u2p_old;
	Dune::PDELab::interpolate( ic_er, gfs_er, uer_old );
	uer_new = uer_old;

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
	using LOP_2P = LocalOperatorTwoPhaseFlow<GV,Parameters,GFS_2P,U2P,GFS_Er,UEr,WBC,NBC>;
	LOP_2P lop_2p( gv, parameter,
				   gfs_2p, &u2p_old,
				   gfs_er, &uer_new, &uer_dot,
				   w_bc, n_bc,
				   &time, &dt,
				   intorder );

	typedef LocalOperatorL2Projection< GV,Parameters,GFS_2P,U2P,GFS_Er,UEr> LOP_L2;
	LOP_L2 lop_l2(gv, parameter, gfs_2p, &u2p_new, gfs_er, &uer_new);

	using LOP_Er = LocalOperatorErosion<GV,Parameters,GFS_Er,UEr,GFS_2P,U2P,GFS_L2,UL2,FSBC,WBC,NBC>;
	LOP_Er lop_er( gv, parameter,
				   gfs_er, &uer_old,
			 	   gfs_2p, &u2p_new, &u2p_dot,
				   gfs_l2, &ul2,
				   fs_bc, w_bc, n_bc,
				   &time, &dt,
				   intorder );

	/****************************************************/
	// GRID OPERATOR
	/****************************************************/
	using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
	MBE mbe_2p(20); // Maximal number of nonzeroes per row
	MBE mbe_l2(20); // Maximal number of nonzeroes per row
	MBE mbe_er(20); // Maximal number of nonzeroes per row

	using GOLOP_2P = Dune::PDELab::GridOperator<GFS_2P,GFS_2P,LOP_2P,MBE,double,double,double,C_2P,C_2P>;
	GOLOP_2P go_lop_2p(gfs_2p,c_2p,gfs_2p,c_2p,lop_2p,mbe_2p);
	typename GOLOP_2P::Traits::Jacobian jac_2p(go_lop_2p);
	if(helper.rank()==0){
		std::cout << jac_2p.patternStatistics() << std::endl;
	}

	using GOLOP_L2 = Dune::PDELab::GridOperator< GFS_L2, GFS_L2, LOP_L2, MBE, double, double, double, C_L2, C_L2>;
	GOLOP_L2 go_lop_l2(gfs_l2, c_l2, gfs_l2, c_l2, lop_l2, mbe_l2);

	using GOLOP_Er = Dune::PDELab::GridOperator<GFS_Er,GFS_Er,LOP_Er,MBE,double,double,double,C_Er,C_Er>;
	GOLOP_Er go_lop_er(gfs_er,c_er,gfs_er,c_er,lop_er,mbe_er);
	typename GOLOP_Er::Traits::Jacobian jac_er(go_lop_er);
	if(helper.rank()==0){
		std::cout << jac_er.patternStatistics() << std::endl;
	}

	/****************************************************/
	// LINEAR SOLVER
	/****************************************************/
	using LS_2P = Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GOLOP_2P>;
	LS_2P ls_2p(gfs_2p,1000,1);
	Dune::Amg::Parameters ls_params_2p = ls_2p.parameters();
	ls_params_2p.setCoarsenTarget(1000000);// max DoF at coarsest level
	ls_2p.setParameters(ls_params_2p);

	using LS_L2 = Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GOLOP_L2>;
	LS_L2 ls_l2(gfs_l2,1000,1);
	//Dune::Amg::Parameters ls_params_l2 = ls_l2.parameters();
	//ls_params_l2.setCoarsenTarget(100000);// max DoF at coarsest level
	//ls_l2.setParameters(ls_params_l2);

	using LS_Er = Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GOLOP_Er>;
	LS_Er ls_er(gfs_er,1000,1);
	Dune::Amg::Parameters ls_params_er = ls_er.parameters();
	ls_params_er.setCoarsenTarget(1000000);// max DoF at coarsest level
	ls_er.setParameters(ls_params_er);

	/****************************************************/
	// NON-LINEAR SOLVER
	/****************************************************/
	using PDESOLVER_2P = Dune::PDELab::Newton<GOLOP_2P,LS_2P,U2P>;
	PDESOLVER_2P pdesolver_2p(go_lop_2p,ls_2p);
	pdesolver_2p.setLineSearchStrategy(PDESOLVER_2P::Strategy::noLineSearch);
	pdesolver_2p.setReassembleThreshold(0.0);
	pdesolver_2p.setVerbosityLevel(ptree.get("newton.flow.verbosity",(int)2));
	pdesolver_2p.setReduction(ptree.get("newton.flow.reduction",(double)1e-6));
	pdesolver_2p.setMinLinearReduction(ptree.get("newton.flow.min_lin_reduction",(double)1e-9));
	pdesolver_2p.setMaxIterations(ptree.get("newton.flow.max_iterations",(int)15));
	pdesolver_2p.setForceIteration(ptree.get("newton.flow.force_iteration",(bool)false));
	pdesolver_2p.setAbsoluteLimit(ptree.get("newton.flow.abs_error",(double)1.e-6));

	typedef Dune::PDELab::StationaryLinearProblemSolver<GOLOP_L2,LS_L2,UL2> SLP_L2;
	SLP_L2 slp_l2(go_lop_l2,ls_l2,ul2,1e-10);
	slp_l2.apply(); // get initial values of ul2

	using PDESOLVER_Er = Dune::PDELab::Newton<GOLOP_Er,LS_Er,UEr>;
	PDESOLVER_Er pdesolver_er(go_lop_er,ls_er);
	pdesolver_er.setLineSearchStrategy(PDESOLVER_Er::Strategy::noLineSearch);
	pdesolver_er.setReassembleThreshold(0.0);
	pdesolver_er.setVerbosityLevel(ptree.get("newton.soil.verbosity",(int)2));
	pdesolver_er.setReduction(ptree.get("newton.soil.reduction",(double)1e-6));
	pdesolver_er.setMinLinearReduction(ptree.get("newton.soil.min_lin_reduction",(double)1e-9));
	pdesolver_er.setMaxIterations(ptree.get("newton.soil.max_iterations",(int)15));
	pdesolver_er.setForceIteration(ptree.get("newton.soil.force_iteration",(bool)false));
	pdesolver_er.setAbsoluteLimit(ptree.get("newton.soil.abs_error",(double)1.e-6));

	/****************************************************/
	// VTK
	/****************************************************/
	int subsampling = 1;
	using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
	VTKWRITER vtkwriter(gv,Dune::refinementIntervals(subsampling));
	using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
	VTKSEQUENCEWRITER vtkSequenceWriter( std::make_shared<VTKWRITER>(vtkwriter),fileName,pathName,"");

	using SUBGFS_Pw = typename Dune::PDELab::GridFunctionSubSpace< GFS_2P, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Pw>> >;
	SUBGFS_Pw subgfs_pw(gfs_2p);
	using DGF_Pw = Dune::PDELab::DiscreteGridFunction< SUBGFS_Pw, U2P >;
	DGF_Pw dgf_pw( subgfs_pw, u2p_new );
	vtkSequenceWriter.addCellData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_Pw> >(dgf_pw,"pw"));

	using SUBGFS_Sw = typename Dune::PDELab::GridFunctionSubSpace< GFS_2P, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Sw>> >;
	SUBGFS_Sw subgfs_sw(gfs_2p);
	using DGF_Sw = Dune::PDELab::DiscreteGridFunction< SUBGFS_Sw, U2P >;
	DGF_Sw dgf_sw( subgfs_sw, u2p_new );
	vtkSequenceWriter.addCellData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_Sw> >(dgf_sw,"sw"));

	using SUBGFS_Pw_l2 = typename Dune::PDELab::GridFunctionSubSpace< GFS_L2, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Pw_l2>> >;
	SUBGFS_Pw_l2 subgfs_pw_l2(gfs_l2);
	using DGF_Pw_l2 = Dune::PDELab::DiscreteGridFunction< SUBGFS_Pw_l2, UL2 >;
	DGF_Pw_l2 dgf_pw_l2( subgfs_pw_l2, ul2 );
	vtkSequenceWriter.addVertexData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_Pw_l2> >(dgf_pw_l2,"pw_l2"));

	using SUBGFS_Pn_l2 = typename Dune::PDELab::GridFunctionSubSpace< GFS_L2, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Pw_l2>> >;
	SUBGFS_Pn_l2 subgfs_pn_l2(gfs_l2);
	using DGF_Pn_l2 = Dune::PDELab::DiscreteGridFunction< SUBGFS_Pn_l2, UL2 >;
	DGF_Pn_l2 dgf_pn_l2( subgfs_pn_l2, ul2 );
	vtkSequenceWriter.addVertexData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_Pn_l2> >(dgf_pn_l2,"pn_l2"));

	using SUBGFS_por = typename Dune::PDELab::GridFunctionSubSpace< GFS_Er, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::porosity>> >;
	SUBGFS_por subgfs_por(gfs_er);
	using DGF_por = Dune::PDELab::DiscreteGridFunction< SUBGFS_por, UEr >;
	DGF_por dgf_por( subgfs_por, uer_new );
	vtkSequenceWriter.addCellData(std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF_por> >(dgf_por,"porosity"));

	using SUBGFS_Cf = typename Dune::PDELab::GridFunctionSubSpace< GFS_Er, Dune::TypeTree::HybridTreePath<Dune::index_constant<PrimaryVariables::Cf>> >;
	SUBGFS_Cf subgfs_cf(gfs_er);
	using DGF_Cf = Dune::PDELab::DiscreteGridFunction< SUBGFS_Cf, UEr >;
	DGF_Cf dgf_cf( subgfs_cf, uer_new );
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
	int newton_iterations_2p = 0;
	int newton_iterations_er = 0;
	double fp_abs_error = ptree.get("fixed_point_iteration.abs_error",(double)1.e-9);
	int fp_max_iter  = ptree.get("fixed_point_iteration.max_iterations",(int)20);
	int fp_min_iter = ptree.get("fixed_point_iteration.min_iterations",(int)1);
	double dtF_min = ptree.get("adaptive_time_control.dtF_min",(double)0.25);
	double dtF_max = ptree.get("adaptive_time_control.dtF_max",(double)0.25);

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
			u2p_dot=0.;
			uer_dot=0.;
			newton_iterations_2p=0;
			newton_iterations_er=0;
			for( int fp = 0.; fp<fp_max_iter; fp++ ){
				if(helper.rank()==0){
				std::cout<< "+++++++++++++++"<< std::endl;
				std::cout<< "FP-num: " << fp << std::endl;
				std::cout<< "+++++++++++++++"<< std::endl;
				}

				if(helper.rank()==0){
				std::cout<< '\n' <<"Solving 2-Phase Flow model:" << std::endl;
				}
				pdesolver_2p.apply( u2p_new );

				if(newton_iterations_2p<pdesolver_2p.result().iterations)
					newton_iterations_2p = pdesolver_2p.result().iterations;

				u2p_dot  = u2p_new;
				u2p_dot -= u2p_old;
				u2p_dot *= 1./dt;

				if(helper.rank()==0){
				std::cout<< '\n' << "Projecting pw, pn:" << std::endl;
				}
				slp_l2.apply();

				if(helper.rank()==0){
				std::cout<< '\n' << "Solving Erosion model:" << std::endl;
				}

				pdesolver_er.apply( uer_new );

				if(newton_iterations_er<pdesolver_er.result().iterations)
					newton_iterations_er = pdesolver_er.result().iterations;

//				uer_dot  = uer_new;
//				uer_dot -= uer_old;
//				uer_dot *= 1./dt;

				if( fp>=fp_min_iter and pdesolver_2p.result().first_defect <fp_abs_error and pdesolver_er.result().first_defect<fp_abs_error ){
					fp=fp_max_iter ;
				}
			}

			newton_iterations = std::max(newton_iterations_2p,newton_iterations_er);

			exceptionCaught = false;

		}catch ( Dune::Exception &e ) {
			exceptionCaught = true;
			if( dt > 1e-15 ){

				if(helper.rank()==0){
				std::cout << "Catched Error, Dune reported error: " << e << std::endl;
				}
				u2p_new = u2p_old;
				uer_new = uer_old;
				newton_iterations_2p = 0;
				newton_iterations_er = 0;
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
		u2p_old = u2p_new;
		uer_old = uer_new;
		//2. ADVANCE TIME:
		time += dt;

		if(helper.rank()==0){
		std::cout<<" "<< std::endl;
		std::cout<< " time = " << time ;
		std::cout<< std::flush;
		}

		if( is_adaptive_time_control ){
			if(newton_iterations>maxAllowableIterations){
				dt=std::max(dt*(1.0-dtF_min),dt_min);
			}
			else if(newton_iterations<=minAllowableIterations){
				dt=std::min(dt*(1.0+dtF_max),dt_max);
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

#endif /* PROBLEM09_DRIVER_DECOUPLED_HH_ */
