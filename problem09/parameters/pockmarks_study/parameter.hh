/*
 * parameters.hh
 *
 *  Created on: Apr 6, 2021
 *      Author: sgupta
 */

#ifndef PROBLEM09_PARAMETERS_TEST_2PFLOW_PARAMETERS_HH_
#define PROBLEM09_PARAMETERS_TEST_2PFLOW_PARAMETERS_HH_

template<typename GV,typename PTree>
class PockmarksStudyParameters
{
private:
	  const GV& gv;
	  const PTree& ptree;
	  const static int dim = GV::dimension;
	  double eps = 1.0e-9;
	  using FieldVector = typename Dune::FieldVector<double,dim>;
	  using FieldMatrix = typename Dune::FieldMatrix<double,dim,dim>;

	  // characteristic values
	  double x_ast;
	  double t_ast;
	  double p_ast;

	  //domain
	  double L;//length of domain
	  double H;//heightof domain
	  double l;//center pos of inlet on the bottom boundary
	  double dl;//width of bottom inlet
	  double Hp;//height of porous domain

	  //gravity
	  bool gravity_flag;
	  double g_magnitude;

	  //(w)etting and (n)on-wetting phase properties
	  double compressibility_w;
	  double compressibility_n;
	  double compressibility_s;
	  double density_w;
	  double exp_rhow;
	  double density_n;
	  double density_s;
	  double viscosity_w;
	  double viscosity_n;

	  // Capillary pressure (Brooks-Corey) parameters
	  double pc_max0;
	  double pe_bc0;
	  double lambda_bc0;
	  double beta_kc0;

	  //sediment properties
	  double porosity_0;
	  double K0_min;
	  double K0_max;
	  double KF;
	  double expK0;
	  double A0_min;
	  double A0_max;
	  double a0;
	  double expA0;
	  double dhz;

	  //internal erosion parameters
	  double ler_w0;
	  double ler_n0;
	  double vw_mag_cr0;
	  double vn_mag_cr0;

	  //sediment deposition parameters
	  double ds0;
	  double n0;

public:
  	//! construct from grid view
	  PockmarksStudyParameters (const GV& gv_, const PTree& ptree_)
	: gv( gv_ ) ,
	  ptree(ptree_)
  	{
		  //characteristic values
		  x_ast  = ptree.get("characteristic_values.length"	 ,(double)1.0);//m
		  t_ast  = ptree.get("characteristic_values.time"	 ,(double)1.0);//s
		  p_ast  = ptree.get("characteristic_values.pressure",(double)1.0);//Pa

		  //domain
		  L  = ptree.get("domain.length",(double)1.0);//m
		  L *= 1./x_ast;//ndim
		  H  = ptree.get("domain.height",(double)1.0);//m
		  H *= 1./x_ast;//ndim
		  l  = ptree.get("domain.bottom_inlet.center" ,(double)0.5);//m
		  l *= 1./x_ast;//ndim
		  dl = ptree.get("domain.bottom_inlet.length" ,(double)0.1);//m
		  dl*= 1./x_ast;//ndim
		  Hp  = ptree.get("domain.porous_domain_height",(double)1.0);//m
		  Hp *= 1./x_ast;//ndim

		  //gravity
		  gravity_flag = ptree.get("gravity.flag"		,(bool)false);
		  g_magnitude = ptree.get("gravity.magnitude"	,(double)9.81); //m/s^2
		  g_magnitude *= x_ast/p_ast; //

		  //(w)etting and (n)on-wetting phase properties
		  viscosity_w = ptree.get("phase.wetting.viscosity"		,(double)1.e-3);//Pa.s
		  viscosity_w *= 1./(p_ast*t_ast);//ndim
		  viscosity_n = ptree.get("phase.nonwetting.viscosity"	,(double)1.e-5);//Pa.s
		  viscosity_n *= 1./(p_ast*t_ast);//ndim
		  density_w = ptree.get("phase.wetting.density"			,(double)1000.0);//kg/m^3
		  exp_rhow  = ptree.get("phase.wetting.density_exponent",(double)0.025);//--
		  density_n = ptree.get("phase.nonwetting.density"		,(double)1.0);	 //kg/m^3
		  density_s = ptree.get("phase.soil.density"			,(double)2000.0);//kg/m^3
		  compressibility_w = ptree.get("phase.wetting.compressibility"		,(double)0.0);//1/Pa
		  compressibility_w *= p_ast;//ndim
		  compressibility_n = ptree.get("phase.nonwetting.compressibility"	,(double)0.0);//1/Pa
		  compressibility_n *= p_ast;//ndim
		  compressibility_s = ptree.get("phase.soil.compressibility"		,(double)0.0);//1/Pa
		  compressibility_s *= p_ast;//ndim

		  //van-Genuchten parameters
		  pc_max0    = ptree.get("sediment.capillary_pressure_function.pcmax0"	,(double)1.e8);//Pa
		  pc_max0   *= 1./p_ast;//ndim
		  pe_bc0     = ptree.get("sediment.capillary_pressure_function.pe0"		,(double)5000.0);//Pa
		  pe_bc0    *= 1./p_ast;//ndim
		  lambda_bc0 = ptree.get("sediment.capillary_pressure_function.lambda0"	,(double)2.0);//-
		  beta_kc0   = ptree.get("sediment.capillary_pressure_function.beta0"	,(double)3.0);//-

		  //hydraulic properties
		  porosity_0 	 = ptree.get("sediment.porosity.phi0"	,(double)0.3);//-
		  K0_min = ptree.get("sediment.permeability.K0_sediment" ,(double)1.e-12);//m^2
		  K0_min *= 1./(x_ast*x_ast);//ndim
		  K0_max = ptree.get("sediment.permeability.K0_surface" ,(double)1.e-9);//m^2
		  K0_max *= 1./(x_ast*x_ast);//ndim
		  KF = ptree.get("sediment.permeability.Kxx_by_Kzz" ,(double)1.);//
		  expK0 = ptree.get("sediment.permeability.expK" ,(double)1.);//
		  A0_min = ptree.get("sediment.permeability.A0_min" ,(double)0.0);//-
		  A0_max = ptree.get("sediment.permeability.A0_max" ,(double)0.0);//-
		  a0 = ptree.get("sediment.permeability.a0" ,(double)1.0);//-
		  expA0 = ptree.get("sediment.permeability.expA" ,(double)1.);//
		  dhz = ptree.get("sediment.permeability.transition_depth" ,(double)0.);//

		  //internal erosion parameters
		  ler_w0 = ptree.get("sediment.erosion.wetting.l0" 		   ,(double)0.0);//1/m
		  ler_w0 *= x_ast;//ndim
		  ler_n0 = ptree.get("sediment.erosion.nonwetting.l0" 	   ,(double)0.0);//1/m
		  ler_n0 *= x_ast;//ndim
		  vw_mag_cr0 = ptree.get("sediment.erosion.wetting.vcr0"   ,(double)0.0);//m/s
		  vw_mag_cr0 *= t_ast/x_ast;//ndim
		  vn_mag_cr0 = ptree.get("sediment.erosion.nonwetting.vcr0",(double)0.0);//m/s
		  vn_mag_cr0 *= t_ast/x_ast;//ndim

		  //sediment deposition parameters
		  ds0 = ptree.get("sediment.deposition.ds0",(double)0.0);
		  n0 = ptree.get("sediment.deposition.n0",(double)0.0);
  	}

	/****************************************************
	 *  GRAVITY VECTOR
	 ****************************************************/
	FieldVector g( ) const {
		FieldVector gravity( 0. );
		double g = 0.;
		if(gravity_flag) g = g_magnitude;
		gravity[dim-1] = g;
		return gravity;
	}

	/****************************************************
	 * COMPUTATIONAL DOMAIN
	 ****************************************************/
	bool isLeftBoundary(FieldVector globalpos)const{
		if(globalpos[0]<0.+eps) return true;
		else return false;
	}

	bool isRightBoundary(FieldVector globalpos) const {
		if(globalpos[0]>L-eps) return true;
		else return false;
	}

	#if DIMENSION==3
	bool isFrontBoundary(FieldVector globalpos)const{
		if(globalpos[1]<0.+eps) return true;
		else return false;
	}

	bool isBackBoundary(FieldVector globalpos) const {
		if(globalpos[1]>L-eps) return true;
		else return false;
	}
	#endif

	bool isBottomBoundary(FieldVector globalpos) const {
		if(globalpos[dim-1]<0.+eps) return true;
		else return false;
	}

	bool isTopBoundary(FieldVector globalpos) const {
		if(globalpos[dim-1]>H-eps) return true;
		else return false;
	}

	bool isBottomInletBoundary(FieldVector globalpos) const {
		if( isBottomBoundary(globalpos)
			and globalpos[0]>l-dl-eps
			and globalpos[0]<l+dl+eps
			#if DIMENSION==3
			and globalpos[1]>l-dl-eps
			and globalpos[1]<l+dl+eps
			#endif
			)
			return true;
		else return false;
	}

	bool isPorousDomain(FieldVector globalpos)const{
		if(globalpos[dim-1]<Hp+eps) return true;
		else return false;
	}

	/****************************************************
	 * PHASE PROPERTIES
	 ****************************************************/

	double WettingPhaseDensity( double pw, double cf ) const{
		return density_w*std::exp(exp_rhow*std::max(cf,0.0));
	}

	double NonwettingPhaseDensity( double pn ) const{
		return density_n;
	}

	double SoilPhaseDensity( double peff ) const{
		return density_s;
	}

	double WettingPhaseCompressibility() const{
		return compressibility_w;
	}

	double NonwettingPhaseCompressibility() const{
		return compressibility_n;
	}

	double SoilPhaseCompressibility() const{
		return compressibility_s;
	}

	double WettingPhaseViscosity() const{
		return viscosity_w;
	}

	double NonwettingPhaseViscosity() const{
		return viscosity_n;
	}

	/****************************************************
	 * VOLUMETRIC SOURCE TERMS
	 ****************************************************/

	double WettingVolumetricSource(FieldVector globalpos) const{
		return t_ast*0.;/*ndim*/
	}

	double NonwettingVolumetricSource(FieldVector globalpos) const{
		return t_ast*0.;/*ndim*/
	}

	double ErosionRate(FieldVector globalpos,
					   double vw_mag,
					   double vn_mag,
					   double porosity )const {
		double qs_w=0.;
		if(vw_mag>vw_mag_cr0){
			qs_w = (-1.) * t_ast*ler_w0/*ndim*/ * (1.-porosity) * (vw_mag/*ndim*/-vw_mag_cr0/*ndim*/);
		}
		double qs_n=0.;
		if(vn_mag>vn_mag_cr0){
			qs_n = (-1.) * t_ast*ler_n0/*ndim*/ * (1.-porosity) * (vn_mag/*ndim*/-vn_mag_cr0/*ndim*/);
		}
		return qs_w+qs_n;/*ndim*/
	}

	double DepositionRate(FieldVector globalpos, double Cf, double Sw, double porosity)const{
		double ds=0.;
		ds = t_ast*ds0*std::pow(porosity*Cf/Sw,n0);
		return ds;
	}

	/****************************************************
	 * SEDIMENT AND HYDRAULIC PROPERTIES
	 ****************************************************/

	double InitialPorosity(FieldVector globalpos)const{
		if( isPorousDomain(globalpos) ) return porosity_0;
		else return 1.0;
	}

	FieldVector
	Permeability(FieldVector globalpos, double porosity) const {

		// double KF_max=100.0;
		// double K0xx=std::min( K0_min*KF+(K0_max*KF_max-K0_min*KF)*std::exp(-expK0*(Hp+dhz-globalpos[dim-1])),K0_max*KF_max);
		// double K0zz=std::min( K0_min   +(K0_max       -K0_min   )*std::exp(-expK0*(Hp+dhz-globalpos[dim-1])),K0_max       );
		// double A0=std::min( A0_min+(A0_max-A0_min)*std::exp(-expA0*(Hp+dhz-globalpos[dim-1])),A0_max);
		// double f0=std::exp( a0*A0*(porosity-porosity_0)/(1.-porosity_0));
		// double f1=std::exp(    A0*(porosity-porosity_0)/(1.-porosity_0));
		
		// FieldVector K(0.);
		// K[0] = K0xx*f0;
		// K[dim-1] = K0zz*f1;

		double f0=std::exp( a0*A0_min*(porosity-porosity_0)/(1.-porosity_0));
		double f1=std::exp(    A0_min*(porosity-porosity_0)/(1.-porosity_0));
		double F=std::min( KF+(1000.0-KF)*std::exp(-expK0*(Hp+dhz-globalpos[dim-1])),1000.0);
		FieldVector K(0.);
		K[0] = K0_min*f0*F;
		#if DIMENSION==3
		K[1] = K0_min*f0*F;
		#endif
		K[dim-1] = K0_min*f1;

		return K; //ndim
	}

	FieldMatrix
	PermeabilityTensor(FieldVector globalpos, double porosity)const{
		auto k = Permeability(globalpos,porosity);

		FieldMatrix K(0.);
		for(int d=0; d<dim; d++) K[d][d] = k[d];

		return K; //ndim
	}

	double CapillaryPressure(FieldVector globalpos, double Sw, double porosity)const{
		double f=std::pow( (1.-porosity)/(1.-porosity_0) , beta_kc0 );
		double pc = std::min( pe_bc0 * f * std::pow(Sw,-(1./lambda_bc0)) , pc_max0 );
		return pc;
	}

	double dPc_dSw(FieldVector globalpos, double Sw, double porosity)const{
		double pc = CapillaryPressure(globalpos,Sw,porosity);
		double swmin = std::pow(pe_bc0/pc_max0,lambda_bc0);
		double dpc = 0.;
		if( Sw>swmin and pc>0. ) dpc = -pc/(lambda_bc0*Sw);
		return dpc;
	}

	double dPc_dPor(FieldVector globalpos, double Sw, double porosity)const{
		double pc = CapillaryPressure(globalpos,Sw,porosity);
		double dpc = -(beta_kc0*pc);
		if(porosity<1.0) dpc *= 1./(1.-porosity);
		return dpc;
//		return 0.;
	}

	double WettingRelativePermeability(FieldVector globalpos, double Sw) const{
		double kr = std::pow( std::max(Sw,0.0) , (2.0/lambda_bc0 + 3.0) );
		return std::min(kr,1.0) ;
	}

	double NonwettingRelativePermeability(FieldVector globalpos, double Sw) const{
		double kr = std::pow(1.0-std::min(Sw,1.0) , 2.0)
				  * ( 1.0 - std::pow( std::min(Sw,1.0) , (2.0/lambda_bc0 + 1.0) ) );
		return std::max(kr,0.0) ;
	}

};

template<typename GV,typename PTree,typename Parameters>
class InitialConditions
{
private:
	const GV& gv;
	const PTree& ptree;
	const Parameters& param;
	const static int dim = GV::dimension;
	using FieldVector = typename Dune::FieldVector<double,dim>;

	double x_ast;
	double t_ast;
	double p_ast;

	double pw;
	double sw;
	double cf;

public:
	//! construct from grid view
	InitialConditions ( const GV& gv_ ,
				 	 	const PTree& ptree_,
						const Parameters& param_)
	: gv( gv_ ) ,
	  ptree(ptree_),
	  param(param_)
	{
		//characteristic values
		x_ast  = ptree.get("characteristic_values.length"	,(double)1.0);//m
		t_ast  = ptree.get("characteristic_values.time"		,(double)1.0);//s
		p_ast  = ptree.get("characteristic_values.pressure"	,(double)1.0);//Pa

		pw = ptree.get("initial.Pw"	,(double)1.e5);//Pa
		pw *= 1./p_ast;//ndim
		sw = ptree.get("initial.Sw",(double)0.01);//-
		cf = ptree.get("initial.Cf",(double)0.0);//?
	}

	template<typename E, typename X>
	double Pw (const E& e, const X& x) const {
		auto xglobal = e.geometry().global(x);
		return pw;//ndim
	}

	template<typename E, typename X>
	double Sw (const E& e, const X& x) const {
		auto xglobal = e.geometry().global(x);
		return sw;//-
	}

	template<typename E, typename X>
	double porosity (const E& e, const X& x) const {
		auto xglobal = e.geometry().global(x);
		return param.InitialPorosity(xglobal);//-
	}

	template<typename E, typename X>
	double Cf (const E& e, const X& x) const {
		auto xglobal = e.geometry().global(x);
		return cf;//?
	}

};

struct ConvectionDiffusionBoundaryConditions
{
  enum Type { Dirichlet=1, Neumann=-1, Outflow=-2, None=-3 };
};

template<typename GV,typename PTree,typename Parameters>
class WettingPhaseBoundaryCondition
{
	typedef ConvectionDiffusionBoundaryConditions::Type BCType;

private:
	const GV& gv;
	const PTree& ptree;
	const Parameters& param;
	const static int dim = GV::dimension;
	using FieldVector = typename Dune::FieldVector<double,dim>;

	double x_ast;
	double t_ast;
	double p_ast;

	int bcset;

	double pw_top;
	double sw_top;
	double dpw_inB;
	double sw_inB;
	double cf_top;
	double t_cutoff;
	double H;

public:
	//! construct from grid view
	WettingPhaseBoundaryCondition ( const GV& gv_ ,
						  	  	    const PTree& ptree_,
									const Parameters& param_)
	: gv( gv_ ) , ptree(ptree_), param(param_)
	{
		//characteristic values
		x_ast  = ptree.get("characteristic_values.length"	,(double)1.0);//m
		t_ast  = ptree.get("characteristic_values.time"		,(double)1.0);//s
		p_ast  = ptree.get("characteristic_values.pressure"	,(double)1.0);//Pa

		bcset = ptree.get("boundary.case",(int) 0 );//-

		pw_top = ptree.get("boundary.top.Pw" ,(double)1.e5);//Pa
		pw_top *= 1./p_ast;//ndim
		sw_top  = ptree.get("boundary.top.Sw",(double)1.0);//-

		dpw_inB = ptree.get("boundary.bottom_inlet.dPw",(double)2.e5);//Pa
		dpw_inB *= 1./p_ast;//ndim
		sw_inB = ptree.get("boundary.bottom_inlet.Sw",(double)0.99);//-

		cf_top = ptree.get("boundary.top.Cf",(double)0.);//kg_soil/kg_water

		t_cutoff = ptree.get("boundary.bottom_inlet.cutoff_time",(double)1.e15);//s

		H  = ptree.get("domain.height",(double)1.0);//m
		H *= 1./x_ast;//ndim
	}

	//! boundary condition type function (true = Dirichlet)
	template<typename I, typename X>
	BCType type (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);

		if( bcset == 1 ){
			if( param.isBottomBoundary(xglobal) or param.isTopBoundary(xglobal) )
				return ConvectionDiffusionBoundaryConditions::Dirichlet;
			else return ConvectionDiffusionBoundaryConditions::Neumann;
		}else if( bcset == 2 ){
			if( param.isBottomInletBoundary(xglobal) )
				return ConvectionDiffusionBoundaryConditions::Dirichlet;
			else return ConvectionDiffusionBoundaryConditions::Neumann;
		}else{//default
			if( param.isBottomInletBoundary(xglobal) or param.isTopBoundary(xglobal) )
				return ConvectionDiffusionBoundaryConditions::Dirichlet;
			else return ConvectionDiffusionBoundaryConditions::Neumann;
		}
	}

	//! Dirichlet value
	//! Specifies Pw|_D
	template<typename I, typename X>
	double dirichlet_value (const I& i, const X& x, double& t, double& cf) const {
		auto xglobal = i.geometry().global(x);

		if( bcset == 1 ){
			if( param.isBottomInletBoundary(xglobal) )
				return pw_top + param.WettingPhaseDensity(pw_top,cf)*param.g()[dim-1]*(H-xglobal[dim-1]) + dpw_inB;
			else if(param.isBottomBoundary(xglobal) and !param.isBottomInletBoundary(xglobal))
				return pw_top + param.WettingPhaseDensity(pw_top,cf)*param.g()[dim-1]*(H-xglobal[dim-1]);
			else return pw_top;
		}else if( bcset == 2 ){
			return pw_top + param.WettingPhaseDensity(pw_top,cf)*param.g()[dim-1]*(H-xglobal[dim-1]) + dpw_inB;
		}else{//default
			if( param.isBottomInletBoundary(xglobal) ){
				if( t<t_cutoff ) return pw_top + param.WettingPhaseDensity(pw_top,cf)*param.g()[dim-1]*(H-xglobal[dim-1]) + dpw_inB;
				else return pw_top + param.WettingPhaseDensity(pw_top,cf)*param.g()[dim-1]*(H-xglobal[dim-1]);
			}
			else return pw_top;
		}
	}

	//! Neumann boundary condition
	template<typename I, typename X>
	FieldVector
	neumann_value (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);
		FieldVector jF(0.0);
		return jF;
	}

	//! Outflow boundary condition
	template<typename I, typename X>
	FieldVector
	outflow_value (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);
		FieldVector oF(0.0);
		return oF;
	}

};

template<typename GV,typename PTree,typename Parameters>
class NonwettingPhaseBoundaryCondition
{
	typedef ConvectionDiffusionBoundaryConditions::Type BCType;

private:
	const GV& gv;
	const PTree& ptree;
	const Parameters& param;
	const static int dim = GV::dimension;
	using FieldVector = typename Dune::FieldVector<double,dim>;

	double x_ast;
	double t_ast;
	double p_ast;

	int bcset;

	double sw_inB;
	double sw_top;

public:
	//! construct from grid view
	NonwettingPhaseBoundaryCondition ( const GV& gv_ ,
						  	  	  	   const PTree& ptree_,
									   const Parameters& param_)
	: gv( gv_ ) , ptree(ptree_), param(param_)
	{
		//characteristic values
		x_ast  = ptree.get("characteristic_values.length"	,(double)1.0);//m
		t_ast  = ptree.get("characteristic_values.time"		,(double)1.0);//s
		p_ast  = ptree.get("characteristic_values.pressure"	,(double)1.0);//Pa

		bcset = ptree.get("boundary.case",(int) 0 );//-

		sw_inB = ptree.get("boundary.bottom_inlet.Sw",(double)0.99);//-
		sw_top  = ptree.get("boundary.top.Sw",(double)0.01);//-
	}

	//! boundary condition type function (true = Dirichlet)
	template<typename I, typename X>
	BCType type (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);

		if( bcset == 1 ){
			if( param.isBottomBoundary(xglobal) or param.isTopBoundary(xglobal) )
				return ConvectionDiffusionBoundaryConditions::Dirichlet;
			else return ConvectionDiffusionBoundaryConditions::Neumann;
		}else if( bcset == 2 ){
			if( param.isBottomInletBoundary(xglobal) )
				return ConvectionDiffusionBoundaryConditions::Dirichlet;
			else return ConvectionDiffusionBoundaryConditions::Neumann;
		}else{//default
			if( param.isBottomInletBoundary(xglobal) or param.isTopBoundary(xglobal) )
				return ConvectionDiffusionBoundaryConditions::Dirichlet;
			else return ConvectionDiffusionBoundaryConditions::Neumann;
		}

	}

	//! Dirichlet value
	//! Specifies Sw|_D
	template<typename I, typename X>
	double dirichlet_value (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);

		if( bcset == 1 ){
			if( param.isBottomInletBoundary(xglobal) )
				return sw_inB;
			else return sw_top;
		}else if( bcset == 2 ){
			return sw_inB;
		}else{//default
			if( param.isBottomInletBoundary(xglobal) ) return sw_inB;
			else return sw_top;
		}
	}

	//! Neumann boundary condition
	template<typename I, typename X>
	FieldVector
	neumann_value (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);
		FieldVector jF(0.0);
		return jF;
	}

	//! Outflow boundary condition
	template<typename I, typename X>
	FieldVector
	outflow_value (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);
		FieldVector oF(0.0);
		return oF;
	}
};

template<typename GV,typename PTree,typename Parameters>
class FluidizedSoilBoundaryCondition
{
	typedef ConvectionDiffusionBoundaryConditions::Type BCType;

private:
	const GV& gv;
	const PTree& ptree;
	const Parameters& param;
	const static int dim = GV::dimension;
	using FieldVector = typename Dune::FieldVector<double,dim>;

	double x_ast;
	double t_ast;
	double p_ast;

	int bcset;

	double cf_top;
	double cf_inB;

public:
	//! construct from grid view
	FluidizedSoilBoundaryCondition ( const GV& gv_ ,
						  	  	  	   const PTree& ptree_,
									   const Parameters& param_)
	: gv( gv_ ) , ptree(ptree_), param(param_)
	{
		//characteristic values
		x_ast  = ptree.get("characteristic_values.length"	,(double)1.0);//m
		t_ast  = ptree.get("characteristic_values.time"		,(double)1.0);//s
		p_ast  = ptree.get("characteristic_values.pressure"	,(double)1.0);//Pa

		bcset = ptree.get("boundary.case",(int) 0 );//-

		cf_top = ptree.get("boundary.top.Cf",(double)0.);//?
		cf_inB  = ptree.get("boundary.bottom_inlet.Cf",(double)0.0);//?
	}

	//! boundary condition type function (true = Dirichlet)
	template<typename I, typename X>
	BCType type (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);
		if( param.isTopBoundary(xglobal) )
			return ConvectionDiffusionBoundaryConditions::Dirichlet;
		else return ConvectionDiffusionBoundaryConditions::Neumann;
	}

	//! Dirichlet value
	//! Specifies Cf|_D
	template<typename I, typename X>
	double dirichlet_value (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);
		return cf_top;
	}

	//! Neumann boundary condition
	template<typename I, typename X>
	FieldVector
	neumann_value (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);
		FieldVector jF(0.0);
		return jF;
	}

	//! Outflow boundary condition
	template<typename I, typename X>
	FieldVector
	outflow_value (const I& i, const X& x, double& t) const {
		auto xglobal = i.geometry().global(x);
		FieldVector oF(0.0);
		return oF;
	}
};

#endif /* PROBLEM09_PARAMETERS_TEST_2PFLOW_PARAMETERS_HH_ */
