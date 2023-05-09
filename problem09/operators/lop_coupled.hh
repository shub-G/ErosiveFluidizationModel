/*
 * lop_coupled.hh
 *
 *  Created on: Apr 7, 2021
 *      Author: sgupta
 */

#ifndef PROBLEM09_OPERATORS_LOP_COUPLED_HH_
#define PROBLEM09_OPERATORS_LOP_COUPLED_HH_

template <class GV,class Params,class GFS,class U,class WBC,class NBC,class FSBC>
class LocalOperatorCoupled :
  public Dune::PDELab::NumericalJacobianApplyVolume		< LocalOperatorCoupled<GV,Params,GFS,U,WBC,NBC,FSBC> >,
  public Dune::PDELab::NumericalJacobianVolume			< LocalOperatorCoupled<GV,Params,GFS,U,WBC,NBC,FSBC> >,
  public Dune::PDELab::NumericalJacobianApplySkeleton	< LocalOperatorCoupled<GV,Params,GFS,U,WBC,NBC,FSBC> >,
  public Dune::PDELab::NumericalJacobianSkeleton		< LocalOperatorCoupled<GV,Params,GFS,U,WBC,NBC,FSBC> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary	< LocalOperatorCoupled<GV,Params,GFS,U,WBC,NBC,FSBC> >,
  public Dune::PDELab::NumericalJacobianBoundary		< LocalOperatorCoupled<GV,Params,GFS,U,WBC,NBC,FSBC> >,
  public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
    const GV& gv;
    const Params& param;
    GFS gfs;
    U *u_old;
	const WBC& w_bc;
	const NBC& n_bc;
	const FSBC& fs_bc;
	double *time;
	double *dt;
	unsigned int intorder;
	double epsilon;

public:
	  // pattern assembly flags
	  enum { doPatternVolume	= true };
	  enum { doPatternSkeleton	= true };

	  // residual assembly flags
	  enum { doAlphaVolume  	= true };
	  enum { doAlphaSkeleton	= true };
	  enum { doAlphaBoundary	= true };

	  using LFS = Dune::PDELab::LocalFunctionSpace<GFS>;
	  using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;
	  using VectorView = typename U::template LocalView<LFSCache>;
	  using RF = typename LFS::template Child<PrimaryVariables::Pw_q1>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
	  using RangeType = typename LFS::template Child<PrimaryVariables::Pw_q1>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
	  using JacobianType = typename LFS::template Child<PrimaryVariables::Pw_q1>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;
	  using size_type = typename LFS::template Child<PrimaryVariables::Pw_q1>::Type::Traits::SizeType;

	  // constructor stores parameters
	  LocalOperatorCoupled( const GV&	gv_,
						 	const Params& param_,
							GFS gfs_,
							U *u_old_,
							const WBC& w_bc_,
							const NBC& n_bc_,
							const FSBC& fs_bc_,
							double *time_,
							double *dt_,
							unsigned int intorder_ = 4,
							double 	epsilon_ = 1.e-6)
	  : gv(gv_),
		param(param_),
		gfs(gfs_),
		u_old(u_old_),
		w_bc(w_bc_),
		n_bc(n_bc_),
		fs_bc(fs_bc_),
		time(time_),
		dt(dt_),
		intorder( intorder_ ),
		epsilon( epsilon_ )
	  {}

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
		  // Reference to cell
		  const auto& cell = eg.entity();

		  // get geometry
		  auto geo = eg.geometry();

		  // dimension
		  const auto dim = gv.dimension;

		  // cell geometry
		  auto ref_el = referenceElement(geo);
		  auto cell_center_local = ref_el.position(0,0);
		  auto cell_center_global = cell.geometry().global(cell_center_local);
		  auto cell_volume = geo.volume();

		  // compute PVs at local center
		  //-----------------------------------------
		  // FV primary variables Pw, Sw, porosity, cf
		  //-----------------------------------------
		  double pw = x(lfsu.template child<PrimaryVariables::Pw>(),0);
		  double sw = x(lfsu.template child<PrimaryVariables::Sw>(),0);
		  double por= x(lfsu.template child<PrimaryVariables::porosity>(),0);
		  double cf = x(lfsu.template child<PrimaryVariables::Cf>(),0);
		  //-----------------------------------------
		  // primary variables Pw_old, Sw_old,
		  //-----------------------------------------
		  LFS lfs(gfs);
		  LFSCache lfscache(lfs);
		  VectorView u_old_view((*u_old));
		  lfs.bind(cell);
		  lfscache.update();
		  u_old_view.bind(lfscache);
		  std::vector<double> ul_old(lfs.size());
		  u_old_view.read( ul_old );
		  double pw_old = ul_old[ lfs.template child<PrimaryVariables::Pw>().localIndex(0)];
		  double sw_old = ul_old[ lfs.template child<PrimaryVariables::Sw>().localIndex(0)];
		  double por_old= ul_old[ lfs.template child<PrimaryVariables::porosity>().localIndex(0)];
		  double cf_old = ul_old[ lfs.template child<PrimaryVariables::Cf>().localIndex(0)];

		  //-----------------------------------------
		  // PROPERTIES
		  //-----------------------------------------
		  auto pc = param.CapillaryPressure(cell_center_global,sw,por);
		  auto pc_old = param.CapillaryPressure(cell_center_global,sw_old,por_old);
		  auto Cw = param.WettingPhaseCompressibility();
		  auto Cn = param.NonwettingPhaseCompressibility();
		  auto Cs = param.SoilPhaseCompressibility();
          auto rhow = param.WettingPhaseDensity(pw);
          auto rhon = param.NonwettingPhaseDensity(pw+pc);
          auto rhos = param.SoilPhaseDensity(pw+(1.-sw)*pc);
          auto muw = param.WettingPhaseViscosity();
          auto mun = param.NonwettingPhaseViscosity();
          auto K = param.PermeabilityTensor(cell_center_global,por);
          auto krw = param.WettingRelativePermeability(cell_center_global,sw);
		  auto krn = param.NonwettingRelativePermeability(cell_center_global,sw);
		  auto g = param.g() ;

		  //-----------------------------------------
		  // SOIL EROSION TERM
		  //-----------------------------------------
		  // Transformation matrix
		  typename EG::Geometry::JacobianInverseTransposed jac;
		  jac = geo.jacobianInverseTransposed(cell_center_local);

		  //evaluate vw_mag
		  std::vector<JacobianType> js_pw(lfsu.template child<PrimaryVariables::Pw_q1>().size());
		  lfsu.template child<PrimaryVariables::Pw_q1>().finiteElement().localBasis().evaluateJacobian(cell_center_local,js_pw);
		  std::vector<Dune::FieldVector<double,dim> > gradphi_pw(lfsu.template child<PrimaryVariables::Pw_q1>().size());
          for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pw_q1>().size(); i++){
              gradphi_pw[i] = 0.0;
              jac.umv(js_pw[i][0],gradphi_pw[i]);
          }
          Dune::FieldVector<double,dim> grad_pw(0.0);
          for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pw_q1>().size(); i++)
            grad_pw.axpy(x(lfsu.template child<PrimaryVariables::Pw_q1>(),i),gradphi_pw[i]);
          //KgradPw
		  Dune::FieldVector<double,dim> Kgrad_pw(0.0);
		  K.mv(grad_pw, Kgrad_pw);
		  //K*rho*g
		  Dune::FieldVector<double,dim> KFw(0.0);
		  auto Fw = g;
		  Fw *= rhow;
		  K.mv(Fw, KFw);
		  // magnitude of vw
		  double vw_mag_sq = 0.;
		  for(int i=0; i<dim; i++){
			  vw_mag_sq += ( (-1.)*(krw/muw)*(Kgrad_pw[i]+KFw[i]) )*( (-1.)*(krw/muw)*(Kgrad_pw[i]+KFw[i]) );
		  }
		  auto vw_mag = sqrt(vw_mag_sq);

		  //evaluate vn_mag
		  std::vector<JacobianType> js_pn(lfsu.template child<PrimaryVariables::Pn_q1>().size());
		  lfsu.template child<PrimaryVariables::Pn_q1>().finiteElement().localBasis().evaluateJacobian(cell_center_local,js_pn);
		  std::vector<Dune::FieldVector<double,dim> > gradphi_pn(lfsu.template child<PrimaryVariables::Pn_q1>().size());
          for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pn_q1>().size(); i++){
              gradphi_pn[i] = 0.0;
              jac.umv(js_pn[i][0],gradphi_pn[i]);
          }
          Dune::FieldVector<double,dim> grad_pn(0.0);
          for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pn_q1>().size(); i++)
            grad_pn.axpy(x(lfsu.template child<PrimaryVariables::Pn_q1>(),i),gradphi_pn[i]);
		  // KgradPn
		  Dune::FieldVector<double,dim> Kgrad_pn(0.0);
		  K.mv(grad_pn, Kgrad_pn);
		  Dune::FieldVector<double,dim> KFn(0.0);
		  auto Fn = g;
		  Fn *= rhon;
		  K.mv(Fn, KFn);
		  double vn_mag_sq = 0.;
		  for(int i=0; i<dim; i++){
			  vn_mag_sq += ( (-1.)*(krn/mun)*(Kgrad_pn[i]+KFn[i]) )
					      *( (-1.)*(krn/mun)*(Kgrad_pn[i]+KFn[i]) );
		  }
		  auto vn_mag = sqrt(vn_mag_sq);

		  //-----------------------------------------
		  // SOURCE TERMS
		  auto qw /*ndim*/ = param.WettingVolumetricSource(cell_center_global);
		  auto qn /*ndim*/ = param.NonwettingVolumetricSource(cell_center_global);
		  auto q_er /*ndim*/ = param.ErosionRate(cell_center_global,vw_mag,vn_mag,por);
		  //-----------------------------------------

		  // RESIDUALS
		  //-----------------------------------------
		  // wetting-phase mass balance
		  //-----------------------------------------
		  double term_w = (sw*por - sw_old*por_old)/(*dt)
						+ por * sw * Cw * (pw-pw_old)/(*dt)
						- qw;
		  r.accumulate(lfsu.template child<PrimaryVariables::Pw>(), 0., term_w*cell_volume );

		  //-----------------------------------------
		  // nonwetting-phase mass balance
		  //-----------------------------------------
		  double term_n = ( (1.-sw)*por - (1.-sw_old)*por_old)/(*dt)
						+ por * (1.-sw) * Cn * ( (pw-pw_old)/(*dt) + (pc-pc_old)/(*dt) )
						- qn;
		  r.accumulate(lfsu.template child<PrimaryVariables::Sw>(), 0., term_n*cell_volume );

		  //-----------------------------------------
		  // soil-phase mass balance
		  //-----------------------------------------
		  double term_soil = (-1.)*(por-por_old)/(*dt)
				  	  	   + (1.-por)*Cs*(pw-pw_old)/(*dt)
						   + (1.-por)*Cs*((1.-sw)*pc-(1.-sw_old)*pc_old)/(*dt)
						   - q_er;
		  r.accumulate(lfsu.template child<PrimaryVariables::porosity>(), 0., term_soil*cell_volume );

		  //-----------------------------------------
		  // fluidized-soil mass balance
		  //-----------------------------------------
		  double term_fs = (cf*sw*por-cf_old*sw_old*por_old)/(*dt)
						 + cf*sw*por*Cw*(pw-pw_old)/(*dt)
						 + q_er*(rhos/rhow);
		  r.accumulate(lfsu.template child<PrimaryVariables::Cf>(), 0., term_fs*cell_volume );

		  //-----------------------------------------
		  // L2-projection
		  //-----------------------------------------
		  // loop over quadrature points
		  for (const auto& ip : quadratureRule(geo,intorder))
		  {
			  auto ip_global = geo.global(ip.position());

			  // evaluate basis functions
			  std::vector<RangeType> phi_pw(lfsu.template child<PrimaryVariables::Pw_q1>().size());
			  lfsu.template child<PrimaryVariables::Pw_q1>().finiteElement().localBasis().evaluateFunction(ip.position(),phi_pw);

			  std::vector<RangeType> psi_pw(lfsv.template child<PrimaryVariables::Pw_q1>().size());
			  lfsv.template child<PrimaryVariables::Pw_q1>().finiteElement().localBasis().evaluateFunction(ip.position(),psi_pw);

			  std::vector<RangeType> phi_pn(lfsu.template child<PrimaryVariables::Pn_q1>().size());
			  lfsu.template child<PrimaryVariables::Pn_q1>().finiteElement().localBasis().evaluateFunction(ip.position(),phi_pn);

			  std::vector<RangeType> psi_pn(lfsv.template child<PrimaryVariables::Pn_q1>().size());
			  lfsv.template child<PrimaryVariables::Pn_q1>().finiteElement().localBasis().evaluateFunction(ip.position(),psi_pn);

			  // evaluate u (u->phase pressures)
			  RF u_pw=0.0;
			  for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pw_q1>().size(); i++)
				  u_pw += x(lfsu.template child<PrimaryVariables::Pw_q1>(),i)*phi_pw[i];

			  RF u_pn=0.0;
			  for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pn_q1>().size(); i++)
				  u_pn += x(lfsu.template child<PrimaryVariables::Pn_q1>(),i)*phi_pn[i];

			  // evaluate gradient of basis functions
			  std::vector<JacobianType> jsu_pw(lfsu.template child<PrimaryVariables::Pw_q1>().size());
			  lfsu.template child<PrimaryVariables::Pw_q1>().finiteElement().localBasis().evaluateJacobian(ip.position(),jsu_pw);

			  std::vector<JacobianType> jsv_pw(lfsv.template child<PrimaryVariables::Pw_q1>().size());
			  lfsv.template child<PrimaryVariables::Pw_q1>().finiteElement().localBasis().evaluateJacobian(ip.position(),jsv_pw);

			  std::vector<JacobianType> jsu_pn(lfsu.template child<PrimaryVariables::Pn_q1>().size());
			  lfsu.template child<PrimaryVariables::Pn_q1>().finiteElement().localBasis().evaluateJacobian(ip.position(),jsu_pn);

			  std::vector<JacobianType> jsv_pn(lfsv.template child<PrimaryVariables::Pn_q1>().size());
			  lfsv.template child<PrimaryVariables::Pn_q1>().finiteElement().localBasis().evaluateJacobian(ip.position(),jsv_pn);

			  // transform gradients of shape functions to real element
			  jac = geo.jacobianInverseTransposed(ip.position());

			  // evaluade gradients of shape fncs.
			  std::vector<Dune::FieldVector<RF,dim> > gradphi_pw(lfsu.template child<PrimaryVariables::Pw_q1>().size());
			  for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pw_q1>().size(); i++){
				  gradphi_pw[i] = 0.0;
				  jac.umv(jsu_pw[i][0],gradphi_pw[i]);
			  }

			  std::vector<Dune::FieldVector<RF,dim> > gradpsi_pw(lfsv.template child<PrimaryVariables::Pw_q1>().size());
			  for (unsigned int i=0; i<lfsv.template child<PrimaryVariables::Pw_q1>().size(); i++){
				  gradpsi_pw[i] = 0.0;
				  jac.umv(jsv_pw[i][0],gradpsi_pw[i]);
			  }

			  std::vector<Dune::FieldVector<RF,dim> > gradphi_pn(lfsu.template child<PrimaryVariables::Pn_q1>().size());
			  for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pn_q1>().size(); i++){
				  gradphi_pn[i] = 0.0;
				  jac.umv(jsu_pn[i][0],gradphi_pn[i]);
			  }

			  std::vector<Dune::FieldVector<RF,dim> > gradpsi_pn(lfsv.template child<PrimaryVariables::Pn_q1>().size());
			  for (unsigned int i=0; i<lfsv.template child<PrimaryVariables::Pn_q1>().size(); i++){
				  gradpsi_pn[i] = 0.0;
				  jac.umv(jsv_pn[i][0],gradpsi_pn[i]);
			  }

			  // compute gradient of u_pw, u_pn
			  Dune::FieldVector<RF,dim> gradu_pw(0.0);
			  for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pw_q1>().size(); i++)
				  gradu_pw.axpy(x(lfsu.template child<PrimaryVariables::Pw_q1>(),i),gradphi_pw[i]);

			  Dune::FieldVector<RF,dim> gradu_pn(0.0);
			  for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pn_q1>().size(); i++)
				  gradu_pn.axpy(x(lfsu.template child<PrimaryVariables::Pn_q1>(),i),gradphi_pn[i]);

			  // integrate
			  // FV --> FEM
			  //|| u_FV - u_FE || --> Min
			  RF factor = ip.weight() * geo.integrationElement(ip.position());
			  double tmp = 0.;
			  for (int j=0; j<lfsv.template child<PrimaryVariables::Pw_q1>().size(); j++){
				  tmp = epsilon*( gradu_pw*gradpsi_pw[j] ) + u_pw*psi_pw[j] - pw*psi_pw[j] ;
				  r.accumulate(lfsv.template child<PrimaryVariables::Pw_q1>(),j,tmp*factor );
			  }

			  tmp = 0.;
			  for (int j=0; j<lfsv.template child<PrimaryVariables::Pn_q1>().size(); j++){
				  tmp = epsilon*( gradu_pn*gradpsi_pn[j] ) + u_pn*psi_pn[j] - (pw+pc)*psi_pn[j] ;
				  r.accumulate(lfsv.template child<PrimaryVariables::Pn_q1>(),j,tmp*factor );
			  }

		  }//End Quadrature Rule

	  }//END:alpha_volume


	  // skeleton integral depending on test and ansatz functions
	  // each face is only visited ONCE!
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_skeleton (const IG& ig,
			  	  	  	   const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
						   const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
						   R& r_s, R& r_n) const
	  {
	        // References to inside and outside cells
	        const auto& cell_inside  = ig.inside();
	        const auto& cell_outside = ig.outside();

	        // get geometries
	        auto geo = ig.geometry();
	        const auto dim = gv.dimension;
	        auto geo_inside  = cell_inside.geometry();
	        auto geo_outside = cell_outside.geometry();

	        // cell geometries
	        auto ref_el_inside 	= referenceElement(geo_inside);
	        auto ref_el_outside = referenceElement(geo_outside);
	        auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	        auto outside_cell_center_local 	= ref_el_outside.position(0,0);
	        auto inside_cell_center_global 	= geo_inside.center();
	        auto outside_cell_center_global = geo_outside.center();

	        // distance of cell centers
	        auto d = outside_cell_center_global;
	        d -= inside_cell_center_global;
	        auto distance = d.two_norm();

	        // face geometry
	        auto ref_el = referenceElement(geo);
	        auto face_center_local = ref_el.position(0,0);
	        auto face_center_global = geo.center();
	        auto face_volume = geo.volume();
	        auto normal = ig.unitOuterNormal(face_center_local);

			// compute primary vars at local self and neighbour centers
	        //-------------------------------------------------
	        // pw_s, pw_n, sw_s, sw_n, por_s, por_n, cf_s, cf_n
	        //-------------------------------------------------
	        double pw_s = x_s(lfsu_s.template child<PrimaryVariables::Pw>(),0);
	        double sw_s = x_s(lfsu_s.template child<PrimaryVariables::Sw>(),0);
	        double por_s= x_s(lfsu_s.template child<PrimaryVariables::porosity>(),0);
	        double cf_s = x_s(lfsu_s.template child<PrimaryVariables::Cf>(),0);
	        double pw_n = x_n(lfsu_n.template child<PrimaryVariables::Pw>(),0);
	        double sw_n = x_n(lfsu_n.template child<PrimaryVariables::Sw>(),0);
	        double por_n= x_n(lfsu_n.template child<PrimaryVariables::porosity>(),0);
	        double cf_n = x_n(lfsu_n.template child<PrimaryVariables::Cf>(),0);

	        // evaluate properties
	        // at self center
	        auto pc_s = param.CapillaryPressure(inside_cell_center_global,sw_s,por_s);
	        auto rhow_s = param.WettingPhaseDensity(pw_s);
	        auto rhon_s = param.NonwettingPhaseDensity(pw_s+pc_s);
	        auto muw_s = param.WettingPhaseViscosity();
	        auto mun_s = param.NonwettingPhaseViscosity();
	        auto K_s = std::abs(param.Permeability(inside_cell_center_global,por_s) * normal);
	        auto krw_s = param.WettingRelativePermeability(inside_cell_center_global,sw_s);
	        auto krn_s = param.NonwettingRelativePermeability(inside_cell_center_global,sw_s);
	        // at neighbour center
	        auto pc_n = param.CapillaryPressure(outside_cell_center_global,sw_n,por_n);
	        auto rhow_n = param.WettingPhaseDensity(pw_n);
	        auto rhon_n = param.NonwettingPhaseDensity(pw_n+pc_n);
	        auto muw_n = param.WettingPhaseViscosity();
	        auto mun_n = param.NonwettingPhaseViscosity();
	        auto K_n = std::abs(param.Permeability(outside_cell_center_global,por_n) * normal);
	        auto krw_n = param.WettingRelativePermeability(outside_cell_center_global,sw_n);
	        auto krn_n = param.NonwettingRelativePermeability(outside_cell_center_global,sw_n);
	        // at interface
	        double K_int = 0.;
			if(K_s*K_n>0.) K_int = 2.*K_s*K_n/(K_s+K_n);
			auto gravity = param.g() * ig.unitOuterNormal(face_center_local) ;

			// upwinding
			// wrt w-phase velocity
			auto potential_w = (pw_n - pw_s)/distance + 0.5*(rhow_s+rhow_n) * gravity ;
			double upw_s = 0., upw_n = 0.;
			if( potential_w>0.){
				upw_s = 0.;
				upw_n = 1.;
			}else{
				upw_s = 1.;
				upw_n = 0.;
			}
			//wrt n-phase velocity
			auto potential_n = ( (pw_n+pc_n) - (pw_s+pc_s) )/distance + 0.5*(rhon_s+rhon_n) * gravity ;
			double upn_s = 0., upn_n = 0.;
			if( potential_n>0.){
				upn_s = 0.;
				upn_n = 1.;
			}else{
				upn_s = 1.;
				upn_n = 0.;
			}

			// RESIDUALS
			//-----------------------------------------
			// wetting-phase mass balance
			//-----------------------------------------
			double term_w = (-1.)*K_int*( upw_s*krw_s/muw_s + upw_n*krw_n/muw_n ) * potential_w;
//			std::cout<< inside_cell_center_global[0] << "," << inside_cell_center_global[1] << '\t' << term_w << std::endl;
			r_s.accumulate(lfsu_s.template child<PrimaryVariables::Pw>() , 0,  term_w*face_volume);
			r_n.accumulate(lfsu_n.template child<PrimaryVariables::Pw>() , 0, -term_w*face_volume);
			//-----------------------------------------
			// nonwetting-phase mass balance
			//-----------------------------------------
			double term_n = (-1.)*K_int*( upn_s*krn_s/mun_s + upn_n*krn_n/mun_n ) * potential_n;
//			std::cout<< outside_cell_center_global[0] << "," << outside_cell_center_global[1] << '\t' << term_n << std::endl;
			r_s.accumulate(lfsu_s.template child<PrimaryVariables::Sw>() , 0,  term_n*face_volume);
			r_n.accumulate(lfsu_n.template child<PrimaryVariables::Sw>() , 0, -term_n*face_volume);
			//-----------------------------------------
			// fluidized-soil mass balance
			//-----------------------------------------
			double term_fs = (-1.)*K_int*( upw_s*cf_s*krw_s/muw_s + upw_n*cf_n*krw_n/muw_n ) * potential_w;
			r_s.accumulate(lfsu_s.template child<PrimaryVariables::Cf>() , 0,  term_fs*face_volume);
			r_n.accumulate(lfsu_n.template child<PrimaryVariables::Cf>() , 0, -term_fs*face_volume);

	  }//END:alpha_skeleton


	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary ( const IG& ig,
			  	  	  	  	const LFSU& lfsu, const X& x, const LFSV& lfsv,
							R& r ) const
	  {
	        // References to inside and outside cells
	        const auto& cell_inside  = ig.inside();

	        // get geometries
	        auto geo = ig.geometry();
	        const auto dim = gv.dimension;
	        auto geo_inside = cell_inside.geometry();

	        // cell geometries
	        auto ref_el_inside 	= referenceElement(geo_inside);
	        auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	        auto inside_cell_center_global 	= geo_inside.center();

	        // face geometry
	        auto ref_el = referenceElement(geo);
	        auto face_center_local = ref_el.position(0,0);
	        auto face_center_global = geo.center();
	        auto face_volume = geo.volume();
	        auto normal = ig.unitOuterNormal(face_center_local);

	        // distance of cell centers
	        auto d = geo.global(face_center_local);
	        d -= inside_cell_center_global;
	        auto distance = d.two_norm();

	        // boundary conditions
	        //-----------------------------------------
	        // Boundary types
	        //-----------------------------------------
	        auto t_new = (*time)+(*dt);
	        auto bct_w = w_bc.type( ig,face_center_local,t_new);
			auto bct_n = n_bc.type( ig,face_center_local,t_new);
			auto bct_fs= fs_bc.type(ig,face_center_local,t_new);
			//-----------------------------------------
			// dirichlet boundary values
			//-----------------------------------------
			auto pw_D = w_bc.dirichlet_value( ig,face_center_local,t_new);
			auto sw_D = n_bc.dirichlet_value( ig,face_center_local,t_new);
			auto cf_D = fs_bc.dirichlet_value(ig,face_center_local,t_new);
			//-----------------------------------------
			// neumann boundary values
			//-----------------------------------------
			auto jw_N = w_bc.neumann_value( ig,face_center_local,t_new);
			auto jn_N = n_bc.neumann_value( ig,face_center_local,t_new);
			auto jfs_N= fs_bc.neumann_value(ig,face_center_local,t_new);
			//-----------------------------------------
			// outflow boundary values
			//-----------------------------------------
			auto jw_O = w_bc.outflow_value( ig,face_center_local,t_new);
			auto jn_O = n_bc.outflow_value( ig,face_center_local,t_new);
			auto jfs_O= fs_bc.outflow_value(ig,face_center_local,t_new);

			// compute primary vars at local self and neighbour centers
	        //-----------------------------------------
	        // pw_s, pw_n, sw_s, sw_n
	        //-----------------------------------------
	        double pw_s = x(lfsu.template child<PrimaryVariables::Pw>(),0);
	        double sw_s = x(lfsu.template child<PrimaryVariables::Sw>(),0);
	        double por_s= x(lfsu.template child<PrimaryVariables::porosity>(),0);
	        double cf_s = x(lfsu.template child<PrimaryVariables::Cf>(),0);
	        double pw_n = pw_s;
	        if( bct_w==ConvectionDiffusionBoundaryConditions::Dirichlet ) pw_n = pw_D;
	        double sw_n = sw_s;
	        if( bct_n==ConvectionDiffusionBoundaryConditions::Dirichlet ) sw_n = sw_D;
	        double por_n = por_s;
	        double cf_n = cf_s;
	        if( bct_fs==ConvectionDiffusionBoundaryConditions::Dirichlet) cf_n = cf_D;

	        // evaluate properties
	        // at self center
	        auto pc_s = param.CapillaryPressure(inside_cell_center_global,sw_s,por_s);
	        auto rhow_s = param.WettingPhaseDensity(pw_s);
	        auto rhon_s = param.NonwettingPhaseDensity(pw_s+pc_s);
	        auto muw_s = param.WettingPhaseViscosity();
	        auto mun_s = param.NonwettingPhaseViscosity();
	        auto K_s = std::abs(param.Permeability(inside_cell_center_global,por_s) * normal);
	        auto krw_s = param.WettingRelativePermeability(inside_cell_center_global,sw_s);
	        auto krn_s = param.NonwettingRelativePermeability(inside_cell_center_global,sw_s);
	        // at neighbour center
	        auto pc_n = param.CapillaryPressure(face_center_global,sw_n,por_n);
	        auto rhow_n = param.WettingPhaseDensity(pw_n);
	        auto rhon_n = param.NonwettingPhaseDensity(pw_n+pc_n);
	        auto muw_n = param.WettingPhaseViscosity();
	        auto mun_n = param.NonwettingPhaseViscosity();
	        auto K_n = std::abs(param.Permeability(face_center_global,por_n) * normal);
	        auto krw_n = param.WettingRelativePermeability(face_center_global,sw_n);
	        auto krn_n = param.NonwettingRelativePermeability(face_center_global,sw_n);
	        // at interface
	        double K_int = 0.;
			if(K_s*K_n>0.) K_int = 2.*K_s*K_n/(K_s+K_n);
			auto gravity = param.g() * ig.unitOuterNormal(face_center_local) ;

			// RESIDUALS
			//-----------------------------------------
			// fluidized-soil mass balance (NEUMANN case)
			//-----------------------------------------
			if( bct_fs==ConvectionDiffusionBoundaryConditions::Neumann ){
				double term_fs = jfs_N * normal;
				r.accumulate(lfsu.template child<PrimaryVariables::Cf>() , 0,  term_fs*face_volume);
			}

			//-----------------------------------------
			// wetting-phase mass balance
			//-----------------------------------------
			if( bct_w==ConvectionDiffusionBoundaryConditions::Neumann ){
				double term_w = jw_N * normal;
				r.accumulate(lfsu.template child<PrimaryVariables::Pw>() , 0,  term_w*face_volume);

				// fluidized-soil mass balance (NEUMANN case)
				if( bct_fs==!ConvectionDiffusionBoundaryConditions::Neumann ){
					double upw_s = 0., upw_n = 0.;
					if( term_w>0.){
						upw_s = 0.;
						upw_n = 1.;
					}else{
						upw_s = 1.;
						upw_n = 0.;
					}
					double term_fs = (upw_s*cf_s+upw_n*cf_n)*term_w;
					if( bct_fs==ConvectionDiffusionBoundaryConditions::Outflow ){
						term_fs += jfs_O * normal;
					}
					r.accumulate(lfsu.template child<PrimaryVariables::Cf>() , 0,  term_fs*face_volume);
				}
			}else{
				// upwinding wrt w-phase velocity
				auto potential_w = (pw_n - pw_s)/distance + 0.5*(rhow_s+rhow_n) * gravity ;
				double upw_s = 0., upw_n = 0.;
				if( potential_w>0.){
					upw_s = 0.;
					upw_n = 1.;
				}else{
					upw_s = 1.;
					upw_n = 0.;
				}
				double term_w = (-1.)*K_int*( upw_s*krw_s/muw_s + upw_n*krw_n/muw_n ) * potential_w;
				if( bct_w==ConvectionDiffusionBoundaryConditions::Outflow ){
					term_w += jw_O * normal;
				}
				r.accumulate(lfsu.template child<PrimaryVariables::Pw>() , 0,  term_w*face_volume);

				double term_fs = (-1.)*K_int*( upw_s*cf_s*krw_s/muw_s + upw_n*cf_n*krw_n/muw_n ) * potential_w;
				if( bct_fs==ConvectionDiffusionBoundaryConditions::Outflow ){
					term_fs += jfs_O * normal;
				}
				r.accumulate(lfsu.template child<PrimaryVariables::Cf>() , 0,  term_fs*face_volume);
			}

			//-----------------------------------------
			// nonwetting-phase mass balance
			//-----------------------------------------
			if( bct_n==ConvectionDiffusionBoundaryConditions::Neumann ){
				double term_n = jn_N * normal;
				r.accumulate(lfsu.template child<PrimaryVariables::Sw>() , 0,  term_n*face_volume);
			}else{
				// upwinding wrt n-phase velocity
				auto potential_n = ( (pw_n+pc_n) - (pw_s+pc_s) )/distance + 0.5*(rhon_s+rhon_n) * gravity ;
				double upn_s = 0., upn_n = 0.;
				if( potential_n>0.){
					upn_s = 0.;
					upn_n = 1.;
				}else{
					upn_s = 1.;
					upn_n = 0.;
				}
				double term_n = (-1.)*K_int*( upn_s*krn_s/mun_s + upn_n*krn_n/mun_n ) * potential_n;
				if( bct_n==ConvectionDiffusionBoundaryConditions::Outflow ){
					term_n += jn_O * normal;
				}
				r.accumulate(lfsu.template child<PrimaryVariables::Sw>() , 0,  term_n*face_volume);
			}

	  }//END:alpha_boundary

};

#endif /* PROBLEM09_OPERATORS_LOP_COUPLED_HH_ */
