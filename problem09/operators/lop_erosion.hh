/*
 * lop_erosion.hh
 *
 *  Created on: Apr 6, 2021
 *      Author: sgupta
 */

#ifndef PROBLEM09_OPERATORS_LOP_EROSION_HH_
#define PROBLEM09_OPERATORS_LOP_EROSION_HH_

template <class GV,class Params,class GFSEr,class UEr,class GFS2P,class U2P,class GFSL2,class UL2,class FSBC,class WBC,class NBC>
class LocalOperatorErosion :
  public Dune::PDELab::NumericalJacobianApplyVolume		< LocalOperatorErosion<GV,Params,GFSEr,UEr,GFS2P,U2P,GFSL2,UL2,FSBC,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianVolume			< LocalOperatorErosion<GV,Params,GFSEr,UEr,GFS2P,U2P,GFSL2,UL2,FSBC,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianApplySkeleton	< LocalOperatorErosion<GV,Params,GFSEr,UEr,GFS2P,U2P,GFSL2,UL2,FSBC,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianSkeleton		< LocalOperatorErosion<GV,Params,GFSEr,UEr,GFS2P,U2P,GFSL2,UL2,FSBC,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary	< LocalOperatorErosion<GV,Params,GFSEr,UEr,GFS2P,U2P,GFSL2,UL2,FSBC,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianBoundary		< LocalOperatorErosion<GV,Params,GFSEr,UEr,GFS2P,U2P,GFSL2,UL2,FSBC,WBC,NBC> >,
  public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
    const GV& gv;
    const Params& param;
    GFS2P gfs_2p;
    U2P *u2p;
    U2P *u2p_dot;
    GFSEr gfs_er;
    UEr *uer_old;
    GFSL2 gfs_l2;
    UL2 *ul2;
	const WBC& w_bc;
	const NBC& n_bc;
	const FSBC& fs_bc;
	double *time;
	double *dt;
	unsigned int intorder ;

public:
	  // pattern assembly flags
	  enum { doPatternVolume	= true };
	  enum { doPatternSkeleton	= true };

	  // residual assembly flags
	  enum { doAlphaVolume  	= true };
	  enum { doAlphaSkeleton	= true };
	  enum { doAlphaBoundary	= true };

	  using LFS2P = Dune::PDELab::LocalFunctionSpace<GFS2P>;
	  using LFS2PCache = Dune::PDELab::LFSIndexCache<LFS2P>;
	  using VectorView2P = typename U2P::template LocalView<LFS2PCache>;

	  using LFSEr = Dune::PDELab::LocalFunctionSpace<GFSEr>;
	  using LFSErCache = Dune::PDELab::LFSIndexCache<LFSEr>;
	  using VectorViewEr = typename UEr::template LocalView<LFSErCache>;

	  using LFSL2 = Dune::PDELab::LocalFunctionSpace<GFSL2>;
	  using LFSL2Cache = Dune::PDELab::LFSIndexCache<LFSL2>;
	  using VectorViewL2 = typename UL2::template LocalView<LFSL2Cache>;
	  using JacobianType = typename LFSL2::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;
	  using size_type = typename LFSL2::template Child<0>::Type::Traits::SizeType;

	  // constructor stores parameters
	  LocalOperatorErosion( const GV&	gv_,
							const Params& param_,
							GFSEr gfs_er_,
							UEr *uer_old_,
							GFS2P gfs_2p_,
							U2P *u2p_,
							U2P *u2p_dot_,
							GFSL2 gfs_l2_,
							UL2 *ul2_,
							const FSBC& fs_bc_,
							const WBC& w_bc_,
							const NBC& n_bc_,
							double *time_,
							double *dt_	,
							unsigned int  intorder_ = 4)
	  : gv(gv_),
		param(param_),
		gfs_er(gfs_er_),
		uer_old(uer_old_),
		gfs_2p(gfs_2p_),
		u2p(u2p_),
		u2p_dot(u2p_dot_),
		gfs_l2(gfs_l2_),
		ul2(ul2_),
		fs_bc(fs_bc_),
		w_bc(w_bc_),
		n_bc(n_bc_),
		time(time_),
		dt(dt_),
		intorder( intorder_ )
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
		  // primary variables Pw and Sw
		  //-----------------------------------------
		  double porosity = x(lfsu.template child<PrimaryVariables::porosity>(),0);
		  double cf = x(lfsu.template child<PrimaryVariables::Cf>(),0);
		  //-----------------------------------------
		  // primary variables porosity_old and Cf_old
		  //-----------------------------------------
		  LFSEr lfser(gfs_er);
		  LFSErCache lfsercache(lfser);
		  VectorViewEr uer_old_view((*uer_old));
		  lfser.bind(cell);
		  lfsercache.update();
		  uer_old_view.bind(lfsercache);
		  std::vector<double> uler_old(lfser.size());
		  uer_old_view.read( uler_old );
		  double porosity_old = uler_old[ lfser.template child<PrimaryVariables::porosity>().localIndex(0)];
		  double cf_old = uler_old[ lfser.template child<PrimaryVariables::Cf>().localIndex(0)];
		  //-----------------------------------------
		  // variables Pw, Sw, dt_Pw, dt_Sw
		  //-----------------------------------------
		  LFS2P lfs2p(gfs_2p);
		  LFS2PCache lfs2pcache(lfs2p);
		  VectorView2P u2p_view((*u2p));
		  VectorView2P u2p_dot_view((*u2p_dot));
		  lfs2p.bind(cell);
		  lfs2pcache.update();
		  u2p_view.bind(lfs2pcache);
		  u2p_dot_view.bind(lfs2pcache);
		  std::vector<double> ul2p(lfs2p.size());
		  u2p_view.read( ul2p );
		  std::vector<double> ul2p_dot(lfs2p.size());
		  u2p_dot_view.read( ul2p_dot );
		  double pw = ul2p[ lfs2p.template child<PrimaryVariables::Pw>().localIndex(0)];
		  double sw = ul2p[ lfs2p.template child<PrimaryVariables::Sw>().localIndex(0)];
		  double dt_pw = ul2p_dot[ lfs2p.template child<PrimaryVariables::Pw>().localIndex(0)];
		  double dt_sw = ul2p_dot[ lfs2p.template child<PrimaryVariables::Sw>().localIndex(0)];
		  //-----------------------------------------
		  // variables Pw_l2, Pn_l2
		  //-----------------------------------------
		  LFSL2 lfsl2(gfs_l2);
		  LFSL2Cache lfsl2cache(lfsl2);
		  VectorViewL2 ul2_view((*ul2));
		  lfsl2.bind(cell);
		  lfsl2cache.update();
		  ul2_view.bind(lfsl2cache);
		  std::vector<double> ull2(lfsl2.size());

		  //-----------------------------------------
		  // SOIL EROSION TERM
		  //-----------------------------------------
		  //evaluate vw, vn
		  typename EG::Geometry::JacobianInverseTransposed jac;
		  jac = geo.jacobianInverseTransposed(cell_center_local);

		  std::vector<JacobianType> js_pw(lfsl2.template child<PrimaryVariables::Pw_l2>().size());
		  lfsl2.template child<PrimaryVariables::Pw_l2>().finiteElement().localBasis().evaluateJacobian(cell_center_local,js_pw);
		  std::vector<Dune::FieldVector<double,dim> > gradphi_pw(lfsl2.template child<PrimaryVariables::Pw_l2>().size());
          for (unsigned int i=0; i<lfsl2.template child<PrimaryVariables::Pw_l2>().size(); i++){
              gradphi_pw[i] = 0.0;
              jac.umv(js_pw[i][0],gradphi_pw[i]);
          }
          Dune::FieldVector<double,dim> grad_pw(0.0);
          for (unsigned int i=0; i<lfsl2.template child<PrimaryVariables::Pw_l2>().size(); i++)
            grad_pw.axpy(ull2[lfsl2.template child<PrimaryVariables::Pw_l2>().localIndex(i)],gradphi_pw[i]);

		  std::vector<JacobianType> js_pn(lfsl2.template child<PrimaryVariables::Pn_l2>().size());
		  lfsl2.template child<PrimaryVariables::Pn_l2>().finiteElement().localBasis().evaluateJacobian(cell_center_local,js_pn);
		  std::vector<Dune::FieldVector<double,dim> > gradphi_pn(lfsl2.template child<PrimaryVariables::Pn_l2>().size());
          for (unsigned int i=0; i<lfsl2.template child<PrimaryVariables::Pn_l2>().size(); i++){
              gradphi_pn[i] = 0.0;
              jac.umv(js_pn[i][0],gradphi_pn[i]);
          }
          Dune::FieldVector<double,dim> grad_pn(0.0);
          for (unsigned int i=0; i<lfsl2.template child<PrimaryVariables::Pn_l2>().size(); i++)
            grad_pn.axpy(ull2[lfsl2.template child<PrimaryVariables::Pn_l2>().localIndex(i)],gradphi_pn[i]);

          auto pc = param.CapillaryPressure(cell_center_global,sw,porosity);
          auto rhow = param.WettingPhaseDensity(pw,cf);
          auto rhon = param.NonwettingPhaseDensity(pw+pc);
          auto rhos = param.SoilPhaseDensity(pw+(1.-sw)*pc);
          auto muw = param.WettingPhaseViscosity();
          auto mun = param.NonwettingPhaseViscosity();
          auto K = param.PermeabilityTensor(cell_center_global,porosity);
          auto krw = param.WettingRelativePermeability(cell_center_global,sw);
		  auto krn = param.NonwettingRelativePermeability(cell_center_global,sw);
		  auto g = param.g() ;
		  // KgradP
		  Dune::FieldVector<double,dim> Kgrad_pw(0.0);
		  K.mv(grad_pw, Kgrad_pw);
		  // KgradP
		  Dune::FieldVector<double,dim> Kgrad_pn(0.0);
		  K.mv(grad_pn, Kgrad_pn);
		  // K*rho*g
		  Dune::FieldVector<double,dim> KFw(0.0);
		  auto Fw = g;
		  Fw *= rhow;
		  K.mv(Fw, KFw);
		  // K*rho*g
		  Dune::FieldVector<double,dim> KFn(0.0);
		  auto Fn = g;
		  Fn *= rhon;
		  K.mv(Fn, KFn);

		  // magnitudes of vw and vn
		  double vw_mag_sq = 0.;
		  for(int i=0; i<dim; i++){
			  vw_mag_sq += ( (-1.)*(krw/muw)*(Kgrad_pw[i]+KFw[i]) )*( (-1.)*(krw/muw)*(Kgrad_pw[i]+KFw[i]) );
		  }
		  auto vw_mag = sqrt(vw_mag_sq);

		  double vn_mag_sq = 0.;
		  for(int i=0; i<dim; i++){
			  vn_mag_sq += ( (-1.)*(krn/mun)*(Kgrad_pn[i]+KFn[i]) )*( (-1.)*(krn/mun)*(Kgrad_pn[i]+KFn[i]) );
		  }
		  auto vn_mag = sqrt(vn_mag_sq);
		  //-----------------------------------------
		  // compute erosion rate
		  auto q_er /*ndim*/ = param.ErosionRate(cell_center_global,vw_mag,vn_mag,porosity);
		  // compute deposition rate
		  auto q_sed/*ndim*/ = param.DepositionRate(cell_center_global,cf,sw,porosity);
		  //-----------------------------------------

		  // phase compressibilities
		  auto Cw = param.WettingPhaseCompressibility();
		  auto Cn = param.NonwettingPhaseCompressibility();
		  auto Cs = param.SoilPhaseCompressibility();

		  // RESIDUALS
		  //-----------------------------------------
		  // soil-phase mass balance
		  //-----------------------------------------
		  auto dpc_dpor= param.dPc_dPor(cell_center_global,sw,porosity);
		  auto dpc_dsw = param.dPc_dSw( cell_center_global,sw,porosity);
		  double term_soil = ( Cs*(1.-porosity)*(1.-sw)*dpc_dpor - 1. ) * (porosity-porosity_old)/(*dt)
				  	  	   + (1.-porosity)*Cs*dt_pw
						   + Cs*(1.-porosity)*((1.-sw)*dpc_dsw - pc ) * dt_sw
						   - q_er
						   - q_sed;
		  r.accumulate(lfsu.template child<PrimaryVariables::porosity>(), 0., term_soil*cell_volume );

		  //-----------------------------------------
		  // fluidized-soil mass balance
		  //-----------------------------------------
		  double term_fs = sw*(cf*porosity-cf_old*porosity_old)/(*dt)
						 + porosity*cf*dt_sw
						 + porosity*cf*sw*Cw*dt_pw
						 + q_er*(rhos/rhow)
						 + q_sed*(rhos/rhow);
		  r.accumulate(lfsu.template child<PrimaryVariables::Cf>(), 0., term_fs*cell_volume );

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
	        //-----------------------------------------
	        // por_s, por_n, cf_s, cf_n
	        //-----------------------------------------
	        double por_s = x_s(lfsu_s.template child<PrimaryVariables::porosity>(),0);
	        double cf_s  = x_s(lfsu_s.template child<PrimaryVariables::Cf>(),0);
	        double por_n = x_n(lfsu_n.template child<PrimaryVariables::porosity>(),0);
	        double cf_n  = x_n(lfsu_n.template child<PrimaryVariables::Cf>(),0);
	        //-----------------------------------------
	        // variables Pw_s, Sw_s, Pw_n, Sw_n
	        //-----------------------------------------
	        LFS2P lfs2p_s(gfs_2p);
			LFS2PCache lfs2pcache_s(lfs2p_s);
			VectorView2P u2p_view_s((*u2p));
			VectorView2P u2p_dot_view_s((*u2p_dot));
			lfs2p_s.bind(cell_inside);
			lfs2pcache_s.update();
			u2p_view_s.bind(lfs2pcache_s);
			std::vector<double> ul2p_s(lfs2p_s.size());
			u2p_view_s.read( ul2p_s );
			double pw_s = ul2p_s[ lfs2p_s.template child<PrimaryVariables::Pw>().localIndex(0)];
			double sw_s = ul2p_s[ lfs2p_s.template child<PrimaryVariables::Sw>().localIndex(0)];

	        LFS2P lfs2p_n(gfs_2p);
			LFS2PCache lfs2pcache_n(lfs2p_n);
			VectorView2P u2p_view_n((*u2p));
			VectorView2P u2p_dot_view_n((*u2p_dot));
			lfs2p_n.bind(cell_outside);
			lfs2pcache_n.update();
			u2p_view_n.bind(lfs2pcache_n);
			std::vector<double> ul2p_n(lfs2p_n.size());
			u2p_view_n.read( ul2p_n );
			double pw_n = ul2p_n[ lfs2p_n.template child<PrimaryVariables::Pw>().localIndex(0)];
			double sw_n = ul2p_n[ lfs2p_n.template child<PrimaryVariables::Sw>().localIndex(0)];


	        // evaluate properties
	        // at self center
	        auto rhow_s = param.WettingPhaseDensity(pw_s,cf_s);
	        auto muw_s = param.WettingPhaseViscosity();
	        auto K_s = std::abs(param.Permeability(inside_cell_center_global,por_s) * normal);
	        auto krw_s = param.WettingRelativePermeability(inside_cell_center_global,sw_s);
	        // at neighbour center
	        auto rhow_n = param.WettingPhaseDensity(pw_n,cf_n);
	        auto muw_n = param.WettingPhaseViscosity();
	        auto K_n = std::abs(param.Permeability(outside_cell_center_global,por_n) * normal);
	        auto krw_n = param.WettingRelativePermeability(outside_cell_center_global,sw_n);
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

			// RESIDUALS
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

			// compute primary vars at local self and neighbour centers
	        //-----------------------------------------
	        // por_s, por_n, cf_s, cf_n
	        //-----------------------------------------
	        double por_s= x(lfsu.template child<PrimaryVariables::porosity>(),0);
	        double por_n = por_s;
	        double cf_s = x(lfsu.template child<PrimaryVariables::Cf>(),0);
	        double cf_n = cf_s;
				        //-----------------------------------------
	        // variables Pw_s, Sw_s, Pw_n, Sw_n
	        //-----------------------------------------
	        LFS2P lfs2p_s(gfs_2p);
			LFS2PCache lfs2pcache_s(lfs2p_s);
			VectorView2P u2p_view_s((*u2p));
			VectorView2P u2p_dot_view_s((*u2p_dot));
			lfs2p_s.bind(cell_inside);
			lfs2pcache_s.update();
			u2p_view_s.bind(lfs2pcache_s);
			std::vector<double> ul2p_s(lfs2p_s.size());
			u2p_view_s.read( ul2p_s );
			double pw_s = ul2p_s[ lfs2p_s.template child<PrimaryVariables::Pw>().localIndex(0)];
			double pw_n = pw_s;
			double sw_s = ul2p_s[ lfs2p_s.template child<PrimaryVariables::Sw>().localIndex(0)];
			double sw_n = sw_s;

	        // boundary conditions
	        //-----------------------------------------
	        // Boundary types
	        //-----------------------------------------
	        auto t_new = (*time)+(*dt);
	        auto bct_w = w_bc.type(ig,face_center_local,t_new);
	        auto bct_n = n_bc.type(ig,face_center_local,t_new);
			auto bct_fs = fs_bc.type(ig,face_center_local,t_new);
			//-----------------------------------------
			// dirichlet boundary values
			//-----------------------------------------
			auto cf_D = fs_bc.dirichlet_value(ig,face_center_local,t_new);
			auto pw_D = w_bc.dirichlet_value(ig,face_center_local,t_new,cf_D);
			auto sw_D = n_bc.dirichlet_value(ig,face_center_local,t_new);
			if( bct_fs==ConvectionDiffusionBoundaryConditions::Dirichlet ) cf_n = cf_D;
			if( bct_w==ConvectionDiffusionBoundaryConditions::Dirichlet ) pw_n = pw_D;
			if( bct_n==ConvectionDiffusionBoundaryConditions::Dirichlet ) sw_n = sw_D;
			//-----------------------------------------
			// neumann boundary values
			//-----------------------------------------
			auto jw_N = w_bc.neumann_value(ig,face_center_local,t_new);
			auto jfs_N = fs_bc.neumann_value(ig,face_center_local,t_new);
			//-----------------------------------------
			// outflow boundary values
			//-----------------------------------------
			auto jw_O = w_bc.outflow_value(ig,face_center_local,t_new);

	        // evaluate properties
	        // at self center
	        auto rhow_s = param.WettingPhaseDensity(pw_s,cf_s);
	        auto muw_s = param.WettingPhaseViscosity();
	        auto K_s = std::abs(param.Permeability(inside_cell_center_global,por_s) * normal);
	        auto krw_s = param.WettingRelativePermeability(inside_cell_center_global,sw_s);
	        // at neighbour center
	        auto rhow_n = param.WettingPhaseDensity(pw_n,cf_n);
	        auto muw_n = param.WettingPhaseViscosity();
	        auto K_n = std::abs(param.Permeability(face_center_global,por_n) * normal);
	        auto krw_n = param.WettingRelativePermeability(face_center_global,sw_n);
	        // at interface
	        double K_int = 0.;
			if(K_s*K_n>0.) K_int = 2.*K_s*K_n/(K_s+K_n);
			auto gravity = param.g() * ig.unitOuterNormal(face_center_local) ;

			// RESIDUALS
			//-----------------------------------------
			// fluidized-soil mass balance
			//-----------------------------------------
			if( bct_w==ConvectionDiffusionBoundaryConditions::Neumann ){
				double term_w = jw_N * normal;
				double upw_s = 0., upw_n = 0.;
				if( term_w>0.){
					upw_s = 0.;
					upw_n = 1.;
				}else{
					upw_s = 1.;
					upw_n = 0.;
				}
				double term_fs = (upw_s*cf_s+upw_n*cf_n)*term_w;
//				if( bct_fs==ConvectionDiffusionBoundaryConditions::Neumann ){
//					term_fs = jfs_N * normal;
//				}
				r.accumulate(lfsu.template child<PrimaryVariables::Cf>() , 0,  term_fs*face_volume);
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
//				double term_fs = (-1.)*K_int*( upw_s*cf_s*krw_s/muw_s + upw_n*cf_n*krw_n/muw_n ) * potential_w;
				double term_fs = (-1.)*cf_n*K_int*( upw_s*krw_s/muw_s + upw_n*krw_n/muw_n ) * potential_w;
//				if( bct_w==ConvectionDiffusionBoundaryConditions::Outflow ){
//					term_fs += (upw_s*cf_s+upw_n*cf_n)*jw_O * normal;
//				}
//				if( bct_fs==ConvectionDiffusionBoundaryConditions::Neumann ){
//					term_fs = jfs_N * normal;
//				}
				r.accumulate(lfsu.template child<PrimaryVariables::Cf>() , 0,  term_fs*face_volume);
			}

	  }//END:alpha_boundary

};

#endif /* PROBLEM09_OPERATORS_LOP_EROSION_HH_ */
