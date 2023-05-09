/*
 * lop_2pflow.hh
 *
 *  Created on: Apr 6, 2021
 *      Author: sgupta
 */

#ifndef PROBLEM09_OPERATORS_LOP_2PFLOW_HH_
#define PROBLEM09_OPERATORS_LOP_2PFLOW_HH_

template <class GV,class Params,class GFS2P,class U2P,class GFSEr,class UEr,class WBC,class NBC>
class LocalOperatorTwoPhaseFlow :
  public Dune::PDELab::NumericalJacobianApplyVolume		< LocalOperatorTwoPhaseFlow<GV,Params,GFS2P,U2P,GFSEr,UEr,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianVolume			< LocalOperatorTwoPhaseFlow<GV,Params,GFS2P,U2P,GFSEr,UEr,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianApplySkeleton	< LocalOperatorTwoPhaseFlow<GV,Params,GFS2P,U2P,GFSEr,UEr,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianSkeleton		< LocalOperatorTwoPhaseFlow<GV,Params,GFS2P,U2P,GFSEr,UEr,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary	< LocalOperatorTwoPhaseFlow<GV,Params,GFS2P,U2P,GFSEr,UEr,WBC,NBC> >,
  public Dune::PDELab::NumericalJacobianBoundary		< LocalOperatorTwoPhaseFlow<GV,Params,GFS2P,U2P,GFSEr,UEr,WBC,NBC> >,
  public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
    const GV& gv;
    const Params& param;
    GFS2P gfs_2p;
    U2P *u2p_old;
    GFSEr gfs_er;
    UEr *uer;
    UEr *uer_dot;
	const WBC& w_bc;
	const NBC& n_bc;
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

	  // constructor stores parameters
	  LocalOperatorTwoPhaseFlow( const GV&	gv_,
			  	  	  	  	  	 const Params& param_,
			  	  	 	 	 	 GFS2P gfs_2p_,
								 U2P *u2p_old_,
								 GFSEr gfs_er_,
								 UEr *uer_,
								 UEr *uer_dot_,
								 const WBC& w_bc_	,
								 const NBC& n_bc_	,
								 double *time_	,
								 double *dt_	,
								 unsigned int  intorder_ = 4)
	  : gv(gv_),
		param(param_),
		gfs_2p(gfs_2p_),
		u2p_old(u2p_old_),
		gfs_er(gfs_er_),
		uer(uer_),
		uer_dot(uer_dot_),
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
		  double pw = x(lfsu.template child<PrimaryVariables::Pw>(),0);
		  double sw = x(lfsu.template child<PrimaryVariables::Sw>(),0);
		  //-----------------------------------------
		  // primary variables Pw_old and Sw_old
		  //-----------------------------------------
		  LFS2P lfs2p(gfs_2p);
		  LFS2PCache lfs2pcache(lfs2p);
		  VectorView2P u2p_old_view((*u2p_old));
		  lfs2p.bind(cell);
		  lfs2pcache.update();
		  u2p_old_view.bind(lfs2pcache);
		  std::vector<double> ul2p_old(lfs2p.size());
		  u2p_old_view.read( ul2p_old );
		  double pw_old = ul2p_old[ lfs2p.template child<PrimaryVariables::Pw>().localIndex(0)];
		  double sw_old = ul2p_old[ lfs2p.template child<PrimaryVariables::Sw>().localIndex(0)];
		  //-----------------------------------------
		  // variables porosity, cf, and dt_por
		  //-----------------------------------------
		  LFSEr lfser(gfs_er);
		  LFSErCache lfsercache(lfser);
		  VectorViewEr uer_view((*uer));
		  VectorViewEr uer_dot_view((*uer_dot));
		  lfser.bind(cell);
		  lfsercache.update();
		  uer_view.bind(lfsercache);
		  uer_dot_view.bind(lfsercache);
		  std::vector<double> uler(lfser.size());
		  uer_view.read( uler );
		  std::vector<double> uler_dot(lfser.size());
		  uer_dot_view.read( uler_dot );
		  double porosity = uler[ lfser.template child<PrimaryVariables::porosity>().localIndex(0)];
		  double cf = uler[ lfser.template child<PrimaryVariables::Cf>().localIndex(0)];
		  double dt_porosity = uler_dot[ lfser.template child<PrimaryVariables::porosity>().localIndex(0)];

		  // evaluate properties at cell center
		  auto pc = param.CapillaryPressure(cell_center_global,sw,porosity);
		  auto Cw = param.WettingPhaseCompressibility();
		  auto Cn = param.NonwettingPhaseCompressibility();
		  auto qw = param.WettingVolumetricSource(cell_center_global);
		  auto qn = param.NonwettingVolumetricSource(cell_center_global);

		  // RESIDUALS
		  //-----------------------------------------
		  // wetting-phase mass balance
		  //-----------------------------------------
		  double term_w = sw * dt_porosity
				  	    + porosity * (sw-sw_old)/(*dt)
						+ porosity * sw * Cw * (pw-pw_old)/(*dt)
						- qw;
		  r.accumulate(lfsu.template child<PrimaryVariables::Pw>(), 0., term_w*cell_volume );

		  //-----------------------------------------
		  // nonwetting-phase mass balance
		  //-----------------------------------------
		  double term_n = (1.-sw) * dt_porosity
				  	  	- porosity * (sw-sw_old)/(*dt)
						+ porosity * (1.-sw) * Cn * (pw-pw_old)/(*dt)
						+ porosity * (1.-sw) * Cn * param.dPc_dSw(cell_center_global,sw,porosity) * (sw-sw_old)/(*dt)
						+ porosity * (1.-sw) * Cn * param.dPc_dPor(cell_center_global,sw,porosity) * dt_porosity
						- qn;
		  r.accumulate(lfsu.template child<PrimaryVariables::Sw>(), 0., term_n*cell_volume );

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
	        // pw_s, pw_n, sw_s, sw_n
	        //-----------------------------------------
	        double pw_s = x_s(lfsu_s.template child<PrimaryVariables::Pw>(),0);
	        double sw_s = x_s(lfsu_s.template child<PrimaryVariables::Sw>(),0);
	        double pw_n = x_n(lfsu_n.template child<PrimaryVariables::Pw>(),0);
	        double sw_n = x_n(lfsu_n.template child<PrimaryVariables::Sw>(),0);
	        //-----------------------------------------
	        // variable porosity
	        //-----------------------------------------
	        LFSEr lfser_s(gfs_er);
	        LFSEr lfser_n(gfs_er);
	        LFSErCache lfsercache_s(lfser_s);
	        LFSErCache lfsercache_n(lfser_n);
	        VectorViewEr uer_view_s((*uer));
	        VectorViewEr uer_view_n((*uer));
			lfser_s.bind(cell_inside);
			lfser_n.bind(cell_outside);
			lfsercache_s.update();
			lfsercache_n.update();
			uer_view_s.bind(lfsercache_s);
			uer_view_n.bind(lfsercache_n);
			std::vector<double> uler_s(lfser_s.size());
			uer_view_s.read( uler_s );
			std::vector<double> uler_n(lfser_n.size());
			uer_view_n.read( uler_n );
			double por_s = uler_s[ lfser_s.template child<PrimaryVariables::porosity>().localIndex(0)];
			double por_n = uler_n[ lfser_n.template child<PrimaryVariables::porosity>().localIndex(0)];
			double cf_s = uler_s[ lfser_s.template child<PrimaryVariables::Cf>().localIndex(0)];
			double cf_n = uler_n[ lfser_n.template child<PrimaryVariables::Cf>().localIndex(0)];

	        // evaluate properties
	        // at self center
	        auto pc_s = param.CapillaryPressure(inside_cell_center_global,sw_s,por_s);
	        auto rhow_s = param.WettingPhaseDensity(pw_s,cf_s);
	        auto rhon_s = param.NonwettingPhaseDensity(pw_s+pc_s);
	        auto muw_s = param.WettingPhaseViscosity();
	        auto mun_s = param.NonwettingPhaseViscosity();
	        auto K_s = std::abs(param.Permeability(inside_cell_center_global,por_s) * normal);
	        auto krw_s = param.WettingRelativePermeability(inside_cell_center_global,sw_s);
	        auto krn_s = param.NonwettingRelativePermeability(inside_cell_center_global,sw_s);
	        // at neighbour center
	        auto pc_n = param.CapillaryPressure(outside_cell_center_global,sw_n,por_n);
	        auto rhow_n = param.WettingPhaseDensity(pw_n,cf_n);
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

	        //-----------------------------------------
	        // variables porosity, cf
	        //-----------------------------------------
	        LFSEr lfser_s(gfs_er);
	        LFSErCache lfsercache_s(lfser_s);
	        VectorViewEr uer_view_s((*uer));
			lfser_s.bind(cell_inside);
			lfsercache_s.update();
			uer_view_s.bind(lfsercache_s);
			std::vector<double> uler_s(lfser_s.size());
			uer_view_s.read( uler_s );
			double por_s = uler_s[ lfser_s.template child<PrimaryVariables::porosity>().localIndex(0)];
			double por_n = por_s;
			double cf_s = uler_s[ lfser_s.template child<PrimaryVariables::Cf>().localIndex(0)];
			double cf_n = cf_s;

				        // boundary conditions
	        //-----------------------------------------
	        // Boundary types
	        //-----------------------------------------
	        auto t_new = (*time)+(*dt);
	        auto bct_w = w_bc.type(ig,face_center_local,t_new);
			auto bct_n = n_bc.type(ig,face_center_local,t_new);
			//-----------------------------------------
			// dirichlet boundary values
			//-----------------------------------------
			auto pw_D = w_bc.dirichlet_value(ig,face_center_local,t_new,cf_n);
			auto sw_D = n_bc.dirichlet_value(ig,face_center_local,t_new);
			//-----------------------------------------
			// neumann boundary values
			//-----------------------------------------
			auto jw_N = w_bc.neumann_value(ig,face_center_local,t_new);
			auto jn_N = n_bc.neumann_value(ig,face_center_local,t_new);
			//-----------------------------------------
			// outflow boundary values
			//-----------------------------------------
			auto jw_O = w_bc.outflow_value(ig,face_center_local,t_new);
			auto jn_O = n_bc.outflow_value(ig,face_center_local,t_new);

			// compute primary vars at local self and neighbour centers
	        //-----------------------------------------
	        // pw_s, pw_n, sw_s, sw_n
	        //-----------------------------------------
	        double pw_s = x(lfsu.template child<PrimaryVariables::Pw>(),0);
	        double sw_s = x(lfsu.template child<PrimaryVariables::Sw>(),0);
	        double pw_n = pw_s;
	        if( bct_w==ConvectionDiffusionBoundaryConditions::Dirichlet ) pw_n = pw_D;
	        double sw_n = sw_s;
	        if( bct_n==ConvectionDiffusionBoundaryConditions::Dirichlet ) sw_n = sw_D;

	        // evaluate properties
	        // at self center
	        auto pc_s = param.CapillaryPressure(inside_cell_center_global,sw_s,por_s);
	        auto rhow_s = param.WettingPhaseDensity(pw_s,cf_s);
	        auto rhon_s = param.NonwettingPhaseDensity(pw_s+pc_s);
	        auto muw_s = param.WettingPhaseViscosity();
	        auto mun_s = param.NonwettingPhaseViscosity();
	        auto K_s = std::abs(param.Permeability(inside_cell_center_global,por_s) * normal);
	        auto krw_s = param.WettingRelativePermeability(inside_cell_center_global,sw_s);
	        auto krn_s = param.NonwettingRelativePermeability(inside_cell_center_global,sw_s);
	        // at neighbour center
	        auto pc_n = param.CapillaryPressure(face_center_global,sw_n,por_n);
	        auto rhow_n = param.WettingPhaseDensity(pw_n,cf_n);
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
			// wetting-phase mass balance
			//-----------------------------------------
			if( bct_w==ConvectionDiffusionBoundaryConditions::Neumann ){
				double term_w = jw_N * normal;
				r.accumulate(lfsu.template child<PrimaryVariables::Pw>() , 0,  term_w*face_volume);
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

#endif /* PROBLEM09_OPERATORS_LOP_2PFLOW_HH_ */
