/*
 * l2projection.hh
 *
 *  Created on: Apr 6, 2021
 *      Author: sgupta
 */

#ifndef PROBLEM09_OPERATORS_LOP_L2PROJECTION_HH_
#define PROBLEM09_OPERATORS_LOP_L2PROJECTION_HH_

template <class GV, class PARAMS, class GFS2P, class U2P, class GFSEr, class UEr>
class LocalOperatorL2Projection :
	public Dune::PDELab::NumericalJacobianApplyVolume		<LocalOperatorL2Projection<GV,PARAMS,GFS2P,U2P,GFSEr,UEr> >,
	public Dune::PDELab::NumericalJacobianVolume			<LocalOperatorL2Projection<GV,PARAMS,GFS2P,U2P,GFSEr,UEr> >,
	public Dune::PDELab::NumericalJacobianApplySkeleton		<LocalOperatorL2Projection<GV,PARAMS,GFS2P,U2P,GFSEr,UEr> >,
	public Dune::PDELab::NumericalJacobianSkeleton			<LocalOperatorL2Projection<GV,PARAMS,GFS2P,U2P,GFSEr,UEr> >,
	public Dune::PDELab::NumericalJacobianApplyBoundary		<LocalOperatorL2Projection<GV,PARAMS,GFS2P,U2P,GFSEr,UEr> >,
	public Dune::PDELab::NumericalJacobianBoundary			<LocalOperatorL2Projection<GV,PARAMS,GFS2P,U2P,GFSEr,UEr> >,
	public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
	public Dune::PDELab::FullVolumePattern,
	public Dune::PDELab::LocalOperatorDefaultFlags,
	public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV& gv;
	const PARAMS& param;
	GFS2P gfs_2p;
	U2P *u2p;
	GFSEr gfs_er;
	UEr *uer;
	unsigned int intorder ;
	double epsilon;

public:
	// pattern assembly flags
	enum { doPatternVolume 	 = true  };
	enum { doPatternSkeleton = false };

	// residual assembly flags
	enum { doAlphaVolume  	= true 	};
	enum { doLambdaVolume  	= true 	};
	enum { doAlphaSkeleton  = false };
	enum { doAlphaBoundary  = false };
	enum { doLambdaBoundary = false };

	typedef Dune::PDELab::LocalFunctionSpace<GFS2P> LFS2P;
	typedef Dune::PDELab::LFSIndexCache<LFS2P> LFS2PCache;
	typedef typename U2P::template LocalView<LFS2PCache> VectorView2P;

	typedef Dune::PDELab::LocalFunctionSpace<GFSEr> LFSEr;
	typedef Dune::PDELab::LFSIndexCache<LFSEr> LFSErCache;
	typedef typename UEr::template LocalView<LFSErCache> VectorViewEr;

	// constructor stores parameters
	LocalOperatorL2Projection(	const GV& gv_,
								const PARAMS& param_,
								GFS2P gfs_2p_,
								U2P	*u2p_,
								GFSEr gfs_er_,
								UEr	*uer_,
								unsigned int intorder_ = 6,
								double 		 epsilon_ = 1.e-6)
	: gv(gv_),
	  param(param_),
	  gfs_2p(gfs_2p_),
	  u2p(u2p_),
	  gfs_er(gfs_er_),
	  uer(uer_),
	  intorder( intorder_ ),
	  epsilon( epsilon_ )
	{}

	// volume integral depending on test and ansatz functions
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	{
		// define types
		using RF = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
		using RangeType = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
		using JacobianType = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;

		// dimensions
		const int dim = EG::Entity::dimension;

		// Get cell
		const auto& cell = eg.entity();

		// Get geometry
		auto geo = eg.geometry();

		// evaluate diffusion tensor at cell center, assume it is constant over elements
		auto ref_el = referenceElement(geo);
		auto localcenter = ref_el.position(0,0);

		// Transformation matrix
		typename EG::Geometry::JacobianInverseTransposed jac;

		// loop over quadrature points
		for (const auto& ip : quadratureRule(geo,intorder))
		{
			auto ip_global = geo.global(ip.position());

			// evaluate basis functions
			std::vector<RangeType> phi_pw(lfsu.template child<PrimaryVariables::Pw_l2>().size());
			lfsu.template child<PrimaryVariables::Pw_l2>().finiteElement().localBasis().evaluateFunction(ip.position(),phi_pw);

			std::vector<RangeType> psi_pw(lfsv.template child<PrimaryVariables::Pw_l2>().size());
			lfsv.template child<PrimaryVariables::Pw_l2>().finiteElement().localBasis().evaluateFunction(ip.position(),psi_pw);

			std::vector<RangeType> phi_pn(lfsu.template child<PrimaryVariables::Pn_l2>().size());
			lfsu.template child<PrimaryVariables::Pn_l2>().finiteElement().localBasis().evaluateFunction(ip.position(),phi_pn);

			std::vector<RangeType> psi_pn(lfsv.template child<PrimaryVariables::Pn_l2>().size());
			lfsv.template child<PrimaryVariables::Pn_l2>().finiteElement().localBasis().evaluateFunction(ip.position(),psi_pn);

			// evaluate u (u->phase pressures)
            RF u_pw=0.0;
            for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pw_l2>().size(); i++)
              u_pw += x(lfsu.template child<PrimaryVariables::Pw_l2>(),i)*phi_pw[i];

            RF u_pn=0.0;
            for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pn_l2>().size(); i++)
              u_pn += x(lfsu.template child<PrimaryVariables::Pn_l2>(),i)*phi_pn[i];

			// evaluate gradient of basis functions
			std::vector<JacobianType> jsu_pw(lfsu.template child<PrimaryVariables::Pw_l2>().size());
			lfsu.template child<PrimaryVariables::Pw_l2>().finiteElement().localBasis().evaluateJacobian(ip.position(),jsu_pw);

			std::vector<JacobianType> jsv_pw(lfsv.template child<PrimaryVariables::Pw_l2>().size());
			lfsv.template child<PrimaryVariables::Pw_l2>().finiteElement().localBasis().evaluateJacobian(ip.position(),jsv_pw);

			std::vector<JacobianType> jsu_pn(lfsu.template child<PrimaryVariables::Pn_l2>().size());
			lfsu.template child<PrimaryVariables::Pn_l2>().finiteElement().localBasis().evaluateJacobian(ip.position(),jsu_pn);

			std::vector<JacobianType> jsv_pn(lfsv.template child<PrimaryVariables::Pn_l2>().size());
			lfsv.template child<PrimaryVariables::Pn_l2>().finiteElement().localBasis().evaluateJacobian(ip.position(),jsv_pn);

			// transform gradients of shape functions to real element
			jac = geo.jacobianInverseTransposed(ip.position());

			// evaluade gradients of shape fncs.
			std::vector<Dune::FieldVector<RF,dim> > gradphi_pw(lfsu.template child<PrimaryVariables::Pw_l2>().size());
            for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pw_l2>().size(); i++){
                gradphi_pw[i] = 0.0;
                jac.umv(jsu_pw[i][0],gradphi_pw[i]);
            }

			std::vector<Dune::FieldVector<RF,dim> > gradpsi_pw(lfsv.template child<PrimaryVariables::Pw_l2>().size());
            for (unsigned int i=0; i<lfsv.template child<PrimaryVariables::Pw_l2>().size(); i++){
                gradpsi_pw[i] = 0.0;
                jac.umv(jsv_pw[i][0],gradpsi_pw[i]);
            }

			std::vector<Dune::FieldVector<RF,dim> > gradphi_pn(lfsu.template child<PrimaryVariables::Pn_l2>().size());
            for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pn_l2>().size(); i++){
                gradphi_pn[i] = 0.0;
                jac.umv(jsu_pn[i][0],gradphi_pn[i]);
            }

			std::vector<Dune::FieldVector<RF,dim> > gradpsi_pn(lfsv.template child<PrimaryVariables::Pn_l2>().size());
            for (unsigned int i=0; i<lfsv.template child<PrimaryVariables::Pn_l2>().size(); i++){
                gradpsi_pn[i] = 0.0;
                jac.umv(jsv_pn[i][0],gradpsi_pn[i]);
            }

			// compute gradient of u_pw, u_pn
            Dune::FieldVector<RF,dim> gradu_pw(0.0);
            for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pw_l2>().size(); i++)
              gradu_pw.axpy(x(lfsu.template child<PrimaryVariables::Pw_l2>(),i),gradphi_pw[i]);

            Dune::FieldVector<RF,dim> gradu_pn(0.0);
            for (unsigned int i=0; i<lfsu.template child<PrimaryVariables::Pn_l2>().size(); i++)
              gradu_pn.axpy(x(lfsu.template child<PrimaryVariables::Pn_l2>(),i),gradphi_pn[i]);

			// integrate
			// FV --> FEM
			//|| u_FV - u_FE || --> Min
			RF factor = ip.weight() * geo.integrationElement(ip.position());
			double tmp = 0.;
			for (int j=0; j<lfsv.template child<PrimaryVariables::Pw_l2>().size(); j++){
				tmp = epsilon*( gradu_pw*gradpsi_pw[j] ) + u_pw*psi_pw[j] ;
				r.accumulate(lfsv.template child<PrimaryVariables::Pw_l2>(),j,tmp*factor );
			}

			tmp = 0.;
			for (int j=0; j<lfsv.template child<PrimaryVariables::Pn_l2>().size(); j++){
				tmp = epsilon*( gradu_pn*gradpsi_pn[j] ) + u_pn*psi_pn[j] ;
				r.accumulate(lfsv.template child<PrimaryVariables::Pn_l2>(),j,tmp*factor );
			}

		}//End Quadrature Rule

	}


	// volume integral depending only on test functions
	template<typename EG, typename LFSV, typename R>
	void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
	{
		// define types
		using RF = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
		using RangeType = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
		using JacobianType = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;

		// dimensions
		const int dim = EG::Entity::dimension;

		// Get cell
		const auto& cell = eg.entity();

		// Get geometry
		auto geo = eg.geometry();

		// evaluate diffusion tensor at cell center, assume it is constant over elements
		auto ref_el = referenceElement(geo);
		auto localcenter = ref_el.position(0,0);
		auto globalcenter = cell.geometry().global(localcenter);

		LFS2P lfs_2p( gfs_2p );
		LFS2PCache lfscache_2p( lfs_2p );
		VectorView2P u2p_view( (*u2p) );
		lfs_2p.bind( (cell) );
		lfscache_2p.update();
		u2p_view.bind( lfscache_2p );
		std::vector<RF> ul2p(lfs_2p.size());
		u2p_view.read( ul2p );

		LFSEr lfs_er( gfs_er );
		LFSErCache lfscache_er( lfs_er );
		VectorViewEr uer_view( (*uer) );
		lfs_er.bind( (cell) );
		lfscache_er.update();
		uer_view.bind( lfscache_er );
		std::vector<RF> uler(lfs_er.size());
		uer_view.read( uler );

		double Pw = ul2p[lfs_2p.template child<PrimaryVariables::Pw>().localIndex(0)];
		double Sw = ul2p[lfs_2p.template child<PrimaryVariables::Sw>().localIndex(0)];
		double por= uler[lfs_er.template child<PrimaryVariables::porosity>().localIndex(0)];
		double Pc = param.CapillaryPressure(globalcenter,Sw,por);
		double Pn=Pw+Pc;

		// loop over quadrature points
		for (const auto& ip : quadratureRule(geo,intorder))
		{
			// evaluate shape functions
			std::vector<RangeType> psi_pw(lfsv.template child<PrimaryVariables::Pw_l2>().size());
			lfsv.template child<PrimaryVariables::Pw_l2>().finiteElement().localBasis().evaluateFunction(ip.position(),psi_pw);

			std::vector<RangeType> psi_pn(lfsv.template child<PrimaryVariables::Pn_l2>().size());
			lfsv.template child<PrimaryVariables::Pn_l2>().finiteElement().localBasis().evaluateFunction(ip.position(),psi_pn);

			// integrate
			RF factor = ip.weight() * geo.integrationElement(ip.position());

			for (int i=0; i<lfsv.template child<PrimaryVariables::Pw_l2>().size(); i++){
				r.accumulate(lfsv.template child<PrimaryVariables::Pw_l2>(),i,(-Pw*psi_pw[i])*factor);
			}

			for (int i=0; i<lfsv.template child<PrimaryVariables::Pn_l2>().size(); i++){
				r.accumulate(lfsv.template child<PrimaryVariables::Pn_l2>(),i,(-Pn*psi_pn[i])*factor);
			}

		}//END: quadrature rule

	}//END: lambda_volume

};

#endif /* PROBLEM09_OPERATORS_LOP_L2PROJECTION_HH_ */
