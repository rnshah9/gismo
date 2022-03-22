/** @file -

@brief -

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): -
*/

#include <gismo.h>

#ifdef GISMO_WITH_IPOPT
#include <gsIpopt/gsOptProblem.h>
#endif

using namespace gismo;

#ifdef GISMO_WITH_IPOPT

template <typename T>
class gsParameterization2DExample
//class gsParameterization2DExample : public gsOptProblem<T>  // When I use this one, I will get some link errors
{
//public:
	//typedef typename std::function<gsSparseMatrix<real_t>(gsVector<real_t> const &)>    Jacobian_t;
	//typedef typename std::function<gsVector<real_t>(gsVector<real_t> const &) >         Residual_t;

public:

	gsParameterization2DExample(const gsDofMapper & mapper, const gsMultiPatch<T> & bRep,
		const index_t & method, const bool & plot_init)
		:
		m_mapper(mapper),
		m_bRep(bRep),
		//m_mp(mp),
		m_eps(1e-8), // in our paper, we set this parameter 
					// always equals to 0.05 * Area
		m_method(method),
		m_plot(false),
		m_plot_init(plot_init)
	{
		// TODO: assign m_mp by using initialization
		// m_mp = ;
		//initialization();

		// parameters related with design variables
		// Number of design variables 
		m_numDesignVars = m_mapper.freeSize();

		m_desLowerBounds.resize(m_numDesignVars);
		m_desUpperBounds.resize(m_numDesignVars);

		// parameters related with constraints
		// There is no constraint in our method
		m_numConstraints = 0;

		m_conLowerBounds.resize(m_numConstraints);
		m_conUpperBounds.resize(m_numConstraints);
		//m_conLowerBounds[0] = m_conUpperBounds[0] = 0;

		// parameters related with nonzero entries in the Constraint Jacobian 
		// this->currentDesign_into(mp, m_curDesign);
		m_numConJacNonZero = m_numDesignVars;
		m_conJacRows.resize(m_numConJacNonZero);
		m_conJacCols.resize(m_numConJacNonZero);
	}


public:

	void enablePlot() { m_plot = true; }
	void enableStress() { m_stress = true; }

	void disablePlot() { m_plot = false; }
	void disableStress() { m_stress = false; }

	/*void currentDesign_into(const gsMultiPatch<T> & mp, gsMatrix<T> & currentDesign)
	{
		currentDesign.resize(m_numDesignVars, 1);
		for (size_t i = 0; i != mp.nPatches(); i++)
			for (size_t c = 0; c != mp.targetDim(); c++)
				for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)
					if (m_mapper.is_free(k, i, c))
						currentDesign(m_mapper.index(k, i, c), 0) = mp.patch(i).coefs()(k, c);
	}

	void setCurrentDesign(const gsMultiPatch<T> & mp)
	{
		this->currentDesign_into(mp, m_curDesign);
	}*/

	//void initialization()
	//{
	//	// TO DO: different initialization methods: 
	//	// 1. discrete Coons 2. Smoothness energy 3. Spring model etc.
	//	switch (m_method)
	//	{
	//	case 1:
	//	{
	//		// Spring model method
	//		gsInfo << "Using spring patch construction.\n";
	//		gsSpringPatch<T> spring(m_bRep);
	//		gsInfo << "Created a " << spring.compute() << "\n";
	//		//if (save) gsWrite(spring.result(), "result_patch");
	//		m_mp.addPatch(spring.result());
	//		if (m_plot_init)
	//			gsWriteParaview(m_mp, "mp_init_spring", 1000, true, true);
	//		break;
	//	}
	//	case 2:
	//	{
	//		// Cross Approximation patch method
	//		// it does not work???
	//		gsInfo << "Using cross approximation construction.\n";
	//		gsCrossApPatch<T> cross(m_bRep);
	//		gsDebug << "xxxxxxxxxxxxxxxxxxx" << "\n";
	//		gsInfo << "Created a " << cross.compute() << "\n";
	//		//if (save) gsWrite(spring.result(), "result_patch");
	//		m_mp.addPatch(cross.result());
	//		if (m_plot_init)
	//			gsWriteParaview(m_mp, "mp_init_cross", 1000, true, true);
	//		break;
	//	}
	//	case 3:
	//	{
	//		// 
	//		break;
	//	}
	//	case 4:
	//	{
	//		// Smoothness energy method
	//		break;
	//	}
	//	case 0:
	//	default:
	//		// discrete Coons method
	//		gsInfo << "Using Coons' patch construction.\n";
	//		gsCoonsPatch<T> coons(m_bRep);
	//		gsInfo << "Created a " << coons.compute() << "\n";
	//		//if (save) gsWrite(coons.result(), "result_patch");
	//		m_mp.addPatch(coons.result());
	//		if (m_plot_init)
	//			gsWriteParaview(m_mp, "mp_init_coons", 1000, true, true);
	//		break;
	//	}

	//	/*if (save)
	//		gsInfo << "Result saved to result_patch.xml\n";
	//	else
	//		gsInfo << "Done. No output created, re-run with --save to get xml "
	//		"file containing the data.\n";*/

	//}

	void initialization()
	{
		// TO DO: different initialization methods: 
		// 1. discrete Coons 2. Smoothness energy 3. Spring model etc.
		switch (m_method)
		{
		case 1:
		{
			// Spring model method
			gsInfo << "Using spring patch construction.\n";
			gsSpringPatch<T> spring(m_bRep);
			gsInfo << "Created a " << spring.compute() << "\n";
			//if (save) gsWrite(spring.result(), "result_patch");
			m_mp.addPatch(spring.result());
			if (m_plot_init)
				gsWriteParaview(m_mp, "mp_init_spring", 1000, true, true);
			break;
		}
		case 2:
		{
			// Cross Approximation patch method
			// Question: only works for 3D surfaces? i.e., for mapping: x: \mathbb{R}^2 --> \mathbb{R}^3
			// I am not sure, ask Hugo or Matthias later.
			gsInfo << "Using cross approximation construction.\n";
			gsCrossApPatch<T> cross(m_bRep);
			gsDebug << "xxxxxxxxxxxxxxxxxxx" << "\n";
			gsInfo << "Created a " << cross.compute() << "\n";
			//if (save) gsWrite(spring.result(), "result_patch");
			m_mp.addPatch(cross.result());
			if (m_plot_init)
				gsWriteParaview(m_mp, "mp_init_cross", 1000, true, true);
			break;
		}
		case 3:
		{
			// consturt a parameterization with the inner control points all equal to (0, 0)
			// TODO: make the following step easier to handle
			// TODO: make this method dimensional-independent

			// must type conversion? any other alternative solution?
			// it seems not good here.. uhhh...
			// i.e., how can I get the knot vector of patch(i)?
			gsBSplineBasis<real_t> westBoundary = static_cast<gsBSplineBasis<real_t>&> (m_bRep.basis(0));
			gsBSplineBasis<real_t> northBoundary = static_cast<gsBSplineBasis<real_t>&> (m_bRep.basis(3));

			// 1. construction of a knot vector for each direction
			gsKnotVector<> uKnot = northBoundary.knots();
			gsKnotVector<> vKnot = westBoundary.knots();

			// 2. construction of a tensor basis
			gsTensorBSplineBasis<2, real_t> basis(uKnot, vKnot);

			// 3. construction of a coefficients
			index_t noCtrPts = northBoundary.size() * westBoundary.size();
			gsMatrix<> coefs(noCtrPts, m_bRep.targetDim());

			for (index_t col = 0; col != westBoundary.size(); col++)
			{
				for (index_t row = 0; row != northBoundary.size(); row++)
				{
					size_t idx = row + col * northBoundary.size();
					if (col == 0)
					{
						// South boundary case
						coefs(idx, 0) = m_bRep.patch(3).coefs()(row, 0);
						coefs(idx, 1) = m_bRep.patch(3).coefs()(row, 1);
						continue;
					}
					else if (col == westBoundary.size() - 1)
					{
						// North boundary case
						coefs(idx, 0) = m_bRep.patch(2).coefs()(row, 0);
						coefs(idx, 1) = m_bRep.patch(2).coefs()(row, 1);
						continue;
					}
					else if (row == 0)
					{
						// West boundary case
						coefs(idx, 0) = m_bRep.patch(0).coefs()(col, 0);
						coefs(idx, 1) = m_bRep.patch(0).coefs()(col, 1);
						continue;
					}
					else if (row == northBoundary.size() - 1)
					{
						// East boundary case
						coefs(idx, 0) = m_bRep.patch(1).coefs()(col, 0);
						coefs(idx, 1) = m_bRep.patch(1).coefs()(col, 1);
						continue;
					}
					else
					{
						// inner control points
						coefs(idx, 0) = 0;
						coefs(idx, 1) = 0;
					}
				}
			}

			// 4. putting basis and coefficients toghether
			gsTensorBSpline<2, real_t> samePoint(basis, coefs);
			m_mp.addPatch(samePoint);
			if (m_plot_init)
				gsWriteParaview(m_mp, "mp_init_samePoint", 1000, true, true);
			break;
		}
		case 4:
		{
			// Smoothness energy method
			// TODO: need to implement later...
			// However, the results seems the same as Spring model method

			break;
		}
		case 0:
		default:
			// discrete Coons method
			gsInfo << "Using Coons' patch construction.\n";
			gsCoonsPatch<T> coons(m_bRep);
			gsInfo << "Created a " << coons.compute() << "\n";
			//if (save) gsWrite(coons.result(), "result_patch");
			m_mp.addPatch(coons.result());
			if (m_plot_init)
				gsWriteParaview(m_mp, "mp_init_coons", 1000, true, true);

			//gsDebug << convert_mp_to_vec() << "\n";
			break;
		}

		/*if (save)
		gsInfo << "Result saved to result_patch.xml\n";
		else
		gsInfo << "Done. No output created, re-run with --save to get xml "
		"file containing the data.\n";*/

	}

	//void makeMapper()
	//{
	//	// Now, we set all the inner control points as optimization variables
	//	// It is possible to set only a part of them as optimization variables later

	//	////! [Make mapper for the design DoFs]
	//	m_mapper.init(gsMultiBasis<>(m_mp), m_mp.targetDim());
	//	//gsDofMapper m_mapper(gsMultiBasis<>(m_mp), m_mp.targetDim());
	//	// 1. Mark the vertical displacements of the control points (except the corners) as design variables
	//	//      a. Eliminate boundary control points
	//	
	//	gsMatrix<index_t> idx;
	//	for (size_t p = 0; p != m_mp.nPatches(); p++)
	//	{
	//		idx = m_mp.basis(p).allBoundary(); // if it need to compute all basis or not? 
	//										   // if YES, is there any more efficient way?

	//		for (size_t c = 0; c != idx.size(); c++)
	//		{
	//			for (size_t d = 0; d != m_mp.targetDim(); d++)
	//			{
	//				m_mapper.eliminateDof(idx(c), p, d);
	//			}
	//		}
	//	}
	//	m_mapper.finalize();

	//	/*gsDebug << "#Numb of free  variables is " << m_mapper.freeSize() << "\n";
	//	gsDebug << "#Numb of fixed variables is " << m_mapper.boundarySize() << "\n";
	//	gsDebug << "#Numb of total variables is " << m_mapper.size() << "\n";

	//	for (size_t p = 0; p != m_mp.nPatches(); p++)
	//	for (size_t k = 0; k != m_mp.basis(p).size(); k++)
	//	for (size_t d = 0; d != m_mp.targetDim(); d++)
	//	gsDebug << "p=" << p << "; k=" << k << "; d=" << d <<
	//	(m_mapper.is_free(k, p, d) ? " is free" : " is eliminated") << "\n";*/
	//}

	T computeArea()
	{
		// compute the area of the computational domain by B-Rep
		// using the following Green's formulation:
		// S = \int_{\Omega} 1 d \Omega
		//   = \oint_{\partial \Omega} x(t) y'(t) dt
		//   = \sum_{1}^4 \int_0^1 x_i(t) y'_i(t) dt
		// Here, one must take care of the orientation of the boundary curves

		T result = 0.;

		gsBSplineBasis<real_t> westBoundary = static_cast< gsBSplineBasis<real_t>& > (m_bRep.basis(0));
		gsBSplineBasis<real_t> northBoundary = static_cast< gsBSplineBasis<real_t>& > (m_bRep.basis(3));

		// line integral along the South and the North boundary

		// 1. get Gauss points and weights
		//northBoundary.eval_into();
		gsOptionList legendreOpts;
		legendreOpts.addInt("quRule", "Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule", gsQuadrature::GaussLegendre);
		legendreOpts.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0);
		legendreOpts.addInt("quB", "Number of quadrature points: quA*deg + quB", 1);
		legendreOpts.addSwitch("overInt", "Apply over-integration or not?", false);
		gsQuadRule<real_t>::uPtr legendre = gsQuadrature::getPtr(northBoundary, legendreOpts);

		gsMatrix<> points;
		gsVector<> weights;

		gsBasis<real_t>::domainIter domIt = northBoundary.makeDomainIterator();

		gsMatrix<> gaussPts(northBoundary.dim(), 0);
		gsMatrix<> gaussWts(1, 0); // should be a gsVector? but how to operate gsVector?

								   // TODO: here need to be optimized! especially when different
								   //       intergration scheme is used!!!
		legendre->mapTo(domIt->lowerCorner(), domIt->upperCorner(),
			points, weights);
		gsMatrix<index_t> gaussIdx(points.cols(), 0);

		index_t start;
		for (; domIt->good(); domIt->next())
		{
			//---------------------------------------------------------------------------
			// Gauss-Legendre rule (w/o over-integration)
			legendre->mapTo(domIt->lowerCorner(), domIt->upperCorner(),
				points, weights);

			start = gaussPts.cols();
			gaussPts.conservativeResize(Eigen::NoChange, gaussPts.cols() + points.cols());
			gaussPts.block(0, start, gaussPts.rows(), points.cols()) = points;

			gaussWts.conservativeResize(Eigen::NoChange, gaussWts.cols() + points.cols());
			gaussWts.block(0, start, gaussWts.rows(), weights.rows()) = weights.transpose();

			gsVector<index_t> localIdx = northBoundary.active(domIt->lowerCorner());

			gaussIdx.conservativeResize(Eigen::NoChange, gaussIdx.cols() + points.cols());
			gaussIdx.block(0, start, gaussIdx.rows(), points.cols()) = localIdx.replicate(1, points.cols());
		}

		// 2. Perform Gauss integration
		gsMatrix<> basisValAlongBoundaryNS;
		gsMatrix<> basis1stDervsAlongBoundaryNS;

		northBoundary.eval_into(gaussPts, basisValAlongBoundaryNS);

		//gsMatrix<real_t> allBasis1stDervs;
		northBoundary.deriv_into(gaussPts, basis1stDervsAlongBoundaryNS);

		/*gsDebug << "basis funtion: " << basisValAlongBoundaryNS << "\n";
		gsDebug << "Size of basis: " << basisValAlongBoundaryNS.rows() << " x " << basis1stDervsAlongBoundaryNS.cols() << "\n";
		gsDebug << "basis dervs: " << basis1stDervsAlongBoundaryNS << "\n";
		gsDebug << "Size of dervs: " << basis1stDervsAlongBoundaryNS.rows() << " x " << basis1stDervsAlongBoundaryNS.cols() << "\n";
		gsDebug << "index: " << gaussIdx << "\n";
		gsDebug << "Size of index: " << gaussIdx.rows() << " x " << gaussIdx.cols() << "\n";*/

		// in MATLAB, we can use localBasis = basisValAlongBoundaryNS(:,gp)
		// is there similar method in GiSmo?
		for (index_t gp = 0; gp != basisValAlongBoundaryNS.cols(); gp++)
		{
			gsVector<T> localBasis(basisValAlongBoundaryNS.rows());
			gsVector<T> localBasis1stDervs(basisValAlongBoundaryNS.rows());
			gsVector<T> localCoordXS(basisValAlongBoundaryNS.rows());
			gsVector<T> localCoordYS(basisValAlongBoundaryNS.rows());
			gsVector<T> localCoordXN(basisValAlongBoundaryNS.rows());
			gsVector<T> localCoordYN(basisValAlongBoundaryNS.rows());

			for (index_t i = 0; i != basisValAlongBoundaryNS.rows(); i++)
			{
				localBasis(i) = basisValAlongBoundaryNS(i, gp);
				localBasis1stDervs(i) = basis1stDervsAlongBoundaryNS(i, gp);

				localCoordXS(i) = m_bRep.patch(3).coefs()(gaussIdx(i, gp), 0);
				localCoordYS(i) = m_bRep.patch(3).coefs()(gaussIdx(i, gp), 1);

				localCoordXN(i) = m_bRep.patch(2).coefs()(gaussIdx(i, gp), 0);
				localCoordYN(i) = m_bRep.patch(2).coefs()(gaussIdx(i, gp), 1);
			}

			/*real_t xS = localBasis * localCoordXS;
			real_t ySDervs = localBasis1stDervs * localCoordYS;
			real_t xN = localBasis * localCoordXN;
			real_t yNDervs = localBasis1stDervs * localCoordYN;*/

			//result += (xS * ySDervs - xN * yNDervs) * gaussWts(gp);

			//gsDebug << typeid(localBasis.dot(localCoordXS)).name() << "\n";
			result += (localBasis.dot(localCoordXS) * localBasis1stDervs.dot(localCoordYS) -
				localBasis.dot(localCoordXN) * localBasis1stDervs.dot(localCoordYN)) * gaussWts(gp);

			/*gsDebug << localBasis * localCoordXS << "\n";
			gsDebug << localBasis1stDervs * localCoordYS << "\n";
			gsDebug << localBasis * localCoordXN << "\n";
			gsDebug << localBasis1stDervs * localCoordYN << "\n";*/
		}

		// line integral along the West and the East boundary

		// 1. get Gauss points and weights
		//westBoundary.eval_into();

		// 2. Perform Gauss integration

		legendre = gsQuadrature::getPtr(westBoundary, legendreOpts);

		/*gsMatrix<> points;
		gsVector<> weights;*/

		//gsBasis<real_t>::domainIter domIt = westBoundary.makeDomainIterator();
		domIt = westBoundary.makeDomainIterator();

		gaussPts.resize(westBoundary.dim(), 0);
		gaussWts.resize(1, 0);
		gaussIdx.resize(points.cols(), 0);

		legendre->mapTo(domIt->lowerCorner(), domIt->upperCorner(),
			points, weights);
		gaussIdx.resize(points.cols(), 0);

		//index_t start;
		for (; domIt->good(); domIt->next())
		{
			//---------------------------------------------------------------------------
			// Gauss-Legendre rule (w/o over-integration)
			legendre->mapTo(domIt->lowerCorner(), domIt->upperCorner(),
				points, weights);

			start = gaussPts.cols();
			gaussPts.conservativeResize(Eigen::NoChange, gaussPts.cols() + points.cols());
			gaussPts.block(0, start, gaussPts.rows(), points.cols()) = points;

			gaussWts.conservativeResize(Eigen::NoChange, gaussWts.cols() + points.cols());
			gaussWts.block(0, start, gaussWts.rows(), weights.rows()) = weights.transpose();

			gsVector<index_t> localIdx = westBoundary.active(domIt->lowerCorner());

			gaussIdx.conservativeResize(Eigen::NoChange, gaussIdx.cols() + points.cols());
			gaussIdx.block(0, start, gaussIdx.rows(), points.cols()) = localIdx.replicate(1, points.cols());
		}

		// 2. Perform Gauss integration
		gsMatrix<> basisValAlongBoundaryWE;
		gsMatrix<> basis1stDervsAlongBoundaryWE;

		westBoundary.eval_into(gaussPts, basisValAlongBoundaryWE);

		//gsMatrix<real_t> allBasis1stDervs;
		westBoundary.deriv_into(gaussPts, basis1stDervsAlongBoundaryWE);

		/*gsDebug << "basis funtion: " << basisValAlongBoundaryNS << "\n";
		gsDebug << "Size of basis: " << basisValAlongBoundaryNS.rows() << " x " << basis1stDervsAlongBoundaryNS.cols() << "\n";
		gsDebug << "basis dervs: " << basis1stDervsAlongBoundaryNS << "\n";
		gsDebug << "Size of dervs: " << basis1stDervsAlongBoundaryNS.rows() << " x " << basis1stDervsAlongBoundaryNS.cols() << "\n";
		gsDebug << "index: " << gaussIdx << "\n";
		gsDebug << "Size of index: " << gaussIdx.rows() << " x " << gaussIdx.cols() << "\n";*/

		// in MATLAB, we can use localBasis = basisValAlongBoundaryNS(:,gp)
		// is there similar method in GiSmo?
		for (index_t gp = 0; gp != basisValAlongBoundaryWE.cols(); gp++)
		{
			gsVector<T> localBasis(basisValAlongBoundaryWE.rows());
			gsVector<T> localBasis1stDervs(basisValAlongBoundaryWE.rows());
			gsVector<T> localCoordXW(basisValAlongBoundaryWE.rows());
			gsVector<T> localCoordYW(basisValAlongBoundaryWE.rows());
			gsVector<T> localCoordXE(basisValAlongBoundaryWE.rows());
			gsVector<T> localCoordYE(basisValAlongBoundaryWE.rows());

			for (index_t i = 0; i != basisValAlongBoundaryWE.rows(); i++)
			{
				localBasis(i) = basisValAlongBoundaryWE(i, gp);
				localBasis1stDervs(i) = basis1stDervsAlongBoundaryWE(i, gp);

				localCoordXW(i) = m_bRep.patch(0).coefs()(gaussIdx(i, gp), 0);
				localCoordYW(i) = m_bRep.patch(0).coefs()(gaussIdx(i, gp), 1);

				localCoordXE(i) = m_bRep.patch(1).coefs()(gaussIdx(i, gp), 0);
				localCoordYE(i) = m_bRep.patch(1).coefs()(gaussIdx(i, gp), 1);
			}

			/*real_t xS = localBasis * localCoordXS;
			real_t ySDervs = localBasis1stDervs * localCoordYS;
			real_t xN = localBasis * localCoordXN;
			real_t yNDervs = localBasis1stDervs * localCoordYN;*/

			//result += (xS * ySDervs - xN * yNDervs) * gaussWts(gp);

			//gsDebug << typeid(localBasis.dot(localCoordXS)).name() << "\n";
			result += (localBasis.dot(localCoordXE) * localBasis1stDervs.dot(localCoordYE) -
				localBasis.dot(localCoordXW) * localBasis1stDervs.dot(localCoordYW)) * gaussWts(gp);

			/*gsDebug << localBasis * localCoordXS << "\n";
			gsDebug << localBasis1stDervs * localCoordYS << "\n";
			gsDebug << localBasis * localCoordXN << "\n";
			gsDebug << localBasis1stDervs * localCoordYN << "\n";*/
		}

		return result;
	}

	void preComputeBasis()
	{
		// pre-compute basis funtion values at Gauss points to further accelerate 

		// TODO: To reduce computational memory, only univariate information is stored in memory, 
		// and the multivariate information is computed on the fly whenever needed

		// TODO: with different quadrature rules? 
		// Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule

		// TODO: now, only one patch case is considered; multi-patch case will be implemented later

		/// Step 1: Constructs a quadrature rule
		// Gauss Legendre
		gsOptionList legendreOpts;
		legendreOpts.addInt("quRule", "Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule", gsQuadrature::GaussLegendre);
		legendreOpts.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0);
		legendreOpts.addInt("quB", "Number of quadrature points: quA*deg + quB", 1);
		legendreOpts.addSwitch("overInt", "Apply over-integration or not?", false);
		gsQuadRule<real_t>::uPtr legendre = gsQuadrature::getPtr(m_mp.basis(0), legendreOpts);

		gsMatrix<> points;
		gsVector<> weights;

		gsBasis<real_t>::domainIter domIt = m_mp.basis(0).makeDomainIterator();

		gsMatrix<> gaussPts(m_mp.basis(0).dim(), 0);
		gsMatrix<> gaussWts(1, 0); // should be a gsVector? but how to operate gsVector?
		
		// TODO: here need to be optimized! especially when different
		//       intergration scheme is used!!!
		legendre->mapTo(domIt->lowerCorner(), domIt->upperCorner(),
			points, weights);
		gsMatrix<index_t> gaussIdx(points.cols(), 0);

		index_t start;
		for (; domIt->good(); domIt->next())
		{
			//---------------------------------------------------------------------------
			// Gauss-Legendre rule (w/o over-integration)
			legendre->mapTo(domIt->lowerCorner(), domIt->upperCorner(),
				points, weights);
			/*gsInfo << "* \t Gauss-Legendre\n"
				<< "- points:\n" << points << "\n"
				<< "- weights:\n" << weights << "\n";*/
			start = gaussPts.cols();
			gaussPts.conservativeResize(Eigen::NoChange, gaussPts.cols() + points.cols());
			gaussPts.block(0, start, gaussPts.rows(), points.cols()) = points;

			gaussWts.conservativeResize(Eigen::NoChange, gaussWts.cols() + points.cols());
			gaussWts.block(0, start, gaussWts.rows(), weights.rows()) = weights.transpose();

			/*gsDebug << "lowerCorner = " << domIt->lowerCorner() << "\n";
			gsDebug << "upperCorner = " << domIt->upperCorner() << "\n";*/

			/*gsVector<T> coord = domIt->lowerCorner();
			m_mp.basis(0).active(coord);
			gsDebug << "active dof = " << m_mp.basis(0).active(domIt->lowerCorner()) << "\n";*/
			
			gsVector<index_t> localIdx = m_mp.basis(0).active(domIt->lowerCorner());

			gaussIdx.conservativeResize(Eigen::NoChange, gaussIdx.cols() + points.cols());
			gaussIdx.block(0, start, gaussIdx.rows(), points.cols()) = localIdx.replicate(1, points.cols());

			/*gsDebug << "active dof = " << localIdx << "\n";
			gsDebug << gaussIdx << "\n";
			gsDebug << "Size of gaussIdx = " << gaussIdx.rows() << "x" << gaussIdx.cols() << "\n";
			*/

			/*m_mp.patch(0).active_into(coord, idx);
			gsDebug << "active dof = " << m_mp.patch(0).active(coord) << "\n";*/

			//gsDebug << "active dof =" << m_mp.active(domIt->lowerCorner()) << "\n";
			//gsInfo << "check weights: " << weights.sum() << "\n";
			//gsInfo << weights.sum() << "\n";
			//---------------------------------------------------------------------------
		}

		/*gsInfo << gaussPts << "\n";
		gsDebug << "Size of gaussPts = " << gaussPts.rows() << "x" << gaussPts.cols() << "\n";*/

		/*gsInfo << gaussWts << "\n";
		gsDebug << "Size of gaussWts = " << gaussWts.rows() << "x" << gaussWts.cols() << "\n";

		gsDebug << gaussIdx << "\n";
		gsDebug << "Size of gaussIdx = " << gaussIdx.rows() << "x" << gaussIdx.cols() << "\n";*/

		m_gaussWts = gaussWts;
		m_gaussIdx = gaussIdx;
		/// Step 2: evalute basis functions with its 1st and 2nd derivates at all integration points

		//gsMatrix<real_t> allBasisVal;
		m_mp.basis(0).eval_into(gaussPts, m_allBasisVal);

		/*gsDebug << "xxxxxxxxxxxxxxxxxxxxxxx" << "\n";
		gsDebug << m_allBasisVal << "\n";
		gsDebug << "Size of allBasisVal = " << m_allBasisVal.rows() << "x" << m_allBasisVal.cols() << "\n";*/

		//gsMatrix<real_t> allBasis1stDervs;
		m_mp.basis(0).deriv_into(gaussPts, m_allBasis1stDervs);

		/*gsDebug << "xxxxxxxxxxxxxxxxxxxxxxx" << "\n";
		gsDebug << m_allBasis1stDervs << "\n";
		gsDebug << "Size of allBasis1stDervs = " << m_allBasis1stDervs.rows() << "x" << m_allBasis1stDervs.cols() << "\n";*/

		//gsMatrix<real_t> allBasis2ndDervs;
		m_mp.basis(0).deriv2_into(gaussPts, m_allBasis2ndDervs);

		/*gsDebug << "xxxxxxxxxxxxxxxxxxxxxxx" << "\n";
		gsDebug << m_allBasis2ndDervs << "\n";
		gsDebug << "Size of allBasis2ndDervs = " << m_allBasis2ndDervs.rows() << "x" << m_allBasis2ndDervs.cols() << "\n";*/

	}

	gsVector<T> convert_mp_to_vec()
	{
		// TODO: Now, it just set all free variables into the vector u
		// I will integrate it with the initialization step

		gsVector<T> u(m_mapper.freeSize());

		index_t idx;
		for (size_t i = 0; i != m_mp.nPatches(); i++)
		{
			for (size_t c = 0; c != m_mp.targetDim(); c++)
			{
				for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)
					// if it is possible to just loop over the free index
					// since this function is called very often during optimization
					if (m_mapper.is_free(k, i, c))
					{
						idx = m_mapper.index(k, i, c);
						u(idx) = m_mp.patch(i).coefs()(k, c);
					}
			}
		}

		return u;
	}
	
	void convert_vec_to_mp(const gsVector<T> & u)
	{
		// Make the geometry
		index_t idx;
		for (size_t i = 0; i != m_mp.nPatches(); i++)
		{
			for (size_t c = 0; c != m_mp.targetDim(); c++)
			{
				for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)
					// if it is possible to just loop over the free index
					// since this function is called very often during optimization
					if (m_mapper.is_free(k, i, c))
					{
						idx = m_mapper.index(k, i, c);
						m_mp.patch(i).coefs()(k, c) = u(idx);
					}
			}
		}
	}

	// What are the differences between gsAsConstVector<T> and gsVector<T>??
	//T evalObj(const gsAsConstVector<T> & u) const
	T evalObj(const gsVector<T> & u) const
	{
		//m_mp = m_mp0;

		// but why this line does not work???
		//convert_vec_to_mp(u);

		// Make the geometry
		index_t idx;
		for (size_t i = 0; i != m_mp.nPatches(); i++)
		{
			for (size_t c = 0; c != m_mp.targetDim(); c++)
			{
				for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)
					// if it is possible to just loop over the free index
					// since this function is called very often during optimization
					if (m_mapper.is_free(k, i, c))
					{
						idx = m_mapper.index(k, i, c);
						m_mp.patch(i).coefs()(k, c) = u(idx);
					}
			}
		}

		T result = 0;
		// For ease of understanding, we use $\hat{G}$ and localCoord here,
		// I will optimize this part later
		gsMatrix<T> hatG(m_mp.targetDim(), m_allBasisVal.rows());
		gsMatrix<T> localCoord(m_allBasisVal.rows(), m_mp.targetDim());

		gsMatrix<T> localIdx();
		for (index_t gp = 0; gp != m_allBasisVal.cols(); gp++)
		{
			for (index_t i = 0; i != m_allBasisVal.rows(); i++)
			{
				hatG(0, i) = m_allBasis1stDervs(2*i, gp);
				hatG(1, i) = m_allBasis1stDervs(2*i+1, gp);

				localCoord(i, 0) = m_mp.patch(0).coefs()(m_gaussIdx(i, gp), 0);
				localCoord(i, 1) = m_mp.patch(0).coefs()(m_gaussIdx(i, gp), 1);
			}

			gsMatrix<T> jacMat = hatG * localCoord;
			real_t detJac = jacMat.determinant();

			if (detJac < m_eps)
				result += (m_eps - detJac)*m_gaussWts(gp);
		}

		// 
		// gsDebug << m_mp.patch(0).coefs() << "\n";
		return result;
	}

	//// TO DO: implement analytical gradient
	void gradObj_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
	{
		GISMO_NO_IMPLEMENTATION;
	}

	void evalCon_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
	{
		GISMO_NO_IMPLEMENTATION;
	}

	void jacobCon_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
	{
		 GISMO_NO_IMPLEMENTATION;
	}

private:

	/// Number of design variables
	int m_numDesignVars;

	/// Number of constraints
	int m_numConstraints;

	/// Number of nonzero entries in the Constraint Jacobian
	int m_numConJacNonZero;

	/// Lower bounds for the design variables
	gsVector<T> m_desLowerBounds;

	/// Upper bounds for the design variables
	gsVector<T> m_desUpperBounds;

	/// Lower bounds for the constraints
	gsVector<T> m_conLowerBounds;

	/// Upper bounds for the constraints
	gsVector<T> m_conUpperBounds;

	/// Constraint Jacobian non-zero entries rows
	std::vector<index_t> m_conJacRows;

	/// Constraint Jacobian non-zero entries columns
	std::vector<index_t> m_conJacCols;

	/// Current design variables (and starting point )
	gsMatrix<T> m_curDesign;

	// IpIpt???
	//using gsOptProblem<T>::m_numDesignVars;				// Number of design variables
	//using gsOptProblem<T>::m_desLowerBounds;			// Lower bounds for the design variables
	//using gsOptProblem<T>::m_desUpperBounds;			// Upper bounds for the design variables

	//using gsOptProblem<T>::m_numConstraints;			// Number of constraints
	//using gsOptProblem<T>::m_conLowerBounds;			// Lower bounds for the constraints
	//using gsOptProblem<T>::m_conUpperBounds;			// Upper bounds for the constraints

	//using gsOptProblem<T>::m_numConJacNonZero;			// Number of nonzero entries in the Constraint Jacobian
	//using gsOptProblem<T>::m_conJacRows;				// Constraint Jacobian non-zero entries rows
	//using gsOptProblem<T>::m_conJacCols;				// Constraint Jacobian non-zero entries columns

	//using gsOptProblem<T>::m_curDesign;					// Current design variables (and starting point )


	//gsThinShellAssemblerBase<T> * m_assembler; 
	//const gsMatrix<T> m_refPoint;
	//std::vector<gsMatrix<bool>> m_freeIndices;
	gsDofMapper m_mapper;

	const gsMultiPatch<T> m_bRep;
	mutable gsMultiPatch<T> m_mp;

	mutable gsMultiBasis<T> m_computeBases;

	mutable gsMatrix<T> m_allBasisVal;
	mutable gsMatrix<T> m_allBasis1stDervs;
	mutable gsMatrix<T> m_allBasis2ndDervs;
	mutable gsMatrix<T> m_gaussWts;
	mutable gsMatrix<index_t> m_gaussIdx;

	T m_eps;
	/*mutable gsMultiPatch<T> m_computePatches;
	mutable gsMultiBasis<T> m_computeBases;

	Residual_t m_Residual;
	Jacobian_t m_Jacobian;

	index_t m_numRefine, m_numElevate;*/

	index_t m_method;
	bool m_plot, m_plot_init;
};
#endif

int main(int argc, char* argv[])
{

	//////////////////// STEP 1: read the file of the boundary represention /////////////////

	bool save = false;
	index_t method = 0;
	real_t tol = 1e-4;
	bool plot_init = false;

	// Load XML file containing the boundary represention 
	// TODO: give some different cases as defalut benchmarks, 
	//       otherwise deal with the input filename
	std::string filename_input("duck_boundary.xml");
	std::string filename_output("results");

	// Read input from command line arguments
	//! [Parse command line]
	gsCmdLine cmd("Hi, give me a file (eg: .xml) containing boundary representation (B-Rep)"
		"and I will try to parameterize it!");

	cmd.addPlainString("input", "Name of the input file containing boundary data", filename_input);
	cmd.addString("o", "output", "Name of the output file", filename_output);
	cmd.addInt("m", "method", "Method: 0 Coons' patch (default), 1 Spring patch, 2: Cross-Ap. patch", method);
	cmd.addReal("t", "tolerance", "Tolerance for identifing patch interfaces", tol);
	cmd.addSwitch("save", "Save result in XML format", save);
	cmd.addSwitch("plotInit", "Plot resulting initaialization for Paraview", plot_init);
	try { cmd.getValues(argc, argv); }
	catch (int rv) { return rv; }
	//! [Parse command line]

	// Load XML file

	gsMultiPatch<> bRep;
	gsReadFile<>(filename_input, bRep);
	GISMO_ENSURE(!bRep.empty(), "The gsMultiPatch is empty - maybe file is missing or corrupt.");
	gsInfo << "Got " << bRep << "\n";
	bRep.computeTopology(tol);
	GISMO_ENSURE(bRep.isClosed(), "The boundary is not closed, adjust tolerance.");
	bRep.closeGaps(tol);

	gsWriteParaview(bRep, "bRep", 1000, true, true);

	////////////////////////////////////////////////////////////////////////////////////////////////////

	//-------------------------------------------------------------------------------------------------
	//////////////////////// STEP 1.2: constructing a fake initialization///////////////////////////////
	// Question: If it is possible to avoid this step? I will ask Hugo or Matthias to make this part more elegant!!
	// consturt a parameterization with the inner control points all equal to (0, 0)
	// to make the following step easier to handle

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// must type conversion? any other alternative solution?
	// it seems not good here.. uhhh...
	// i.e., how can I get the knot vector of patch(i)?
	gsBSplineBasis<real_t> westBoundary = static_cast< gsBSplineBasis<real_t>& > (bRep.basis(0));
	gsBSplineBasis<real_t> northBoundary = static_cast< gsBSplineBasis<real_t>& > (bRep.basis(3));

	/*gsDebug << westBoundary.knots() << "\n";
	gsDebug << northBoundary.knots() << "\n";*/

	// 1.2.1. construction of a knot vector for each direction
	gsKnotVector<> uKnot = northBoundary.knots();
	gsKnotVector<> vKnot = westBoundary.knots();

	/*index_t noPtsX = northBoundary.size();  // # of control points in $\xi$-direction
	index_t noPtsY = westBoundary.size(); //# of control points in $\eta$-direction
	gsDebug << "Size of control points = " << noPtsX
		<< " x " << noPtsY << "\n";*/

	// 1.2.2. construction of a basis
	gsTensorBSplineBasis<2, real_t> basis(uKnot, vKnot);

	// 1.2.3. construction of a coefficients
	index_t noCtrPts = northBoundary.size() * westBoundary.size();
	gsMatrix<> coefs(noCtrPts, bRep.targetDim());

	// 1.2.4. putting basis and coefficients toghether
	gsTensorBSpline<2, real_t>  surface(basis, coefs);

	gsMultiPatch<> mp;
	mp.addPatch(surface);
	
	//////////////////////// STEP 2: constructing analysis-suitable parameterization ///////////////////
	//// eliminate foldovers first

	/// for test purpose
	//! [Refine and elevate]
	//index_t numRefine = 1;
	//index_t numElevate = 1;

	//// p-refine
	//if (numElevate != 0)
	//	mp.degreeElevate(numElevate);

	//// h-refine
	//for (int r = 0; r < numRefine; ++r)
	//	mp.uniformRefine();

	//gsWriteParaview<>(mp.basis(0), "basis");
	//! [Refine and elevate]
	
	////! [Make mapper for the design DoFs]
	gsDofMapper mapper(gsMultiBasis<>(mp), mp.targetDim());
	// 1. Mark the vertical displacements of the control points (except the corners) as design variables
	//      a. Eliminate boundary control points
	gsMatrix<index_t> idx;
	for (size_t p = 0; p != mp.nPatches(); p++)
	{

		idx = mp.basis(p).allBoundary(); // if it need to compute all basis or not? 
										 // if YES, is there any more efficient way?

		for (size_t c = 0; c != idx.size(); c++)
		{
			for (size_t d = 0; d != mp.targetDim(); d++)
			{
				mapper.eliminateDof(idx(c), p, d);
			}
		}
	}
	mapper.finalize();

	gsDebug << "#Numb of free  variables is "  << mapper.freeSize() << "\n";
	gsDebug << "#Numb of fixed variables is " << mapper.boundarySize() << "\n";
	gsDebug << "#Numb of total variables is " << mapper.size() << "\n";

	/*for (size_t p = 0; p != mp.nPatches(); p++)
		for (size_t k = 0; k != mp.basis(p).size(); k++)
			for (size_t d = 0; d != mp.targetDim(); d++)
				gsDebug << "p=" << p << "; k=" << k << "; d=" << d <<
				(mapper.is_free(k, p, d) ? " is free" : " is eliminated") << "\n";*/

	gsParameterization2DExample<real_t> opt(mapper, bRep, method, plot_init);

	gsDebug << "Area of computational domain is " 
		<< opt.computeArea() << "\n";

	opt.initialization();
	//opt.makeMapper();
	opt.preComputeBasis();

	gsVector<real_t> initU = opt.convert_mp_to_vec();

	//gsInfo << "inital guess = " << initU << "\n";
	gsInfo << "value of objective function is: " << opt.evalObj(initU) << "\n";

	//gsDebug << opt.makeInitGuess() << "\n";
	//gsInfo << "control points are:" << mp.patch(0).coefs() << "\n";





	// Make initialization guess

	//gsTensorBSplineBasis<2, real_t> m_computeBases = gsTensorBSplineBasis<real_t>(mp);

	//gsMultiBasis<real_t> m_computeBases = gsMultiBasis<real_t>(mp);
	//gsDebug << "#Numb of control points is: " << m_computeBases.size() << "\n";

	//gsDebug << "The dimensions is: " << m_computeBases.dim() << "\n";

	///// step 2.1: Gauss integration rule
	//// Gauss Legendre
	//gsOptionList legendreOpts;
	//legendreOpts.addInt("quRule", "Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule", gsQuadrature::GaussLegendre);
	//legendreOpts.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0);
	//legendreOpts.addInt("quB", "Number of quadrature points: quA*deg + quB", 1);
	//legendreOpts.addSwitch("overInt", "Apply over-integration or not?", false);
	//gsQuadRule<real_t>::uPtr legendre = gsQuadrature::getPtr(mp.basis(0), legendreOpts);

	//gsMatrix<> points;
	//gsVector<> weights;

	//gsBasis<real_t>::domainIter domIt = mp.basis(0).makeDomainIterator();

	//gsMatrix<> gaussPts(mp.basis(0).dim(), 0);
	//gsMatrix<> gaussWts(1, 0);

	//index_t start;
	//for (; domIt->good(); domIt->next())
	//{
	//	//---------------------------------------------------------------------------
	//	// Gauss-Legendre rule (w/o over-integration)
	//	legendre->mapTo(domIt->lowerCorner(), domIt->upperCorner(),
	//		points, weights);
	//	gsInfo << "* \t Gauss-Legendre\n"
	//		<< "- points:\n" << points << "\n"
	//		<< "- weights:\n" << weights << "\n";
	//	start = gaussPts.cols();
	//	gaussPts.conservativeResize(Eigen::NoChange, gaussPts.cols() + points.cols());
	//	gaussPts.block(0, start, gaussPts.rows(), points.cols()) = points;

	//	gaussWts.conservativeResize(Eigen::NoChange, gaussWts.cols() + points.cols());
	//	gaussWts.block(0, start, gaussWts.rows(), weights.rows()) = weights.transpose();

	//	gsInfo << "check weights: " << weights.sum() << "\n";
	//	//gsInfo << weights.sum() << "\n";
	//	//---------------------------------------------------------------------------
	//}

	//gsInfo << gaussPts << "\n";
	//gsDebug << "Size of gaussPts = " << gaussPts.rows() << "x" << gaussPts.cols() << "\n";

	//gsInfo << gaussWts << "\n";
	//gsDebug << "Size of gaussWts = " << gaussWts.rows() << "x" << gaussWts.cols() << "\n";


	///// step 2.2 evalute basis functions and their derivs

	//gsMatrix<real_t> allBasisVal;
	//mp.basis(0).eval_into(gaussPts, allBasisVal);

	//gsDebug << "xxxxxxxxxxxxxxxxxxxxxxx" << "\n";
	//gsDebug << allBasisVal << "\n";
	//gsDebug << "Size of allBasisVal = " << allBasisVal.rows() << "x" << allBasisVal.cols() << "\n";

	//gsMatrix<real_t> allBasis1stDervs;
	//mp.basis(0).deriv_into(gaussPts, allBasis1stDervs);

	//gsDebug << "xxxxxxxxxxxxxxxxxxxxxxx" << "\n";
	//gsDebug << allBasis1stDervs << "\n";
	//gsDebug << "Size of allBasis1stDervs = " << allBasis1stDervs.rows() << "x" << allBasis1stDervs.cols() << "\n";

	//gsMatrix<real_t> allBasis2ndDervs;
	//mp.basis(0).deriv2_into(gaussPts, allBasis2ndDervs);

	//gsDebug << "xxxxxxxxxxxxxxxxxxxxxxx" << "\n";
	//gsDebug << allBasis2ndDervs << "\n";
	//gsDebug << "Size of allBasis2ndDervs = " << allBasis2ndDervs.rows() << "x" << allBasis2ndDervs.cols() << "\n";





	/*for (size_t i=0; i!=oooo.col; i++)
		gsInfo << oooo(i) << "\n";*/
	
	////m_computeBases.eval_into(ssss, oooo);

	//gsDebug << oooo << "\n";

	//gsMatrix<real_t> ssss(3,1);
	//gsFuncData<real_t> oooo;

	//ssss << 0.1, 0.2, 0.3;

	///*for (size_t i = 0; i != 3; i++)
	//	ssss[i,1] = 0.3;*/

	//gsDebug << ssss << "\n";

	//m_computeBases.compute(ssss, oooo);

	//gsDebug << oooo << "\n";

	//m_computeBases.compute(in, out);

	//gsDebug << out << "\n";

	/*in(0.0);

	gsDebug << in[0] << "\n";*/

	//m_computeBases.compute(in, out);

	//gsDebug << out(0) << "\n";

	//opt.initialization(); // different initialization methods 1. discrete Coons 2. Smoothness energy 3. Spring model etc.

	//gsDebugVar(gsAsVector<real_t>(u));
	//gsDebugVar(opt.evalObj(u));

	//opt.solve();

	////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////// STEP 3: improving the quality of parameterization ///////////////////////////
	//// Another optimization problem


	////////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////// STEP 4: output the resulting parameterization ////////////////////////
	// writing the resulting parameterization to a G+Smo .xml file
	// filename_output is a string. The extention .xml is added automatically

	// TODO: other formats? To make it easy to visulization

	gsFileData<> fd;
	fd << mp;
	fd.save(filename_output);
	gsInfo << "Wrote G+Smo file:     " << filename_output << ".xml \n";

	gsWriteParaview(mp, "mp", 1000, true, true);

	//gsWrite(mp, filename_output);
	////////////////////////////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////// STEP 5: VISUALIZATION ///////////////////////////////////////
	// TODO: visualization, in MATLAB first? then ?

	////////////////////////////////////////////////////////////////////////////////////////////////////

	system("pause");
	return EXIT_SUCCESS;
}

