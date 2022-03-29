/** @file -

@brief - A reference implementation of the following paper:
		Ji, Y., Yu, Y. Y., Wang, M. Y., & Zhu, C. G. (2021). 
		Constructing high-quality planar NURBS parameterization for 
		isogeometric analysis by adjustment control points and weights. 
		Journal of Computational and Applied Mathematics, 396, 113615.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): Ye Ji (jiye@mail.dlut.edu.cn)
*/

#include <gismo.h>

// for using HLBFGS
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif
#include "HLBFGS.h"
#include "Lite_Sparse_Matrix.h"

#include <chrono>

//#include <ultimaille/all.h>
#include <HLBFGS_wrapper.h>
using namespace UM;

// For gsOptProblem class, the default solver is IpOpt, but I was failed to install it 
// under Windows, so I use HLBFGS solver instead. This solver is lighter, just a few .h
// files, not with too many dependencies
#ifdef GISMO_WITH_IPOPT
#include <gsIpopt/gsOptProblem.h>
#endif

using namespace gismo;

#ifdef GISMO_WITH_IPOPT

template <typename T>
class gsParameterization2DExample
	//class gsParameterization2DExample : public gsOptProblem<T>  // When I use this one, I will get some link errors
{

public:

	gsParameterization2DExample(const gsMultiPatch<T> & bRep,
		const index_t & method, const bool & plot_init)
		:
		//m_mapper(mapper),
		m_bRep(bRep),
		//m_mp(mp),
		m_method(method),
		m_plot(false),
		m_plot_init(plot_init)
	{
		// TODO: assign m_mp by using initialization
		// m_mp = ;

		m_area = computeArea();
		m_eps = 0.05*m_area;
		initialization();
		makeMapper();

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

	void disablePlot() { m_plot = false; }

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
			//// Cross Approximation patch method
			//// Question: only works for 3D surfaces? i.e., for mapping: x: \mathbb{R}^2 --> \mathbb{R}^3
			//// I am not sure, ask Hugo or Matthias later.

			//gsInfo << "Using cross approximation construction.\n";
			//gsCrossApPatch<T> cross(m_bRep);
			//gsDebug << "xxxxxxxxxxxxxxxxxxx" << "\n";
			//gsInfo << "Created a " << cross.compute() << "\n";
			////if (save) gsWrite(spring.result(), "result_patch");
			//m_mp.addPatch(cross.result());
			//if (m_plot_init)
			//	gsWriteParaview(m_mp, "mp_init_cross", 1000, true, true);
			//break;
		}
		case 3:
		{
			// consturt a parameterization with the inner control points all equal to (0, 0)
			// TODO: make the following step easier to handle
			// TODO: make this method dimensional-independent

			// Practice: write a class named same point to put all the inner CtrPts to its barycenter,
			// refer to gsCoonsPatch class

			// first, get a coons patcj; and then put the inner CtrPts to (0,0)
			// We will implement it later like gsCoonsPatch class

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
			// However, the results seems the same as Spring model method?

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
			break;
		}

		/*if (save)
			gsInfo << "Result saved to result_patch.xml\n";
		else
			gsInfo << "Done. No output created, re-run with --save to get xml "
			"file containing the data.\n";*/
	}

	void makeMapper()
	{
		// Now, we set all the inner control points as optimization variables
		// It is possible to set only a part of them as optimization variables later

		////! [Make mapper for the design DoFs]
		m_mapper.init(gsMultiBasis<>(m_mp), m_mp.targetDim());
		//gsDofMapper m_mapper(gsMultiBasis<>(m_mp), m_mp.targetDim());
		// 1. Mark the vertical displacements of the control points (except the corners) 
		//    as design variables
		//      a. Eliminate boundary control points

		gsMatrix<index_t> idx;
		for (size_t p = 0; p != m_mp.nPatches(); p++)
		{
			idx = m_mp.basis(p).allBoundary(); // if it need to compute all basis or not? 
											   // if YES, is there any more efficient way?

			for (size_t c = 0; c != idx.size(); c++)
			{
				for (size_t d = 0; d != m_mp.targetDim(); d++)
				{
					m_mapper.eliminateDof(idx(c), p, d);
				}
			}
		}
		m_mapper.finalize();

		gsDebug << "#Numb of free  variables is " << m_mapper.freeSize() << "\n";
		gsDebug << "#Numb of fixed variables is " << m_mapper.boundarySize() << "\n";
		gsDebug << "#Numb of total variables is " << m_mapper.size() << "\n";

		/*for (size_t p = 0; p != m_mp.nPatches(); p++)
			for (size_t k = 0; k != m_mp.basis(p).size(); k++)
				for (size_t d = 0; d != m_mp.targetDim(); d++)
					gsDebug << "p=" << p << "; k=" << k << "; d=" << d <<
					(m_mapper.is_free(k, p, d) ? " is free" : " is eliminated") << "\n";*/
	}

	T computeArea()
	{
		// compute the area of the computational domain by B-Rep
		// using the following Green's formulation:
		// S = \int_{\Omega} 1 d \Omega
		//   = \oint_{\partial \Omega} x(t) y'(t) dt
		//   = \sum_{1}^4 \int_0^1 x_i(t) y'_i(t) dt

		// Here, one must take care of the orientation of the boundary curves
		// I found there exist some files get negative values of area
		// We should make this part more robust!! make it indenpent of boundary orientation
		// make some pre-check? or just get the absolute value of the result?

		// Or if the opposite (like NS and WE) boundary curves always with the same direction?

		T result = 0.;

		// Take care: here, we assume the opposite boundary are successive and original orientation
		//			  is positive!!
		// FIX IT LATER!! for the case that the order of boundary curves is randam
		// gsCoonPatch.hpp might be helpful!
		// line integral along the West and the East boundary
		index_t idxBdry = 0;
		gsBSplineBasis<T> westBoundary = static_cast<gsBSplineBasis<T>&> (m_bRep.basis(idxBdry));
		T resultWE = oppoBdryIntegral(westBoundary, idxBdry);

		// line integral along the North and the South boundary
		idxBdry = 2;
		gsBSplineBasis<T> northBoundary = static_cast<gsBSplineBasis<T>&> (m_bRep.basis(idxBdry));
		T resultNS = oppoBdryIntegral(northBoundary, idxBdry);

		return result = resultWE + resultNS;
	}

	T oppoBdryIntegral(gsBSplineBasis<T> bdryOnedirection, index_t idxBrdy) const
	{
		T area = 0.;
		// line integral along the opposite boundaries along one parametric direction

		// 1. get Gauss points and weights
		gsOptionList legendreOpts;
		legendreOpts.addInt("quRule", "Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule", gsQuadrature::GaussLegendre);
		legendreOpts.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0);
		legendreOpts.addInt("quB", "Number of quadrature points: quA*deg + quB", 1);
		legendreOpts.addSwitch("overInt", "Apply over-integration or not?", false);
		gsQuadRule<T>::uPtr legendre = gsQuadrature::getPtr(bdryOnedirection, legendreOpts);

		gsMatrix<> points;
		gsVector<> weights;

		gsBasis<T>::domainIter domIt = bdryOnedirection.makeDomainIterator();

		gsMatrix<> gaussPts(bdryOnedirection.dim(), 0);
		gsMatrix<> gaussWts(1, 0); // should be a gsVector? but how to operate gsVector?

		// TODO: here need to be optimized! especially when different
		//       intergration scheme is used!!!
		legendre->mapTo(domIt->lowerCorner(), domIt->upperCorner(), points, weights);
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

			gsVector<index_t> localIdx = bdryOnedirection.active(domIt->lowerCorner());

			gaussIdx.conservativeResize(Eigen::NoChange, gaussIdx.cols() + points.cols());
			gaussIdx.block(0, start, gaussIdx.rows(), points.cols()) = localIdx.replicate(1, points.cols());
			//---------------------------------------------------------------------------
		}

		// 2. Perform Gauss integration
		gsMatrix<> basisVal;
		gsMatrix<> basis1stDervs;

		bdryOnedirection.eval_into(gaussPts, basisVal);
		bdryOnedirection.deriv_into(gaussPts, basis1stDervs);

		// in MATLAB, we can use localBasis = basisValAlongBoundaryNS(:,gp)
		// is there similar method in GiSmo?
		for (index_t gp = 0; gp != basisVal.cols(); gp++)
		{
			gsVector<T> localBasis(basisVal.rows());
			gsVector<T> localBasis1stDervs(basisVal.rows());
			gsVector<T> localCoordXS(basisVal.rows());
			gsVector<T> localCoordYS(basisVal.rows());
			gsVector<T> localCoordXN(basisVal.rows());
			gsVector<T> localCoordYN(basisVal.rows());

			for (index_t i = 0; i != basisVal.rows(); i++)
			{
				localBasis(i) = basisVal(i, gp);
				localBasis1stDervs(i) = basis1stDervs(i, gp);

				localCoordXS(i) = m_bRep.patch(idxBrdy + 1).coefs()(gaussIdx(i, gp), 0);
				localCoordYS(i) = m_bRep.patch(idxBrdy + 1).coefs()(gaussIdx(i, gp), 1);

				localCoordXN(i) = m_bRep.patch(idxBrdy).coefs()(gaussIdx(i, gp), 0);
				localCoordYN(i) = m_bRep.patch(idxBrdy).coefs()(gaussIdx(i, gp), 1);
			}

			area += (localBasis.dot(localCoordXS) * localBasis1stDervs.dot(localCoordYS) -
				localBasis.dot(localCoordXN) * localBasis1stDervs.dot(localCoordYN)) * gaussWts(gp);

		}
		return area;
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

			start = gaussPts.cols();
			gaussPts.conservativeResize(Eigen::NoChange, gaussPts.cols() + points.cols());
			gaussPts.block(0, start, gaussPts.rows(), points.cols()) = points;

			gaussWts.conservativeResize(Eigen::NoChange, gaussWts.cols() + points.cols());
			gaussWts.block(0, start, gaussWts.rows(), weights.rows()) = weights.transpose();

			gsVector<index_t> localIdx = m_mp.basis(0).active(domIt->lowerCorner());

			gaussIdx.conservativeResize(Eigen::NoChange, gaussIdx.cols() + points.cols());
			gaussIdx.block(0, start, gaussIdx.rows(), points.cols()) = localIdx.replicate(1, points.cols());
			//---------------------------------------------------------------------------
		}

		m_gaussWts = gaussWts;
		m_gaussIdx = gaussIdx;

		/// Step 2: evalute basis functions with its 1st and 2nd derivates at all integration points
		m_mp.basis(0).eval_into(gaussPts, m_allBasisVal);

		/*gsDebug << "xxxxxxxxxxxxxxxxxxxxxxx" << "\n";
		gsDebug << m_allBasisVal << "\n";
		gsDebug << "Size of allBasisVal = " << m_allBasisVal.rows() << "x" << m_allBasisVal.cols() << "\n";*/

		m_mp.basis(0).deriv_into(gaussPts, m_allBasis1stDervs);

		/*gsDebug << "xxxxxxxxxxxxxxxxxxxxxxx" << "\n";
		gsDebug << m_allBasis1stDervs << "\n";
		gsDebug << "Size of allBasis1stDervs = " << m_allBasis1stDervs.rows() << "x" << m_allBasis1stDervs.cols() << "\n";*/

		m_mp.basis(0).deriv2_into(gaussPts, m_allBasis2ndDervs);

		/*gsDebug << "xxxxxxxxxxxxxxxxxxxxxxx" << "\n";
		gsDebug << m_allBasis2ndDervs << "\n";
		gsDebug << "Size of allBasis2ndDervs = " << m_allBasis2ndDervs.rows() << "x" << m_allBasis2ndDervs.cols() << "\n";*/

	}

	gsVector<T> convert_mp_to_gsvec() const
	{
		// TODO: Now, it just set all free variables into the vector u
		// I will integrate it with the initialization step

		gsVector<T> currentDesign(m_mapper.freeSize());

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
						currentDesign(idx) = m_mp.patch(i).coefs()(k, c);
					}
			}
		}

		return currentDesign;
	}

	void convert_gsVec_to_mp(const gsVector<T> & currentDesign, gsMultiPatch<T> & mp) const
	{
		// Make the geometry
		index_t idx;
		for (size_t i = 0; i != mp.nPatches(); i++)
		{
			for (size_t c = 0; c != mp.targetDim(); c++)
			{
				for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)

					// if it is possible to just loop over the free index
					// since this function is called very often during optimization
					if (m_mapper.is_free(k, i, c))
					{
						idx = m_mapper.index(k, i, c);
						mp.patch(i).coefs()(k, c) = currentDesign(idx);
					}
			}
		}
	}

	void convert_vec_to_mp(const std::vector<T>& currentDesign, gsMultiPatch<T> & mp) const
	{
		// Make the geometry
		index_t idx;
		for (size_t i = 0; i != mp.nPatches(); i++)
		{
			for (size_t c = 0; c != mp.targetDim(); c++)
			{
				for (index_t k = 0; k != m_mapper.patchSize(i, c); k++)

					// if it is possible to just loop over the free index
					// since this function is called very often during optimization
					if (m_mapper.is_free(k, i, c))
					{
						idx = m_mapper.index(k, i, c);
						mp.patch(i).coefs()(k, c) = currentDesign[idx];
					}
			}
		}
	}

	// What are the differences between gsAsConstVector<T> and gsVector<T>??
	//T evalObj(const gsAsConstVector<T> & u) const
	void objInterFreeFunc_into(const std::vector<double>& X, double& F, std::vector<double>& G) const
	{
		// Make the geometry
		convert_vec_to_mp(X, m_mp);

		gsVector<index_t> freeIdx = m_mapper.findFree(0);

		index_t nVarsCtrPts = freeIdx.size();
		index_t nAllCtrPts = m_mp.basis(0).size();
		std::vector<T> allG(2 * nAllCtrPts, 0);

		//gsVector<T> allG(2 * nAllCtrPts);
		// For ease of understanding, we use $\hat{G}$ and localCoord here,
		// I will optimize this part later
		gsMatrix<T> hatG(m_mp.targetDim(), m_allBasisVal.rows());
		gsMatrix<T> localCoord(m_allBasisVal.rows(), m_mp.targetDim());

		for (index_t gp = 0; gp != m_allBasisVal.cols(); gp++)
		{
			for (index_t i = 0; i != m_allBasisVal.rows(); i++)
			{
				hatG(0, i) = m_allBasis1stDervs(2 * i, gp);
				hatG(1, i) = m_allBasis1stDervs(2 * i + 1, gp);

				localCoord(i, 0) = m_mp.patch(0).coefs()(m_gaussIdx(i, gp), 0);
				localCoord(i, 1) = m_mp.patch(0).coefs()(m_gaussIdx(i, gp), 1);
			}

			gsMatrix<T> jacMat = hatG * localCoord;
			real_t detJac = jacMat.determinant();

			// TODO: need to optimize!!
			if (detJac < m_eps)
			{
				//gsDebug << "detJac < m_eps" << "\n";
				F += (m_eps - detJac) * m_gaussWts(gp);

				/*parDers_objFun(:, id) = parDers_objFun(:, id) - ...
					(det_Jac * (JacobiMat\G_cap)) * wtJ2(i, 1);*/

				for (index_t i = 0; i != hatG.cols(); i++)
				{
					
					allG[m_gaussIdx(i, gp)] += (jacMat(0, 1)*hatG(1, i) - jacMat(1, 1)*hatG(0, i)) * m_gaussWts(gp);
					allG[m_gaussIdx(i, gp) + nAllCtrPts] += (jacMat(1, 0)*hatG(0, i) - jacMat(0, 0)*hatG(1, i))*m_gaussWts(gp);
					/*allG[m_gaussIdx(i,gp)] += parDers_objFun(0,i);
					allG[m_gaussIdx(i,gp) + nAllCtrPts] += parDers_objFun(1,i);*/
				}
			}
		}

		for (index_t i = 0; i != nVarsCtrPts; i++)
		{
			G[i] = allG[freeIdx(i)];
			G[i + nVarsCtrPts] = allG[freeIdx(i) + nAllCtrPts];
		}
	}

	void objImprovePts_into(const std::vector<double>& X, double& F, std::vector<double>& G) const
	{
		// Make the geometry
		convert_vec_to_mp(X, m_mp);

		gsVector<index_t> freeIdx = m_mapper.findFree(0);

		index_t nVarsCtrPts = freeIdx.size();
		index_t nAllCtrPts = m_mp.basis(0).size();
		std::vector<T> allG(2 * nAllCtrPts, 0);

		//gsVector<T> allG(2 * nAllCtrPts);
		// For ease of understanding, we use $\hat{G}$ and localCoord here,
		// I will optimize this part later
		gsMatrix<T> hatG(m_mp.targetDim(), m_allBasisVal.rows());
		gsMatrix<T> localCoord(m_allBasisVal.rows(), m_mp.targetDim());

		for (index_t gp = 0; gp != m_allBasisVal.cols(); gp++)
		{
			for (index_t i = 0; i != m_allBasisVal.rows(); i++)
			{
				hatG(0, i) = m_allBasis1stDervs(2 * i, gp);
				hatG(1, i) = m_allBasis1stDervs(2 * i + 1, gp);

				localCoord(i, 0) = m_mp.patch(0).coefs()(m_gaussIdx(i, gp), 0);
				localCoord(i, 1) = m_mp.patch(0).coefs()(m_gaussIdx(i, gp), 1);
			}

			gsMatrix<T> jacMat = hatG * localCoord;
			real_t detJac = jacMat.determinant();

			// FIX IT!
			real_t jacMatFroNorm = jacMat.norm() * jacMat.inverse().norm();

			//real_t sigma1 = jacMat.singularValues()(0);
			//real_t sigma2 = jacMat.singularValues()(1); // make it more robust: 1 = svd.singularValues().size() - 1
			//real_t jacMatFroNorm = sigma1 / sigma2 + sigma2 / sigma1;
			//real_t detJac = sigma1 * sigma2;

			gsMatrix<T> parDers_objFun;
			if (detJac > m_eps)
			{
				F += (jacMatFroNorm + pow(detJac/m_area-1,2)) * m_gaussWts(gp);

				// TODO: need to optimize it!
				parDers_objFun = (2 / detJac * jacMat.transpose() * hatG 
					- jacMatFroNorm * jacMat.inverse() * hatG 
					+ 2 / m_area * (detJac / m_area - 1) * (detJac*jacMat.inverse()*hatG)) * m_gaussWts(gp);

				/*parDers_objFun(:, id) = parDers_objFun(:, id) + ...
					(2 * JacobiMat' * G_cap / det_Jac - energy_FNorm * (JacobiMat\G_cap)...
						+ 2 / area * (det_Jac / area - 1) * (det_Jac * (JacobiMat\G_cap))) * wtJ2(i, 1);*/
			}
			else
			{
				// TODO: make it simple, just make F to 1e10 and allG to 0? TEST IT!
				//F += (1e10 * pow(jacMat.norm(),2) + pow(detJac / m_area - 1, 2)) * m_gaussWts(gp);

				F = 1e15;
				return;
				/*parDers_objFun(:, id) = parDers_objFun(:, id) + ...
					(10 ^ 10 * 2 * sum(JacobiMat, 'all') * G_cap...
						+ 2 / area * (det_Jac / area - 1) * (det_Jac * (JacobiMat\G_cap))) * wtJ2(i, 1);*/

				// TODO: need to optimize it!
				/*parDers_objFun = (2e10 * pow(jacMat.norm(),2) *hatG
					+ 2 / m_area * (detJac / m_area - 1) * (detJac*jacMat.inverse()*hatG)) * m_gaussWts(gp);*/
			}
			for (index_t i = 0; i != hatG.cols(); i++)
			{
				allG[m_gaussIdx(i, gp)] += parDers_objFun(0, i);
				allG[m_gaussIdx(i, gp) + nAllCtrPts] += parDers_objFun(1, i);
			}
		}

		for (index_t i = 0; i != nVarsCtrPts; i++)
		{
			G[i] = allG[freeIdx(i)];
			G[i + nVarsCtrPts] = allG[freeIdx(i) + nAllCtrPts];
		}
	}

	T evalObj(const gsVector<T> & currentDesign) const
	{
		// Make the geometry
		convert_gsVec_to_mp(currentDesign, m_mp);

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
				hatG(0, i) = m_allBasis1stDervs(2 * i, gp);
				hatG(1, i) = m_allBasis1stDervs(2 * i + 1, gp);

				localCoord(i, 0) = m_mp.patch(0).coefs()(m_gaussIdx(i, gp), 0);
				localCoord(i, 1) = m_mp.patch(0).coefs()(m_gaussIdx(i, gp), 1);
			}

			gsMatrix<T> jacMat = hatG * localCoord;
			real_t detJac = jacMat.determinant();

			if (detJac < m_eps)
				result += (m_eps - detJac)*m_gaussWts(gp);
		}

		return result;
	}

	//// TO DO: implement analytical gradient
	//void gradObj_into(const gsAsConstVector<T> & currentDesign, gsAsVector<T> & result) const
	void gradObj_into(const gsVector<T> & currentDesign, gsVector<T> & result) const
	{
		//GISMO_NO_IMPLEMENTATION;

		// Now, it is computed numerically, and this part will be replaced by 
		// analytically gradient computation.
		const T h = T(1e-8);

		T currVal = evalObj(currentDesign);
		for (index_t i = 0; i != currentDesign.size(); i++)
		{
			gsVector<T> nextDesign = currentDesign;
			nextDesign(i) = nextDesign(i) + h;
			T nextVal = evalObj(nextDesign);
			result(i) = (nextVal - currVal) / h;
		}
	}

	void ElimFoldovers()
	{
		// Step 2. Eliminating foldovers
		// Section 4.2 in our paper

		int debug = 1;
		const LBFGS_Optimizer::func_grad_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
			F = 0;
			std::fill(G.begin(), G.end(), 0);
			objInterFreeFunc_into(X, F, G);
		};

		/*double E_prev, E;
		std::vector<double> trash(X.size());
		func(X, E_prev, trash);*/

		LBFGS_Optimizer opt(func);
		/*opt.gtol = bfgs_threshold;
		opt.maxiter = bfgs_maxiter;*/
		opt.gtol = 1e-6;
		opt.maxiter = 100;
		opt.run(X);

		//gsInfo << X << "\n";
		//func(X, E, trash);
		//if (debug > 0) gsInfo << "E: " << E << "\n";
	}

	void ImprovePts()
	{
		// Step 2. Eliminating foldovers
		// Section 4.2 in our paper

		int debug = 1;
		const LBFGS_Optimizer::func_grad_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
			F = 0;
			std::fill(G.begin(), G.end(), 0);
			objImprovePts_into(X, F, G);
		};

		LBFGS_Optimizer opt(func);
		/*opt.gtol = bfgs_threshold;
		opt.maxiter = bfgs_maxiter;*/
		opt.gtol = 1e-6;
		opt.maxiter = 1000;
		opt.run(X);
	}

	void compute()
	{
		// Eliminate foldovers
		ElimFoldovers();
		
		// Further impove the parameterization quality
		ImprovePts();

		// Push the result into the resulting parameterization (multi-patch)
		convert_vec_to_mp(X, m_mp);
		gsWriteParaview(m_mp, "result", 1000, false, false);
	}

	std::vector<T> convert_mp_to_vec() const
	{
		// TODO: Now, it just set all free variables into the vector u
		// I will integrate it with the initialization step

		std::vector<T> currentDesign(m_mapper.freeSize());

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
						currentDesign[idx] = m_mp.patch(i).coefs()(k, c);
					}
			}
		}

		return currentDesign;
	}

	void evalCon_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
	{
		GISMO_NO_IMPLEMENTATION;
	}

	void jacobCon_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
	{
		GISMO_NO_IMPLEMENTATION;
	}

	gsMultiPatch<> outputResult() const
	{
		return m_mp;
	}

public:
	std::vector<double> X;     // current geometry

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

	// IpOpt???
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

	T m_eps; // need to handle later, set m_eps = 0.05*S

	T m_area; // area of computational domain

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
	std::string filename_input("breps/duck_boundary.xml");
	//std::string filename_input("breps/filedata/rotor_bdry.xml");
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

	// make it robust 
	// to whatever the order of b-rep curves
	// gsWrite(bRep, "bRep.xml");

	GISMO_ENSURE(bRep.isClosed(), "The boundary is not closed, adjust tolerance.");
	bRep.closeGaps(tol);

	//gsWriteParaview(bRep, "bRep", 1000, true, true);

	//--------------------------------------------------------------------------------------------------


	/////////////////////////////////////// STEP 2: eliminate foldovers ////////////////////////////////

	gsParameterization2DExample<real_t> opt(bRep, method, plot_init);

	gsDebug << "Area of computational domain is "
		<< opt.computeArea() << "\n";

	opt.preComputeBasis();

	gsVector<real_t> initU = opt.convert_mp_to_gsvec();

	//gsInfo << "inital guess = " << initU << "\n";
	gsDebug << "value of objective function is: " << opt.evalObj(initU) << "\n";

	opt.X = opt.convert_mp_to_vec();

	//gsInfo << "initU1 (std::vector) = " << "\n";
	/*for (index_t i = 0; i != opt.X.size(); i++)
	{
		gsInfo << opt.X[i] << "\n";
	}*/

	//gsDebug << opt.X.size() << initU.size() << "\n";
	/*gsVector<real_t> grads(initU.size());
	opt.gradObj_into(initU, grads);
	gsDebug << "gradient = " << grads << "\n";
	gsDebug << "grad.size = " << grads.size() << "\n";*/

	//opt.compute();

	//for (index_t i = 0; i != initU.size(); i++)
		//opt.X[i] = initU(i);

	/*gsVector<real_t> resultNum(initU.size());
	opt.gradObj_into(initU, resultNum);

	real_t fVal = opt.evalObj(initU);

	double E_prev=0;
	std::vector<double> trash(opt.X.size());
	std::fill(trash.begin(), trash.end(), 0);
	opt.objInterFreeFunc_into(opt.X, E_prev, trash);

	gsDebug << fVal << "  " << E_prev << "\n";*/
	//gsDebug << initU.size() << "  " << opt.X.size() << "\n";

	//gsDebug << resultNum.size() << "  " << trash.size() << "\n";

	/*for (index_t i = 0; i != resultNum.size(); i++)
		gsDebug << resultNum(i) << "  " << trash[i] << "\n";*/

	auto t1 = std::chrono::high_resolution_clock::now();
	//bool success = opt.compute();
	opt.compute();
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time = t2 - t1;
	gsDebug << "running time : " << time.count() << "\n";

	/*gsInfo << "inital guess = " << initU << "\n";
	gsInfo << "optimal X = " << "\n";
	for (index_t i = 0; i != opt.X.size(); i++)
		gsInfo << opt.X[i] << "\n";*/

	/*for (index_t i = 0; i != opt.X.size(); i++)
		gsInfo << opt.X[i]-initU(i) << "\n";*/

	/*if (success)
		std::cerr << "SUCCESS; running time: " << time.count() << " s; min det J = " << opt.detmin << std::endl;
	else
		std::cerr << "FAIL TO UNTANGLE!" << std::endl;*/


	//--------------------------------------------------------------------------------------------------

	///////////////////// STEP 3: improving the quality of parameterization ///////////////////////////
	//// Another optimization problem

	//--------------------------------------------------------------------------------------------------

    ///////////////////////////// STEP 4: output the resulting parameterization ////////////////////////
	// writing the resulting parameterization to a G+Smo .xml file
	// filename_output is a string. The extention .xml is added automatically

	// TODO: other formats? To make it easy to visulization

	gsMultiPatch<> mp = opt.outputResult();

	gsFileData<> fd;
	fd << mp;
	fd.save(filename_output);
	gsInfo << "Wrote G+Smo file:     " << filename_output << ".xml \n";

	gsWriteParaview(mp, "mp", 1000, true, true);

	//gsWrite(mp, filename_output);
	//--------------------------------------------------------------------------------------------------


	////////////////////////////////////// STEP 5: VISUALIZATION ///////////////////////////////////////
	// TODO: visualization, in MATLAB first? then ?

	// GNUPLOT, looks good; a MATLAB plot style

	//--------------------------------------------------------------------------------------------------

	system("pause");
	return EXIT_SUCCESS;
}

