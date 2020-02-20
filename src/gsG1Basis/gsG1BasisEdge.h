/** @file gsG1BasisEdge.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include <gsG1Basis/gsGluingData.h>
#include <gsG1Basis/gsVisitorG1BasisEdge.h>


namespace gismo
{
template<class T, class bhVisitor = gsVisitorG1BasisEdge<T>>
class gsG1BasisEdge : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsG1BasisEdge(gsMultiPatch<T> const & mp,
                 gsMultiBasis<T> & mb,
                 gsOptionList & optionList)
        : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {
        // Computing the gluing data
        gsGluingData<T> gluingData(m_mp,m_mb,m_optionList);
        m_gD = gluingData;

        // Computing the G1 - basis function at the edge
        // Spaces for computing the g1 basis
        index_t m_r = m_optionList.getInt("regularity");

        gsBSplineBasis<> basis_1 = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(0).component(1)); // 0 -> v
        gsBSplineBasis<> basis_2 = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(1).component(0)); // 1 -> u

        index_t m_p; // Minimum degree at the interface
        if (basis_1.degree() >= basis_2.degree())
            m_p = basis_2.maxDegree();
        else
            m_p = basis_1.maxDegree();

        // first,last,interior,mult_ends,mult_interior
        gsKnotVector<T> kv_plus(0,1,0,m_p+1,m_p-1-m_r); // p,r+1 //-1 bc r+1
        gsBSplineBasis<> basis_plus(kv_plus);

        if (basis_1.numElements() <= basis_2.numElements()) //
            for (size_t i = basis_1.degree()+1; i < basis_1.knots().size() - (basis_1.degree()+1); i = i+(basis_1.degree()-m_r))
                basis_plus.insertKnot(basis_1.knot(i),m_p-1-m_r);
        else
            for (size_t i = basis_2.degree()+1; i < basis_2.knots().size() - (basis_2.degree()+1); i = i+(basis_2.degree()-m_r))
                basis_plus.insertKnot(basis_2.knot(i),m_p-1-m_r);

        m_basis_plus = basis_plus;
        n_plus = m_basis_plus.size();

        gsKnotVector<T> kv_minus(0,1,0,m_p+1-1,m_p-1-m_r); // p-1,r //-1 bc p-1
        gsBSplineBasis<> basis_minus(kv_minus);

        if (basis_1.numElements() <= basis_2.numElements()) //
            for (size_t i = basis_1.degree()+1; i < basis_1.knots().size() - (basis_1.degree()+1); i = i+(basis_1.degree()-m_r))
                basis_minus.insertKnot(basis_1.knot(i),m_p-1-m_r);
        else
            for (size_t i = basis_2.degree()+1; i < basis_2.knots().size() - (basis_2.degree()+1); i = i+(basis_2.degree()-m_r))
                basis_minus.insertKnot(basis_2.knot(i),m_p-1-m_r);

        m_basis_minus = basis_minus;
        n_minus = m_basis_minus.size();

        // Basis for the G1 basis
        m_basis_g1 = m_mb;

        refresh();
        assemble();
        solve();
    }


    void refresh();
    void assemble();
    inline void apply(bhVisitor & visitor, int patchIndex, boxSide side = boundary::none);
    void solve();

    void constructSolution(gsMultiPatch<T> & result_L, gsMultiPatch<T> & result_R);

    index_t get_n_plus() { return n_plus; }
    index_t get_n_minus() { return n_minus; }

    void plotG1Basis(gsMultiPatch<T> & basisG1_L, gsMultiPatch<T> & basisG1_R, std::string baseName)
    {

        const std::string baseName1(baseName + "_0");
        gsParaviewCollection collection1(baseName1);

        const std::string baseName2(baseName + "_1");
        gsParaviewCollection collection2(baseName2);

        std::string fileName, fileName2;
        for (unsigned i = 0; i < basisG1_L.nPatches(); i++)
        {

            fileName = baseName1 + "_" + util::to_string(i);

            gsField<> temp_field_L(m_mp.patch(0),basisG1_L.patch(i));
            gsWriteParaview(temp_field_L,fileName,5000);
            collection1.addTimestep(fileName,i,"0.vts");

        }
        for (unsigned i = 0; i < basisG1_R.nPatches(); i++)
        {

            fileName2 = baseName2 + "_" + util::to_string(i);

            gsField<> temp_field_R(m_mp.patch(1),basisG1_R.patch(i));
            gsWriteParaview(temp_field_R,fileName2,5000);
            collection2.addTimestep(fileName2,i,"0.vts");

        }
        collection1.save();
        collection2.save();
    }

    void g1Condition()
    {
        gsMatrix<> points(1,1000);
        points.setRandom();
        points = points.array().abs();

        gsMatrix<> points2d_L(2, 1000);
        gsMatrix<> points2d_R(2, 1000);

        points2d_L.setZero();
        points2d_R.setZero();
        points2d_L.row(1) = points; // v
        points2d_R.row(0) = points; // u

        real_t g1Error = 0;
        for (size_t i = 0; i < g1Basis_L.nPatches(); i++)
        {
            gsMatrix<> temp;
            temp = m_gD.get_alpha_tilde_1().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).topRows(1))
                + m_gD.get_alpha_tilde_0().eval(points).cwiseProduct(g1Basis_R.patch(i).deriv(points2d_R).bottomRows(1))
                + m_gD.get_beta_bar().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).bottomRows(1));

            if (temp.array().abs().maxCoeff() > g1Error)
                g1Error = temp.array().abs().maxCoeff();
        }

        gsInfo << "Conditiontest G1 continuity: \n" << g1Error << "\n\n";


        //gsInfo << "\nConditiontest G1 continuity: \n" << g1Basis_L.patch(0).coefs() << "\n\n";

    }


protected:

    // Input
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsOptionList m_optionList;

    // Gluing data
    gsGluingData<T> m_gD;

    // Basis for getting the G1 Basis
    gsBSplineBasis<> m_basis_plus;
    gsBSplineBasis<> m_basis_minus;

    // Basis for the G1 Basis
    gsMultiBasis<> m_basis_g1;

    // Size of the basis
    index_t n_plus, n_minus;

    // System
    std::vector<gsSparseSystem<T> > m_f_0, m_f_1;
    // For Dirichlet boundary
    using Base::m_ddof;

    std::vector<gsMatrix<>> solVec_t, solVec_b;
    gsMultiPatch<T> g1Basis_L, g1Basis_R;

}; // class gsG1BasisEdge


template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::constructSolution(gsMultiPatch<T> & result_L, gsMultiPatch<T> & result_R)
{

    result_L.clear();
    result_R.clear();

    // Dim is the same for all basis functions
    const index_t dim = ( 0!=solVec_t.at(0).cols() ? solVec_t.at(0).cols() :  m_ddof[0].cols() );

    gsMatrix<T> coeffs;
    for (size_t patchNum = 0; patchNum < 2; patchNum++)
    {
        for (index_t p = 0; p < m_basis_plus.size(); ++p)
        {
            const gsDofMapper & mapper = m_f_0.at(p).colMapper(0); // unknown = 0

            // Reconstruct solution coefficients on patch p
            index_t sz;
            sz = m_basis_g1.basis(patchNum).size();


            coeffs.resize(sz, dim);

            for (index_t i = 0; i < sz; ++i)
            {
                if (mapper.is_free(i, patchNum)) // DoF value is in the solVector // 0 = unitPatch
                {
                    coeffs.row(i) = solVec_t.at(p).row(mapper.index(i, patchNum));
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                    coeffs.row(i) = m_ddof[0].row( mapper.bindex(i, patchNum) ).head(dim); // = 0
                }
            }

            if (patchNum == 0)
                result_L.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));
            if (patchNum == 1)
                result_R.addPatch(m_basis_g1.basis(1).makeGeometry(give(coeffs)));
        }

        for (index_t p = 0; p < m_basis_minus.size(); ++p)
        {
            const gsDofMapper & mapper = m_f_1.at(p).colMapper(0); // unknown = 0

            // Reconstruct solution coefficients on patch p
            index_t sz;
            sz = m_basis_g1.basis(patchNum).size();

            coeffs.resize(sz, dim);

            for (index_t i = 0; i < sz; ++i)
            {
                if (mapper.is_free(i, patchNum)) // DoF value is in the solVector // 0 = unitPatch
                {
                    //gsInfo << "mapper index: " << mapper.index(i, p) << "\n";
                    coeffs.row(i) = solVec_b.at(p).row(mapper.index(i, patchNum));
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                    coeffs.row(i) = m_ddof[0].row( mapper.bindex(i, patchNum) ).head(dim); // = 0
                }
            }

            if (patchNum == 0)
                result_L.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));
            if (patchNum == 1)
                result_R.addPatch(m_basis_g1.basis(1).makeGeometry(give(coeffs)));

        }
    }

    g1Basis_L = result_L;
    g1Basis_R = result_R;

    //g1Condition();
}

template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsBasis<> & B_0 = m_mb.basis(0); // Patch 0
    gsBasis<> & B_1 = m_mb.basis(1); // Patch 1

    gsVector<index_t> sz(2);
    sz(0) = B_0.size();
    sz(1) = B_1.size();

    gsDofMapper map(sz);

    gsMatrix<unsigned> act1, act0;

    for (index_t i = 2; i < B_0.component(0).size(); i++) // only the first two u-columns are Dofs
    {
        act0 = B_0.boundaryOffset(1, i); // WEST
        map.markBoundary(0, act0); // Patch 0
    }
    for (index_t i = 2; i < B_1.component(1).size(); i++) // only the first two v-columns are Dofs
    {
        act1 = B_1.boundaryOffset(3, i); // SOUTH
        map.markBoundary(1, act1); // Patch 1
    }

    map.finalize();
    //gsInfo << "map : " << map.asVector() << "\n";
    //map.print();

    // 2. Create the sparse system
    gsSparseSystem<T> m_system = gsSparseSystem<T>(map);
    for (index_t i = 0; i < m_basis_plus.size(); i++)
        m_f_0.push_back(m_system);
    for (index_t i = 0; i < m_basis_minus.size(); i++)
        m_f_1.push_back(m_system);


} // refresh()

template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::assemble()
{
    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_mb[0],2,1,0.333333);
    for (index_t i = 0; i < m_basis_plus.size(); i++)
        m_f_0.at(i).reserve(nz, 1);

    for (index_t i = 0; i < m_basis_minus.size(); i++)
        m_f_1.at(i).reserve(nz, 1);

    if(m_ddof.size()==0)
        m_ddof.resize(2); // 0,1

    const gsDofMapper & map_0 = m_f_0.at(0).colMapper(0); // Map same for every 0
    const gsDofMapper & map_1 = m_f_1.at(0).colMapper(0); // Map same for every 1

    m_ddof[0].setZero(map_0.boundarySize(), 1 ); // plus
    m_ddof[1].setZero(map_1.boundarySize(), 1 ); // minus

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor,0); // patch 0
    apply(visitor,1); // patch 1

    for (index_t i = 0; i < m_basis_plus.size(); i++)
        m_f_0.at(i).matrix().makeCompressed();
    for (index_t i = 0; i < m_basis_minus.size(); i++)
        m_f_1.at(i).matrix().makeCompressed();

} // assemble()

template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::apply(bhVisitor & visitor, int patchIndex, boxSide side)
{
#pragma omp parallel
    {

        gsQuadRule<T> quRule ; // Quadrature rule
        gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
        gsVector<T> quWeights; // Temp variable for mapped weights

        bhVisitor
#ifdef _OPENMP
        // Create thread-private visitor
    visitor_(visitor);
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#else
            &visitor_ = visitor;
#endif

        gsBasis<T> & basis_g1 = m_basis_g1.basis(patchIndex); // basis for construction

        // Same for all patches
        gsBasis<T> & basis_geo = m_mb.basis(patchIndex).component(patchIndex);
        gsBasis<T> & basis_plus = m_basis_plus;
        gsBasis<T> & basis_minus = m_basis_minus;

        // Initialize reference quadrature rule and visitor data
        visitor_.initialize(basis_g1, basis_plus, basis_minus, quRule);

        const gsGeometry<T> & patch = m_mp.patch(patchIndex);

        // Initialize domain element iterator
        typename gsBasis<T>::domainIter domIt = basis_g1.makeDomainIterator(boundary::none);

#ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
        for (; domIt->good(); domIt->next() )
#endif
        {

            // Map the Quadrature rule to the element
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(basis_g1, basis_geo, basis_plus, basis_minus, patch, quNodes, m_gD, m_optionList);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_f_0, m_f_1); // omp_locks inside
        }
    }//omp parallel
} // apply

template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::solve()
{
    gsSparseSolver<real_t>::CGDiagonal solver;

    for (index_t i = 0; i < m_basis_plus.size(); i++) // Tilde
    {
        solver.compute(m_f_0.at(i).matrix());
        solVec_t.push_back(solver.solve(m_f_0.at(i).rhs()));
    }
    for (index_t i = 0; i < m_basis_minus.size(); i++)
    {
        solver.compute(m_f_1.at(i).matrix());
        solVec_b.push_back(solver.solve(m_f_1.at(i).rhs()));
    }
} // solve

} // namespace gismo