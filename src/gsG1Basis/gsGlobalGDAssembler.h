/** @file gsGluingDataAssembler.h

    @brief Provides assembler for the gluing data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

#include <gsG1Basis/gsVisitorGlobalGD.h>

namespace gismo
{

template <class T, class bhVisitor = gsVisitorGlobalGD<T> >
class gsGlobalGDAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsGlobalGDAssembler(gsBasis<T> const & basis,
                           gsMultiPatch<T> const & mp,
                           index_t const & gamma)
        : m_mp(mp), m_gamma(gamma)
    {

        m_basis.push_back(basis); // Basis for alpha and beta

        refresh();
    }

    void refresh();

    void assemble();

    inline void apply(bhVisitor & visitor,
                      int patchIndex = 0,
                      boxSide side = boundary::none);

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix_alpha_0() const { return m_system_alpha_0.matrix(); }
    const gsMatrix<T> & rhs_alpha_0() const { return m_system_alpha_0.rhs(); }

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix_alpha_1() const { return m_system_alpha_1.matrix(); }
    const gsMatrix<T> & rhs_alpha_1() const { return m_system_alpha_1.rhs(); }

    const gsSparseMatrix<T> & matrix_beta_0() const { return m_system_beta_0.matrix(); }
    const gsMatrix<T> & rhs_beta_0() const { return m_system_beta_0.rhs(); }

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix_beta_1() const { return m_system_beta_1.matrix(); }
    const gsMatrix<T> & rhs_beta_1() const { return m_system_beta_1.rhs(); }

protected:
    // interface + geometry for computing alpha and beta
    gsMultiPatch<T> m_mp;

    index_t m_gamma;

    // Space for phi_0,i, phi_1,j
    std::vector< gsMultiBasis<T> > m_basis;

    //using Base::m_system;
    using Base::m_ddof;

    gsSparseSystem<T> m_system_alpha_0;
    gsSparseSystem<T> m_system_alpha_1;
    gsSparseSystem<T> m_system_beta_0;
    gsSparseSystem<T> m_system_beta_1;

}; // class gsG1BasisAssembler


template <class T, class bhVisitor>
void gsGlobalGDAssembler<T, bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map_alpha_L(m_basis[0]);
    gsDofMapper map_alpha_R(m_basis[0]);

    gsDofMapper map_beta_L(m_basis[0]);
    gsDofMapper map_beta_R(m_basis[0]);

    map_alpha_L.finalize();
    map_alpha_R.finalize();
    map_beta_L.finalize();
    map_beta_R.finalize();

    // 2. Create the sparse system
    m_system_alpha_0 = gsSparseSystem<T>(map_alpha_L);
    m_system_alpha_1 = gsSparseSystem<T>(map_alpha_R);
    m_system_beta_0 = gsSparseSystem<T>(map_beta_L);
    m_system_beta_1 = gsSparseSystem<T>(map_beta_R);

} // refresh()

template <class T, class bhVisitor>
void gsGlobalGDAssembler<T, bhVisitor>::assemble()
{
    //GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis[0][0],2,1,0.333333);
    m_system_alpha_0.reserve(nz, 1);
    m_system_alpha_1.reserve(nz, 1);
    m_system_beta_0.reserve(nz, 1);
    m_system_beta_1.reserve(nz, 1);


    if(m_ddof.size()==0)
        m_ddof.resize(4); // One for L, one for R, x2 for beta, alpha

    const gsDofMapper & mapper_L = m_system_alpha_0.colMapper(0);
    const gsDofMapper & mapper_R = m_system_alpha_1.colMapper(0);

    m_ddof[0].setZero(mapper_L.boundarySize(), 1 ); // tilde
    m_ddof[1].setZero(mapper_R.boundarySize(), 1 ); // bar
    m_ddof[2].setZero(mapper_L.boundarySize(), 1 ); // tilde
    m_ddof[3].setZero(mapper_R.boundarySize(), 1 ); // bar

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor);

    m_system_alpha_0.matrix().makeCompressed();
    m_system_alpha_1.matrix().makeCompressed();
    m_system_beta_0.matrix().makeCompressed();
    m_system_beta_1.matrix().makeCompressed();

}

template <class T, class bhVisitor>
inline void gsGlobalGDAssembler<T, bhVisitor>::apply(bhVisitor & visitor,
                                                        int patchIndex,
                                                        boxSide side)
{
#pragma omp parallel
    {
        bhVisitor
#ifdef _OPENMP
        // Create thread-private visitor
        visitor_(visitor);
        const int tid = omp_get_thread_num();
        const int nt  = omp_get_num_threads();
#else
            &visitor_ = visitor;
#endif

        gsQuadRule<T> quRule ; // Quadrature rule
        gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
        gsVector<T> quWeights; // Temp variable for mapped weights

        gsBasis<T> & basis = m_basis[0].basis(patchIndex); // = 0

        // Initialize reference quadrature rule and visitor data
        visitor_.initialize(basis,quRule);

        //const gsGeometry<T> & patch = m_geo.patch(patchIndex); // 0 = patchindex

        // Initialize domain element iterator -- using unknown 0
        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(boundary::none);

#ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
        for (; domIt->good(); domIt->next() )
#endif
        {

            // Map the Quadrature rule to the element
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(basis, quNodes, m_mp, m_gamma);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_system_alpha_0, m_system_alpha_1,
                                   m_system_beta_0, m_system_beta_1); // omp_locks inside

        }

    }//omp parallel
}


} // namespace gismo