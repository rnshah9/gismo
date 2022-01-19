/** @file gsApproxC1Edge.hpp

    @brief Creates the approx C1 space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/gsApproxC1Edge.h>

#include <gsUnstructuredSplines/gsPatchReparameterized.h>

#include <gsUnstructuredSplines/gsApproxGluingData.h>

namespace gismo
{
    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::compute(std::vector<patchSide> & sidesContainer) {

        // Compute GLuing data
        gsApproxGluingData<d, T> approxGluingData(m_auxPatches, m_optionList, sidesContainer);
/*
        gsMultiBasis<T> basis;
        index_t dir;
        if (sidesContainer.size() == 2)
        {
            gsMultiBasis<T> basis_11(m_auxPatches[0].getBasisRotated().piece(9));
            gsMultiBasis<T> basis_22(m_auxPatches[1].getBasisRotated().piece(9));
            index_t dir_1 = 1;
            index_t dir_2 = 0;

            if (basis_11.basis(0).component(dir_1).numElements() > basis_22.basis(0).component(dir_2).numElements())
            {
                basis = basis_22;
                dir = dir_2;
            }
            else
            {
                basis = basis_11;
                dir = dir_1;
            }
        }
        else
        {
            gsMultiBasis<T> basis_11(m_auxPatches[0].getBasisRotated().piece(9));
            basis = basis_11;
            dir = 1;
        }
*/
        //! [Problem setup]
        basisEdgeResult.clear();
        for (size_t patchID = 0; patchID < sidesContainer.size(); patchID++) {
            gsMultiPatch<T> result;

            index_t dir = patchID == 0 ? 1 : 0;

            gsBSplineBasis<T> basis_plus, basis_minus;

            gsMultiBasis<T> initSpace(m_auxPatches[patchID].getBasisRotated().piece(9));
            createPlusSpace(m_auxPatches[0].getPatchRotated(), initSpace.basis(0), dir, basis_plus);
            createMinusSpace(m_auxPatches[0].getPatchRotated(), initSpace.basis(0), dir, basis_minus);

            gsGeometry<T> &geo = m_auxPatches[patchID].getPatchRotated();

            gsBSpline<T> beta, alpha;
            bool bdy = true;
            if (sidesContainer.size() == 2)
            {
                bdy = false;
                beta = approxGluingData.betaS(dir);
                alpha = approxGluingData.alphaS(dir);
            }

            // [!The same setup for each bf!]
            gsSparseSolver<real_t>::SimplicialLDLT solver;
            gsExprAssembler<> A(1, 1);

            // Elements used for numerical integration
            gsMultiBasis<T> edgeSpace(
                    m_auxPatches[patchID].getBasisRotated().piece(sidesContainer[patchID]));

            A.setIntegrationElements(edgeSpace);
            gsExprEvaluator<> ev(A);

            // Set the discretization space
            auto u = A.getSpace(edgeSpace);

            // Create Mapper
            gsDofMapper map(edgeSpace);
            gsMatrix<index_t> act;
            for (index_t i = 2; i < edgeSpace[0].component(1 - dir).size();
                 i++) // only the first two u/v-columns are Dofs (0/1)
            {
                act = edgeSpace[0].boundaryOffset(dir == 0 ? 3 : 1, i); // WEST
                map.markBoundary(0, act); // Patch 0
            }
            map.finalize();

            u.setupMapper(map);

            gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
            fixedDofs.setZero( u.mapper().boundarySize(), 1 );

            A.initSystem();
            A.assemble(u * u.tr()); // The Matrix is the same for each bf
            solver.compute(A.matrix());
            // [!The same setup for each bf!]

            index_t n_plus = basis_plus.size();
            index_t n_minus = basis_minus.size();

            index_t bfID_init = 3;
            for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
            {
                A.initVector(); // Just the rhs

                gsTraceBasis<real_t> traceBasis(geo, beta, basis_plus, initSpace.basis(0), bdy, bfID, dir);
                auto aa = A.getCoeff(traceBasis);

                A.assemble(u * aa);

                gsMatrix<> solVector = solver.solve(A.rhs());

                auto u_sol = A.getSolution(u, solVector);
                gsMatrix<> sol;
                u_sol.extract(sol);

                result.addPatch(edgeSpace.basis(0).makeGeometry(give(sol)));
            }

            bfID_init = 2;
            for (index_t bfID = bfID_init; bfID < n_minus - bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
            {
                A.initVector(); // Just the rhs

                gsNormalDerivBasis<real_t> normalDerivBasis(geo, alpha, basis_minus, initSpace.basis(0), bdy, bfID, dir);
                auto aa = A.getCoeff(normalDerivBasis);

                A.assemble(u * aa);

                gsMatrix<> solVector = solver.solve(A.rhs());

                auto u_sol = A.getSolution(u, solVector);
                gsMatrix<> sol;
                u_sol.extract(sol);

                result.addPatch(edgeSpace.basis(0).makeGeometry(give(sol)));
            }

            // parametrizeBasisBack
            m_auxPatches[patchID].parametrizeBasisBack(result);

            basisEdgeResult.push_back(result);
        }
}


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::computeAuxTopology()
    {
        for(unsigned i = 0; i <  m_auxPatches.size(); i++)
        {
            if(m_auxPatches[i].getPatchRotated().orientation() == -1)
                m_auxPatches[i].swapAxis();
        }
    }


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::reparametrizeInterfacePatches()
    {
        computeAuxTopology();

        gsMultiPatch<> temp_mp;
        for(unsigned i = 0; i <  m_auxPatches.size(); i++)
            temp_mp.addPatch(m_auxPatches[i].getPatchRotated());

        temp_mp.computeTopology();

        // Right patch along the interface. Patch 0 -> v coordinate. Edge west along interface
        switch (temp_mp.interfaces()[0].second().side().index())
        {
            case 1:
                //gsInfo << "Global patch: " << patch_2 << "\tLocal patch: " << temp_mp.interfaces()[0].second().patch << " not rotated\n";
                break;
            case 4: m_auxPatches[0].rotateParamClock();
                //gsInfo << "Global patch: " << patch_2 <<"\tLocal patch: " << temp_mp.interfaces()[0].second().patch << " rotated clockwise\n";
                break;
            case 3: m_auxPatches[0].rotateParamAntiClock();
                //gsInfo << "Global patch: " << patch_2 <<"\tLocal patch: " << temp_mp.interfaces()[0].second().patch << " rotated anticlockwise\n";
                break;
            case 2: m_auxPatches[0].rotateParamAntiClockTwice();
                //gsInfo << "Global patch: " << patch_2 <<"\tLocal patch: " << temp_mp.interfaces()[0].second().patch << " rotated twice anticlockwise\n";
                break;
            default:
                break;
        }

        // Left patch along the interface. Patch 1 -> u coordinate. Edge south along interface
        switch (temp_mp.interfaces()[0].first().side().index())
        {
            case 3:
                //gsInfo << "Global patch: " << patch_1 <<"\tLocal patch: " << temp_mp.interfaces()[0].first().patch << " not rotated\n";
                break;
            case 4: m_auxPatches[1].rotateParamAntiClockTwice();
                //gsInfo << "Global patch: " << patch_1 <<"\tLocal patch: " << temp_mp.interfaces()[0].first().patch << " rotated twice anticlockwise\n";
                break;
            case 2: m_auxPatches[1].rotateParamAntiClock();
                //gsInfo << "Global patch: " << patch_1 <<"\tLocal patch: " << temp_mp.interfaces()[0].first().patch << " rotated anticlockwise\n";
                break;
            case 1: m_auxPatches[1].rotateParamClock();
                //gsInfo << "Global patch: " << patch_1 <<"\tLocal patch: " << temp_mp.interfaces()[0].first().patch << " rotated clockwise\n";
                break;
            default:
                break;
        }
    } // reparametrizeInterfacePatches


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::reparametrizeSinglePatch(index_t side)
    {
        computeAuxTopology();

        if(m_auxPatches[0].getOrient())
        {
            switch (side)
            {
                case 3:
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " not rotated\n";
                    break;
                case 2:
                    m_auxPatches[0].rotateParamClock();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated clockwise\n";
                    break;
                case 4:
                    m_auxPatches[0].rotateParamAntiClockTwice();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated twice anticlockwise\n";
                    break;
                case 1:
                    m_auxPatches[0].rotateParamAntiClock();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated anticlockwise\n";
                    break;
            }
        }
        else
        {
            switch (side)
            {
                case 1:
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " not rotated\n";
                    break;
                case 4:
                    m_auxPatches[0].rotateParamClock();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated clockwise\n";
                    break;
                case 2:
                    m_auxPatches[0].rotateParamAntiClockTwice();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated twice anticlockwise\n";
                    break;
                case 3:
                    m_auxPatches[0].rotateParamAntiClock();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated anticlockwise\n";
                    break;
            }
        }
    } // reparametrizeSinglePatch

} // namespace gismo
