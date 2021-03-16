/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

# include <gsArgyris/gsGluingData/gsApproxGluingDataAssembler.h>


namespace gismo
{

template<short_t d, class T>
class gsApproxGluingData
{
private:
    typedef typename std::vector<gsC1ArgyrisAuxiliaryPatch<d,T>> ArgyrisAuxPatchContainer;

public:
    gsApproxGluingData()
    { }


    gsApproxGluingData(ArgyrisAuxPatchContainer const & auxPatchContainer,
                       gsOptionList const & optionList,
                       std::vector<index_t> sidesContainer = std::vector<index_t>{})
        : m_auxPatches(auxPatchContainer), m_optionList(optionList)
    {

        if (m_auxPatches.size() == 2) // Interface
        {
            setGlobalGluingData(0,m_auxPatches[0].side(), 1); // v Order is important!!!
            setGlobalGluingData(1,m_auxPatches[1].side(), 0); // u
        }
        else if (m_auxPatches.size() == 1) // Vertex
        {

            for (size_t dir = 0; dir < sidesContainer.size(); dir++)
            {
                if(auxPatchContainer[0].getArygrisBasisRotated().isInterface(sidesContainer[dir])) // West
                    setGlobalGluingData(0, sidesContainer[dir], dir);
                else
                {
                    gsBSpline<> alpha;
                    gsBSpline<> beta;
                    alphaSContainer.push_back(alpha);
                    betaSContainer.push_back(beta);
                }
            }
        }
        else
            gsInfo << "Something went wrong \n";

    }

    // Computed the gluing data globally
    void setGlobalGluingData(index_t patchID = 0,  index_t side = 1, index_t dir = 1);

    gsBSpline<T> & alphaS(index_t patchID) { return alphaSContainer[patchID]; }
    gsBSpline<T> & betaS(index_t patchID) { return betaSContainer[patchID]; }

protected:

    // Spline space for the gluing data + multiPatch
    ArgyrisAuxPatchContainer m_auxPatches;

    const gsOptionList m_optionList;

    // Result
    std::vector<gsBSpline<T>> alphaSContainer, betaSContainer;

}; // class gsGluingData


template<short_t d, class T>
void gsApproxGluingData<d, T>::setGlobalGluingData(index_t patchID, index_t side, index_t dir)
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsBSplineBasis<T> bsp_gD = m_auxPatches[patchID].getArygrisBasisRotated().getBasisGluingData(side);

    gsApproxGluingDataAssembler<T> approxGluingDataAssembler(m_auxPatches[patchID].getPatch(), bsp_gD, dir, m_optionList);
    alphaSContainer.push_back(approxGluingDataAssembler.getAlphaS());
    betaSContainer.push_back(approxGluingDataAssembler.getBetaS());

    if (patchID == 0)
        gsWriteParaview(approxGluingDataAssembler.getAlphaS(), "alpha_R", 1000);
    if (patchID == 1)
        gsWriteParaview(approxGluingDataAssembler.getAlphaS(), "alpha_L", 1000);

    if (patchID == 0)
        gsWriteParaview(approxGluingDataAssembler.getBetaS(), "beta_R", 1000);
    if (patchID == 1)
        gsWriteParaview(approxGluingDataAssembler.getBetaS(), "beta_L", 1000);

} // setGlobalGluingData


} // namespace gismo

