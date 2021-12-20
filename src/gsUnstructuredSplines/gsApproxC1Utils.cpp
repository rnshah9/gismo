/** @file gsApproxC1Utils.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#include <gsNurbs/gsTensorBSplineBasis.h>

#include <gsUnstructuredSplines/gsApproxC1Utils.h>


namespace gismo
{
    void createGluingDataSpace(const gsGeometry<real_t> & patch, const gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & result)
    {
        const gsBSplineBasis<real_t> basis_1 = dynamic_cast<const gsBSplineBasis<real_t> &>(basis.component(dir));

        index_t p_tilde = basis_1.degree() - 1;
        index_t r_tilde = p_tilde - 1;
        gsKnotVector<real_t> kv_gluingData(basis_1.knots().unique(), p_tilde, r_tilde);

        result = gsBSplineBasis<real_t>(kv_gluingData); // S(\tilde{p},\tilde{r},h)
/*
        // Geometry check
        gsBSplineBasis<real_t> basis_geo_1 = dynamic_cast<const gsBSplineBasis<real_t> &>(patch.basis().component(dir));
        index_t p_geo = basis_geo_1.degree();
        std::vector<index_t> vec_1 = basis_geo_1.knots().multiplicities();
        gsAsVector<index_t> mult_1(vec_1);
        for (index_t i = 1; i < mult_1.size()-1; i++) // First and last are ignored
        {
            real_t knot = basis_geo_1.knots().uValue(i);
            index_t r_geo = p_geo - mult_1(i);
            if (r_tilde > r_geo)
                result.insertKnot(knot,r_tilde-r_geo);
        }
*/
    }

    void createPlusSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & result)
    {
        gsBSplineBasis<real_t> basis_1 = dynamic_cast<gsBSplineBasis<real_t> &>(basis.component(dir));

        index_t m, p;
        p = basis_1.degree();
        m = basis_1.knots().multiplicityIndex(p+1);

        result = gsBSplineBasis<real_t>(basis_1);

        if (m != 1)
            result.elevateContinuity(1);

/*
        gsKnotVector<real_t> basis_knots1 = basis_1.knots();
        gsKnotVector<real_t> geo_knots1 = patch.basis().knots(dir);
        reduceContinuityPreserveGeo(basis_knots1, geo_knots1);
        result = gsBSplineBasis<real_t>(basis_knots1);
*/
    }

    void createMinusSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & result)
    {
        gsBSplineBasis<real_t> basis_1 = dynamic_cast<gsBSplineBasis<real_t> &>(basis.component(dir));

        index_t m, p;
        p = basis_1.degree();
        m = basis_1.knots().multiplicityIndex(p+1);

        result = gsBSplineBasis<real_t>(basis_1);

        result.degreeDecrease(1);
        if (m != 1)
            result.reduceContinuity(1);
/*
        gsKnotVector<real_t> basis_knots1 = basis_1.knots();
        basis_knots1.degreeDecrease(1);
        gsKnotVector<real_t> geo_knots1 = patch.basis().knots(dir);
        reduceContinuityPreserveGeo(basis_knots1, geo_knots1, 1, 1);
        result = gsBSplineBasis<real_t>(basis_knots1);
*/
    }

    void createEdgeSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & basis_plus,
                        gsBSplineBasis<real_t> & basis_minus, gsBSplineBasis<real_t> & basis_gluingData,
                        gsTensorBSplineBasis<2, real_t> & result)
    {
        gsTensorBSplineBasis<2, real_t> basis_edge = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(basis);

        basis_edge.component(dir).setDegreePreservingMultiplicity(basis_plus.degree()+basis_gluingData.degree()-1);

        index_t r_plus, r_minus, r_edge, r_gD, r;
        r_plus = basis_plus.degree() - basis_plus.knots().multiplicityIndex(basis_plus.degree()+1); // p+1, since c++ starts at 0
        r_minus = basis_minus.degree() - basis_minus.knots().multiplicityIndex(basis_minus.degree()+1);
        r_gD = basis_gluingData.degree() - basis_gluingData.knots().multiplicityIndex(basis_gluingData.degree()+1);
        r_edge = basis_edge.degree(dir) - basis_edge.knots(dir).multiplicityIndex(basis_edge.degree(dir)+1);

        r = math::min(r_gD, math::min(r_plus, r_minus));
        if (r_edge > r)
            basis_edge.component(dir).reduceContinuity(r_edge - r);
        else if (r_edge < r)
            basis_edge.component(dir).elevateContinuity(r - r_edge);

        result = gsTensorBSplineBasis<2, real_t>(basis_edge);
    }

    void createEdgeSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & basis_plus,
                        gsBSplineBasis<real_t> & basis_minus, gsTensorBSplineBasis<2, real_t> & result)
    {
        gsTensorBSplineBasis<2, real_t> basis_edge = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(basis);

        gsKnotVector<real_t> basis_edge_knots1 = basis_edge.knots(dir);
        gsKnotVector<real_t> basis_plus_knots = basis_plus.knots();
        gsKnotVector<real_t> basis_minus_knots = basis_minus.knots();
        
        basis_edge.component(dir).setDegreePreservingMultiplicity(basis_plus.degree());

        index_t r_plus, r_minus, r_edge, r;
        r_plus = basis_plus.degree() - basis_plus.knots().multiplicityIndex(basis_plus.degree()+1); // p+1, since c++ starts at 0
        r_minus = basis_minus.degree() - basis_minus.knots().multiplicityIndex(basis_minus.degree()+1);
        r_edge = basis_edge.degree(dir) - basis_edge.knots(dir).multiplicityIndex(basis_edge.degree(dir)+1);

        r = math::min(r_plus, r_minus);
        if (r_edge > r)
            basis_edge.component(dir).reduceContinuity(r_edge - r);
        else if (r_edge < r)
            basis_edge.component(dir).elevateContinuity(r - r_edge);

        result = gsTensorBSplineBasis<2, real_t>(basis_edge);

/*
        std::vector<real_t> knots_container;
        std::vector<index_t> mult_edge = basis_edge.knots(dir).multiplicities();
        //std::vector<index_t> mult_plus = basis_plus.knots().multiplicities();
        std::vector<index_t> mult_minus = basis_minus.knots().multiplicities();
        for (size_t i = 1; i < mult_minus.size()-1; i++)
        {
            if (mult_edge.at(i) < mult_minus.at(i))
            {
                real_t knot = basis_edge.knots(dir).unique().at(i);
                basis_edge.knots(dir).insert(knot, mult_minus.at(i)-mult_edge.at(i));
            }
        }
        result = gsTensorBSplineBasis<2, real_t>(basis_edge);
*/
    }

    // Corner Vertex
    void createVertexSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, bool isInterface_1, bool isInterface_2, gsTensorBSplineBasis<2, real_t> & result)
    {
        gsTensorBSplineBasis<2, real_t> basis_vertex = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(basis);

        for (index_t dir = 0; dir < basis.domainDim(); dir ++)
        {
            if (dir == 0 ? isInterface_1 : isInterface_2) // If edge is an interface
            {
                gsBSplineBasis<real_t> basis_gluingData, basis_plus, basis_minus;
                createGluingDataSpace(patch, basis, dir, basis_gluingData);
                createPlusSpace(patch, basis, dir, basis_plus);
                createMinusSpace(patch, basis, dir, basis_minus);

                basis_vertex.component(dir).setDegreePreservingMultiplicity(basis_plus.degree()+basis_gluingData.degree()-1);

                index_t r_plus, r_minus, r_edge, r_gD, r;
                r_plus = basis_plus.degree() - basis_plus.knots().multiplicityIndex(basis_plus.degree()+1); // p+1, since c++ starts at 0
                r_minus = basis_minus.degree() - basis_minus.knots().multiplicityIndex(basis_minus.degree()+1);
                r_gD = basis_gluingData.degree() - basis_gluingData.knots().multiplicityIndex(basis_gluingData.degree()+1);
                r_edge = basis_vertex.degree(dir) - basis_vertex.knots(dir).multiplicityIndex(basis_vertex.degree(dir)+1);

                r = math::min(r_gD, math::min(r_plus, r_minus));
                if (r_edge > r)
                    basis_vertex.component(dir).reduceContinuity(r_edge - r);
                else if (r_edge < r)
                    basis_vertex.component(dir).elevateContinuity(r - r_edge);
            }
            else
            {
                index_t r_12, p_12;
                p_12 = basis.degree(dir);
                r_12 = p_12 - basis.knots(dir).multiplicityIndex(p_12+1);
                if (r_12 == p_12 - 1) // == basis_vertex_1.degree(1)
                    basis_vertex.component(dir).reduceContinuity(1); // In the case for the max. smoothness
            }
        }

        result = gsTensorBSplineBasis<2, real_t>(basis_vertex);
    }

    void reduceContinuityPreserveGeo(gsKnotVector<real_t> & basis_knots, gsKnotVector<real_t> & geo_knots,
                                     index_t mult_basis_reduce, index_t mult_geo_reduce)
    {
        index_t p = basis_knots.degree();

        size_t j = 1; // iterator for geo
        for (size_t i = 1; i < basis_knots.multiplicities().size()-1; ++i) // for all interior knots
        {
            index_t m = basis_knots.multiplicities().at(i);
            if ( m <= p - mult_basis_reduce )
            {
                real_t knot = basis_knots.unique().at(i);
                if (knot == geo_knots.unique().at(j) && j<geo_knots.unique().size()-1)
                {
                    index_t p_geo = geo_knots.degree();
                    index_t m_geo = geo_knots.multiplicities().at(j);
                    if (p-m > p_geo-m_geo-mult_geo_reduce)
                        basis_knots.insert(knot, (p-m)-(p_geo-m_geo-mult_geo_reduce));
                    j++;
                }
                else
                    basis_knots.insert(knot, mult_basis_reduce); // Reduce multiplicity
            }
        }
    }

}



