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
    void createGluingDataSpace(gsTensorBSplineBasis<2, real_t> basis, index_t dir, gsBSplineBasis<real_t> & result)
    {
        gsBSplineBasis<real_t> basis_1 = dynamic_cast<gsBSplineBasis<real_t> &>(basis.component(dir));

        index_t p_tilde = basis_1.degree() - 1;
        index_t r_tilde = p_tilde - 1;
        gsKnotVector<real_t> kv_gluingData(basis_1.knots().unique(), p_tilde, r_tilde);

        result = gsBSplineBasis<real_t>(kv_gluingData); // S(\tilde{p},\tilde{r},h)
    }

    void createPlusSpace(gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & result)
    {
        gsBSplineBasis<real_t> basis_1 = dynamic_cast<gsBSplineBasis<real_t> &>(basis.component(dir));

        index_t m, p;
        p = basis_1.degree();
        m = basis_1.knots().multiplicityIndex(p+1);

        result = gsBSplineBasis<real_t>(basis_1);

        if (m != 1)
            result.elevateContinuity(1);

    }

    void createMinusSpace(gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & result)
    {
        gsBSplineBasis<real_t> basis_1 = dynamic_cast<gsBSplineBasis<real_t> &>(basis.component(dir));

        index_t m, p;
        p = basis_1.degree();
        m = basis_1.knots().multiplicityIndex(p+1);

        result = gsBSplineBasis<real_t>(basis_1);

        result.degreeDecrease(1);
        if (m != 1)
            result.reduceContinuity(1);

    }

}



