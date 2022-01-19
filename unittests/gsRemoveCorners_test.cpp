/** @file gsRemoveCorners_test.cpp

    @brief Tests for the removeCorners function in gsPatchPreconditionersCreator.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Schneckenleitner, S. Takacs
*/

#include "gismo_unittest.h"

SUITE(gsRemoveCorners_test)
{

    TEST(MassOp)
    {
        gsLinearOperator<>::uPtr massInvWithoutCorners;
        gsSparseMatrix<> massWithCorners;
        std::vector<index_t> corners;

        index_t numRefine = 2;
        index_t degree = 3;

        gsFunctionExpr<> bc("0.0",2);

        gsMultiPatch<> square(*gsNurbsCreator<>::BSplineSquare());
        square.computeTopology();

        gsMultiBasis<> mb(square);

        gsSparseEntries<> seReducedMatrix;
        gsSparseMatrix<> reducedMatrix;

        for (index_t i = 0; i < numRefine; ++i) {
            mb.uniformRefine();
        }
        mb[0].component(1).uniformRefine();

        mb[0].setDegreePreservingMultiplicity(degree);

        gsBoundaryConditions<> bcInfo;
        for (gsMultiPatch<real_t>::biterator it = square.bBegin();  it != square.bEnd(); ++it)
            bcInfo.addCondition(*it, condition_type::neumann, bc);

        gsOptionList opt = gsAssembler<>::defaultOptions();
        massWithCorners = gsPatchPreconditionersCreator<>::massMatrix(mb.basis(0), bcInfo, opt);

        corners.resize(3);
        for(index_t c = 2; c <= 1<<mb.dim(); c++) {
            corners[c-2] = mb[0].functionAtCorner(c);
        }

        index_t row = 0, col = 0;
        for (index_t i = 0; i < massWithCorners.rows(); i++)
        {
            if ( std::find(corners.begin(), corners.end(), i) == corners.end() )
            {
                for (index_t j = 0; j < massWithCorners.cols(); j++)
                {
                    if( std::find(corners.begin(), corners.end(), j) == corners.end() )
                    {
                        seReducedMatrix.add(row, col, massWithCorners(i,j));
                        col++;
                    }
                }
                col=0;
                row++;
            }
        }

        for(index_t i = 2; i <= 1<<mb.dim(); i++)
            bcInfo.addCornerValue(i, 0, 0);

        reducedMatrix.resize(massWithCorners.rows() - corners.size(), massWithCorners.cols()-corners.size());
        reducedMatrix.setFrom(seReducedMatrix);
        reducedMatrix.makeCompressed();

        massInvWithoutCorners = gsPatchPreconditionersCreator<>::massMatrixInvOp(mb.basis(0), bcInfo, opt);

        gsMatrix<> result;
        gsMatrix<> massInvMass(reducedMatrix.rows(), reducedMatrix.cols());
        for (index_t i = 0; i < reducedMatrix.cols(); ++i) {
            massInvWithoutCorners->apply(reducedMatrix.toDense().col(i), result);
            massInvMass.col(i) = result;
        }

        CHECK((massInvMass-gsMatrix<>::Identity(massInvMass.rows(),massInvMass.cols())).norm() < 1.e-11);
        CHECK(massInvWithoutCorners->rows() + (index_t)corners.size() == massWithCorners.rows() && massInvWithoutCorners->cols() + (index_t)corners.size() == massWithCorners.cols());
    }

    TEST(MassOp2)
    {
        gsLinearOperator<>::uPtr massInvWithoutCorners;
        gsSparseMatrix<> massWithCorners;
        std::vector<index_t> corners;

        index_t numRefine = 3;
        index_t degree = 3;

        gsFunctionExpr<> bc("0.0",2);

        gsMultiPatch<> square(*gsNurbsCreator<>::BSplineSquare());
        square.computeTopology();

        gsMultiBasis<> mb(square);

        for (index_t i = 0; i < numRefine; ++i) {
            mb.uniformRefine();
        }
        mb[0].setDegreePreservingMultiplicity(degree);

        gsBoundaryConditions<> bcInfo;
        for (gsMultiPatch<real_t>::biterator it = square.bBegin();  it != square.bEnd(); ++it)
            bcInfo.addCondition(*it, condition_type::dirichlet, bc);

        gsOptionList opt = gsAssembler<>::defaultOptions();
        massWithCorners = gsPatchPreconditionersCreator<>::massMatrix(mb.basis(0), bcInfo, opt);

        for(index_t i = 1; i <= 1<<mb.dim(); i++)
            bcInfo.addCornerValue(i, 0, 0);

        massInvWithoutCorners = gsPatchPreconditionersCreator<>::massMatrixInvOp(mb.basis(0), bcInfo, opt);

        gsMatrix<> result;
        gsMatrix<> massInvMass(massWithCorners.rows(), massWithCorners.cols());
        for (index_t i = 0; i < massWithCorners.cols(); ++i) {
            massInvWithoutCorners->apply(massWithCorners.toDense().col(i), result);
            massInvMass.col(i) = result;
        }

        CHECK((massInvMass-gsMatrix<>::Identity(massInvMass.rows(),massInvMass.cols())).norm() < 1.e-11);
        CHECK(massInvWithoutCorners->rows() == massWithCorners.rows() && massInvWithoutCorners->cols() == massWithCorners.cols());
    }

    TEST(StiffnessOp)
    {
        gsLinearOperator<>::uPtr stiffnessInvWithoutCorners;
        gsSparseMatrix<> stiffnessWithCorners;
        std::vector<index_t> corners;

        index_t numRefine = 3;
        index_t degree = 2;

        gsFunctionExpr<> bc("0.0",2);

        gsMultiPatch<> square(*gsNurbsCreator<>::BSplineSquare());
        square.computeTopology();

        gsMultiBasis<> mb(square);

        gsSparseEntries<> seReducedMatrix;
        gsSparseMatrix<> reducedMatrix;

        for (index_t i = 0; i < numRefine; ++i)
            mb.uniformRefine();

        mb[0].component(0).setDegreePreservingMultiplicity(degree);
        mb[0].component(1).setDegreePreservingMultiplicity(degree+1);

        gsBoundaryConditions<> bcInfo;
        bcInfo.addCondition(0, 1, condition_type::dirichlet, bc);
        bcInfo.addCondition(0, 2, condition_type::neumann, bc);
        bcInfo.addCondition(0, 3, condition_type::neumann, bc);
        bcInfo.addCondition(0, 4, condition_type::neumann, bc);

        gsOptionList opt = gsAssembler<>::defaultOptions();
        stiffnessWithCorners = gsPatchPreconditionersCreator<>::stiffnessMatrix(mb.basis(0), bcInfo, opt);

        corners.resize(2);
        corners[0] = mb[0].functionAtCorner(2) - 1;
        corners[1] = mb[0].functionAtCorner(4) - mb[0].boundary(1).rows();

        index_t row = 0, col = 0;
        for (index_t i = 0; i < stiffnessWithCorners.rows(); i++)
        {
            if ( std::find(corners.begin(), corners.end(), i) == corners.end() )
            {
                for (index_t j = 0; j < stiffnessWithCorners.cols(); j++)
                {
                    if( std::find(corners.begin(), corners.end(), j) == corners.end() )
                    {
                        seReducedMatrix.add(row, col, stiffnessWithCorners(i,j));
                        col++;
                    }
                }
                col=0;
                row++;
            }
        }

        reducedMatrix.resize(stiffnessWithCorners.rows()-corners.size(), stiffnessWithCorners.cols()-corners.size());
        reducedMatrix.setFrom(seReducedMatrix);
        reducedMatrix.makeCompressed();

        bcInfo.addCornerValue(2, 0, 0);
        bcInfo.addCornerValue(4, 0, 0);

        stiffnessInvWithoutCorners = gsPatchPreconditionersCreator<>::fastDiagonalizationOp(mb.basis(0), bcInfo, opt);

        gsMatrix<> result;
        gsMatrix<> stiffnessInvStiffness(reducedMatrix.rows(), reducedMatrix.cols());
        for (index_t i = 0; i < reducedMatrix.cols(); ++i) {
            stiffnessInvWithoutCorners->apply(reducedMatrix.toDense().col(i), result);
            stiffnessInvStiffness.col(i) = result;
        }

        CHECK((stiffnessInvStiffness-gsMatrix<>::Identity(stiffnessInvStiffness.rows(),stiffnessInvStiffness.cols())).norm() < 1.e-11);
        CHECK(stiffnessInvWithoutCorners->rows() + (index_t)corners.size() == stiffnessWithCorners.rows() && stiffnessInvWithoutCorners->cols() + (index_t)corners.size() == stiffnessWithCorners.cols());
    }

    TEST(StiffnessOp2)
    {
        gsLinearOperator<>::uPtr stiffnessInvWithoutCorners;
        gsSparseMatrix<> stiffnessWithCorners;
        std::vector<index_t> corners;

        index_t numRefine = 3;
        index_t degree = 5;

        gsFunctionExpr<> bc("0.0",2);

        gsMultiPatch<> square(*gsNurbsCreator<>::BSplineSquare());
        square.computeTopology();

        gsMultiBasis<> mb(square);

        gsSparseEntries<> seReducedMatrix;
        gsSparseMatrix<> reducedMatrix;

        for (index_t i = 0; i < numRefine; ++i)
            mb.uniformRefine();

        mb[0].component(0).uniformRefine();

        mb[0].component(0).setDegreePreservingMultiplicity(degree);
        mb[0].component(1).setDegreePreservingMultiplicity(degree+1);

        gsBoundaryConditions<> bcInfo;
        bcInfo.addCondition(0, 1, condition_type::dirichlet, bc);
        bcInfo.addCondition(0, 2, condition_type::neumann, bc);
        bcInfo.addCondition(0, 3, condition_type::dirichlet, bc);
        bcInfo.addCondition(0, 4, condition_type::neumann, bc);

        gsOptionList opt = gsAssembler<>::defaultOptions();
        stiffnessWithCorners = gsPatchPreconditionersCreator<>::stiffnessMatrix(mb.basis(0), bcInfo, opt);

        corners.resize(1);
        corners[0] = stiffnessWithCorners.rows()-1;

        index_t row = 0, col = 0;
        for (index_t i = 0; i < stiffnessWithCorners.rows(); i++)
        {
            if ( std::find(corners.begin(), corners.end(), i) == corners.end() )
            {
                for (index_t j = 0; j < stiffnessWithCorners.cols(); j++)
                {
                    if( std::find(corners.begin(), corners.end(), j) == corners.end() )
                    {
                        seReducedMatrix.add(row, col, stiffnessWithCorners(i,j));
                        col++;
                    }
                }
                col=0;
                row++;
            }
        }

        reducedMatrix.resize(stiffnessWithCorners.rows()-corners.size(), stiffnessWithCorners.cols()-corners.size());
        reducedMatrix.setFrom(seReducedMatrix);
        reducedMatrix.makeCompressed();

        bcInfo.addCornerValue(2, 0, 0);
        bcInfo.addCornerValue(3, 0, 0);
        bcInfo.addCornerValue(4, 0, 0);

        stiffnessInvWithoutCorners = gsPatchPreconditionersCreator<>::fastDiagonalizationOp(mb.basis(0), bcInfo, opt);

        gsMatrix<> result;
        gsMatrix<> stiffnessInvStiffness(reducedMatrix.rows(), reducedMatrix.cols());
        for (index_t i = 0; i < reducedMatrix.cols(); ++i) {
            stiffnessInvWithoutCorners->apply(reducedMatrix.toDense().col(i), result);
            stiffnessInvStiffness.col(i) = result;
        }

        CHECK((stiffnessInvStiffness-gsMatrix<>::Identity(stiffnessInvStiffness.rows(),stiffnessInvStiffness.cols())).norm() < 1.e-11);
        CHECK(stiffnessInvWithoutCorners->rows() + (index_t)corners.size() == stiffnessWithCorners.rows() && stiffnessInvWithoutCorners->cols() + (index_t)corners.size() == stiffnessWithCorners.cols());
    }

    TEST(StiffnessOp3)
    {
        gsLinearOperator<>::uPtr matrixInvWithoutCorners;
        gsSparseMatrix<> stiffnessWithCorners, massWithCorners;
        std::vector<index_t> corners;
        real_t alpha = 1.0;

        index_t numRefine = 1;
        index_t degree = 2;

        gsFunctionExpr<> bc("0.0",2);

        gsMultiPatch<> square(*gsNurbsCreator<>::BSplineSquare());
        square.computeTopology();

        gsMultiBasis<> mb(square);

        gsSparseEntries<> seReducedMatrix;
        gsSparseMatrix<> reducedMatrix;

        for (index_t i = 0; i < numRefine; ++i)
            mb.uniformRefine();

        //mb[0].component(0).uniformRefine();

        mb[0].component(0).setDegreePreservingMultiplicity(degree);
        mb[0].component(1).setDegreePreservingMultiplicity(degree+1);

        gsBoundaryConditions<> bcInfo;
        bcInfo.addCondition(0, 1, condition_type::neumann, bc);
        bcInfo.addCondition(0, 2, condition_type::dirichlet, bc);
        bcInfo.addCondition(0, 3, condition_type::neumann, bc);
        bcInfo.addCondition(0, 4, condition_type::dirichlet, bc);

        gsOptionList opt = gsAssembler<>::defaultOptions();
        stiffnessWithCorners = gsPatchPreconditionersCreator<>::stiffnessMatrix(mb.basis(0), bcInfo, opt);
        massWithCorners = gsPatchPreconditionersCreator<>::massMatrix(mb.basis(0), bcInfo, opt);

        corners.resize(1);
        corners[0] = 0;

        index_t row = 0, col = 0;
        for (index_t i = 0; i < stiffnessWithCorners.rows(); i++)
        {
            if ( std::find(corners.begin(), corners.end(), i) == corners.end() )
            {
                for (index_t j = 0; j < stiffnessWithCorners.cols(); j++)
                {
                    if( std::find(corners.begin(), corners.end(), j) == corners.end() )
                    {
                        seReducedMatrix.add(row, col, stiffnessWithCorners(i,j) + alpha * massWithCorners(i,j));
                        col++;                    }
                }
                col=0;
                row++;
            }
        }

        // reduce the stiffness matrix
        reducedMatrix.resize(stiffnessWithCorners.rows()-corners.size(), stiffnessWithCorners.cols()-corners.size());
        reducedMatrix.setFrom(seReducedMatrix);
        reducedMatrix.makeCompressed();

        gsInfo << "reduced: \n"<<(stiffnessWithCorners + alpha * massWithCorners).toDense()<<"\n";

        bcInfo.addCornerValue(1, 0, 0);
        bcInfo.addCornerValue(2, 0, 0);
        bcInfo.addCornerValue(3, 0, 0);
        bcInfo.addCornerValue(4, 0, 0);

        matrixInvWithoutCorners = gsPatchPreconditionersCreator<>::fastDiagonalizationOp(mb.basis(0), bcInfo, opt, alpha);

        gsMatrix<> result;
        gsMatrix<> matrixInvMatrix(reducedMatrix.rows(), reducedMatrix.cols());
        for (index_t i = 0; i < reducedMatrix.cols(); ++i) {
            matrixInvWithoutCorners->apply(reducedMatrix.toDense().col(i), result);
            matrixInvMatrix.col(i) = result;
        }

        CHECK((matrixInvMatrix-gsMatrix<>::Identity(matrixInvMatrix.rows(),matrixInvMatrix.cols())).norm() < 1.e-11);
        CHECK(matrixInvWithoutCorners->rows() + (index_t)corners.size() == stiffnessWithCorners.rows() && matrixInvWithoutCorners->cols() + (index_t)corners.size() == stiffnessWithCorners.cols());
    }
}
