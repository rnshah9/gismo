#pragma once

#include <gismo.h>
#include <Eigen/Dense>
#include <gsModeling/gsFitting.h>
#include <gsModeling/gsFitting.hpp>
#include <gsSolver/gsSolverUtils.h>

namespace gismo {

template <class T>
class gsFastFitting : public gsFitting<T>
{
public:
    gsFastFitting(const gsMatrix<T> & param_values,
                  const gsMatrix<T> & points,
                  gsBasis<T> & basis,
                  const gsMatrix<T> & ugrid,
                  const gsMatrix<T> & vgrid)
        : gsFitting<T>(param_values, points, basis), m_ugrid(ugrid), m_vgrid(vgrid), m_weights(ugrid.cols(),vgrid.cols())
    {
        m_weights.setZero();
    }
    index_t FindGridPoint1Dfast(const gsMatrix<T> &grid, const T& curr_param, index_t beg, index_t end);
    void FindGridPoint1D(const gsMatrix<T> &grid, const T& curr_param, index_t& k2) const;
    void FindGridPoint(const gsMatrix<T>& curr_param, index_t &k1, index_t &k2);
    void FindGridPointfast(const gsMatrix<T>& curr_param, index_t& k1, index_t& k2);
    gsMatrix<T> ProjectParamToGrid(const gsMatrix<T> & uv);
    void GridProjection();
    gsMatrix<index_t> GridProjectionNEW();
    void assembleSystem(gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B);
    void assembleSystemNEW(gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B, gsMatrix<index_t>& gridindex);
    void compute(const bool condcheck);
    void computeNEW(const bool condcheck);
    void computeAllProjectedErrors(const gsMatrix<T>& uv, const gsMatrix<T>& xyz);
    void computeProjectedAverageErrors();
    void plotErrors(const std::string & fname) const;
    void computeHaussdorfErrors(const gsMatrix<T>& uv, const gsMatrix<T>& xyz, const bool grid);
    gsMatrix<T> ParameterCorrection(const gsMatrix<T>& uv, const gsMatrix<T>& xyz);
    T getL2error() const;
    void computeErrors(const gsMatrix<T>& uv, const gsMatrix<T>& xyz);
    using gsFitting<T>::computeErrors;
protected:
    void changeParam(const gsMatrix<T> &uv, const gsMatrix<T> &xyz);
    gsMatrix<T> BuildLookupTable(gsMatrix<index_t> & uactives, gsMatrix<T>& uvalues);
    gsMatrix<T> BuildLookupTableNEW(gsMatrix<index_t> & uactives, gsMatrix<T>& uvalues, gsMatrix<index_t>& gridindex);
protected:
    /// u parameter for grid in u-direction for grid projection / fast fitting
    gsMatrix<T> m_ugrid;
    /// v parameter for grid in v-direction for grid projection / fast fitting
    gsMatrix<T> m_vgrid;
    /// weights, index how many data get projected to one grid point
    gsMatrix<index_t> m_weights;
};


template<class T>
void gsFastFitting<T>::changeParam(const gsMatrix<T> &uv, const gsMatrix<T> &xyz)
{
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  "Wrong input");
    this -> m_param_values = uv;
    this -> m_points = xyz.transpose();
}


template<class T>
index_t gsFastFitting<T>::FindGridPoint1Dfast(const gsMatrix<T> &grid, const T& curr_param, index_t beg, index_t end)
{
    //if (beg == end)
    if (end-beg == 1)
    {
        //return beg;
        if (curr_param < (grid(beg)+grid(end))/2.0 )
            return beg;
        else
            return end;
    }
    else
    {
        index_t mid = std::ceil(beg+(end-beg)/2.0);

        if (curr_param < (grid(mid-1)+grid(mid))/2.0 )
        {
            return FindGridPoint1Dfast(grid,curr_param,beg,mid);
        }
        else
            return FindGridPoint1Dfast(grid,curr_param,mid,end);

    }
}


template<class T>
void gsFastFitting<T>::FindGridPoint1D(const gsMatrix<T> &grid, const T& curr_param, index_t& k2) const
{
    // Could be improved with a faster search method
    if ((grid(0)+grid(1))/2 >= curr_param)
    {
        k2 = 0;
        return;
    }

    for (index_t i= 1; i< grid.cols()-1 ; i++)
    {
        if ((grid(i-1)+grid(i))/2< curr_param && (grid(i)+grid(i+1))/2 >= curr_param)
        {
            k2 = i;
            return;
        }
    }

    k2 = grid.cols()-1;
}

template<class T>
void gsFastFitting<T>::FindGridPoint(const gsMatrix<T>& curr_param, index_t& k1, index_t& k2)
{
    // Search function, k is the global index of the nearest neighbour of curr_param, with k=(k1,k2)
    FindGridPoint1D(m_ugrid,curr_param(0),k1);
    FindGridPoint1D(m_vgrid,curr_param(1),k2);

//    index_t k1fast = FindGridPoint1Dfast(m_ugrid,curr_param(0),0,(m_ugrid.cols()));
//    if (k1 != k1fast || k2 != FindGridPoint1Dfast(m_vgrid,curr_param(1),0,m_vgrid.cols()))
//        gsInfo << "Wrong FindGridPoint1Dfast" << std::endl;
}

template<class T>
void gsFastFitting<T>::FindGridPointfast(const gsMatrix<T>& curr_param, index_t& k1, index_t& k2)
{
    // Search function, k is the global index of the nearest neighbour of curr_param, with k=(k1,k2)
    k1=FindGridPoint1Dfast(m_ugrid,curr_param(0),0,m_ugrid.cols()-1);
    k2=FindGridPoint1Dfast(m_vgrid,curr_param(1),0,m_vgrid.cols()-1);
}


template<class T>
gsMatrix<T> gsFastFitting<T>::ProjectParamToGrid(const gsMatrix<T> & uv)
{
    gsMatrix<T> uv_temp (2,uv.cols());
    index_t k1,k2;
    for (index_t k=0; k<uv.cols(); k++)
    {
        FindGridPointfast(uv.col(k),k1,k2);
        uv_temp(0,k) = m_ugrid(k1);
        uv_temp(1,k) = m_vgrid(k2);
    }

    return uv_temp;
}


template<class T>
void gsFastFitting<T>::GridProjection()
{
    gsInfo << "Begin GridProjection" << std::endl;
    // Initializing some variables
    index_t n1,n2;
    n1 = m_weights.rows();
    n2 = m_weights.cols();
    m_weights.setZero(n1,n2);

    const int dimension = this->m_points.cols();
    const int numpoints = this->m_points.rows();

    /*gsSparseMatrix<T> A_mat(num_basis, num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here, to do: improve
    int nonZerosPerCol = 1;
    for (int i = 0; i < this->m_basis->dim(); ++i)
        nonZerosPerCol *= ( 2 * this->m_basis->degree(i) + 1 ) * 4;
    A_mat.reservePerColumn( nonZerosPerCol );*/

    // Initialize new m_params, m_points
    /*gsMatrix<T> params_temp, points_temp;
    params_temp.setZero(2,n1*n2);
    points_temp.setZero(n1*n2,dimension);*/
    gsMatrix<T> curr_param,curr_point;

    //gsMatrix<T> params_temp(2,n1*n2);
    gsSparseMatrix<T> points_temp(n1*n2,dimension);
    points_temp.reservePerColumn(numpoints);

    for (index_t i=0; i<numpoints; i++)
    {
        // we have to check which points belong to which grid point. New variables:
        // m_params: ugrid x vgrid as tensor product structure -- Done in the next step
        // m_points: init with 0, add p_i to corresponding col i, where gridpoint (u_i1,v_i2) (of ugrid, vgrid) are the nearest neighbour of m_params
        // weights: increase weights(i) by one for each m_points added to i

        curr_param = this->m_param_values.col(i);
        curr_point = this->m_points.row(i);
        index_t k,k1=0,k2=0;
        //FindGridPoint(curr_param,k1,k2);
        //k = k2 * n1 + k1;

        //index_t k1fast=0,k2fast=0;
        FindGridPointfast(curr_param,k1,k2);
        k = k2 * n1 + k1;

        //if (k1 != k1fast)
        //    gsInfo << "Wrong fast gridprojection" << std::endl;

        points_temp(k,0) += curr_point(0);
        points_temp(k,1) += curr_point(1);
        points_temp(k,2) += curr_point(2);

        m_weights(k1,k2) += 1;
    }

    // Build m_param as mesh of parameters
    /*for (index_t i2= 0; i2<n2; i2++)
    {
        for (index_t i1= 0; i1<n1; i1++)
        {
            index_t iglobal = i2 * n1 + i1;
            params_temp(0,iglobal)=m_ugrid(i1);
            params_temp(1,iglobal)=m_vgrid(i2);
        }
    }*/

    // Write new parameter and point data into member variables
    //this -> m_param_values = params_temp;
    this -> m_points = points_temp;
}

template<class T>
gsMatrix<index_t> gsFastFitting<T>::GridProjectionNEW()
{
    index_t n1;
    n1 = m_weights.rows();
    const int numpoints = this->m_points.rows();
    gsMatrix<T> curr_param;
    gsMatrix<index_t> gridindex(numpoints,3);

    for (index_t i=0; i<numpoints; i++)
    {
        curr_param = this->m_param_values.col(i);
        index_t k1=0,k2=0;
        FindGridPointfast(curr_param,k1,k2);

        gridindex(i,0)= k2 * n1 + k1;
        gridindex(i,1)= k1;
        gridindex(i,2)= k2;
    }
    return gridindex;    // This now only contains the global index and u-, v-indices of the according grid point with respect to m_param_values
}

template<class T>
gsMatrix<T> gsFastFitting<T>::BuildLookupTable(gsMatrix<index_t> & uactives, gsMatrix<T>& uvalues)
{
    gsInfo<< "Begin Build Lookup Table" << std::endl;
    index_t N1 = this->m_basis->component(0).size();    // number of basis functions in u direction
    int p = this->m_basis->degree(0);   // degree for u
    int n1,n2;
    n1 = m_vgrid.cols();  // gridpoints in u-direction
    n2 = m_ugrid.cols();  // gridpoints in v-direction
    gsMatrix<T> Table;
    Table.setZero(N1*N1,n2);
    gsInfo<< "SetZero Build Lookup Table" << std::endl;

    // Build look-up table, i.e. elements h_{k2,i2,j2} in Table
    index_t flops = 0;
    for (index_t k1=0; k1<n1; k1++)
    {
        gsVector<T> bcol = uvalues.col(k1);
        gsVector<index_t> icol = uactives.col(k1);
        for (index_t i=0; i<p+1; i++)
        {
            T bi = bcol(i);//uvalues(i,k1);
            index_t i1 = icol(i);//uactives(i,k1);
            for (index_t j=0; j<p+1; j++)
            {
                T bj = bcol(j);//uvalues(j,k1);
                index_t j1 = icol(j);//uactives(j,k1);
                index_t ind_temp = j1*N1+i1;   // lexicographical running index
                for (index_t k2=0; k2<n2; k2++)
                {
                    index_t weight = m_weights(k1, k2);
                    if(weight!=0)
                    {
                        // Judge i < j and use symmetry of bi*bj -> more efficiently?
                        Table(ind_temp,k2) += bi*bj*weight;
                        flops += 3;
                    }
                }
            }
        }
    }
    // Output number flops:
    gsInfo << "Counted flops fast fitting--PART 1: " << flops << std::endl;

    // How can this be improved? -> parallelization?
    return Table;
}

template<class T>
gsMatrix<T> gsFastFitting<T>::BuildLookupTableNEW(gsMatrix<index_t> & uactives, gsMatrix<T>& uvalues, gsMatrix<index_t>& gridindex)
{
    index_t N1 = this->m_basis->component(0).size();    // number of basis functions in u direction
    int p = this->m_basis->degree(0);                   // degree for u
    int n2 = m_ugrid.cols();                            // gridpoints in v-direction
    //gsSparseMatrix<T> Table(N1*N1,n2);
    //Table.reservePerColumn((p+1)*(p+1));
    gsMatrix<T> Table(N1*N1,n2);
    index_t flops = 0;

    for (index_t k=0; k<gridindex.rows(); k++)
    {
        index_t k1 = gridindex(k,1);
        index_t k2 = gridindex(k,2);
        for (index_t i=0; i<p+1; i++)
        {
            T bi = uvalues(i,k1);
            index_t i1 = uactives(i,k1);
            for (index_t j=0; j<p+1; j++)
            {
                T bj = uvalues(j,k1);
                index_t j1 = uactives(j,k1);
                index_t ind_temp = j1*N1+i1;   // lexicographical running index
                Table(ind_temp,k2) += bi*bj;
                flops += 3;
            }
        }
    }

    // Output number flops:
    gsInfo << "Counted flops fast fitting--PART 1: " << flops << std::endl;
    return Table;
}

template<class T>
void gsFastFitting<T>::assembleSystem( gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B)
{
    //Initialize look-up Table
    const int num_basis=this->m_basis->size();
    //const int dimension=this->m_points.cols();

    // Get 1D basis in v direction:
    gsBasis<T>* vbasis_tmp;
    vbasis_tmp = &(this->m_basis->component(1));
    gsBSplineBasis<T> vbasis = *(static_cast<gsBSplineBasis<T>*>(vbasis_tmp));

    index_t N1;
    index_t N2;
    N2 = vbasis.size();    // number basis v direction
    N1 = num_basis/N2;    // number basis u direction

    int p;
    p = this->m_basis->degree(1);      // degree v direction
    int n1;
    n1 = m_ugrid.cols();  // gridpoints in v-direction
    int n2;
    n2 = m_vgrid.cols();  // gridpoints in u-direction
    //A_mat.setZero();

    // pre-evaluate vbasis in grid-points
    gsMatrix<index_t> vactives;
    gsMatrix<T> vvalues;
    vbasis.active_into(m_vgrid,vactives);
    vbasis.eval_into(m_vgrid,vvalues);

    // Get 1D basis in u direction:
    gsBasis<T>* ubasis_tmp;
    ubasis_tmp = &(this->m_basis->component(0));
    gsBSplineBasis<T> ubasis = *(static_cast<gsBSplineBasis<T>*>(ubasis_tmp));

    // pre-evaluate ubasis in grid-points
    gsMatrix<index_t> uactives;
    gsMatrix<T> uvalues;
    ubasis.active_into(m_ugrid,uactives);
    ubasis.eval_into(m_ugrid,uvalues);

    // Compute look-up table
    gsMatrix<T> Table;
    Table = BuildLookupTable(uactives,uvalues);

    // tensorbasis used in order to find lexicographical running index in inner slope
    //gsTensorBSplineBasis<2,T> *tensorbasis = (static_cast<gsTensorBSplineBasis<2,T>*>(this->m_basis));
    // Assemble matix A, right hand side b
    index_t flops = 0;
    for (index_t k2=0; k2<n2; k2++)    // grid points v
    {
        for (index_t i2=0; i2<p+1; i2++)     // active basis functions of v basis at grid point, approximation
        {
            T bi = vvalues(i2,k2);
            index_t ui = vactives(i2,k2);
            for (index_t j2=0; j2<p+1; j2++)      // active basis functions of v basis at grid point, test function
            {
                T bj = vvalues(j2,k2);
                index_t uj = vactives(j2,k2);
                for (index_t i1=0; i1<N1; i1++)       // all basis functions of u basis, approximation
                {
                    for (index_t j1=std::max(0,i1-p); j1<std::min(i1+p+1,N1); j1++)       // neighboring basis functions of u basis, test function
                    {
                        //index_t i_lexiglobal = tensorbasis->index(i1,ui);
                        index_t i_lexiglobal = ui * N1 + i1;
                        //index_t j_lexiglobal = tensorbasis->index(j1,uj);
                        index_t j_lexiglobal = uj * N1 + j1;
                        //if (i_lexiglobal!= i_lexiglobal2)
                        //    gsInfo << "Wrong lexicographical numbering" << std::endl;

                        index_t ind_temp = j1*N1+i1;
                        if (Table(ind_temp,k2) != 0)  // Check if Table is non-zero
                        {
                            A_mat(i_lexiglobal,j_lexiglobal) += bi*bj*Table(ind_temp,k2);   // Symmetry of bi*bj -> could be done more efficiently?
                            flops += 3;
                        }
                    }
                }
            }
        }
    }
    // Info for the flop count:
    gsInfo << "Counted flops fast fitting--PART 2: " << flops << std::endl;

    for(index_t k = 0; k < n1*n2; ++k)     // all grid points
    {
        // Check if weight(k) = 0: Find local k1,k2 from k
        index_t k1 = k % n1;
        index_t k2 = (k-k1)/n1;

        if (m_weights(k1,k2)!=0)
        {
            for (index_t i1=0; i1<p+1; i1++)
            {
                index_t iu = uactives(i1,k1);
                real_t bu = uvalues(i1,k1);
                for (index_t i2=0; i2<p+1; i2++)
                {
                    index_t iv = vactives(i2,k2);
                    real_t bv = vvalues(i2,k2);
                    index_t ilexi = iv*N1+iu;
                    m_B.row(ilexi) += bu*bv* this->m_points.row(k);
                }
            }
        }
    }

    //Old version: new version is more efficient
    /*gsMatrix<T> value, curr_point;
    gsMatrix<index_t> active;
    for(index_t k = 0; k < n1*n2; ++k)
    {
        // Check if weight(k) = 0: Find local k1,k2 from k
        index_t k1 = k % n1;
        index_t k2 = (k-k1)/n1;

        if (m_weights(k1,k2)!=0)
        {
            curr_point = this->m_param_values.col(k);

            //computing the values of the basis functions at the current point
            this->m_basis->eval_into(curr_point, value);

            // which functions have been computed i.e. which are active 500, 500. 1.
            this->m_basis->active_into(curr_point, active);

            const index_t numActive = active.rows();

            for (index_t i = 0; i != numActive; ++i)
            {
                const int ii = active.at(i);
                m_B.row(ii) += value.at(i) * this->m_points.row(k);
            }
        }
    }*/

    /*// For testing purpose to compare A_mat:
    gsFileData<> fd;
    gsMatrix<T> C=A_mat.toDense();
    fd << C;
    fd.dump("FastFittingMatrix");
    gsFileData<> fb;
    fb << m_B;
    fb.dump("FastFittingb");*/
}

template<class T>
void gsFastFitting<T>::assembleSystemNEW( gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B, gsMatrix<index_t>& gridindex)
{
    //Initializing
    const int num_basis=this->m_basis->size();
    const int numpoints = this->m_points.rows();
    // Get 1D basis in v direction:
    gsBasis<T>* vbasis_tmp;
    vbasis_tmp = &(this->m_basis->component(1));
    gsBSplineBasis<T> vbasis = *(static_cast<gsBSplineBasis<T>*>(vbasis_tmp));
    index_t N2 = vbasis.size();    // number basis v direction
    index_t N1 = num_basis/N2;    // number basis u direction
    int p = this->m_basis->degree(1);      // degree v direction
    //A_mat.setZero();

    // pre-evaluate vbasis in grid-points
    gsMatrix<index_t> vactives;
    gsMatrix<T> vvalues;
    vbasis.active_into(m_vgrid,vactives);
    vbasis.eval_into(m_vgrid,vvalues);

    // Get 1D basis in u direction:
    gsBasis<T>* ubasis_tmp;
    ubasis_tmp = &(this->m_basis->component(0));
    gsBSplineBasis<T> ubasis = *(static_cast<gsBSplineBasis<T>*>(ubasis_tmp));

    // pre-evaluate ubasis in grid-points
    gsMatrix<index_t> uactives;
    gsMatrix<T> uvalues;
    ubasis.active_into(m_ugrid,uactives);
    ubasis.eval_into(m_ugrid,uvalues);

    // Compute look-up table
    gsMatrix<T> Table;//(N1*N1,n2);
    //Table.reservePerColumn((p+1)*(p+1));
    Table = BuildLookupTableNEW(uactives,uvalues,gridindex);
    index_t flops = 0;

    std::set<index_t> k2set;
    for (index_t i=0; i< numpoints; i++)
        k2set.insert(gridindex(i,2));

    for ( std::set<index_t>::const_iterator k2=k2set.begin() ; k2!=k2set.end(); ++k2 )
    {
        for (index_t i2=0; i2<p+1; i2++)     // active basis functions of v basis at grid point, approximation
        {
            T bi = vvalues(i2,*k2);
            index_t ui = vactives(i2,*k2);
            for (index_t j2=0; j2<p+1; j2++)      // active basis functions of v basis at grid point, test function
            {
                T bj = vvalues(j2,*k2);
                index_t uj = vactives(j2,*k2);
                for (index_t i1=0; i1<N1; i1++)       // all basis functions of u basis, approximation
                {
                    for (index_t j1=std::max(0,i1-p); j1<std::min(i1+p+1,N1); j1++)       // neighboring basis functions of u basis, test function
                    {
                        index_t i_lexiglobal = ui * N1 + i1;
                        index_t j_lexiglobal = uj * N1 + j1;
                        index_t ind_temp = j1*N1+i1;
                        real_t table = Table(ind_temp,*k2);
                        if (table != 0)
                        {
                            A_mat(i_lexiglobal,j_lexiglobal) += bi*bj*table;   // could be done more efficiently?
                            flops += 3;
                        }
                    }
                }
            }
        }
    }
    // Info for the flop count:
    gsInfo << "Counted flops fast fitting--PART 2: " << flops << std::endl;

    for(index_t k = 0; k < numpoints; ++k)     // all grid points
    {
        // Check if weight(k) = 0: Find local k1,k2 from k
        index_t k1 = gridindex(k,1);
        index_t k2 = gridindex(k,2);

        for (index_t i1=0; i1<p+1; i1++)
        {
            index_t iu = uactives(i1,k1);
            real_t bu = uvalues(i1,k1);
            for (index_t i2=0; i2<p+1; i2++)
            {
                index_t iv = vactives(i2,k2);
                real_t bv = vvalues(i2,k2);
                index_t ilexi = iv*N1+iu;
                m_B.row(ilexi) += bu*bv* this->m_points.row(k);
            }
        }
    }
}


template<class T>
void gsFastFitting<T>::compute(const bool condcheck)
{
    // Wipe out previous result
    if ( this->m_result )
        delete this->m_result;

    // Initialization of variables num_basis, dimension
    const int num_basis=this->m_basis->size();
    const int dimension=this->m_points.cols();

    //left side matrix
    gsSparseMatrix<T> A_mat(num_basis, num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here, to do: improve
    int nonZerosPerCol = 1;
    for (int i = 0; i < this->m_basis->dim(); ++i)
        nonZerosPerCol *= ( 2 * this->m_basis->degree(i) + 1 ) * 4;
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis, dimension);
    m_B.setZero(); // ensure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    // Initialize weights and look-up table
    //gsMatrix<index_t> weights (m_ugrid.cols(),m_vgrid.cols());
    //gsMatrix<T> Table;
    GridProjection();
    //weights.setOnes();

    // stopping time for complete assembly (including right-hand side)
    gsStopwatch time;
    time.restart();
    assembleSystem( A_mat, m_B);
    time.stop();
    gsInfo<<"Assembly time                     : "<< time <<"\n";
    // To do: include regularization later after projection step (not needed at this point)
    //applySmoothing(lambda, A_mat, dreg);

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    //gsDebugVar( A_mat.nonZerosPerCol().maxCoeff() );
    //gsDebugVar( A_mat.nonZerosPerCol().minCoeff() );
    A_mat.makeCompressed();

    typename gsSparseSolver<T>:: BiCGSTABILUT solver( A_mat );
    gsMatrix<T> x;
    if (condcheck)
    {
        Eigen::SparseMatrix<double> I(A_mat.rows(),A_mat.cols());
        I.setIdentity();
        auto A_inv = solver.solve(I);
        gsInfo << "Condition number: " << A_mat.norm() * A_inv.norm() << std::endl;
        x = A_inv * m_B;
    }

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        std::cerr<<  "The preconditioner failed. Aborting.";// << std::endl;
        this->m_result = NULL;
        return;
    }

    // Solve for x
    x = solver.solve(m_B); //toDense()

    // Generate the B-spline curve / surface
    this->m_result = this->m_basis->makeGeometry( give(x) ).release();

    // Compute average errors
    //computeProjectedAverageErrors();
}

template<class T>
void gsFastFitting<T>::computeNEW(const bool condcheck)
{
    // Wipe out previous result
    if ( this->m_result )
        delete this->m_result;

    // Initialization of variables num_basis, dimension
    const int num_basis=this->m_basis->size();
    const int dimension=this->m_points.cols();
    const int numpoints=this->m_points.rows();

    //left side matrix
    gsSparseMatrix<T> A_mat(num_basis, num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here, to do: improve
    int nonZerosPerCol = 1;
    for (int i = 0; i < this->m_basis->dim(); ++i)
        nonZerosPerCol *= ( 2 * this->m_basis->degree(i) + 1 ) * 4;
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis, dimension);
    m_B.setZero(); // ensure that all entries are zero in the beginning

    gsMatrix<index_t> gridindex (numpoints,3);
    gridindex = GridProjectionNEW();

    // stopping time for complete assembly (including right-hand side)
    gsStopwatch time;
    time.restart();
    assembleSystemNEW( A_mat, m_B, gridindex);
    time.stop();
    gsInfo<<"Assembly time NEW                 : "<< time <<"\n";
    A_mat.makeCompressed();

    typename gsSparseSolver<T>:: BiCGSTABILUT solver( A_mat );
    gsMatrix<T> x;
    if (condcheck)
    {
        Eigen::SparseMatrix<double> I(A_mat.rows(),A_mat.cols());
        I.setIdentity();
        auto A_inv = solver.solve(I);
        gsInfo << "Condition number: " << A_mat.norm() * A_inv.norm() << std::endl;
        x = A_inv * m_B;
    }

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        std::cerr<<  "The preconditioner failed. Aborting.";// << std::endl;
        this->m_result = NULL;
        return;
    }

    // Solve for x
    x = solver.solve(m_B); //toDense()

    // Generate the B-spline curve / surface
    this->m_result = this->m_basis->makeGeometry( give(x) ).release();

    // Compute average errors
    //computeProjectedAverageErrors();
}


template<class T>
void gsFastFitting<T>::computeAllProjectedErrors(const gsMatrix<T>& uv, const gsMatrix<T>& xyz)
{
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  "Wrong input");

    // since m_param_values and m_points were adjusted we need to look at the original values uv and xyz:
    this->m_pointErrors.clear();

    // Leave xyz the same but project each uv to the grid:
    gsMatrix<T> uv_temp (2,uv.cols());
    uv_temp = ProjectParamToGrid(uv);

    computeErrors(uv_temp,xyz.transpose());

/*    index_t k1,k2;
    for (index_t k=0; k<uv.cols(); k++)
    {
        FindGridPoint(uv.col(k),k1,k2);
        uv_temp(0,k) = m_ugrid(k1);
        uv_temp(1,k) = m_vgrid(k2);
    }

    // We change the member variables and compute the error:
    //changeParam(uv_temp,xyz);
    computeErrors(uv_temp,xyz);

    // change back param and points?*/
}


template<class T>
void gsFastFitting<T>::computeProjectedAverageErrors()
{
    this->m_pointErrors.clear();

    gsMatrix<T> val_i;
    bool firsterrorfound = false;
    this->m_result->eval_into(this->m_param_values, val_i);

    index_t k1,k2;
    FindGridPoint(this->m_param_values.col(0),k1,k2);
    if (m_weights(k1,k2) != 0)
    {
        this->m_pointErrors.push_back( (this->m_points.row(0)/m_weights(k1,k2) - val_i.col(0).transpose()).norm() );
        this->m_max_error = this->m_min_error = this->m_pointErrors.back();
        firsterrorfound = true;
    }

    for (index_t i = 1; i < this->m_points.rows(); i++)
    {
        FindGridPoint(this->m_param_values.col(i),k1,k2);   // Find grid point for the current parameter
        if (m_weights(k1,k2) != 0)
        {
            const T err = (this->m_points.row(i)/m_weights(k1,k2) - val_i.col(i).transpose()).norm() ;    // this averages all point error contributions of the data projected to this grid point (k1,k2)
            this->m_pointErrors.push_back(err);

            if (firsterrorfound==false)
            {
                this->m_max_error = this->m_min_error = this->m_pointErrors.back();
                firsterrorfound = true;
            }

            if ( err > this->m_max_error ) this->m_max_error = err;
            if ( err < this->m_min_error ) this->m_min_error = err;
        }
    }
}


template<class T>
void gsFastFitting<T>::plotErrors(const std::string & fname) const
{
    // Plot Errors
    const std::vector<real_t>& eval_field = this->pointWiseErrors();
    gsMatrix<real_t> bigMatrix(4, eval_field.size());
    for(size_t i=0; i<eval_field.size(); i++)
    {
        bigMatrix(0,i) = this->m_param_values(0,i);
        bigMatrix(1,i) = this->m_param_values(1,i);
        bigMatrix(2,i) = 0;
        bigMatrix(3,i) = eval_field[i];
    }
    gsWriteParaviewPoints(bigMatrix, fname);
}


template<class T>
gsMatrix<T> gsFastFitting<T>::ParameterCorrection(const gsMatrix<T>& uv, const gsMatrix<T>& xyz)
{
    // We need to calculate the approximate parameter with minimal distance to the result.
    // Input needs to be in the same form as m_params, m_points
    GISMO_ASSERT( uv.cols() == xyz.rows() && uv.rows() == 2 && xyz.cols() == 3,
                  "Wrong input");

    gsMatrix<T> delta_uv(2,uv.cols());
    // Evaluation
    gsMatrix<T> derivs, evals;
    this->result()->deriv_into(uv,derivs);
    this->result()->eval_into(uv,evals);

    for (index_t i =0; i != uv.cols(); i++)
    {
        // Find 2x2 system for parameter correction
        gsMatrix<T> A_sys(2,2);
        gsMatrix<T> b_sys(2,1);
        A_sys(0,0) = derivs(0,i)*derivs(0,i)+derivs(2,i)*derivs(2,i)+derivs(4,i)*derivs(4,i);
        A_sys(1,0) = derivs(0,i)*derivs(1,i)+derivs(2,i)*derivs(3,i)+derivs(4,i)*derivs(5,i);
        A_sys(0,1) = A_sys(1,0);
        A_sys(1,1) = derivs(1,i)*derivs(1,i)+derivs(3,i)*derivs(3,i)+derivs(5,i)*derivs(5,i);

        b_sys(0,0) = (xyz(i,0)-evals(0,i))*derivs(0,i)+(xyz(i,1)-evals(1,i))*derivs(2,i)+(xyz(i,2)-evals(2,i))*derivs(4,i);
        b_sys(1,0) = (xyz(i,0)-evals(0,i))*derivs(1,i)+(xyz(i,1)-evals(1,i))*derivs(3,i)+(xyz(i,2)-evals(2,i))*derivs(5,i);

        // Solve system
        delta_uv.col(i) = A_sys.fullPivLu().solve(b_sys);
    }
    return delta_uv;
}


template<class T>
void gsFastFitting<T>::computeHaussdorfErrors(const gsMatrix<T>& uv, const gsMatrix<T>& xyz, const bool grid)
{
    // We need to calculate the error according to the approx. minimal distance to the result.
    // We use parameter correction function in order to calculate this:
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  "Wrong input");

    gsMatrix<T> uv_temp (2,uv.cols());
    if (grid)
        uv_temp = ProjectParamToGrid(uv);
    else
        uv_temp = uv;

    gsMatrix<T> delta_uv(2,uv.cols());
    delta_uv = ParameterCorrection(uv_temp,xyz.transpose());

    // Change Parameters:
    gsMatrix<T> new_param(2,uv.cols());
    real_t damp = 0.5;
    for (index_t i = 0; i< uv.cols(); i++)
    {
        new_param(0,i) = std::min(1.0, std::max(uv_temp(0,i)+damp*delta_uv(0,i),0.0));
        new_param(1,i) = std::min(1.0, std::max(uv_temp(1,i)+damp*delta_uv(1,i),0.0));
    }
    //new_param = std::min(1.0, std::max( uv_temp + damp*delta_uv, 0.0));

    computeErrors(new_param, xyz.transpose());
}


template<class T>
T gsFastFitting<T>::getL2error() const
{
    // TODO: gismo assert d = 2.
    // type : L^2 error if h is uniform & the same in both directions
    real_t result = 0;

    for (size_t i=0; i < this->m_pointErrors.size(); i++)
    {
         result += this->m_pointErrors[i]*this->m_pointErrors[i];
    }
    return sqrt(result / this->m_pointErrors.size());
}


template<class T>
void gsFastFitting<T>::computeErrors(const gsMatrix<T>& uv, const gsMatrix<T>& xyz)
{
    GISMO_ASSERT( uv.cols() == xyz.rows() && uv.rows() == 2 && xyz.cols() == 3,
                  "Wrong input");

    this->m_pointErrors.clear();

    gsMatrix<T> val_i;
    //m_result->eval_into(m_param_values.col(0), val_i);
    this->m_result->eval_into(uv, val_i);
    this->m_pointErrors.push_back( (xyz.row(0) - val_i.col(0).transpose()).norm() );
    this->m_max_error = this->m_min_error = this->m_pointErrors.back();

    for (index_t i = 1; i < xyz.rows(); i++)
    {
        //m_result->eval_into(m_param_values.col(i), val_i);

        const T err = (xyz.row(i) - val_i.col(i).transpose()).norm() ;

        this->m_pointErrors.push_back(err);

        if ( err > this->m_max_error ) this->m_max_error = err;
        if ( err < this->m_min_error ) this->m_min_error = err;
    }
}

} //namespace gismo
