#include <gismo.h>
#include <stdlib.h>

namespace gismo
{

    template<short_t d, class T>
    class gsVarianceFitting : gsHFitting<d, T>
    {
	
    public:

	gsVarianceFitting(const gsMatrix<T>& params,
			  const gsMatrix<T>& points,
			  gsTHBSplineBasis<d, T>& basis,
			  T lambda)
	    : gsHFitting<d, T>(params, points, basis, 0.1, std::vector<unsigned>(2, 2), lambda),
	    m_ref(0.5), m_tol(1e-5)
	{
	    // Reduce m_params and m_point_values to 90% randomly.
	    srand(time(NULL));

	    std::vector<index_t> chosen, notChosen;
	    for(index_t i=0; i<params.cols(); i++)
	    {
		if(rand() % 10 == 9)
		    notChosen.push_back(i);
		else
		    chosen.push_back(i);
	    }

	    size_t numChosen = chosen.size();

	    m_param_values.resize(Eigen::NoChange, numChosen);
	    m_points.resize      (Eigen::NoChange, numChosen);

	    for(size_t i=0; i<chosen.size(); i++)
	    {
		m_param_values.col(i) = params.col(chosen[i]);
		m_points.col(i)       = points.col(chosen[i]);
	    }

	    // Save the remaining 10% for verification.
	    m_testParams = gsMatrix<T>(params.rows(), notChosen.size());
	    m_testPoints = gsMatrix<T>(points.rows(), notChosen.size());
	    for(size_t i=0; i<notChosen.size(); i++)
	    {
		m_testParams.col(i) = params.col(notChosen[i]);
		m_testPoints.col(i) = points.col(notChosen[i]);
	    }

	    //gsInfo << m_param_values << std::endl << std::endl << m_testParams << std::endl;
	    gsInfo << "tol:    " << m_tol << std::endl;
	    gsInfo << "lambda: " << m_lambda << std::endl;
	    gsInfo << "chosen: " << numChosen << std::endl;
	    gsInfo << "test:   " << notChosen.size() << std::endl;
	    gsInfo << "check:  " << numChosen + notChosen.size() << " ?= " << params.cols() << std::endl;
	}

	void compute();

    protected:
	T perc()
	{
	    this->computeErrors();
	    index_t numBelow = this->numPointsBelow(m_tol);
	    index_t numError = m_pointErrors.size();
	    gsInfo << "below: " << numBelow << std::endl;
	    gsInfo << "error: " << numError << std::endl;
	    return (100.0 * this->numPointsBelow(m_tol)) / m_pointErrors.size();
	}

	T validate()
	{
	    gsMatrix<T> values;
	    m_result->eval_into(m_testParams, values);

	    index_t numBelow = 0;
	    for(index_t i=0; i<m_testParams.rows(); i++)
	    {
		const T err = (m_testPoints.row(i) - values.col(i).transpose()).norm();
		if(err > m_tol)
		    numBelow++;
	    }
	    return (100.0 * numBelow)/values.cols();
	}

    protected:
	const real_t m_ref, m_tol;

	gsMatrix<T> m_testParams, m_testPoints;

	using gsFitting<T>::m_param_values;
	using gsFitting<T>::m_points;
	using gsFitting<T>::m_basis;
	using gsFitting<T>::m_result;
	using gsFitting<T>::m_pointErrors;

	using gsHFitting<d, T>::m_lambda;
    };

    template <short_t d, class T>
    void gsVarianceFitting<d, T>::compute()
    {
	if(m_pointErrors.size() == 0)
	{
	    gsHFitting<d, T>::compute(m_lambda);
	    gsHFitting<d, T>::computeErrors();
	}
	for(index_t iter=0; iter < 10; iter++)
	{
	    gsHFitting<d, T>::nextIteration(m_tol, m_tol);
	    gsHFitting<d, T>::computeErrors();
	    gsInfo << "iter: " << iter << std::endl;
	    //gsInfo << "below: " << this->numPointsBelow(m_tol) << std::endl;
	    gsInfo << "full: " << perc() << "%, test: " << validate() << "%" << std::endl;
	    // When I remove the call to maxPointError, we are suddenly at 100%. Why?!
	    //gsInfo << "max: "  << this->maxPointError() << std::endl;
	    //gsInfo << *m_basis << std::endl;
	}
    }

} // namespace gismo

using namespace gismo;

int main()
{
    
    std::string fn = "fitting/deepdrawingC.xml";

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);

    gsKnotVector<real_t> knots(-1.0, 1.0, 0, 4);
    gsTensorBSplineBasis<2> basis(knots, knots);
    gsTHBSplineBasis<2> thb(basis);
    gsVarianceFitting<2, real_t> fitting(uv, xyz, thb, 1e-6);
    fitting.compute();

    //gsWriteParaview(*fitting.result(), "result", 10000, false, true);

    return 0;
}
