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
                          T ref, index_t ext,
			  T lambda, T tol, T tst)
	    : gsHFitting<d, T>(params, points, basis, ref, std::vector<unsigned>(2, ext), lambda),
              m_tol(tol)
	{
	    // Reduce m_params and m_point_values to tst randomly.
	    srand(time(NULL));
            index_t testPerc = 100.0 * tst; // tst expressed in %
            gsInfo << "testPerc: " << testPerc << std::endl;

	    std::vector<index_t> chosen, notChosen;
	    for(index_t i=0; i<params.cols(); i++)
	    {
		if(rand() % 100 < testPerc)
		    notChosen.push_back(i);
		else
		    chosen.push_back(i);
	    }

	    size_t numChosen = chosen.size();

            // Warning: for whatever reason, the parameters are stored
            // as cols but the fitting points as rows!
	    m_param_values.resize(Eigen::NoChange, numChosen);
	    m_points.resize      (numChosen, Eigen::NoChange);

	    for(size_t i=0; i<chosen.size(); i++)
	    {
		m_param_values.col(i) = params.col(chosen[i]);
		m_points.row(i)       = points.col(chosen[i]);
	    }

	    // Save the remaining 10% for verification.
	    m_testParams = gsMatrix<T>(params.rows(), notChosen.size());
	    m_testPoints = gsMatrix<T>(notChosen.size(), points.rows());
	    for(size_t i=0; i<notChosen.size(); i++)
	    {
		m_testParams.col(i) = params.col(notChosen[i]);
		m_testPoints.row(i) = points.col(notChosen[i]);
	    }

            gsInfo << "Setting aside " << m_testParams.cols() << "/" << m_param_values.cols() + m_testParams.cols() << " points for testing." << std::endl;
	}

	void compute(index_t maxIter);

    protected:
	T perc()
	{
	    return (100.0 * this->numPointsBelow(m_tol)) / m_pointErrors.size();
	}

	T validate()
	{
	    gsMatrix<T> values;
	    m_result->eval_into(m_testParams, values);

	    index_t numBelow = 0;
	    for(index_t i=0; i<m_testParams.cols(); i++)
	    {
		const T err = (m_testPoints.row(i) - values.col(i).transpose()).norm();
		if(err < m_tol)
		    numBelow++;
	    }
            //gsInfo << "below: " << numBelow << "/" << values.cols() << std::endl;
	    return (100.0 * numBelow)/values.cols();
	}

    protected:
	const real_t m_tol;

	gsMatrix<T> m_testParams, m_testPoints;

	using gsFitting<T>::m_param_values;
	using gsFitting<T>::m_points;
	using gsFitting<T>::m_basis;
	using gsFitting<T>::m_result;
	using gsFitting<T>::m_pointErrors;

	using gsHFitting<d, T>::m_lambda;
    };

    template <short_t d, class T>
    void gsVarianceFitting<d, T>::compute(index_t maxIter)
    {
	if(m_pointErrors.size() == 0)
	{
	    gsHFitting<d, T>::compute(m_lambda);
	    gsHFitting<d, T>::computeErrors();
	}
	for(index_t iter=0; iter < maxIter; iter++)
	{
	    gsHFitting<d, T>::nextIteration(m_tol, -1); // The -1 means we use m_ref.
	    gsHFitting<d, T>::computeErrors();
	    gsInfo << "iter: " << iter << std::endl;
            gsInfo << "basis: " << *m_basis << std::endl;
	    gsInfo << "full: " << perc() << "%, test: " << validate() << "%" << std::endl;
	    gsWriteParaview(*m_result, "iter-" + util::to_string(iter), 14400, false, true);
	}

        //gsWriteParaview(*m_result, "result", 640000, false, true);
        gsMatrix<T> data = m_points.transpose();
        gsWriteParaviewPoints(data, "data");
        gsMatrix<T> test = m_testPoints.transpose();
        gsWriteParaviewPoints(test, "test");
    }

} // namespace gismo

using namespace gismo;

int main(int argc, char *argv[])
{
    
    std::string fn = "fitting/deepdrawingC.xml";
    index_t iter = 5;
    real_t tol = 1e-5;
    index_t ext = 2;
    real_t ref = 0.1;
    real_t lambda = 1e-9;
    real_t tst = 0.1;

    gsCmdLine cmd("Variance fitting.");
    cmd.addInt("i", "iter", "number of iterations", iter);
    cmd.addReal("s", "lambda", "smoothing coefficient", lambda);
    cmd.addReal("p", "refPercent", "percentage of points to refine in each iteration", ref);
    cmd.addInt("q", "extension", "extension size", ext);
    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tol);
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addReal("t", "tst", "ratio of points set aside for validation", tst);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);

    gsKnotVector<real_t> uKnots(0.0, 1.0, 1, 4);
    gsKnotVector<real_t> vKnots(0.0, 1.0, 0, 4);
    gsTensorBSplineBasis<2> basis(uKnots, vKnots);
    gsTHBSplineBasis<2> thb(basis);
    gsVarianceFitting<2, real_t> fitting(uv, xyz, thb, ref, ext, lambda, tol, tst);
    fitting.compute(iter);

    return 0;
}
