#ifndef SPARSE_LIB_CG_SOLVER_H
#define SPARSE_LIB_CG_SOLVER_H

namespace sparse_lib
{

template <typename M, typename V1, typename V2, typename P, typename Iter>
bool CG_Solver(const M &A, V1& x, const V2& b, P& preConditioner, Iter &It, bool setPreconditioner = false)
{
	typedef DenseVector< typename M::value_type > VectorType;
	std::string oldId = It.getId();
	It.setId(oldId + " - CG");
	if(setPreconditioner)
		preConditioner.set(A);
	VectorType r(x.size()), z(x.size()), p(x.size()), q(x.size());
	It.setBNorm(nrm2(b));
	It.setStep(0);
	typename M::value_type alpha = 0, beta = 0, rho = 0, rho_1 = 0;
	r = b - A*x;
	//std::cout << "nrm2(r): " << nrm2(r) << std::endl;
	//std::cout << "nrm2(b): " << nrm2(b) << std::endl;
	while(!It.finished(nrm2(r))) {
		z = preConditioner.apply(r, z);
		//std::cout << "nrm2(z): " << nrm2(z) << std::endl;
		rho = r*z;
		if(It.getStep()==0)
			p = z;
		else {
			beta = rho / rho_1;
			p*=beta;
			p+=z;
		}
		q = A*p;
		alpha = rho / (p*q);
		x+=(alpha*p);
		r-=(alpha*q);
		rho_1 = rho;
		++It;
	}
	It.setId(oldId);
	return It.hasConverged();
}

} /* sparse_lib */ 


#endif /* end of include guard: CG_SOLVER_H */
