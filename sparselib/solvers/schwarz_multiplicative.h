#ifndef SCHWARZ_MULTIPLICATIVE_H
#define SCHWARZ_MULTIPLICATIVE_H

namespace sparse_lib
{

template <typename M, typename V1, typename V2, typename P, typename Iter>
bool MultiplicativeSchwarz_Solver(const M &A, V1& x, const V2& b, P& preConditioner, Iter &It, bool setPreconditioner = false)
{
	typedef DenseVector<typename M::value_type>		VectorType;
	typedef DenseVectorView<typename M::value_type>	VectorViewType;
	const std::vector<SchwarzLinearSystemBase<typename M::value_type>*> blocks = A.engine().get_blocks();
	const std::vector<int> block_ptr = A.engine().get_block_ptr();
	
	VectorType r(x.size()), rhs(x.size());
	It.setBNorm(nrm2(b));
	
	r=b-A*x;
	rhs=b;
	while(!It.finished(nrm2(r))) {
		for(int i=0; i<blocks.size(); ++i) {
			VectorViewType xv = x.range(block_ptr[i], block_ptr[i+1]-1);
			VectorViewType rhsv = rhs.range(block_ptr[i], block_ptr[i+1]-1);
			if(!blocks[i]->solve(xv,rhsv)){
				abort();
			}
			rhs = b - A.engine().get_background()*x;
		}
		r=b-A*x;
		++It;
	}
	
	return It.hasConverged();
}

} /* sparse_lib */ 

#endif /* end of include guard: SCHWARZ_MULTIPLICATIVE_H */
