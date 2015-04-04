#ifndef SPARSE_LIB_GMRES_SOLVER_H
#define SPARSE_LIB_GMRES_SOLVER_H

namespace sparse_lib
{

template < typename V1, typename M, typename V2, typename V3 >
void Update(V1& x, int k, M& h, V2& s, const V3 & V)
{
	V2 y(s.size());
	y=s;
	// Backsolve:
	for (int i = k; i >= 0; --i) {
		y(i) /= h(i,i);
		for (int j = i - 1; j >= 0; --j)
			y(j) -= h(j,i) * y(i);
	}
	
	for (int j = 0; j <= k; ++j) {
		x +=  y(j) * V[j];
	}
}

template < typename V1, typename M, typename V2, typename V3, typename PC >
void Update(V1& x, int k, M& h, V2& s, const V3 & V, PC& P)
{
	V2 y(s.size());
	V1 Mvy(x.size());
	V1 vy(x.size());
	y=s;
	// Backsolve:
	for (int i = k; i >= 0; --i) {
		y(i) /= h(i,i);
		for (int j = i - 1; j >= 0; --j)
			y(j) -= h(j,i) * y(i);
	}
	
	for (int j = 0; j <= k; ++j) {
		vy = y(j) * V[j];
		Mvy = P.apply(vy, Mvy);
		x +=  Mvy;
	}
}


template<typename T>
inline void GeneratePlaneRotation(const T &dx, const T &dy, T &cs, T &sn)
{
	/*if (dy == 0.0) {
		cs = 1.0;
		sn = (T)0.0;
	} else if (dx == (T)0.0) {
		cs = 0.0;
		sn = (T)1.0;
	}
	else {   //dy is always larger or equal than dx 
		const double rel = std::pow(dx/dy,2.0);
		const double ynorm = 1.0 / std::sqrt(1.0 + rel);
		cs = std::sqrt(rel) * ynorm;
		sn = cs * dy / dx;
	}*/
	/*if (dy == 0.0) {
		cs = (T)1.0;
		sn = (T)0.0;
	} else if(dx == 0.0) {
		cs = (T)0.0;
		sn = dy/std::abs(dy);
	} else if (std::abs(dy) > std::abs(dx)) {
		T tmp = dx/dy;
		sn = (T)1.0 / std::sqrt( (T)1.0 + std::pow(tmp, 2));
		cs = tmp * sn;
	} else {
		T tmp = dy/dx;
		cs = (T)1.0/sqrt((T)1.0 + std::pow(tmp, 2));
		sn = tmp * cs;
	}*/
	if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (std::abs(dy) > std::abs(dx)) {
    T temp = dx / dy;
    sn = 1.0 / std::sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    T temp = dy / dx;
    cs = 1.0 / std::sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}

template<typename T>
inline void ApplyPlaneRotation(T &dx, T &dy, const T &cs, const T &sn)
{
	T temp  =  cs * dx + sn * dy;

	dy = -sn * dx + cs * dy ;
	dx = temp;
}

template <typename M, typename V1, typename V2, typename P, typename Iter>
bool GMRES_Solver_RPC(const M &A, V1& x, const V2& b, P& preConditioner, Iter &It_out, int restart, bool setPreconditioner = false, bool output = false)
{
	typedef DenseVector< typename M::value_type > VectorType;
	typedef Matrix<MatrixEngine<typename M::value_type, GENERAL, ASYMMETRIC, DENSE_FULL> > DenseMatrixType;
	if(setPreconditioner)
		preConditioner.set(A);
	const int dim = A.rows();
	if (restart == 0)
		restart = dim/8.0+It_out.getMaxIter()/4.0;
	if (restart > It_out.getMaxIter())
		restart = It_out.getMaxIter();
	//restart_ = std::min(restart_, dim*2);
	VectorType s(restart+1), sn(restart+1), w(dim);
	VectorType cs(restart+1);
	DenseMatrixType H(restart+1,restart);
	
	std::vector<givens_rotation<typename M::value_type> > rotations(restart+1);
	
	typename M::value_type bnorm = nrm2(b);
	
	It_out.setStep(0);
	It_out.setBNorm(bnorm);
	
	if(bnorm == 0.0) {
		x = 0.0;
		return true;
	}
	
	std::string originalId = It_out.getId();
	It_out.setId(originalId + " - GMRES Outer");
	if (It_out.getBNorm() < 1e-16)
		It_out.setBNorm(1e-16);
	
	VectorType r(b.size());
	r = b-A*x;
	
	typename M::value_type beta = nrm2(r);
	
	if(It_out.finished(beta)) { // x passed to the solver was already within tolerance
		return true;
	}
	
	std::vector<VectorType> V(restart+1);
	
	while(!It_out.finished(beta)) {
		V[0] = r * (1.0/beta);
		s = 0.0;
		s(0) = beta;
		Iter It_in(It_out);
		It_in.setMaxIter(restart-1);
		It_in.setStep(0);
		It_in.setId(originalId + " - GMRES Inner");
		if(It_out.getOutputLevel() == ITER_SUMMARY)
			It_in.silence();
		int i = 0;
		do {
			VectorType Mv(w.size());
			Mv = preConditioner.apply(V[i], Mv);
			if(output) std::cout << "Calling w = preConditioner.apply(V[i],Mv); " << std::endl;
			w = A*Mv;
			for(int k=0; k<=i; k++) {
				H(k,i) = w*V[k];
				w -= H(k,i) * V[k];
			}
			const typename M::value_type normw = nrm2(w);
			H(i+1,i) = normw;
			V[i+1] = w * (1.0 / H(i+1,i));
			
			for (int k = 0; k < i; k++)							// [REF]: MED 2008-12-17- pp 268 Step 11a
				rotations[k].scalar_apply(H(k,i), H(k+1,i));
			rotations[i] = givens_rotation<typename M::value_type>(H(i,i), H(i+1,i));
			rotations[i].scalar_apply(H(i,i), H(i+1,i));
			rotations[i].scalar_apply(s(i), s(i+1));
			
			It_in++; It_out++; i++;
			//if(output) std::cout << "S(" << i << ") = " << s(i) << std::endl;
		} while(!It_in.finished(std::fabs(s(i))));
		beta = std::abs(s(i));
		if(output) std::cout << "Calling Update with " << i-1 << std::endl;
		if(output) std::cout << "nrm2(w): " << nrm2(w) << std::endl;
		Update(x, i-1, H, s, V, preConditioner); 								// [REF]: MED 2008-12-17- pp 268 Step 12
		if(!It_in.hasConverged()) {
			r=b-A*x;
			beta = nrm2(r);  								// [REF]: MED 2008-12-17- pp 268 Step 1b for next iteration
		}
	}
}

template <typename M, typename V1, typename V2, typename P, typename Iter>
bool FGMRES_Solver(const M &A, V1& x, const V2& b, P& preConditioner, Iter &It_out, int restart, bool setPreconditioner = false, bool output = false)
{
	typedef DenseVector< typename M::value_type > VectorType;
	typedef Matrix<MatrixEngine<typename M::value_type, GENERAL, ASYMMETRIC, DENSE_FULL> > DenseMatrixType;
	if(setPreconditioner)
		preConditioner.set(A);
	const int dim = A.rows();
	if (restart == 0)
		restart = dim/8.0+It_out.getMaxIter()/4.0;
	if (restart > It_out.getMaxIter())
		restart = It_out.getMaxIter();
	//restart_ = std::min(restart_, dim*2);
	VectorType s(restart+1), sn(restart+1), w(dim);
	VectorType cs(restart+1);
	DenseMatrixType H(restart+1,restart);
	
	std::vector<givens_rotation<typename M::value_type> > rotations(restart+1);
	
	typename M::value_type bnorm = nrm2(b);
	
	It_out.setStep(0);
	It_out.setBNorm(bnorm);
	
	if(bnorm == 0.0) {
		x = 0.0;
		return true;
	}
	
	std::string originalId = It_out.getId();
	It_out.setId(originalId + " - FGMRES Outer");
	if (It_out.getBNorm() < 1e-16)
		It_out.setBNorm(1e-16);
	
	VectorType r(b.size());
	r = b-A*x;
	
	typename M::value_type beta = nrm2(r);
	
	if(It_out.finished(beta)) { // x passed to the solver was already within tolerance
		return true;
	}
	
	std::vector<VectorType> V(restart+1);
	std::vector<VectorType> Z(restart+1);
	
	while(!It_out.finished(beta)) {
		V[0] = r * (1.0/beta);
		s = 0.0;
		s(0) = beta;
		Iter It_in(It_out);
		It_in.setMaxIter(restart-1);
		It_in.setStep(0);
		It_in.setId(originalId + " - FGMRES Inner");
		if(It_out.getOutputLevel() == ITER_SUMMARY)
			It_in.silence();
		int i = 0;
		do {
			VectorType Mv(w.size());
			Mv = preConditioner.apply(V[i], Mv);
			Z[i] = Mv;
			if(output) std::cout << "Calling w = preConditioner.apply(V[i],Mv); " << std::endl;
			w = A*Mv;
			for(int k=0; k<=i; k++) {
				H(k,i) = w*V[k];
				w -= H(k,i) * V[k];
			}
			const typename M::value_type normw = nrm2(w);
			H(i+1,i) = normw;
			V[i+1] = w * (1.0 / H(i+1,i));
			
			for (int k = 0; k < i; k++)							// [REF]: MED 2008-12-17- pp 268 Step 11a
				rotations[k].scalar_apply(H(k,i), H(k+1,i));
			rotations[i] = givens_rotation<typename M::value_type>(H(i,i), H(i+1,i));
			rotations[i].scalar_apply(H(i,i), H(i+1,i));
			rotations[i].scalar_apply(s(i), s(i+1));
			
			It_in++; It_out++; i++;
			//if(output) std::cout << "S(" << i << ") = " << s(i) << std::endl;
		} while(!It_in.finished(std::fabs(s(i))));
		beta = std::abs(s(i));
		if(output) std::cout << "Calling Update with " << i-1 << std::endl;
		if(output) std::cout << "nrm2(w): " << nrm2(w) << std::endl;
		Update(x, i-1, H, s, Z); 								// [REF]: MED 2008-12-17- pp 268 Step 12
		if(!It_in.hasConverged()) {
			r=b-A*x;
			beta = nrm2(r);  								// [REF]: MED 2008-12-17- pp 268 Step 1b for next iteration
		}
	}
}


template <typename M, typename V1, typename V2, typename P, typename Iter>
bool GMRES_Solver(const M &A, V1& x, const V2& b, P& preConditioner, Iter &It_out, int restart, bool setPreconditioner = false, bool output = false)
{
	typedef DenseVector< typename M::value_type > VectorType;
	typedef Matrix<MatrixEngine<typename M::value_type, GENERAL, ASYMMETRIC, DENSE_FULL> > DenseMatrixType;
	if(setPreconditioner)
		preConditioner.set(A);
	const int dim = A.rows();
	if (restart == 0)
		restart = dim/8.0+It_out.getMaxIter()/4.0;
	if (restart > It_out.getMaxIter())
		restart = It_out.getMaxIter();
	//restart_ = std::min(restart_, dim*2);
	VectorType s(restart+1), sn(restart+1), w(dim);
	VectorType cs(restart+1);
	DenseMatrixType H(restart+1,restart);
	
	std::vector<givens_rotation<typename M::value_type> > rotations(restart+1);
	
	VectorType Mb(b.size());
	if(output) std::cout << "Calling Mb = preConditioner.apply(b,Mb); " << std::endl;
	Mb = preConditioner.apply(b, Mb);
	It_out.setStep(0);
	typename M::value_type mbnorm = nrm2(Mb);
	It_out.setBNorm(mbnorm);
	if(It_out.getBNorm() == 0.0) {
		x = 0.0;
		return true;
	}
	std::string originalId = It_out.getId();
	It_out.setId(originalId + " - GMRES Outer");
	if (It_out.getBNorm() < 1e-16)
		It_out.setBNorm(1e-16);
	
	VectorType resid(b.size());
	VectorType r(b.size());
	r = b-A*x;
	if(output) std::cout << "Calling resid = preConditioner.apply(r,resid); " << std::endl;
	resid = preConditioner.apply(r, resid);
	typename M::value_type beta = nrm2(resid);
	if(It_out.finished(beta)) {
		return true;
	}
	std::vector<VectorType> V(restart+1);
	while(!It_out.finished(beta)) {
		V[0] = resid * (1.0/beta);		
		s=0.0;
		s(0) = beta;
		Iter It_in(It_out);
		It_in.setMaxIter(restart-1);
		It_in.setStep(0);
		It_in.setId(originalId + " - GMRES Inner");
		if(It_out.getOutputLevel() == ITER_SUMMARY)
			It_in.silence();
		int i = 0;
		do {
			VectorType AV(w.size());
			AV=A*V[i];
			if(output) std::cout << "Calling w = preConditioner.apply(AV,w); " << std::endl;
			w = preConditioner.apply(AV, w);
			for(int k=0; k<=i; k++) {
				H(k,i) = w*V[k];
				w -= H(k,i) * V[k];
			}
			const typename M::value_type normw = nrm2(w);
			H(i+1,i) = normw;
			V[i+1] = w * (1.0 / H(i+1,i));
			
			for (int k = 0; k < i; k++)							// [REF]: MED 2008-12-17- pp 268 Step 11a
				rotations[k].scalar_apply(H(k,i), H(k+1,i));
			rotations[i] = givens_rotation<typename M::value_type>(H(i,i), H(i+1,i));
			rotations[i].scalar_apply(H(i,i), H(i+1,i));
			rotations[i].scalar_apply(s(i), s(i+1));
			
			It_in++; It_out++; i++;
			//if(output) std::cout << "S(" << i << ") = " << s(i) << std::endl;
		} while(!It_in.finished(std::fabs(s(i))));
		beta = std::abs(s(i));
		if(output) std::cout << "Calling Update with " << i-1 << std::endl;
		if(output) std::cout << "nrm2(w): " << nrm2(w) << std::endl;
		Update(x, i-1, H, s, V); 								// [REF]: MED 2008-12-17- pp 268 Step 12
		if(!It_in.hasConverged()) {
			r=b-A*x;
			resid = preConditioner.apply(r, resid);  	// [REF]: MED 2008-12-17- pp 268 Step 1a for next iteration
			beta = nrm2(resid);  								// [REF]: MED 2008-12-17- pp 268 Step 1b for next iteration
		}
	}

	return It_out.hasConverged();
}



} /* sparse_lib */ 


#endif /* end of include guard: GMRES_SOLVER_H */
