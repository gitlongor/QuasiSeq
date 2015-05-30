#include <R.h>
#include <Rinternals.h>
#include <map>

template <typename T>
class rowVec {
	public:
		T * x; 		// pointer to the first element
		int len;    // length of vector
		int shift;  // index shift between adjacent elements
		friend inline bool operator< (const rowVec& lhs, const rowVec& rhs){
			// elementwise comparison of two vectors from the end
			// assuming operator< and operator> defined for type T 
			T L, R;
			for(int i=lhs.len-1; i>=0; i--){
				L= *(lhs.x+lhs.shift*i) ;
				R= *(rhs.x+rhs.shift*i) ;
				if(L > R) return(false);
				else if(L < R) return(true);
			}
			return(false);
		}
		friend inline bool operator> (const rowVec& lhs, const rowVec& rhs){return rhs < lhs;}
		friend inline bool operator<=(const rowVec& lhs, const rowVec& rhs){return !(lhs > rhs);}
		friend inline bool operator>=(const rowVec& lhs, const rowVec& rhs){return !(lhs < rhs);}
};

template <typename T>
class vecMap {  // a map with key being rowVec type; values are not used.
	public:
		typedef std::map<rowVec<T>, int> tMap;
		rowVec<T> arow; 
		typename tMap::iterator it;
		tMap rowMap;
		std::pair<typename tMap::iterator, bool> returnPair;
		
		void dupMat		(const T* x, const int* nrow, const int*ncol, int* const out, bool const fromLast=false);
};

template <typename T>
void vecMap<T>::dupMat (const T* x, const int* nrow, const int*ncol, int* const out, bool const fromLast)
{
/* put a logical vector of duplicated rows of numeric matrix x into out */
	int i;	
	arow.shift = (int)(*nrow);
	arow.len = (int)(*ncol);
	if (fromLast) {
		arow.x=const_cast<T*>(x) + (*nrow)-1;
		for(i=*nrow-1; i>=0; --i, --(arow.x))
			out[i] = (int) !(rowMap.insert( std::pair<rowVec<T>, int>(arow, 0) ).second);
	}else {
		arow.x=const_cast<T*>(x);
		for(i=0; i<*nrow; ++i, ++(arow.x))
			out[i] = (int) !(rowMap.insert( std::pair<rowVec<T>, int>(arow, 0) ).second);
	}
	rowMap.clear();
}

// instantiation of global objects:
vecMap<int> 	intVecMap;
vecMap<double> 	doubleVecMap;

extern "C" {

SEXP dupNumMat(SEXP x, SEXP fromLast)
{/* returns a logical vector of duplicated rows of numeric matrix x */
	SEXP out;
	int* dim;
	dim=INTEGER(getAttrib(x, R_DimSymbol));
	out = PROTECT(allocVector(LGLSXP, *dim));
	
	if (TYPEOF(x) == REALSXP) {
		doubleVecMap.dupMat	(REAL(x), dim, dim+1,  LOGICAL(out), (bool)(*(LOGICAL(fromLast))) );
	}else if (isInteger(x)){
		intVecMap.dupMat	(INTEGER(x), dim, dim+1,  LOGICAL(out), (bool)(*(LOGICAL(fromLast))) );
	}else if  (TYPEOF(x) == LGLSXP) {
		intVecMap.dupMat	(LOGICAL(x), dim, dim+1,  LOGICAL(out), (bool)(*(LOGICAL(fromLast))) );
	}else{
		error("C function 'dumNumMat' only accepts REALSXP, LGLSXP and INTSXP");
	}
	
	UNPROTECT(1);
	return out;
}

}
