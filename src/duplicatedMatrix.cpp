#include <R.h>
#include <Rinternals.h>
#include <map>

/*
bool operator==(const char* lhs, const char* rhs)
{
	return(strcmp(lhs, rhs) == 0);
}
bool operator<(const char* lhs, const char* rhs)
{
	return(strcmp(lhs, rhs) < 0);
}
bool operator>(const char* lhs, const char* rhs)
{
	return(strcmp(lhs, rhs) > 0);
}
*/

template <typename T>
class rcVec {		// a row vec or a col vec
	public:
		T * x; 		// pointer to the first element
		int len;    // length of vector: ncol for row vec; nrow for col vec
		int eltShift;  // index shift between adjacent elements: nrow for row vec; 1 for col vec
		int vecShift;  // index shift between adjacent vectors: 1 for row vec; nrow for col vec
		friend inline bool operator< (const rcVec& lhs, const rcVec& rhs){
			// elementwise comparison of two vectors from the end
			// assuming operator== and operator< defined for type T 
			T L, R;
			for(int i=lhs.len-1; i>=0; i--){
				if ( (L= *(lhs.x+lhs.eltShift*i)) == (R= *(rhs.x+rhs.eltShift*i)) ) continue;
				return(L < R);
			}
			return(false);
		}
		friend inline bool operator> (const rcVec& lhs, const rcVec& rhs){return rhs < lhs;}
		friend inline bool operator<=(const rcVec& lhs, const rcVec& rhs){return !(lhs > rhs);}
		friend inline bool operator>=(const rcVec& lhs, const rcVec& rhs){return !(lhs < rhs);}
};

template <typename T>
class vecMap {  // a map with key being rcVec type; values are not used.
	public:
		typedef std::map<rcVec<T>, int> tMap;
		rcVec<T> aRC; 
		typename tMap::iterator it;
		tMap rowMap;
		std::pair<typename tMap::iterator, bool> returnPair;
		
		void dupMat		(const T* x, const int* nrow, const int* ncol, int* const out, bool const byRow=true, bool const fromLast=false);
};

template <typename T>
void vecMap<T>::dupMat (const T* x, const int* nrow, const int* ncol, int* const out, bool const byRow, bool const fromLast)
{
/* put a logical vector of duplicated rows of numeric matrix x into out */
	int i;	
	if(byRow){
		aRC.eltShift = (int)(*nrow);
		aRC.vecShift = 1;
		aRC.len = (int)(*ncol);
	}else{
		aRC.eltShift = 1;
		aRC.vecShift = (int)(*nrow);
		aRC.len = (int)(*nrow);
	}
	if (fromLast) {
		aRC.x=const_cast<T*>(x) + ( byRow ? (*nrow)-1 : ((*ncol)-1)*(*nrow) ); 
		for(i=aRC.len-1; i>=0; aRC.x -= aRC.vecShift)
			out[i--] = (int) !(rowMap.insert( std::pair<rcVec<T>, int>(aRC, 0) ).second);
	}else {
		aRC.x=const_cast<T*>(x);
		for(i=0; i<aRC.len; aRC.x += aRC.vecShift)
			out[i++] = (int) !(rowMap.insert( std::pair<rcVec<T>, int>(aRC, 0) ).second);
	}
	rowMap.clear();
}

// instantiation of global objects:
vecMap<int> 	intVecMap;
vecMap<double> 	doubleVecMap;


extern "C" {

SEXP dupNumMat(SEXP x, SEXP MARGIN, SEXP fromLast)
{/* returns a logical vector of duplicated rows of numeric matrix x */
	SEXP out;
	int* dim;
	dim=INTEGER(getAttrib(x, R_DimSymbol));
	out = PROTECT(allocVector(LGLSXP, dim[*INTEGER(MARGIN)-1]));
	
	if (TYPEOF(x) == REALSXP) {
		doubleVecMap.dupMat	(REAL(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
	}else if (isInteger(x)){
		intVecMap.dupMat	(INTEGER(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
	}else if  (TYPEOF(x) == LGLSXP) {
		intVecMap.dupMat	(LOGICAL(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
	}else if (TYPEOF(x) == CHARSXP) {
	}else{
		error("C function 'dumNumMat' only accepts REALSXP, LGLSXP and INTSXP");
	}
	
	UNPROTECT(1);
	return out;
}

}
