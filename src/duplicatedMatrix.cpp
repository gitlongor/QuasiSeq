#include <R.h>
#include <Rinternals.h>
#include <map>

#include "CharSEXP.h"

template <typename T>
class rcVec {		// a row vec or a col vec from a column-major order matrix
	public:
		T * x; 		// pointer to the first element
		int len;    // length of vector: ncol for row vec; nrow for col vec
		int eltShift;  // index shift between adjacent elements: nrow for row vec; 1 for col vec
		int vecShift;  // index shift between adjacent vectors: 1 for row vec; nrow for col vec
		int nVec;		// number of vectors: nrow for row vec; ncol for col vec
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
class vecMap {  // a map with key being rcVec type; values (0) are not used.
	public:
		typedef std::map<rcVec<T>, int> tMap;
		rcVec<T> aRC; 
		typename tMap::iterator it;
		tMap rcvMap;
		std::pair<typename tMap::iterator, bool> returnPair;
		
		void duplicatedMat		(const T* x, const int* nrow, const int* ncol, int* const out, bool const byRow=true, bool const fromLast=false);
};

template <typename T>
void vecMap<T>::duplicatedMat (const T* x, const int* nrow, const int* ncol, int* const out, bool const byRow, bool const fromLast)
{
/* put a logical vector of duplicated rows of numeric matrix x into out */
	if(byRow){
		aRC.eltShift = aRC.nVec = (int)(*nrow);
		aRC.vecShift = 1;
		aRC.len = (int)(*ncol);
	}else{
		aRC.eltShift = 1;
		aRC.vecShift = aRC.len = (int)(*nrow);
		aRC.nVec = (int)(*ncol);
	}
	// map insert: if not previously inserted, the .second of returned pair is true; otherwise false. the .first is an iterator for the (previously) inserted element, which is not used. 
	if (fromLast) {
		aRC.x=const_cast<T*>(x) + ( byRow ? (*nrow)-1 : ((*ncol)-1)*(*nrow) ); 
		for(int i=aRC.nVec-1; i>=0; aRC.x -= aRC.vecShift)
			out[i--] = (int) !(rcvMap.insert( std::pair<rcVec<T>, int>(aRC, 0) ).second);
	}else {
		aRC.x=const_cast<T*>(x);
		for(int i=0; i<aRC.nVec; aRC.x += aRC.vecShift) 
			out[i++] = (int) !(rcvMap.insert( std::pair<rcVec<T>, int>(aRC, 0) ).second);
	}
	rcvMap.clear();
}

// instantiation of global objects:
vecMap<int> 		intVecMap;
vecMap<double> 		doubleVecMap;
vecMap<CharSEXP>	charsexpVecMap; 


extern "C" {

SEXP dupNumMat(SEXP x, SEXP MARGIN, SEXP fromLast)
{/* returns a logical vector of duplicated rows of numeric matrix x */
	SEXP out;
	int* dim;
	dim=INTEGER(getAttrib(x, R_DimSymbol));
	out = PROTECT(allocVector(LGLSXP, dim[*INTEGER(MARGIN)-1]));
	
	switch (TYPEOF(x)) {
		case REALSXP:
			doubleVecMap.duplicatedMat	(REAL(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
			break;
		case INTSXP:  // factor type is also covered here
			// if(!inherits(x, "factor"))
				intVecMap.duplicatedMat	(INTEGER(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
			// else {;} 
			break;
		case LGLSXP:
			intVecMap.duplicatedMat	(LOGICAL(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
			break;
		case STRSXP: {
			CharSEXP* charSexpPtr = new CharSEXP [ dim[0]*dim[1] ];
			for(int i=dim[0]*dim[1]-1; i>=0; --i)
				charSexpPtr[i].sexp = STRING_ELT(x, i);
			
			charsexpVecMap.duplicatedMat	(charSexpPtr, dim, dim+1, LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
			
			delete[] charSexpPtr;
			break;
		}
		default:
			error("C function 'dumNumMat' only accepts REALSXP, LGLSXP, INTSXP and STRSXP");
	}
	
	UNPROTECT(1);
	return out;
}

}
