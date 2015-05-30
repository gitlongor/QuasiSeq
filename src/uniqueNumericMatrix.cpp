#include <R.h>
#include <Rinternals.h>
#include <map>

typedef struct {
	double * x; 
	int len;
	int shift;
} tDoubleRowVec;

bool fncomp (tDoubleRowVec lhs, tDoubleRowVec rhs) 
{	double tmp;
	for(int i=lhs.len-1; i>=0; i--){
		// Rprintf("lhs[%d]=%.0f\trhs[%d]=%.0f\n", i, *(lhs.x+lhs.shift*i), i, *(rhs.x+rhs.shift*i));
		tmp =(*(lhs.x+lhs.shift*i) - *(rhs.x+rhs.shift*i)) ;
		if (tmp>0.0){
			// Rprintf("\tReturn false\n");
			return(false);
		}else if (tmp<0.0) {
			// Rprintf("\tReturn true\n");
			return(true);
		}
	}
	// Rprintf("\tReturn false\n");
	return(false);
}
bool(*fn_pt)(tDoubleRowVec,tDoubleRowVec) = fncomp;

typedef std::map<tDoubleRowVec,int,bool(*)(tDoubleRowVec,tDoubleRowVec)> tdmap;

tDoubleRowVec arow;
tdmap::iterator it;
tdmap rowMap(fn_pt);
// std::pair<tdmap::iterator, bool> returnPair;

extern "C" {

void duplicatedRowsNumericMatrix(const double* x, const int* nrow, const int* ncol, int* const out)
{/* put a logical vector of duplicated rows of numeric matrix x into out */
	int i;	
	arow.shift = (int)(*nrow);
	arow.len = (int)(*ncol);
	arow.x=(double *)x;
	for(i=0; i<*nrow; ++i, ++(arow.x))
		out[i] = (int) !(rowMap.insert( std::pair<tDoubleRowVec, int>(arow, 0) ).second);
	rowMap.clear();
}

void duplicatedRowsNumericMatrixFromLast(const double* x, const int* nrow, const int* ncol, int* const out)
{/* put a logical vector of duplicated rows of numeric matrix x into out */
 /* out is assumed to be initialized to FALSE */
	int i;	
	arow.shift = (int)(*nrow);
	arow.len = (int)(*ncol);
	arow.x=(double *)x + (*nrow)-1;
	for(i=*nrow-1; i>=0; --i, --(arow.x))
		out[i] = (int) !(rowMap.insert( std::pair<tDoubleRowVec, int>(arow, 0) ).second);
	rowMap.clear();
}

SEXP dupRowNumMat(SEXP x, SEXP fromLast)
{/* returns a logical vector of duplicated rows of numeric matrix x */
	SEXP out;
	int* dim;
	dim=INTEGER(getAttrib(x, R_DimSymbol));
	out = PROTECT(allocVector(LGLSXP, *dim));
	if(*(LOGICAL(fromLast))) duplicatedRowsNumericMatrixFromLast(REAL(x), dim, dim+1,  LOGICAL(out));
	else					 duplicatedRowsNumericMatrix(REAL(x), dim, dim+1,  LOGICAL(out));
	UNPROTECT(1);
	return out;
}




}
