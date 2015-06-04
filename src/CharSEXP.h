#ifndef R_R_H
#include <R.h>
#endif

#ifndef R_INTERNALS_H_
#include <Rinternals.h>
#endif

#include <cstring>

/* 
	NOTE: R_NaString is a different SEXP than mkChar("NA"), but holding the same string "NA". 
		  We will treat R_NaString to be smaller than every usual string, including mkChar("NA"). 
*/

class CharSEXP{
	public:
		SEXP sexp;
		inline bool valid() {return( TYPEOF(sexp) == CHARSXP );}
		
		CharSEXP(SEXP );
		CharSEXP();
		
		friend inline bool operator< (const CharSEXP& lhs, const CharSEXP& rhs) 
		{
			if (lhs.sexp == R_NaString) return( rhs.sexp != R_NaString );
			if (rhs.sexp == R_NaString) return(false);
			return( 
				strcmp(const_cast<const char*>(CHAR(lhs.sexp)), const_cast<const char*>(CHAR(rhs.sexp)) )<0
			); 
		}
					
		friend inline bool operator> (const CharSEXP& lhs, const CharSEXP& rhs) 
		{
			if (rhs.sexp == R_NaString) return( lhs.sexp != R_NaString );
			if (lhs.sexp == R_NaString) return(false);
			return( 
				strcmp(const_cast<const char*>(CHAR(lhs.sexp)), const_cast<const char*>(CHAR(rhs.sexp)) )>0
			); 
		}
			
		 
		friend inline bool operator== (const CharSEXP& lhs, const CharSEXP& rhs) 
		{	
			return (lhs.sexp == rhs.sexp);  // R CHARSXP objects are cached (only one copy per string)
			/* return( 
				strcmp(const_cast<const char*>(CHAR(lhs.sexp)), const_cast<const char*>(CHAR(rhs.sexp)) )==0
			); */
		}
					
		friend inline bool operator!= (const CharSEXP& lhs, const CharSEXP& rhs) 
		{ 
			return( lhs.sexp != rhs.sexp);
			/* return( 
				strcmp(const_cast<const char*>(CHAR(lhs.sexp)), const_cast<const char*>(CHAR(rhs.sexp)) )!=0
			); */
		}
};

CharSEXP::CharSEXP(SEXP x)
{
	if (TYPEOF(x) == CHARSXP) sexp = x;
	else error("CharSEXP should be initialized with a CHARSXP type object");
}
		
CharSEXP::CharSEXP()
{
	sexp = R_NaString; 
}
	

