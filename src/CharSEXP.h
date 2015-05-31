#ifndef R_R_H
#include <R.h>
#endif

#ifndef R_INTERNALS_H_
#include <Rinternals.h>
#endif

#include <cstring>

class CharSEXP{
	public:
		SEXP sexp;
		inline bool valid() {return( TYPEOF(sexp) == CHARSXP );}
		
		CharSEXP(SEXP );
		CharSEXP();
		
		friend inline bool operator< (const CharSEXP& lhs, const CharSEXP& rhs) 
		{return( 
			strcmp(const_cast<const char*>(CHAR(lhs.sexp)), const_cast<const char*>(CHAR(rhs.sexp)) )<0
		); }
					
		friend inline bool operator> (const CharSEXP& lhs, const CharSEXP& rhs) 
		{return( 
			strcmp(const_cast<const char*>(CHAR(lhs.sexp)), const_cast<const char*>(CHAR(rhs.sexp)) )>0
		); }
					
		friend inline bool operator== (const CharSEXP& lhs, const CharSEXP& rhs) 
		{return( 
			strcmp(const_cast<const char*>(CHAR(lhs.sexp)), const_cast<const char*>(CHAR(rhs.sexp)) )==0
		); }
					
		friend inline bool operator!= (const CharSEXP& lhs, const CharSEXP& rhs) 
		{return( 
			strcmp(const_cast<const char*>(CHAR(lhs.sexp)), const_cast<const char*>(CHAR(rhs.sexp)) )!=0
		); }
};

CharSEXP::CharSEXP(SEXP x)
{
	if (TYPEOF(x) == CHARSXP) sexp = x;
	else error("CharSEXP should be initialized with a CHARSXP type object");
}
		
CharSEXP::CharSEXP()
{
	sexp = mkChar("");
}
	

