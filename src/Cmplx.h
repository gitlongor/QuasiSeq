#ifndef R_COMPLEX_H
#include <Comlex.h>
#endif

inline bool operator< (const Rcomplex& lhs, const Rcomplex& rhs)
{
	if (lhs.r == rhs.r) return (lhs.i < rhs.i);
	return (lhs.r < rhs.r);
}
inline bool operator> (const Rcomplex& lhs, const Rcomplex& rhs){return rhs < lhs;}
inline bool operator<=(const Rcomplex& lhs, const Rcomplex& rhs){return !(lhs > rhs);}
inline bool operator>=(const Rcomplex& lhs, const Rcomplex& rhs){return !(lhs < rhs);}
inline bool operator==(const Rcomplex& lhs, const Rcomplex& rhs)
{
	return (lhs.r == rhs.r && lhs.i == rhs.i);
}
inline bool operator!=(const Rcomplex& lhs, const Rcomplex& rhs){return !(lhs == rhs);}