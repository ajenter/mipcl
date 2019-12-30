#include <cfloat>
#include <cmath>
#include <limits>

double DBL_Zero=0.0001;//sqrt(DBL_EPSILON);

template<typename T> inline T min(T x, T y)
{
  return (x <= y)? x: y;
}

template<typename T> inline T abs(T x)
{
  return (x < 0)? -x: x;
}

template<typename T> inline bool lessOrEqual(T x, T y);
template<typename T> inline bool isZero(T x);
template<typename T> inline bool isPositive(T x);
template<typename T> inline bool isNegative(T x);
template<typename T> inline bool isNonNegative(T x);

inline void setDblZero(double zero)
{DBL_Zero=zero;}

inline double getDblZero()
{return DBL_Zero;}

template<> inline bool lessOrEqual<int>(int x, int y)
{
	return (x <= y)? true: false;
}

template<> inline bool lessOrEqual<double>(double x, double y)
{
	return (x <= y + DBL_Zero)? true: false;
}

template<> inline bool isZero<int>(int x)
{
	return (x == 0)? true: false;
}

template<> inline bool isZero<double>(double x)
{
	return (x <= DBL_Zero && -x <= DBL_Zero)? true: false;
}

template<> inline bool isPositive<int>(int x)
{
	return (x > 0)? true: false;
}

template<> inline bool isPositive<double>(double x)
{
	return (x > DBL_Zero)? true: false;
}

template<> inline bool isNegative<int>(int x)
{
	return (x < 0)? true: false;
}

template<> inline bool isNegative<double>(double x)
{
	return (-x > DBL_Zero)? true: false;
}


template<> inline bool isNonNegative<int>(int x)
{
	return (x >= 0)? true: false;
}

template<> inline bool isNonNegative<double>(double x)
{
	return (-x <= DBL_Zero)? true: false;
}

