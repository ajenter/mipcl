///////////////////////////////////////////////////////////////
/**
 * \file Sort.h Function for sorting arrays of different types.
 *
 * |  __Author__  | N.N. Pisaruk                              |
 * |-------------:|:------------------------------------------|
 * |  __e-mail__  | nicolaipisaruk@gmail.com                  |
 * | __home page__| wwww.mipcl-cpp.appspot.com                |
 *
 *   \copyright __2015 Nicolai N. Pisaruk__
 */

 /*  This file is part of the Mixed Integer Class Library (MIPCL).
 *
 *  MIPCL is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  MIPCL is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with MIPCL; if not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __SORT__H
#define __SORT__H

#ifdef _WINDOWS
#ifndef MIP_API
#ifdef MIP_EXPORTS
#define MIP_API __declspec(dllexport)
#else
#define MIP_API __declspec(dllimport)
#endif // MIP_EXPORTS
#endif // MIP_API
typedef __int64 __LONG; ///< For Visual C++, `long` means 32 bit integer; so we redefine it.
#else
#ifndef MIP_API
/// WINDOWS specifics.
/**
 * All files within this DLL are compiled with the MIP_EXPORTS
 * symbol defined on the command line. This symbol should not be defined on any project
 * that uses this DLL. This way any other project whose source files include this file see
 * MIP_API functions as being imported from a DLL, whereas this DLL sees symbols
 * defined with this macro as being exported.
 */
#define MIP_API
#endif // MIP_API
typedef long long __LONG; ///< For GNU C++, `long` means 64 bit integer.
#endif // _WINDOWS

/// Namespace `SORT` contains function for sorting arrays of different types.
/**
 * Given a set of values \f$a_1,a_2,\dots,a_m\f$, and a subset _S_ of _n_ indices from {1,2,...,_m_},
 * list indices in _S_, \f$i_1,i_2,\ldots,i_n\f$, so that
 *   - either \f$a_{i_1}\ge a_{i_2}\ge\dots\ge a_{i_n}\f$,
 *   - or \f$a_{i_1}\le a_{i_2}\le\dots\le a_{i_n}\f$.
 */
namespace SORT {	
	MIP_API void incSortInt(int n, int* ipInd, const int* ipVal); ///< lists indices `ipInd[i]` in non-decreasing order of values `ipVal[ipInd[i]]`.
	MIP_API void decSortInt(int n, int* ipInd, const int* ipVal); ///< lists indices `ipInd[i]` in non-increasing order of values `ipVal[ipInd[i]]`.
	MIP_API void incSortLong(int n, int* ipInd, const __LONG* lpVal); ///< lists indices `ipInd[i]` in non-decreasing order of values `lpVal[ipInd[i]]`.
	MIP_API void incSortDouble(int n, int* ipInd, const double* dpVal); ///< lists indices `ipInd[i]` in non-decreasing order of values `dpVal[ipInd[i]]`.
	MIP_API void decSortDouble(int n, int* ipInd, const double* dpVal); ///< lists indices `ipInd[i]` in non-increasing order of values `dpVal[ipInd[i]]`.
	MIP_API void incSortPairs(int n, int* ipInd, const int* ipVal); ///< lists indices `ipInd[i]` in non-increasing lexicographic order of pair-values `(ipVal[ipInd[i]<<1],ipVal[(ipInd[i]<<1)+1])`.
	MIP_API void incSortPairs(int n, double* dpVal); ///< lists indices `ipInd[i]` in non-increasing lexicographic order of pair-values `(dpVal[ipInd[i]<<1],dpVal[(ipInd[i]<<1)+1])`.
	MIP_API void minK(int k, int n, int* ipInd, const double* dpVal); ///< lists indices `ipInd[i]` so that each of the first `k` values `dpVal[ipInd[i]]` is not greater than any other value.
	MIP_API void maxK(int k, int n, int* ipInd, const double* dpVal); ///< lists indices `ipInd[i]` so that each of the first `k` values `dpVal[ipInd[i]]` is not greater than any other values.
} // end of namespace SORT

#endif // #ifndef __SORT__H
