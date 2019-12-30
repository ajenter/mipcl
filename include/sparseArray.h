///////////////////////////////////////////////////////////////
/**
 * \file sparseArray.h Interface for `SparseArray` structure
 *
 * |  __Author__  | N.N. Pisaruk                              |
 * |-------------:|:------------------------------------------|
 * |  __e-mail__  | nicolaipisaruk@gmail.com                  |
 * | __home page__| wwww.mipcl-cpp.appspot.com                |
 *
 *   \copyright __2019 Nicolai N. Pisaruk__
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
#ifndef __SPARSEARRAY_H_
#define __SPARSEARRAY_H_

#ifdef _WINDOWS
#ifndef MIP_API
/// WINDOWS specifics.
/**
 * All files within this DLL are compiled with the MIP_EXPORTS
 * symbol defined on the command line. This symbol should not be defined on any project
 * that uses this DLL. This way any other project whose source files include this file see
 * MIP_API functions as being imported from a DLL, whereas this DLL sees symbols
 * defined with this macro as being exported.
 */
#ifdef MIP_EXPORTS
#define MIP_API __declspec(dllexport)
#else
#define MIP_API __declspec(dllimport)
#endif
#endif
#ifndef __ALWAYS_INLINE
#define __ALWAYS_INLINE __forceinline
#endif
#else
#ifndef MIP_API
#define MIP_API  ///< WINDOWS specifics, not used in other systems.
#endif
#ifndef __ALWAYS_INLINE
/**
 * All functions preceded with this attribute are always inlined.
 */
#define __ALWAYS_INLINE __attribute__((always_inline)) inline
#endif
#endif

 ///  This structure has been designed for storing sparse arrays.
 /**
  * This structure implements a simple but very efficient structure for storing sparse arrays.
  * It has been designed just for using within `MIPCL`.
  * For efficiency reasons, `SparseArray` lacks a lot what is needed for safely using in external applications.
  *
  * It is up to the user to allocate memory for three internal plain arrays.
  */
struct SparseArray {
	int sz; ///< Number of non-zero entries.
	double *Val; ///< Non-zero array entries.
	int *ipInd; ///< Indices of non-zero entries.
	int *ipPos; ///< Position of non-zero entries: for `p=ipPos[i]`, if `p >= 0 && p < sz`, and if `ipInd[p]==i`, then the value of the `i`-th array element is `Val[p]`; otherwise, the value of the `i`-th array element is `0`.

//	/**
//	 *  The function returns the position of a given elemen.
//	 *
//	 * \return position of element `ind` if it is in the array, or `-1` otherwise.
//	 */
//	__ALWAYS_INLINE int getPos(int ind)
//	{
//		int pos{ipPos[ind]};
//		if (pos >=0 && pos < sz) {
//			if (ipInd[pos] != ind)
//				pos=-1;
//		}
//		else pos=-1;
//		return pos;
//	};


	/**
	 *  This function appends to the end of this array a new element.
	 *
	 * \param[in] ind, w `w` is value of element indexed by `ind`.
	 */
	__ALWAYS_INLINE void append(int ind, double w)
	{
		ipPos[ipInd[sz]=ind]=sz;
		Val[sz++]=w;
	};

	/**
	 *  The function sets a new value of the given element.
	 *
	 * \param[in] ind, w `w` is new value of element indexed by `ind`.
	 */
	__ALWAYS_INLINE void set(int ind, double w)
	{
		int pos{ipPos[ind]};
		bool flag{false};
		if (pos >=0 && pos < sz)
			if (ipInd[pos] == ind)
				flag=true;
		if (flag)
			Val[pos]=w;
		else {
			ipPos[ipInd[sz]=ind]=sz;
			Val[sz++]=w;
		}
	};

	/**
	 * \param[in] ind index of array element.
	 * \return the value of the element indexed by `ind`.
	 */
	__ALWAYS_INLINE double get(int ind)
	{
		double w{0.0};
		int pos{ipPos[ind]};
		if (pos >=0 && pos < sz)
			if (ipInd[pos] == ind)
				w=Val[pos];
		return w;
	};

	/**
	 *  The function adds a given value to a given element.
	 *
	 * \param[in] ind,w element indexed by `ind` is incremented by `w`.
	 */
	__ALWAYS_INLINE void add(int ind, double w)
	{
		int pos{ipPos[ind]};
		bool flag{false};
		if (pos >=0 && pos < sz)
			if (ipInd[pos] == ind)
				flag=true;
		if (flag)
			Val[pos]+=w;
		else {
			ipPos[ipInd[sz]=ind]=sz;
			Val[sz++]=w;
		}
	};

	/**
	 *  The function swaps (exchanges) indices (in this array) of two given elements.
	 *
	 * \param[in] ind1,ind2 indices of two distinct elements.
	 */
	__ALWAYS_INLINE void swap(int ind1, int ind2)
	{
		int pos1{ipPos[ind1]}, pos2{ipPos[ind2]};
		bool flag1{false}, flag2{false};
		if (pos1 >=0 && pos1 < sz)
			if (ipInd[pos1] == ind1)
				flag1=true;
		if (pos2 >=0 && pos2 < sz)
			if (ipInd[pos2] == ind2)
				flag2=true;
		if (flag1) {
			if (flag2) {
				double w{Val[pos1]}; Val[pos1]=Val[pos2]; Val[pos2]=w;
			}
			else
				ipPos[ipInd[pos1]=ind2]=pos1;
		}
		else if (flag2)
			ipPos[ipInd[pos2]=ind1]=pos2;
	};
};

#endif /* __SPARSEARRAY_H_ */
