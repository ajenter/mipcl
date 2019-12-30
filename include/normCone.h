///////////////////////////////////////////////////////////////
/**
 * \file normCone.h interface for `CNormCone` class
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
 *  modify it under the terms of the GNU Lesser General Public Licens
 *  ///////////////////////////////////////////////////////////////
/**
 * \file normCone.h interface for `CRecord` class
 * |  __Author__  | N.N. Pisaruk                              |
 * |-------------:|:------------------------------------------|
 * |  __e-mail__  | nicolaipisaruk@gmail.com                  |
 * | __home page__| wwww.mipcl-cpp.appspot.com                |
 *
 *   \copyright __2018 Nicolai N. Pisaruk__
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
#ifndef _NORMCONE_H_
#define _NORMCONE_H_

class CLP;

/**
 * This class supports inclusion into LP formulations norm constraints:
 *      `||x^J|| <= t`, where `J` is a subset of variable indices, `t` is a real variable.
 */
class CNormCone
{
private:
	int *m_ipBuf; ///< Memory buffer to store constraints.
	char *m_sMsg; ///< Pointer to `CLP::m_sWarningMsg`.
	int m_iCtrNum; ///< Number of norm constraints.
	int m_iBufSize; ///< Size of memory buffer (in integer words).
	int m_iSz; ///< Size (in integer words) of used buffer part.
public:
	/** Constructor.
	 *
	 * \param[in] maxCtrNum maximum number of norm constraints;
	 * \param[in] avCtrSize average norm constraint size;
	 * \param[in] msg pointer to message string (usually `CLP::m_sWarningMsg`).
	 */
	CNormCone(int maxCtrNum, int avCtrSize, char *msg);

	/// Destructor.
	virtual ~CNormCone();

// Implementation
	/**
	 * The procedure adds a new norm constraint.
	 *
	 * \param[in] t right-hand side variable;
	 * \param[in] sz number of variables involved into left-hand side;
	 * \param[in] ipVars list of  `sz` variables involved into left-hand side;
	 * \param[in] tol this norm constraint is considered violated at given point
	 * if violation value is greater than `tol`.
	 */
	int addCtr(int t, int sz, int *ipVars, double tol);

	/**
	 * Given a solution to a node LP, the procedure separates this solution point
	 * from the second-order cones given by the cones constraints.
	 *
	 * \param[in] fromLP if `true`, this procedure has been called from an `CLP` object, and otherwise from `CMIP` object;
	 * \param[in] lp pointer to calling `CLP` object;
	 * \param[in] dpX point to be separated;
	 * \param[in] ipCol,dpA memory buffers used for adding separating inequalities to the calling object.
	 * \param[in] genCuts if `true` separating inequalities (cuts) are added to the calling object;
	 * otherwise, the procedure returns a non-zero value if tested point violates at least one norm constraint,
	 * and it returns `0` if all norm constraints are satisfied at tested point;
	 * \param[in] maxCutNum if `k` norm constraints are violated, procedure add min{k,maxCutNum}` cuts.
	 * \return number of violated norm constraints.
	 *
	 */
	int separate(bool fromLP, CLP *lp, const double *dpX, int* &ipCol, double* &dpA, bool genCuts, int maxCutNum=1000000);

	/**
	 * During preprocessing, indices of variables involved into a norm-constraint may be changed.
	 * This procedure is used to update norm-constraint index set.
	 *
	 * \param[in] ipNewCol list of new indices.
	 */
	void updateCtrs(int *ipNewCol);

#ifdef __PYTHON_MODULE_
	/**
	 * It is not easy to pass arrays from a Phyton program to The MIPCL library,
	 * therefore, the variables involved into a norm constraint
	 *
	 * \param[in] t right-hand side variable;
	 * \param[in] sz number of variables in left-hand side;
	 * \param[in] tol constraint tolerance.
	 * \sa addCtr(), addVar().
	 */
	int startCtr(int t, int sz, double tol);

	/**
	 * The procedure adds a variable to the list of variables
	 * of the previously started norm constraint.
	 *
	 * \param[in] k number of variables in norm constraint list;
	 * \param[in] var index of added variable.
	 * \return `k+1`;
	 * \sa startCtr().
	 */
	int addVar(int k, int var);
#endif


};

#endif /*_NORMCONE_H_*/
