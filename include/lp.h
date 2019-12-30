///////////////////////////////////////////////////////////////
/**
 * \file lp.h Interface for `CLP` class
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
#ifndef __LP__H
#define __LP__H

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
typedef __int64 __LONG; ///< For Visual C++, `long` means 32 bit integer; so we redefine it.
#else
#ifndef MIP_API
#define MIP_API
#endif
typedef long long __LONG; ///< For GNU C++, `long` means 64 bit integer.
#endif

#include <fstream>
#include "sparseArray.h"

class CException;
class CPrepStack;
class CLU;
class CNormCone;

/// This class has been designed for solving Linear Programs (LPs)
/**
 * \f{align*}{
 *     & c^Tx \to \max,\\
 *     & b_1 \le Ax \le b_2,\\
 *     & l \le  x \le u,
 * \f}
 * by both, prime and dual, simplex methods.
 */
class MIP_API CLP 
{
    friend class CLU;
    friend class CMIP;
    friend class CNormCone;
    bool m_bCLP; ///< `true` for LPs, and `false` for MIPs.
public:
    /**
     * If preprocessing is switched on (default option)
     * or columns are generated, the solver do not preserve the order
     * in which the variables were added to the matrix.
     * The handles are primarily used for mapping the variables in the solver memory onto user's variables.
     */
    typedef int tagHANDLE; ///< program type for handles of variables and constraints

    /// Matrix alignment.
    /**
     * The non-zero matrix entries may be listed either _column-wise_ or _row-wise_.
     * It is also possible that the entries are not ordered at all.
     * Furthermore, even if the non-zero entries are ordered before the solution procedure starts,
     * after adding to the matrix new rows and/or columns, the previous alignment will be lost.
     * In the latter case, we say that the matrix is _approximately aligned_ either column-wise or row-wise.
     */
    enum enAlign {
        ALIGN_NONE        = 0x00000000, ///< no alignment
        ALIGN_COLUMN_APPR = 0x00000001, ///< mainly by columns
        ALIGN_COLUMN      = 0x00000011, ///< by columns
        ALIGN_ROW_APPR    = 0x00000100, ///< mainly by rows
        ALIGN_ROW         = 0x00001100  ///< by rows
    };

	/// LP types of constraints.
	/**
	 * The type of a constraint is the bitwise OR of the members of `enCtrType`.
	 */
    enum enCtrType {
        CTR_ATTACHED     = 0x00000001, ///< constrains of this type cannot be removed from matrix
        CTR_LEFT         = 0x00000002, ///< constraints is _left bounded_ (its left hand side is finite)
        CTR_RIGHT        = 0x00000004, ///< constraints is _right bounded_ (its right hand side is finite)
        CTR_EQ           = 0x00000008, ///< constraint is equation
        CTR_REF          = 0x00000010, ///< flag is _privately_ used in `CLP`
		CTR_NOT_STABLE   = 0x00000020, ///< With too big or too small coefficients, or has been derived from non-stable constraints.
		CTR_STRONG_CUT   = 0x80000000  ///< cutting planes (cuts)  with this flag are treated as ordinary inequalities
    };

	/// LP types of variables.
	/**
	 * The type of a variable is the bitwise OR of the members of `enVarType`.
	 */
    enum enVarType {
        VAR_LEFT        = 0x10000000,  ///< variable is _lower bounded_ (its lower bound is finite)
        VAR_RIGHT       = 0x20000000,  ///< variable is _upper bounded_ (its upper bound is finite)
        VAR_FX          = 0x40000000,  ///< variable is _fixed_
        VAR_REF         = 0x80000000,  ///< is _privately_ used in `CLP`
        VAR_NOT_MOD     = 0x04000000,  ///< a variable with this flag set cannot be modified in any way (deleted or scaled)
		VAR_INT         = 0x00001000, ///< integral variable, meaningless in `CLP`.
   };

    /// This is the list of all possible ways to scale matrix.
    enum enScaling {
    	/**
    	 * A row (resp., column) scaling factor is the power of 2 that is closest to
    	 *  the average value of the minimum and maximum absolute values of coefficients in that row (resp., column).
    	 *  The min-max scaling procedure 10 times, in turn, first computes row scaling factors and multiplies all rows by their factors,
    	 *  and then computes column scaling factors  and multiplies all the columns by their factors.
    	 */
        SCL_MINMAX, ///< Min-max scaling.
        /**
         * The _ideal scaling_ procedure applies a complex algorithm to compute scaling factors such that, after scaling,
         *  the difference between the maximum and minimum of the absolute values of matrix coefficients is minimal.
         */
        SCL_IDEAL, ///< Ideal scaling.
        SCL_NO, ///< Switches off scaling.
        /**
         * A row (resp., column) scaling factor is the power of 2 that is closest to
         *  the geometric mean value of the absolute values of coefficients in that row (resp., column).
         *  The _geometric-mean scaling_ procedure 10 times, in turn, first computes column scaling factors and multiplies all columns by their factors,
         *  and then computes rows scaling factors  and multiplies all rows by their factors.
         */
        SCL_GM_ROWS, ///< Geometric-mean scaling, rows are scaled after columns.
        SCL_GM_COLUMNS, ///< Geometric-mean scaling, rows are scaled before columns (for details, see `SCL_GM_ROWS`).
        /**
         * A row (resp., column) scaling factor is the power of 2 that is closest to
         *  the maximum absolute values of coefficients in that row (resp., column).
         *  The _ max-value scaling_ procedure 10 times, in turn, first computes column scaling factors and multiplies all columns by their factors,
         *  and then computes rows scaling factors  and multiplies all rows by their factors.
         *
         */
        SCL_MAX_ROWS, ///< Max-value scaling, rows are scaled after columns.
        SCL_MAX_COLUMNS, ///< Max-value scaling, rows are scaled before columns (for details, see `SCL_MAX_ROWS`).
        NOT_SCALED=127, ///< When adding to the matrix a new constraint (or column), the value of scaling factor set to `NOT_SCALED` means that this constraint (or column) is to be scaled by the solver.
		SCL_MIN_EXP=-19 ///< When scaling, entries with exponents less than `SCL_MIN_EXP` are treated as zeroes.
    };

    /// LP methods.
    enum enLPmethod {
        AUTO_DETECT, ///< It is left to the solver to decide what method to use for solving LPs.
        PRIME_SIMPLEX, ///< Prime simplex is used for solving LPs.
        DUAL_SIMPLEX ///< Dual simplex is used for solving LPs.
    };

    ///Separation rules.
    enum enSepRule {
        SEP_MOST_VIOLATED,  ///< Most violated constraint is chosen.
        SEP_ONLY_EQUATIONS, ///< Violated constraint is chosen only among equations.
		/**
		 * Steepest edge rule, i.e.,
		 * violation of each constraint is divided by the norm of corresponding row,
         * and then constraint with minimum ratio is chosen.
		 */
        SEP_STEEPEST_EDGE
    };

    /// Pricing rules.
    enum enPricingRule {
        PRC_MOST_NEGATIVE, ///< Variable with most negative reduced cost is chosen.
        PRC_STEEPEST_EDGE ///< Reduced cost of each variable is divided by the norm of corresponding column, and then variable with minimum ratio is chosen.
    };

    /// Problem states.
    /**
     * Any current state of the problem being solved is composed as the bitwise OR of the members of `enProbState`.
     */
    enum enProbState {
        PROB_IN_MEMORY    = 0x00000001, ///< Memory necessary to solve problem has been allocated.
        PROB_PREPARED     = 0x00000002, ///< Ready (memory allocated, matrix preprocessed and scaled) for optimizing.
        PROB_SOLVED       = 0x00000004, ///< Problem has been solved.
        PROB_SOLUTION     = 0x00000008, ///< Solution has been found.
        PROB_OPTIMAL      = 0x00000010, ///< Optimality has been proven.
        PROB_INFEASIBLE   = 0x00000020, ///< Constraint inconsistency has been proven.
        PROB_UNBOUNDED    = 0x00000040, ///< Objective function is not bounded.
        PROB_TIME_LIMIT   = 0x00000080, ///< Time limit exceeded.
        PROB_IT_LIMIT    = 0x00000100, ///< Limit on number of iterates exceeded.
        PROB_GAP_LIMIT    = 0x00000200, ///< Required integrality gap attained.
        PROB_SOLVER_FLAGS = PROB_SOLVED|PROB_SOLUTION|PROB_OPTIMAL|PROB_INFEASIBLE|PROB_UNBOUNDED|PROB_TIME_LIMIT|PROB_GAP_LIMIT ///< Bitwise OR of all `enProbState` masks.
    };

    /**
     * A bitwise ORing of the members of `enRowColGenRule1` stored in `m_uRowColGenRule` determines whether rows and/or columns are generated.
     */
    enum enRowColGenRule1 {
    	ROW_GEN     	= 0x1, ///< Rows/cuts can be generated by overloading `separate()`, `genCut1()`, or cut generating procedures implemented in `CMIP`.
    	COL_GEN     	= 0x2,  ///< Columns can be generated by overloading `generateColumns()`.
    	SEP_PROC	    = 0x4 ///< If `separate()` has been overloaded.
    };

// Constants
    static const double INF; ///< default value of infinity.
    static const double VAR_INF; ///< default value of infinity for variables.
    static const double MAX_BIG_M; ///< default value for big `M` used  in the one phase implementation of the prime simplex algorithm.
    static const int SHIFT; ///< equals to (std::numeric_limits<int>::max() >> 1).
    static const char msgNoSolution[]; ///< constant string for no-solution message.

    /**
     * \warning Each message must contain no more than `255` characters.
     */
    char m_sWarningMsg[256]; ///< Character array for using as message strings.

protected:
	int m_iThread; ///< Index of a thread in a multi-threaded MIP-application.
    char m_strProblemName[32]; ///< Stores problem name.

    std::ostream* m_pLogStream; ///< LOG file.

//////  start M A T R I X
    bool m_bSense; ///< Objective sense: if `true` - maximize, if `false` - minimize.
    int m_iN; ///< Number of variables in the currently processed LP.
    int m_iN0; ///< Number of variables in the root LP; if columns are generated, `m_iN0 <= m_iN`.

    /*
     *  If an extra variable (above this limit) is added, memory for the matrix is reallocated.
     */
    int m_iNmax; ///< Maximum number of variables.
    int m_iM;  ///< Number of constraints in the currently processed LP.
    int m_iM0; ///< Number of constraints in the root LP; if rows (cuts) are generated, `m_iM0 <= m_iM`.

    /**
     *  If an extra constraint (above this limit) is added, memory for the matrix is reallocated.
     */
    int m_iMmax; ///< Maximum  number of constraints.
    int m_iNZ;  ///< Number of nonzero entries in the currently processed matrix.
    int m_iNZ0; ///< Number of nonzero entries in the root LP matrix.

    /**
     *  If an extra entry (above this limit) is added, memory for the matrix is reallocated.
     */
    int m_iNZmax; ///< Maximum number of entries.

    double* m_dpC; ///< Cost vector, array of size `m_iNmax`.
    double* m_dpC0; ///< Pointer to stored objective vector. To prevent cycling, `dualSimplex()` perturbs the objective; `m_dpC0` is used to restore the original objective.

    /**
     * Cost vector `m_dpC` is multiplied by 2^{m_iObjScaleExp}.
     */
    int m_iObjScaleExp; ///< objective scaling factor.

    /**
     * Array of size `2*m_iMmax`, `m_dpB[i<<1]` and `m_dpB[(i<<)+1]` are left and right hand sides of constraint `i`.
     */
    double* m_dpB; ///< left and right hand sides of constraints.

    /*
     *  Array of size `2*m_iNmax`, `m_dpD[j<<1]` and `m_dpD[(j<<1)+1]` are lower and upper  bound of variable `j`.
     */
    double* m_dpD; ///< lower and upper bounds of variables.

private:
    char* m_cpMatrix; ///< pointer to the memory allocated for the matrix.

    CNormCone *m_pNormCone; ///< pointer to the norm-constraint storage.

protected:
    int* m_ipPrevColEntry; ///< array of size `m_iNZmax`, `m_ipPrevColEntry[j]` is the entry in column `m_ipCol[j]` preceding entry `j`.
    int* m_ipPrevRowEntry; ///< array of size `m_iNZmax`, `m_ipPrevRowEntry[j]` is the entry in row `m_ipRow[j]` preceding entry `j`.
    int* m_ipRow; ///<  array of size `m_iNZmax`, `m_ipRow[j]` is row index of entry `j`.
    int* m_ipCol; ///<  array of size `m_iNZmax`, `m_ipCol[j]` is column index of entry `j`.
    double* m_dpVal; ///< array of size `m_iNZmax`, `m_dpVal[j]` is coefficient for entry `j`.
// Row lists
    int* m_ipLastRowEntry; ///< array of size `m_iMmax`, `m_ipLastRowEntry[i]` starts the linked list of entries in row `i`.
    int* m_ipRowSize; ///< array of size `m_iMmax`, `m_ipRowSize[i]` is number of entries in row `i`.
// Column Lists
    int* m_ipLastColEntry; ///< array of size `m_iNmax`, `m_ipLastColEntry[i]` starts the linked list of entries in column `i`.
    int* m_ipColSize; ///< array of size `m_iNZmax`, `m_ipColSize[i]` is number of entries in column `i`.

    /**
     * Array of size `m_iMmax`, `s=m_cpRowScale[i]` is a scale factor of row `i`, i.e., row `i` is multiplied by 2^s.
     */
    signed char* m_cpRowScale; ///< row factors.

    /**
     *  Array of size `m_iNmax`, `s=m_cpColScale[j]` is a scale factor of column `j`, i.e., column `j` is multiplied by 2^s.
     */
    signed char* m_cpColScale; ///< column factors.

    enAlign m_iAlign; ///< specifies how matrix entries are listed.

    /**
     *  Array of size `m_iMmax`,
     * `m_ipRowHd[i]` is _handle_ of row (or constraint) `i`,
     *     - if `m_ipRowHd[i] < 0`, constraint `i` is also stored in the pool;
     *     - if `m_ipRowHd[i] >= 0`, it is up to the user to assign a value to the handle, say,
     *       this may be an index specifying the location of a data structure describing constraint `i`.
     *
     * \attention each row handle must be unique.
     */
    tagHANDLE* m_ipRowHd; ///< row handles.

    /**
     *  Array of size `m_iNmax`,
     * `m_ipColHd[j]`     CNormCone *m_pNormCone;

     * is _handle_ of column (or variable) `j`,
     *     - if `m_ipColHd[j] < 0`, column `i` is also stored in the pool;
     *     - if `m_ipColHd[j] >= 0`, it is up to the user to assign a value to the handle, say,
     *       this may be an index specifying the location of a data structure describing column `j`.
     *
     * \attention each column handle must be unique.
     */
    tagHANDLE* m_ipColHd; ///< column handles.

    /**
     *  Array of size `m_iNmax`, `m_ipVarType[j]` is a type of variable `j`,
     *  which is the bitwise OR of the members of two enumerations CLP::enVarType and CMIP::enVarType.
     */
    unsigned int* m_ipVarType;  ///< types of variables.

    /**
     * Array of size `m_iMmax`, `m_ipCtrType[i]` is a type of constraint `i`,
   	 * which is the bitwise OR of the members of two enumerations `CLP::enCtrType` and `CMIP::enCtrType`.
     */
    unsigned int* m_ipCtrType; ///< constraint types.

//////  end of M A T R I X

// Tolerance parameters

    double m_dZeroEntry; ///< all numbers that are less than `m_dZeroEntry` are treated as zeroes.
    double m_dVarInf; ///< `-m_dVarInf <= x[j] <= m_dVarInf`, default value is `CLP::VAR_INF=1.0e+12`.
    double m_dCtrTol; ///< maximum violation allowed for constraints.
    double m_dVarTol; ///< maximum violation allowed for variables.
    double m_dShadowPriceTol; ///< maximum violation allowed for shadow prices (dual variables).
    double m_dRedCostTol; ///< maximum violation allowed for reduced costs.

private:
    double m_dPivotTol; ///< minimum value for pivots.
    double m_dRelPivotTol; ///< pivot cannot be less than maximum absolute value in the pivot row times `m_dRelPivotTol`.
    double m_dGoodPivot; ///< minimum value for "good" pivots.
    double m_dPivotErrTol; ///< if the difference of two pivot values, one from the pivot column and the other one from the pivot row, are greater than `m_dPivotErrTol`, then the basic matrix is refactored.
    double m_dRelPivotErrTol; ///< if the ratio of two pivot values, one from the pivot column and the other one from the pivot row, are greater than `m_dRelPivotErrTol`, then the basic matrix is refactored.
    double m_dTolPrimeDeg; ///< if the increase of the objective in some iterate of the `primalSimplex()` is less than `m_dTolPrimeDeg`, then this iterate is treated as degenerate.
    double m_dTolDualDeg; ///< if the increase of the objective in some iterate of the `dualSimplex()` is less than `m_dTolDualDeg`, then this iterate is treated as degenerate.
    double m_dBigM; ///< maximum value for big Ms.

    /**
     * If a bound for a variable is violated by more than `m_dVarViolatThreshold`,
     * the separation procedure will not verify any non-trivial inequalities.
     */
    double m_dVarViolatThreshold;

    /**
     * If a dual variable is less than `-m_dShadowPriceViolatThreshold`,
     * the corresponding column is chosen as a pivot column;
 	 * the latter is not true for the steepest edge strategy.
     */
    double m_dShadowPriceViolatThreshold;

    // Memory for the matrix
    char *m_cpColMatrMem; ///< memory for storing arrays `m_dpC`,  `m_dpD`, `m_ipLastColEntry`, `m_ipColSize`, `m_ipVarType`, `m_ipColHd`, `m_cpColScale`.
    char *m_cpRowMatrMem; ///< memory for storing arrays `m_dpB`, `m_ipLastRowEntry`, `m_ipRowSize`, `m_ipCtrType`, `m_ipRowHd`, `m_cpRowScale`.

    char *m_cpColSimplexMem; ///< memory for storing arrays `m_dpPi`, `m_dpUd`, `m_ipBasicColumn`, `m_ipColMap`.
    char *m_cpRowSimplexMem; ///< memory for storing arrays `m_dpBeta`, `m_dpX`, `m_dpUb`, `m_ipRowMap`, `m_ipBasicRow`.

protected:
    double m_dFixedCost; ///< total cost of all variables fixed during preprocessing.
    double m_dObjFactor; ///< original objective was divided by m_dObjFactor;
    double m_dBasicVarCost; ///< total cost of basic variables.
    double m_dObjVal; ///< current value of the objective
    double* m_dpX; ///< prime solution; array of size `m_iMmax.
    double* m_dpPi; ///< dual solution (potentials), array of size `m_iNmax`.

    /**
     * `m_iState` is the bitwise OR of the members of the enumeration `enProbState`.
     */
    int m_iState; ///< state of the problem.

    /**
     * Problem is
     *  - infeasible if m_iS != 0 && m_bDualFeasible,
	 *  - unbounded if m_iS >= 0 && m_bPrimeFeasible
     */
    int m_iS; ///< Defines a state of the last LP solved.

// The following 5 auxiliary arrays are extensively used as temporal memory in different ways.
    double *m_dpUd; ///< auxiliary array of size `m_iNmax`.
    double *m_dpUb; ///< auxiliary array of size `m_iMmax`.
    int* m_ipArray; ///< auxiliary array of size max{`m_iNmax,m_iMmax`}
    double* m_dpArray; ///< auxiliary array of size `m_iMmax+m_iNmax`
    double* m_dpW; ///< auxiliary array of size `2*m_iMmax`
	double *m_dpFd; ///< Auxiliary array of size max{`2*m_iNmax`}, this array may be freely used in user defined cut or column generating procedures.
	int *m_ipFd; ///< Auxiliary array that shares memory with `m_dpFd'.
	double *m_dpFb; ///< Auxiliary array of size `max{2*m_iMax,m_iNmax}`, this array may be freely used in user defined cut or column generating procedures.

    /**
     * Array of size max{`m_iNmax,m_iMmax+m_iMaxLuUpdateNum`},
     *  may be freely used in user cut and column generating procedures.
     */
    double *m_dpNorm; ///< row or column norms in the current basis.

private:
    SparseArray m_spUb;  ///< Sparse analog for `m_dpUb`.
    SparseArray m_spUd;  ///< Sparse analog for `m_dpUd`.
    SparseArray m_spGamma; ///< Sparse analog for `m_dpGamma`.

///////////////////////////////////////////////////////////////
// F L A G S
/////////////
public:
    bool m_bPrimeFeasible; ///< if `true`, the current basic solution is prime feasible.
    bool m_bDualFeasible; ///< if `true`, the current basic solution is dual feasible.
    bool m_bRowNorm;  ///< if `true`, `m_dpNorm[i]` is the norm of `m_ipBasicRow[i]` in the current basis
    bool m_bColNorm; ///< if `true`, `m_dpNorm[i]` is the norm of `i`-th non-basic row.
private:
    enScaling m_eScaling; ///< scaling strategy.

    /**
     * If `m_iInfoMsgFlag & 0x1`, run time CLP-messages are blocked;
     * if `m_iInfoMsgFlag & 0x2`, run time branch-and-cut CMIP-messages are blocked;
     * if `m_iInfoMsgFlag & 0x4`, cut-statistics CMIP-messages are blocked;
     * if `m_iInfoMsgFlag & 0x10`, all run time messages are blocked.
     *
     * \sa `setInfoMsgState()` and `beSilent()`.
     */
    int m_iInfoMsgFlag;

    int m_iRowNormNum; ///< row norms are computed for all rows plus for those columns that are in current row basis.
    int m_iColNormNum; ///< number of computed column norms.

    /**
     * This constant is used when implementing the steepest edge pricing and separating rules.
     * If the recomputed norm of a pivot row or column exceeds its actual value
     * by more than `BAD_NORM_FACTOR` times, all row or column norms are recomputed.
     */
    static const double BAD_NORM_FACTOR;

    /**
     *memory to store the objective and the  row and variable bounds
     *before applying perturbations to prevent cycling
     */
    char *m_cpPrtbMem; ///< only a pointer, memory is allocated in case of necessity.

    /**
     * If an approximate norm (square of the norm) is BAD_NORM_FACTOR times bigger or smaller of
 	 * the exact norm, m_iBadNormNum is incremented
 	 * if `m_iBadNormNum > MAX_BAD_NORM_NUM`, all norms are recomputed.
     */
    int m_iBadNormNum; ///< Number of error encountered among the norms.

    /**
     * If m_iNormUpdatesNum > MAX_UPDATE_NORM_NUM, switch to Dancig rule.
     */
    int m_iNormUpdatesNum; ///< number of norm updates.
//
    double m_dPivot; ///< value of pivot element when pivot operation is performed.
    double m_dObjInc; ///< increase of prime or dual objective.
    double m_dPrimalStep; ///< step value when prime solution is updated.
    double m_dDualStep; ///< step value when dual solution is updated.

//////////////////////////////////////////////////////////////
// BASIS
    int m_iNoOfTrivial; ///< in current basis, there are `m_iNoOfTrivial` trivial (representing bounds for variables) rows among basic rows.

    int m_iBasicRowSize; ///< total number of non-zeroes in basic rows.
    int m_iBasicColSize; ///< total number of non-zeroes in basic columns.
protected:
    int m_iBasisSize; ///< number of nontrivial (b1 <= Ax <= b2) inequalities in the basis.

    /**
     * If `(r=m_ipBasicRow[i]) >= 0`, row `r` is the `i`-th row in the basis;
     * if `(r=m_ipBasicRow[i]) < 0`, row e_{-r-1} is the `i`-th row in the basis.
     */
    int *m_ipBasicRow;  ///< list of basic rows.

    /**
     * If `(k=m_ipRowMap[i]) < SHIFT`,  then `m_ipBasicRow[abs(k)-1]=i` and the value of the inequality is equal to its
     *   - left-hand-side if `k < 0`;
     *   - right-hand-side if `k > 0`.
     *
     *   If `k >= SHIFT`, then `m_ipBasicRow[k-SHIFT]=i`.
     *  \remark The pair (`m_ipRowMap`,`m_ipColMap`) is called a _short basis_.
     * \sa m_ipBasicRow.
     */
    int *m_ipRowMap;

    int *m_ipBasicColumn; ///< list of basic columns.

    /**
     * If `(k=m_ipColMap[j]) < SHIFT`,  then `m_ipBasicColumn[abs(k)-1]=j` and the value of the value of variable is equal to its
     *   - lower bound if `k < 0`;
     *   - upper bound if `k > 0`.
     *m_dRelPivotTol
     *   If `k >= SHIFT`, then `m_ipBasicColumn[k-SHIFT]=j`.
     *  \remark The pair (`m_ipRowMap`,`m_ipColMap`) is called a _short basis_.
     * \sa `m_ipBasicColumn`.
     */
    int *m_ipColMap;

    double *m_dpBeta; ///< values of non-basic rows, the value of row `i` is A_ix.

    enLPmethod m_eDefLPmethod; ///< LP method to solve the root LP.
    enLPmethod m_eLPmethod; ///< currently used LP method.

    // Parameters
    unsigned m_uRowColGenRule; ///< bitwise OR of the members of `enRowColGenRule`; by default, all flags are set to `1`.

    enSepRule m_eSepRule; ///< currently used separation rule.
    enPricingRule m_ePricingRule; ///< currently used pricing rule.
    enSepRule m_eLpSepRule; ///< separation rule used for solving pure LPs as well as root LPs when solving MIPs.
    enPricingRule m_eLpPricingRule; ///< pricing rule used for solving pure LPs as well as root LPs when solving MIPs.

protected:
    int m_iPreproc; ///< if `true` preprocessing is on; default value is `true`.
    CPrepStack* m_pPrepStack; ///< pointer to the preprocessing stack.

    /**
     * Sometimes we need to solve an LP that is a slightly modified version of  previously another LP.
     * If those modification does not affect the basic submatrix and the objective,
     * one can apply the dual simplex algorithm to re-optimizing the LP starting from an optimal basis of the original LP.
     * But if the bounds of some non-basic variables were modified, we need to update an initial basic solution as well.
     */
    int m_iUpdatePrimeSol; ///< number of modified bounds of non-basic variables.

// Factoring the basis
    CLU* m_pLU; ///< pointer to an object of class CLU, that solves linear systems.
private:
    int m_iMaxLuUpdateNum; ///< maximal number of multiplicators allowed.

// STATISTICS
    int m_iItNum; ///< number of iterates performed by `primeSimplex()` or `dualSimplex()`.
    int m_iDegItNum; ///< number of degenerate  iterates performed by `primeSimplex()` or `dualSimplex()`.
    int m_iLpItToInform; ///< number of LP iterations to print an info message.
    int m_iPrtbNum; ///< number of objective(dual)/bounds(primal) perturbations.
    int m_iPartitionNum; ///< number of times the basis was refactored during a call to `primeSimplex()` or `dualSimplex()`.

protected:
    __LONG m_lStartTime; ///< start time of solving problem.
    __LONG m_lSolTime; ///< solution time.

//////////////////////////////////////////////////////////////////////////
//                        Functions of general use                      //
//////////////////////////////////////////////////////////////////////////
private:
    /**
     * This function resets all tolerance parameters to their default values.
     */
    void resetToleranceParameters();

    void closeLogStream(); ///< closes LOG-file.

public:
	__LONG getStartTime() const
		{return m_lStartTime;} ///< \return start time (given in seconds since epoch).

	__LONG getSolTime() const
		{return m_lSolTime;} ///< \return solution time (in seconds).

    /**
     * \return `true` if this function is called from an object of type `CLP`,
     *          and `false` if it is called from an object of type `CMIP`.
     */
    bool isCLP()
    	{return m_bCLP;}

    /**
     * The procedure opens a log file, which name is the value of `name` appended by `_mipcl.log`.
     *
     * \param[in] name log file name; if `name=0`, `m_strProblemName` is used as log file name.
     */
    void openLogStream(const char* name=0);

    /**
     * If the lower bound of a variable is less than `-val`, then the variable is assumed not bounded from below.
     * Similarly, if the upper bound of a variable is greater than `val`, then the variable is assumed not bounded from above.
     * \param[in] inf new value of "infinity" for variables.
     * \sa getVarInf().
     */
    void setVarInf(double inf)
        {m_dVarInf=inf;}

    /**
     * \return currently used value of "infinity" for variables.
     * \sa setVarInf().
     */
    double getVarInf() const
        {return m_dVarInf;}

    /**
     * __MIPCL__ writes to the LOG file, which name is the concatenation of the problem name and the extension ".log"
     * quite a lot of information describing all phases of the solution process. After having been analyzed the contents
     * of the LOG file, the user may adjust the parameters of the solver to improve its efficiency.
     * The user can also write to this LOG file any messages by calling `writeStrToLogStream()`.
     * \param[in] str message string.
     */
    void writeStrToLogStream(const char* str);


///////////////////////////////////////////////////////////////////////////////
    /**
     * Sets problem name.
     * \param[in] name problem name.
     * \remark The name is truncated to have at most 31 symbol.
     */
    void setProblemName(const char* name);

    /**
     * \param[out] name problem name.
     * \attention allocate for `name` at least 32 bytes.
     */
    void getProblemName(char* name) const;

//    std::ostream* getLogStreamPtr()
//        {return m_pLogStream;}
//    void setLogStreamPtr(std::ostream* pLogStream)
//        {m_pLogStream=pLogStream;}

/////////////////////////////////////
// Getting and setting tolerance variables
    /**
   	 * All numeric values that are less than `m_dZeroEntry` are considered as zeroes.
   	 * \param[in] zero new value of "zero".
   	 * \remark The default value is  1.0e-13.
   */
    void setZero(double zero)
        {m_dZeroEntry=zero;}

    /**
     * \return currently used value of "zero".
     */
    double getZero() const
        {return m_dZeroEntry;}

    /**
     * The absolute value of any pivot element must always be greater than `m_dPivotTol` (an internal MICL parameter),
     * and it also must be greater than either `m_dGoodPivot` or `m_dRelPivotTol * gamma`,
     * where gamma is the maximum absolute value of a nonzero entry in the pivot column
     * for the prime simplex method, or in the pivot row for the dual simplex method.
     * The other two parameters, `m_dPivotErrTol` and `m_dRelPivotErrTol`
     * are used to estimate absolute and relative errors when computing pivot elements.
     *
     * Default values:
     *     - 1.0e-9 for `m_dPivotTol`,
     *     - 1.0e-3 for `m_dGoodPivot`,
     *     - 1.0e-10 for `m_dPivotErrTol`
     *     - 1.0e-7 for `m_dRelPivotTol` and `m_dRelPivotErrTol`.
     *
     * \param[in] tolPiv new value for `m_dPivotTol`;
     * \param[in] relTolPiv new value for `m_dRelPivotTol`;
     * \param[in] goodPiv new value for `m_dGoodPivot`;
     * \param[in] tolPivErr new value for `m_dPivotErrTol`;
     * \param[in] relTolPivErr new value for `m_dRelPivotErrTol`;
     * \attention All these tolerance parameters are interdependent; therefore, tuning them is a subtle task.
     * \sa getPivTol(), getGoodPiv(), getRelPivTol(), getPivErrTol(), getRelPivErrTol()
     */
    void setPivTol(double tolPiv, double relTolPiv, double goodPiv,
                   double tolPivErr, double relTolPivErr);

    /**
     * \return minimum absolute value of pivots.
     * \sa setPivTol().
     */
    double getPivTol() const
        {return m_dPivotTol;}

    /**
     * \return minimum absolute value of pivots that are considered as "good".
     * \sa setPivTol().
     */
   double getGoodPiv() const
        {return m_dGoodPivot;}

   /**
    * \return minimum relative pivot value in pivot row and column.
     * \sa setPivTol().
    */
    double getRelPivTol() const
        {return m_dRelPivotTol;}

    /**
     * \return tolerance for relative pivot errors.
      * \sa setPivTol().
     */
    double getPivErrTol() const
        {return m_dPivotErrTol;}

    /**
     * \return tolerance for pivot errors.
      * \sa setPivTol().
     */
    double getRelPivErrTol() const
        {return m_dRelPivotErrTol;}

    /**
     * \param[in] tol any variable may maximum violation allowed.
     * \sa getVarTol().
     */
   void setVarTol(double tol);

   /**
    * \return maximum violation allowed for bounds on variables.
    */
   double getVarTol() const
        {return m_dVarTol;}

    /**
     * \param[in] tol maximum violation allowed for regular constraints (not cuts).
     * \sa getCtrTol().
     */
    virtual void setCtrTol(double tol);

    /**
     * \return maximum violation allowed for regular constraints (not cuts).
     */
    double getCtrTol() const
        {return m_dCtrTol;}

    /**
     * _Shadow prices_ (dual variables) are considered as non-negative if their values are not less than `-m_dShadowPriceTol`.
     * \param[in] tol new value for `m_dShadowPriceTol`.
     * \sa getShadowPriceTol().
    */
    void setShadowPriceTol(double tol);

    /**
     * \return currently used value for shadow price tolerance.
     * \sa setShadowPriceTol().
    */
    double getShadowPriceTol() const
        {return m_dShadowPriceTol;}

    /**
     * _Reduced costs_ of variables are considered as non-negative if their values are not less than `-m_dRedCostTol`.
     * \param[in] tol new value for `m_dRedCostTol`.
     * \sa getRedCostTol().
    */
    void setRedCostTol(double tol);

    /**
     * \return currently used value for reduced cost tolerance.
     * \sa setRedCostTol().
    */
    double getRedCostTol() const
        {return m_dRedCostTol;}

// setting/getting parameters
    /**
     * The function sets _degeneracy tolerance_ values for both, prime and dual` simplex procedures.
     * \param[in] primeDegTol  if increase of objective after performing any iterate of `primeSimplex()`
     *  is less than `primeDegTol` this iterate is considered as _degenerate_.
     * \param[in] dualDegTol  if decrease of objective after performing any iterate of `dualSimplex()`
     *  is less than `dualDegTol` this iterate is considered as _degenerate_.
     *  \remark If one of two input parameters has a negative value, the corresponding tolerance value remains unaffected.
     *  \sa getPrimeDegTol(), getDualDegTol().
     */
    void setDegenTol(double primeDegTol, double dualDegTol);

    /**
     * \return degeneracy tolerance values for prime simplex algorithm.
     * \sa setDegenTol().
     */
    double getPrimeDegTol() const
        {return m_dTolPrimeDeg;}

    /**
     * \return degeneracy tolerance values for dual simplex algorithm.
     * \sa setDegenTol().
     */
    double getDualDegTol() const
        {return m_dTolDualDeg;}
    
    /**
     * If a bound for a variable is violated by more than a given value,
     * the separation procedure will not verify any non-trivial inequalities.
     * \param[in] th tolerance value.
     * \sa getVarViolatThreshold().
     */
    void setVarViolatThreshold(double th)
    	{m_dVarViolatThreshold=th;}

    /**
     * If a bound for a variable is violated by more than a given value, called variable violation threshold,
     * the separation procedure will not verify any non-trivial inequalities.
     * \return variable violation threshold.
     * \sa setVarViolatThreshold().
     */

    double getVarViolatThreshold() const
    	{return m_dVarViolatThreshold;}

    /**
     * If a dual variable is less than a given value,
     * the corresponding column is chosen as a pivot column;
 	 * the latter is not true for the steepest edge strategy.
     * \param[in] th tolerance value.
     */
    void setShadowPriceToleranceThreshold(double th)
    	{m_dShadowPriceViolatThreshold=th;}

    /**
     * \param[in] scalingMethod new scaling method.
     * \sa enScaling.
     */
    void setScaling(enScaling scalingMethod)
        {m_eScaling=scalingMethod;}

    /**
     * Call this function to switch off all preprocessing actions, excluding scaling the matrix.
     */
    void preprocOff()
        {m_iPreproc=0;}

    /**
     * \return `true` if CLP run time messages are to be printed.
     * \sa switchLpInfoMsg()
     */
    bool lpInfoMsg() const
    {return (!(m_iInfoMsgFlag & 0x11))? true: false;}

    /**
	 * The function is called to switch on or off printing CLP run time messages.
     * \param[in] flag if set to `true`, run time CLP-messages will be printing; otherwise, not,
     */
    void switchLpInfoMsg(bool flag)
    {
    	if (!flag) m_iInfoMsgFlag|=0x1; else m_iInfoMsgFlag&=~0x1;
    }

    /**
     * \param[in] flag if set to `true`, all run time messages are blocked.
     */
    void beSilent(bool flag=true)
    {if (flag) m_iInfoMsgFlag|=0x10; else m_iInfoMsgFlag&=~0x10;}

    bool isSilent() const
    {return (m_iInfoMsgFlag & 0x10)? true: false;} ///< \return `true` if silent mode has been set.

    /**
     * \param[in] fr CLP-messages will be printed after those LP-iterates which index is a multiple of `fr`.
     * \sa switchPrintingInfoMsg().
     */
    void setFrequencyForInfoMsg(int fr)
        {if (fr > 0) m_iLpItToInform=fr;}

//    /// \param[out] str running time.
//    void getRunTime(char* str);

    /**
     * Some modification of the matrix may make the currently optimal basis either prime or dual infeasible.
     * Therefore, after having done any such modifications, the user must call this function to
     * change the status of the basis.
     * \param[in] primeFeasible if `false`, the basis becomes prime infeasible;
     * \param[in] dualFeasible if `false, the basis becomes dual infeasible;
     */
    void setOptFlags(bool primeFeasible, bool dualFeasible)
        {m_bPrimeFeasible=primeFeasible; m_bDualFeasible=dualFeasible;}

    /**
     * \param[in] dpC array of size `m_iN` of objective coefficients;
     * \param[in] sense  if `true`, the goal is to _maximize_ the objective function; otherwise, to _minimize_;
     * \param[in] scale if `true`, the objective is given in the scaled variables; otherwise, in non-scaled.
     */
    void setObjective(const double *dpC, bool sense=true, bool scale=true);

    /**
     * \param[in] i if `i < SHIFT`,  the objective is row `i`; otherwise,
     *  if `i >= SHIFT`, the objective is to optimize the (i-SHIFT)-th variable;
     * \param[in] sense  if `true`, the goal is to _maximize_ the objective function; otherwise, to _minimize_.
     */
    void setObjective(int i, bool sense=true);

    /**
     * \param[in] sense  if `true`, the goal is to _maximize_ the objective function; otherwise, to _minimize_.
     */
    void setObjSense(bool sense)
    	{m_bSense=sense;}

    /**
     * \return `true` if the goal is to _maximize_ the objective function,
     *  or `false` if the goal is to _minimize__ the objective.
     */
    bool getObjSense() const
    {return m_bSense;}

//////////////////////////////////////////////////////////////////////////
//                           Matrix attributes                          //
//////////////////////////////////////////////////////////////////////////

    /// \return number of variables.
    int getVarNum() const
        {return  m_iN;}

    /// \return number of constraints.
    int getCtrNum() const
        {return m_iM;}

    /// \return number of nonzero entries in the matrix.
    int getNonZerosNum() const
        {return m_iNZ;}

// Separation/Pricing rule
    /** \return `true` if rows can be generated.
     * \sa `enRowGenRule`.
     */
    bool isRowGen() const
    {return (m_uRowColGenRule & ROW_GEN)? true: false;}

    /** \return `true` if columns can be generated.
     * \sa `enRowGenRule`.
     */
   bool isColGen()
    {return (m_uRowColGenRule & COL_GEN)? true: false;}

protected:

   /**
    * \return tolerance value for a given constraint; for LPs, this is always `m_dCtrTol`.
    */
   virtual double getThisCtrTol(int /*i*/) const
   	   {return m_dCtrTol;}

    /**
     * \return `true` if the problem is an LP, otherwise, `false`.
     * \remark Default implementation always returns `true`.
     */
    virtual bool isPureLP() const
        {return true;}
//////////////////////////////////////////////////////////////////////////
//                       C O N S T R U C T O R S                        //
//////////////////////////////////////////////////////////////////////////
//private:
//    void init(); ///< Initializes class members.
public:
    /**
     *  The constructor initializes an "empty" `CLP` object and sets the problem name to `name`.
     *  \param[in] name problem name of no more than '30' characters; if `name==0`, problem name will be a featureless "LP".
     */
    CLP(const char* name);
#ifndef __ONE_THREAD_
    /// Clone constructor.
    /**
     * This constructor is used in multithreaded MIP applications.
     * \param[in] other reference to `CLP` object to be cloned;
     * \param[in] thread thread index (positive integer).
     * \note The clone constructors in any derived class first must call the clone constructor of the base class.
     * \sa `CMIP::CMIP(const CMIP &other, int thread)`, `CMIP::clone()`.
     */
    CLP(const CLP &other, int thread);
#endif
    virtual ~CLP(); ///< The destructor.

//    void* operator new(size_t iSize) throw(); ///< Operator `new`.
//    void operator delete(void* pPtr) throw(); ///< Operator `delete`.
//	void* operator new[](size_t iSize) throw(); ///< Operator `new[]`.
//    void operator delete[](void* pPtr) throw(); ///< Operator `delete[]`.

/////////////////////////////////////////////////
//             Operations with matrix
/////////////////////////////////////////////////

    /**
     * Allocates memory for an LP problem of the required  size.
     * \param[in] m number of rows;
     * \param[in] n number of columns;
     * \param[in] nz number of non-zeroes in the matrix;
     * \param[in] bRowGen `true` if generation of rows is assumed;
     * \param[in] bColGen `true` if generation of columns is assumed;
     * \param[in] mMax maximum number of rows;
     * \param[in] nMax maximum number of columns;
     * \param[in] nzMax maximum number of non-zeroes in the matrix.
     * \remark If `bRowGen=true` and `mMax=0`, `openMatrix()` will estimate `mMax`;
     *   if `bColGen=true` and `nMax=0`, `openMatrix()` will estimate `nMax`.
     * \throws CMemoryException lack of memory.
     */
    virtual void openMatrix(int m, int n, int nz,
                     bool bRowGen=true, bool bColGen=false, 
                     int mMax=0, int nMax=0, int nzMax=0);

    /**
     * Does all preparations for solving the problem: allocates memory, does preprocessing, scales the matrix.
     */
    virtual void closeMatrix();

    /**
     * The function adds a constraint (empty row) to the matrix.
     * \param[in] hd handle of the constraint;
	 * \param[in] type type of the constraint;
	 * \param[in] lhs,rhs  respectively, left (LHS) and right (RHS) hand sides of constraint.
 	 * \return constraint index.
 	 * \throws CMemoryException lack of memory.
     * \sa addRow(), safeAddRow().
     */
    int addCtr(tagHANDLE hd, unsigned type, double lhs, double rhs);

    /**
     * This function is used to add a row (constraint) to the matrix.
     * \param[in] hd handle of row (constraint);
	 * \param[in] type type of constraint;
	 * \param[in] lhs,rhs respectively, left (LHS) and right (RHS) hand sides of constraint;
 	 * \param[in] sz number of non-zero entries;
	 * \param[in] dpVal, ipCol references to arrays of size `sz`, `dpVal[i]` is the coefficient at variable `ipCol[i]`;
	 * \param[in] bSort if `true`, row entries are sorted in increasing order of their column indices.
 	 * \return constraint index.
 	 * \throws CMemoryException lack of memory.
 	 * \sa addCtr(), safeAddRow(), getRow().
     */
    int addRow(tagHANDLE hd, unsigned type, double lhs, double rhs,int sz, double* dpVal, int* ipCol, bool bSort=true);

    /**
     * This function is a safer version of addRow().
     * If `dpVal` and  `ipCol` are references to internal MIPCL arrays such as `m_dpFb` or `m_ipArray`,
     * memory for such arrays may be reallocated during the call to `addRow()`,
     * and, as a consequence, the pointers `dpVal` and `ipCol` becomes not valid.
     * \param[in] hd,type,lhs,rhs,sz,dpVal,ipCol,bSort have the same meanings as parameters of `addRow()`.
 	 * \return constraint index.
     * \throws CMemoryException lack of memory.
     * \sa addCtr(), addRow(), getRow().
    */
    int safeAddRow(tagHANDLE hd, unsigned type, double lhs, double rhs, int sz, double* &dpVal, int* &ipCol, bool bSort=true);

    /**
     * The function adds a variable (empty column) to the matrix.
     * \param[in] hd handle of variable; if hd < 0, then you must overload
     *  `getColumn()` (version with 11 parameters);
     * \param[in] type type of variable;
     * \param[in] cost objective coefficient;
     * \param[in] l,u lower and upper bounds of variable;
     * \return index of added column (variable).
     * \throws CMemoryException lack of memory.
     */
    int addVar(tagHANDLE hd, unsigned type, double cost,double l, double u);

    /**
     * The function adds a column to the matrix.
     * \param[in] hd handle of the variable; if hd >= 0, then you must overload
     *  `getColumn()` (version with 11 parameters);
     * \param[in] type type of variable;
     * \param[in] cost objective coefficient;
     * \param[in] l,u lower and upper bounds of variable;
     * \param[in] sz number of nonzero entries;
     * \param[in] dpVal,ipRow arrays of size `sz`; `dpVal[i]` is coefficient in row `ipRow[i]`;
 	 * \param[in] bSort if `true`, column entries are sorted in increasing order of their row indices.
     * \return index of added column (variable).
     * \throws CMemoryException lack of memory.
     */
    int addColumn(tagHANDLE hd, unsigned type, double cost,double l, double u,
        int sz, double* dpVal, int* ipRow, bool bSort=true);

    /**
     * This function is a safer version of addColumn().
     * If `dpVal` and  `ipCol` are references to internal MIPCL arrays such as `m_dpFb` or `m_ipArray`,
     * memory for such arrays may be reallocated during the call to `addColumn()`,
     * and, as a consequence, the pointers `dpVal` and `ipRow` becomes not valid.
     * \param[in] hd,type,cost,l,u,sz,dpVal,ipRow,bSort have the same meanings as parameters of `addColumn()`.
     * \return index of added column (variable).
     * \throws CMemoryException lack of memory.
     * \sa addVar(), addColumn(), getColumn().
    */
    int safeAddColumn(tagHANDLE hd, unsigned type, double cost,double l, double u,
        int sz, double* &dpVal, int* &ipRow, bool bSort=true);

    /**
     * The functions adds to the matrix the entry of value `val` into row `i` and column `j`;
     * if `i < 0`, the function sets the cost of variable `j`, `m_dpC[j]=val`.
     * \param[in] val entry value;
     * \param[in] i row index;
     * \param[in] j column index.
     * \throws CMemoryException lack of memory.
     */
    void addEntry(double val, int i, int j);

    /**
     * The function changes the matrix entry in row `i` and column `j`;
     * if `i < 0`, the function changes the cost of variable `j`, setting `m_dpC[j]=val`.
     * \param[in] val new entry value;
     * \param[in] i row index;
     * \param[in] j column index.
     * \throws CMemoryException lack of memory.
     */
    void changeEntry(double val, int i, int j);
    
    /**
     * \return number of non-zero entries in row `i`.
     */
    int getRowSize(int i) const
    	{return m_ipRowSize[i];}

    /**
     * \return number of non-zero entries in column `j`.
     */
    int getColumnSize(int j) const
    	{return m_ipColSize[j];}

    /**
     * \param[in] i constraint index.
     * \return left hand side of constraint `i`.
     * \sa  getRHS().
     */
    double getLHS(int i) const;

    /**
     * \param[in] i constraint index.
     * \return right hand side of constraint `i`.
     * \sa  getLHS().
     */
    double getRHS(int i) const;

    /**
     * \param[in] j index of a variable.
     * \return lower bound of variable `j`.
     * \sa  getVarUpBound().
     */
    double getVarLoBound(int j) const;

    /**
     * \param[in] j index of a variable.
     * \return upper bound of variable `j`.
     * \sa  getVarLoBound().
     */
    double getVarUpBound(int j) const;

    /**
     * \param[in] j index of a variable.
     * \return cost of variable `j`.
     * \sa setObjCoeff().
     */
    double getObjCoeff(int j) const;

    /**
     * The function sets (or changes) the cost of a variable.
     * \param[in] j index of a variable.
     * \param[in] val cost of variable `j`.
     * \sa getObjCoeff(), addEntry(), changeEntry().
     */
    void setObjCoeff(int j, double val);

protected:
    /**
     * The function multiplies the constraint coefficients by a given factor.
     * \param[in] i constraint index;
     * \param[in] factor multiplier.
     */
    void multiplyCtr(int i, double factor);

    /**
     * The functions extends the type of variable `j`
     *  by  bitwise ORing its current type with the flags stored in the parameter `type`.
     * \param[in] j index of variable;
     * \param[in] type is bitwise OR of members of `enVarType` and `CMIP::enVarType`.
     * \sa enVarType, CMIP::enVarType.
     */
    virtual void extendVarType(int j, unsigned type)
        {m_ipVarType[j]|=type;}

    /**
     * The functions extends the type of constraint `i`
     *  by  bitwise ORing its current type with the flags stored in the parameter`type`.
     * \param[in] i index of constraint;
     * \param[in] type is bitwise OR of members of `enCtrType` and `CMIP::enCtrType`.
     * \sa enCtrType, CMIP::enCtrType.
     */
    virtual void extendCtrType(int i, unsigned type)
        {m_ipCtrType[i]|=type;}

   	/**
   	 * \param[in] j variable index.
   	 * \return `true` if variable `j` is integer-valued.
   	 */
    bool isVarIntegral(int j) const
    {
        return (m_ipVarType[j] & VAR_INT)? true: false;
    }

   	/**
   	 * \param[in] j variable index.
   	 * \return `true` if variable `j` is can be used for branching on it.
   	 * \sa CMIP::isVarUsedForBranching().
   	 */
    virtual bool isVarUsedForBranching(int j) const {return false;}

    /**
     * To be overloaded in `CMIP`.
     * \sa CMIP::isVarStrongIntegral().
     */
    virtual bool isVarStrongIntegral(int /*j*/) const
        {return false;}

    /**
     * To be overloaded in `CMIP`.
     * \sa CMIP::isVarBinary().
     */
    virtual bool isVarBinary(int /*j*/) const
        {return false;}

    /**
     * To be overloaded in `CMIP`.
     * \sa CMIP::isVarScalable().
     */
    virtual bool isVarScalable(int j) const
        {return (m_ipVarType[j] & VAR_NOT_MOD)? false: true;}

    /**
     * To be overloaded in `CMIP`.
     * \sa CMIP::isVarDeletable().
     */
    virtual bool isVarDeletable(int j) const
        {return (m_ipVarType[j] & VAR_NOT_MOD)? false: true;}

    /**
     * To be overloaded in `CMIP`.
     * \sa CMIP::isVarModifyable().
     */
    virtual bool isCtrModifyable(int /*i*/) const
        {return true;}

    /**
     * \param[in] row row index;
     * \param[out] dpVal,ipCol arrays of size `k`, where `k` is return value;
     * `dpVal[i]` is coefficient in column `ipCol[i]`;
     * \param[in] bScaled if `true`, row of scaled matrix is returned; otherwise, row of non-scaled matrix.
     * \return number of non-zero coefficients in row `row`.
     */
    int getRow(int row, double* dpVal, int* ipCol, bool bScaled=false);

    /**
     * \param[in] col column index;
     * \param[out] dpVal,ipRow arrays of size `k`, where `k` is return value;
     * `dpVal[i]` is the coefficient in row `ipRow[i]`;
     * \param[in]  bScaled if `true`, column of scaled matrix is returned; otherwise, column of non-scaled matrix.
     * \return number of non-zero coefficients in column `col`.
     */
    int getColumn(int col, double* dpVal, int* ipRow, bool bScaled=false);

    /**
     * \param[in] recomputeBasicVarsCost if `true`,
     * total cost of non-basic variables stored in `m_dBasicVarCost` is recomputed.
     * \return cost of of current basic solution.
     */
    double computeObjValue(bool recomputeBasicVarsCost=false);

////////////////////////////////////////
//       Preprocessing
///////////////////////////////////////
protected:
    /**
     * \param[in] a,b two integers (although written as double reals).
     * \return global common divisor of input integers.
     */
    static double GCD(double a, double b);

////////
private:
    /**
     * The procedure preprocess a given constraint.
     *
     * \param[in] ind constraint index;
     * \param[in,out] last index of last non-processed constraint in the circular list of non-processed constraints;
     * \param[in,out] lhs, rhs left and right hand sides of constraint to be processed;
     * \param[in,out] loBd,upBd minimum and maximum constraint value.
     * \return `false` if it was detected that the constraint being processed is inconsistent; otherwise, the return value is `true`.
     */
    bool preprocessCtr(int ind, int& last,double& lhs, double& rhs, double& loBd, double& upBd);

    /**
     *  The procedure deletes redundant (with INFINITY bounds) and singleton rows.
     *
     *  \return number of deleted rows.
     */
    int deleteRedundantRows();

    /**
     * If a free variable belongs to an equation, its deleted from the matrix.
     *
     * \return number of deleted free variables.
     */
    int deleteFreeVars();

    /**
     * Sometimes, we can fix a free variable if it admits a value for which
     * all the constraints involving this free variable become free,
     * i.e., they are valid for all feasible values of the other variables.
     *
     * \param[in] j index of free variable;
     */
    void fixFreeVariable(int j);

    void deleteFixedVar(int col); ///< deletes fixed variable indexed by `col`.
    int deleteFixedVars(); ///< iteratively calls `deleteFixedVar()`.

    /**
     * The function is called from `delFixedVars()`.
     * \return  `true` if variable `col` has been fixed by duality reasoning.
     */
    bool dualFixVar(int col);

    int processParallelRows(); ///< combines any pairs of parallel constraints into one two-sided constraint.

    /**
     * The procedure compares two rows: a row is less than another one if the value of the first row is less than the value of the second.
     * \param[in] r1,r2 row indices.
     * \return number of fixed variables.
     * \remark Each row `r` is associated with a flag `cpRowPassiv[r]` non-zero value of which indicates that
     *  row `r` is dominated (less or equal) by another row. So, if `compareRows()` establishes that
     *  row `r1` is less or equal than row `r2` or vice versa, it assigns the value of `1` to either `cpRowPassiv[r]` or `cpRowPassiv[r2]'.
     */
    int compareRows(int r1, int r2);

    /**
     *  The procedure eliminates all constraints that dominate (or are dominated by) a given constrain;
     *  \param[in] row constraint index.
     *  \return number of fixed variables plus number of removed rows.
     */
    int removeDominatedRows(int row);

    /**
     *  The procedure eliminates all dominated constraints; does its job by calling `removeDominatedRows()`.
     *  \param[in] timeLimit time limit.
     *  \return number of fixed variables plus number of removed rows.
     */
    int processDominatedRows(__LONG timeLimit);

    /**
     *  The procedure deletes a given constraint with single variable, and simultaneously updates the bounds on that variable.
     *  \param[in] row constraint index.
     */
    void deleteSingletonRow(int row);

    /**
     * The function, in the equality indexed by `row`, seeks for a variable to be substituted for
     * the linear expression induced from equality `row` so that the number of fill-ins (additional non-zero entries in the matrix)
     * does not exceeds some predefined value.
     * \param[in] row index of some equality constraint;
     * \param[in] flag if `flag=0`, only real variables are considered; if `flag=1`, only integer variables are considered,
     * and if `flag=2`, only binary variables are considered.
     * \returns index of variable to be substituted, or `-1` if no variable is good for substituting.
     */
    int estimateFillIns(int row, int flag);

    /**
     * Only small size equalities are used when substituting variables.
     * \param[in] row - row index of equality constraint `\sum_{i\in I} a_{row,i} x_i = b_{row}`
     * \return `-1` if no variable has been chosen to be substituted;
     *   otherwise, variable `j` to be substituted for expression
     *       `b_{row} - sum(i in I\ j) a_{row,i} x_i`.
     */
    int getBestFillInColumn(int row);

    /**
     * The procedure substitute variable `x_j` for the expression
     *         `b_{row} - sum(i in I\ j} a_{row,j} x_j`.
     * \param[in] row - row index of equality constraint `\sum_{i in I} a_{row,i} x_i = b_{row}`;
     * \param[in] j  index from `I`.
     */
    void substitute(int row, int j);

    /**
     * The procedure rearranges the matrix so that
     * the non-zeroes of any line (column or row) `i` are followed by the non-zeroes of line `i+1`.
     *
     * \param[in] align if `align==ALIGN_COLUMN`, lines are columns; otherwise, lines are rows.
     */
    void alignMatrix(enAlign align);

    /**
     * After the matrix has been closed, any new constraint
     * `b1 <= sum(i=0,...,sz-1) dpVal[i]*x[ipCol[i]] <= b2`,
     * before adding to the matrix, must be preprocessed.
     * \param[in] n number of variables in initial (non-preprocessed) matrix;
     * \param[in,out] b1,b2 left and right hand sides;
     * \param[in,out] sz number of non-zero entries (size);
     * \param[in,out] dpVal array of coefficients (of size `sz`);
     * \param[in,out] ipCol array of column indices (of size `sz`).
     */
    void preprocNewCtr(int n, double& b1, double& b2, int& sz, double* dpVal, int* ipCol);

    /**
     * The adds all constrains containing variable `x[col]` to the list of _active_ constraints.
     *
     * \param[in] col column index;
     * \param[in,out] last position in list of active constraints.
     */
    void labelActiveCtrs(int col, int& last);

    /**
     *  The procedure applies standard preprocessing techniques to a given list of constraints.
     *  \param[in] last number of constraint indices written in array `ipCtrLst=m_ipArray+m_iN`.
     *  \return `true` if no contradictions have been detected; otherwise, `false`.
     *
     */
    bool preprocessMatrix(int last);

    /**
     * The procedure does some preparations for calling `preprocessMatrix(int&)`.
     */
    bool preprocessMatrix();

    /**
     * The procedure compares pairs of rows to remove from the formulation
     * dominated (or equivalent) inequalities.
     *
     * \param[in] maxInclRowSize only rows of size less than `maxInclRowSize` are processed.
     * \return number of deleted rows.
     */
    int seekIncludedCtrs(int maxInclRowSize);

protected:
    /**
     * The procedures replaces bounds `[d1[j],d2[j]]` for `[0,d2[j]-d1[j]]`) for those variables, `j`, which lower bounds `d1[j] <!= 0`.
     * \return number of changed bounds.
     */
    int shiftBounds();

    /**
     * Just a skeleton. To be overloaded in `CMIP`.
     * \return `true`.
     * \remark Default implementations does nothing.
     */
    virtual bool preprocessInit() {return true;}

    /**
     * Just a skeleton. To be overloaded in `CMIP`.
     * \return `true`.
     * \remark Default implementations does nothing.
     */
    virtual bool preprocessPlus() {return true;}

    /**
     * Preprocessing without probing.
     * \param[in] bDominant if `true`, dominated constraints are detected and processed.
     * \return total number of variables fixed and constraints eliminated.
     */
    int basicPreprocess(bool bDominant=true);

    /**
     * The procedure reformulate an LP (or MIP) by calling, in turn,
     *      - `preprocessInit()`,
     *      - `basicPreprocess()`,
     *      - `preprocessPlus()`.
     * \return `false` if inconsistency has been detected; otherwise, `true`.
     */
    bool preprocess();

////////////////////////////////////////
// S C A L I N G
private:
    void SCL_ShiftScaleFactors(); ///< The procedure is used only in `SCL_scaleMatrix()`.

    /**
     * \param[out] minExp minimum exponent among matrix entries.
     *  \return the maximum absolute value among exponents of matrix entries.
     */
    int getMaxEntryExponent(int &minExp);

    /**
     * The procedure compute a scaling factor for a given row.
     *
     * \param[in] row index of row to be scaled;
     * \param[in] GMflag if `true`, geometric-mean rule be used; otherwise, arithmetic-mean rule be used.
     * \return exponent of the scale factor.
     */
    int SCL_HScaleRow(int row, bool GMflag);

    /**
     * The procedure compute a scaling factor for a given column.
     *
     * \param[in] col index of column to be scaled;
     * \param[in] GMflag if `true`, geometric-mean rule be used; otherwise, arithmetic-mean rule be used.
     * \return exponent of the scale factor.
     */
    int SCL_HScaleColumn(int col, bool GMflag);

    /**
     * The procedure scales all matrix rows.
     *
     * \param[in] GMflag if `true`, geometric-mean rule be used; otherwise, arithmetic-mean rule be used.
     */
    void SCL_HScaleRows(bool GMflag);

    /**
     * The procedure scales all matrix columns.
     *
     * \param[in] GMflag if `true`, geometric-mean rule be used; otherwise, arithmetic-mean rule be used.
     */
    void SCL_HScaleColumns(bool GMflag);

    /**
     * The procedure iteratively scales rows an columns processing each of them `roundNum` times.
     * \param[in] roundNum number of iterates (rounds);
     * \param[in] bGMflag  if `true`, geometric mean scaling is applied; otherwise, min-max scaling.
     * \param[in] bRowAlign if true, rows are scaled after columns.
     */
    void SCL_HScale(int roundNum, bool bGMflag, bool bRowAlign);

    /**
     * The procedure is used internally by `idealScaling()`,
     * which implements a scaling algorithm that iteratively solves
     * a number of the mean-cost cycle problems defined on subgraphs
     * of the bipartite graph induced by the constraint matrix.
     *
     *
     * \param[in] dMaxExp maximum scaling exponent;
     * \param[in] iNodeNum,ipNodeSep,ipNode,ipRowNode,ipColNode parameters that describe underlying graph at current iterate;
     * \param[out] ipNewNodeSep,ipNewNode,ipNewRowNode,ipNewColNode parameters that describe underlying graph for next iterate;
     * \param[in] ipNum,ipLevel,ipLowLink,ipStack,ipParent,ipCurRow,ipCurCol,ipCurNode these arrays are used to implement depth-search algorithm
     * that computes node levels;
     * \param[in,out] dpRowScale, dpColScale arrays that store row and column scaling exponents.
     * \return `true` if the procedure has decreased the maximum scaling exponent.
     */
    bool SCL_computeLevels(double dMaxExp, int& iNodeNum,
        int* &ipNodeSep,int* &ipNode, int* &ipRowNode, int* &ipColNode, 
        int* &ipNewNodeSep,int* &ipNewNode, int* &ipNewRowNode, int* &ipNewColNode, 
        int* ipNum, int* ipLevel, int* ipLowLink, int* ipStack, int* ipParent,
        int *ipCurRow, int* ipCurCol, int* ipCurNode,
        double *dpRowScale, double *dpColScale);

    /**
     * The procedure implements the "_ideal scaling_", when it is not possible to
     * improve scaling of any line (row or column) without deterioration of scalings the other lines.
     *
     * \param[in] maxExp maximum exponent among matrix entries.
     */
    void idealScaling(int maxExp);

    /**
     * The procedure is used internally in `minMaxScaling()`.
     * \remark set `flag` to `true` only for the last call to the procedure.
     */
    bool SCL_estimateExponent(int q, int *ipS, int* ipT, bool flag=false);

    /**
     * The procedure implements min-max scaling, i.e.,
     * in each line (row or column) the mean value of the minimum and maximum absolute coefficient values are computed,
     * and then the line is multiplied by the power of 2 nearest to that mean value.

     * \param[in] maxExp maximum exponent among matrix entries.
    */
    void minMaxScaling(int maxExp);

protected:
    void scaleObj(); ///<  Scales objective with previously computed column scale factors.
    void SCL_scaleMatrix(); ///< Scales matrix with previously computed scale factors.

    /**
     * \remark The method used for computing factors is determined by the current value of `m_eScaling`,
     *  which is set by `setScaling()`.
     */
    void scaleMatrix(); ///< Compute factors and then scales matrix.
    void unscaleMatrix(); ///< Recovers matrix to initial state.

/////////////////////////////////////////////////////////////////////

    /**
     * The function deletes variable (column) `j` from the matrix and adds its expression
     *       x(j) = b +sum(i=0,...,sz-1} dpVal[i]*x(ipCol[i])
     * to the preprocessing stack.
     * \param[in] j index of variable;
     * \param[in] hd handle of deleted variable, or `hd in {-1,-2}`;
     *   if `hd = -1`, then `sz = 1` and variable `x[ipCol[0]]` is substituted for `x[ipCol[0]]*dpVal[0]+b`;
     *   if `hd = 2`, then `sz = 1`, and, for `c=reinterpret_cast<int*>(dpVal)[0]`,
     *   if `x[ipCol[0]] >= 0.0`, then `x[c]=0.0`; otherwise, `x[c]=-x[ipCol[0]]` and `x[ipCol[0]]=0.0`;
     * \param[in] b free term of expression;
     * \param[in] sz number of variables in expression
     * \param[in] ipCol list of `sz` variables (columns) in expression;
     * \param[in] dpVal list of `sz` coefficients in expression.
     */
    virtual void deleteVariable(int j, int hd,  double b, int sz, int *ipCol, double *dpVal);

    /**
     * The procedure deletes from the matrix all the rows and columns that were deleted during preprocessing.
     * This procedure is also used in `CMIP` to delete from the matrix those cuts that are not "tight" for optimal node LP-solutions.
     * \param[in] m0 initial row to start search for unneeded rows;
     * \param[in] n0 initial column to start search for unneeded columns;
     * \param[in] prepFlag if `true`, expressions in preprocessing stack are also updated;
     * \param[in] bEntry is used internally in `CMIP`.
     * \attention Never set both, `prepFlag` and `bEntry`, flags to `false` in your applications.
     */
    void compressMatrix(int m0=0, int n0=0, bool prepFlag=true, bool bEntry=true);

    /**
     * The function is called when a set of constraints is deleted from the matrix.
     * Default `CLP`-implementation does nothing.
     * \sa CMIP::setCtrsInactive().
     */
    virtual void setCtrsInactive(int /*sz*/, int* /*ipCtr*/) {};

    /**
     * The function is called when a set of columns is deleted from the matrix.
     * Default `CLP`-implementation does nothing.
     * \sa CMIP::setColumnsInactive().
     */
    virtual void setColumnsInactive(int /*sz*/, int* /*ipCtr*/) {};

    /**
     * The procedure deletes all "_non-tight_" inequalities and columns previously generated during the solution process,
     * and which indices `>=` m for rows, and `>= n` for columns.
     * \param[in] m only rows which index is not less than `m` are considered as candidates for deletion;
     * \param[in] n only columns which index is not less than `n` are considered as candidates for deletion;
     * \param[in] bFull if `true`, matrix must be compressed and solution be updated;
     * \param[in] tight inequality (resp, column) is "_tight_" if its _slack_ (resp., _reduced cost_)
     *  is less than value of `tight`.
     */
    void deleteNonBasicLines(int m, int n, bool bFull, double tight);

    /**
     * The function is used to express constraints in scaled variables.
     * \param[in] sz size (number of variables) of constraint;
     * \param[in,out] dpVal,ipCol arrays of size `sz`; in constraint being scaled,
     *  `dpVal[i]` is coefficient at variable `ipCol[i]`.
     *  \sa scaleCtr().
     */
	void scaleRow(int sz, double* dpVal, int* ipCol);

	/**
	 * When implementing a separation, or cut generation procedure,
	 *  the user knows nothing about transformations of the matrix
	 * that will be done by the solver. Therefore, in user defined procedures inequalities (cuts) are generated
	 * assuming that the matrix had not been changed.
	 *  Before adding such inequalities to the matrix, they must be previously scaled.
	 *
     * Given a constraint expressed in scaled variables,
     *  the function computes _scaling factor_ and multiplies that constraint by this factor.
     * \param[in,out] lhs,rhs left and right sides of constraint being scaled;
     * \param[in] sz size (number of variables) of constraint;
     * \param[in,out] dpVal,ipCol arrays of size `sz`; in constraint being scaled,
     *  `dpVal[i]` is coefficient at variable `ipCol[i]`.
     * \return scaling factor.
     *  \sa scaleRow().
     */
    int scaleCtr(double& lhs, double& rhs, int sz, double* dpVal, int* ipCol);

    /**
     * The function multiplies each coefficient (including the cost) in the column
     * by the factor of that row the coefficient belongs to.
     * \param[in,out] cost cost of column being scaled;
     * \param[in] sz size (number of rows with non-zero coefficients) of column;
     * \param[in] dpVal,ipCol arrays of size `sz`; in column being scaled,
     *  `dpVal[i]` is coefficient in row `ipRow[i]`.
     *  \sa scaleColumn().
     */
	void scaleColumn(double &cost, int sz, double* dpVal, int* ipCol);

	/**
	 * When implementing a column generating procedure, the user knows nothing about transformations of the matrix
	 * that will be done by the solver. Therefore, in user defined procedures columns are generated assuming
	 * that the matrix had not been changed. Before adding such columns to the matrix, they must be previously scaled.
	 *
     * Given a variable (cost coefficient, column of the matrix, and lower and upper bounds),
     * `scaleVar()` computes _scaling factor_ and multiplies that column by this factor.
     * \param[in,out] cost cost of column being scaled;
     * \param[in,out] l,u lower and upper bounds that corresponds to column being scaled;
     * \param[in] sz size (number of variables) of constraint;
     * \param[in,out] dpVal,ipRow arrays of size `sz`; in column being scaled,
     *  `dpVal[i]` is coefficient in row `ipRow[i]`.
     * \return scaling factor.
     *  \sa scaleColumn().
     */
    int scaleVar(double &cost, double& l, double& u, int sz, double* dpVal, int* ipRow);

//////////////////////////////////////////////////////////////////////////
//                   Initialization/memory allocation                   //
//////////////////////////////////////////////////////////////////////////
    /**
     * The procedure allocates memory needed for both, prime and dual, simplex algorithms.
     * \throws CMemoryException lack of memory
     */
    void allocMemForSimplex();

    /**
     * The function is used internally in `CMIP`.
     * It builds linked lists of both, rows and columns.
     */
    void buildRowColLists();

private:
    void allocColMemForSimplex(); ///< The procedure allocates memory for "column" arrays used in simplex procedures.
    void allocRowMemForSimplex(); ///< The procedure allocates memory for "row" arrays used in simplex procedures.

    /**
     * The procedure reallocates memory for "column" arrays used in simplex procedures.
     *
     * If the number of columns (`m_iN`) approaches its limit (`m_iNmax`),
     * this limit is increased 1/3-d of its initial value, and then memory is reallocated.
     *
     */
    void reallocColMemForSimplex();

    /**
     * The procedure reallocates memory for "row" arrays used in simplex procedures.
     *
     * If the number of rows (`m_iM`) approaches its limit (`m_iMmax`),
     * this limit is increased 1/3-d of its initial value, and then memory is reallocated.
     *
     */
    void reallocRowMemForSimplex();

    void allocMemForNorms(); ///< The procedure allocates memory for arrays storing row or column norms.
    void reallocMemForNorms(); ///< The procedure reallocates memory for arrays storing row or column norms.

    void allocMemForEntries(); ///< The procedure allocates memory for the matrix.
protected:
    /**
     * The function reallocates memory for the matrix.
     *
     * If the number of nonzero entries (`m_iNZ`) approaches its limit (`m_iNZmax`),
     * this limit is increased 1/3-d of its initial value, and then memory is reallocated.
     *
     * \param[in] nz maximum number of entries.
     * \remark This function is overloaded in CMIP to make it thread safe.
     * \attention The procedure should not be called in user applications.
     */
    virtual void reallocMemForEntries(int nz);

    /**
     * The function is called when the number of rows (constraints) exceeds `m_iMmax`.
     * \attention `incMaxRowNumber()` should not be called in user applications.
     * \remark This function is overloaded in CMIP to make it thread safe.
     */
    virtual void incMaxRowNumber();

    /**
     * The function is called when the number of columns (variables) exceeds `m_iNmax`.
     * \attention `incMaxColumnNumber()` should not be called in user applications.
     * \remark This function is overloaded in CMIP to make it thread safe.
     */
    virtual void incMaxColumnNumber();

    /**
     * The function allocates memory for three auxiliary arrays,
     *  `m_dpArray` and `m_ipArray`,  both of size `m_iNmax+m_iMmax`,
     *  `m_dpArray` and `m_ipArray`,  both of size `m_iNmax+m_iMmax`,
     *  `dpW` of size `2*m_iMmax`, `m_dpFd` of size 2*m_iNmax`, and `m_dpFb` of size `max{2*m_iMmax,m_iNmax}`.
     *
     * \attention `allocMemForAuxArrays()` should not be called in user applications.
     * \sa reallocMemForAuxArrays().
     */
    virtual void allocMemForAuxArrays();

    /**
     * The function reallocates memory for five auxiliary arrays,
     *  `m_dpArray` and `m_ipArray`,  both of size `m_iNmax+m_iMmax`,
     *  `dpW` of size `2*m_iMmax`, `m_dpFd` of size 2*m_iNmax`, and `m_dpFb` of size `max{2*m_iMmax,m_iNmax}`.
     *  This function is called when either the number of columns
     *  exceeds `m_iNmax`, or the number of rows exceeds `m_iMmax`.
     *
     * \param[in] rowMem if `true`, memory for arrays which size depends on `m_iMmax` is reallocated;
     * \param[in] colMem if `true`, memory for arrays which size depends on `m_iNmax` is reallocated.
     * \attention `reallocMemForAuxArrays()` should not be called in user applications.
     * \sa allocMemForAuxArrays().
     */
    void reallocMemForAuxArrays(bool rowMem, bool colMem);

private:
    void allocMemForProblem(); ///<  \throws CMemoryException
    void allocMemForColumns(); ///<  \throws CMemoryException
    void reallocMemForColumns(); ///<  \throws CMemoryException
    void allocMemForRows(); ///<  \throws CMemoryException
    void reallocMemForRows(); ///<  \throws CMemoryException

    /**
     * For small LPs, which number of rows and number of columns do not exceed `MAX_SHORT_SIZE`,
     * the procedure stores in the given memory buffer the current basis.
     *
     *\param[in] mem memory buffer of sufficient size.
     * \sa saveIntBasis(), saveBasisBitwise().
     */
    int saveShortBasis(int* mem);

    /**
     *  For LPs of medium size, the procedure stores in the given memory buffer the current basis.
     *
     * \param[in] mem memory buffer of sufficient size.
     * \sa saveShortBasis(), saveIntBasis().
     */
    int saveIntBasis(int* mem);

    /**
     * For LPs of medium size, the procedure stores in the given memory buffer the current basis.
     *
     * \param[in] mem memory buffer of sufficient size.
     * \sa saveShortBasis(), saveBasisBitwise().
     */
    int saveBasisBitwise(unsigned* mem);
///
    /**
     * Restores the basis previously stored by `saveShortBasis()`.
     *
     * \param[in] mem memory buffer that stores basis.
     */
    void restoreShortBasis(int* mem);

    /**
     * Restores the basis previously stored by `saveIntBasis()`.
     *
     * \param[in] mem memory buffer that stores basis.
     */
    void restoreIntBasis(int* mem);

    /**
     * Restores the basis previously stored by `saveBasisBitwise()`.
     *
     * \param[in] mem memory buffer that stores basis.
     */
    void restoreBasisBitwise(unsigned* mem);
protected:
    /**
     * The function packs the basis into a memory buffer (pointed to by `mem`).
     * \param[out] mem pointer to memory array of size at least `m_iN`.
     * \return size (in integer words) of stored data.
     * \sa restoreBasis().
     */
    int saveBasis(int* mem);

    /**
     * The function restores basis previously stored by `saveBasis()`.
     * \param[in] mem pointer to memory array with stored basis.
     * \sa saveBasis().
     */
    void restoreBasis(int* mem);

    /**
     * The function restores full basis from its short version given by two arrays `m_ipRowMap` and `m_ipColMap`.
     * \param[in] ipRowMap,ipColMap pointer to memory array with stored basis  in short form;
     *  usually, `ipRowMap` and `ipColMap` are copies of `m_ipRowMap` and `m_ipColMap`.
     *  \remark The function is used internally in `CMIP`.
     * \sa `m_ipRowMap` and `m_ipColMap`.
     */
    void restoreBasis(int* ipRowMap, int* ipColMap);
/////////////////////////////////////////////////////////////

// SEPARATION

    /**
     * \param[in] i constraint index;
     * \param[in] scaled if `true`, slack for scaled matrix is returned; otherwise, for not-scaled.
     * \return slack value of right part of constraint `i` for basic solution currently stored in memory.
     */
	double getCtrRightSlack(const int i, bool scaled=true) const;

    /**
     * \param[in] i constraint index;
     * \param[in] scaled if `true`, slack for scaled matrix is returned; otherwise, for not-scaled.
     * \return slack value of left part of constraint `i` for basic solution currently stored in memory.
     */
	double getCtrLeftSlack(const int i, bool scaled=true) const;

private:
	/**
	 * Multiplies row `i` by vector stored in `X`, assuming that components in `X` are listed so as in `m_ipBasicColumn`.
	 * \param[in] i row index;
	 * \param[in] X array of size `m_iN`.
	 * \return   sum(a(j,m_iBasicColumn[j])*X[j] for j=0,...,sz-1).
	 */
    double getBasicLeftHand(int i, double* X) const;

    /**
     * The procedure verifies whether the current basis satisfies or violates a given inequality.
     *
     * \param[in] i constraint index;
     * \param[out] beta violation value;
     * \param[out] side if `true`, right-hand side inequality is violated;
     * otherwise, left-hand side inequality is violated.
     * \return `true` if given inequality is not violated, and `false` otherwise.
     */
    bool checkInEq(int i, double& beta, bool& side) const;

    /**
     * The procedure looks for an equation that is violated by the current basic solution.
     *
     * \param[out] beta violation value;
     * \param[out] side if `true`, right-hand side inequality is violated;
     * otherwise, left-hand side inequality is violated.
     * \return index of violated inequality.
     */
    int getViolatedEquation(double& beta, bool& side) const;

    /**
     * The procedure changes the separation rule.
     * \param[in] eSepRule new separation rule;
     * \param[in] bUpdateFlag if `true`, current basis must be updated.
     */
    void switchSeparation(enSepRule eSepRule, bool bUpdateFlag=false);

    /**
     * This procedure is to find a violated constrain according to the the steepest edge rule.
     *
     * \param[out] r row index of violated inequality;
     * \param[out] side if `true`, right hand side inequality is violated, otherwise, left hand side inequality is violated;
     * \param[in,out] edg only those inequalities are processed which `edge` is less than input `edg` vaqlue.
     * \return number of processed inequalities.
     */
    int rowSeparate(int& r, bool& side, double& edg);

    /**
     * The procedure separates the current basic solution stored in `m_dpX`.
     * \param[out] side `true` if left side of inequality `r`, where `r` is return value, is violated; otherwise, `false`.
     * \return index of violated inequality.
     */
    int innerSeparate(bool& side);

/////////////////////////////////////////////////////////////
// PRICING
    /**
     * The function sets new pricing rule.
     *
     * \param[in] ePricingRule new pricing rule;
     * \param[in] bUpdateFlag if `true`, LU-partition of basic matrix is updated.
     * \sa enPricingRule.
     */
    void switchPricing(enPricingRule ePricingRule, bool bUpdateFlag=false);

    /**
     * This procedure is called in `pricing()` to detect a column to enter the basis.
     *
     * \param[out] r if non-negative, column `m_ipBasicColumn[r]` is to enter basis;
     * \param[in] edg only columns which norm is less than `edg` are considered.
     * \return 0 if no required column has been detected (`r < 0`),
     *  and non-zero otherwise (`r >= 0`).
     */
    int pricingColumns(int& r, double edg);

    /**
     * The procedure performs the pricing operation of the primal simplex method.
     * It detects a row (constraint) that leaves the basis, or a column (variable) that enters the basis.
     *
     * \return non-negative integer `t`; if `t < m_iBasisSize`, row `m_ipBasicRow[t]` leaves the basis,
     * and if `t >= m_iBasisSize`, column `m_ipBasicColumn[t]` enters the basis.
     */
    int pricing();
//////////////////////////////////////////////////////////////////////////
//                          Pivoting Functions                          //
//////////////////////////////////////////////////////////////////////////
public:
    /**
     * The functions sets the pricing rule to be used for solving LPs and root LPs of MIPs.
     * \param[in] pricingRule new pricing rule.
     * \sa enPricingRule, getPricingRule().
     */
    void setLpPricingRule(enPricingRule pricingRule)
    {
        m_eLpPricingRule=pricingRule;
    }

    /**
     * \return currently used pricing rule.
     */
    enPricingRule getPricingRule() const
    {return m_ePricingRule;}

    /**
     * The function sets the __separaion rule__ applied when solving LPs,
     *  and is also used when solving root LPs of MIPs.
     *  \param[in] sepRule new separation rule.
     * \sa enSepRule, getSepRule(), setMipSepRule().
     */
    void setLpSepRule(enSepRule sepRule)
    {
    	if (sepRule != SEP_ONLY_EQUATIONS)
    		m_eLpSepRule=sepRule;
    }

    /**
     * \return currently used _separation rule_.
     */
    enSepRule getSepRule() const
    {return m_eSepRule;}

    /**
     * At each iterate both, prime and dual, simplex algorithms solves at least two linear systems, `B x = b` and `B^Ty=c`.
     * Here $B$ denotes basic matrix, vector $b$ is composed from the basic components of the right or left hand side vectors,
     * and vector $c$ is composed from the basic components of the upper or lower bound vectors.
     * To solve the systems `B x = b` and `B^Ty=c`, `CLP` _factors_ the basic matrix, i.e., computes an LU-partition: `B=LU`,
     * where `L` is lower, and `U` upper triangular matrices. Since computing LU-partition is time-consuming operation,
     * it is computed only in those simplex-iterations that are multiple of some predefined number (say, 50).
     * In all other iterations, `CLP` updates the partition of the previous iteration.
     * Although updating LU-partitions is less costly than re-factoring the basic matrix,
     * such updates introduce additional errors in the solutions. To eliminate these the solver has to re-factor the basis
     * after carrying out a predefined number of updates, or in case the solver detects computational instability.
     *
     * `setMaxLuUpdateNum()` sets a new limit on the number of consecutive updates of the LU-partition.
     * \param[in] maxUpdateNum maximum number of LU-updates.
     */
    void setMaxLuUpdateNum(int maxUpdateNum)
        {m_iMaxLuUpdateNum=maxUpdateNum;}

private:
//         Steepest Edge
//////////////////////////////////////////////////////////////////////////
    /**
     * The function verifies whether a given constraint is in the reference framework.
     *
     * \param[in] i constraint index.
     * \return `true` if constraint `i` is in the reference framework.
     */
    bool isCtrRef(int i) const
    {return (m_ipCtrType[i] & CTR_REF)? true: false;}

    /**
     * The function verifies whether a given variable is in the reference framework.
     *
     * \param[in] i variable index.
     * \return `true` if variable `x[i]` is in the reference framework.
     */
    bool isVarRef(int i) const ///<
    {return (m_ipVarType[i] & VAR_REF)? true: false;}

    /**
     * The function adds a given constraint to the reference framework.
     *
     * \param[in] i constraint index.
     */
    void ctr2Ref(int i)
    {m_ipCtrType[i] |= CTR_REF;}

    /**
     * The function adds a given variable to the reference framework.
     *
     * \param[in] i index of variable.
     */
    void var2Ref(int i)
    {m_ipVarType[i] |= VAR_REF;}

    /**
     * The function removes a given constraint to the reference framework.
     *
     * \param[in] i constraint index.
     */
    void clearCtrRef(int i)
    {m_ipCtrType[i] &= ~CTR_REF;}

    /**
     * The function removes a given variable to the reference framework.
     *
     * \param[in] i index of variable.
     */
    void clearVarRef(int i)
    {m_ipVarType[i] &= ~VAR_REF;}

    /**
     *  The procedure reinitializes the row norms and the reference framework.
     */
    void DUAL_rectify();

    /**
     * The procedure computes the norm of a given column expressed in the current basis.
     *
     * \param[in] i norm of column `m_ipBasicColumn[i]` be computed.
     * \return norm value.
     */
    double computeColumnNorm(int i);

    /**
     *  The procedure computes norms for columns `m_ipBasicColumn[k1},...,m_ipBasicColumn[k2-1]`.
     *
     *  \param[in] k1,k2 indices of the first and last columns which norms to be computed.
     */
    void computeColumnNorms(int k1=0, int k2=0);

    /**
     * The procedure updates all column norms after the basis has been changed.
     *
     * \param[in] t index of row that left the basis.
     */
    void updateColumnNorms(int t);

    /**
     * The function computes the norm of a given row.
     *
     * \param[in] i row index.
     * \return norm value.
     */
    double computeRowNorm(int i);

    /**
     * The procedure computes the norms of all rows.
     */
    void computeRowNorms();

    /**
     * The procedure updates all row norms after the basis has been changed.
     *
     * \param[in] s `m_ipBasicRow[s]` is index of row that entered the basis.
     * \param[in] ur pivot value;
     * \param[in] tau norm of row that entered the basis.
     *
     */
    void updateRowNorms(int s, double ur, double tau);

//////////////////////////////////////////////////////////////////////////
    /**
     * The function is called only in `primeRatioTest()` to compute `Ub`.
     *
     * \param[in] t if `t < bs`, basic row `m_ipBasicRow[t]` is leaving basis; otherwise, column `m_ipBasicColumn[t]` enters basis.
     * \return number of non-zero entries in the basic part of `Ub`.
     */
    int computeUb(int t); ///< is called only in `primeRatioTest()`

    /**
     *  The procedure is called only in `primeRatioTest()` to compute the dual direction stored in `m_dpUd`,
     *  as well as the updates for column norms stored in `m_dpArray`.
     *
     *  \param[in] s is `s > 0`, row `s-1` enters basis; otherwise, column `-s-1` leaves basis;
     *  \param[in] t if `t < m_iBasisSize`, row `m_ipBasicRow[t]` leaves basis; otherwise, column `m_ipBasicColumn[t]` enters basis;
     *  \return `false` if column norms must be recomputed; otherwise, `true`.
     */
    bool computeUdAndColNormUpdates(int s, int t);

    /**
     * The procedure is called when the primal simplex method detects a basic variable
     *  which current value violates one of the bounds (left or right) of this variable.
     *
     *  \param[in] j index of variable;
     *  \param[in] side if `true`, right bound is violated; otherwise, left bound is violated;
     *  \param[in] delta violation value.
     */
    void makePrimeVarFeasible(int j, bool side, double delta); ///< is called when variable `j` is infeasible

    /**
     * The procedure is called when the primal simplex method detects a constraint
     * that are violated by the current basic solution.
     *
     *  \param[in] j constraint index;
     *  \param[in] side if `true`, right-hand side is violated; otherwise, left-hand side is violated;
     *  \param[in] delta violation value.
     */
    void makePrimeCtrFeasible(int j, bool side, double delta);

    /**
     * The procedure performs the ratio test of the primal simplex method.
     *
     * \param[in] t row that leaves basis;
     * \param[out] side if `true`, right-hand inequality of row entering basis is tight;
     *   otherwise, left-hand inequality is tight;
     * \param[out] bDegen if `true`, current simplex iteration is degenerate;
     * \param[in] bSafeMode if `true`, this function was called after some numeric instability had been detected.
     * \return row index of row to enter the basis.
     */
    int primeRatioTest(int t, bool& side, bool& bDegen, bool bSafeMode);

    /**
     * The procedure performs first the primal pricing to pick up a row that will leave the basis,
     * and then the primal ratio test to detect a row that enters the basis.
     *
     * \param[out] s if `s < m_iBasisSize`, row `m_iBasicRow[s]` enters the basis;
     *  if `s >= m_iBasisSize`, variable `x[m_ipBasicColumn[s]]` leaves the basis.
     * \param[out] t if `t < m_iBasisSize`, row `m_ipBasicRow[t]` leaves the basis;
     *  if `t >= m_iBasisSize`, variable `x[m_ipBasicColumn[t]]` enters the basis;
     *  \param[out] side if `true`, right-hand side of row `m_iBasicRow[s]` is tight (`s < m_iBasisSize`),
     *  or variable `x[m_ipBasicColumn[s]]` must take its upper-bound value (`s >= m_iBasisSize`);
     *  otherwise, left-hand side of row `m_iBasicRow[s]` is tight (`s < m_iBasisSize`),
     *  or variable `x[m_ipBasicColumn[s]]` must take its lower-bound value (`s >= m_iBasisSize`);
     * \param[out] bDegen if upcoming simplex iteration will be degenerate.
     */
    void primeSearchPivot(int& s, int& t, bool& side, bool& bDegen);

    /**
     *  The procedure is used only in `dualRatioTest()` to compute the dual direction stored in `m_dpUd`.
     *
     *  \param[in] s is `s > 0`, row `s-1` leaves basis, otherwise, column `-s-1` leaves basis;
     *  \param[out] dMaxVal maximum absolute value of coefficients in `m_dpUd`.
     *  \return always `0`.
     */
    int computeUd(int s, double& dMaxVal);

    /**
     *  The procedure is called only in `dualRatioTest()` to compute the dual direction stored in `m_dpUd`,
     *  as well as the updates for row norms stored in `m_dpArray`.
     *
     *  \param[in] s is `s > 0`, row `s-1` leaves basis, otherwise, column `-s-1` leaves basis;
     *  \param[out] dMaxVal maximum absolute value of coefficients in `m_dpUd`.
     *  \return `SHIFT+1` if row norms must be updated; otherwise, `0`.
     */
    int computeUdAndRowNormUpdates(int s, double& dMaxVal);

    /**
     * The procedure is called from `dualSeekPivot()` to correct the current basis so that to make it dual feasible.
     *
     * \param[in] i index of dual variable (shadow price or reduced cost) such that `m_dpPi[i]` of the wrong sign.
     */
    void makeBasisDualFeasible(int i);

    /**
     * \param[in] s specifies index of violated (by current basic solution) constraint:
     *     if `s > 0`, then inequality `s` is violated; otherwise, if `s < 0`, then one bound of variable `-s-1` is not satisfied;
     *  \param[in] side if `true`, then inequality right hand side or variable upper bound is not satisfied;
     *  otherwise, inequality left hand side or variable lower bound is violated;
     *  \param[out] degen if `true`, then this simplex iterate is degenerate.
     *
     *  \return in case of success, index of basic constraint or variable to be removed from the basis;
     *   `SHIFT+1` if norm of the pivot row is not correct (a signal to recompute row norms).
     */
    int dualRatioTest(int s, bool side, bool& degen);
 
    /**
     * The function is called only in `dualSeekPivot()` to compute `Ub`.
     *
     * \param[in] s if `s > 0`, row `s-1` enters basis; if `s < 0`, column `-s-1` leaves basis;
     * \param[in] t if `t < m_iBasisSize`, row `m_ipBasicRow[t]` leaves basis;
     *  if `t >= m_iBasisSize`, column `m_ipBasicColumn[t]` enters basis.
     * \return pivot value, which is computed twice: first when computing `Ud`, and next when computing `Ub`;
     * if these two values are different, the basis is refactored.
     */
    double dualComputeUb(int s, int t);

    /**
     * The procedure first calls `innerSeparate()` to find an inequality violated by the current basic solution,
     * or a basic variable which current value violates one of its bounds.
     * Next, finds a row to leave the basis, or a column to enter the basis.
     *
     * \param[out] s if `s < m_iBasisSize`, then row `m_ipBasicRow[s]` is to enter basis;
     * if `s >= m_iBasisSize`, then column `m_ipBasicColumn[s]` is to leave basis;
     * \param[out] t if `t < m_iBasisSize`, then row `m_ipBasicRow[t]` is to leave basis;
     * if `t >= m_iBasisSize`, then column `m_ipBasicColumn[t]` is to enter basis;
     * \param[out] side if `true`, right-side inequality of constraint `m_ipBasicRow[s]` (`s < m_iBasisSize`)
     *  or upper bound of variable `m_ipBasicColumn[s]` (`s >= m_iBasisSize`) is used;
     *  otherwise, left-side inequality or lower bound is used;
     * \param[out] degen if `true`, upcoming dual-simplex iteration is degenerate, and not degenerate otherwise;
     * \param[in] needSol if set to `false`, in case of numeric instability,
     * `dualSeekPivot()` immediately terminates without trying to rectify it.
     *
     */
    void dualSeekPivot(int& s, int& t, bool& side, bool& degen, bool needSol=true);

    /**
     * The procedure extracts from the (sparse) matrix a given row.
     *
     * \param[in] i row index;
     * \param[out] d array to store output row.
     * \sa unpackRowForBasis().
     */
    void unpackRow(int i, double* d);

    /**
     * The procedure extracts from the (sparse) matrix a given row,
     * the row entries are stored in the order given by `m_ipBasicColumn`.
     *
     * \param[in] s row index;
     * \param[out] U array to store output row.
     * \sa unpackRow().
     */
    void unpackRowForBasis(int s, double* U);

    /**
     * The procedure extracts from the (sparse) matrix the basic components of a given column,
     * the column entries are stored in the order given by `m_ipBasicRow`.
     *
     * \param[in] s column index;
     * \param[out] U array to store output column.
     * \param[in] ExtraRow, v if non-negative and `v != 0`, entry from row `ExtraRow` is assigned to `*v`.
     * \sa unpackRow().
     */
    void unpackColumnForBasis(int s, double* U, int ExtraRow=-1, double* v=0);

    /**
     * if bSafeMode == true, LU::Factor is called in the safe mode,
     * i.e. with (full) column permutations to reduce rounding off errors
    */
    void updatePartition(bool bSafeMode=false);

    /**
     * This procedure updates current dual solution.
     *
     * \param[in] r pivot row;
     * \param[in] step value of step: `m_dpPi = m_dpPi + step * m_dpUd` and `mPdpPi[r]=step`.
     */
    void updateDualSolution(int r, double step);

    /**
     * This procedure updates current prime solution.
     *
     * \param[in] step value of step: `m_dpX = m_dpX + step * m_dpUb`.
     */
    void updatePrimalSolution(double step);

    /**
     * \param[in,out] y on input, if not zero, vector of basic components of dual solution; on output, full (including reduced costs) dual solution;
     * \param[in] how =1 for column base extension, =-1 for row base extension, = 0 to be decided in the procedure;
     * \param[in] bWithFixed if `false`, the reduced costs of fixed variables are not computed.
     */
    void extendDualSolution(double* y=0, int how=0, bool bWithFixed=false);

    void extendPrimeSolution(double* x=0, int col=-1, int how=0); ///< Set how = 1 in case of partial separation


    /**
     * The procedure performs the pivot iteration when a row is substituted for a row.
     *
     * \param[in] r row `m_ipBasicRow[r]` leaves basis;
     * \param[in] s row that leaves basis;
     * \param[in] side if `true`, right-hand side inequality of constraint `s` is used, and left-hand side is used otherwise.
     *
     * \return
     *      - `0` in rare simple cases when the right-hand side inequality of constraint `s` is substituted for its left-hand side inequality or vice versa;
     *      - '1` if it is not needed to refactor the basic matrix;
     *      - `-1` if it is needed to refactor the basic matrix;
     *      - `-2` if it is needed to enter the safe mode and repeat the simplex iteration.
     */
    int substituteRowForRow(int r, int s, bool side);

    /**
     * The procedure performs the pivot iteration when a row is substituted for a column.
     *
     * \param[in] r row `m_ipBasicRow[r]` leaves basis;
     * \param[in] s column that leaves basis;
     * \param[in] side if `true`, `x[s]` is to take value of its upper bound, and lower bound otherwise.
     *
     * \return
     *      - `0` in rare simple cases when the value of variable `x[s]` is changed from its lover to its upper bound or vice versa;
     *      - '1` if it is not needed to refactor the basic matrix;
     *      - `-1` if it is needed to refactor the basic matrix;
     *      - `-2` if it is needed to enter the safe mode and repeat the simplex iteration.
     */
    int substituteRowForColumn(int r, int s, bool side);

    /**
     * The procedure swaps two non-basic columns.
     *
     * \param[in] c1,c2 column indices.
     */
    void swapNonBasicColumns(int c1, int c2);

    /**
     * The procedure swaps two non-basic rows.
     *
     * \param[in] r1,r2 row indices.
     */
    void swapNonBasicRows(int r1, int r2);

    /**
     * The procedure performs the pivot operation when a column is substituted for a row.
     * \param[in,out] r column `c=m_ipBasicColumn[r]` is to be replaced by row `s`;
     *  since `substituteColumnForRow()` swaps columns, on return, `r` gives a new position of `c`
     *  in array `m_ipBasicColumn`;
     * \param[in] s row to enter basis;
     * \param[in] side if `true`, right hand side constraint given by `row` enters basis;
     *   otherwise, right hand side constraint enters basis.
     * \return
     *      - '1` if it is not needed to refactor the basic matrix;
     *      - `-1` if it is needed to refactor the basic matrix;
     *      - `-2` if it is needed to enter the safe mode and repeat the simplex iteration.
     */
    int substituteColumnForRow(int& r, int s, bool side);

    /**
     * The procedure performs the pivot operation when a column is substituted for a column.
     * \param[in,out] r column `c=m_ipBasicColumn[r]` is to be replaced by column `s`;
     * \param[in] s column to enter basis;
     * \param[in] side if `true`,  variable `x[c]` is to take its upper bound value;
     *   otherwise, x[c]` is to take its lower bound value.
     * \return
     *      - '1` if it is not needed to refactor the basic matrix;
     *      - `-1` if it is needed to refactor the basic matrix;
     *      - `-2` if it is needed to enter the safe mode and repeat the simplex iteration.
     */
    int substituteColumnForColumn(int r, int s, bool side);

    /**
     * The procedure performs the pivot operation when auxiliary arrays are dense.
     *
     * \param[in] r if `r < m_iBasisSize`, row `m_ipBasicRow[r]` leaves basis;
     * if `r >= m_iBasisSize`, column `m_ipBasicColumn[r]` enters basis;
     * \param[in] s if `s > 0`, row `s-1` enters basis;
     * if `s < 0`, column `-s-1` leaves basis;
     * \param[in] side if `side == true`,
     *     - in case `s > 0`, right-hand side inequality `s-1` is tight,
     *     - in case `s < 0`, variable `x[-s-1]` takes its upper bound value;
     * if `side==false`,
     *     - in case `s > 0`, left-hand side inequality `s-1` is tight,
     *     - in case `s < 0`, variable `x[-s-1]` takes its lower bound value.
     *
     * \sa sparseDoPivot().
     */
    void doPivot(int s, int r, bool side);

    /**
     * The procedure performs the pivot operation when auxiliary arrays are sparse.
     * It has the same parameters as `doPivot()`.
     */
    void sparseDoPivot(int s, int r, bool side);

protected:
    /**
     * The function transforms a non-scaled basic solution for the LP in the memory
     * into a solution of the original LP or MIP. To do its job, the function uses expressions
     * in the preprocessing stack to recover the values of those variables that were eliminated
     *  from the matrix during preprocessing.
     * \param[in] dpBasicX basic solution;
     * \param[out] dpX array of size `m_iN` that sores extended solution;
     * `dpX[i]` is value of variable with handle `m_ipColHd[i]`.
     * \return
     */
    int getPrimeSolution(double* dpBasicX, double* dpX);

    /**
     * \param[in] j index of variable;
     * \return value of variable `j`.
     */
    double getVarValue(int j) const;

    /**
     * \param[in] i row index;
     * \return value of row `i`.
     */
    double getRowValue(int i) const;

public:
    /**
     * \param[in] i row index;
     * \param[in] dpX array of size `m_iN`.
     * \return the scalar product of the constraint vector of row `i` and the vector stored in `dpX`.
     */
    double getRowValue(int i, double * dpX) const;

protected:

    /**
     * The function computes the expression of a given vector in the current basis.
     * If `col >= 0`, then the expression (in current basis) of the unit vector with 1 in position `col` is computed.
     * \param[in,out] x array of size `m_iN` storing vector.
     * \param[in] col column index.
     * \remark This function is used internally in `CLP` and `CMIP`.
     * \sa computeY().
     */
    void computeX(double* x=0, int col=-1);

    /**
     * The function computes the expression of a given dual vector in the current basis.
     * \param[in,out] y array of size `m_iN` storing vector.
     * \param[in] withFixedVars if `false` reduced costs of fixed variables are not computed.
     * \remark This function is used internally in `CLP` and `CMIP`.
     * \sa computeX().
     */
    void computeY(double* y=0, bool withFixedVars=false);

    /**
     * The function re-factors the basis, and then recomputes both, prime and dual, basic solutions.
     */
    void updateSolution();

    /**
     * The function computes the vector weighted sum of all non-basic columns,
     *  i.e. each column is first multiplied by the value of corresponding variable.
     *  \sa incColumnSum().
     */
    void computeBasicColumnSum();

    /**
     * The function adds a given matrix column, `col`, multiplied by `delta` to the sum of basic columns.
     * \param[in] col column index;
     * \param[in] delta multiplier.
     */
    void incColumnSum(int col, double delta);

private:
//////////////////////////////////////////////////////////////////////////
//                      Prime and Dual Simplex Methods                  //
//////////////////////////////////////////////////////////////////////////
    /**
     * \param[in] col column index;
     * \param[in] dCost column cost,
 	 * \param[in] dpPi vector of size `m_iBasisSize`, row basic potentials.
     * \return reduced cost of variable indexed with `col`.
     */
    double getReducedCost(int col, double dCost, double* dpPi);

    /**
     * The procedure deletes from the matrix all rows and columns for previously deleted constraints and variables.
     */
    void compressBasis();

    /**
     * The procedure deletes from the matrix all columns corresponding to slack variables.
     */
    void deleteSlackVars();

    /**
     * The procedure extends the basis to make it prime feasible by adding slack variables.
     * \param[in] bigM big `M` value;
     * \param[in] dpX array of size `m_iN`, basic solution to be extended.
     */
    void addSlackVars(double bigM, double* dpX);


//////////////////////////////////////////////////////////////////////////
// Perturbation of the objective, row and column bounds to prevent Cycling
/////////////////////
    void storeObj(); ///< The procedure stores the objective coefficients before perturbing their values.

    /**
     * The procedure restores objective coefficients previously stored  by `storeObj()`.
     *
     * \param[in] updateY if `true`, dual solution is recomputed.
     */
    void restoreObj(bool updateY);

    void prtbObj(); ///< The function randomly stirs up objective coefficients of non-basic variables.

    /**
     * The procedure stores the bounds of all variables before perturbing their values.
     *
     * \return `true` in case of success.
     */
    bool storeBounds();

    void prtbBounds(); ///< The function randomly stirs up bounds of non-basic variables.

    /**
     * The procedure restores the bounds of variables previously stored  by `storeBounds()`.
     *
     * \param[in] updateX if `true`, primal solution is recomputed.
     */
    void restoreBounds(bool updateX);
///////////////////////////////////////////////////////////////////////////
public:
    /**
     * The function sets the solution algorithm to be used for solving LPs and root LPs of MIPs.
     * \param[in] method solution algorithm.
     * \sa enLPmethod, getCurrentLPmethod().
     */
    void setLPmethod(enLPmethod method)
        {m_eDefLPmethod=method;}

    /**
     * \return currently used LP algorithm.
     */
    enLPmethod getCurrentLPmethod() const
        {return m_eLPmethod;}

    /**
     * This procedure does a lot of useful things: allocates memory, applies different preprocessing techniques,
     * scales the matrix, and so on. As a results, in many cases, we can get substantially simple (for solving)
     * problem, which is still equivalent to the original (user) problem.
     * The latter means that, given an optimal solution to the transformed problem, we can
     * easily compute an optimal solution to the original problem.
     * \return  `false` if inconsistency has been detected; otherwise, `true`.
     * \throws CMemoryException lack of memory.
    */
    virtual bool prepare();

    /**
     * The procedure implements a prime simplex algorithm.
     * \param[in] timeToStop procedure will stop at time `timeToStop` (given in seconds since epoch); `timeToStop==0l` means no time limitation;
     *  if you want `primeSimplex()` to stop running in `timeLimit` seconds, set `timeToStop=getStartTime()/1000l + timeLimit`;
     * \param[in] upperBound procedure stops if objective value exceeds `upperBound`;
     * \param[in] maxItNum  procedure stops after carrying out at most `maxItNum` iterations;
     * \param[in] degCheckInterval,maxDegPrc degeneracy test is activated after `degCheckInt` consecutive iterations;
     * if more than `maxDegPrc`% iterations among the last `degCheckInterval` ones were degenerated, bounds on variables are perturbed.
     * \return
     *   - `0` solution found;
     *   - `-1` objective function upper bound, `upperBound`, exceeded;
     *   - `2` limit on number of iterates,  `maxItNum`, exceeded;
     *   - `3` time limit, `timeLimit`, exceeded.
     *   \throws CMemoryException lack of memory.
     *   \sa dualSimplex().
     */
    int primeSimplex(__LONG timeToStop=0l, double upperBound=INF, int maxItNum=1000000000,
    		int degCheckInterval=250, int maxDegPrc=90);

    /**
     * The procedure implements a dual simplex algorithm.
     * \param[in] timeToStop procedure will stop running at time `timeToStop` (given in seconds since epoch); `timeToStop==0l` means no time limitation;
     *  if you want `dualSimplex()` to stop running in `timeLimit` seconds, set `timeToStop=getStartTime()/1000l + timeLimit`;
     * \param[in] lowerBound procedure stops if objective value becomes less than `lowerBound`
     * \param[in] maxItNum  procedure stops after carrying out at most `maxItNum` iterates;
     * \param[in] needSol if set to `true` and if solution has been found, on termination basic matrix is re-factored to eliminate computational errors;
     * \param[in] inconsistSertificate if set to `true`, on termination, even if problem has been proven to be infeasible, basic matrix is re-factored to eliminate computational errors;
     * \param[in] degCheckInterval,maxDegPrc degeneracy test is activated after `degCheckInterval` consecutive iterations;
     * if more than `maxDdegPrc`% iterations among the last `degCheckInterval` ones were degenerated, objective is perturbed.
     * \return
     *   - `0` solution found;
     *   - `1` inconsistency detected;
     *   - `2` limit on number of iterates,  `maxItNum`, exceeded;
     *   - `3` objective function is less than `lowerBound`;
     *   - `4` time limit, `timeLimit`, exceeded.
     *   \throws CMemoryException lack of memory.
     *   \sa primeSimplex().
     */
    int dualSimplex(__LONG timeToStop=0l, double lowerBound = -INF, int maxItNum=1000000000,
             bool needSol=true, bool inconsistSertificate=false,
             int degCheckInterval=250, int maxDegPrc=90);

private:
    /**
     * The repeatedly calls `dualSimplex()` and `primeSimplex()` until prime and dual feasibility of the basis is achieved.
     * \param[in] timeToStop procedure will stop running at time `timeToStop` (given in seconds since epoch);
     * \param[in] bRowGen if `true`, row (cut) generation is allowed;
     * \param[in] bColGen if `true`, column generation is allowed.
     * \throws CMemoryException lack of memory
     */
    int solveProblem(__LONG timeToStop, bool bRowGen=false, bool bColGen=false);

    /**
     * 1) calls `initPrimalBasis()` to get an initial primal feasible basis;
     * 2)   solve LP by the prime simplex method;
     * 3) deletes all slack variables introduced by ``initPrimalBasis()`.
     */
    int PRIME_solve(__LONG timeToStop);

///////////////////////////////////////////////////////
// Initial Basis
///////////////////////////////////////////////////////

    /**
     * The procedure estimates whether changing the value of a given
     * non-basic variable increases the objective, and if so,
     * the value of this variable is set to its opposite bound value
     * (from lower to upper and vice versa).
     *
     * \param[in] col index of non-basic variable.
     * \return `true` if the basis has been changed, and `false` otherwise.
     */
    bool INIT_changeSide(int col);

    /**
     * This function computes an optimal solution to the LP with only one inequality
     * corresponding to a given problem constraint.
     * \param[in] row constraint index;
     * \param[in] side is `true` for right hand side constraint, and `false` for left hand side constraint;
     * \param[out] pivot pivot value.
     * \return index of a variable to be included in the column basic set.
     */
    int lpGreedy(int row, bool side, double &pivot);

    /**
     * The procedure tries to move to the dual initial basis all unbounded variables.
     * This procedure is called only from `initDualBasis()`.
     *
     * \param[out] freeRowNum number of non-basic rows that do not contain basic variables.
     * \return `-1` in case the problem has been proven to be unbounded, and `0` otherwise.
     */
    int unboundedVarsToBasis(int &freeRowNum);

    /**
     * The procedure is called from `initDualBasis()` to build an initial LP dual feasible basis.
     *
     * \param[in] freeRowNum number of rows without entries in basic columns.
     */
    void crackBasis(int freeRowNum);

    virtual void initPrimalBasis(); ///< Builds an initial LP primal feasible basis by introducing slack variables.

    /**
     * The procedure builds an initial LP dual feasible basis.
     *
     * \return `-1` in case the problem has been proven to be unbounded, and `0` otherwise.
     */
    virtual int initDualBasis();

protected:
    /**
     * The procedure first computes an initial basis, and then implements an iterative procedure
     * that, until the current basic solution is prime and dual feasible,
     * applies either `primeSimplex()` or `dualSimplex()`.
     * \param[in] timeLimit procedure will stop running in at most `timeLimit` seconds; `timeLimit==0l` means no time limitation;
     * \param[in] genFlag if `false`, both, inequalities and columns, are not generated.
     * \return 0 in case a solution has been found; otherwise, the return value is that of
     * `primeSimplex()` or `dualSimplex()` depending which of them returned a non-zero value.
     * \sa optimize().
     */
    int solveLP(__LONG timeLimit=0, bool genFlag=true);
///////////////////////////////////////////////////////
//       S E P A R A T I O N
///////////////////////////////////////////////////////
public:
    /**
     * Call this function to prevent generation of any inequalities.
     * \attention The function should be called before `openMatrix()`, which allocates memory.
     * Otherwise, unnecessary additional amount of memory will be allocated.
     */
    void switchOffRowGen()
        {m_uRowColGenRule&=COL_GEN;}

    void switchOnColGen()
        {m_uRowColGenRule&=~COL_GEN;} ///<  Call this function to allow generation of new columns.

protected:
    /**
     * This function is used to add new constraints to the closed matrix.
     * \param[in] hd handle of constraint;
	 * \param[in] type type of constraint;
	 * \param[in] b1,b2  respectively, left (LHS) and right (RHS) hand sides of constraint;
 	 * \param[in] sz number of non-zero entries;
	 * \param[in] dpVal, ipCol references to arrays of size `sz`, `dpVal[i]` is coefficient at variable `ipCol[i]`;
	 * \param[in] bVarScaled  if `true`, constraint is expressed in scaled variables;
 	 * \param[in] factor scale factor; `NOT_SCALED` value means that constraint is to be scaled by __MIPCL__ itself;
 	 * \param[in] n if non-zero, number of columns in original (non-preprocessed) problem;
 	 * in such a case the expressions in the preprocessing stack are used to substitute the variables `m_iN,...,n-1`;
 	 * \param[in] toBasis if `true`, data structures representing basis are slightly modified;
 	 * in rare cases when the matrix has been closed already, but optimization has not been started yet,
 	 *  `  toBasis` must be set to `false`.
 	 *  \return index of newly created constraint.
 	 * \throws CMemoryException lack of memory.
 	 * \attention The procedure adds the constraint to the matrix but not to the MIP pool.
 	 * \sa addRow(), getRow(), CMIP::genCut1(), CMIP::genCut2().
     */
    int addNewRow(tagHANDLE hd, unsigned type, double b1, double b2,
                  int sz, double* &dpVal, int* &ipCol,
                  bool bVarScaled=true, int factor=NOT_SCALED, int n=0,
                  bool toBasis=true);

    /**
     * The function may be useful for adding to the matrix several very similar constraints.
     */
    int dublicateRow(int row); ///< Adds to the matrix another copy of row `row`.

    /**
     *
     * \param[in] hd handle of the variable; if `hd >= 0`, then you must overload
     *  `getColumn()` (version with 11 parameters);
     * \param[in] type type of variable;
     * \param[in] cost objective coefficient;
     * \param[in] l,u lower and upper bounds of variable;
     * \param[in] sz number of nonzero entries;
     * \param[in] dpVal,ipRow arrays of size `sz`; `dpVal[i]` is coefficient in row `ipRow[i]`;
     * \param[in] side  if `true`, then in current basis variable takes value of `u`;
     *           otherwise, variable takes value of `l`;
     * \param[in] scaled if `true`, column is given in scaled variables;
     * \param[in] factor   - scale factor, set `factor=NOT_SCALED` if column must be scaled;
     * \param[in] flag set it to `true`, only if `addNewColumn()` is called within `generateColumns()`;
     * \return index of newly created column (variable).
     * \throws CMemoryException lack of memory.
     */
    int addNewColumn(tagHANDLE hd, unsigned type, double cost, double l, double u,
                      int sz, double* &dpVal, int* &ipRow, 
                      bool side, bool scaled=false, int factor=0, bool flag=false);

    /**
     * This function is helpful in branch-and-price algorithms when pricing columns that are not in the matrix.
     * \param[in] sz number of non-zero coefficients in column;
     * \param[in] dpVal,ipRow arrays of size `sz`; `dpVal[i]` is coefficient in row `ipRow[i]`.
     * \return if `c` is the cost coefficient of corresponding variable and `q` denotes return value, then
     * `c-q` is reduced cost of that variable.
     */
    double estimateCol(int sz, double* dpVal, int* ipRow);

    /**
     * Usually, this function is used for generating (adding to the matrix) strong inequalities, i.e.,
     * those that are a part of problem formulation, or those that are facet-defining cuts, and so on.
     * \param[in] n number of variables;
     * \param[in] dpX,ipColHd arrays of size `n`; `dpX[i]` is value of variable having handle `ipColHd[i]`;
     * \param[in] genFlag if `true`, all generated inequalities are added to matrix by calling `addNewRow()`
     *  (or `CMIP::addCut()` in MIP applications); otherwise,
     * none inequality should be added;
     * \return
     *   - `genFlag==true`: `true` if at least one inequality were generated; otherwise, `false`;
     *   - `genFlag==false`: `false` (!!!) if solution presented in `dpX,ipColHd`  is feasible,
 	 * or `true` if it is infeasible.
 	 * \throws CMemoryException lack of memory.
     * \remark Default implementation does nothing.
 	 * \sa addNewRow(), CMIP::genCut1(), CMIP::genCut2().
     */
    virtual bool separate(int n, const double* dpX, const tagHANDLE* ipColHd, bool genFlag);
//////////////////////////////////////////////////////////////////////////
//                           Column Generation                          //
//////////////////////////////////////////////////////////////////////////
    /**
     * This is used when working with preprocessed matrix;
     * otherwise, `getShadowPrices()` should be used.
     * \param[in] z vector of potentials;
     * \param[out] y array of size `m_iM`, `y[i]` dual variable for constraint `i`;
     * \param[in] bScaled if `true`, output `y` is computed for scaled problem; otherwise, for not-scaled.
     * \sa getShadowPrices().
     */
    void getDualRowVars(double* y, double* z, bool bScaled);

    /**
     * This function must be overloaded in any derived class that generates columns and stores them
     * in its own pool. If a new column is added to the matrix with a positive handle,
     * later, the solver may call `getColumn()` to reload this column to the matrix.
     * \param[in] hd handle of column to be returned;
     * \param[in] m number of rows in matrix;
     * \param[in] ipHd array of size `n`; `ipHd[i]` is handle of row `i`;
     * \param[out] type type of corresponding variable;
     * \param[out] cost objective coefficient;
     * \param[out] l,u  lower and upper bounds of corresponding variable;
     * \param[out] sz number of non-zero coefficients in column;
     * \param[out] dpVal, ipRow arrays of size `sz`; `dpVal[i]` is coefficient in row `ipRow[i]`.
      * \return `true` if column with handle `hd` was successfully restored; otherwise, `false`.
     * \sa addNewColumn().
     */
    virtual bool getColumn(tagHANDLE hd, int m, const tagHANDLE* ipHd,
                unsigned& type, double& cost, double& l, double& u,
                int& sz, double* dpVal, int* ipRow)
        {return false;}

    /**
     * This function must be overloaded in any derived class that generates rows (cuts) and stores them
     * in its own pool. If a new inequality is added to the matrix with a positive handle,
     * later, the solver may call `getRow()` to reload this inequality to the matrix.
     * \param[in] hd handle of constraint (row) to be returned;
     * \param[in] n number of columns in matrix;
     * \param[in] ipHd array of size `n`; `ipHd[j]` is handle of column `j`;
     * \param[out] type type of constraint;
     * \param[out] lhs,rhs  LHS and RHS;
     * \param[out] sz number of variables in constraint;
     * \param[out] dpVal, ipCol arrays of size `sz`; `dpVal[i]` is coefficient in column `ipCol[i]`.
     * \param[out] scaled if `true` the constraint (row) is returned in scaled variables; otherwise, `false`.
     * \return `true` if constraint with handle `hd` was successfully restored; otherwise, `false`.
     * \sa addNewRow().
     */
    virtual bool getRow(tagHANDLE hd, int n, const tagHANDLE* ipHd,
                        unsigned& type, double& lhs, double& rhs,
                        int& sz, double* dpVal, int* ipCol, bool &scaled)
        {return false;}

    /**
     * The function must be overloaded in user applications generating columns.
     * The functions is to produce a number of columns (presently not in the matrix)
     * with negative reduced costs defined to be the cost of the column minus the value
     * returned by `estimateColumn()`.
     * \param[in] m number of rows in matrix;
     * \param[in] ipRowHd,dpY arrays of size `m`; `dpY[i] ` is shadow price (value of dual variable)
     *  of constraint (row) having handle `ipRowHd`.
     * \return `true` if at least one column was generated; otherwise, `false`.
     * \attention The procedure must add columns to the matrix by calling `addNewColumn()`
     *  with parameter `flag` set to `false`.
     */
    virtual bool generateColumns(int m, const tagHANDLE* ipRowHd, const double* dpY);

    tagHANDLE getVarHandle(int j) const
        {return m_ipColHd[j];} ///< \return handle of variable indexed by `j`.

    tagHANDLE getCtrHandle(int i) const
        {return m_ipRowHd[i];} ///< \return handle of constraint indexed by `i`.

    /**
     * Overload this function to allow `CLP` to use symbolic constraint names when storing solutions.
     * \param[in] hd constraint (row) handle;
     * \param[out] name memory to store returned name;
     * \return pointer to `name` parameter.
     * \remark Default implementation composes the names as the concatenation
     *  of string "ctr_" with the value of `hd`.
     * \sa getVarName().
     */
    virtual char* getCtrName(tagHANDLE hd, char* name) const;

    /**
     * Overload this function to allow `CLP` to use symbolic names of variables when storing solutions.
     * \param[in] hd handle of variable (column);
     * \param[out] name memory to store returned name;
     * \return pointer to `name` parameter.
     * \remark Default implementation represents the name as the string "x(q)",
     *  where `q` is the value of `hd`.
     * \sa getCtrName().
     */
    virtual char* getVarName(tagHANDLE hd, char* name) const;

public:
    /**
     * \param[in] j index (not handle) of variable;
     * \param[in] l,u new lower and upper bounds for variable `j`.
     * \attention Changing bounds for variables may make the basis prime or dual infeasible.
     */
    void setVarBounds(int j, double l, double u);

    /**
     * \param[in] j index of a variable;
     * \param[in] l new lower bound for variable `j`.
     * \attention Changing bounds for variables may make the basis prime or dual infeasible.
     */
    void setVarLoBound(int j, double l);

    /**
     * \param[in] j index of a variable;
     * \param[in] u new upper bound for variable `j`.
     * \attention Changing bounds for variables may make the basis prime or dual infeasible.
     */
    void setVarUpBound(int j, double u);

    /**
     * The procedure removes both, lower and upper, bounds.
     * Sets the lower bound to minus infinity, and the upper to infinity.
     * \attention Removing  bounds for variables may make the basis dual infeasible.
     * \sa setVarBounds().
     */
    void setVarFree(int j)
        {m_dpD[j<<1]=-(m_dpD[(j<<1)+1]=m_dVarInf); m_ipVarType[j]&=~(VAR_LEFT | VAR_RIGHT | VAR_FX);}

    /**
     * The procedure sets both, lower and upper, bounds for a given constraint.
     * \param[in] i constraint index (not handle);
     * \param[in] lhs,rhs new left and right hand sides for constraint `i`.
     * \attention Changing bounds for constraints may make the basis prime or dual infeasible.
     */
    void setCtrBounds(int i, double lhs, double rhs);

    /**
     * The procedure sets the lower (left hand side) bound for a given constraint.
     * \param[in] i constraint index (not handle);
     * \param[in] lhs new left hand side for constraint `i`.
     * \attention Changing bounds for constraints may make the basis prime or dual infeasible.
     */
    void setLHS(int i, double lhs);

    /**
     * The procedure sets the upper (right hand side) bound for a given constraint.
     * \param[in] i constraint index (not handle);
     * \param[in] rhs new right hand side for constraint `i`.
     * \attention Changing bounds for constraints may make the basis prime or dual infeasible.
     */
    void setRHS(int i, double rhs);

    /**
     * The procedure sets the left hand side to minus infinity, and the right hand side to infinity.
     * \attention Removing  bounds for constraints may make the basis dual infeasible.
     * \sa setCtrBounds().
     */
    void setCtrFree(int i)
        {m_dpB[i<<1]=-(m_dpB[(i<<1)+1]=INF); m_ipCtrType[i]&=~(CTR_LEFT | CTR_RIGHT | CTR_EQ);}

protected:
    /**
     * \param[in] j index of variable;
     * \return `true` if variable `j` is _fixed_ (its lower and upper bounds are equal), otherwise, `false`.
     */
    bool isVarFixed(int j) const
        {return (m_ipVarType[j] & VAR_FX)? true: false;}

    /**
     * \param[in] j index of variable;
     * \return `true` if variable `j` is _free_, otherwise, `false`.
     */
    bool isVarFree(int j) const
        {return (m_ipVarType[j] & (VAR_LEFT | VAR_RIGHT))? false: true;}

    /**
     * \param[in] j index of variable;
     * \return `true` if variable `j` is _bounded_ (both its bounds, lower and upper, are finite), otherwise, `false`.
     */
    bool isVarBounded(int j) const
        {return ((m_ipVarType[j] & VAR_LEFT) && (m_ipVarType[j] & VAR_RIGHT))? true: false;}

    /**
     * The function is used to verify whether a variable is lower bounded.
     * \param[in] j index of variable;
     * \return `true` if variable `j` is _lower bounded_ (its lower bound is finite), otherwise, `false`.
     */
    bool isVarLoBounded(int j) const
        {return (m_ipVarType[j] & VAR_LEFT)? true: false;}  

    /**
     * The function is used to verify whether a variable is upper bounded.
     * \param[in] j index of variable;
     * \return `true` if variable `j` is _upper bounded_ (its upper bound is finite), otherwise, `false`.
     */
    bool isVarUpBounded(int j) const
        {return (m_ipVarType[j] & VAR_RIGHT)? true: false;}

/////////////////////////////////////////////////////////////////////////////
    /**
     * The function is used to verify whether a constraint is an equality constraint.
     * \param[in] i index of constraint;
     * \return `true` if constrain `i` is _equality_ (its left  and right hand are equal), otherwise, `false`.
     */
   bool isCtrEq(int i) const
        {return (m_ipCtrType[i] & CTR_EQ)? true: false;}

   /**
    * The function is used to verify whether a constraint is _free_.
    * \param[in] i index of constraint;
    * \return `true` if constrain `i` is _free_, otherwise, `false`.
    */
    bool isCtrFree(int i) const
        {return (m_ipCtrType[i] & (CTR_LEFT | CTR_RIGHT))? false: true;}

    /**
     * The function is used to verify whether a constraint is two side bounded.
     * \param[in] i index of constraint;
     * \return `true` if constrain `i` is _bounded_ (its left  and right hand are finite), otherwise, `false`.
     */
    bool isCtrBounded(int i) const
        {return ((m_ipCtrType[i] & CTR_LEFT) && (m_ipCtrType[i] & CTR_RIGHT))? true: false;}

    /**
     * The function is used to verify whether a constraint is lower bounded.
     * \param[in] i index of constraint;
     * \return `true` if constrain `i` is _lower bounded_ (its left hand is finite), otherwise, `false`.
     */
    bool isCtrLoBounded(int i) const
        {return (m_ipCtrType[i] & CTR_LEFT)? true: false;}

    /**
     * The function is used to verify whether a constraint is upper bounded.
     * \param[in] i index of constraint;
     * \return `true` if constrain `i` is _upper bounded_ (its right hand is finite), otherwise, `false`.
     */
    bool isCtrUpBounded(int i) const
        {return (m_ipCtrType[i] & CTR_RIGHT)? true: false;}

/////////////////////////////////////////////////
//   Access the solution
/////////////////////////////////////////////////
public:
    /**
     * The function resumes the state in which the problem was before the solution procedure has started.
     */
    virtual void reset();

    /**
     * The procedure solves LPs.
     * \param[in] timeLimit limit on solution time (in seconds);
     * \param[in] gap integrality gap;
     * \param[in] solFileName pointer to string with file name for storing intermediate solutions.
     * \remark When solving LPs, input parameters are not used.
     * \sa CMIP::optimize().
     */
    virtual void optimize(__LONG timeLimit=1000000l, double gap=0.0, const char *solFileName=0);

    /**
     * \return `true` if LP (MIP) in memory has been prepared for optimization.
     */
    bool isPrepared() const
    {return (m_iState & PROB_PREPARED)? true: false;}

    /**
     * \return `true` if LP (or MIP) in memory has been solved.
     */
    bool isSolved() const
    {return (m_iState & PROB_SOLVED)? true: false;}
    
    /**
     * \return `true` if LP in memory has been solved and optimal solution found.
     * \sa CMIP::isSolution().
     */
    virtual bool isSolution() const
        {return (m_iState & PROB_SOLVED) && m_bPrimeFeasible && m_bDualFeasible;}

    /**
     * \return `true` if LP in memory has been solved and proven to be _infeasible_.
     * \sa whyLpInfeasible(), showWhyLpInfeasible().
     */
    bool isLpInfeasible() const
    {return (m_iS != 0) && m_bDualFeasible;}

    /**
     * \return `true` if this LP or MIP is infeasible (has no solution).
     */
	bool isInfeasible() const
	{return (m_iState & PROB_INFEASIBLE)? true: false;}


    /**
     * \return `true` if LP in memory has been solved and proven to be _unbounded_.
     * \sa whyLpUnbounded(), showWhyLpUnbounded().
     */
    bool isLpUnbounded() const
        {return (m_iBasisSize == -1) || (m_iS >= 0 && m_bPrimeFeasible);}

    /**
     * This function is usually called when the problem has been solved already.
     * \return optimal objective value of solved LP.
     */
    double getObjVal() const;

    /**
     * The function returns two pointers to the internal MIPCL arrays storing the LP solution.
     * \param[out] dpX,ipHd `dpX[j]` is value of variable whose handle is `ipHd[j]`, `j=1,...,n`,
     *    where `n` is return value.
     * \return number of variables.
     * \attention Do not modify the values of both arrays  `dpX` and `ipHd`.
     * \sa CMIP::getSolution().
     */
    int getSolution(double* &dpX, int* &ipHd);

    /**
     * The function returns _reduced costs_ of variables for the solved LP.
     * \param[out]  dpC,ipHd - arrays of size `n`, where `n` denotes return value;
     * `dpC[i]` is reduced cost of variable with handle `ipHd[i]`.
     * \return number of variables.
     * \remark If `getReducedCosts()` is called with `dpC` and/or `ipHd` set to zero,
     * on return `dpC` and/or `ipHd` are pointers to internal `CLP` arrays.
     */
    int getReducedCosts(double* &dpC, int* &ipHd);

    /**
     * The function returns constraint _shadow prices_ (optimal values of dual variables).
     * \param[out] dpP,ipHd arrays of size `m`, where `m` is return value;
     *  `dpP[i]` is shadow price of constraint with handle `ipHd[i].
     * \return number of constraints.
     */
    int getShadowPrices(double* &dpP, int* &ipHd);

    /**
     * According to Farkas Lemma, a system of inequalities is infeasible
     * if can get the wrong inequality `0 < -1` by multiplying these inequalities by some numbers (multiplies)
     * --- of course,if we multiply an inequality by a negative number, we also have to change the sign of that inequality) ---
     * an then sum up the results. A collection of such multipliers is known as a _sertificate_ of inconsistency.
     *
     * `whyLpInfeasible()` returns a certificate of inconsistency.
     *
     * \param[out] m number of constraints;
     * \param[out] dpYctr,ipRowHd arrays of size `m`; if `dpYctr[i] > 0` (`dpYctr[i] < 0`),
     * then right (left) part of constraint with handle `ipRowHd[i]`
     * belongs to the set of contradicting constraints;
     * \param[out] n number of variables;
     * \param[out] dpYbd,ipColHd arrays of size `n`;
     * if `dpYbd[i] > 0` (`dpYbd[i] < 0`), then right (left) bound of variable
     * with handle `ipColHd[i]` belongs to the set of contradicting constraints.
     * \sa isLpInfeasible(), showWhyLpInfeasible().
     * \attention 1.  Call `whyLpInfeasible()` only if `isLpInfeasible()` returns `true`.
     *
     * 2. If the value of any input parameter is `0`,
     * then the output value of this parameter points to an internal MIPCL buffer.
     * Therefore, calling two function in sequence, keep in mind that the second call
     * may override the values returned by the first call.
     * This is not true for the handle-parameters, you can always
     * set initial values for these parameters to `0`.
     *
     * If you do not allocate memory to some input parameter and its value is not `0`,
     * then the behavior of the calling function is not predicted!!!
     */
    void whyLpInfeasible(int &m, int* &ipRowHd, double* &dpYctr,
                          int &n, int* &ipColHd, double* &dpYbd);

    /**
     * The function simply writes to the output stream the result returned by `whyLpInfeasible()`.
     * \param[in] out output file name.
     * \sa whyLpInfeasible(), showWhyLpInfeasible(const char* fileName).
     * \attention Call `showWhyLpInfeasible()` only if `isLpInfeasible()` returns `true`.
     */
    void showWhyLpInfeasible(std::ostream &out);

    /**
     * The function simply writes to the file the result returned by `whyLpInfeasible()`.
     * \param[in] fileName output file name.
     * \sa whyLpInfeasible(), showWhyLpInfeasible(std::ostream &out).
     * \attention Call `showWhyLpInfeasible()` only if `isLpInfeasible()` returns `true`.
     */
    void showWhyLpInfeasible(const char* fileName);

    /**
     * An LP is unbounded if the there exists a ray such that all its points are
     *  feasible solutions to the system of LP inequalities.
     * This procedure returns such a ray.
     * \param[out] n number of columns;
     * \param[out] dpX,dpRay arrays of size `n`;
     * `dpX + lambda * dpRay` is feasible point for all `lambda >= 0`;
     * \param[out] ipColHd arrays of size `n`; variable with handle `ipColHd[i]` has index `i`.
     * \return `true` if the problem is feasible but unbounded (`dpX` contains a feasible point),
     *   and `false` if the problem is either infeasible or unbounded; in the latter case,
     *   any feasible point has not been found, and, therefore, `dpX` contains nothing.
     *
     * \sa whyLpInfeasible().
     * \attention Call `whyLpUnbounded()` only if `isLpUnbounded()` returns `true`.
     */
    bool whyLpUnbounded(int &n, double* &dpX, double* &dpRay, int* &ipColHd);

    /**
     * The function simply writes to the file the result returned by `whyLpUnbounded()`.
     * \param[in] fileName output file name.
     * \sa whyLpUnbounded().
     * \attention Call `showWhyLpInfeasible()` only if `isLpUnbounded()` returns `true`.
     */
    void showWhyLpUnbounded(const char* fileName);

	/**
	 * The function writes LP solutions into the file.
	 * The user can overload this function to store solutions in an appropriate way.
	 * \param[in] fileName name of the file to store solutions; if `fileName=`0`,
	 * the solver makes up the file name by appending the extension ".sol" to the name of the problem being solved.
	 * \throws CFileException.
	 */
    virtual void printSolution(const char *fileName=0);

protected:
    /**
     * \return `-INF`.
     * \sa CMIP::getLastLowerBound().
     */
    virtual double getLastLowerBound()
    	{return -INF;}

    /**
     * Usually, the function is used in MIP applications.
     * \param[in] val objective value of some solution in original (not-scaled) variables.
     * \return objective value of that solution in scaled variables.
     * \sa getNotScaledObjVal().
     */
    double getScaledObjVal(double val) const;

    /**
     * Usually, the function is used in MIP applications.
     * \param[in] val objective value of some solution in scaled variable.
     * \return objective value of that solution in original (not-scaled) variables..
     * \sa getScaledObjVal().
     */
    double getNotScaledObjVal(double val) const;

    /**
     * The function is used privately in CMIP.
     * \return the number of times the basic matrix has been refactored
     *  during last call to `primeSimplex()` or `dualSimplexthis()`.
     */
    int getPartitionNum() const
    	{return m_iPartitionNum;}


///////////////////////////////////////////////////
//  Serialization
///////////////////////////////////////////////////
private:
    /**
     * The function stores into (or restores from) the stream `ar` the current basis.
     *
     * \param[in] ar reference to a stream;
     * \param[in] is_storing if `true`, the basis is written to the stream;
     *  otherwise, the basis is restored from the stream.
     */
    void serializeBasis(std::fstream& ar, bool is_storing);

    /**
     * The function stores into (or restores from) the stream `ar` the constraint matrix.
     *
     * \param[in] ar reference to a stream;
     * \param[in] is_storing if `true`, the matrix is written to the stream;
     *  otherwise, the matrix is restored from the stream.
     */
    virtual void serializeMatrix(std::fstream& ar, bool is_storing);

    /**
     * The function stores into (or restores from) the stream `ar`
     * the values of all tolerance parameters.
     *
     * \param[in] ar reference to a stream;
     * \param[in] is_storing if `true`, the parameter values are written to the stream;
     *  otherwise, the parameter values are restored from the stream.
     */
    virtual void serializeTolVars(std::fstream& ar, bool is_storing);

    /**
     * The function stores into (or restores from) the stream `ar`
     * the values of all parameters (flags)
     * that affect the performance of the simplex algorithms (primal and dual).
     *
     * \param[in] ar reference to a stream;
     * \param[in] is_storing if `true`, the parameter values are written to the stream;
     *  otherwise, the parameter values are restored from the stream.
     */
    virtual void serializeFlags(std::fstream& ar, bool is_storing);

protected:
    /**
     * The function stores into (or restores from) the stream `ar`
     * `CLP` objects (all its member storing permanent data).
     * Derived classes may overload this function to store additional information.
     * In such a case, call first `serialize()` of the base class.
     * \param[in] ar reference to a stream;
     * \param[in] is_storing if `true`, the object is written to the stream;
     *  otherwise, the object is restored from the stream.
     *  \sa CMIPP::serialize().
     */
    virtual void serialize(std::fstream& ar, bool is_storing);

/////////////////////////////////////////////////////
// Statistics
/////////////////////////////////////////////////////
    /**
     * \return number of iterates performed by `primeSimplex()` or `dualSimplex()`.
     */
    int getLPItNum() const
        {return m_iItNum;}

    /**
     * \return number of degenerate iterates performed by `primeSimplex()` or `dualSimplex()`.
     */
    int getDegItNum() const
        {return m_iDegItNum;}

    /**
     * The function is to compose a string describing the problem being solved.
     * \param[out] str memory for output string;
     * \return pointer to output string.
     * \remark If `str=0` when calling `getProbStatStr()`, on return `str` points to
     * an internal `CLP` array.
     */
    virtual char* getProbStatStr(char *str=0);

    /**
     * The function prints messages to the standard output stream.
     * This function should be overloaded in applications with graphic user interfaces (GUI).
     * \param[in] msg message string;
     * \param[in] level if set to 1) 0 -  info message; 2) 1 - warning message; 3) 2 - error message.
     * \remark Warning and error messages are also printed to LOG stream.
     */
    virtual void infoMessage(const char* msg, int level=0);

    /**
     * The function prints into the standard output stream a message string
     * which describes the current state of the solution process when
     * running `primeSimplex()` or `dualSimplex()`.
     * When developing an application with a GUI interface,
     * the user may wish to overload this function.
     * \param[in] method string describing simplex method that is currently running;
     * \param[in] time string representation of time elapsed since application start;
     * \param[in] itNum number of simplex-iterates;
     * \param[in] degItNum number of degenerate simplex-iterates;
     * \param[in] objVal current objective value.
     */
    virtual void lpInfo(const char *method, const char *time, int itNum, int degItNum, double objVal);

private:
    /**
     * The function prepares parameter for calling `lpInfo()` and then calls it.
     */
    void __lpInfo();

public:

#ifdef __PYTHON_MODULE_
    /**
     * To make communication between __Python__ and __MIPCL__
     * easier, `setSolution()` lists components of the prime and dual solutions
     * in order of their handles (initial indices).
     */
    void setSolution();

    /**
     *  When the problem has been proven to be infeasible,
     * this procedure is called to set an inconsistency certificate.
     *
     */
    void setInconsistCertificate();

    /**
     * This function must be called only after calling `setSolution()`.
     * \param[in] ind index (handle) of variable.
     * \return value of variable indexed by `ind`.
     */
    double getOptVarValue(int ind) const;

    /**
     * This function must be called only after calling `setSolution()`.
     * \param[in] ind constraint index.
     * \return shadow price of constraint indexed by `ind`.
     */
    double getShadowPrice(int ind) const;

#endif


// Debugging tools

	/**
	 * The function is useful for debugging. It prints a given row of the matrix to the standard error stream.
	 * \param[in] i row index;
	 * \param[in] scaled if `true`, row of scaled matrix is printed; otherwise, not-scaled;
	 * \param[in] varValues if `true`, current values of variables that occur in row are also printed.
	 */
   void printRow(int i, bool scaled=true, bool varValues=false);

   /**
	* The function is used for debugging. It prints the constraint `sum(i=0,..,sz-1) dpVal[i]*x(ipCol[i]) <= (>=) b`
	*  given given by its argument to the standard error stream.
	* \param[in] sz size of arrays `dpVal` and `ipCol`;
	* \param[in] dpVal array of size `sz` storing coefficients;
	* \param[in] ipCol array of size `sz` storing indices of variables;
	* \param[in] b depending on value of `side`, right-hand-side or left-hand-side of inequality to be printed;
	* \param[in] side if `true`, the inequality sign is `<=`; otherwise, the sign is `>=`.
	*/
   void printCtr(int sz, double* dpVal, int* ipCol, double b,  bool side=true);

   /**
	* The function is useful for debugging. It prints a given column of the matrix to the standard error stream.
	* \param[in] j column index.
	*/
	void printColumn(int j);

	/**
	 * The function is useful for debugging. Prints a the matrix to the file.
	 * \param[in] fileName file name;
	 * \param[in] scaled if `true`, scaled matrix is printed; otherwise, not-scaled;
	 * \throws CFileException.
	 */
	void printMatrix(const char* fileName, bool scaled=false);

/////////////////////////////////////////////////////
private:
	/**
	 * The procedure performs the pivot iteration when a row is substituted for a row.
	 * It does the same as `substituteRowForRow()` does, but without affecting LU-partition
	 * of the basic matrix.
	 */
	int STRBR_substituteRowForRow(int r, int s, bool side);

	/**
	 * The procedure performs the pivot iteration when a row is substituted for a column.
	 * It does the same as `substituteRowForColumn()` does, but without affecting LU-partition
	 * of the basic matrix.
	 */
	int STRBR_substituteRowForColumn(int r, int s, bool side);

	/**
	 * The procedure performs the pivot iteration when a column is substituted for a row.
	 * It does the same as `substituteColumnForRow()`, but without affecting LU-partition
	 * of the basic matrix.
	 */
	int STRBR_substituteColumnForRow(int& r, int s, bool side);

	/**
	 * The procedure performs the pivot iteration when a column is substituted for a Column.
	 * It does the same as `substituteColumnForColumn()`, but without affecting LU-partition
	 * of the basic matrix.
	 */
	int STRBR_substituteColumnForColumn(int r, int s, bool side);

	/**
	 * This procedure is a modification of `doPivot()` to speed up `CMIP::strongBranching()`.
	 */

	bool STRBR_doPivot(int s, int r, bool side);

	/**
	 * This procedure is a modification of `computeUd()` to speed up `CMIP::strongBranching()`.
	 */
	void STRBR_computeUd(int s, double& dMaxVal);

	/**
	 * This procedure is a modification of `dualRatioTest()` to speed up `CMIP::strongBranching()`.
	 */
	int STRBR_dualRatioTest(int s, bool side);

	/**
	 * This procedure is a modification of `dualGetPivotVal()` to speed up `CMIP::strongBranching()`.
	 */
	double STRBR_dualGetPivotVal(int s, int t, double* U);

	/**
	 * This procedure is a modification of `dualSeekPivot()` to speed up `CMIP::strongBranching()`.
	 */
	void STRBR_dualSeekPivot(int& s, int& t, bool& side);
protected:
	/**
	 * This procedure implements a very restricted variation of the dual simplex algorithm
	 * in order to fulfill just a few iterates (without re-factoring the basic matrix).
	 * Such a behavior is needed to efficiently implement a strong branching procedure.
	 *
	 * \param[in] lowerBound lower bound on optimal objective value;
	 * \param[in] maxItNum maximum number of iterates to be accomplished.
	 * \return 1: in case of success;
	 *     2: exceeded maximum number of iterations;
	 *     3: exceeded given lower bound on optimal objective value;
	 *     6: any numeric instability.
	 */
	int STRBR_estimateObjDecrease(double lowerBound, int maxItNum);

///////////////// Norm constraints ||x|| <= t

public:
	/**
	 * Using quadratic constraints in MIPs is not common.
	 * Therefore, __MIPCL__ treats norm-cone constraints of the form \f$\|x\|\le t\f$ as a supplementary feature.
	 * Despite of apparent simplicity of such norm-constraints, we can use them to model
	 * many restrictions of practical importance.
	 * To use such restrictions in your __MIPCL__ application,
	 * call `allowNormCtrs()` before `closeMatrix()`.
	 * In most case (if not in all) you must also disable the preprocessing,
	 * the application of which to an incomplete formulation will almost certainly result in an error.
	 *
	 * \param[in] maxCtrNum maximum number of norm constraints;
	 * \param[in] avCtrSize average size of norm constraints to be added.
     * \throws CMemoryException lack of memory.
	 */
    void allowNormCtrs(int maxCtrNum, int avCtrSize);

    /**
     * This function adds a norm constraint.
     *
     * \param[in] t right-hand-side variable;
     * \param[in] sz number of variables involved in left-hand-side vector;
     * \param[in] ipVars list (of size `sz`) of variables involved in left-hand-side vector;
     * \param[in] tol tolerance value.
     * \throws CMemoryException lack of memory.
     * \sa allowNormCtrs(), separate()
     * \attention do not call `addNormCtr()` without calling `allowNormCtrs()`.
     */
	int addNormCtr(int t, int sz, int *ipVars, double tol=0.0001);

#ifdef __PYTHON_MODULE_
    /**
     * When __MIPCL__ is called from __Python__, any norm constraint is posed in two steps:
     * first `startNormCtr()` does almost all job except storing components of the left-hand-side vector;
     * then `addNormVar()` is repeatedly called to add the components of the right-hand-side vector.
     *
     * \param[in] t right-hand-side variable;
     * \param[in] sz size of left-hand-side vector;
     * \param[in] tol tolerance value.
     * \return position (in the norm constraint buffer) of the first component of the left-hand-side vector.
     */
	int startNormCtr(int t, int sz, double tol);

    /**
     * This function must be called a required number of times to complete posing
     * a norm constraint started by `startNormCtr()`.
     * \param[in] k position (in the norm constraint buffer) of added variable;
     * \param[in] var index of added variable.
     * \return position (in the norm constraint buffer) of the next component of the left-hand-side vector.
     *
     * \sa startNormCtr()
     */
	int addNormVar(int k, int var);
#endif

};

#endif // #ifndef __LP__H
