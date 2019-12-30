///////////////////////////////////////////////////////////////
/**
 * \file cmip.h Interface for `CMIP` class
 * |  __Author__  | N.N. Pisaruk                              |
 * |-------------:|:------------------------------------------|
 * |  __e-mail__  | nicolaipisaruk@gmail.com                  |
 * | __home page__| wwww.mipcl-cpp.appspot.com                |
 *
 * \copyright __2019 Nicolai N. Pisaruk__
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

#ifndef __CMIP__H
#define __CMIP__H

#include "lp.h"
#include "thread.h"

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
#endif // MIP_EXPORTS
#endif // MIP_API
//typedef __int64 __LONG; ///< For Visual C++, `long` means 32 bit integer; so we redefine it.
#else
#ifndef MIP_API
#define MIP_API
#endif // MIP_API
//typedef long long __LONG; ///< For GNU C++, `long` means 64 bit integer.
#endif // _WINDOWS

#ifndef __PSEUDOCOSTS_
#define __PSEUDOCOSTS_ ///< Switches on computing pseudocosts.
#endif

class CPool;
class CNode;
class CTree;
class CRecord;
class CImpl;
class CAUT;

///  This structure specifies how cuts of a particular type be generated.
/**
 * Any record of type `cutProps` affects the way the cuts of a particular type are generated.
 */
struct cutProps {
    double maxSlack; ///< Only constraints which slack is less than `maxSlack` are used for generating cuts of this type.
    double cutTol; ///< Tolerance value for cuts of this type.
    int maxCutRounds; ///< The solver calls at most `maxCutRounds` times cut generating procedure for cuts of this type.
    int	maxCutsPerRound; ///< At most `m_ipMaxCutsPerRound[i]` cuts are generated during each call of cut generating procedure for cuts of this type.
    int	maxCutSize; ///< For cuts of this type,  all cuts of size greater than `maxCutSize` are rejected.
    int cutNodes; ///< Number of nodes for which cuts of this type are generated.
    int cutNodeHeight; ///< Minimum height of a node for which cuts of this type are generated.
  	int maxCutPrc;  ///< For cuts of this type, only cuts of size not greater than `(maxCutPrc*m_iN)/100` are generated.
  	int maxCutExp; ///< Cuts with biggest coefficient exponent greater than `maxCutExp` are rejected.
};

 /// Class `CMIP` has been designed for solving Mixed Integer Programs (MIPs)
 /**
 * \f{align*}{
 *     & c^Tx \to \max,\\
 *     & b_1 \le Ax \le b_2,\\
 *     & l \le  x \le u,\\
 *     & x_i\in\mathbb{Z}\quad\forall\; i\in I,
 * \f}
 * by the branch-and-cut or branch-and-price method.
 *
 * __Examples__ of how to use `CMIP` class:
 *    - \ref primer;
 *    - \ref fcnf;
 *    - \ref infLP;
 *    - \ref subgraph;
 *    - \ref MSched;
 *    - \ref gap;
 *    - \ref tsp
 */
class MIP_API CMIP: public CLP
{
	friend class CTree;
	friend class CImpl;
	friend class CAUT;

public:
	/// MIP types of constraints.
	/**
	 * The type of a constraint is the bitwise OR of the members of two enumerations CLP::enCtrType and CMIP::enCtrType.
	 * A constraint is of a particular type if it can be transformed to that type by complementing some binary variables.
	 */
	enum enCtrType {
		CTR_LOCAL            = 0x00000040, ///< The constraint is _local_ for node being processed.
		CTR_INT_VARS         = 0x00000080, ///< All variables and coefficients are integral in this constraint.
		CTR_INT_COEFF        = 0x00000100, ///< All coefficients are integral in this constraint.
		CTR_INT              = CTR_INT_VARS|CTR_INT_COEFF, ///< The constraint with all variables and coefficients being integral.
		CTR_VAR              = 0x00000200, ///< _Variable bound_, i.e. \f$z_i \le (\ge) a_jx_j +b,\quad z_i\in\mathbb{R},\; x_j\in \{0,1\}\f$.
		CTR_BINPACK          = 0x00000400, ///< _Bin packing_, i.e., \f$\sum_{i\in I} a_i x_i + a_j x_j \le a_j,\quad j\not\in I, x_j\in\{0,1\},\; x_i\in\{0,1\}\;\forall\: i\in I\f$.
		CTR_KNAPSACK         = 0x00000800, ///< \f$\sum_{i\in I} a_i x_i \le b,\quad x_i\in\{0,1\}\;\forall\: i\in I\f$.
		CTR_MX_KNAPSACK      = 0x00001000, ///< \f$\sum_{i\in I} a_i x_i + \sum_{j\in J} q_j z_j \le b,\quad x_i\in\{0,1\}\;\forall\: i\in I\f$.
		CTR_PACKING          = 0x00002000, ///< \f$\sum_{i\in I} x_i \le 1,\quad x_i\in\{0,1\}\;\forall\: i\in I\f$.
		CTR_INV_KNAPSACK     = 0x00004000, ///< \f$\sum_{i\in I} x_i \le b,\quad x_i\in \{0,1\}\;\forall\: i\in I\f$.
		CTR_COVERING         = CTR_INV_KNAPSACK, ///< \f$\sum_{i\in I} x_i \ge 1,\quad x_i\in\{0,1\}\;\forall\: i\in I\f$.
		CTR_CARDINALITY      = 0x00008000, ///< \f$\sum_{i\in I} x_i \le k,\quad x_i\in\{0,1\}\;\forall\: i\in I\f$.
		CTR_GUB              = CTR_PACKING|CTR_CARDINALITY, ///< \f$\sum_{i\in I} x_i \le 1,\quad x_i\in\{0,1\}\;\forall\: i\in I\f$.
		/**
		 * A __SOS1__-constraint is a __GUB__-constraint (of type `CTR_GUB`) such that
		 *  a) all its variables are present in exactly the same other constraints,
		 *  and b) in any of those constraints, if we list the coefficients of SOS1-variables in the order
		 *  these variables are present in the base SOS1-constraint, we get an increasing or decreasing sequence.
		 */
		CTR_SOS1             = 0x00010000, ///< \f$\sum_{i\in I} z_i \le (=) 1,\quad z_i \ge 0\;\forall\: i\in I \text{ and only one } z_i \text{ takes non-zero value} \f$.
		CTR_SOS2             = 0x00020000, ///< \f$\sum_{i\in I} z_i \le (=) 1,\quad z_i \ge 0\;\forall\: i\in I,\text{ and if } z_i> 0,\: z_j>0, \text{ then } |i-j| = 1\f$.
		CTR_01FLOW           = 0x00040000, ///< \f$\sum_{i\in I} a_i z_i \le (=) b,\quad 0\le z_i \le \alpha_ix_i +\beta_i \text{ with } x_i\in \{0,1\},\; i\in I\f$.
		CTR_WITH_VAR_BOUNDS  = 0x00080000, ///< \f$\sum_{i\in I} a_i z_i \le (=) b,\quad 0\le z_i \le a_ix_i +b_i \text{ with } x_i\in\{0,1\},\; i\in J\subset I\f$.
		CTR_WITH_UNIQUE      = 0x00100000, ///< Any constraint with only one non-binary variable.
		CTR_PARITY           = 0x00200000, ///< Any integer constraint with more than two odd coefficients.
		CTR_MX_01            = 0x00400000, ///< The constraint has binary and non-binary variables.
		CTR_MX_INT           = 0x00800000, ///< The constraint has integer and non-integer variables.
		CTR_GEN              = 0x00800000, ///< All other constraints (currently not used).
		CTR_NOT_INV          = 0x01000000, ///< Not all absolute coefficient values are ones.
		CTR_BRANCHING_INV    = 0x02000000, ///< Invariant knapsack used in balancing branching.
		CTR_WITH_DEP_BINS    = 0x04000000, ///< There are implication relations, \f$x_i=a\; \Rightarrow\; x_j=b,\; (a,b \in \{0,1\}\f$), between their binary variables.
		CTR_IN_POOL          = 0x08000000, ///< The constraint is stored in pool.
		CTR_LB_UNBOUNDED     = 0x10000000, ///< The flag is _privately_ used in preprocessing subroutines.
    	CTR_UB_UNBOUNDED     = 0x20000000, ///<The flag is _privately_ used in preprocessing subroutines.
//		CTR_DISJUNCTION      = 0x40000000, ///< The constraint contains binary variables setting which to one of its bounds (0 or 1) makes this constraint _free_.
		CTR_OBJ              = 0x40000000, ///< The constraint represents the objective, its lhs and rhs are lower and upper bounds on the optimal objective value.
// flag `0x80000000` is used in `CLP` for `CTR_STRONG_CUT`
		CTR_BINARY           = CTR_BINPACK | CTR_KNAPSACK | CTR_INV_KNAPSACK | CTR_PACKING | CTR_COVERING | CTR_CARDINALITY, ///< All variables are binaries.
		CTR_WITH_INT_VARS    = CTR_INT_VARS | CTR_MX_01 | CTR_MX_INT | CTR_MX_KNAPSACK, ///The constraint includes integer variables.
		CTR_MIR              = CTR_MX_INT|CTR_MX_01|CTR_MX_KNAPSACK|CTR_WITH_VAR_BOUNDS|CTR_01FLOW|CTR_KNAPSACK|CTR_BINPACK,//|CTR_WITH_UNIQUE, ///< Mixed Integer Rounding (MIR) constraint
		CTR_FULLY_CLASSIFIED = CTR_VAR|CTR_BINPACK|CTR_KNAPSACK |CTR_INV_KNAPSACK|CTR_PACKING|CTR_COVERING|CTR_PARITY|CTR_SOS1|CTR_SOS2, ///< Mask of all fully classified constraints
		CTR_CLASSIFIED = CTR_MIR|CTR_01FLOW|CTR_INT|CTR_VAR|CTR_INV_KNAPSACK |CTR_PACKING|CTR_COVERING|CTR_CARDINALITY|CTR_PARITY|CTR_BINPACK|CTR_BRANCHING_INV|CTR_WITH_DEP_BINS|CTR_NOT_INV ///< Mask of automatically assigning flags
	};

	/// MIP types of variables.
	/**
	 * The type of a variable is the bitwise OR of the members of two enumerations, `CLP::enVarType` and `CMIP::enVarType`.
	 */
	enum enVarType {
//		VAR_INT           = 0x00001000, ///< Integral variable.
		VAR_BIN           = 0x00002000, ///< Binary variable.
		VAR_4_INT         = 0x00004000, ///< Both bounds, lower and upper, are in [0,15].
		VAR_8_INT         = 0x00008000, ///< Both bounds, lower and upper,  are in [0,255].
		VAR_IN_VAR_LB     = 0x00010000, ///< The non-binary variable occurs in some lower variable bound constraint.
		VAR_0_IN_VAR_CTR  = VAR_IN_VAR_LB, ///< If binary variable is set to `0`, then some variable bound constraint becomes a restriction.
		VAR_IN_VAR_UB     = 0x00020000, ///< The variable occurs in some upper variable bound constraint.
		VAR_1_IN_VAR_CTR  = VAR_IN_VAR_UB, ///< If the binary variable is set to `1`, then some variable bound constraint becomes a restriction.
		VAR_IN_VAR_CTR    = VAR_IN_VAR_UB | VAR_IN_VAR_LB, ///< The variable occurs in either lower, or upper variable bound constraint.  \sa CTR_VAR.
		VAR_IN_GUB        = 0x00040000, ///< The variable occurs in _GUB_ constraints, which are cardinality packing constraints. \sa CTR_PACKING and CTR_CARDINALITY.
		VAR_IN_PACKING    = 0x00080000, ///< The variable occurs in _packing_ constraints. \sa CTR_PACKING.
		VAR_UNIQUE        = 0x00100000, ///< The non-binary variable that occurs in constraints which all other variables are binary.
		VAR_MON_UP        = 0x00200000, ///< The variable can be increased up to its upper bound without violating feasibility.
		VAR_MON_DOWN      = 0x00400000, ///< The variable can be decreased up to its lower bound without violating feasibility.
		VAR_MONOTONE      = VAR_MON_UP | VAR_MON_DOWN, ///< The variable is monotone down or up.
		VAR_BINPACK       = 0x00000000, ///< The variable occurs in _bin-packing_ constraint (currently is not used). \sa CTR_BINPACK.
		VAR_IN_POOL       = 0x00800000, ///< The column corresponding to this variable is stored in the pool.
		VAR_SOS           = 0x08000000, ///< The variable occurs in a _SOS1_ or _SOS2_ constraint. \sa CTR_SOS1, CTR_SOS2.
		VAR_CLASSIFIED    = VAR_IN_VAR_CTR|VAR_IN_PACKING|VAR_IN_GUB|VAR_MONOTONE|VAR_UNIQUE|VAR_IN_POOL|VAR_BINPACK, ///< mask of automatically assigned flags.
		VAR_PRI_MIN       = -50, ///< The priority of any variable is an integer from `VAR_PRI_MIN` to `VAR_PRI_MAX`.
		VAR_PRI_MAX       = 50 ///< The priority of any variable is an integer from `VAR_MIN_PRI` to `VAR_MAX_PRI`.
	};

	/// Branching rules.
	/**
	 * For each integer-valued variable \f$j\f$ with fractional value \f$x_j\f$, we estimate from below
	 *    - \f$\lambda^{\mathrm{down}}_j\f$: decrease  of the objective if we branch by this variable down;
	 *    - \f$\lambda^{\mathrm{up}}_j\f$: decrease  of the objective if we branch by this variable up.
	 *
	 * Then we compute the score \f$\lambda_j\f$ as a function of \f$\lambda^{\mathrm{down}}_j\text{ and } \lambda^{\mathrm{up}}_j\f$,
	 * say, \f$\lambda_j=\lambda^{\mathrm{down}}_j \cdot \lambda^{\mathrm{up}}_j\f$, or \f$\lambda_j=\min{\lambda^{\mathrm{down}}_j,\lambda^{\mathrm{up}}_j}\f$.
	 * A variable \f$j^*\f$ with maximum score is selected for doing branching.
	 */
	enum enBranchRule {
		/// Strong branching (default strategy).
		/**
		 * A candidate set of fractional variables is selected;
		 * for every candidate \f$x_j\f$, for each of 2 branches a fixed number (approximately 10)
		 * of dual iterates are performed and the decrease of the objective value is fixed;
		 * the product of these two decreases is the score.
		 */
		STRONG_BR,
		/// Score branching.
		/**
		 * Let\f$\lambda^{\mathrm{down}}_j\f$ (resp., \f$\lambda^{\mathrm{up}}_j\f$) be the decrease of the objective
		 * after performing one dual iterate; the minimum of these two decreases is the score.
		 */
		MAX_SCORE
	};

	/// Rounding rules.
	/** Specifies all possible values of m_eRoundType, which determines
	 * how non-integral solutions are rounded.
	 * if m_eRoundType != ROUND_NONE, at each node, an optimal solution is rounding
	 * in hope to get a feasible solution.
	 */
	enum enRoundType {
		ROUND_NONE, ///< None of rounding procedures is called.
		ROUND_OFF, ///< Round off values of all integral variables.
		ROUND_DOWN, ///< Round down values of all integral variables.
		ROUND_UP, ///< Round up values of all integral variables.
		ROUND_USER ///< User defined procedure is called for rounding integer variables.
	};

	/// Cut types.
	enum enCutType {
		CUT_TYPE_NUM           =13, ///< Number of cut types (classes).
		_CLICK                  =0, ///< Click cuts.
		_KNAPSACK               =1, ///< Cuts generated for knapsack constraints.
		_MX_KNAPSACK            =2, ///< Generalized cover cuts generated for mixed-knapsacks.
		 _MIR                   =3, ///< Mixed-Integer-Rounding (MIR) cuts.
		_FLOW_COVER             =4, ///< Flow-cover cuts generated by 0,1-flow constraints.
		_SPARSE_MOD2            =5, ///< Sparse (0,1/2)-Chvatal-Gomory (mod2) cuts.
		_DENSE_MOD2             =6, ///< Dense (0,1/2)-Chvatal-Gomory (mod2) cuts.
		_SPARSE_GOMORY          =7, ///< Sparse mixed-integer Gomory cuts.
		_DENSE_GOMORY           =8, ///< Dense mixed-integer Gomory cuts.
		_PARITY                 =9, ///< Parity cuts.
		_VAR_BOUND             =10, ///< Implied variable bounds.
		_SIMPLE_DJ             =11, ///< Simple disjunctions.
		_USER_DEF              =12, ///< Cuts generated in user procedures `saparate()`, genCut1()`, and genCut2()`.
		CUTS_TO_POOL	       = 0  ///< The mask of cuts send to the pool.
	};

	/// Extension of `enRowColGenRule`.
    enum enRowColGenRule2 {
    	GEN1_PROC	    = 0x8, ///< If `genCut1()` has been overloaded.
    	GEN2_PROC	    = 0x10, ///< If `genCut2()` has been overloaded.
    	STOP_AUTO_CUTS  = 0x20, ///< If the flag is set, generation of  auto cuts is stopped.
    	WITH_AUTO_CUTS  = 0x40, ///< `autoCuts()` was not called at the previous call to `generateCuts()`.
    	NO_SOLVER_DECISIONS = 0x80 ///< If user resets any parameter affecting cut generation, _MIPCL_ does not change standard settings of all parameter affecting cut generation.
    };


/////////////////////////////////////////////////////////////
private:
	static const char m_cpCutName[CUT_TYPE_NUM][16]; ///< names of cuts

    enSepRule m_eNodeSepRule; ///< Separation rule used when solving not-root node LPs.
    enPricingRule m_eNodePricingRule; ///<< Pricing rule used when solving not-root node LPs..

    int m_iVarLbNum;///< Number of variable lower bounds.
    int m_iVarUbNum; ///< Number of variable upper bounds.
    int m_iIntNum;  ///< Number of integer variables.
    int m_iBinNum; ///< Number of binary variables.
    int m_iFracVarNum; ///< Number of fractional components in the current basic solution.
    int m_iParityNum; ///< Number of parity constraints.

    int m_iProbingDepth; ///< Depth of probing search, default value is `PROBING_DEPTH`.
    int m_iProbingDepthAtNodes; ///< Default value `PROBING_DEPTH_AT_NODES`.

    CImpl *m_pImpl; ///< Class that supports logical implications.
////////////////////////////////////////////////////
    /**
     * Array of size `2*m_iN`,
     * `m_dpD0[j<<1]` and `m_dpD0[(j<<1)+1]` are lower and upper bound of variable `j`.
     */
	double *m_dpD0; ///< Bounds of variables at the root node.

//////// Branching ///////////////////////////////////////////////////
    enBranchRule m_eBranchRule; ///< Currently used branching rule.

    int m_iBrCol; ///< Variable selected for branching.
    double m_dBrVarVal; ///< Value of variable selected for branching.

    int m_iBrRow; ///< _GUB_-constraint selected for branching.
	int m_iFxNum; ///< Number of variables fixed in `strongBranching`.
	int *m_ipFxVar; ///< Array of size 4*STRONG_BR_LOW_NODE_CAND_NUM used in `GUBbranching` and `strongBranching`.

#ifdef __PSEUDOCOSTS_
    /**
     * `m_dpBrDec[j<<1]` (resp., `m_dpBrDec[(j<<1)+1]`) is total relative objective decrease
     * when branching by variable `j` rounding it down (resp., up)
     */
	double *m_dpBrDec;

	/**
	 * An array of size `3*m_iN`;
     * For `j=0,...,m_iN-1`,
     *     - `m_ipBrNum[3*j]`: two lower (upper) bytes store the number of consistent (node LP remains consistent) branchings when variable `j` was rounded down (resp., up);
     *     - `m_ipBrNum[3*j+1]`: two lower (upper) bytes store the number of successful branchings (node LP remains consistent and objective has increased)
     *     - `m_ipBrNum[3*j+2]`: two lower (upper) bytes store the number inconsistent branchings (node LP became inconsistent) when variable `j` was rounded down (resp., up);
     *       (node LPs were inconsistent or the objective increased when variable `j` was rounded down (resp., up).
	 */
	unsigned *m_ipBrNum;
#endif

////////  Diving   ///////////////////////////////////////////////////
	int m_iDivingCalls; ///< Number of times `dive()` has been called so far.
	int m_iToughDivingCalls; ///< If positive, number of times `dive()` exceeds time limit; negative value indicates `dive()` is called rarely.
//////////////////////////////////////////////////////////////////////////
//                        POOL to store generated cuts                  //
//////////////////////////////////////////////////////////////////////////
    CPool* m_pPool; ///< Pointer to the pool that stores cuts.
    bool m_bPool; ///< If set to _false_, memory for the pool is not allocated.
	double m_dPoolCtrTol; ///< Maximum violation allowed for constraints in the pool.
	double m_dCtrCutTol; ///< Maximum violation allowed for cuts that are currently in the constraint matrix.

// ////////////////////////////////////////////////////////////////////////
//                           Branch And Cut TREE                         //
// ////////////////////////////////////////////////////////////////////////
    CTree* m_pTree; ///< Pointer to the root node of the search tree.
    int m_iParentNode; ///< Index of the parent node.
    int m_iNode; ///< Index of the currently processed node.
    int m_iHeight; ///< Height of the currently processed node.

    int m_iObjSize; ///< Number of non-zero coefficients in the objective.
	bool m_bObjInt; ///< `true` if optimal objective value is integral.
	/**
	 * If `0 <= m_iObjRow < m_iM`, then row `m_iObjRow` represents the objective,
	 * and LHS and RHS of constrains `m_iObjRow` are lower and upper bounds for the optimal objective value.
	 * If `m_iObjRow >= SHIFT`, then the objective is to maximize variable `m_iObjRow-SHIFT`.
	 * If `m_iObjRow=-1`, then the objective is not represented by any constraint.
	 */
	int m_iObjRow;
	/**
	 * If `m_dObjVal < m_dMipObjVal + m_dAbsObjTol`, the node is pruned.
     */
    double m_dAbsObjTol; ///< Tolerance value for MIP objective
    /**
     * If the difference between some value and its nearest integer is less than `m_dIntTol`,
     * then this value is treated as an integer.
     */
    double m_dIntTol;

	CRecord *m_pRecord; ///< Pointer to structure that stores the best MIP solution found so far.
    double m_dLoBound; ///< equals `m_pRecord->m_dLoBound` when the currently processed node (subproblem) was restored.

    char *m_strSolFileName; ///< Pointer to string with file name for storing intermediate solutions; if `m_strSolFileName=0`, intermediate solutions are not stored.

    enRoundType m_eRoundType; ///< Specifies how non-integer solutions are rounded.

    __LONG m_lTimeToStop; ///< `BranchAndCut()` procedure will stop running at time `m_lTimeToStop` (given in seconds since epoch); `m_lTimeToStop==0l` means no time limitation.
    double m_dDualGap; ///< Branch-and-cut stops if the dual gap is less than `m_dDualGap`.

///// STATISTICS ////////////////////////////////////////
    unsigned int m_iBranchAndCutNodes; ///< Number of branch-and-cut nodes processed.
    unsigned int m_iDifficultNodes; ///< Number of nodes where numerical difficulties were met.

///// CUT GENERATION ///////////////////////////////////
	int m_iM1; ///< Number of rows when a new round of cut generation starts.
    int m_iLastKnapsack0; ///< Index of the last knapsack constraint in the root-node matrix.
    int m_iLastKnapsack; ///< Starts linked list of knapsack constraints.
    int m_iLast01Flow0; ///< Index of the last flow constraint in the root-node matrix.
    int m_iLast01Flow; ///< Starts linked list of variable bound constraints in the currently processed matrix.
    int m_iLastVarBound0; ///< Number of variable bound constraints in the root matrix.
    int m_iLastVarBound;  ///< Starts linked list of variable bounds in the currently processed matrix.
    int m_iLastGUB;  ///< Starts linked list of GUB-constraints in the currently processed matrix.
    int m_iLastGUB0; ///< Index of the last GUB-constraint in the root-node matrix.
    int m_iLastPacking0; ///< Index of the last packing constraint in the root-node matrix.
    int m_iLastPacking;  ///< Starts linked list of packing constraints.
    int m_iLastMixedKnapsack0; ///< Index of the last mixed-knapsack constraint in the root-node matrix.
    int m_iLastMixedKnapsack;  ///< Starts linked list of mixed-knapsack constraints.
    int m_iLastParity0; ///< Index of the last parity constraint in the root-node matrix.
    int m_iLastParity;  ///< Starts linked list of parity constraints.
//    int m_iLastSimpleDJ0; ///< Index of the last simple-disjunction constraint in the root-node matrix.
//    int m_iLastSimpleDJ;  ///< Starts linked list of simple-disjunction constraints.
//    int m_iLastBranchingInv0; ///< Index in the root-node matrix of the last invariant knapsacks, which all variables are in variable bound.
//    int m_iLastBranchingInv;  ///< Starts linked list of invariant knapsacks, which all variables are in variable bound; such invariant knapsacks are used in balanced branching.

    int* m_ipCtrLink; ///< Linked lists of constraint classes.

// constraint statistics
    int m_iAutoCutRound; ///< Current round of generating cuts.
    int m_iMaxAutoCutRounds; ///< Maximum number of rounds (iterates) during one call to `generateCuts()` for currently processed node.
    int m_iMinAutoCutRounds; ///< Minimum number of rounds (iterates) during one call to `generateCuts()` for currently processed node.
    int m_iMaxAutoCutRoundsAtNodes; ///< Number of rounds for generating cuts per node.
    int m_iMaxAutoCutRoundsAtRoot; ///< Number of rounds for generating cuts at root node.
    int m_iMaxCutSize; ///< Maximum cut size; it is set to either `m_iMaxCutSizeAtRoot` or `m_iMaxCutSizeAtNodes`.
    int m_iMaxCutSizeAtRoot; ///< Maximum cut size at root nodes.
    int m_iMaxCutSizeAtNodes; ///< Maximum cut size at non-root nodes.

    int m_iUserCutPattern; ///< Determines cuts which patterns cannot be changed by the solver.

    int m_iAutoCutNodes; ///< Cuts are automatically generated for first `m_iAutoCutNodes` processed nodes.
    int m_iAutoCutNodeHeight; ///< Cuts are automatically generated for all nodes of height less or equal than `m_iAutoCutNodeHeight`.

    /**
     * At the root node, if after one-round-cut generation the objective decreases less than
 	 * by `m_dCutRoundPrcAtRoot` times its value before cuts were added,
 	 * then the cut generation procedure terminates.
     */
	double m_dCutRoundPrcAtRoot;

    /**
     * At non-root nodes, if after one-round-cut generation the objective decreases less than
 	 * by `m_dCutRoundPrcAtRoot` times its value before cuts were added,
 	 * then the cut generation procedure terminates.
     */
	double m_dCutRoundPrc;

	int m_iMinAutoCutRoundsAtRoot; ///< Minimum round to generate cuts at root node.
	int m_iMinAutoCutRoundsAtNodes;  ///< Minimum round to generate cuts at not-root nodes.
	int m_iMinDenseCutSize; ///< Cuts of size not greater than `m_iMinDenseCutSize` are never rejected as dense ones.

	bool m_bCutsForThisNode; ///< If `true`, at least one cut has been generated for node being processed.

	unsigned m_uiCutsToPool; ///< See `sendCutsToPool()` and `areCutsSentToPool()` for the explanation of using this field.

    int m_ipCutRound[CUT_TYPE_NUM]; ///< `m_ipCutRound[i]` is maximum number of times `generateCuts()` call procedure that generates cuts of type `i`.
	cutProps m_pCutPropsBuf[CUT_TYPE_NUM-1]; ///< Memory buffer for `m_pCutProps`.
	cutProps *m_pCutProps; ///< The record `m_pCutProps[i]` describes properties of cuts of type `i`;

// Statistics of using cuts
  	int m_ipCurrentRoundCutNum[CUT_TYPE_NUM-1]; ///< `m_ipCurrentRoundCutNum[i]` is number of cuts of type `i` produced at the current cut-generation round.
  	int m_ipTotalCutNumBuf[CUT_TYPE_NUM<<1]; ///< Memory buffer for `m_ipTotalCutNum`.
  	int *m_ipTotalCutNum; ///< `m_ipTotalCutNum[i<<1]` is total (generated) number and `m_ipTotalCutNum[(i<<1)+1]` is used (present in the basis just after re-optimization) number of cuts of type `i`.

    int m_iCurrentMxTemplate; ///< used internally to smoothly use templates when generating  MIR cuts.

////// SOS constraints ////////////////////////////
    int m_iSOS1num; ///< Number of SOS1 constraints.
    int m_iSOS2num; ///< Number of SOS2 constraints.
    int m_iLastSOS1;  ///< Starts linked list of SOS1-constraints in the currently processed matrix.
    int m_iLastSOS10; ///< Index of the last SOS1-constraint in the root-node matrix.
    int m_iLastSOS2; ///< Starts list of SOS2 constraints.

////// SYMMETRY Breaking
    CAUT *m_pAut; ///< Automorphism group.
//////////////////////////////////////////////////
#ifndef __ONE_THREAD_
// Multithreading
    int m_iCoreNum; ///< Number of physical cores.
    int m_iThreadNum; ///< Thread number used when solving MIPs,; by default, `m_iThreadNum=m_iCoreNum`.
	_RWLOCK *m_rwStatLock; ///< Locks changing statistics attributes: `m_ipTotalCutNum`.
	_MUTEX *m_pMemMutex; ///< Locks any memory reallocations.

#ifdef __PSEUDOCOSTS_
	_RWLOCK *m_rwBrLock; ///< Locks changing arrays `m_dpBrDec` and `m_ipBrNum` used in score branching.
#endif

static int getNumberOfCores(); ///<The function  returns number of physical cores.

	/**
	 * The function starts a new thread.
	 * \param[in] param pointer to `CMIP` objects.
	 * \return always `0`.
	 */
#ifdef _WIN32
	static unsigned int __stdcall startThread(void* param);
#else
	static void* startThread(void* param);
#endif

protected:
	/// The clone constructor.
	/**
	 * The function clones an instance of `CMIP` class.
	 * \param[in] pMip pointer to an instance of CMIP class to be cloned;
	 * \param[in] thread index from {1,2,...,__MAX_THREAD_NUM}.
	*/
	virtual CMIP* clone(const CMIP *pMip, int thread);

	int getThreadNum() const
	{return m_iThreadNum;} ///< \return number of threads currently used.

	/**
	 * In multithreaded application a separate CMIP object is created for each thread.
	 * The threads share common memory, which is allocated and released by the object of thread 0.
	 * \return thread index of calling object.
	 * \see CLP::m_iThread
	 */
	int getThreadIndex() const
	{return m_iThread;}
public:
	/**
	 * If `threadNum > 0`  and `threadNum < __MAX_THREAD_NUM`,
	 * the function sets the number of threads to threadNum;
	 *  otherwise, the function does nothing.
	 *  \param[in] threadNum number of threads.
	 */
	void setThreadNum(int threadNum);
#endif
//////////////////////////////////////////
//     Functions of general use
//////////////////////////////////////////
public:
    /**
     * Call this function to prevent MIPCL from searching for symmetries.
     */
	void skipSymmetrySearch()
	{m_iPreproc&=~0x2;}

    /**
     * The function sets the _separaion rule_ used when not-root node LPs of MIPs.
     *  \param[in] sepRule new separation rule.
     * \sa enSepRule, getSepRule(), CMIP::setMipSepRule().
     */
    void setNodeSepRule(enSepRule sepRule)
    {
    	if (sepRule == SEP_MOST_VIOLATED || sepRule == SEP_STEEPEST_EDGE)
    		m_eSepRule=m_eNodeSepRule=sepRule;
    }

    /**
     * The functions sets the _pricing rule_ to be used for solving not-root node LPs.
     * \param[in] pricingRule new pricing rule.
     * \sa enPricingRule, getPricingRule(), setMipPricingRule().
     */
    void setNodePricingRule(enPricingRule pricingRule)
    {
    	if (pricingRule == PRC_MOST_NEGATIVE || pricingRule == PRC_STEEPEST_EDGE)
    		m_ePricingRule=m_eNodePricingRule=pricingRule;
    }

    /**
     *  /return `true` if the problem being solved is an LP (has no integer variable).
     */
    bool isPureLP() const final;

    /**
     * \param[in] tol new tolerance value for pool cuts;
     *     only if an appropriate cut in the pool is violated by more than `tol`,
     *     then it is extracted to the constraint matrix.
     */
    void setPoolCutTol(double tol)
    	{m_dPoolCtrTol=tol;}

    /**
     * \param[in] tol new tolerance value for cuts currently present in the matrix;
     *     only if a non-basic cut-inequality is violated by more than `tol`,
     *     then the corresponding row can be chosen as the pivot row.
     */
    void setCtrCutTol(double tol)
    	{m_dCtrCutTol=tol;}

private:
    /**
     * Overrides the base class, `CLP`, function.
     * \param[in] i row (constraint) index.
     * \return tolerance value for constraint `i`.
     */
    double getThisCtrTol(int i) const final
    {return (i <= m_iM0 || m_ipRowHd[i] >= 0 || (m_ipCtrType[i] & (CTR_ATTACHED|CTR_STRONG_CUT)))? m_dCtrTol: m_dCtrCutTol;}

public:
    /**
     * The function sets the tolerance value for regular (not cuts) constraints.
     * \param[in] ctrTol maximum violation allowed for all constraints (excluding cuts).
     * \sa CLP::setCtrTol(), CLP::getCtrTol().
     */
    void setCtrTol(double ctrTol) final;

    /**
     *  \param[in] i constraint index.
     *  \return `true` if constraint `i` is global; otherwise, `false`.
     */
    bool isCtrGlobal(int i) const
    {
        return (m_ipCtrType[i] & CTR_LOCAL)? false: true;
    }

    /**
     * \return number of integer variables.
     */
    int getIntegerVarNum() const
        {return m_iIntNum;}

    /**
     * \return number of binary variables.
     */
    int getBinaryVarNum() const
        {return m_iBinNum;}

    /**
     * \return number of real valued variables.
     */
    int getRealVarNum() const
        {return m_iN-m_iIntNum;}

 /**
  * \return maximum number of rounds when probing binary variables.
  * \sa setProbingDepth(), setProbingRoundNum().
  */
    int getProbingDepth() const
        {return m_iProbingDepth;}

/**
 * \param[in] depth new maximum number of rounds when probing binary variables.
 * \sa getProbingDepth(), setProbingRoundNum().
 */
    void setProbingDepth(int depth)
        {m_iProbingDepth=depth;}

/** Integer variables can be assigned priorities, i.e., integer numbers from `VAR_PRI_MIN` to `VAR_PRI_MAX`,
 * which are used when selecting a variable for branching.
 * The higher priority the more chances for a fractional variable to be chosen.
 * \param[in] j index of variable which priority to be changed
 * \param[in] pr new priority
 * \sa incVarPriority(), getVarPriority().
*/
    void setVarPriority(int j, int pr);

 /**
  * The function increases the priority of a given variable by a given value.
  * \param[in] j variable index;
  * \param[in] inc value of increase.
  * \sa setVarPriority(), getVarPriority().
  */
    void incVarPriority(int j, int inc);

 /**
  * \param[in] j variable index.
  *  \return priority of variable `j`.
  *  \sa setVarPriority(), incVarPriority().
  */
    int getVarPriority(int j) const;

    /**
     * The function sets the objective tolerance value:
     * If the difference between the upper and lower bounds on the optimal objective value
     *  is less than the objective tolerance value, the solver returns the best solution found
     *  as a good approximate solution.
     * \param[in] tol tolerance value.
     * \sa getAbsObjTol().
     */
    void setAbsObjTol(double tol)
        {if (tol >= 0.0) m_dAbsObjTol=tol;}

    /**
     * \return objective tolerance value.
     * \sa setAbsObjTol().
     */
    double getAbsObjTol() const
        {return m_dAbsObjTol;}

protected:
    /**
     * This function overload the function of the base class.
     * The functions extends the type of variable `j`
     *  by  bitwise ORing its current type with the flags stored in the parameter `type`.
     * \param[in] j index of variable;
     * \param[in] type bitwise OR of members of `enVarType` and `CMIP::enVarType`.
     * \sa CLP::enVarType, CMIP::enVarType, and `CLP::extendVarType()`.
     */
    void extendVarType(int j, unsigned type) final;

///////////////////////////////////////////
//   Overloaded matrix functions
///////////////////////////////////////////
public:
    /**
     * The function overloads the function of the base class `CLP`.
     * \sa CLP::openMatrix().
     */
	void openMatrix(int m, int n, int nz,
			bool bRowGen=true, bool bColGen=false,
			int mMax=0, int nMax=0, int nzMax=0) final;

    /**
     * The function overloads the function of the base class `CLP`.
     * \sa CLP::closeMatrix().
     */
    void closeMatrix() final;

    /**
     * The function adds a new cut to the matrix.
     * \param[in] hd handle of the constraint,
     *           - if `hd >= 0`, you must overload getRow function;
     *           - if `hd == -1`, the constraint is added to the pool immediately;
     *           - if `hd < -1` (most common option), the constraint is added to the pool if
     *            it is tight (holds as equality) for an optimal LP solution;
     * \param[in] type type of the constraint;
     * \param[in] b1,b2 left hand side (LHS) and right hand side (RHS), respectively;
     * \param[in] sz number of nonzero entries;
     * \param[in] dpVal,ipCol arrays of size `sz`, `dpVal[i]` is coefficient in column `ipCol[i]`;
     * \param[in] bVarScaled is `true` if cut is written in scaled variables;
     * \param[in] factor the constraint has been already scaled (multiplied) by `2^{factor}`,
     *         if `factor==NOT_SCALED`, __MIPCL__ will scale this constraint.
     * \param[in] n if `n > m_iN`, then `n` is the number of variables in the original (non-preprocessed)
     * problem; in this case the constraint must be preprocessed;
     * \return index of the newly created constraint.
     * \throws CMemoryException lack of memory.
     *
     */
    int addCut(tagHANDLE hd, unsigned type, double b1, double b2,
             int sz, double* dpVal, int* ipCol,
             bool bVarScaled=true, int factor=NOT_SCALED, int n=0);

    /**
     * This function is a safer version of `addCut()`.
     * If `dpVal` and  `ipCol` are references to __ MIPCL__ internal arrays such as `m_dpFb` or `m_ipArray`,
     * memory for such arrays may be reallocated during the call to `addCut()`,
     * and, as a consequence, the pointers `dpVal` and `ipCol` becomes not valid.
     * \param[in] hd handle of the constraint,
     *           - if `hd >= 0`, you must overload getRow function;
     *           - if `hd == -1`, the constraint is added to the pool immediately;
     *           - if `hd < -1` (most common option), the constraint is added to the pool if
     *            it is tight (holds as equality) for an optimal LP solution;
     * \param[in] type inequality type (bitwise OR of `enVarType` members);
     * \param[in] b1,b2 left hand side (LHS) and right hand side (RHS), respectively;
     * \param[in] sz number of nonzero entries;
     * \param[in] dpVal,ipCol arrays of size `sz`, `dpVal[i]` is coefficient in column `ipCol[i]`;
     * \param[in] n if `n > m_iN`, then `n` is the number of variables in the original (non-preprocessed)
     * problem; in this case the constraint must be preprocessed;
     * \param[in] bVarScaled is `true` if cut is written in scaled variables;
     * \param[in] factor the constraint has been already scaled (multiplied) by `2^{factor}`,
     *         if `factor==NOT_SCALED`, __MIPCL__ will scale this constraint.
     * \return index of the newly created constraint.
     * \throws CMemoryException lack of memory.
     * \sa addCut().
    */
    int safeAddCut(tagHANDLE hd, unsigned type, double b1, double b2,
             int sz, double* &dpVal, int* &ipCol,
             bool bVarScaled=true, int factor=NOT_SCALED, int n=0);

private:
    /**
     * The function compute a hash value `(hash1,hash2)` for a given constraint
     *     `b1 <= sum(i in 0..sz-1) dpVal[i]*x(ipCol[i]) <= b2`.
     * This function is used to prevent adding identical cuts.
     *
     * \param[in] sz number of entries;
     * \param[in] ipCol column indices, list of size 'sz';
     * \param[in] dpVal coefficients, list of size 'sz';
     * \param[in] b1,b2 left and right hand sides;
     * \param[in] factor input inequality is to be multiplied by `2^{factor}`.
     * \return 'true' if cut is new; otherwise, 'false'.
     * \remark 'm_dpW' is used to store hash values.
     */
    bool isNewCut(int sz, const int* ipCol, const double* dpVal, double b1, double b2, int factor=0);

//////////////////////////////////////////////////////////////////////////
//     C O N S T R U C T O R S / D E S T R U C T O R                    //
//////////////////////////////////////////////////////////////////////////
    void resetCutProperties(); ///< This procedure resets to default values the parameters that affect cut generation.
public:
    /**
     *  The constructor initializes an "empty" `CMIP` object and sets the problem name to `name`.
     *  \param[in] name problem name of no more than '30' characters; if `name==0`, problem name will be a featureless "LP".
     */
    CMIP(const char* name);

#ifndef __ONE_THREAD_
    /// Clone constructor.
	/**
	 * It is  used in multithreaded applications.
	 * The cloned object shares many members with the object being cloned.
	 * \param[in] other reference to a `CMIP` object being cloned;
	 * \param[in]  thread thread-index.
     * \throws `CMemoryException` lack of memory.
     */
    CMIP(const CMIP &other, int thread);
#endif


    virtual ~CMIP(); ///< The destructor.

// ////////////////////////////////////////////////////////////////////////
//                   Initialization/memory allocation                   //
// ////////////////////////////////////////////////////////////////////////
    /**
     * Call this function only in case when cuts of all types are not supposed to be generated.
     * In such a case you can still overload separate(), genCut1(), and genCut2() procedures
     * to generate problem specific constraints. In such a case, you have to overload
     * getRow() function from the base class LP.
     */
    void doNotUsePool()
        {m_bPool=false;}

private:
    void allocMemForBC(); ///< allocates memory for the search tree, pool, and for some working arrays.
    /**
     * The function increases the maximum number of columns in the matrix,
     * and then reallocates memory for all structures which size depend on the number of columns.
     * \throws CMemoryException lack of memory.
     * \sa CLP::incMaxColumnNumber().
     */
    void incMaxColumnNumber() final;

    /**
     * The function increases the maximum number of rows in the matrix,
     * and then reallocates memory for all structures which size depend on the number of rows.
     * \throws CMemoryException lack of memory.
     * \sa CLP::incMaxColumnNumber().
     */
    void incMaxRowNumber() final;

   	void reallocMemForEntries(int nz) final; ///< thread safe version of the base-class-function.

protected:
    /**
   	 * \param[in] j variable index.
   	 * \return `true` if variable `j` is binary.
     */
    bool isVarBinary(int j) const final
    {
        return (m_ipVarType[j] & VAR_BIN)? true: false;
    }

   	/**
   	 * Using different preprocessing techniques one can deduce that some variables initially not declared as integer-valued
   	 * may take only integer values. The solver declares such variables as integral
   	 * (this may help the solver in generating stronger cuts) but assigns them the lowest priority to avoid branching by such variables.
   	 * It worth noting, that, for technical reasons, __MIPCL__ assigns to variables of the lowest priority the maximum possible priority value.
   	 * \param[in] j variable index.
   	 * \return `true` if variable `j` is integer-valued.
   	 * \sa isVarUsedForBranching().
   	 */
    bool isVarStrongIntegral(int j) const final
    {
        return ((m_ipVarType[j] & VAR_INT) && (getVarPriority(j) > VAR_PRI_MIN))? true: false;
    }

    /**
     * We say that a variable is _used for branching_ if its lower and upper bounds
     *  may be changed during the branching process.
     * Let us also note that real-valued variables (say, those in SOS2-constraints) can be fixed during the branching process.
     * So, they are used for branching as well as strong integral variables.
   	 * \param[in] j variable index.
   	 * \return `true` if variable `j` can be chosen for branching.
   	 * \sa isVarUsedForBranching().
     * \sa isVarStrongIntegral(), isVarSOS(), startBranching().
     */
    bool isVarUsedForBranching(int j) const final;

   /**
    * \param[in] j variable index.
   	 * \return `true` if variable `j` appears in a SOS1 or SOS2-constraint.
   	 * \sa isVarUsedForBranching().
    */
    bool isVarSOS(int j) const
    {
        return (m_ipVarType[j] & VAR_SOS)? true: false;
    }

    /**
     * Strong integral and SOS variables are _not scalable_.
    * \param[in] j variable index.
   	 * \return `true` if variable `j` is _scalable_.
   	 * \sa isVarStrongIntegral(), isVarSOS().
    */
    bool isVarScalable(int j) const final
    {
        return (m_ipVarType[j] & (VAR_NOT_MOD | VAR_INT | VAR_SOS))? false: true;
    }

    /**
     *  Constraint is declared as being _integral_ if all its variables are integer-valued as well as all its coefficients are integers.
     *  The slack of an inequality (defined as the difference between its right and left hand sides) takes only integral value;
     *  this property may be useful in cut generating procedures.
     *  \param[in] i constraint index.
     *  \return `true` if constraint `i` is integral; otherwise, `false`.
     */
    bool isCtrIntegral(int i) const
    {
        return ((m_ipCtrType[i] & CTR_INT_VARS) && (m_ipCtrType[i] & CTR_INT_COEFF))? true: false;
    }

    /**
     * The objective is _integral_ if its optimal value is an integer.
     * \return `true` if the objective is integral; otherwise, `false`.
     */
    bool isObjIntegral() const
    {return m_bObjInt;}

    /**
     * For a not-root node, a lower variable bound for a variable is _local_ if its value differs from that at the root node.
     * \param[in] j variable index.
     * \return `true` if the lower bound of variable `j` is local; otherwise, `false`.
     * \sa isVarLoBoundLocal().
     */
    bool isVarUpBoundLocal(int j) const;

    /**
     * For a not-root node, an upper  bound for a variable is _local_ if its value differs from that at the root node.
     * \param[in] j variable index.
     * \return `true` if the lower bound of variable `j` is local; otherwise, `false`.
     * \sa isVarUpBoundLocal().
     */
    bool isVarLoBoundLocal(int j) const;

	/**
	 * A variable is _monotone up_ if, given any solution of the relaxation LP,
	 *   we can increase the value of the considered variable without violating any constraint other than
	 *   the upper bound of this variable.
	 * \param[in] j variable index.
	 * \return `true` if variable `j` is monotone up; otherwise, `false`.
	 * \sa isVarMonotoneDown(), isVarMonotone().
	 */
    bool isVarMonotoneUp(int j) const
        {return (m_ipVarType[j] & VAR_MON_UP)? true: false;}

	/**
	 * A variable is _monotone down_ if, given any solution of the relaxation LP,
	 *   we can decrease the value of the considered variable without violating any constraint other than
	 *   the lower bound of this variable.
	 * \param[in] j variable index.
	 * \return `true` if variable `j` is monotone down; otherwise, `false`.
	 * \sa isVarMonotoneUp(), isVarMonotone().
	 */
    bool isVarMonotoneDown(int j) const
        {return (m_ipVarType[j] & VAR_MON_DOWN)? true: false;}

	/**
	 * A variable is _monotone_ if it is either monotone up or monotone down.
	 * Usually, a monotone variable is among the last candidates for choosing it as the branching variable
	 * since it is always possible to round off such a variable.
	 * \param[in] j variable index.
	 * \return `true` if variable `j` is monotone; otherwise, `false`.
	 * \sa isVarMonotoneUp(), isVarMonotone().
	 */
    bool isVarMonotone(int j) const
        {return (m_ipVarType[j] & (VAR_MON_UP | VAR_MON_DOWN))? true: false;}

private:
    /**
     * \param[in] j index of variable.
     * \return `true` if variable `x[j]` is in a variable bound constraint.
     */
    bool isVarInVarBound(int j) const
    {return (!(m_ipVarType[j] & VAR_BIN) && (m_ipVarType[j] & VAR_IN_VAR_CTR))? true: false;}

    /**
     * \param[in] j index of variable.
     * \return `true` if variable `x[j]` is in a variable upper bound constraint.
     */
    bool isVarInVarUpperBound(int j) const
    {return (!(m_ipVarType[j] & VAR_BIN) && (m_ipVarType[j] & VAR_IN_VAR_UB))? true: false;}

    /**
     * \param[in] j index of variable.
     * \return `true` if variable `x[j]` is in a variable lower bound constraint.
     */
    bool isVarInVarLowerBound(int j) const
    {return (!(m_ipVarType[j] & VAR_BIN) && (m_ipVarType[j] & VAR_IN_VAR_LB))? true: false;}

    /**
     * \param[in] i index of constraint.
     * \return `true` if constraint indexed by `i`  is a variable upper bound.
     */
	bool isCtrVarUb(int i) const
	{return  ((m_ipCtrType[i] & CTR_VAR) && (m_ipCtrType[i] & CTR_RIGHT))? true: false;}

    /**
     * \param[in] i index of constraint.
     * \return `true` if constraint indexed by `i`  is a variable lower bound.
     */
	bool isCtrVarLb(int i) const
	{return  ((m_ipCtrType[i] & CTR_VAR) && (m_ipCtrType[i] & CTR_LEFT))? true: false;}

	/**
	 *  All but SOS1 and SOS2-variables (i.e., those included in SOS1 SOS2-constraints) can be deleted from the matrix by the solver.
	 *  \param[in] j index of variable.
     * \return `true` if variable `j` can be deleted.
	 */
    bool isVarDeletable(int j) const final
    {
        return (m_ipVarType[j] & (VAR_NOT_MOD | VAR_SOS))? false: true;
    }

    /**
     *  All but SOS1 and SOS2-constraints can be modified by the solver.
     * \param[in] i constraint index.
     * \return `true` if constraint `i` can be modified.
     */
    bool isCtrModifyable(int i) const final
    {
    	return (m_ipCtrType[i] & (CTR_SOS1 | CTR_SOS2))? false: true;
    }

    /**
     * \param[in] i constraint index.
     * \return `true` if constraint `i` is a GUB-constraint.
     */
    bool isCtrGUB(int i) const
    {
    	return ((m_ipCtrType[i] & CTR_PACKING) && (m_ipCtrType[i] & CTR_CARDINALITY))? true: false;
    }

///////////////////////////////////////////////////////
//  Branch and Cut attributes
//////////////////////////////////////////////////////
public:
    /**
     * Any MIP is transformed into the form when the objective is maximized.
     * The function returns the transformed objective value of the currently best solution.
     * It is usually called when the problem is being solved.
     * \return objective value of the best MIP solution found.
     * \sa getObjVal(), safeGetObjLowerBound().
     */
    double getObjLowerBound() const;

    /**
     * A thread safe version of getObjLowerBound().
     * \return `m_dAbsObjTol` plus the objective value of the best MIP solution found.
     * \sa getObjVal(), safeGetObjLowerBound().
     * \remark this function locks for reading `m_pRecord`.
     */
    double safeGetObjLowerBound();
private:
    /**
     * This function is called from `CLP::dualSimplex()`, to update the new value of `m_dLoBound`.
     * If some thread finds a new best solution and the current objective value exceeds this new lower bound,
     * `CLP::dualSimplex()` stops running.
     */
    double getLastLowerBound() final;

public:
    /**
     * This function is usually called when the problem has been solved already.
     * \return original objective value of the best MIP solution found.
     * \sa getObjLowerBound(), safeGetObjLowerBound().
     */
    double getObjVal() const;

    /**
     * \return number of MIP solutions found so far.
     * \sa safeGetSolNum().
     */
    int getSolNum() const;

    /**
     * Thread safe version getSolNum().
     * \return number of MIP solutions found so far.
     * \sa  getSolNum().
     */
    int safeGetSolNum();

    /**
     * \return index of the currently processed node.
     * \sa getCurrentNodeHeight().
     */
    int getCurrentNode() const
        {return m_iNode;}

    /**
     * \return height of the currently processed node.
     * \sa getCurrentNode().
     */
    int getCurrentNodeHeight() const;

    /**
     * An integer variable is assumed to be integral if the absolute difference
     * between its value and the nearest integer to this value is less than the integrality tolerance.
     * \param[in] intTol new integrality tolerance.
     * \sa getIntTol().
     */
    void setIntTol(double intTol)
        {m_dIntTol=intTol;}

    /**
     * \return integrality tolerance.
     * \sa setIntTol().
     */
    double getIntTol() const
        {return m_dIntTol;}

   /**
    * This function sets a new search rule.
    * \param[in] rule search rule to for low-height nodes;
    * \sa enBranchRule.
    */
    void setBranchingRule(enBranchRule rule)
    {
    	m_eBranchRule=rule;
    }

    /**
     * The solver automatically generates cuts for first `m_iAutoCutNodes` processed nodes,
     * and also for nodes of height less or equal than `m_iAutoCutNodeHeight`.
     * \param[in] nodeNum number of nodes for which cuts are generated;
     * \param[in] height cuts are generated for nodes of height `<= height`.
     * \sa setCutTypePattern(), setAutoCutRounds().
     */
    void setAutoCutPattern(int nodeNum, int height)
	{
		m_iAutoCutNodes=nodeNum;
		m_iAutoCutNodeHeight=height;
	}

    /**
     * Cuts of size greater than some predefined limit are rejected.
     * This procedure sets new limits.
     *
     * \param[in] sizeAtRoot maximum cut size at root nodes;
     *    if `sizeAtRoot` is negative, then the existing limit remains unchanged;
     * \param[in] sizeAtNodes maximum cut size at non-root nodes;
     *    if `sizeAtNodes` is negative, then the existing limit remains unchanged.
     */
    void setMaxCutSize(int sizeAtRoot, int sizeAtNodes)
    {
    	if (sizeAtRoot >= 0)
    		m_iMaxCutSizeAtRoot=sizeAtRoot;
    	if (sizeAtNodes >= 0)
    		m_iMaxCutSizeAtNodes=sizeAtNodes;
    	if (m_iMaxCutSizeAtRoot < m_iMaxCutSizeAtNodes)
    		m_iMaxCutSizeAtRoot=m_iMaxCutSizeAtNodes;
	}

    /**
     * The solver automatically generates cuts of type `cutType` only for first `m_ipCutNodes[cutType]` processed nodes.
     * \param[in] cutType cut type;
     * \param[in] nodeNum number of nodes for which cuts of type `cutType` are generated;
     * \param[in] height cuts of type `i` are generated for nodes of height `<= height`.
     * \sa enCutType, setAutoCutPattern();
    */
    void setCutTypePattern(enCutType cutType, int nodeNum, int height)
    {
    	m_iUserCutPattern|= (1 << cutType);
    	m_pCutProps[cutType].cutNodes=nodeNum;
    	m_pCutProps[cutType].cutNodeHeight=height;
    }

    /**
     * At any node the solver calls cut generating procedures a restricted number of times.
     * \param[in] atRoot at root node, solver repeatedly calls its cut generating procedure at most `atRoot` times.
     * \param[in] atNodes at not-root nodes, solver repeatedly calls its cut generating procedure at most `atNodes` times.
     * \sa setAutoCutPattern().
     */
    void setAutoCutRounds(int atRoot, int atNodes)
        {m_iMaxAutoCutRoundsAtRoot=atRoot; m_iMaxAutoCutRoundsAtNodes=atNodes;}

    /**
     * The solver adds to the matrix as well as to the pool cuts of restricted size.
     * \param[in] size new value of maximum size.
     */
    void setMaxCutSize(int size)
    {m_iMaxCutSize=size;}

    /**
     * The rules that are used by the solver
     *  to choose a right moment to stop generating cuts are complicated.
     *  The progress in decreasing the upper bound
     *  on the optimal objective value is the key factor affecting the whole cut generation process.
     *  setObjDecFracPerCutRoundAtRoot() sets a new minimum relative decrease of the upper bound (objective value of the relaxation LP)
     *  per cut round for the root node.
     * \param[in] prc minimum relative decrease of the upper bound per cut round for the root node.
     * \sa getRelObjDecPerCutRoundAtRoot(), setRelObjDecPerCutRoundAtNodes().
     */
	void  setRelObjDecPerCutRoundAtRoot(double prc)
	{
		if (prc > 1.0e-6) m_dCutRoundPrcAtRoot=prc;
	}

	/**
     * \return prc minimum relative decrease of the upper bound per cut round at the root node.
     * \sa setRelObjDecPerCutRoundAtRoot().
	 */
    double  getRelObjDecPerCutRoundAtRoot() const
	{return m_dCutRoundPrcAtRoot;}

    /**
     *  This function sets a new minimum relative decrease of the upper bound (objective value of the relaxation LP)
     *  per cut round for not-root nodes.
     * \param[in] prc minimum relative decrease of the upper bound per cut round for not-root nodes.
     * \sa getRelObjDecPerCutRoundAtNodes(), setRelObjDecPerCutRoundAtRoot().
     */
	void  setRelObjDecPerCutRoundAtNodes(double prc)
	{
		if (prc > 1.0e-6)
			m_dCutRoundPrc=prc;
	}

	/**
     * \return prc minimum relative decrease of the upper bound per cut round for not-root nodes.
     * \sa setRelObjDecPerCutRoundAtNodes().
	 */
	double getRelObjDecPerCutRoundAtNodes() const
	{return m_dCutRoundPrc;}

	/**
	 * The procedure sets minimum numbers of rounds to generate cuts at the root and not-root nodes.
	 * \param[in] atRoot minimum number of rounds to generate cuts at root node;
	 * \param[in] atNodes minimum number of rounds to generate cuts at not-root nodes.
	 */
	void setMinCutRounds(int atRoot, int atNodes)
		{m_iMinAutoCutRoundsAtRoot=atRoot; m_iMinAutoCutRoundsAtNodes=atNodes;}

	/**
	 * The solver calls the cut generating procedure for a particular cut type a limited number of times.
	 * This function reestablishes that limit.
     * \param[in] cutType cut type;
     * \param[in] roundNum number of rounds.
     * \sa enCutType, setMaxCutsPerRound().
	 */
    void setMaxCutRoundNum(enCutType cutType, int roundNum)
        {m_pCutProps[cutType].maxCutRounds=roundNum;}

	/**
	 * At each cut generation round, the cut generating procedure for a particular cut type, `cutType`,
	 * generates only those cuts which size is not greater than the value of `m_ipMaxCutSize[cutType]`.
	 * This function sets a new limit.
     * \param[in] cutType cut type;
     * \param[in] size maximum cut size.
     * \sa enCutType.
	 */
    void setMaxCutSize(enCutType cutType, int size)
    	{m_pCutProps[cutType].maxCutSize=size;}

	/**
	 * At each cut generation round, the cut generating procedure for a particular cut type
	 * generates no more than a predefined number of cuts.
	 * This function reestablishes that limit.
     * \param[in] cutType cut type;
     * \param[in] cutNum number of cuts.
     * \sa enCutType, setMaxCutRoundNum().
	 */
    void setMaxCutsPerRound(enCutType cutType, int cutNum)
    	{m_pCutProps[cutType].maxCutsPerRound=cutNum;}

    /**
     * \param[in] cutType cut type;
     * \param[in] tol cut inequality is assumed as being violated, and then added to the matrix,
     *  if the inequality is violated by at least `tol`.
     * \sa enCutType, getToleranseForCut().
     */
    void setToleranceForCut(enCutType cutType, double tol)
        {m_pCutProps[cutType].cutTol=tol;}

    /**
     * \param[in] cutType cut type.
     * \return tolerance value for the cuts of type `cutType`.
     * \sa enCutType, setToleranseForCut().
     */
    double getToleranceForCut(enCutType cutType) const
        {return m_pCutProps[cutType].cutTol;}

    /**
     * Usually, cuts are sent to the pool
     * if they are in the matrix when the processed subproblem is added to the search tree as a new leaf.
     * Call this function to immediately send to the pool all newly generated cuts of type `cutType`.
     * \param[in] cutType cut type.
     * \sa enCutType, areCutsSentToPool.
     */
    void sendCutsToPool(enCutType cutType)
    	{m_uiCutsToPool|=(1 << cutType);}

    /**
     * \param[in] cutType cut type.
     * \return `true` if cuts of type `cutType` are immediately sent to the pool.
     * \sa enCutType, sendCutsToPool.
     */
    bool areCutsSentToPool(enCutType cutType) const
    	{return (m_uiCutsToPool & (1 << cutType))? true: false;}

    /**
     * Mod2-cuts (or (0,1/2)-Chvatal-Gomory cuts) are considered as sparse
     *  if the number of variables in them is less than `sparsePrc*m_iN/100`
     * (let us remember that `m_iN` is the number of a` variables),
     * and as dense if the number of variables in them is less than `sdensePrc*m_iN`.
     * Any cut of size greater than `sdensePrc*m_iN/100` is considered as super-dense and is rejected (not generated).
     * \param[in] sparsePrc percentage (relative to the total number of variables in the problem) of variable in sparse cuts;
     * \param[in] densePrc percentage of variable in dense cuts.
     */
    void setMod2CutDensity(int sparsePrc, int densePrc)
        {m_pCutProps[_SPARSE_MOD2].maxCutPrc=sparsePrc; m_pCutProps[_DENSE_MOD2].maxCutPrc=densePrc;}

    /**
     * Gomory cuts are considered as sparse if the number of variables in them is less than `sparsePrc*m_iN/100`
     * (let us remember that `m_iN` is the number of all variables),
     * and as dense if the number of variables in them is less than `sdensePrc*m_iN`.
     * Any cut of size greater than `sdensePrc*m_iN/100` is considered as super-dense and is rejected (not generated).
     * \param[in] sparsePrc percentage (relative to the total number of variables in the problem) of variable in sparse cuts;
     * \param[in] densePrc percentage of variable in dense cuts.
     */
    void setGomoryCutDensity(int sparsePrc, int densePrc)
        {m_pCutProps[_SPARSE_GOMORY].maxCutPrc=sparsePrc; m_pCutProps[_DENSE_GOMORY].maxCutPrc=densePrc;}
private:
  	/// \return number of generated cuts of type `type`.
  	int getCutGenerated(int type) {return m_ipTotalCutNum[type<<1];}

  	/// \return number of cuts of type `type` used.
  	int getCutUsed(int type) {return m_ipTotalCutNum[(type<<1)+1];}
///////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
//  Classification Reformulation, and preprocessing   //
////////////////////////////////////////////////////////
    /**
     *  The function marks all variables present in SOS1 and SOS2-constraints as SOS-variables.
     *  \see `CTR_SOS1`, `CTR_SOS2` and `VAR_SOS`.
     */
    void selectSOSvars();

    /**
     *  The function labels integer-valued variables
     *   with lower and upper bounds, respectively, `0` and `1` as _binary_.
     */
    void selectBinaries();

    /**
     *  The procedure seeks for flow structures to classify flow variables as integer ones.
     */
	void seekFlowStructures();

    /**
     *  The procedure looks for flow structures to classify flow variables as integer ones.
     */
	void checkTrFlowStructure();

	/**
	 *  The procedure looks for identical columns.
	 */
	int seekIdenticalCols();

	/**
	 * The procedure, for each invariant knapsack equation, seeks for
	 * the constraints that contain all variables from this equation,
	 * and then the invariant knapsack equation multiplied by an appropriate number
	 * is added to each found constraint.
	 * The goal is to simplify the constraint matrix by reducing the number of non-zero entries.
	 *
	 */
	int seekIncludedInvKnapsacks();

// Constraint classification
	/**
	 * The procedure detects SOS2-constraints.
	 *
	 * \return number of SOS2-constraits detected.
	 * \sa isSOS2().
	 */
    int seekSOS2();

    /**
     * The functions verifies whether a given row (together with some other rows) represents a __SOS2__ constraint.
     * \param[in] r row index.
     * \return `true` if row `r` represents a __SOS2__ constraint.
     * \remark the function has not the `const` attribute since it may alter internal `MIPCL` arrays.
     */
    bool isSOS2(int r);

    /**
     * \param[in] r row index.
     * \return `true` if row `r` represents a _parity_ constraint.
     */
    bool isParityCtr(int r) const;

    /**
     * The procedure tries to multiply a given constraint, in turn, by `min{a(r,j): j in supp(A_r)}` and `600`,
     * to get all coefficients integral. If such a multiplier exists, let us denote it by `a`,
     * and the constraint contains at most one non-integral variable,
     * then this constraint is multiplied by `a`.
     * \param[in] row constraint index.
     */
    void makeCtrCoeffIntegral(int row);

    /**
     * The procedure does for cuts the same as `makeCtrCoeffIntegral()` does for matrix constraints.
     * \param[in] sz size of arrays `dpVal` and `ipCol`;
     * \param[in,out] dpVal,ipCol `sum(i=0,..,sz-1) dpVal[i]*x[ipCol[i]] <= b` is input inequality;
     * \param[in,out] b right hand side of input inequality.
     * \sa makeCtrCoeffIntegral().
     */
    void makeCutCoeffIntegral(int sz, double *dpVal, int *ipCol, double &b);

    /**
     * The procedure establishes additional properties of constraints that are not fully classified.
     * \param[in] ind constraint index.
     * \return bitwise OR of `enVarType` members.
     */
    unsigned classifyCtrPlus(int ind);

    /**
     * The function classifies constraint `row`.
     * \param[in] row index of constraint.
     */
    void classifyCtr(int row);

    /**
     * The procedure classifies constraints and variables.
     * \param[in] bPlus when applied within preprocessing step with `bPlus` set to `false`,
     *  some properties of constraints, e.g. those that are used only when cuts are generated,
     *  are skipped to be detected.
     */
    void classify(bool bPlus=true);

    /**
     * The procedure builds linked lists of constraint indices
     *  that are written in the array `m_ipCtrLink`.
     * For example, the linked list of all knapsack constraints
     *  starts with the constraint indexed by `i(0)=m_iLastKnapsack`,
     *  the next knapsack constraints are `i(1)=m_ipCtrLink[i(0)],...,i(k)=m_ipCtrLink[i(k-1)]`
     *  with `m_ipCtrLink[i(k)=-1` indicating the last constraint in the list.
     *  If `m_iLastKnapsack=-1`, then there is no knapsack constraint.
     *  \param[in] bUpdate if set to `true`, currently existing list is updated;
     *  otherwise, the list is build from scratch.
     */
    void buildCtrLists(bool bUpdate=false);

///// PREPROCESSING //////////////////////////////////////

    /**
     * Substitutes variable `x[col]` with the expression `b + sum(j in 0..sz-1) dpVal[j]*x[ipCol[j]]`.
     * \param[in] j variable to be deleted;
     * \param[in] hd handle of deleted variable, or `hd in {-1,-2}`;
     *   if `hd = -1`, then `sz = 1` and variable `x[ipCol[0]]` is substituted for `x[ipCol[0]]*dpVal[0]+b`;
     *   if `hd = 2`, then `sz = 1`, and, for `c=reinterpret_cast<int*>(dpVal)[0]`,
     *   if `x[ipCol[0]] >= 0.0`, then `x[c]=0.0`; otherwise, `x[c]=-x[ipCol[0]]` and `x[ipCol[0]]=0.0`;
     * \param[in] b free term in the substitution
     * \param[sz] sz number of variables in the substitution;
     * \param[in] ipCol list of `sz` variables;
     * \param[in] dpVal list of `sz` coefficients.
     */
    void deleteVariable(int j, int hd, double b, int sz, int *ipCol, double *dpVal) final;

    /**
     * The function calls `reduceCoefficients()` to tighten the constraint given by its arguments,
     * and then classifies the resulting inequality.
     * \param[in,out] sz number of variables in constraint to be classified; on output, `sz` may be smaller because some coefficients may be zeroes;
     * \param[in,out] dpVal,ipCol arrays of size `sz`, `dpVal[i]` is coefficient for variable `ipCol[i]`;
     * \param[in,out] b left or right hand side of constraint;
     * \param[in] side if `true` (resp., `false`), `b` is right (resp., left) hand side;
     * \param[in] makeIntCoeff if `true`, `makeCutCoeffIntegral()` is called.
     * \return constraint type.
     */
	int classifyAndReduceCoeff(int &sz, double* dpVal, int* ipCol, double&  b,  bool side, bool makeIntCoeff=false);

	/**
	 * The function tries to make integral the coefficients of the cut given by first four arguments:
	 *     `sum(i in 0..sz) dpVal[i]*x[ipCol[i]] <= (>=) b`.
	 *
     * \param[in] size number of variables in constraint to be classified;
     * \param[in,out] dpVal,ipCol arrays of size `sz`, `dpVal[i]` is coefficient for variable `ipCol[i]`;
     * \param[in,out] b right hand side of constraint;
     * \param[in] side if `true`, the sign in of the input inequality is `<=`; otherwise, the sign is `>=`.
     * \return constraint type (bitwise OR of `enVarType` members).
	 */
	unsigned simplifyCut(int &size, double *dpVal, int *ipCol, double &b, bool side=true);

//////////////// P R O B I N G
	/**
	 * To speed up probing of binary variables, not all constraint
	 * are processed at probing rounds. Usually, only constraints of small size
	 * and some structured constraints - say, packing constraints - are useful
	 * when probing binary variables.
	 * \param[in] r row index;
	 * \param[in] dpBd for row `r`, `dpBd[r<<1]` and `dpBd[(r<<1)+1]` are, respectively,
	 *            lower and upper bounds on left-hand-side value of constraint `r`.
	 * \return `true` is the constraint indexed by `r` is to be probed.
	 */
		bool isCtrUsedInProbing(int r, double * dpBd);

	/**
	 * When a binary variable is fixed, this procedure is called in probing subroutines
	 * to derive implications from the implication matrix.
	 *
	 * \param[in] row row-index in implication matrix;
	 * \param[in,out] last last constraint in list (queue) of active and non-processed constraints;
	 * \param[in,out] ipCtrLst list of active and non-processed constraints;
	 * \param[in,out] cpActiveCtr  for row `r`, `cpActiveCtr[r]` is `true` if `r` is in `ipCtrLst`;
	 * \param[in,out] dpBd for row `r`, `dpBd[r<<1]` and `dpBd[(r<<1)+1]` are, respectively,
	 *            lower and upper bounds on left-hand-side value of constraint `r`;
	 * \param[in,out] dpD for column `j`, `dpD[j<<1]` and `dpD[(j<<1)+1]` are, respectively,
	 *            lower and upper bounds for variable `j`;
	 * \param[in] round probing round;
	 * \param[in,out] cpRound for column `j`, `cpRound[j<<1]` and `cpRounf[(j<<1)+1]` are
	 *       last probing rounds when, respectively, lower and upper bounds for variable `j` changed.
	 * \return `false` if propagation results in any inconsistency; otherwise, `true`.
	 */
	bool propagateImplications(int row, int& last, int *ipCtrLst, char *cpActiveCtr,
			double* dpBd, double* dpD, int round, unsigned char *cpRound);

	/**
	 * The procedure is called only within `PROBE_preprocessCtr()` when the lower or upper bound on variable `col` was changed.
	 * \param[in] m number of rows in preprocessed submatrix;
	 * \param[in] col index of variables whose bound was changed;
	 * \param[in] side1 if `true`, lower bound of variable `col` was changed;
	 * \param[in] side2 if `true`, upper bound of variable `col` was changed;
	 * \param[in] delta1 value by which lower bound of `col` was changed;
	 * \param[in] delta2 value by which upper bound of `col` was changed;
	 * \param[in] round iterate (round) of `probe()`;
	 * \param[in,out] last index of last constraint in list of active (to be preprocessed) constraints;
	 * \param[in,out] ipCtrLst list (of maximum size `m`) of active constraints;
	 * \param[in,out] cpActiveCtr constraint `r` is active if `cpActiveCtr[r]==1`;
	 * \param[in,out] dpBd array of size `m`; `dpBd[r<<1]` and `dpBd[(r<<1)+1]` are respectively lower and upper values of constraint `r`.
	 * \return `false` if any inconsistency has been detected; otherwise, `true`.
	 */
	bool correctBounds(int m, int col, bool side1, bool side2, double delta1, double delta2,
			 int round, int& last, int *ipCtrLst, char *cpActiveCtr, double* dpBd);

	/**
	 * The procedure is called from `probe()` to preprocess a given constraint, i.e.,
	 * to improve bounds on variables and the constraint being processed.
	 * \param[in] ind index of constraint;
	 * \param[in,out] last constraint in queue (given by array `ipCtrLst`) of  constraints to be preprocessed;
	 * \param[in,out] ipCtrLst queue of  constraints to be preprocessed;
	 * \param[in,out] cpActiveCtr if `cpActiveCtr[r]`, then inequality `r` belongs to
	 *  queue of  constraints to be preprocessed;
	 * \param[in,out] dpBd `dpBd[r<<1]` and `dpBd[(r<<1)+1]` are, respectively, lower and upper bound on value of constraint `r`;
	 * \param[in,out] dpD `dpD[j<<1]` and `dpD[(j<<1)+1]` are, respectively, lower and upper bound on value of variable `j`;
	 * \param[in] round preprocessing round;
	 * \param[in,out] cpRound value of `cpRound[j<<1]` (resp., `cpRound[(j<<1)+1]`) is last round at which
	 * lower (resp., upper) bound on variable `j` has been changed.
	 * \return `false` if constraint `ind` is infeasible subject to the bounds on variables written in `dpD`;
	 * otherwise, the procedure returns `true`.
	 */
	bool PROBE_preprocessCtr(int ind, int& last, int *ipCtrLst, char *cpActiveCtr,
        double* dpBd, double* dpD, int round, unsigned char *cpRound);

	/**
	 * The procedure at each preprocessing _round_ repeatedly calls `PROBE_preprocessCtr()`
	 *  to preprocess constraints from the queue given by its three first parameters.
	 *  `PROBE_preprocessCtr()` adds to the queue new constraints to be preprocessed in the next processing round.
	 * \param[in,out] last constraint in queue (given by array `ipCtrLst`) of  constraints to be preprocessed;
	 * \param[in,out] ipCtrLst queue of  constraints to be preprocessed;
	 * \param[in,out] cpActiveCtr if `cpActiveCtr[r]`, then inequality `r` belongs to
	 *  queue of  constraints to be preprocessed;
	 * \param[in,out] dpBd `dpBd[r<<1]` and `dpBd[(r<<1)+1]` are, respectively, lower and upper bound on value of constraint `r`;
	 * \param[in,out] dpD `dpD[j<<1]` and `dpD[(j<<1)+1]` are, respectively, lower and upper bound on value of variable `j`;
	 * \param[in,out] cpRound value of `cpRound[j<<1]` (resp., `cpRound[(j<<1)+1]`) is last round at which
	 * lower (resp., upper) bound on variable `j` has been changed;
	 * \param[in] implRow row in implication matrix.
	 * \return `false` if constraint `ind` is infeasible subject to the bounds on variables written in `dpD`;
	 * otherwise, the procedure returns `true`.
	 */
	bool probe(int& last, int *ipCtrLst, char *cpActiveCtr,
		double *dpD, double *dpBd, unsigned char *cpRound, int implRow=-1);

	/**
	 * The procedure prepares a `CImpl` object for storing implications between binary variables,
	 * as well as variable bounds.
	 * \param[in] impl pointer to `CImpl` object.
	 */
	void probeInit(CImpl *impl);

	/**
	 * The procedure updates the lower and upper values of all constraints
	 * when a bound (lower or upper) of a given variable were changed.
	 *
	 * \param[in] c column index;
	 * \param[in] side,delta if `side==true`, upper bound of variable `x[c]` is decreased by `delta`;
	 * otherwise, lower bound of `x[c]` is increased by `delta`;
	 * \param[in,out] dpBd array of size `2*m_iM`, `dpBd[i<<1]` and `dpBd[(i<<1)+1]`
	 * are respectively lower and upper values of constraint `i`.
	 */
	void updateRowBds(int c, bool side, double delta, double *dpBd);

	/**
	 * This procedure is called after probing a binary variable to tighten constraints.
	 *
	 * \param[in] c column index of binary variable having been probed;
	 * \param[in] side if `true`, variable `c` was fixed to `1`; otherwise, to `0`.
	 * \return number of tightened constraints
	 */
	int tightenCtrs(int c, bool side);

	/**
	 * After probing a binary variable `c1`, the procedure adds a new lower variable bound, `x[c2]-a*x[c1] >= b`, on variable `c2`.
	 *
	 * \param[in] c1 variable having been probed;
	 * \param[in] c2 dependent variable whose lower bound has been increased when probing variable `c1`;
	 * \param[in] side `true` if `x[c1]=1` when probing; otherwise, value of `side` is `false`.
	 */
	void addNewLoVarBound(int c1, int c2, bool side);

	/**
	 * After probing a binary variable `c1`, the procedure adds a new upper variable bound, `x[c2]-a*x[c1] <= b`, on variable `c2`.
	 *
	 * \param[in] c1 variable having been probed;
	 * \param[in] c2 dependent variable whose upper bound has been increased when probing variable `c1`;
	 * \param[in] side `true` if `x[c1]=1` when probing; otherwise, value of `side` is `false`.
	 */
	void addNewUpVarBound(int c1, int c2, bool side);

	/**
	 * This procedure adds to the matrix inequalities for implications derived when probing a binary variable.
	 *
	 * \param[in] c column index of binary variable having been probed;
	 * \param[in] side if `true`, variable `c` was fixed to `1`; otherwise, to `0`;
	 * \param[out] varBdAdded  number of variable bounds added to constraint matrix;
	 * \param[out] implAdded number of implications generated.
	 * \remark the procedure is used only in rare cases when storing implications in a `CImpl` object requires to much memory.
	 */
	void processImplications(int c, bool side, int &varBdAdded, int &implAdded);

	/**
	 * This procedure is called to further process implications derived when probing a binary variable.
	 *
	 * \param[in] c column index of binary variable having been probed;
	 * \param[in] side if `true`, variable `c` was fixed to `1`; otherwise, to `0`;
	 * \param[in] impl pointer to `CImpl` object that stores implications.
	 * \param[out] varBdAdded  number of variable bounds added to implication matrix;
	 * \param[out] implAdded number of implications added to implication matrix.
	 * \return  number of variable bounds added.
	 */
	void processImplications(int c, bool side, CImpl *impl, int &varBdAdded, int &implAdded);

	/**
	 * The procedure implements probing down or up a given variable.
	 *
	 * \param[in] c column index of binary variable having been probed;
	 * \param[in] side if `true`, variable `c` was fixed to `1`; otherwise, to `0`;
	 * \param[out] varBdAdded number of variable bounds added;
	 * \param[out] ctrTightened number of constraints having been tightened;
	 * \param[out] implAdded number of implications added to implication matrix.
	 * \return `false` if some inconsistencies were detected (constraints are infeasible); otherwise, `true`.
	 */
	bool probeVar(int c, bool side, int &varBdAdded, int &ctrTightened, int &implAdded);

	/**
	 * The procedure substitutes variable `x[k]` for the expression `a*x[j]+b`.
	 *
	 * \param[in] k,j indices of variables;
	 * \param[in] a,b coefficients in expression `a*x[j]+b`.
	 */
	void substituteVar(int k, int j, double a, double b);

	/**
	 * The procedure implements an iterative process of probing all binary (and some integer) variables.
	 *
	 * \param[in] round current preprocessing iterate;
	 * \param[in] timeLimit maximum running time (in seconds);
	 * \param[out] probeVarNum number of variables having been probed;
	 * \param[out] varFixed number of fixed variables;
	 * \param[out] ctrTightened number of modified constraints;
	 * \param[out] varBdAdded number of variable bounds added;
	 * \param[out] implAdded number of implications added to implication matrix.
	 * \return `false` if problem domain is empty.
	 */
	bool probeVars(int round, __LONG timeLimit, int &probeVarNum,
			int &varFixed, int &ctrTightened, int &varBdAdded, int &implAdded);

	/**
	 * The procedure selects a list of variables to be probed.
	 * \param[in] round current preprocessing iterate;
	 * \param[out] ipVars list of `sz` variables to be preprocessed, where `sz` denotes return value.
	 * \return number of variables to be preprocessed.
	 */
	int selectVarsForProbing(int round, int *ipVars);

	/**
	 * The procedure calls `selectVarsForProbing()` to select a list of variables to be probed,
	 * and then calls `probeVars()` to probe variables from that list.
	 * \param[in] timeLimit maximum running time (in seconds);
	 * \param[in] timeLimitPerRound maximum running time (in seconds) for each probing round.
	 * \return `false` if some inconsistencies were detected (constraints are infeasible); otherwise, `true`.
	 */
	bool probing(__LONG timeLimit, __LONG timeLimitPerRound);

protected:
	/**
	 * The procedure displays preprocessing statistics.
	 * \param[in] timeStr string representation of time elapsed since MIPCL have been started;
	 * \param[in] round probing round; if `round < 0`, round field should not be displayed;
	 *            if `round==0`, only header is displayed;
	 * \param[in] probeVarNum number of variables having been probed;
	 * \param[in] varFixed number of variables having been fixed;
	 * \param[in] ctrTightened number of constraints having been tightened;
	 * \param[in] varBdAdded number of new variable bounds having been added;
	 * \param[in] implications number of entries in implication matrix.
	 */
	virtual void probingInfo(char* timeStr, int round,
			int probeVarNum,int varFixed,int ctrTightened,int varBdAdded, int implications);

private:
///////////////////////////////////////////
	/**
	 * Columns are parallel (opposite) if their entries are identical (of different signs)
	 * in all the rows but variable bounds.
	 * \param[in] col1,col2 column indices.
	 * \return `1` if columns are parallel; `-1` if columns are opposite; otherwise, `0`.
	 * \see detectParallelColumns().
	 */
	int areColumnsParallel(int col1, int col2);

	/**
	 * The procedure computes hash values for every matrix column,
	 * and then, for each pair of columns with equal hash values,
	 * calls `areColumnsParallel()` to detect parallel columns.
	 * \return number of columns (variables) deleted.
	 */
	int detectParallelColumns();

	/**
	 * For each pairs of upper bounded mixed-integer inequalities, which sum is an integer inequality,
	 * the procedure preprocesses this sum inequality to tighten the bounds of its (integer) variables.
	 * \return number of bounds having been tightened.
	 */
	int detectRowsWithParallelFracParts();

	/**
	 * The procedure computes implied bounds for monotone variables.
	 * \return number of changed bounds.
	 */
	int processMonotoneVars();

	/**
	 * The procedure _disaggregates_ a given knapsack inequality:
	 *      - sum(i in 0..ctrSize-1) dpVal[i]*x[ipCol[i]] <= b.
	 * In other words, the input inequality is substitutes for a number of inequalities
	 * that together imply the input inequality.
	 *
	 * For example, let us consider the knapsack inequality
	 *     - x_1 + 3x_2 + 4x_3 + 2x_4 <= 8
	 * If we divide it by `k=2` (`k` must be a divider of `b=8`),
	 * we see that the sum of the fractional parts of sll coefficients is not greater than `1`;
	 * therefore, the input inequality is a convex combination of `k` valid constraints
	 *     - x_1 + x_2 + 2x_3 + x_4 <= 4 and 2x_2 + 2x_3 + x_4 <= 4.
	 * Further, the latter inequality is simplified to x_2 + x_3 + x_4 <= 2.
	 *
	 * \param[in] ctrSize number of entries in input inequality;
	 * \param[in] dpVal,ipCol arrays of size `ctrSize`;
	 * \param[in] b right hand side;
	 * \param[in] cpReversed if `cpReversed[j]==1`, than binary variable `x[j]` was complemented.
	 * \return `true` if input constraint has been disaggregated.
	 */
	bool disaggregateKnapsack(int ctrSize, double *dpVal, int *ipCol, double b, char *cpReversed);

	/**
	 * The procedure just calls `disaggregateKnapsack()` for each knapsack constraint.
	 * \return number of _disaggregated_ knapsack constraints.
	 */
	int disaggregateKnapsacks();

	/**
	 * The procedure uses special kind of disjunctions to simplify the constraint matrix.
	 *
	 * \return number of removed non-zero entries.
	 */
	int coverDisjunctions();

	/**
	 * For each integer constraint of size less than 10, the procedure
	 * computes the GCD of all coefficients, and if this GCD does not divide
	 * the right (or left) hand side, the constraint is divided by this GCD and
	 * the right (left) hand size if rounded down (up).
	 * \return number of constraints having been tightened.
	 */
	int roundToMakeStronger();

/////////////////////////////////////////////////////////////
////// Symmetry detection
	/**
	 * The function estimates "connection rates" of all rows with a given set of columns.
	 *
	 * \param[in] sz size of column block;
	 * \param[in] S list (of size `sz`) of columns;
	 * \param[out] d array (of size `m_iM`) to store `row dependencies` with column set given by `S` and `sz`.
	 * \return number of rows `r` in column `v` such that `S[r]=true`.
	 */
	void rowDegFromColSet(int sz, int *S, __LONG *d);

	/**
	 * The function estimates "connection rates" of all columns with a given set of rows.
	 *
	 * \param[in] sz size of rows block;
	 * \param[in] S list (of size `sz`) of rows;
	 * \param[out] d array (of size `m_iM`) to store `column dependencies` with row set given by `S` and `sz`.
	 * \return number of rows `r` in column `v` such that `S[r]=true`.
	 */
	void colDegFromRowSet(int sz, int *S, __LONG *d);

	/**
	 * \param[in] rowPi,colPi permutations of rows and columns.
	 * \return `true` if `pi` is a formulation automorphism.
	 */
	bool isAutomorphism(int *rowPi, int *colPi);

	/**
	 * The procedure partitions column and rows so that
	 * each partition block contains (hopefully) an orbit from the orbit partition
	 * of the formulation automorphism group.
	 *
	 * \param[out] colBlocks number of blocks in column partition;
	 * \param[out] rowBlocks number of blocks in column partition;
	 */
	void buildInitPartition(int &colBlocks, int &rowBlocks);

	/**
	 * The procedure builds a set of generators of the formulation automorphism group.
	 *
	 * \return procedure running time.
	 */
	__LONG detectSymmetry();

	/**
	 *
	 * \param[in] sz number of generators;
	 * \param[in] genList list of `sz` generators;
	 * \param[out] pi,sep orbit partition; for `i=0,...,k-1` (`k` denotes return value),
	 *  `pi[sep[i]],...,pi[sep[i+1]-1]` is orbit `i`;
	 * \param[in] buf memory buffer of size `3*m_iN`.
	 * \return number of orbits.
	 */
	int getOrbitPartition(int sz, int *genList[], int *pi, int *sep, int *buf);

	/**
	 * \return number of variables fixed.
	 */
	int processOrbits(CAUT &aut);

	// Orbital branching is subsumed by orbital fixing !!!?????
	// In case it is used, modify memory allocation for `m_ipFxVar`
	/**
	 * The procedure performs _orbital branching_.
	 * It is applied for the current search tree node, and
	 *    - computes the stabilizer subgroup of the _global symmetry group_;
	 *    - chooses an orbit (of this computed subgroup) that contains a binary variable taking a fractional value;
	 *      the orbit of size `m_iFxNum` is stored in `m_ipFxVar` array, and `m_iBrCol` stores an orbit variable.
	 *
	 * \return `true` if required orbit and variable in it do exist.
	 *
	 */
	bool orbitalBranching();


	/**
	 * This procedure implements _orbital fixing_, which works as follows:
	 *     - select the set of generators that are compatible with the lower bounds at a particular search tree node;
	 *     - if this set is not empty, build the orbit partition for the subgroup induced by these generators;
	 *     - the upper bounds of all variables in each orbit are set to the minimum upper bound value of all variables in this orbit.
	 * It is also worth noting, that orbital fixing subsumes _orbital branching_.
	 *
	 * \return number of changed bounds.
	 */
	int orbitalFixing();

////////////////////////////////

	/**
	 * The procedure just calls `makeCoeffIntegral()` for each row.
	 * \return always `true`.
	 */
	bool preprocessInit() final;

	/**
	 * The procedure reformulate a MIP by applying a number of preprocessing techniques (including probing).
	 * \return `false` in case the problem being solved has been proven to be inconsistent.
	 */
	bool preprocessPlus() final;

	void BranchAndCutWorker(); ///< This procedure implements a branch-and-cut algorithm.

	/**int CMIP::isEntry(int col, int row)
	 *
	 * When `BranchAndCut()` terminates, it calls `fixSolutionState()`
	 * to assign a value to `m_iState`, which is the bitwise OR of the members of the enumeration `enProbState`.
	 */
    void fixSolutionState();
	//////
public:
    /**
     * This procedure does a lot of useful things: allocates memory, applies different preprocessing techniques,
     * scales the matrix, and so on. In many cases, we can get substantially simpler (for solving)
     * problem, which is still equivalent to the original (user) problem.
     * The latter means that, given an optimal solution to the transformed problem, we can
     * easily compute an optimal solution to the original problem.
     * \return `false` if inconsistency has been detected; otherwise, `true`.
     * \throws CMemoryException lack of memory.
    */
    bool prepare() final;

    /**
     * The function resumes the state in which the problem was before the solution procedure has started.
     */
    void reset() final;

    /**
     * The procedure creates a number of threads each of which calls `BranchAndCutWorker()`
     * that implements the branch-and-cut/price algorithm.
     * Each instance of BranchAndCutWorker() stops if duality (integrality) gap
     * (i.e., difference between upper and lower bounds on the optimal objective value) is less than the value of `gap`.
     * \param[in] timeLimit limit on solution time (in seconds); `timeLimit==0l` means no time limitation;
     * \param[in] gap duality gap.
     * \sa optimize(), setDualGap().
     */
    void BranchAndCut(__LONG timeLimit=0l, double gap=0.0);

    /**
     * The procedure solves the root LP, and then, if its solution is not integral, calls `BranchAndCut()`.
     * \param[in] timeLimit limit on solution time (in seconds); `timeLimit==0l` means no time limitation;
     * \param[in] gap integrality gap;
     * \param[in] solFileName pointer to string with file name for storing intermediate solutions; if `solFileName=0`, intermediate solutions are not stored.
     * \sa CLP::optimize(), BranchAndCut(), setDualGap().
     */
    void optimize(__LONG timeLimit=0l, double gap=0.0, const char *solFileName=0) final;

    /**
     * The function returns two pointers to the internal __MIPCL__ arrays storing the solution found so far.
     * \param[out] dpX,ipHd `dpX[j]` is the value of variable whose handle is `ipHd[j]`, `j=1,...,n`,
     *    where `n` is the return value. Do not modify the values of both arrays  `dpX` and `ipHd`.
     * \return number of variables.
     * \sa CLP::getSolution().
     */
    int getSolution(double* &dpX, int* &ipHd);

    /**
     * \return `true` if at least one feasible solution has been found so far.
     */
    bool isSolution() const final;

    /**
     * First, `isSolution()` is called; if it returns `true`, then
     * we can call this function to verify whether
     * the solution found by __MIPCL__ is optimal.
     * \return `true`, if an optimal solution has been found.
     */
    bool isSolutionOptimal() const
	{return (m_iState & PROB_OPTIMAL)? true: false;}

    /**
     * \return `true` if solution procedure stopped after exceeding given time limit.
     */
	bool timeLimitStop() const
	{return (m_iState & PROB_TIME_LIMIT)? true: false;}

    /**
     * \return currently upper bound on the objective value defined to be
     * the largest value of the nodes LPs among all search tree leaves.
     */
    double getObjBound() const;

	/**
	 * The procedure writes MIP solutions into the file.
	 * The user can overload this function to store solutions in an appropriate way.
	 * \param[in] fileName name of the file to store solutions; if `fileName=`0`,
	 * the solver makes up the file name by appending the extension ".sol" to the name of the problem being solved.
	 */
    virtual void printSolution(const char* fileName=0) override;

protected:
    /**
     * To generate problem specific cuts, __MIPCL__ allows the user to overload any of three virtual functions:
     * `separate()` (a member of CLP class), `genCut1()`, and `genCut2()`.
     * The solver in turn calls `separate()`, `genCut1()`,
     * its internal cut generating procedure, and then `genCut2()`.
     * Both, `genCut1()` and `genCut2()` are called only if their predecessors failed to produce any cut.
     *
     * Rule of Thumb:
     *  - `separate()` produces all the inequalities that constitute a part of the problem formulation
     *    (for example, TSP cut inequalities);
     *  - strong (facet defining and etc.) cuts are generated by `genCut1()`
     *    (for example, TSP comb inequalities);
     *  - less stronger cuts must be generated by `genCut2()`.
     *
     * \param[in] n number of variables;
     * \param[in] dpX,ipColHd arrays of size `n` that represent a solution to cut off;
     *     `dpX[i]` is the value of the variable with handle `ipColHd[i]`.
     * \return `true` if at least one cut has been detected
     * \sa  `genCut2()`, `CLP::separate()`.
     */
    virtual bool genCut1(int n, const double* dpX, const tagHANDLE* ipColHd);

    /**
     * Similar to `genCut1()`, except that `genCut2()` is called only if all the other cut generating
     * procedures failed to produce any cut.
     * \sa genCut1().
     */
    virtual bool genCut2(int n, const double* dpX, const tagHANDLE* ipColHd);

    /**
     * This function must be overloaded in any derived class
     *  that implements its own pool for storing generated inequalities.
     * When a new node is added to the search tree,
     *  `lockCtr()` is called for any constraint
     *  generated by the derived class cut generating procedures (i.e. `separate()`, `genCut1()`, and `genCut2()`)
     *  unless that constraint was sent to the __MIPCL__ pool.
     *  Properly implemented two functions, `lockCtr()` and `unlockCtr()`, allow derived classes
     *  to delete, if necessary, unused cuts from their pools.
     *  \param[in] hd constraint handle.
     *  \sa unlockCtr(), delNodeLocalCtrs().
     */
    virtual void lockCtr(tagHANDLE hd)
        {}

    /**
     * This function must be overloaded in any derived class
     *  that implements its own pool for storing generated inequalities.
     * When a new node is extracted from the search tree,
     *  `unlockCtr()` is called for any constraint
     *  generated by the derived class cut generating procedures (i.e. `separate()`, `genCut1()`, and `genCut2()`)
     *  unless that constraint was sent to the __MIPCL__ pool.
     *  Properly implemented two functions, `lockCtr()`, `unlockCtr()`, and `delNodeLocalCtrs()` allow derived classes
     *  to delete, if necessary, unused cuts from their pools.
     *  \param[in] hd constraint handle.
     *  \sa lockCtr(), delNodeLocalCtrs().
     */
	virtual void unlockCtr(tagHANDLE hd)
        {}

    /**
     * This function must be overloaded in any derived class
     *  that implements its own pool for storing generated columns.
     * When a new node is added to the search tree,
     *  `lockColumn()` is called for any column
     *  generated by the derived class column generating procedure (i.e., generateColumns()
     *  unless that column was sent to the __MIPCL__ pool.
     *  Properly implemented two functions, `lockColumn()`, `unlockColumn()`, and `delNodeLocalColumns()`
     *  allow derived classes to delete, if necessary, unused columns from their pools.
     *  \param[in] hd column handle.
     *  \sa unlockColumn(), delNodeLocalColumns().
     */
    virtual void lockColumn(tagHANDLE hd)
        {}

    /**
     * This function must be overloaded in any derived class
     *  that implements its own pool for storing generated columns.
     * When a new node is extracted from the search tree,
     *  `unlockColumn()` is called for any column
     *  generated by the derived class column generating procedure (i.e., `generateColumns()`
     *  unless that column was sent to the __MIPCL__ pool.
     *  Properly implemented two functions, `lockColumn()` and `unlockColumn()`, and `delNodeLocalColumns()`
     *  allow derived classes to delete, if necessary, unused columns from their pools.
     *  \param[in] hd column handle.
     *  \sa lockColumn(), delNodeLocalColumns().
     */
    virtual void unlockColumn(tagHANDLE hd)
        {}

    /**
     * This function must be overloaded in any derived class
     *  that implements its own pool for storing generated inequalities.
     * When a node is pruned, `delNodeLocalCtrs()` is called to allow the derived class
     * to eliminate, if necessary, all the constraints that are local for node `nd`.
     * \param[in] nd pruned node.
     * \sa lockCtr(), unlockCtr().
     */
	virtual void delNodeLocalCtrs(int nd) {};

	/**
     * This function must be overloaded in any derived class
     *  that implements its own pool for storing generated columns.
     * When a node is pruned, `delNodeLocalColumns()` is called to allow the derived class
     * to eliminate, if necessary, all the columns that are local for node `nd`.
     * \param[in] nd pruned node.
     * \sa lockColumn(), unlockColumn().
     */
	virtual void delNodeLocalColumns(int nd) {};

//////////////////////////////////////////////////////////////////////////
//                    SOS feasibility and branching                     //
//////////////////////////////////////////////////////////////////////////
private:
    /**
     * \param[out] row violated SOS1-constraint;
     * \param[out] l position in row `row`; in left branch all variables in positions preceding `l` are set to `0`,
     * while in right branch  all variables in positions following `l` are set to `0`.
     * \return `true` if violated SOS1-constraint has been detected.
     */
    bool SOS1_Branching(int& row, int& l);

    /**
     * The procedure decides on which SOS2-constraint to perform branching.
     * \param[out] l if return value, say `r` is non-negative,
     *  `l` is position in row `r`; in left branch all variables in positions preceding `l` are set to `0`,
     * while in right branch  all variables in positions following `l` are set to `0`.
     * \return if non-negative, index of violated SOS2-constraint; otherwise, all SOS2-constraints are satisfied.
     */
    int SOS2_Branching(int& l);

    /**
     * \param[in] X current node LP solution;
     *    it is assumed that `X[j]` is the value of variable `j` in the preprocessed problem.
     * \returns  `true` if all SOS1 and SOS2-constraints are satisfied, and `false` otherwise.
     */
    bool isSOSFeasible(double* X);

    /**
     * The function picks up a free node in the search tree.
     * \param[out] height of node being returned.
     * \return node index.
     */
    int getNode(int &height);

protected:
    /**
     * The procedure decides how to branch at the node being processed.
     * It is called by the solver when it starts processing a new node.
     * Overloading `startBranching()` together with updateBranch(), allows the user
     * to take full control on the branching process.
     * \param[in] nodeHeight height (in the search tree) of the currently processed node.
     * \return number of branches.
     * \sa updateBranch().
     */
    virtual int startBranching(int nodeHeight);

    /**
     * The function sets new bounds for integer variables or/and adds new constraints to the matrix
     * to transform the parent node problem into the problem of branch `i`.
     * \param[in] i branch index, when `i` is nonnegative integer that is less than the value returned by `startBranching()`.
     * \return `false` if inconsistency has been detected; otherwise, `true`.
     * \sa startBranching().
     */
    virtual bool updateBranch(int i);

    /**
     * There are practical applications with some requirements that are difficult to be expressed
     * by means of inequalities. As alternative, we can implement
     * a special branching strategy that eliminates all the solutions not valid for those difficult constraints.
     * In such a case, we need to overload three functions: `startBranching()`, `updateBranch()`, and `isFeasible()`.
     *
     * Any time when a new solution has been computed that satisfies all the constraints known to the solver,
     * `isFeasible()` is called to verify whether that solution is really feasible for the problem being solved.
     * \param[in] n number of variables;
     * \param[in] dpX,ipColHd for `j=1,...,n`, `dpX[j]` is the value of variable with handle `ipColHd[j]`.
     * \return `true` if the solution given by `dpX` and `ipColHd` is feasible; otherwise, `false`.
     * \sa startBranching(), updateBranch().
     */
    virtual bool isFeasible(int n, const double* dpX, const tagHANDLE* ipColHd)
        {return true;}

private:
    /**
     * The function first calls `updateBranch()` and then tries to propagate changes done by `updateBranch()`.
     * \param[in] i branch number (numbering starts from `0`).
     */
    bool _updateBranch(int i);

#ifdef __PSEUDOCOSTS_
    /**
     * The procedure is called from `strongBranching()` to build a list of integral variables with fractional values.
     *     - If there are variables on which branchings have been done more that 4 times, and if during those branchings
     *       either the objective had increased or inconsistency had been detected,
     *       one of such variables with maximum score (computed based on pseudocosts) is selected for branching.
     *     - If for all variables on which branchings have been done more that 4 times,
     *       neither the objective had increased nor inconsistency had been detected,
     *            + a list of variables on which branchings have been done less that 4 times
     *              is generated;
     *            + if for all variables branchings have been done more that 4 times,
     *              one of variables with maximum score (computed based on the current basis) is selected for branching.
     * \param[in] maxSize maximum number of variables in return list;
     * \param[out] ipLst list of `k` integral variables taking fractional value, where `k` is return value.
     * \returns size of returned list.
     */
    int selectFracListUsingPseudocosts(int maxSize, int *ipLst);

    /**
     * The procedure does not use pseudocosts to build a list of fractional variables.
     * In particular, this procedure is used in place of `selectFracListUsingPseudocosts()`
     * when columns are generated.
     * \param[in] maxSize maximum number of variables in return list;
     * \param[out] ipLst list of `k` integral variables taking fractional value, where `k` is return value.
     * \returns size of returned list.
     */
    int selectFracListWithoutPseudocosts(int maxSize, int *ipLst);
#endif

    /**
     * The procedure is used in `strongBranching()` to build a list of fractional variables.
     * \param[in] maxSize maximum number of variables in return list;
     * \param[out] ipLst list of `k` integral variables taking fractional value, where `k` is return value.
     * \returns size of returned list.
     */
	int selectFracList(int maxSize, int *ipLst);

	/**
	 * The procedure implements strong branching.
	 *
	 * First, a subset of variables with biggest estimates (scores) are selected by calling `selectFracList()`.
	 * Then, the values of these selected variables are rounded down and up
	 * to compute a weighted decreases of the objective after a prescribed number of dual iterates;
	 * a variable with the biggest decrease is chosen to do branching.
	 *
	 * \param[in] nodeHeight height of processed node.
	 * \return `0` if everything was OK.
	 */
	int strongBranching(int nodeHeight);

#ifdef __PSEUDOCOSTS_
	/**
	 * An implementation of `getFractional()` when pseudocosts are used.
	 * \param[out] d current value of variable `j` which is return value.
	 * \return variable `j` to split its domain into two parts, `[m_dpD[j<<1],d]` and `[d+1,m_dpD[(j<<1)+1]]`
     * if `j` is negative, then no fractional variable has been detected.
	 */
	int getFractionalUsingPseudocosts(double& d);

	/**
	 * An implementation of `getFractional()` when pseudocosts are not used.
	 * \param[out] d current value of variable `j` which is return value.
	 * \return variable `j` to split its domain into two parts, `[m_dpD[j<<1],d]` and `[d+1,m_dpD[(j<<1)+1]]`
     * if `j` is negative, then no fractional variable has been detected.
	 */
	int getFractionalWithoutPseudocosts(double& d);
#endif

protected:
	/**
	 * The procedure selects a variable to branch on it.
	 * \param[out] d current value of variable `j` which is return value.
	 * \return variable `j` to split its domain into two parts, `[m_dpD[j<<1],d]` and `[d+1,m_dpD[(j<<1)+1]]`
     * if `j` is negative, then no fractional variable has been detected.
	 */
    virtual int getFractional(double& d);

    ///< User defined prime heuristic.
    /**
     * The function is used to round a non-integer solution stored in `dpX` and `ipColHd`.
     * \param[in] n size of `dpX` and `ipColHd` (number of variables);
     * \param[in,out] dpX `dpX[i]` is the value of the variable with handle `ipColHd[i]`;
     * \param[in] ipColHd stores handles of the variables;
     * \param[in,out] objVal on input, upper (for maximization problem) or lower (minimization) bound on objective value;
     * \return `true` if the procedure succeeded in constructing
     * a feasible solution ( stored in `dpX` and `ipColHd`) with objective value better than input value of `objVal`.
     */
    virtual bool roundSolution(double& objVal, int n, double* dpX, const int* ipColHd)
    {return false;}

public:
    /**
     * The function specifies a way of rounding optimal LP solutions.
     * \param[in] dir direction of rounding fractional components;
     * if `dir==ROUND_USER`, `roundSolution()` with four argument has to
     * be overloaded to implement problem specific rounding strategy.
     * \sa enRoundType, getRoundingType().
     */
    void setRoundingType(enRoundType dir)
        {m_eRoundType=dir;}

    /**
     * \return type (direction )of rounding fractional components.
     * \sa enRoundType, setRoundingType().
     */
    enRoundType getRoundingType() const
        {return m_eRoundType;}
private:
    /**
     * The function is called in` roundSolution()` to round the LP solution of currently processed node.
     * \param[out] n size of solution (+ variables from propagation stack) stored in `m_dpUd`.
     * \param[out] objVal objective value of new solution found.
     * \returns `true` if a new record solution has been  found .
     */
    bool roundX(int& n, double& objVal);

    /**
     * The function calls a user defined `roundSolution()` (with four arguments)
     * if `m_eRoundType == ROUND_USER`, otherwise, `roundX()` is called.
     * \param[out] objVal objective value of solution found (if any);
     * \param[out] n number of variables.
     * \sa `enRoundType`, `setRoundingDir()`, `roundX()`.
     */
    bool roundSolution(double& objVal, int& n);
    
    /**
     * The procedure uses dual variables to change bounds of the prime variables.
     * \return number of fixed binary variables.
     */
    int propagateByDualVars();

    /**
     * The procedure uses reduced costs and shadow prices to reduce lower and upper bounds
     * for variables and constraints when an upper bound on the objective is known.
     * The procedure is used only in `propagate()`.
     * \param[in,out] dpD array of size `2*m_iN`, where `dpD[j<<1], dpD[(j<<1)+1` are, respectively,
     * lower and upper bounds on variable `j`;
     * \param[in,out] dpB array of size `2*m_iM`, where `dpD[r<<1], dpD[(r<<1)+1` are, respectively,
     * lower and upper bounds for constraint `r`.
     *
     */
    void propagateByDualVars(double *dpD, double *dpB);
    
protected:
    /**
     * The user can overload this function to implement some problem specific propagation algorithm.
     * To access variable bounds, one can call the following functions:
     *  `getVarNum()`, `getVarHd()`, `getVarLoBound()`, `getVarUpBound()`.
     * Do not change variable bounds directly correcting values in `m_dpD`, instead,
     * use the functions `setVarLoBound()`, `setVarUpBound()`, or `setVarBounds()`.
     * \return `false` if inconsistency has been detected; otherwise, `true`.
     * \sa CLP::setVarLoBound(), CLP::setVarUpBound(), CLP::setVarBounds().
     */
    virtual bool propagate();
      
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////
// D I V I N G
///////////////////////
private:
    /**
     * The procedure is called in `dive()` to return a violated SOS1-constraint to be rounded.
     *
     * \param[in] maxRows maximum number of SOS1-constraints to return;
     * \param[out] ipRow list of SOS1-constraints to be rounded off;
     * \param[in] dpVal memory buffer of size `m_iM`;
     * \return number of SOS1-constraints to round off (size of array `ipRow`).
     */
    int getSOS1toRoundOff(int maxRows, int* ipRow, double* dpVal);

    /**
     * The procedure is called in `dive()` for setting to zero all but 1 variables in a given SOS1-constraint.
     * \param[in] row index of SOS1-constraint.
     * \return number of variables fixed.
     */
    int SOS1_roundOff(int row);

    /**
     * The procedure is called in `dive()` to return a violated SOS2-constraint to be rounded.
     *
     * \param[in] maxRows maximum number of SOS2-constraints to return;
     * \param[out] ipRow list of SOS2-constraints to be rounded off;
     * \param[in] dpVal memory buffer of size `m_iM`;
     * \return number of SOS1-constraints to round off (size of array `ipRow`).
     */
    int getSOS2toRoundOff(int maxRows, int *ipRow, double *dpVal);

    /**
     * The procedure is called in `dive()` for setting to zero all but 2 successive variables in a given SOS2-constraint.
     * \param[in] row index of SOS2-constraint.
     * \return number of variables fixed.
     */
    int SOS2_roundOff(int row);

    /**
     * The procedure is called in `dive()` when inconsistency is detected.
     * It tries to prevent inconsistency  detected in `dive()`
     * by reverting the value of a single fixed binary variable.
     * \param[in] maxObjDec only variables having reduced costs not greater than `maxObjDec` are considered;
     * \param[in,out] ipChanged array of size `m_iN`, where `ipChanged[i]` is number of times the value of binary `i` has been changed.
     * \return `true` in case of success.
     */
    bool DIVE_backTrace(double maxObjDec, int *ipChanged);

    /**
     * This procedure is similar to `DIVE_backTrace()` but its goal is
     * to change fixings of binary variable to reduce the objective of the LP being solved
     * at least by a given value.
     * \param[in] gap value by which objective be reduced;
     * \param[in,out] ipChanged array of size `m_iN`, where `ipChanged[i]` is number of times the value of binary `i` has been changed.
     * \return `true` in case of success.
     */
    bool DIVE_worthObjBackTrace(double gap, int *ipChanged);

    /**
     * For an objective of small size (<= 150), adds the constraint `loBound <= c^Tx <= upBound` to the matrix.
     * \param[in] type type of constraint added;
     * \param[in] loBound lower objective bound;
     * \param[in] upBound upper objective bound.
     * \return index of added row.
     */
    int addObjectiveCtr(int type, double loBound, double upBound);

    /**
     * The procedure is called within `dive()` to fix the values of those integral
     * variables that are closed to their lower or upper bounds.
     * \param[in,out] ipChanged array of size `m_iN`, where `ipChanged[i]` is number of times the value of binary `i` has been changed.
     * \return if non-negative, number of variables fixed, or `-1` in case of failure.
     * \remark procedure is used only in `dive()`.
     */
    int DIVE_fixVars(int *ipChanged);

    /**
     * The procedure is iteratively fixing an increasingly large number of variables until either
     *      - a new, best integral solution is found, a new incumbent,
     *      - or the fixings make the result LP infeasible,
     *      - or an objective value worse than the current incumbent.
     * \param[in] timeLimit limit on running time (in seconds).
     * \return `false` if the node was proven to be infeasible; otherwise, `true`.
     */
    bool dive(__LONG timeLimit);
// end of D I V I N G ///////////////////////////////

protected:
    /**
     * The function is called to prepare information for storing in the node of the search tree.
     * This function can be overloaded to exploit problem specifics and store node data more efficiently.
     * \param[out] ipVal array to write down node data.
     * \return number of bytes stored in `ipVal`.
     * \sa restoreNodeData().
     */
    virtual int storeNodeData(int* ipVal);

    /**
     * Taking as input the data stored by `storeNodeData()`,
     * the function must restore the node to be processed.
     * \param[in] ipVal data needed to restore the node.
     */
    virtual void restoreNodeData(int* ipVal);

    /**
     * The procedure is called any time a group of constraints (inequalities) are removed from the matrix.
     * This procedure must be overloaded in any derived class that generates inequalities and stores them in its own pool.
     * \param[in] sz number of constraints;
     * \param[in] ipCtr list of `sz` columns.
     */
    virtual void setCtrsInactive(int sz, int* ipCtr);

    /**
     * The procedure is called any time a group of columns are removed from the matrix.
     * This procedure must be overloaded in any derived class that generates columns and stores them in its own pool.
     * \param[in] sz number of columns;
     * \param[in] ipCtr list of `sz` constraints.
     */
    virtual void setColumnsInactive(int sz, int* ipCtr);

private:
    /**
     * When a new node starts processing,
     * the procedure adds constraints from the pool to the constant part (comprising first `m_iM0` rows) of the matrix.
     * \throws CMemoryException in multithreaded applications.
     */
    void restoreMatrix();

	/**
	 * The procedure is used only in `preprocessObjective()`.
	 * \param[in] j variable index;
	 * \param[in,out] dpD array of size `2*m_iN`, where `dpD[i<<1]` and `dpD[(i<<1)+1]`
	 * are lower and upper bounds on variable indexed by `i`;
	 * \param[in,out] objUpBound upper bound on objective value.
	 */
	void propagateVarBound(int j, double *dpD, double &objUpBound);

    /**
     * The procedure preprocesses the constraint `loBd <= c^Tx` (defined by the objective)
     * to reduce lower and upper bounds of variables.
     * \param[in] loBd current lower bound value.
     */
    void preprocessObjective(double loBd);

    /**
     * The procedure restores the LP for a given node in the search tree.
     * \param[in] iNode search-tree node to be restored;
     * \param[in] flag if `false`, only basis is restored.
     */
    void restoreSubproblem(int iNode, bool flag=true);

    /**
     * For a given node in the search tree, the procedure "unlocks" cuts and columns from the pool
     * that were used by the node. Unlocking cut (or column) means decrementing the counter of users (that are nodes)
     * of this cut (or column).
     * \param[in] ind index of node in search tree.
     */
	void releaseNode(int ind);

	/**
	 * The procedure deletes from the pool all the cuts that are local for a given search tree node.
	 * \param[in] ind index of search tree node.
	 */
	void delNodeLocals(int ind);

	/**
	 * The procedure removes the currently processed node from the search tree
	 * and deletes all local cuts generated at this node.
	 */
	void freeCurrentNode();

    /**
     * The procedure resets some parameters for a new node of the search tree.
     */
    void resetNodeParam();

    /**
     * This procedure adds a given constraint to the pool.
     * \param[in] i constraint index.
     */
    void addCtrToPool(int i);

    /**
     * This procedure adds a given column to the pool.
     * \param[in] col column index.
     */
    void addColumnToPool(int col);

    /**
     *  The procedure is a truncated version of `solveNodeLP()` without cut generation and deleting rows.
     *  \param[in] loBd procedure stops running if objective value becomes less or equal `loBd`.
     *  \return return code, see `CLP::dualSimplex()`.
     */
    int resolveLP(double loBd);

    /**
     * The procedure is called in `processNode()` to solve the LP for the processing node.
     * \return return code, see `CLP::dualSimplex()`.
     */
    int solveNodeLP();

#ifdef __PSEUDOCOSTS_
    /**
     * The procedure is called in `processNode()` and `_updateBranch()` to update `m_dpBrDec` and `m_ipBrNum` arrays.
     * \param[in] j column (variable) index;
     * \param[in] side if `true`, variable `j` was rounded up; otherwise, variable `j` was rounded down;
     * \param[in] bFeasible if `true`, node LP remained consistent after branching; otherwise, node LP became inconsistent;
     * \param[in] objInc objective increase.
     */
    void updateBranchingStats(int j, bool side, bool bFeasible, double objInc);
#endif

    /**
     * This procedure solves node LPs.
     */
    void processNode();

    /**
     * This procedure solves the root node LP.
     */
    void processRootNode();
    
    /**
     * The procedure verifies whether a given solution is integral.
     * \param[in] dpX array of size `m_iN`; `dpX[i]` is value of variable having handle `m_ipColHd[i]`.
     * \return `true` if solution `dpX` is integral.
     */
    bool isIntegral(double* dpX);

    /**
     * \return number of fractional components in the LP solution of the currently processed node.
     */
    int getFracVarNum() const;

    /**
     * The function does the same as `CLP:getPrimeSolution()` does, plus rounds off integer components.
     * \param[in] dpBasicX `dpBasicX[m_ipBasicColumn[i]]` value of basic variable `m_ipBasicCol[i]`, `i=0,...,m_BasisSize-1`;
     * \param[out] dpX array of size `n`, where `n` is return value, and `dpX[j]` is value of variable with handle `m_ipHd[j]`.
     * \return number of variables in the original (non-preprocessed) problem.
     */
    int getLpSolution(double* dpBasicX, double* dpX);

	/**
	 * The function is a wrapper to call `changeRecord()`, which may be overloaded in derived classes.
	 * \sa changeRecordVal().
	 */
	void _changeRecord(double objVal,
				int n, const double* dpX, const tagHANDLE* ipHd);
public:
	/**
	 * When a new solution is found which is better than the record solution,
	 *  `changeRecord()` procedure is called to replace the record solution for a new one.
	 *  Sometimes, e.g. when an initial solution is known, it is needed to change only the bound
	 *  on the objective value. In such cases, `changeObjBound()` should be called.
	 * \param[in] objVal new lower (for maximization problem) or upper (for minimization problem) bound,
	 *  which is usually objective value of best solution found so far.
	 * \sa changeRecord().
	 */
	void changeObjBound(double objVal);

	/**
	 * The best solution found so far is called a _record solution_, and its objective value is known as a _record_,
	 * which is a _lower_ or _upper_ bound on the optimal objective value.
	 * The procedure is called when a new best solution is found to store it as the record solution.
	 *  When solving structured problems (such as TSP), the user can overload `changeRecord()`
	 *  to store record solutions more efficiently.
	 * \param[in]  objVal objective value;
	 * \param[in] n number of variables;
	 * \param[in] dpX,ipHd solution, `dpX[j]` is value of variable with handle `ipHd[j]`, `j=1,...,n`.
	 * \sa changeObjBound().
	 */
	virtual void changeRecord(double objVal,
		int n, const double* dpX, const tagHANDLE* ipHd);

///////////////////////////////////////////////////////////////
//               CUT Generation
///////////////////////////////////////////////////////////////
protected:
	/**
	 * The procedure reduces size of the cut by removing some variables taking their bound values;
	 * this may make the cut local.
	 * \param[in] maxSize maximally allowed cut size;
	 * \param[in,out] sz size of arrays `ipCol` and `dpVal` (size of cut);
	 * \param[in,out] ipCol,dpVal,b inequality `sum(i=0..sz-1) dpVal[i]*x(ipCol[i]) <= b` represents both, input and output, cuts;
	 * \param[out] iFactor `sum(i=0..sz-1) r(i)`, where `frexp(fabs(m_dpVal[i]),r(i))`;
	 * \param[in,out] bLocal is `true` if resulting cut is local.
	 */
	void shrinkCut(int maxSize, int &sz, int *ipCol, double *dpVal, double &b, int &iFactor, bool &bLocal);

private:
	/**
	 * The procedure preprocesses a constraint previously generated as a cut.
	 * \param[in] row constraint index.
	 * \return number of fixed variables.
	 */
	int preprocSingleCtr(int row);

	/**
	 * The procedure preprocesses all the constraints previously generated as cuts.
	 * \param[in] r1 index of first constraint to be preprocessed;
	 * \param[in] r2 index of last constraint to be preprocessed;
	 * \return number of fixed variables.
	 */
	int preprocCtrs(int r1, int r2);

	/**
	 * The procedure extracts from the pool all the cuts that are valid for the currently processed node,
	 * and that are violated by the optimal solution found for this node.
	 *
	 * \return number of cuts extracted.
	 */
    int separateFromPool();

    /**
     * The procedure generate click cuts based on the implication matrix (a member of the `CImpl` class)
     * that is build when calling `probing()`.
     * \return number of cuts generated.
     */
	int genImpliedClicks();

///////////////////////////////////////////////////////////////////
// Variable Bounds
//////////////////////
    /**
     * The procedure generate variable lower and upper bound cuts stored in `m_pImpl`.
     * \return number of cuts generated.
     */
	int genImplVarBounds();

    ////////////////////////////////////////////////////////////////////
//  K N A P S A C K  cuts
//////////////////////////
// !!!m_dpArray and m_ipArray are used when building cuts
    /**
     * The procedure implements an algorithm for lifting cover inequalities.
     * For the idea of lifting, see the description of `LCI()`.
     */
	void liftCover(int& sz, int n, int *pi,
		double b, double* dpA, int* ipAlpha);

	/**
	 * To find a cover inequality that separates a given point from a given knapsack set,
	 * `LCI()` solves a 0,1-knapsack problem; this computation may be time consuming.
	 * For knapsack sets of big size, `simpleLCI()` solves that knapsack problem approximately
	 * but much faster.
	 * \param[in] sz,dpA,ipCol,b input knapsack inequality `sum(i=0,...,sz-1) dpA[i]*x[ipCol[i]] <= b`;
	 * \param[in] dpX `dpX[j]` is value of variable `j`;
	 * \param[in] bLocal if `true`, input inequality is local;
	 * \param[out] rhs right hand side of returned inequality;
	 * \param[out] type type of returned inequality (bitwise OR of `enVarType` members);
	 * \return size of generated cut (let us denote it by 'cutSize)`, or '0' if no cut has been found.
	 * \remark The returned inequality (cut) is given as
	 *     'sum(i in 0..cutSize-1) m_dpArray[i]*x(m_ipArray[i]) <= rhs'.
	 */
	int simpleLCI(int sz, double* dpA, int* ipCol, double* dpX, double b, bool bLocal, double &rhs, unsigned &type);

	/**
	 * The procedure builds LCI cut based on the mixed knapsack
	 *      'sum(i in 0..sz-1) dpA[i]*x(ipCol[i]) <= b'.
	 * \param[in] sz size of mixed knapsack;
	 * \param[in] dpA constraint coefficients, array of size 'sz';
	 * \param[in] ipCol column indices, array of size 'sz';
	 * \param[in] dpX array of size `sz`, `dpX[i]` is value of variable `ipCol[i]`;
	 * \param[in] bLocal 'true` if input constraint is local; otherwise, 'false';
	 * \param[in,out] b right hand side value;
	 * \param[out] type type of returned constraint (bitwise OR of `enVarType` members).
	 * \return cut size, or '0' if no cut has been found.
	 */
	int mxSimpleLCI(int sz, double* dpA, int* ipCol, double* dpX, bool bLocal, double &b, unsigned &type);

	/**
	 * Given a knapsack inequality, the procedure tries to build
	 * a lifted gub-cover inequality that is violated by a given point.
	 * \param[in] n number of variables;
	 * \param[in] n0 number of variables taking value of `0`;
	 * \param[in] n1 number of variables taking value of `1`;
	 * \param[in] nf number of variables taking value from `(0,1)`;
	 * \param[in] A vector of coefficients of input knapsack inequality;
	 * \param[in] X point to be separated;
	 * \param[in] b right-hand side of input inequality;
	 * \param[out] iBeta right-hand side of output inequality;
	 * \param[out] ipAlpha vector of coefficients of output inequality.
	 * \return `true` if separating inequality has been found, and `false` otherwise.
	 */
	bool LGCI(int n, int n0, int n1, int nf,
			double* A, double* X, double b, int &iBeta, int *ipAlpha);

	/**
	 * For a given knapsack inequality, the procedure builds a lifted cover inequality.
	 *
	 * For the knapsack set `X={x in {0,1}^n: sum(j=1,...,n) a(j)*x(j) <= b}` with all `a(j)` greater than zero,
	 * a subset `I` from `{1,...,n}` is a _cover_ if `sum(j \in I) a(j) > b`;
	 * then the _cover_ inequality `sum(j in I) x(j) <= |I|-1` is valid for `X`.
	 * One can strengthen this cover inequality by _lifting_ it.
	 *
	 * The idea of lifting is as follows.
	 * Given an inequality `sum(j=1,...,n) u(j)*x(j) <= q` (some `u(j)` may be zeros) that are valid for `X`,
	 * if `sum(j=1,...,n) u(j)*x(j) <= s < q` for all `x in X` such that `x(j)=1`,
	 * then we can increase `u(j)` by `q-s`.
	 *
	 * A lifting algorithm is implemented in `liftCover()`.
	 *
	 * \param[in] i constraint index;
	 * \param[in] side if `true`, right hand side inequality of constraint `i` is used; otherwise,
	 * left hand side inequality is used;
	 * \param[out] rhs right hand side of returned inequality;
	 * \param[out] type type of returned inequality (bitwise OR of `enVarType` members).
	 * \return size of generated cut (let us denote it by 'cutSize)`, or '0' if no cut has been found.
	 * \remark The returned inequality (cut) is given as
	 *     'sum(i in 0..cutSize-1) m_dpArray[i]*x(m_ipArray[i]) <= rhs'.
	 * \sa liftCover(), GUB_LCI().
	 */
	int LCI(int i, bool side, double &rhs, unsigned &type);

	/**
	 * For a given constraint, the procedure builds a lifted GUB-cover inequality.
	 *
	 * Given a `X={x in {0,1}^n: sum(j=1,...,n) a(j)*x(j) <= b; sum(j\in J_i) x(j) <=1 for i=1,...,k}`,
	 * where all `a(j)` are greater than zero and `J_1 U J_2 U ... U J_k = {1,...,n}`.
	 * A subset `I` from `{1,...,n}` is a _GUB-cover_ if `sum(i=1,...,k) q(i) > b`,
	 * where `q(i)=max{a(j): j in I and j in J_i} if the intersection of `I` and `J_i` is not empty;
	 * otherwise, `q(i)=0`. Then the GUB-cover inequality
	 *     `sum(j in I) x_j <= (sum(i=1,...,k) sign(q(i))) -1`
	 * is valid for `X`.
	 *
	 * \param[in] i constraint index;
	 * \param[in] side if `true`, right hand side inequality of constraint `i` is used; otherwise,
	 * left hand side inequality is used;
	 * \param[out] rhs right hand side of returned inequality;
	 * \param[out] type type of returned inequality (bitwise OR of `enVarType` members).
	 * \return size of generated cut (let us denote it by 'cutSize)`, or '0' if no cut has been found.
	 * \remark The returned inequality (cut) is given as
	 *     'sum(i in 0..cutSize-1) m_dpArray[i]*x(m_ipArray[i]) <= rhs'.
	 * \sa LCI().
	 */
	int GUB_LCI(int i, bool side, double &rhs, unsigned &type);

	/**
	 * The procedure implements an iterator through the list of knapsack constraints.
	 * The iterator is used in `knapsackCuts()`.
	 * \param[in] row index of currently used knapsack constraint;
	 * to start search, set `row=-1`;
	 * \param[in,out] family if `1`, return value is index of knapsack constraint
	 * that is a parity constraint as well;
	 *  otherwise, return value is `0`.
	 * \return index of knapsack constraint that follows constraint with index `row`,
	 * negative return value indicates that `row` is the last knapsack constraint.
	 */
	int getNextKnapsack(int row, int &family);

	/**
	 * For each knapsack inequality, the procedure calls `LCI()` to generate a lifted cover inequality;
	 * in case of failure (`LCI()` returns `false`), `GUB_LCI()` is called to generate a GUB-cover inequality.
	 * \return number of cuts generated.
	 */
    int knapsackCuts();

    /**
     * Given an inequality
     *     `sum(A[j]*x[j] for j=0,...,n-1) <= b`
     * and a point `X` with all components from `[0,1]`,
     * find a cover `C` with positive excess value `lambda`,
     *      `lambda = b- sum(A[j]*X[j] for j in C) > 0`.
     * To find a cover, this procedure implements the LP-heuristic.
     *
     * \param[in] n number of binary variables;
     * \param[in] A vector of `n` coefficients;
     * \param[in] b right-hand side value;
     * \param[in] X vector of size `n` with all components from `[0,1]`;
     * \param[out] ipCov array of size `k` (`k` denotes return value),
     * which represents cover `C`;
     * \param[out] lambda excess value of output cover.
     * \return size of cover if non-zero.
     */
	int lpCover(int n, double* A, double b, double* X, int* ipCov, double& lambda);

    /**
     * The procedure generate mixed knapsack cuts for sets
     * \f$X =\{(y,s)\in \{0,1\}^N\times R^2_+: \sum_{j\in N} a_j y_j - s_1 + s_2\le b\}\f$,
     * with \f$a_j > 0\text{ for }j\in N\f$.
     *
     * Let _C_ be a cover, \f$\sum_{j\in C} a_j = b + \lambda\f$,
     * \f$\bar{a} > \lambda > 0\f$, where  \f$\bar{a} = \max_{j\in C} a_j\f$.
     * Let \f$E(C) = C \cup \{j\in N\setminus C: a_j > \bar{a}\}\f$ be the extension of _C_.
     * Then the _mixed cover inequality_
     * \f$\sum_{j\in E(C)} \min(a_j,\lambda) y_j - s_1 \le -lambda + \sum{j\in C} \min(a_j,\lambda)\f$
     * is valid for \f$\mathrm{conv}(X)\f$.
     *
     * \param[in] szA number of variables in template constraint;
     * \param[in] szP number of rows used to derive template constraint;
     * \param[in,out] rhs right-hand side of template constraint, or returned cut;
     * \param[out] type type of cut inequality (bitwise OR of `enVarType` members);
     * \param[in] bLocal if `true`, template constraint is local; otherwise, it is global;
     * \return `true` size of cut, or '0' if no cut has been found.
     */
    int mixedKnapsackCut(int szA, int szP, double &rhs, unsigned &type, bool bLocal);

////////////////////////////////
/// Parity cuts
    /**
     * The procedure generate parity cuts by calling `parityCut()` for every parity constraint.
 	 * A _parity constraint_ is any constraint from which we can deduce
 	 * a _parity condition_ that the sum \f$\sum_{j \in S} x_j\f$ of binary variables is odd (even).
 	 * This parity condition is expressed by the following family of inequalities:
     *  \f{align*}{
     *     \sum_{j\in J} x_j - \sum_{j\in S\setminus J} x_j \le |J|-1,\quad \forall\; J\subseteq S,\; |J|\text{ is even (odd)}.
     * \f}
     *
     * \return number of cuts generated.
     */
    int parityCuts();

    /**
     * \param[in] gap initial gap;
     * \param[in] sz size of arrays `ipCol` and `dpVal`;
     * \param[in,out] ipCol,dpVal represent left part of either input inequality or output cut:
     *  `sum(i=0,...,sz-1) dpVal[i]*x(ipCol[i]) <= b`;
     * \param[in,out] b right-hand side of either input inequality or output cut.
     * \param[in,out] parity on input, parity of constraint (1 - even, 0 - odd);
     *  on output, type of cut generated.
     *  \return `true` if cut has been generated.
     */
    bool parityIneq(unsigned &parity, double gap, int sz, int *ipCol, double *dpVal, int &b);

    /**
     * The procedure tries to produce a parity cut for a given parity inequality.
     *
     * \param[in] tol cut is accepted only if it is violated by more than `tol`;
     * \param[in] sz size of arrays `ipCol` and `dpVal`;
     * \param[in,out] ipCol,dpVal represent left part of either input inequality or output cut:
     *  `sum(i=0,...,sz-1) dpVal[i]*x(ipCol[i]) <= rhs`;
     * \param[in,out] rhs right-hand side of either input inequality or output cut;
     * \param[out] cutType type of returned inequality (bitwise OR of `enVarType` members).
     *  \return number of cuts generated.
     */
    int parityCut(double tol, int sz, int *ipCol, double *dpVal, double &rhs, unsigned &cutType);

////////////////////////////////////////////////////////////////////
//  P A C K I N G  cuts
//////////////////////////
    /**
     * \param[in] sz size (number of variables) of click;
     * \param[in,out] ipClick on input array of size `sz` with list of variables in initial click constraint,
     * on output `ipClick` contains list of variables in extended click of size `k`, where `k` is return value;
     * \param[in] dWeight sum of values of all variables in click;
     * \param[in,out] bGlobal if `true` on input, click has been derived from global constraints
     * \param[in] candNum number of variables written in `m_ipArray` that are candidates for adding to click.
     * \return number of variables in extended click, or `0` if extended click is not violated by current LP solution.
     */
	int extendClick(int sz, int *ipClick, double dWeight, bool &bGlobal, int candNum);

	/**
	 * \param[in] row index of packing constraint;
	 * \param[in] candNum number of previously generated click cuts;
	 * \param[in] cutNum number of click cuts generated earlier;
	 *  any time a new click cut is generated, it is added to matrix only if
	 *  its hash value does not match any hash value of click cut generated earlier.
	 * \return `true` if at least one cut was generated.
	 *
	 */
    bool clickCut(int row, int candNum, int cutNum);

    /**
     * The procedure call `clickCut()` for all tight packing constraints.
     * \returns number of generated cuts.
     */
    int clickCut();
    
///////////////////////////////////////////////////////////////
//  MIR inequalities: Chvatal-Gomory, Residual Capacity cuts
///////////////////////
    /**
     * Given a valid inequality that is a linear combination of problem constraints,
     * the procedure looks for another problem inequality to eliminate from the ininial inequality
     * one variables which current value is far from its bounds.
     * \param[in] round number of inequalities used to obtain initial inequality;
     * \param[in] m0 only first `m0` matrix rows are considered as candidates for mixing;
     * \param[in] sz size of arrays `dpA` and `ipCol`;
     * \param[in] dpA,ipCol `sum(i=0,...,sz-1) dpA[i]*x(ipCol[i])` is right part of initial inequality;
     * \param[in] ipRowPos only if `ipRowPos[i]==0`, inequality `i` can be used for mixing;
     * \param[out] row index of inequality to be mixed with initial one;
     * \param[out] col index of variable to be eliminated;
     * \param[out] pivot inequality indexed by `row` must first be multiplied by `pivot` and then added to initial inequality;
     * \param[in] allBinVars mixing inequality must be binary.
     * \return `true` is appropriate inequality has been found; otherwise, `false`.
     */
    bool getCtrForMixing(int round, int m0, int sz, double* dpA, int* ipCol,
    				int* ipRowPos, int& row, int& col, double& pivot, bool allBinVars);

    /**
     * Given a valid inequality that is a linear combination of problem constraints,
     * the procedure mixes this inequality with that found by `getCtrForMixing()`
     * to eliminate one variable, which current value is far from its bounds, from the initial inequality
     * given as
     *     `sum(i in 0..szA-1) dpA[i] x(ipCol[i]) <= b`.
     *
     * \param[in] row constraint (row) to  be mixed with input inequality;
     * \param[in] pivot when mixing, constraint indexed by `row` is multiplied by `pivot`;
     * \param[in,out] szA size of input and output constraints;
     * \param[in,out] dpA coefficients of input and output constraints, array of size `szA`;
     * \param[in,out] ipCol column (variable) indices of input and output constraints, array of size `szA`;
     * \param[in,out] ipColPos if non-negative, `ipColPos[col]` is position of column `col` in array `ipCol`, i.e.,
     *   `ipCol[ipColPos[col]]=col`; negative, value of `ipColPos[col]` indicates that the value `col`
     *   is not present in `ipCol`;
     * \param[in,out] b right-hand side of input and output inequalities;
     * \param[in,out] szP number of constraints that have been mixed;
     * \param[in,out] dpPivot array of size `szP`, when mixing, constraint `ipRow[i]` is multiplied by `dpPivot[i]`;
     * \param[in,out] ipRow indices of mixing constraints, array of size `szP`;
     * \param[in,out] ipRowPos if non-negative, `ipRowPos[row]` is position of constraint `row` in array `ipRow`, i.e.,
     *   `ipRow[ipRowPos[row]]=row`; negative, value of `ipRowPos[row]` indicates that constraint indexed by `row`
     *   is not mixing constraint.
     * \param[in,out] bLocal if `true` on input (output), then input (output) constraint is local.
     */
    void mixCtrs(int row, double pivot,
        int &szA, double* dpA, int* ipCol, int* ipColPos, double& b,
        int &szP, double* dpPivot, int *ipRow, int* ipRowPos, bool& bLocal);

    /**
     * The procedure does some preparations for subsequent calls of 'mirCut()` and `mixedKnapsackCut()`.
     * In particular, the procedure complements some variables,
     * and then compute an array of factors; the inequality is multiplied, in turn,
     * by each of these factors and the result is "mixes-integer rounded".
     *
     * \param[in] szA size of arrays `dpA` and `ipCol`;
     * \param[in] dpA,ipCol `sum(i in 0..szA-1) dpA[i]*x[ipCol[i]] <= f0` is input inequality;
     * \param[in,out] f0 on output, `f0` is right hand side of the inequality with some variables complemented;
     * \param[in,out] bLocal if `true`, constraint is local.
     * \param[out] dpQ array of factors.
     * \return number of factors stored in `dpQ`.
     */
    int mirPrepareForCuts(int szA, double* dpA, int *ipCol, double &f0, bool &bLocal, double *dpQ);
    
    /**
     * Given an inequality an a scalar, divide the inequality by the scalar
     * and use the resulting inequality to build a MIR-cut that cuts off
     * the optimal solution of the current node LP.
     *
     * \param[in] szA size (number of non-zero coefficients in) of input inequality;
     * \param[in] dpA array of size `szA` of inequality coefficients;
     * \param[in] ipCol array of size `szA` of inequality variables;
     * \param[in,out] szP input inequality is a result of mixing `szP`other inequalities;
     * \param[in] delta input inequality is divided by this scalar;
     * \param[out] beta right-hand side value of output inequality (cut);
     * \param[out] violat amount by which output inequality is violated by optimal solution to current node LP
     * \param[out] type type of output inequality;
     * \param[in] flag if `true`, output cut (if such one exists) is produced; otherwise,
     * non-zero return value inly indicates that a required cut exists;
     * \return size of the output inequality (cut), where this output inequality is stored in
     * the arrays `*dpA0{m_dpUd}` (coefficients), `*ipCol0{reinterpret_cast<int*>(m_dpFd)}`
     *  (indices of involved variables).
     */
    int buildDeltaMirCut(int szA, double* dpA, int* ipCol, int szP, double delta, double& beta,
          double& violat, unsigned &type, bool flag);

    /**
     * Given an inequality, eliminate those non-binary variables that are present
     * in tight variable-bound constraints.
     *
     * \param[in,out] szA size of (number of non-zero coefficients in) input and output inequalities;
     * \param[in,out] dpA list (of size of `szA`) of non-zero coefficients of input and output inequalities;
     * \param[in,out] ipCol list of variables (columns) having non-zero coefficients in input and output inequalities;
     * \param[in,out] b right-hand side value of input and output inequalities;
     * \param[in,out] szP input (resp., output) inequality is a result of mixing `szP`other inequalities;
     * \param[in,out] bLocal if `true`, input (resp., output) inequality is local;
     * \param[in] m0 only first `m0` constraints currently present in matrix can be used
     * for mixing with other constraints;
     * \param[out] bMxKnapsack if `true`, input inequality has not been changed, and,
     * therefore, it can be a starting point for producing a mixed-knapsack cut.
     */
    void substituteVarBounds(int& szA, double* dpA, int *ipCol, double& b,
          int& szP, bool& bLocal, int m0, bool &bMxKnapsack);

    /**
     * Given an inequality, find a collection of disjoint sets of its binary variables,
     * so that any two variables from the same set cannot simultaneously take the value of `1`.
     *
     * \param[in] szA size of input inequality.
     * \param[in] dpA array of size `szA` of inequality coefficients;
     * \param[in] ipCol array of size `szA` of inequality variables;
     * \param[in] ipColPos array of size `m_iN`, if `ipColPos[j] >= 0`, then
     * `ipCol[ipColPos[j]] == j`.
     * \param[out] bLocal `false` if output partition has been derived from local GUB-constraints.
     * \return size (number of sets) of the output partition, which is stored in
     * `*ipGubCol{reinterpret_cast<int*>(m_dpFb)}`.
     */
	int mirGetGUBs(int szA, double* dpA, int *ipCol, int *ipColPos, bool &bLocal);
          
	/**
	 * The procedure builds a MIR-cut starting from the input inequality.
	 *
	 * \param[in] szA,dpA,ipCol valid inequality `sum(i in 0..szA-1) dpA[i]*x(ipCol[i]) <= rhs`;
	 * \param[in] szP number of inequalities (other than variable bounds) mixed to get input inequality.
	 * \param[in,out] rhs right hand side of input and output inequality;
	 * \param[out] cutType cut type (bitwise OR of `enVarType` members);
	 * \param[in] bLocal if `true`, cut is local.
	 * \return size of the cut, or '0` if no cut has been found.
	 */
    int mirCut(int szA, double* dpA, int* ipCol, int szP, double &rhs, unsigned &cutType, bool bLocal);

    /**
     * The procedure tries to generate a _MIR_ cut starting from the inequality indexed by `row`.
     *
     * \param[in] row constraint index;
     * \param[in] side side of constraint `row`;
     * \param[in] m0 first `m0` rows of constraint matrix may be used as mixing inequalities;
     * \param[in] needCut if `needCut & 0x1`, _MIR_ cuts to be generated; if `needCut & 0x2`, _mixed knapsack_ cuts to be generated.
     * \param[in,out] sz cut size;
     * \param[out] cutType cut type (bitwise OR of `enVarType` members);
     * \param[in,out] rhs right hand side.
     * \return `1` if MIR cut was generated, `-1` if mixed knapsack cut was generated, and `0` otherwise.
     */
    int mixedIntCut(int row, bool side, int m0, int needCut, int &sz, unsigned &cutType, double &rhs);

    /**
     * The procedure calls MIR and mixed knapsack cut generators.
     * \param[out] mirCutGenerated number of MIR cuts generated;
     * \param[out] mxKnapsackCutGenerated number of mixed-knapsack cuts generated.
     */
    void mixedIntCuts(int &mirCutGenerated, int &mxKnapsackCutGenerated);
        
//////////////////////////////////////////////////////////////////////
////// GOMORY Cuts                      //////////////////////////////
    /**
     * The procedure generates a _mixed-integer Gomory_ cut
     *  for a weighted sum, `sum(i=0..sz-1) dpVal[i]*x[ipCol[i]]`, of integer variables taking fractional values;
     *  all `dpVal[i]` are integers.
     *  Usually, `sz=1`, `dpVal[0]=1.0`, and a cut is generated for only one variable `x[ipCol[0]]`.
     *
     * \param[in] maxCutSize cuts of size greater than `maxCutsSize` are rejected;
     * \param[in] sz number of variables;
     * \param[in] ipCol list of `sz` variable indices;
     * \param[in] dpVal list of `sz` coefficients;
     * \param[out] rhs right hand side of returned inequality;
     * \param[out] type type of returned inequality (bitwise OR of `enVarType` members).
     * \return size of generated cut, or '0' if not cut has been found.
     * \remark If we denote te returned value by 'sz', then the cut is given as
     *      'sum(i in 0..sz-1) m_dpArray[i]*x(m_ipArray[i]) <= rhs'.
     */
    int GomoryCut(int maxCutSize, int sz, int *ipCol, double *dpVal, double &rhs, unsigned &type);

    /**
     * The procedure generates a family of _mixed-integer Gomory_ cuts.
     *
     * \param[in]  bDense if `true` dense cuts are generated;
     * \param[out] denseCutGenerated number of dense cuts generated;
     * \param[out] sparseCutGenerated number of sparse cuts generated
     */
    void GomoryCuts(bool bDense, int &denseCutGenerated, int &sparseCutGenerated);

//////////////////////////////////////////////////
//
    /**
     * The procedure transforms a binary matrix into the upper triangular form.
     *
     * \param[in] m number of rows;
     * \param[in] n number of columns;
     * \param[in] n1 size of row (in integer words);
     * \param[out] S permutation of rows (array of size `m`);
     * \param[in,out] A binary matrix written row-wise (`n1` words per row),
     *  on output it is transformed into upper triangular form.
     * \return rank of matrix `A`.
     */
	int mod2Basis(int m, int n, int n1, int* S, int *A);

	/**
	 * The procedure decides which rows and columns be included into the matrix
	 * that is used for generating Chvatal-Gomory cuts with coefficients `0,1/p,...,(p-1)/p`.
	 *
	 * \param[in] p integer from `{2,3,5}`;
	 * \param[in] neqMsk bitwise OR of the members of `CMIP::enCtrType` members,
	 *       rows representing inequalities that match this mask are added to matrix;
	 * \param[in] eqMsk bitwise OR of the members of `CMIP::enCtrType` members,
	 *       rows representing equalities that match this mask are added to matrix;
	 * \param[in] maxVarSlack maximum slack value for variable to be added to matrix;
	 * \param[in] maxCtrSlack maximum slack value for integer constraints to be added to matrix;
	 * \param[out] rowNum number of rows in matrix;
	 * \param[out] colNum number of columns in matrix.
	 */
	void buildModPMatrix(int p, int neqMsk, int eqMsk,
			double maxVarSlack, double maxCtrSlack, int &rowNum, int &colNum);

	/**
     * The procedure generates a family of {0,1/2}-Chvatal-Gomory cuts.
     *
     * \param[in] maxVarSlack maximum slack value for variable to be included in binary matrix;
     * \param[in] maxCtrSlack maximum slack value for integer constraints to be included in binary matrix;
     * \param[out] denseCutNum number of dense cuts generated;
     * \param[out] sparseCutNum number of sparse cuts generated.
     */
	void mod2Cuts(double maxVarSlack, double maxCtrSlack, int &denseCutNum, int &sparseCutNum);

	////////////////////////////////////
	// Disjunctions cuts

	/**
	 * Given an inequality `sum(i=0,...,sz-1) dpA[ipCol[i]] dpX[ipCol[i]] >= 1`, where `0 <= dpX[j] <= 1` for all `j`.
	 * Let P denote the solution set of the given inequality, and let `x^*=dpX`.
	 * `computeAlpha_P()` finds the value of `alpha_P` defined to be `min{u^Tx^*: u^Tx \ge 1 for all x\in P, u \ge 0}`.
	 * If `alpha_P < 1`, then, for any `u` such that `u^Tx^*=alpha_P`,
	 * the inequality `u^T x \ge 1` is valid for `P` and cuts off `x^*` from `P`.
	 *
	 * \param[in] sz number of variables with non-zero coefficients;
	 * \param[in] ipCol array of `sz` variable indices;
	 * \param[in] dpA array of `sz` coefficients;
	 * \param[in] dpX point `x^*`;
	 * \param[out] alpha `alpha_P`;
	 * \param[out] q see description of return value;
	 * \return size, `k`, of support of `u`, non-zero coefficients of `u` are `dpA[ipCol[i]]/q`, `i=0,...,k-1`,
	 * where the entries of `ipCol` have been permuted so that `dpX[ipCol[0]] <= dpX[ipCol[1]] <= ... <= dpX[ipCol[sz-1]]`.
	 */
	int computeAlpha_P(int sz, int *ipCol, double * dpA, double *dpX, double &alpha, double &q);

	/**
	 * The procedure computes the most violated inequality (if any)
	 * that separates a given point from the union of two sets
	 *     `S_1 ={x: sum(i=0,...,sz1-1) dpA1[ipCol1[i]]*x(ipCol1[i]) >= 1, 0 <= x(ipCol1[i]) <= 1, i=0,...,sz1-1}` and
	 *     `S_2 ={x: sum(i=0,...,sz2-1) dpA2[ipCol2[i]]*x(ipCol2[i]) >= 1, 0 <= x(ipCol2[i]) <= 1, i=0,...,sz2-1}`.
	 *
	 *  \param[in,out] sz1 size of set `S_1`; on output, `sz1` is number of variables from `S_1` in separating inequality;
	 *  \param[in,out] ipCol1 indices of variables from set `S_1`;
	 *  \param[in] dpA1 array of coefficients of first inequality;
	 *  \param[in] dpX1 stores values `x(ipCol1[i]),...,x(ipCol1[sz1-1])`;
	 *  \param[in,out] sz2 size of set `X_2`; on output, `sz2` is number of variables from `X_2` in separating inequality;
	 *  \param[in,out] ipCol2 indices of variables from set `S_2`;
	 *  \param[in] dpA2 array of coefficients of second inequality;
	 *  \param[in] dpX2 stores values `x(ipCol2[i]),...,x(ipCol2[sz2-1])`;
	 *  \return if non-zero, size `k` of separating inequality that is written as
	 *      `sum(i=0,...,k-1) m_dpArray[i]*x(m_ipArray[i]) >= 1`.
	 */
	int disjunctionCut(int &sz1, int *ipCol1, double * dpA1, double *dpX1,
				int &sz2, int *ipCol2, double * dpA2, double *dpX2);

	/**
	 * For a given non-binary variable, `getTightLoVarBound()` verifies
	 * whether there exists a lower variable bound that is tight for
	 * an optimal solution of the node LP being processed.
	 *
	 * \param[in] j index of non-binary variable;
	 * \param[out] val,lhs if return value `c` is not `-1`, then `x[j] + val*x[c] >= lhs` is required variable bound.
	 * \return index of binary variable, or `-1`.
	 */
	int getTightLoVarBound(int j, double &val, double &lhs);

	/**
	 * For a given non-binary variable, `getTightUpVarBound()` verifies
	 * whether there exists an upper variable bound that is tight for
	 * an optimal solution of the node LP being processed.
	 *
	 * \param[in] j index of non-binary variable;
	 * \param[out] val,rhs if return value `c` is not `-1`, then `x[j] + val*x[c] <= rhs` is required variable bound.
	 * \return index of binary variable, or `-1`.
	 */
	int getTightUpVarBound(int j, double &val, double &rhs);

	/**
	 * The procedure calls `disjunctionCut()` to build a disjunction cut base on a disjunction derived from a given inequality.
	 *
	 * Given an inequality
	 *       `sum(j\in B) a[j]*x[j] + sum(j in N) a[j]*x[j] >= b`,
	 * where
	 *     - `B` is a subset of binary variables,
	 *     - `b` and all `a[j]` are positive,
	 *     - all variables take nonnegative values,
	 * then the following disjunction is valid:
	 *      `sum(j in J) x_j \ge 1  or sum(j in B\J) a[j]*x[j] + sum(j in N) a[j]*x[j] >= b`.
	 *
	 * To cut off a given point `x^*`, we move to `J` some binaries with big coefficients `a[j]` and small values `x^*[j]`.
	 *
	 * \param[in] row constraint index;
	 * \param[in] side if `true`, left hand side inequality is used to derive disjunction;
	 * otherwise, right hand side inequality is used;
	 * \param[out] type type of returned cut (bitwise OR of `enVarType` members);
	 * \param[out] lhs left hand side of returned constraint.
	 * \return cut size, or `0` if no cut has been found.
	 */
	int oneRowDisjunction(int row, bool side, unsigned& type, double& lhs);

	/**
	 * The procedure calls `oneRowDisjunction()`
	 *  for each tight constraint of type `CTR_KNAPSACK|CTR_MX_KNAPSACK|CTR_MX_01`.
	 *  \return number of generated cuts.
	 */
	int oneRowDisjunctions();

////////////////////////////////////////////////////////////////////
//                    Functions for managing cuts                //
////////////////////////////////////////////////////////////////////

	/**
	 * The procedure is called to invoke all internal cut generating procedures.
	 *
	 * \return number of cuts generated.
	 */
    int autoCuts();

    /**
     * \param[in] ipFracNum  array of size `AUTO_CUT_CYCLE`,
     *       where `ipFracNum[m_iAutoCutRound-1) % AUTO_CUT_CYCLE]`
     *        is number of integer variables taking fractional values after previous round of generating cuts;
     * \param[in] dpObjProgress array of size `AUTO_CUT_CYCLE`,
     *       where `dpObjProgress[m_iAutoCutRound-1) % AUTO_CUT_CYCLE]`
     *        is objective value after previous round of generating cuts;
     * \param[out] n number of variables (columns);
     * \param[out] firstAutoCut index of first new cut generated by __MIPCL__.
     * \param[out] ipCutNum three values: `ipCutNum[0]` number of strong cuts generated by user defined `separate()` and `genCut1()`,
     *       `ipCutNum[1]` number of auto cuts generated by __MIPCL__ itself,
     *       `ipCutNum[2]` number of week cuts generated by user defined `genCut2()`.
     */
    void generateCuts(int *ipFracNum, double *dpObjProgress, int& n, int &firstAutoCut, int* ipCutNum);
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
//   SERIALIZATION: functions are used when restarting computations
    /**
     * The procedure reads/writes from/to a given stream the matrix of the problem being solved.
     * \param[in,out] ar input or output stream;
     * \param[in] is_storing if `true`, matrix is storing; otherwise, it is restoring.
     */
    void serializeMatrix(std::fstream& ar, bool is_storing) final;

    /**
     * The procedure reads/writes from/to a given stream the best (record) solution.
     * \param[in,out] ar input or output stream;
     * \param[in] is_storing if `true`, solutionis storing; otherwise, it is restoring.
     */
    void serializeRecSolution(std::fstream& ar, bool is_storing);

    /**
     * The procedure reads/writes from/to a given stream the values of tolerance variables.
     * \param[in,out] ar input or output stream;
     * \param[in] is_storing if `true`, tolerance variables are storing; otherwise, they are restoring.
     */
    void serializeTolVars(std::fstream& ar, bool is_storing) final;

    /**
     * The procedure reads/writes from/to a given stream different parameters that affect the solution process.
     * \param[in,out] ar input or output stream;
     * \param[in] is_storing if `true`, parameters are storing; otherwise, they are restoring.
     */
    void serializeFlags(std::fstream& ar, bool is_storing) final;

protected:
    /**
     * The function stores into (or restores from) the stream `ar`
     * CMIP objects (all its member storing permanent data).
     * \param[in] ar reference to a stream;
     * \param[in] is_storing if `true`, the object is written to the stream;
     *  otherwise, the object is restored from the stream.
     * \remark Derived classes may overload this function to store additional information.
     * In such a case, call first `serialize()` of the base class.
     *  \sa CLP::serialize().
     */
    virtual void serialize(std::fstream& ar, bool is_storing) override;

private:
    void writeMIP(); ///< The procedure writes the current problem state into the file `m_strProblemName.res`


    /**
     * The procedure is used in case of restart
     * to continue solving the problem  from the state saved in the file `m_strProblemName.res`.
     */
    void readMIP();
/////////////////////////////////////////////////////
// Statistics
/////////////////////////////////////////////////////
    char* getProbStatStr(char *str=0) final;

    /**
     * The functions calls `mipInfo()` to print
     * the current state of the solution process.
     * \param[in] nodeNum total number of nodes;
     * \param[in] leafNum number of leaves (non-processed nodes) in the search tree;
     * \param[in] upperBound upper bound on the optimal objective value;
     * \param[in] header if `true`, the header is displayed.
     */
    void __mipInfo(int nodeNum, int leafNum, double upperBound, bool header);
protected:
    /**
     * The function prints into the standard output stream a message string
     * which describes the current state of the solution process.
     *
     * \param[in] timeElapsed string representation of time elapsed since solution process started;
     * \param[in] nodeNum total number of nodes;
     * \param[in] leafNum number of leaves (non-processed nodes) in the search tree;
     * \param[in] bestObjVal best objective value achieved so far;
     * \param[in] objBound lower or upper (`sense=true`) bound on the optimal objective value;
     * \param[in] gap relative (in percents) difference between `bestObjVal` and `objBound`;
     * \param[in] solsFound number of solutions found so far;
     * \param[in] sense if `true`, objective is maximized; otherwise, minimized;
     * \param[in] header if `true`, the header is displayed.
     * \remark When developing an application with a GUI interface, the user may wish to overload this function.
     */
    virtual void mipInfo(char *timeElapsed, int nodeNum, int leafNum,
    		double bestObjVal, double objBound, double gap, int solsFound, bool sense, bool header);

    /**
     * After each _cut generating round_ (i.e., after calling the cut generating procedure),
     * the solver calls `cutInfo()` to inform about the progress in decreasing the objective value of the relaxation LP.
     *
     * \param[in] time time elapsed since the solution process has started;
     * \param[in] round cut generating round;
     * \param[in] objVal objective value of the relaxation LP;
     * \param[in] fracNum number of integral variables taking fractional values in the current LP solution;
     * \param[in] cutNum number of cuts generated at this particular round.
     * \remark When developing an application with a GUI interface, the user may wish to overload this function.
     */
    virtual void cutInfo(__LONG time, int round, double objVal, int fracNum, int cutNum);

    /**
     * After having been processed the root node,
     * the solver prints into the standard output stream
     * how many cuts of each type were generated, and how many of them were present in the matrix
     * at the beginning of the next cut generating round.
     *
     * \remark When developing an application with a GUI interface, the user may wish to overload this function.
     */
	virtual void cutStatistics();

	/**
	 * The function prints into the stream `out` solution statistics such as:
	 * whether the problem is feasible or not, solution time,
	 * whether optimality was proven or not, and so on.
	 * \param[in] MIPCLver __MIPCL__ version;
	 * \param[in] solTime string representation of solution time;
	 * \param[in] timeLimit if `true`, time limit has been reached;
	 * \param[in] nodeNum number of processed nodes;
	 * \param[in] feasible if `false`, problem is infeasible;
	 * \param[in] hasSolution if `true`, solution (not necessary optimal) has been found;
	 * \param[in] objVal objective value of best solution;
	 * \param[in] opt if `true`, optimality has been proven;
	 * \param[in] gap duality gap;
	 * \param[in] gapLimit if `true`, solution process has been terminated since required duality gap had been reached;
	 * \param[in] bound depending on value of `sense`, lower or upper bound on optimal objective value;
	 * \param[in] difficultNodes number of nodes processed with numerical difficulties;
	 * \param[in] out reference to a stream.
	 * \sa \ref MIPCLmsgs
	 */
	virtual void solStatistics(std::ostream &out, const char* MIPCLver,
			const char* solTime, bool timeLimit, int nodeNum,
			bool feasible, bool hasSolution, double objVal,
			bool opt, double gap, bool gapLimit, double bound,
			int difficultNodes);

private:
	/**
	 * The function calls `solStatistics()` to  prints into the stream `out` solution statistics such as:
	 * whether the problem is feasible or not, solution time,
	 * whether optimality was proven or not, and so on.
	 * \param[in] out reference to an output stream.
	 */
	void _solStatistics(std::ostream &out);

public:
	/**
	 * Given a MIP, the difference between the optimal values of the relaxation LP and MIP itself
	 * is known as _duality_ (or _integrality_) _gap_.
	 * Generating cuts and doing branching, the solver decreases this duality gap.
	 * At a particular moment, the duality gap is the difference between the current upper and lower bounds
	 * on the optimal objective value. If this current duality gap does not exceeds some limit
	 * (given by the user or computed by the solver), the solver terminates,
	 * and the best found solution is printed as a solution to the problem.
	 * \param[in] gap duality gap.
	 * \sa optimize(), BranchAndCut().
	 */
	void setDualGap(double gap)
		{m_dDualGap=gap;}

	int getBranchAndCutNodeNumber()
		{return m_iBranchAndCutNodes;} ///< \return total number of branch-and-cut nodes.

	int getNoOfActiveNodes() const; ///< \return number of leaves in the search tree.

	double getUpperBound() const; ///< \return current upper bound on the optimal objective value.

    /**
     * \return `true` if CMIP cut-messages are to be printed.
     * \sa CMIP::switchMipInfoMsg()
     */
	bool mipCutInfoMsg()
	{return (!(m_iInfoMsgFlag & 0x12))? true: false;}

    /**
     * \return `true` if CMIP run time messages are to be printed.
     * \sa CMIP::switchMipInfoMsg()
     */
    bool mipTreeInfoMsg()
    {return (!(m_iInfoMsgFlag & 0x14))? true: false;}

	/**
	 * The function is called to switch on or off printing CMIP run time messages.
	 * \param[in] cutInfo if `true`, the solver will start printing cut statistics;
	 *  otherwise,  will stop printing such messages;
	 * \param[in] treeInfo if `true`, the solver will start printing branch-and-cut statistics;
	 *  otherwise,  will stop printing such messages.
	 */
	 void switchMipInfoMsg(bool cutInfo, bool treeInfo);

#ifdef __PYTHON_MODULE_
/**
 * When __MIPCL__ is called from __Python__, all handles of variables are their indices,
 * i.e., the handle of variable `x[i]` is `i`. To make communication between __Python__ and __MIPCL__
 * easier, `setSolution()` lists solution components in order of their handles (initial indices).
 */
void setSolution();

/**
 * This function must be called only after calling `sortSolutionByHandles()`.
 * \param[in] ind index (= handle) of variable.
 * \return value of the variable indexed by `ind`.
 */
double getOptVarValue(int ind) const;
#endif

};
#endif // #ifndef __CMIP__H
