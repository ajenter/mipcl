///////////////////////////////////////////////////////////////
/**
 * \file except.h Interface for CException class
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
#ifndef __EXCEPT__H
#define __EXCEPT__H

#ifdef _WINDOWS
#ifndef MIP_API
#ifdef MIP_EXPORTS
#define MIP_API __declspec(dllexport)
#else
#define MIP_API __declspec(dllimport)
#endif
#endif
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
#endif
#endif

/// Interface for exception classes used in __MIPCL__ applications.
class MIP_API CException
{
#define MAX_MSG_LEN 256 ///< `MAX_MSG_LEN-1` is maximal length of any message.
protected:
  char m_szErrorMessage[MAX_MSG_LEN]; ///< Buffer for storing messages.

public:

	CException() {}; ///< The default constructor.
	/**
	 * The constructor that is usually used.
	 *
	 * \param[in] msg string describing the cause of throwing this exception.
	 */
	CException(const char* msg);
	virtual ~CException() {}; ///< The destructor.

	/**
	 * \return pointer to the explanatory string.
	 */
	const char* what()
		{return m_szErrorMessage;}
};

/// Exceptions thrown if any problem arises when allocating (reallocating) memory.
class MIP_API CMemoryException: public CException
{
	static const char memError[];  ///< The string identifying a memory allocation error.
public:
	/**
	 * \param[in] msg message string.
	 */
	CMemoryException(const char* msg); ///< The constructor.
	~CMemoryException() {}; ///< The destructor.
};


/// Exceptions thrown if any problem arises when manipulating with files.
class MIP_API CFileException: public CException
{
	static const char fileOpenError[]; ///< The string identifying a file open error.
public:
	/**
	 * \param[in] where name of procedure that has thrown this exception;
	 * \param[in] fileName name of file that caused any problem.
	 */
	CFileException(const char* where, const char* fileName); ///< The constructor.

	~CFileException() {}; ///< The destructor.
};

/// Exceptions thrown when any unexpected computation instability is detected.
class MIP_API CDegenException: public CException
{
public:
	/**
	 * \param[in] msg message string.
	 */
	CDegenException(const char* msg); ///< The constructor.

	~CDegenException() {}; ///< The destructor.
};

/// Exceptions thrown when any crucial error is detected in input data.
class MIP_API CDataException: public CException
{
	static const char dataError[]; ///< The string identifying a data error.
public:
	/**
	 * \param[in] msg message string.
	 */
	CDataException(const char* msg);
	~CDataException() {};
};

#endif // #ifndef __EXCEPT__H
