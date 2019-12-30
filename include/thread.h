///////////////////////////////////////////////////////////////
/**
 * \file thread.h Declares macros to simplify using different multithreding libraries.
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

#ifndef _THREAD_H_
#define _THREAD_H_

#ifndef __ONE_THREAD_
#ifdef _CPP_THREADS
#include <thread>
#include <mutex>
#include <shared_mutex>

typedef std::thread _THREAD;  ///< Alias for `pthread_t`.
typedef std::mutex _MUTEX; ///< Alias for `pthread_mutex_t`.
typedef std::shared_mutex _RWLOCK; ///< Alias for `pthread_rwlock_t`.

#define _THREAD_CREATE(_thread,start_thread,param) \
		_thread = std::thread(start_thread,param); ///< Creates a new thread with start function `start_thread()`, that gets `param` as an argument.
#define _THREAD_JOIN(thread) \
		thread.join(); ///< Attaches `thread` to the calling thread.
#define _THREAD_CLOSE(thread) ///< Empty placeholder to be consistent with windows threads.

#define _RWLOCK_INIT(rwLock) ///< Initializes the read-write lock pointed by `p_rwLock`.
#define _RWLOCK_DESTROY(rwLock) ///< Destroys the read-write-lock pointed by `p_rwLock`.
#define _RWLOCK_RDLOCK(rwLock) \
		(rwLock)->lock_shared(); ///< Applies a read lock to the read-write lock pointed by `rwLock`.
#define _RWLOCK_SAFE_RDLOCK(rwLock) \
	if (rwLock) \
		rwLock->lock_shared(); ///< A safer version of `_RWLOCK_RDLOCK()`.
#define _RWLOCK_WRLOCK(rwLock) \
		(rwLock)->lock(); ///< Applies a write lock to the read-write lock pointed by `rwLock`.
#define _RWLOCK_SAFE_WRLOCK(rwLock) \
	if (rwLock) \
		rwLock->lock();  ///< A safer version of `_RWLOCK_WRLOCK()`.
#define _RWLOCK_UNLOCK_RDLOCK(rwLock) \
		(rwLock)->unlock_shared(); ///< Unlocks the read-write lock pointed by `rwLock`.
#define _RWLOCK_UNLOCK_WRLOCK(rwLock) \
		(rwLock)->unlock(); ///< Unlocks the read-write lock pointed by `rwLock`.
#define _RWLOCK_SAFE_UNLOCK_RDLOCK(rwLock) \
	if (rwLock) \
		rwLock->unlock_shared();  ///< A safer version of `_RWLOCK_UNLOCK()`.
#define _RWLOCK_SAFE_UNLOCK_WRLOCK(rwLock) \
	if (rwLock) \
		rwLock->unlock();  ///< A safer version of `_RWLOCK_UNLOCK()`.
#define _RWLOCK_ASSIGN_P(plock,pvalue) \
	plock=pvalue; ///< Assigns the value of `pvalue` to `_RWLOCK` pointer `plock`.

#define _MUTEX_INIT(mutex) ///< Initializes the mutex represented by `mutex`.
#define _MUTEX_DESTROY(mutex) ///< Destroys the mutex represented by `mutex`.
#define _MUTEX_LOCK(pmutex) \
	pmutex->lock(); ///< Locks the mutex  pointed by `pmutex`.
#define _MUTEX_UNLOCK(pmutex) \
	pmutex->unlock(); ///<Unlocks the mutex  pointed by `pmutex`
#define _MUTEX_SAFE_LOCK(pmutex) \
	if (pmutex) \
		pmutex->lock();  ///< A safer version of `_MUTEX_LOCK()`.
#define _MUTEX_SAFE_UNLOCK(pmutex) \
	if (pmutex) \
		pmutex->unlock(); ///< A safer version of `_MUTEX_UNLOCK()`.
#define _MUTEX_ASSIGN_P(pmutex,value) \
	pmutex=value; ///< Assigns the value of `value` to mutex pointer `pmutex`.

///////////////////////////////// end for C++ 17 threads
#else
#ifdef _WINDOWS
#include <windows.h>
#include <process.h>
#undef max
#undef min

typedef HANDLE _THREAD; ///< Alias for `HANDLE`.
typedef CRITICAL_SECTION _MUTEX; ///< Alias for `CRITICAL_SECTION`.
typedef SRWLOCK _RWLOCK; ///< Alias for `SRWLOCK`.

#define _THREAD_CREATE(thread,start_thread,param) \
		{thread = (HANDLE)_beginthreadex(0,0,&(start_thread),param,0,0);} ///< Creates a new thread with start function `start_thread()`, that gets `param` as an argument.
#define _THREAD_JOIN(thread) \
		WaitForSingleObject(thread, INFINITE);  ///< Attaches `thread` to the calling thread.
#define _THREAD_CLOSE(thread) \
		CloseHandle(thread); ///< Threads created by `_beginthreadex()` must be closed.

// Macros for working with mutexes
#define _MUTEX_INIT(mutex) \
	InitializeCriticalSection(&mutex); ///< Initializes a critical section represented by `mutex`.
#define _MUTEX_DESTROY(mutex) \
	DeleteCriticalSection(&mutex); ///< Destroys the critical section represented by `mutex`.
#define _MUTEX_LOCK(pmutex) \
		EnterCriticalSection(pmutex);  ///< Enters critical section pointed by `pmutex`.
#define _MUTEX_UNLOCK(pmutex) \
	{LeaveCriticalSection(pmutex);}  ///< Leaves critical section pointed by `pmutex`.
#define _MUTEX_SAFE_LOCK(pmutex) \
	if (pmutex) \
		{EnterCriticalSection(pmutex);} // A safer version of `_MUTEX_LOCK()`.
#define _MUTEX_SAFE_UNLOCK(pmutex) \
	if (pmutex) \
		{LeaveCriticalSection(pmutex);}// A safer version of `_MUTEX_UNLOCK()`
#define _MUTEX_ASSIGN_P(pmutex,value) \
	pmutex=value; ///< Assigns the value of `value` to mutex pointer `pmutex`.

#define _RWLOCK_INIT(rwLock) \
	InitializeSRWLock(&rwLock); ///< Initializes the read-write lock given by `rwLock`.
#define _RWLOCK_DESTROY(rwLock) ///< Destroys the read-write-lock given by `rwLock`.

#define _RWLOCK_RDLOCK(rwLock) \
		AcquireSRWLockShared(rwLock); ///< Applies a read lock to the read-write lock pointed by `rwLock`.
#define _RWLOCK_SAFE_RDLOCK(rwLock) \
	if (rwLock) \
		AcquireSRWLockShared(rwLock); ///< A safer version of `_RWLOCK_RDLOCK()`.
#define _RWLOCK_WRLOCK(rwLock) \
		AcquireSRWLockExclusive(rwLock); ///< Applies a write lock to the read-write lock pointed by `rwLock`.
#define _RWLOCK_SAFE_WRLOCK(rwLock) \
	if (rwLock) \
		AcquireSRWLockExclusive(rwLock);  ///< A safer version of `_RWLOCK_WRLOCK()`.
#define _RWLOCK_UNLOCK_RDLOCK(rwLock) \
		ReleaseSRWLockShared(rwLock); ///< Unlocks the read-write lock pointed by `rwLock`.
#define _RWLOCK_UNLOCK_WRLOCK(rwLock) \
		ReleaseSRWLockExclusive(rwLock); ///< Unlocks the read-write lock pointed by `rwLock`.
#define _RWLOCK_SAFE_UNLOCK_RDLOCK(rwLock) \
	if (rwLock) \
	ReleaseSRWLockShared(rwLock);  ///< A safer version of `_RWLOCK_UNLOCK()`.
#define _RWLOCK_SAFE_UNLOCK_WRLOCK(rwLock) \
	if (rwLock) \
	ReleaseSRWLockExclusive(rwLock);  ///< A safer version of `_RWLOCK_UNLOCK()`.
#define _RWLOCK_ASSIGN_P(plock,pvalue) \
	plock=pvalue; ///< Assigns the value of `pvalue` to `_RWLOCK` pointer `plock`.
///////////////////////////////// end for WINDOWS threads
#else
#include <unistd.h>
#include <pthread.h>

typedef pthread_t _THREAD;  ///< Alias for `pthread_t`.
typedef pthread_mutex_t _MUTEX; ///< Alias for `pthread_mutex_t`.
typedef pthread_rwlock_t _RWLOCK; ///< Alias for `pthread_rwlock_t`.

#define _THREAD_CREATE(thread,start_thread,param) \
		pthread_create(&(thread),NULL,start_thread,param); ///< Creates a new thread with start function `start_thread()`, that gets `param` as an argument.
#define _THREAD_JOIN(thread) \
		pthread_join((thread),NULL); ///< Attaches `thread` to the calling thread.
#define _THREAD_CLOSE(thread) ///< Empty placeholder to be consistent with windows threads.
		
#define _RWLOCK_INIT(rwLock) \
	pthread_rwlock_init(&rwLock,NULL); ///< Initializes the read-write lock pointed by `rwLock`.
#define _RWLOCK_DESTROY(rwLock) \
	pthread_rwlock_destroy(&rwLock); ///< Destroys the read-write-lock pointed by `rwLock`.
#define _RWLOCK_RDLOCK(rwLock) \
		pthread_rwlock_rdlock(rwLock); ///< Applies a read lock to the read-write lock pointed by `rwLock`.
#define _RWLOCK_SAFE_RDLOCK(rwLock) \
	if (rwLock) \
		pthread_rwlock_rdlock(rwLock); ///< A safer version of `_RWLOCK_RDLOCK()`.
#define _RWLOCK_WRLOCK(rwLock) \
		pthread_rwlock_wrlock(rwLock); ///< Applies a write lock to the read-write lock pointed by `rwLock`.
#define _RWLOCK_SAFE_WRLOCK(rwLock) \
	if (rwLock) \
		pthread_rwlock_wrlock(rwLock);  ///< A safer version of `_RWLOCK_WRLOCK()`.
#define _RWLOCK_UNLOCK_RDLOCK(rwLock) \
		pthread_rwlock_unlock(rwLock); ///< Unlocks the read-write lock pointed by `rwLock`.
#define _RWLOCK_UNLOCK_WRLOCK(rwLock) \
		pthread_rwlock_unlock(rwLock); ///< Unlocks the read-write lock pointed by `rwLock`.
#define _RWLOCK_SAFE_UNLOCK_RDLOCK(rwLock) \
	if (rwLock) \
		pthread_rwlock_unlock(rwLock);  ///< A safer version of `_RWLOCK_UNLOCK()`.
#define _RWLOCK_SAFE_UNLOCK_WRLOCK(rwLock) \
	if (rwLock) \
		pthread_rwlock_unlock(rwLock);  ///< A safer version of `_RWLOCK_UNLOCK()`.
#define _RWLOCK_ASSIGN_P(plock,pvalue) \
	plock=pvalue; ///< Assigns the value of `pvalue` to `_RWLOCK` pointer `plock`.

#define _MUTEX_INIT(mutex) \
	pthread_mutex_init(&mutex,NULL); ///< Initializes the mutex represented by `mutex`.
#define _MUTEX_DESTROY(mutex) \
  	pthread_mutex_destroy(&mutex); ///< Destroys the mutex represented by `mutex`.
#define _MUTEX_LOCK(pmutex) \
	pthread_mutex_lock(pmutex); ///< Locks the mutex  pointed by `pmutex`.
#define _MUTEX_UNLOCK(pmutex) \
	pthread_mutex_unlock(pmutex); ///<Unlocks the mutex  pointed by `pmutex`
#define _MUTEX_SAFE_LOCK(pmutex) \
	if (pmutex) \
		pthread_mutex_lock(pmutex);  ///< A safer version of `_MUTEX_LOCK()`.
#define _MUTEX_SAFE_UNLOCK(pmutex) \
	if (pmutex) \
		pthread_mutex_unlock(pmutex); ///< A safer version of `_MUTEX_UNLOCK()`.
#define _MUTEX_ASSIGN_P(pmutex,value) \
	pmutex=value; ///< Assigns the value of `value` to mutex pointer `pmutex`.
#endif /* #elif _WINDOWS */
#endif /* #ifdef _CPP_THREADFS */
#else
//#define __MAX_THREAD_NUM 1
#define _RWLOCK_INIT(rwLock)
#define _RWLOCK_RDLOCK(rwLock)
#define _RWLOCK_SAFE_RDLOCK(rwLock)
#define _RWLOCK_WRLOCK(rwLock)
#define _RWLOCK_SAFE_WRLOCK(rwLock)
#define _RWLOCK_UNLOCK_RDLOCK(rwLock)
#define _RWLOCK_SAFE_UNLOCK_RDLOCK(rwLock)
#define _RWLOCK_UNLOCK_WRLOCK(rwLock)
#define _RWLOCK_SAFE_UNLOCK_WRLOCK(rwLock)

#define _MUTEX_ASSIGN_P(pmutex,value)
#define _MUTEX_LOCK(pmutex)
#define _MUTEX_UNLOCK(pmutex)
#define _MUTEX_SAFE_LOCK(pmutex)
#define _MUTEX_SAFE_UNLOCK(pmutex)
#endif /* #ifndef __ONE_THREAD_ */

#endif /*_THREAD_H_*/
