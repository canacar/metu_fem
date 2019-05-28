/* $Id: hptimer.h,v 1.2 2008/01/28 07:36:28 canacar Exp $ */
/* 
 * This file is part of the EMSI Tools Package developed at the
 * Brain Research Laboratory, Middle East Technical University
 * Department of Electrical and Electronics Engineering.
 *
 * Copyright (C) 2008 Zeynep Akalin Acar
 * Copyright (C) 2008 Can Erkin Acar
 * Copyright (C) 2008 Nevzat G. Gencer
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
//---------------------------------------------------------------------------
#ifndef hptimerH
#define hptimerH
//---------------------------------------------------------------------------
#define HPT_NAME_LEN 64
#define HPT_BUF_LEN 32
//---------------------------------------------------------------------------
#ifdef _WIN32_

// These are from winnt.h & windefs.h
#ifndef _WINNT_
#define MAXLONGLONG                      (0x7fffffffffffffff)
typedef struct _LARGE_INTEGER {
    __int64 QuadPart;
} LARGE_INTEGER;
extern "C"{
__stdcall QueryPerformanceCounter(LARGE_INTEGER *lpPerfCnt);
__stdcall QueryPerformanceFrequency(LARGE_INTEGER *lpPerfCnt);
}
#endif
#else
#include <sys/time.h>
#include <unistd.h>

typedef long long int __int64;
#define MAXLONGLONG                      (0x7fffffffffffffffLL)
typedef struct _LARGE_INTEGER {
    __int64 QuadPart;
} LARGE_INTEGER;
extern "C"{
	extern struct timeval _hpt_tval;
	extern struct timezone _hpt_tzone;
	inline int QueryPerformanceCounter(LARGE_INTEGER *lpPerfCnt){
		gettimeofday(&_hpt_tval, &_hpt_tzone);
		lpPerfCnt->QuadPart=(__int64)_hpt_tval.tv_sec*(__int64)1000000+(__int64)_hpt_tval.tv_usec;
		return 1;
	}
	inline int QueryPerformanceFrequency(LARGE_INTEGER *lpPerfCnt){
		lpPerfCnt->QuadPart=1000000;
		return 1;
	}
}
#endif
//---------------------------------------------------------------------------
class HPTimer;

class HPTimerMgr{
public:
        ~HPTimerMgr();

        static HPTimer *getTimer(const char *name);

        static HPTimerMgr *getManager(void);
        static void report(void);
        static void report(char *str, int len);

        inline int isOK(void){return m_freq!=0;}
        inline __int64 getFrequency(void){return m_freq;}
        inline HPTimer *getFirst(void){return m_root;}

        static char *dispFreq(double t);
        static char *dispTime(double t);

private:
        friend class HPTimer;
        static HPTimerMgr *m_instance;
        HPTimerMgr();           // singleton
        void addTimer(HPTimer *t);
        void removeTimer(HPTimer *t); // called by timer destructor

        __int64 m_freq;
        HPTimer *m_root, *m_end; // so that we can add to the end of list
};

class HPTimer{
public:
        ~HPTimer();
        inline void start(void){
                QueryPerformanceCounter(&m_start);
        }
        void stop(void){
                QueryPerformanceCounter(&m_stop);
                m_cur=m_stop.QuadPart-m_start.QuadPart;
                if(m_cur<0) return;
                m_total+=m_cur;
                m_count++;
                if(m_min>m_cur) m_min=m_cur;
                if(m_max<m_cur) m_max=m_cur;
        }
        void reset(void);
        void report(void);
        int report(char *buf, int len);

        inline int getCount(void){return m_count;}
        inline __int64 getTotal(void){return m_total;}
        inline __int64 getMin(void){return m_min;}
        inline __int64 getMax(void){return m_max;}

        inline double getTotalSec(void){
                return m_parent->m_freq ? ((double)m_total/(double)m_parent->m_freq):0;}
        inline double getMinSec(void){
                return m_parent->m_freq ? ((double)m_min/(double)m_parent->m_freq):0;}
        inline double getMaxSec(void){
                return m_parent->m_freq ? ((double)m_max/(double)m_parent->m_freq):0;}

        inline const char *getName(void){return m_name;}
        inline HPTimer *next(void){return m_next;}

private:
        // timer data
        LARGE_INTEGER m_start, m_stop;
        __int64 m_cur, m_min, m_max, m_total;
        int m_count;
        char m_name[HPT_NAME_LEN];

       // managed by timer manager
        friend class HPTimerMgr;
        HPTimer(HPTimerMgr *parent, const char *name);
        HPTimerMgr *m_parent;
        HPTimer *m_next;

};
#endif

