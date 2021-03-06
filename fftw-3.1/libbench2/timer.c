/*
 * Copyright (c) 2001 Matteo Frigo
 * Copyright (c) 2001 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* $Id: timer.c,v 1.12 2006-01-15 21:32:54 athena Exp $ */

#include "bench.h"
#include <stdio.h>

/* 
 * System-dependent timing functions:
 */
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_BSDGETTIMEOFDAY
#ifndef HAVE_GETTIMEOFDAY
#define gettimeofday BSDgettimeofday
#define HAVE_GETTIMEOFDAY 1
#endif
#endif

double time_min;
int time_repeat;

#if !defined(HAVE_TIMER) && (defined(__WIN32__) || defined(_WIN32) || defined(_WINDOWS) || defined(__CYGWIN__))
#include <windows.h>
typedef LARGE_INTEGER mytime;

static mytime get_time(void)
{
     mytime tv;
     QueryPerformanceCounter(&tv);
     return tv;
}

static double elapsed(mytime t1, mytime t0)
{
     LARGE_INTEGER freq;
     QueryPerformanceFrequency(&freq);
     return ((double) (t1.QuadPart - t0.QuadPart)) /
	  ((double) freq.QuadPart);
}

#define HAVE_TIMER
#endif


#if defined(HAVE_GETTIMEOFDAY) && !defined(HAVE_TIMER)
typedef struct timeval mytime;

static mytime get_time(void)
{
     struct timeval tv;
     gettimeofday(&tv, 0);
     return tv;
}

static double elapsed(mytime t1, mytime t0)
{
     return (double)(t1.tv_sec - t0.tv_sec) +
	  (double)(t1.tv_usec - t0.tv_usec) * 1.0E-6;
}

#define HAVE_TIMER
#endif

#ifndef HAVE_TIMER
#error "timer not defined"
#endif

static const double tmax_try = 1.0;    /* seconds */
static const int nmin = 128, nmax = 133;

static double time_one(int n)
{
     float X[16], Y[16];
     int i;
     mytime t0, t1;

     for (i = 0; i < 16; ++i)
	  X[i] = 0;

     t0 = get_time();
     for (i = 0; i < n; ++i)
	  bench_fft8(X, X+1, Y, Y+1, 2, 2);
     t1 = get_time();
     return (elapsed(t1, t0));
}

static double time_n(int n)
{
     int     i;
     double  tmin;

     tmin = time_one(n);
     for (i = 1; i < time_repeat; ++i) {
	  double t = time_one(n);
	  if (t < tmin)
	       tmin = t;
     }
     return tmin;
}

static int good_enough_p(int n, double *tp)
{
     int i;
     double t = 0.0;

     /* vary nmin and see if time scales proportionally */
     for (i = nmin; i < nmax; ++i) {
	  double t1 = time_n(n * i);

	  if (t1 <= 0)
	       return 0; /* not enough resolution */

	  if (t1 >= tmax_try) {
	       t = t1;
	       break;
	  }

	  t = (i == nmin) ? t1 : t;

	  if (t1 >= (t * (i + 0.5) / nmin))
	       return 0;

	  if (t1 <= (t * (i - 0.5) / nmin))
	       return 0;
     }

     *tp = t;
     return 1;
}

static double calibrate(void)
{
     double t = tmax_try;
     int n;

     for (n = 1; n < (1 << 20); n += n) 
	  if (good_enough_p(n, &t))
	       break;

     return t;
}


void timer_init(double tmin, int repeat)
{
     static int inited = 0;

     if (inited)
	  return;
     inited = 1;

     if (!repeat)
	  repeat = 8;
     time_repeat = repeat;

     if (tmin > 0)
	  time_min = tmin;
     else
	  time_min = calibrate();
}

static mytime t0[NTIMERS];

void timer_start(int n)
{
     BENCH_ASSERT(n >= 0 && n < NTIMERS);
     t0[n] = get_time();
}

double timer_stop(int n)
{
     mytime t1;
     BENCH_ASSERT(n >= 0 && n < NTIMERS);
     t1 = get_time();
     return elapsed(t1, t0[n]);
}

