/*
 * Copyright (c) 2003, 2006 Matteo Frigo
 * Copyright (c) 2003, 2006 Massachusetts Institute of Technology
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

/* $Id: altivec.c,v 1.14 2006-01-05 03:04:27 stevenj Exp $ */

#include "ifftw.h"
#include "simd.h"

#if HAVE_ALTIVEC

#if HAVE_SYS_SYSCTL_H
#include <sys/sysctl.h>
#endif

#if HAVE_SYS_SYSCTL_H && HAVE_SYSCTL && defined(CTL_HW) && defined(HW_VECTORUNIT)
/* code for darwin */
static int really_have_altivec(void)
{
     int mib[2], altivecp;
     size_t len;
     mib[0] = CTL_HW;
     mib[1] = HW_VECTORUNIT;
     len = sizeof(altivecp);
     sysctl(mib, 2, &altivecp, &len, NULL, 0);
     return altivecp;
} 
#else /* HAVE_SYS_SYSCTL_H etc. */

#include <signal.h>
#include <setjmp.h>

static jmp_buf jb;

static void sighandler(int x)
{
     longjmp(jb, 1);
}

static int really_have_altivec(void)
{
     void (*oldsig)(int);
     oldsig = signal(SIGILL, sighandler);
     if (setjmp(jb)) {
	  signal(SIGILL, oldsig);
	  return 0;
     } else {
	  __asm__ __volatile__ (".long 0x10000484"); /* vor 0,0,0 */
	  signal(SIGILL, oldsig);
	  return 1;
     }
     return 0;
}
#endif

int RIGHT_CPU(void)
{
     static int init = 0, res;
     if (!init) {
	  res = really_have_altivec();
	  init = 1;
     }
     return res;
}
#endif
