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

/* $Id: sse.c,v 1.15 2006-01-23 20:31:24 athena Exp $ */

#include "ifftw.h"
#include "simd.h"

#if HAVE_SSE

# if defined(__x86_64__) || defined(_M_X64) || defined(_M_AMD64)

  int RIGHT_CPU(void)
  {
       return 1;
  }

# else /* !x86_64 */

# include <signal.h>
# include <setjmp.h>

  static inline int cpuid_edx(int op)
  {
#    ifdef _MSC_VER
       int result;
       _asm {
	    push ebx
	    mov eax,op
	    cpuid
	    mov result,edx
	    pop ebx
       }
       return result;
#    else
       int eax, ecx, edx;

       __asm__("push %%ebx\n\tcpuid\n\tpop %%ebx"
	       : "=a" (eax), "=c" (ecx), "=d" (edx)
	       : "a" (op));
       return edx;
#    endif
  }

  static jmp_buf jb;

  static void sighandler(int x)
  {
       UNUSED(x);
       longjmp(jb, 1);
  }

  static int sse_works(void)
  {
       void (*oldsig)(int);
       oldsig = signal(SIGILL, sighandler);
       if (setjmp(jb)) {
	    signal(SIGILL, oldsig);
	    return 0;
       } else {
#        ifdef _MSC_VER
	    _asm { xorps xmm0,xmm0 }
#        else
	    /* asm volatile ("xorps %xmm0, %xmm0"); */
	    asm volatile (".byte 0x0f; .byte 0x57; .byte 0xc0");
#        endif
	    signal(SIGILL, oldsig);
	    return 1;
       }
  }

  int RIGHT_CPU(void)
  {
       static int init = 0, res;
       extern void X(check_alignment_of_sse_mpmp)(void);

       if (!init) {
	    res = (cpuid_edx(1) & (1 << 25)) && sse_works();
	    init = 1;
	    X(check_alignment_of_sse_mpmp)();
       }
       return res;
  }
# endif /* x86_64 */

#endif /* HAVE_SSE */
