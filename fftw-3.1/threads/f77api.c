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

#include "api.h"

/* if F77_FUNC is not defined, then we don't know how to mangle identifiers
   for the Fortran linker, and we must omit the f77 API. */
#if defined(F77_FUNC) || defined(WINDOWS_F77_MANGLING)

#include "x77.h"

#define F77(a, A) F77x(x77(a), X77(A))

#ifndef WINDOWS_F77_MANGLING

#if defined(F77_FUNC)
#  define F77x(a, A) F77_FUNC(a, A)
#  include "f77funcs.h"
#endif

#if defined(F77_FUNC_) && !defined(F77_FUNC_EQUIV)
#  undef F77x
#  define F77x(a, A) F77_FUNC_(a, A)
#  include "f77funcs.h"
#endif

#else /* WINDOWS_F77_MANGLING */

#  define WINDOWS_F77_FUNC(a, A) a ## __
#  define F77x(a, A) WINDOWS_F77_FUNC(a, A)
#  include "f77funcs.h"

#  undef WINDOWS_F77_FUNC
#  if defined(__GNUC__)
#    define WINDOWS_F77_FUNC(a, A) __attribute__((stdcall)) A
#  elif defined(_MSC_VER) || defined(_ICC) || defined(_STDCALL_SUPPORTED)
#    define WINDOWS_F77_FUNC(a, A) __stdcall A
#  else
#    define WINDOWS_F77_FUNC(a, A) A /* oh well */
#  endif
#  include "f77funcs.h"

#endif /* WINDOWS_F77_MANGLING */

#endif				/* F77_FUNC */
