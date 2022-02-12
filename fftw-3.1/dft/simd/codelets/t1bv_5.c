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

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Fri Jan 27 19:59:28 EST 2006 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -n 5 -name t1bv_5 -include t1b.h -sign 1 */

/*
 * This function contains 20 FP additions, 19 FP multiplications,
 * (or, 11 additions, 10 multiplications, 9 fused multiply/add),
 * 26 stack variables, and 10 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle_c.ml,v 1.13 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t1b.h"

static const R *t1bv_5(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     DVK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DVK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DVK(KP618033988, +0.618033988749894848204586834365638117720309180);
     DVK(KP951056516, +0.951056516295153572116439333379382143405698634);
     INT i;
     R *x;
     x = ii;
     for (i = m; i > 0; i = i - VL, x = x + (VL * dist), W = W + (TWVL * 8), MAKE_VOLATILE_STRIDE(ios)) {
	  V T1, T2, T9, T4, T7;
	  T1 = LD(&(x[0]), dist, &(x[0]));
	  T2 = LD(&(x[WS(ios, 1)]), dist, &(x[WS(ios, 1)]));
	  T9 = LD(&(x[WS(ios, 3)]), dist, &(x[WS(ios, 1)]));
	  T4 = LD(&(x[WS(ios, 4)]), dist, &(x[0]));
	  T7 = LD(&(x[WS(ios, 2)]), dist, &(x[0]));
	  {
	       V T3, Ta, T5, T8;
	       T3 = BYTW(&(W[0]), T2);
	       Ta = BYTW(&(W[TWVL * 4]), T9);
	       T5 = BYTW(&(W[TWVL * 6]), T4);
	       T8 = BYTW(&(W[TWVL * 2]), T7);
	       {
		    V T6, Tg, Tb, Th;
		    T6 = VADD(T3, T5);
		    Tg = VSUB(T3, T5);
		    Tb = VADD(T8, Ta);
		    Th = VSUB(T8, Ta);
		    {
			 V Te, Tc, Tk, Ti, Td, Tj, Tf;
			 Te = VSUB(T6, Tb);
			 Tc = VADD(T6, Tb);
			 Tk = VMUL(LDK(KP951056516), VFNMS(LDK(KP618033988), Tg, Th));
			 Ti = VMUL(LDK(KP951056516), VFMA(LDK(KP618033988), Th, Tg));
			 Td = VFNMS(LDK(KP250000000), Tc, T1);
			 ST(&(x[0]), VADD(T1, Tc), dist, &(x[0]));
			 Tj = VFNMS(LDK(KP559016994), Te, Td);
			 Tf = VFMA(LDK(KP559016994), Te, Td);
			 ST(&(x[WS(ios, 2)]), VFNMSI(Tk, Tj), dist, &(x[0]));
			 ST(&(x[WS(ios, 3)]), VFMAI(Tk, Tj), dist, &(x[WS(ios, 1)]));
			 ST(&(x[WS(ios, 4)]), VFNMSI(Ti, Tf), dist, &(x[0]));
			 ST(&(x[WS(ios, 1)]), VFMAI(Ti, Tf), dist, &(x[WS(ios, 1)]));
		    }
	       }
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     VTW(1),
     VTW(2),
     VTW(3),
     VTW(4),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 5, "t1bv_5", twinstr, &GENUS, {11, 10, 9, 0}, 0, 0, 0 };

void X(codelet_t1bv_5) (planner *p) {
     X(kdft_dit_register) (p, t1bv_5, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle_c -simd -compact -variables 4 -pipeline-latency 8 -n 5 -name t1bv_5 -include t1b.h -sign 1 */

/*
 * This function contains 20 FP additions, 14 FP multiplications,
 * (or, 17 additions, 11 multiplications, 3 fused multiply/add),
 * 20 stack variables, and 10 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle_c.ml,v 1.13 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t1b.h"

static const R *t1bv_5(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     DVK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DVK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DVK(KP587785252, +0.587785252292473129168705954639072768597652438);
     DVK(KP951056516, +0.951056516295153572116439333379382143405698634);
     INT i;
     R *x;
     x = ii;
     for (i = m; i > 0; i = i - VL, x = x + (VL * dist), W = W + (TWVL * 8), MAKE_VOLATILE_STRIDE(ios)) {
	  V Tf, T5, Ta, Tc, Td, Tg;
	  Tf = LD(&(x[0]), dist, &(x[0]));
	  {
	       V T2, T9, T4, T7;
	       {
		    V T1, T8, T3, T6;
		    T1 = LD(&(x[WS(ios, 1)]), dist, &(x[WS(ios, 1)]));
		    T2 = BYTW(&(W[0]), T1);
		    T8 = LD(&(x[WS(ios, 3)]), dist, &(x[WS(ios, 1)]));
		    T9 = BYTW(&(W[TWVL * 4]), T8);
		    T3 = LD(&(x[WS(ios, 4)]), dist, &(x[0]));
		    T4 = BYTW(&(W[TWVL * 6]), T3);
		    T6 = LD(&(x[WS(ios, 2)]), dist, &(x[0]));
		    T7 = BYTW(&(W[TWVL * 2]), T6);
	       }
	       T5 = VSUB(T2, T4);
	       Ta = VSUB(T7, T9);
	       Tc = VADD(T2, T4);
	       Td = VADD(T7, T9);
	       Tg = VADD(Tc, Td);
	  }
	  ST(&(x[0]), VADD(Tf, Tg), dist, &(x[0]));
	  {
	       V Tb, Tj, Ti, Tk, Te, Th;
	       Tb = VBYI(VFMA(LDK(KP951056516), T5, VMUL(LDK(KP587785252), Ta)));
	       Tj = VBYI(VFNMS(LDK(KP951056516), Ta, VMUL(LDK(KP587785252), T5)));
	       Te = VMUL(LDK(KP559016994), VSUB(Tc, Td));
	       Th = VFNMS(LDK(KP250000000), Tg, Tf);
	       Ti = VADD(Te, Th);
	       Tk = VSUB(Th, Te);
	       ST(&(x[WS(ios, 1)]), VADD(Tb, Ti), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 3)]), VSUB(Tk, Tj), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 4)]), VSUB(Ti, Tb), dist, &(x[0]));
	       ST(&(x[WS(ios, 2)]), VADD(Tj, Tk), dist, &(x[0]));
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     VTW(1),
     VTW(2),
     VTW(3),
     VTW(4),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 5, "t1bv_5", twinstr, &GENUS, {17, 11, 3, 0}, 0, 0, 0 };

void X(codelet_t1bv_5) (planner *p) {
     X(kdft_dit_register) (p, t1bv_5, &desc);
}
#endif				/* HAVE_FMA */