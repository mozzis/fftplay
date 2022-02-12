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
/* Generated on Fri Jan 27 20:29:35 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -twiddle-log3 -precompute-twiddles -n 4 -dit -name hf2_4 -include hf.h */

/*
 * This function contains 24 FP additions, 16 FP multiplications,
 * (or, 16 additions, 8 multiplications, 8 fused multiply/add),
 * 31 stack variables, and 16 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hf.h"

static const R *hf2_4(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 4, MAKE_VOLATILE_STRIDE(ios)) {
	  E T6, T3, T1, Tb, T7, Tx, Tc, Tm, Ts, T8, Tf, Th, T2, T5;
	  T2 = W[2];
	  T6 = W[1];
	  T3 = W[0];
	  T5 = W[3];
	  {
	       E Tj, Tl, Tk, Tr, Ta, T4;
	       T1 = rio[0];
	       Ta = T2 * T6;
	       T4 = T2 * T3;
	       Tj = rio[WS(ios, 3)];
	       Tl = iio[0];
	       Tb = FNMS(T5, T3, Ta);
	       T7 = FMA(T5, T6, T4);
	       Tk = T2 * Tj;
	       Tr = T2 * Tl;
	       Tx = iio[-WS(ios, 3)];
	       Tc = iio[-WS(ios, 1)];
	       Tm = FMA(T5, Tl, Tk);
	       Ts = FNMS(T5, Tj, Tr);
	       T8 = rio[WS(ios, 2)];
	       Tf = rio[WS(ios, 1)];
	       Th = iio[-WS(ios, 2)];
	  }
	  {
	       E Tv, T9, Tg, Tp;
	       Tv = Tb * T8;
	       T9 = T7 * T8;
	       Tg = T3 * Tf;
	       Tp = T3 * Th;
	       {
		    E Tw, Td, Ti, Tq;
		    Tw = FMA(T7, Tc, Tv);
		    Td = FNMS(Tb, Tc, T9);
		    Ti = FMA(T6, Th, Tg);
		    Tq = FNMS(T6, Tf, Tp);
		    {
			 E Ty, TA, Te, To;
			 Ty = Tw + Tx;
			 TA = Tx - Tw;
			 Te = T1 + Td;
			 To = T1 - Td;
			 {
			      E Tn, Tz, Tu, Tt;
			      Tn = Ti + Tm;
			      Tz = Tm - Ti;
			      Tu = Tq + Ts;
			      Tt = Tq - Ts;
			      iio[-WS(ios, 1)] = Tz + TA;
			      rio[WS(ios, 3)] = Tz - TA;
			      rio[0] = Te + Tn;
			      iio[-WS(ios, 2)] = Te - Tn;
			      rio[WS(ios, 1)] = To + Tt;
			      iio[-WS(ios, 3)] = To - Tt;
			      iio[0] = Tu + Ty;
			      rio[WS(ios, 2)] = Tu - Ty;
			 }
		    }
	       }
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 0, 1},
     {TW_CEXP, 0, 3},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 4, "hf2_4", twinstr, &GENUS, {16, 8, 8, 0}, 0, 0, 0 };

void X(codelet_hf2_4) (planner *p) {
     X(khc2hc_register) (p, hf2_4, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc -compact -variables 4 -pipeline-latency 4 -twiddle-log3 -precompute-twiddles -n 4 -dit -name hf2_4 -include hf.h */

/*
 * This function contains 24 FP additions, 16 FP multiplications,
 * (or, 16 additions, 8 multiplications, 8 fused multiply/add),
 * 21 stack variables, and 16 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hf.h"

static const R *hf2_4(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 4, MAKE_VOLATILE_STRIDE(ios)) {
	  E T2, T4, T3, T5, T6, T8;
	  T2 = W[2];
	  T4 = W[3];
	  T3 = W[0];
	  T5 = W[1];
	  T6 = FMA(T2, T3, T4 * T5);
	  T8 = FNMS(T4, T3, T2 * T5);
	  {
	       E T1, Th, Tl, Tp, Ta, To, Te, Tk, Tf, Tg;
	       T1 = rio[0];
	       Tf = rio[WS(ios, 3)];
	       Tg = iio[0];
	       Th = FMA(T2, Tf, T4 * Tg);
	       Tl = FNMS(T4, Tf, T2 * Tg);
	       Tp = iio[-WS(ios, 3)];
	       {
		    E T7, T9, Tc, Td;
		    T7 = rio[WS(ios, 2)];
		    T9 = iio[-WS(ios, 1)];
		    Ta = FNMS(T8, T9, T6 * T7);
		    To = FMA(T8, T7, T6 * T9);
		    Tc = rio[WS(ios, 1)];
		    Td = iio[-WS(ios, 2)];
		    Te = FMA(T3, Tc, T5 * Td);
		    Tk = FNMS(T5, Tc, T3 * Td);
	       }
	       {
		    E Tb, Ti, Tn, Tq;
		    Tb = T1 + Ta;
		    Ti = Te + Th;
		    iio[-WS(ios, 2)] = Tb - Ti;
		    rio[0] = Tb + Ti;
		    Tn = Tk + Tl;
		    Tq = To + Tp;
		    rio[WS(ios, 2)] = Tn - Tq;
		    iio[0] = Tn + Tq;
	       }
	       {
		    E Tj, Tm, Tr, Ts;
		    Tj = T1 - Ta;
		    Tm = Tk - Tl;
		    iio[-WS(ios, 3)] = Tj - Tm;
		    rio[WS(ios, 1)] = Tj + Tm;
		    Tr = Th - Te;
		    Ts = Tp - To;
		    rio[WS(ios, 3)] = Tr - Ts;
		    iio[-WS(ios, 1)] = Tr + Ts;
	       }
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 0, 1},
     {TW_CEXP, 0, 3},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 4, "hf2_4", twinstr, &GENUS, {16, 8, 8, 0}, 0, 0, 0 };

void X(codelet_hf2_4) (planner *p) {
     X(khc2hc_register) (p, hf2_4, &desc);
}
#endif				/* HAVE_FMA */
