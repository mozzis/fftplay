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
/* Generated on Fri Jan 27 20:42:11 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 4 -dif -name hb_4 -include hb.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 16 additions, 6 multiplications, 6 fused multiply/add),
 * 25 stack variables, and 16 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hb.h"

static const R *hb_4(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 6, MAKE_VOLATILE_STRIDE(ios)) {
	  E Th, Ta, T7, Ti, T9;
	  {
	       E Tq, Td, T3, Tg, Tx, Tm, T6, Tp;
	       {
		    E Tk, T4, Tl, T5;
		    {
			 E Tb, Tc, T1, T2, Te, Tf;
			 Tb = iio[0];
			 Tc = rio[WS(ios, 2)];
			 T1 = rio[0];
			 T2 = iio[-WS(ios, 2)];
			 Te = iio[-WS(ios, 1)];
			 Tq = Tb + Tc;
			 Td = Tb - Tc;
			 Tf = rio[WS(ios, 3)];
			 Tk = T1 - T2;
			 T3 = T1 + T2;
			 T4 = rio[WS(ios, 1)];
			 Tg = Te - Tf;
			 Tl = Te + Tf;
			 T5 = iio[-WS(ios, 3)];
		    }
		    Tx = Tk + Tl;
		    Tm = Tk - Tl;
		    T6 = T4 + T5;
		    Tp = T4 - T5;
	       }
	       iio[-WS(ios, 3)] = Td + Tg;
	       {
		    E Tu, Tr, T8, Ty, Tv, Tw, Tt;
		    Tt = W[4];
		    Tu = Tq - Tp;
		    Tr = Tp + Tq;
		    rio[0] = T3 + T6;
		    T8 = T3 - T6;
		    Ty = Tt * Tx;
		    Tv = Tt * Tu;
		    Tw = W[5];
		    {
			 E Tj, To, Ts, Tn;
			 Tj = W[0];
			 To = W[1];
			 Th = Td - Tg;
			 rio[WS(ios, 3)] = FNMS(Tw, Tu, Ty);
			 iio[0] = FMA(Tw, Tx, Tv);
			 Ts = Tj * Tr;
			 Tn = Tj * Tm;
			 Ta = W[3];
			 T7 = W[2];
			 iio[-WS(ios, 2)] = FMA(To, Tm, Ts);
			 rio[WS(ios, 1)] = FNMS(To, Tr, Tn);
			 Ti = Ta * T8;
			 T9 = T7 * T8;
		    }
	       }
	  }
	  iio[-WS(ios, 1)] = FMA(T7, Th, Ti);
	  rio[WS(ios, 2)] = FNMS(Ta, Th, T9);
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 4},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 4, "hb_4", twinstr, &GENUS, {16, 6, 6, 0}, 0, 0, 0 };

void X(codelet_hb_4) (planner *p) {
     X(khc2hc_register) (p, hb_4, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc -compact -variables 4 -pipeline-latency 4 -sign 1 -n 4 -dif -name hb_4 -include hb.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 16 additions, 6 multiplications, 6 fused multiply/add),
 * 13 stack variables, and 16 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hb.h"

static const R *hb_4(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 6, MAKE_VOLATILE_STRIDE(ios)) {
	  E T3, Ti, Tc, Tn, T6, Tm, Tf, Tj;
	  {
	       E T1, T2, Ta, Tb;
	       T1 = rio[0];
	       T2 = iio[-WS(ios, 2)];
	       T3 = T1 + T2;
	       Ti = T1 - T2;
	       Ta = iio[0];
	       Tb = rio[WS(ios, 2)];
	       Tc = Ta - Tb;
	       Tn = Ta + Tb;
	  }
	  {
	       E T4, T5, Td, Te;
	       T4 = rio[WS(ios, 1)];
	       T5 = iio[-WS(ios, 3)];
	       T6 = T4 + T5;
	       Tm = T4 - T5;
	       Td = iio[-WS(ios, 1)];
	       Te = rio[WS(ios, 3)];
	       Tf = Td - Te;
	       Tj = Td + Te;
	  }
	  rio[0] = T3 + T6;
	  iio[-WS(ios, 3)] = Tc + Tf;
	  {
	       E Tq, Ts, Tp, Tr;
	       Tq = Tn - Tm;
	       Ts = Ti + Tj;
	       Tp = W[4];
	       Tr = W[5];
	       iio[0] = FMA(Tp, Tq, Tr * Ts);
	       rio[WS(ios, 3)] = FNMS(Tr, Tq, Tp * Ts);
	  }
	  {
	       E T8, Tg, T7, T9;
	       T8 = T3 - T6;
	       Tg = Tc - Tf;
	       T7 = W[2];
	       T9 = W[3];
	       rio[WS(ios, 2)] = FNMS(T9, Tg, T7 * T8);
	       iio[-WS(ios, 1)] = FMA(T9, T8, T7 * Tg);
	  }
	  {
	       E Tk, To, Th, Tl;
	       Tk = Ti - Tj;
	       To = Tm + Tn;
	       Th = W[0];
	       Tl = W[1];
	       rio[WS(ios, 1)] = FNMS(Tl, To, Th * Tk);
	       iio[-WS(ios, 2)] = FMA(Th, To, Tl * Tk);
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 4},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 4, "hb_4", twinstr, &GENUS, {16, 6, 6, 0}, 0, 0, 0 };

void X(codelet_hb_4) (planner *p) {
     X(khc2hc_register) (p, hb_4, &desc);
}
#endif				/* HAVE_FMA */
