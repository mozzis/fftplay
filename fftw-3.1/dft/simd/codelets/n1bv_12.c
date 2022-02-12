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
/* Generated on Fri Jan 27 19:44:00 EST 2006 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_notw_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -sign 1 -n 12 -name n1bv_12 -include n1b.h */

/*
 * This function contains 48 FP additions, 20 FP multiplications,
 * (or, 30 additions, 2 multiplications, 18 fused multiply/add),
 * 49 stack variables, and 24 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_notw_c.ml,v 1.16 2006-01-05 03:04:27 stevenj Exp $
 */

#include "n1b.h"

static void n1bv_12(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     const R *xi;
     R *xo;
     xi = ii;
     xo = io;
     for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  V T1, T6, Tc, Th, Td, Te, Ti, Tz, T4, TA, T9, Tj, Tf, Tw;
	  {
	       V T2, T3, T7, T8;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       T6 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
	       T2 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
	       T3 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
	       T7 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
	       T8 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
	       Tc = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
	       Th = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
	       Td = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
	       Te = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
	       Ti = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
	       Tz = VSUB(T2, T3);
	       T4 = VADD(T2, T3);
	       TA = VSUB(T7, T8);
	       T9 = VADD(T7, T8);
	       Tj = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
	  }
	  Tf = VADD(Td, Te);
	  Tw = VSUB(Td, Te);
	  {
	       V T5, Tp, TJ, TB, Ta, Tq, Tk, Tx, Tg, Ts;
	       T5 = VADD(T1, T4);
	       Tp = VFNMS(LDK(KP500000000), T4, T1);
	       TJ = VSUB(Tz, TA);
	       TB = VADD(Tz, TA);
	       Ta = VADD(T6, T9);
	       Tq = VFNMS(LDK(KP500000000), T9, T6);
	       Tk = VADD(Ti, Tj);
	       Tx = VSUB(Tj, Ti);
	       Tg = VADD(Tc, Tf);
	       Ts = VFNMS(LDK(KP500000000), Tf, Tc);
	       {
		    V Tr, TF, Tb, Tn, TG, Ty, Tl, Tt;
		    Tr = VADD(Tp, Tq);
		    TF = VSUB(Tp, Tq);
		    Tb = VSUB(T5, Ta);
		    Tn = VADD(T5, Ta);
		    TG = VADD(Tw, Tx);
		    Ty = VSUB(Tw, Tx);
		    Tl = VADD(Th, Tk);
		    Tt = VFNMS(LDK(KP500000000), Tk, Th);
		    {
			 V TC, TE, TH, TL, Tu, TI, Tm, To;
			 TC = VMUL(LDK(KP866025403), VSUB(Ty, TB));
			 TE = VMUL(LDK(KP866025403), VADD(TB, Ty));
			 TH = VFNMS(LDK(KP866025403), TG, TF);
			 TL = VFMA(LDK(KP866025403), TG, TF);
			 Tu = VADD(Ts, Tt);
			 TI = VSUB(Ts, Tt);
			 Tm = VSUB(Tg, Tl);
			 To = VADD(Tg, Tl);
			 {
			      V TK, TM, Tv, TD;
			      TK = VFMA(LDK(KP866025403), TJ, TI);
			      TM = VFNMS(LDK(KP866025403), TJ, TI);
			      Tv = VSUB(Tr, Tu);
			      TD = VADD(Tr, Tu);
			      ST(&(xo[0]), VADD(Tn, To), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 6)]), VSUB(Tn, To), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 9)]), VFMAI(Tm, Tb), ovs, &(xo[WS(os, 1)]));
			      ST(&(xo[WS(os, 3)]), VFNMSI(Tm, Tb), ovs, &(xo[WS(os, 1)]));
			      ST(&(xo[WS(os, 5)]), VFMAI(TM, TL), ovs, &(xo[WS(os, 1)]));
			      ST(&(xo[WS(os, 7)]), VFNMSI(TM, TL), ovs, &(xo[WS(os, 1)]));
			      ST(&(xo[WS(os, 11)]), VFNMSI(TK, TH), ovs, &(xo[WS(os, 1)]));
			      ST(&(xo[WS(os, 1)]), VFMAI(TK, TH), ovs, &(xo[WS(os, 1)]));
			      ST(&(xo[WS(os, 8)]), VFNMSI(TE, TD), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 4)]), VFMAI(TE, TD), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 2)]), VFMAI(TC, Tv), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 10)]), VFNMSI(TC, Tv), ovs, &(xo[0]));
			 }
		    }
	       }
	  }
     }
}

static const kdft_desc desc = { 12, "n1bv_12", {30, 2, 18, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1bv_12) (planner *p) {
     X(kdft_register) (p, n1bv_12, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_notw_c -simd -compact -variables 4 -pipeline-latency 8 -sign 1 -n 12 -name n1bv_12 -include n1b.h */

/*
 * This function contains 48 FP additions, 8 FP multiplications,
 * (or, 44 additions, 4 multiplications, 4 fused multiply/add),
 * 27 stack variables, and 24 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_notw_c.ml,v 1.16 2006-01-05 03:04:27 stevenj Exp $
 */

#include "n1b.h"

static void n1bv_12(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     const R *xi;
     R *xo;
     xi = ii;
     xo = io;
     for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	  V T5, Ta, TG, TF, Ty, Tm, Ti, Tp, TJ, TI, Tx, Ts;
	  {
	       V T1, T6, T4, Tk, T9, Tl;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       T6 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
	       {
		    V T2, T3, T7, T8;
		    T2 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
		    T3 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
		    T4 = VADD(T2, T3);
		    Tk = VSUB(T2, T3);
		    T7 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
		    T8 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    T9 = VADD(T7, T8);
		    Tl = VSUB(T7, T8);
	       }
	       T5 = VFNMS(LDK(KP500000000), T4, T1);
	       Ta = VFNMS(LDK(KP500000000), T9, T6);
	       TG = VADD(T6, T9);
	       TF = VADD(T1, T4);
	       Ty = VADD(Tk, Tl);
	       Tm = VMUL(LDK(KP866025403), VSUB(Tk, Tl));
	  }
	  {
	       V Tn, Tq, Te, To, Th, Tr;
	       Tn = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
	       Tq = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
	       {
		    V Tc, Td, Tf, Tg;
		    Tc = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
		    Td = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
		    Te = VSUB(Tc, Td);
		    To = VADD(Tc, Td);
		    Tf = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    Tg = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    Th = VSUB(Tf, Tg);
		    Tr = VADD(Tf, Tg);
	       }
	       Ti = VMUL(LDK(KP866025403), VSUB(Te, Th));
	       Tp = VFNMS(LDK(KP500000000), To, Tn);
	       TJ = VADD(Tq, Tr);
	       TI = VADD(Tn, To);
	       Tx = VADD(Te, Th);
	       Ts = VFNMS(LDK(KP500000000), Tr, Tq);
	  }
	  {
	       V TH, TK, TL, TM;
	       TH = VSUB(TF, TG);
	       TK = VBYI(VSUB(TI, TJ));
	       ST(&(xo[WS(os, 3)]), VSUB(TH, TK), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[WS(os, 9)]), VADD(TH, TK), ovs, &(xo[WS(os, 1)]));
	       TL = VADD(TF, TG);
	       TM = VADD(TI, TJ);
	       ST(&(xo[WS(os, 6)]), VSUB(TL, TM), ovs, &(xo[0]));
	       ST(&(xo[0]), VADD(TL, TM), ovs, &(xo[0]));
	  }
	  {
	       V Tj, Tv, Tu, Tw, Tb, Tt;
	       Tb = VSUB(T5, Ta);
	       Tj = VSUB(Tb, Ti);
	       Tv = VADD(Tb, Ti);
	       Tt = VSUB(Tp, Ts);
	       Tu = VBYI(VADD(Tm, Tt));
	       Tw = VBYI(VSUB(Tt, Tm));
	       ST(&(xo[WS(os, 11)]), VSUB(Tj, Tu), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[WS(os, 5)]), VADD(Tv, Tw), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[WS(os, 1)]), VADD(Tj, Tu), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[WS(os, 7)]), VSUB(Tv, Tw), ovs, &(xo[WS(os, 1)]));
	  }
	  {
	       V Tz, TD, TC, TE, TA, TB;
	       Tz = VBYI(VMUL(LDK(KP866025403), VSUB(Tx, Ty)));
	       TD = VBYI(VMUL(LDK(KP866025403), VADD(Ty, Tx)));
	       TA = VADD(T5, Ta);
	       TB = VADD(Tp, Ts);
	       TC = VSUB(TA, TB);
	       TE = VADD(TA, TB);
	       ST(&(xo[WS(os, 2)]), VADD(Tz, TC), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 8)]), VSUB(TE, TD), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 10)]), VSUB(TC, Tz), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 4)]), VADD(TD, TE), ovs, &(xo[0]));
	  }
     }
}

static const kdft_desc desc = { 12, "n1bv_12", {44, 4, 4, 0}, &GENUS, 0, 0, 0, 0 };
void X(codelet_n1bv_12) (planner *p) {
     X(kdft_register) (p, n1bv_12, &desc);
}

#endif				/* HAVE_FMA */
