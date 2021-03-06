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
/* Generated on Fri Jan 27 19:58:29 EST 2006 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -n 16 -name t2fv_16 -include t2f.h */

/*
 * This function contains 87 FP additions, 64 FP multiplications,
 * (or, 53 additions, 30 multiplications, 34 fused multiply/add),
 * 61 stack variables, and 32 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle_c.ml,v 1.13 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t2f.h"

static const R *t2fv_16(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     R *x;
     x = ri;
     for (i = m; i > 0; i = i - VL, x = x + (VL * dist), W = W + (TWVL * 30), MAKE_VOLATILE_STRIDE(ios)) {
	  V TO, Ta, TJ, TP, T14, Tq, T1i, T10, T1b, T1l, T13, T1c, TR, Tl, T15;
	  V Tv;
	  {
	       V Tc, TW, T4, T19, T9, TD, TI, Tj, TZ, T1a, Te, Th, Tn, Tr, Tu;
	       V Tp;
	       {
		    V T1, T2, T5, T7;
		    T1 = LD(&(x[0]), dist, &(x[0]));
		    T2 = LD(&(x[WS(ios, 8)]), dist, &(x[0]));
		    T5 = LD(&(x[WS(ios, 4)]), dist, &(x[0]));
		    T7 = LD(&(x[WS(ios, 12)]), dist, &(x[0]));
		    {
			 V Tz, TG, TB, TE;
			 Tz = LD(&(x[WS(ios, 14)]), dist, &(x[0]));
			 TG = LD(&(x[WS(ios, 10)]), dist, &(x[0]));
			 TB = LD(&(x[WS(ios, 6)]), dist, &(x[0]));
			 TE = LD(&(x[WS(ios, 2)]), dist, &(x[0]));
			 {
			      V Ti, TY, TX, Td, Tg, Tm, Tt, To;
			      {
				   V T3, T6, T8, TA, TH, TC, TF, Tb;
				   Tb = LD(&(x[WS(ios, 1)]), dist, &(x[WS(ios, 1)]));
				   T3 = BYTWJ(&(W[TWVL * 14]), T2);
				   T6 = BYTWJ(&(W[TWVL * 6]), T5);
				   T8 = BYTWJ(&(W[TWVL * 22]), T7);
				   TA = BYTWJ(&(W[TWVL * 26]), Tz);
				   TH = BYTWJ(&(W[TWVL * 18]), TG);
				   TC = BYTWJ(&(W[TWVL * 10]), TB);
				   TF = BYTWJ(&(W[TWVL * 2]), TE);
				   Tc = BYTWJ(&(W[0]), Tb);
				   TW = VSUB(T1, T3);
				   T4 = VADD(T1, T3);
				   T19 = VSUB(T6, T8);
				   T9 = VADD(T6, T8);
				   Ti = LD(&(x[WS(ios, 13)]), dist, &(x[WS(ios, 1)]));
				   TD = VADD(TA, TC);
				   TY = VSUB(TA, TC);
				   TI = VADD(TF, TH);
				   TX = VSUB(TF, TH);
			      }
			      Td = LD(&(x[WS(ios, 9)]), dist, &(x[WS(ios, 1)]));
			      Tg = LD(&(x[WS(ios, 5)]), dist, &(x[WS(ios, 1)]));
			      Tm = LD(&(x[WS(ios, 15)]), dist, &(x[WS(ios, 1)]));
			      Tj = BYTWJ(&(W[TWVL * 24]), Ti);
			      Tt = LD(&(x[WS(ios, 11)]), dist, &(x[WS(ios, 1)]));
			      To = LD(&(x[WS(ios, 7)]), dist, &(x[WS(ios, 1)]));
			      TZ = VADD(TX, TY);
			      T1a = VSUB(TY, TX);
			      Te = BYTWJ(&(W[TWVL * 16]), Td);
			      Th = BYTWJ(&(W[TWVL * 8]), Tg);
			      Tn = BYTWJ(&(W[TWVL * 28]), Tm);
			      Tr = LD(&(x[WS(ios, 3)]), dist, &(x[WS(ios, 1)]));
			      Tu = BYTWJ(&(W[TWVL * 20]), Tt);
			      Tp = BYTWJ(&(W[TWVL * 12]), To);
			 }
		    }
	       }
	       {
		    V Tf, T11, Tk, T12, Ts;
		    TO = VADD(T4, T9);
		    Ta = VSUB(T4, T9);
		    TJ = VSUB(TD, TI);
		    TP = VADD(TI, TD);
		    Tf = VADD(Tc, Te);
		    T11 = VSUB(Tc, Te);
		    Tk = VADD(Th, Tj);
		    T12 = VSUB(Th, Tj);
		    Ts = BYTWJ(&(W[TWVL * 4]), Tr);
		    T14 = VSUB(Tn, Tp);
		    Tq = VADD(Tn, Tp);
		    T1i = VFNMS(LDK(KP707106781), TZ, TW);
		    T10 = VFMA(LDK(KP707106781), TZ, TW);
		    T1b = VFNMS(LDK(KP707106781), T1a, T19);
		    T1l = VFMA(LDK(KP707106781), T1a, T19);
		    T13 = VFNMS(LDK(KP414213562), T12, T11);
		    T1c = VFMA(LDK(KP414213562), T11, T12);
		    TR = VADD(Tf, Tk);
		    Tl = VSUB(Tf, Tk);
		    T15 = VSUB(Tu, Ts);
		    Tv = VADD(Ts, Tu);
	       }
	  }
	  {
	       V T1d, T16, TS, Tw, TU, TQ;
	       T1d = VFMA(LDK(KP414213562), T14, T15);
	       T16 = VFNMS(LDK(KP414213562), T15, T14);
	       TS = VADD(Tq, Tv);
	       Tw = VSUB(Tq, Tv);
	       TU = VSUB(TO, TP);
	       TQ = VADD(TO, TP);
	       {
		    V T1e, T1j, T17, T1m;
		    T1e = VSUB(T1c, T1d);
		    T1j = VADD(T1c, T1d);
		    T17 = VADD(T13, T16);
		    T1m = VSUB(T16, T13);
		    {
			 V TV, TT, TK, Tx;
			 TV = VSUB(TS, TR);
			 TT = VADD(TR, TS);
			 TK = VSUB(Tw, Tl);
			 Tx = VADD(Tl, Tw);
			 {
			      V T1h, T1f, T1o, T1k;
			      T1h = VFMA(LDK(KP923879532), T1e, T1b);
			      T1f = VFNMS(LDK(KP923879532), T1e, T1b);
			      T1o = VFMA(LDK(KP923879532), T1j, T1i);
			      T1k = VFNMS(LDK(KP923879532), T1j, T1i);
			      {
				   V T1g, T18, T1p, T1n;
				   T1g = VFMA(LDK(KP923879532), T17, T10);
				   T18 = VFNMS(LDK(KP923879532), T17, T10);
				   T1p = VFMA(LDK(KP923879532), T1m, T1l);
				   T1n = VFNMS(LDK(KP923879532), T1m, T1l);
				   ST(&(x[WS(ios, 12)]), VFNMSI(TV, TU), dist, &(x[0]));
				   ST(&(x[WS(ios, 4)]), VFMAI(TV, TU), dist, &(x[0]));
				   ST(&(x[0]), VADD(TQ, TT), dist, &(x[0]));
				   ST(&(x[WS(ios, 8)]), VSUB(TQ, TT), dist, &(x[0]));
				   {
					V TN, TL, TM, Ty;
					TN = VFMA(LDK(KP707106781), TK, TJ);
					TL = VFNMS(LDK(KP707106781), TK, TJ);
					TM = VFMA(LDK(KP707106781), Tx, Ta);
					Ty = VFNMS(LDK(KP707106781), Tx, Ta);
					ST(&(x[WS(ios, 1)]), VFNMSI(T1h, T1g), dist, &(x[WS(ios, 1)]));
					ST(&(x[WS(ios, 15)]), VFMAI(T1h, T1g), dist, &(x[WS(ios, 1)]));
					ST(&(x[WS(ios, 7)]), VFMAI(T1f, T18), dist, &(x[WS(ios, 1)]));
					ST(&(x[WS(ios, 9)]), VFNMSI(T1f, T18), dist, &(x[WS(ios, 1)]));
					ST(&(x[WS(ios, 3)]), VFMAI(T1p, T1o), dist, &(x[WS(ios, 1)]));
					ST(&(x[WS(ios, 13)]), VFNMSI(T1p, T1o), dist, &(x[WS(ios, 1)]));
					ST(&(x[WS(ios, 11)]), VFMAI(T1n, T1k), dist, &(x[WS(ios, 1)]));
					ST(&(x[WS(ios, 5)]), VFNMSI(T1n, T1k), dist, &(x[WS(ios, 1)]));
					ST(&(x[WS(ios, 14)]), VFNMSI(TN, TM), dist, &(x[0]));
					ST(&(x[WS(ios, 2)]), VFMAI(TN, TM), dist, &(x[0]));
					ST(&(x[WS(ios, 10)]), VFMAI(TL, Ty), dist, &(x[0]));
					ST(&(x[WS(ios, 6)]), VFNMSI(TL, Ty), dist, &(x[0]));
				   }
			      }
			 }
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
     VTW(5),
     VTW(6),
     VTW(7),
     VTW(8),
     VTW(9),
     VTW(10),
     VTW(11),
     VTW(12),
     VTW(13),
     VTW(14),
     VTW(15),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 16, "t2fv_16", twinstr, &GENUS, {53, 30, 34, 0}, 0, 0, 0 };

void X(codelet_t2fv_16) (planner *p) {
     X(kdft_dit_register) (p, t2fv_16, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle_c -simd -compact -variables 4 -pipeline-latency 8 -n 16 -name t2fv_16 -include t2f.h */

/*
 * This function contains 87 FP additions, 42 FP multiplications,
 * (or, 83 additions, 38 multiplications, 4 fused multiply/add),
 * 36 stack variables, and 32 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle_c.ml,v 1.13 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t2f.h"

static const R *t2fv_16(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     R *x;
     x = ri;
     for (i = m; i > 0; i = i - VL, x = x + (VL * dist), W = W + (TWVL * 30), MAKE_VOLATILE_STRIDE(ios)) {
	  V TJ, T10, TD, T11, T1b, T1c, Ty, TK, T16, T17, T18, Tb, TN, T13, T14;
	  V T15, Tm, TM, TG, TI, TH;
	  TG = LD(&(x[0]), dist, &(x[0]));
	  TH = LD(&(x[WS(ios, 8)]), dist, &(x[0]));
	  TI = BYTWJ(&(W[TWVL * 14]), TH);
	  TJ = VSUB(TG, TI);
	  T10 = VADD(TG, TI);
	  {
	       V TA, TC, Tz, TB;
	       Tz = LD(&(x[WS(ios, 4)]), dist, &(x[0]));
	       TA = BYTWJ(&(W[TWVL * 6]), Tz);
	       TB = LD(&(x[WS(ios, 12)]), dist, &(x[0]));
	       TC = BYTWJ(&(W[TWVL * 22]), TB);
	       TD = VSUB(TA, TC);
	       T11 = VADD(TA, TC);
	  }
	  {
	       V Tp, Tw, Tr, Tu, Ts, Tx;
	       {
		    V To, Tv, Tq, Tt;
		    To = LD(&(x[WS(ios, 14)]), dist, &(x[0]));
		    Tp = BYTWJ(&(W[TWVL * 26]), To);
		    Tv = LD(&(x[WS(ios, 10)]), dist, &(x[0]));
		    Tw = BYTWJ(&(W[TWVL * 18]), Tv);
		    Tq = LD(&(x[WS(ios, 6)]), dist, &(x[0]));
		    Tr = BYTWJ(&(W[TWVL * 10]), Tq);
		    Tt = LD(&(x[WS(ios, 2)]), dist, &(x[0]));
		    Tu = BYTWJ(&(W[TWVL * 2]), Tt);
	       }
	       T1b = VADD(Tp, Tr);
	       T1c = VADD(Tu, Tw);
	       Ts = VSUB(Tp, Tr);
	       Tx = VSUB(Tu, Tw);
	       Ty = VMUL(LDK(KP707106781), VSUB(Ts, Tx));
	       TK = VMUL(LDK(KP707106781), VADD(Tx, Ts));
	  }
	  {
	       V T2, T9, T4, T7, T5, Ta;
	       {
		    V T1, T8, T3, T6;
		    T1 = LD(&(x[WS(ios, 15)]), dist, &(x[WS(ios, 1)]));
		    T2 = BYTWJ(&(W[TWVL * 28]), T1);
		    T8 = LD(&(x[WS(ios, 11)]), dist, &(x[WS(ios, 1)]));
		    T9 = BYTWJ(&(W[TWVL * 20]), T8);
		    T3 = LD(&(x[WS(ios, 7)]), dist, &(x[WS(ios, 1)]));
		    T4 = BYTWJ(&(W[TWVL * 12]), T3);
		    T6 = LD(&(x[WS(ios, 3)]), dist, &(x[WS(ios, 1)]));
		    T7 = BYTWJ(&(W[TWVL * 4]), T6);
	       }
	       T16 = VADD(T2, T4);
	       T17 = VADD(T7, T9);
	       T18 = VSUB(T16, T17);
	       T5 = VSUB(T2, T4);
	       Ta = VSUB(T7, T9);
	       Tb = VFNMS(LDK(KP923879532), Ta, VMUL(LDK(KP382683432), T5));
	       TN = VFMA(LDK(KP923879532), T5, VMUL(LDK(KP382683432), Ta));
	  }
	  {
	       V Td, Tk, Tf, Ti, Tg, Tl;
	       {
		    V Tc, Tj, Te, Th;
		    Tc = LD(&(x[WS(ios, 1)]), dist, &(x[WS(ios, 1)]));
		    Td = BYTWJ(&(W[0]), Tc);
		    Tj = LD(&(x[WS(ios, 13)]), dist, &(x[WS(ios, 1)]));
		    Tk = BYTWJ(&(W[TWVL * 24]), Tj);
		    Te = LD(&(x[WS(ios, 9)]), dist, &(x[WS(ios, 1)]));
		    Tf = BYTWJ(&(W[TWVL * 16]), Te);
		    Th = LD(&(x[WS(ios, 5)]), dist, &(x[WS(ios, 1)]));
		    Ti = BYTWJ(&(W[TWVL * 8]), Th);
	       }
	       T13 = VADD(Td, Tf);
	       T14 = VADD(Ti, Tk);
	       T15 = VSUB(T13, T14);
	       Tg = VSUB(Td, Tf);
	       Tl = VSUB(Ti, Tk);
	       Tm = VFMA(LDK(KP382683432), Tg, VMUL(LDK(KP923879532), Tl));
	       TM = VFNMS(LDK(KP382683432), Tl, VMUL(LDK(KP923879532), Tg));
	  }
	  {
	       V T1a, T1g, T1f, T1h;
	       {
		    V T12, T19, T1d, T1e;
		    T12 = VSUB(T10, T11);
		    T19 = VMUL(LDK(KP707106781), VADD(T15, T18));
		    T1a = VADD(T12, T19);
		    T1g = VSUB(T12, T19);
		    T1d = VSUB(T1b, T1c);
		    T1e = VMUL(LDK(KP707106781), VSUB(T18, T15));
		    T1f = VBYI(VADD(T1d, T1e));
		    T1h = VBYI(VSUB(T1e, T1d));
	       }
	       ST(&(x[WS(ios, 14)]), VSUB(T1a, T1f), dist, &(x[0]));
	       ST(&(x[WS(ios, 6)]), VADD(T1g, T1h), dist, &(x[0]));
	       ST(&(x[WS(ios, 2)]), VADD(T1a, T1f), dist, &(x[0]));
	       ST(&(x[WS(ios, 10)]), VSUB(T1g, T1h), dist, &(x[0]));
	  }
	  {
	       V T1k, T1o, T1n, T1p;
	       {
		    V T1i, T1j, T1l, T1m;
		    T1i = VADD(T10, T11);
		    T1j = VADD(T1c, T1b);
		    T1k = VADD(T1i, T1j);
		    T1o = VSUB(T1i, T1j);
		    T1l = VADD(T13, T14);
		    T1m = VADD(T16, T17);
		    T1n = VADD(T1l, T1m);
		    T1p = VBYI(VSUB(T1m, T1l));
	       }
	       ST(&(x[WS(ios, 8)]), VSUB(T1k, T1n), dist, &(x[0]));
	       ST(&(x[WS(ios, 4)]), VADD(T1o, T1p), dist, &(x[0]));
	       ST(&(x[0]), VADD(T1k, T1n), dist, &(x[0]));
	       ST(&(x[WS(ios, 12)]), VSUB(T1o, T1p), dist, &(x[0]));
	  }
	  {
	       V TF, TQ, TP, TR;
	       {
		    V Tn, TE, TL, TO;
		    Tn = VSUB(Tb, Tm);
		    TE = VSUB(Ty, TD);
		    TF = VBYI(VSUB(Tn, TE));
		    TQ = VBYI(VADD(TE, Tn));
		    TL = VADD(TJ, TK);
		    TO = VADD(TM, TN);
		    TP = VSUB(TL, TO);
		    TR = VADD(TL, TO);
	       }
	       ST(&(x[WS(ios, 7)]), VADD(TF, TP), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 15)]), VSUB(TR, TQ), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 9)]), VSUB(TP, TF), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 1)]), VADD(TQ, TR), dist, &(x[WS(ios, 1)]));
	  }
	  {
	       V TU, TY, TX, TZ;
	       {
		    V TS, TT, TV, TW;
		    TS = VSUB(TJ, TK);
		    TT = VADD(Tm, Tb);
		    TU = VADD(TS, TT);
		    TY = VSUB(TS, TT);
		    TV = VADD(TD, Ty);
		    TW = VSUB(TN, TM);
		    TX = VBYI(VADD(TV, TW));
		    TZ = VBYI(VSUB(TW, TV));
	       }
	       ST(&(x[WS(ios, 13)]), VSUB(TU, TX), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 5)]), VADD(TY, TZ), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 3)]), VADD(TU, TX), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 11)]), VSUB(TY, TZ), dist, &(x[WS(ios, 1)]));
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     VTW(1),
     VTW(2),
     VTW(3),
     VTW(4),
     VTW(5),
     VTW(6),
     VTW(7),
     VTW(8),
     VTW(9),
     VTW(10),
     VTW(11),
     VTW(12),
     VTW(13),
     VTW(14),
     VTW(15),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 16, "t2fv_16", twinstr, &GENUS, {83, 38, 4, 0}, 0, 0, 0 };

void X(codelet_t2fv_16) (planner *p) {
     X(kdft_dit_register) (p, t2fv_16, &desc);
}
#endif				/* HAVE_FMA */
