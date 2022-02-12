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
/* Generated on Fri Jan 27 20:52:53 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2r -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 32 -name hc2rIII_32 -dft-III -include hc2rIII.h */

/*
 * This function contains 174 FP additions, 100 FP multiplications,
 * (or, 106 additions, 32 multiplications, 68 fused multiply/add),
 * 101 stack variables, and 64 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2rIII.h"

static void hc2rIII_32(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP534511135, +0.534511135950791641089685961295362908582039528);
     DK(KP1_763842528, +1.763842528696710059425513727320776699016885241);
     DK(KP303346683, +0.303346683607342391675883946941299872384187453);
     DK(KP1_913880671, +1.913880671464417729871595773960539938965698411);
     DK(KP098491403, +0.098491403357164253077197521291327432293052451);
     DK(KP1_990369453, +1.990369453344393772489673906218959843150949737);
     DK(KP820678790, +0.820678790828660330972281985331011598767386482);
     DK(KP1_546020906, +1.546020906725473921621813219516939601942082586);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP668178637, +0.668178637919298919997757686523080761552472251);
     DK(KP1_662939224, +1.662939224605090474157576755235811513477121624);
     DK(KP198912367, +0.198912367379658006911597622644676228597850501);
     DK(KP1_961570560, +1.961570560806460898252364472268478073947867462);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E T1N, T1K, T1Q, T1H, T1O, T1P;
	  {
	       E T1I, T1e, T1Z, T7, T2E, T2i, T1x, Tz, Te, T2j, T22, T2F, T1h, T1y, TK;
	       E T1J, Tm, T2B, TX, Tp, T2m, T28, T1M, T1C, T1k, TW, TY, T2a, T14, T15;
	       E Ts, TZ;
	       {
		    E TE, T1g, TJ, T1f;
		    {
			 E T4, Tv, T3, T2g, T1d, T5, Tw, Tx;
			 {
			      E T1, T2, T1b, T1c;
			      T1 = ri[0];
			      T2 = ri[WS(ris, 15)];
			      T1b = ii[0];
			      T1c = ii[WS(iis, 15)];
			      T4 = ri[WS(ris, 8)];
			      Tv = T1 - T2;
			      T3 = T1 + T2;
			      T2g = T1c - T1b;
			      T1d = T1b + T1c;
			      T5 = ri[WS(ris, 7)];
			      Tw = ii[WS(iis, 8)];
			      Tx = ii[WS(iis, 7)];
			 }
			 {
			      E Tb, TA, Ta, T20, TD, Tc, TG, TH;
			      {
				   E T8, T9, TB, TC;
				   T8 = ri[WS(ris, 4)];
				   {
					E T1a, T6, T2h, Ty;
					T1a = T4 - T5;
					T6 = T4 + T5;
					T2h = Tx - Tw;
					Ty = Tw + Tx;
					T1I = T1a - T1d;
					T1e = T1a + T1d;
					T1Z = T3 - T6;
					T7 = T3 + T6;
					T2E = T2h + T2g;
					T2i = T2g - T2h;
					T1x = Tv + Ty;
					Tz = Tv - Ty;
					T9 = ri[WS(ris, 11)];
				   }
				   TB = ii[WS(iis, 4)];
				   TC = ii[WS(iis, 11)];
				   Tb = ri[WS(ris, 3)];
				   TA = T8 - T9;
				   Ta = T8 + T9;
				   T20 = TC - TB;
				   TD = TB + TC;
				   Tc = ri[WS(ris, 12)];
				   TG = ii[WS(iis, 3)];
				   TH = ii[WS(iis, 12)];
			      }
			      {
				   E TF, Td, T21, TI;
				   TE = TA - TD;
				   T1g = TA + TD;
				   TF = Tb - Tc;
				   Td = Tb + Tc;
				   T21 = TG - TH;
				   TI = TG + TH;
				   Te = Ta + Td;
				   T2j = Ta - Td;
				   T22 = T20 - T21;
				   T2F = T20 + T21;
				   TJ = TF - TI;
				   T1f = TF + TI;
			      }
			 }
		    }
		    {
			 E TM, Ti, TN, T25, TU, TR, Tl, TO;
			 {
			      E TS, TT, Tg, Th, Tj, Tk;
			      Tg = ri[WS(ris, 2)];
			      Th = ri[WS(ris, 13)];
			      T1h = T1f - T1g;
			      T1y = T1g + T1f;
			      TK = TE + TJ;
			      T1J = TE - TJ;
			      TM = Tg - Th;
			      Ti = Tg + Th;
			      TS = ii[WS(iis, 2)];
			      TT = ii[WS(iis, 13)];
			      Tj = ri[WS(ris, 10)];
			      Tk = ri[WS(ris, 5)];
			      TN = ii[WS(iis, 10)];
			      T25 = TS - TT;
			      TU = TS + TT;
			      TR = Tj - Tk;
			      Tl = Tj + Tk;
			      TO = ii[WS(iis, 5)];
			 }
			 {
			      E T12, T13, Tq, Tr;
			      {
				   E Tn, T1A, TV, T24, T26, TP, To, T27, T1B, TQ;
				   Tn = ri[WS(ris, 1)];
				   T1A = TR - TU;
				   TV = TR + TU;
				   T24 = Ti - Tl;
				   Tm = Ti + Tl;
				   T26 = TN - TO;
				   TP = TN + TO;
				   To = ri[WS(ris, 14)];
				   T12 = ii[WS(iis, 1)];
				   T27 = T25 - T26;
				   T2B = T26 + T25;
				   T1B = TM + TP;
				   TQ = TM - TP;
				   TX = Tn - To;
				   Tp = Tn + To;
				   T2m = T24 + T27;
				   T28 = T24 - T27;
				   T1M = FNMS(KP414213562, T1A, T1B);
				   T1C = FMA(KP414213562, T1B, T1A);
				   T1k = FMA(KP414213562, TQ, TV);
				   TW = FNMS(KP414213562, TV, TQ);
				   T13 = ii[WS(iis, 14)];
			      }
			      Tq = ri[WS(ris, 6)];
			      Tr = ri[WS(ris, 9)];
			      TY = ii[WS(iis, 6)];
			      T2a = T13 - T12;
			      T14 = T12 + T13;
			      T15 = Tq - Tr;
			      Ts = Tq + Tr;
			      TZ = ii[WS(iis, 9)];
			 }
		    }
	       }
	       {
		    E T1L, T1F, T23, T2n, T2k, T2e, T1p, T1t, T1s, T1i, T1o, T19, T1l, T1q;
		    {
			 E T2z, T2G, T2H, T2C, T1j, T17, T2r, T2s, T2u, T2v, T2K, T2D;
			 {
			      E T2L, T2d, T2l, T2O;
			      {
				   E Tf, T2N, Tu, T2M;
				   {
					E T1D, T16, T29, Tt, T2b, T10;
					T2z = T7 - Te;
					Tf = T7 + Te;
					T1D = T15 + T14;
					T16 = T14 - T15;
					T29 = Tp - Ts;
					Tt = Tp + Ts;
					T2b = TY - TZ;
					T10 = TY + TZ;
					T2N = T2F + T2E;
					T2G = T2E - T2F;
					T2H = Tm - Tt;
					Tu = Tm + Tt;
					{
					     E T2c, T2A, T1E, T11;
					     T2c = T2a - T2b;
					     T2A = T2b + T2a;
					     T1E = TX + T10;
					     T11 = TX - T10;
					     T2L = Tf - Tu;
					     T2d = T29 + T2c;
					     T2l = T29 - T2c;
					     T2C = T2A - T2B;
					     T2M = T2B + T2A;
					     T1L = FMA(KP414213562, T1D, T1E);
					     T1F = FNMS(KP414213562, T1E, T1D);
					     T1j = FMA(KP414213562, T11, T16);
					     T17 = FNMS(KP414213562, T16, T11);
					     T2O = T2M + T2N;
					}
				   }
				   O[0] = KP2_000000000 * (Tf + Tu);
				   O[WS(os, 16)] = KP2_000000000 * (T2N - T2M);
			      }
			      T23 = T1Z + T22;
			      T2r = T1Z - T22;
			      O[WS(os, 24)] = KP1_414213562 * (T2O - T2L);
			      O[WS(os, 8)] = KP1_414213562 * (T2L + T2O);
			      T2s = T2m + T2l;
			      T2n = T2l - T2m;
			      T2k = T2i - T2j;
			      T2u = T2j + T2i;
			      T2v = T28 - T2d;
			      T2e = T28 + T2d;
			 }
			 {
			      E T2y, T2t, T2x, T2w;
			      T2y = FMA(KP707106781, T2s, T2r);
			      T2t = FNMS(KP707106781, T2s, T2r);
			      T2x = FMA(KP707106781, T2v, T2u);
			      T2w = FNMS(KP707106781, T2v, T2u);
			      O[WS(os, 14)] = KP1_961570560 * (FMA(KP198912367, T2y, T2x));
			      O[WS(os, 30)] = -(KP1_961570560 * (FNMS(KP198912367, T2x, T2y)));
			      O[WS(os, 22)] = KP1_662939224 * (FNMS(KP668178637, T2t, T2w));
			      O[WS(os, 6)] = KP1_662939224 * (FMA(KP668178637, T2w, T2t));
			      T2K = T2z - T2C;
			      T2D = T2z + T2C;
			 }
			 {
			      E TL, T18, T2J, T2I;
			      T1p = FNMS(KP707106781, TK, Tz);
			      TL = FMA(KP707106781, TK, Tz);
			      T18 = TW + T17;
			      T1t = TW - T17;
			      T1s = FMA(KP707106781, T1h, T1e);
			      T1i = FNMS(KP707106781, T1h, T1e);
			      T2J = T2H + T2G;
			      T2I = T2G - T2H;
			      T1o = FNMS(KP923879532, T18, TL);
			      T19 = FMA(KP923879532, T18, TL);
			      O[WS(os, 12)] = KP1_847759065 * (FMA(KP414213562, T2K, T2J));
			      O[WS(os, 28)] = -(KP1_847759065 * (FNMS(KP414213562, T2J, T2K)));
			      O[WS(os, 20)] = KP1_847759065 * (FNMS(KP414213562, T2D, T2I));
			      O[WS(os, 4)] = KP1_847759065 * (FMA(KP414213562, T2I, T2D));
			      T1l = T1j - T1k;
			      T1q = T1k + T1j;
			 }
		    }
		    {
			 E T1z, T1U, T1Y, T1T, T1V, T1G;
			 {
			      E T1w, T1r, T1n, T1m;
			      T1n = FMA(KP923879532, T1l, T1i);
			      T1m = FNMS(KP923879532, T1l, T1i);
			      T1w = FMA(KP923879532, T1q, T1p);
			      T1r = FNMS(KP923879532, T1q, T1p);
			      O[WS(os, 9)] = -(KP1_546020906 * (FNMS(KP820678790, T1o, T1n)));
			      O[WS(os, 25)] = -(KP1_546020906 * (FMA(KP820678790, T1n, T1o)));
			      O[WS(os, 17)] = -(KP1_990369453 * (FMA(KP098491403, T19, T1m)));
			      O[WS(os, 1)] = KP1_990369453 * (FNMS(KP098491403, T1m, T19));
			      {
				   E T1R, T1S, T1v, T1u;
				   T1z = FNMS(KP707106781, T1y, T1x);
				   T1R = FMA(KP707106781, T1y, T1x);
				   T1S = T1M + T1L;
				   T1N = T1L - T1M;
				   T1K = FNMS(KP707106781, T1J, T1I);
				   T1U = FMA(KP707106781, T1J, T1I);
				   T1v = FNMS(KP923879532, T1t, T1s);
				   T1u = FMA(KP923879532, T1t, T1s);
				   T1Y = FMA(KP923879532, T1S, T1R);
				   T1T = FNMS(KP923879532, T1S, T1R);
				   O[WS(os, 13)] = -(KP1_913880671 * (FNMS(KP303346683, T1w, T1v)));
				   O[WS(os, 29)] = -(KP1_913880671 * (FMA(KP303346683, T1v, T1w)));
				   O[WS(os, 21)] = -(KP1_763842528 * (FMA(KP534511135, T1r, T1u)));
				   O[WS(os, 5)] = KP1_763842528 * (FNMS(KP534511135, T1u, T1r));
				   T1V = T1C + T1F;
				   T1G = T1C - T1F;
			      }
			 }
			 {
			      E T2q, T2f, T1X, T1W, T2p, T2o;
			      T1X = FMA(KP923879532, T1V, T1U);
			      T1W = FNMS(KP923879532, T1V, T1U);
			      T2q = FNMS(KP707106781, T2e, T23);
			      T2f = FMA(KP707106781, T2e, T23);
			      O[WS(os, 15)] = KP1_990369453 * (FMA(KP098491403, T1Y, T1X));
			      O[WS(os, 31)] = -(KP1_990369453 * (FNMS(KP098491403, T1X, T1Y)));
			      O[WS(os, 23)] = KP1_546020906 * (FNMS(KP820678790, T1T, T1W));
			      O[WS(os, 7)] = KP1_546020906 * (FMA(KP820678790, T1W, T1T));
			      T2p = FNMS(KP707106781, T2n, T2k);
			      T2o = FMA(KP707106781, T2n, T2k);
			      T1Q = FNMS(KP923879532, T1G, T1z);
			      T1H = FMA(KP923879532, T1G, T1z);
			      O[WS(os, 10)] = KP1_662939224 * (FMA(KP668178637, T2q, T2p));
			      O[WS(os, 26)] = -(KP1_662939224 * (FNMS(KP668178637, T2p, T2q)));
			      O[WS(os, 18)] = KP1_961570560 * (FNMS(KP198912367, T2f, T2o));
			      O[WS(os, 2)] = KP1_961570560 * (FMA(KP198912367, T2o, T2f));
			 }
		    }
	       }
	  }
	  T1O = FMA(KP923879532, T1N, T1K);
	  T1P = FNMS(KP923879532, T1N, T1K);
	  O[WS(os, 11)] = KP1_763842528 * (FMA(KP534511135, T1Q, T1P));
	  O[WS(os, 27)] = -(KP1_763842528 * (FNMS(KP534511135, T1P, T1Q)));
	  O[WS(os, 19)] = KP1_913880671 * (FNMS(KP303346683, T1H, T1O));
	  O[WS(os, 3)] = KP1_913880671 * (FMA(KP303346683, T1O, T1H));
     }
}

static const khc2r_desc desc = { 32, "hc2rIII_32", {106, 32, 68, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2rIII_32) (planner *p) {
     X(khc2rIII_register) (p, hc2rIII_32, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2r -compact -variables 4 -pipeline-latency 4 -sign 1 -n 32 -name hc2rIII_32 -dft-III -include hc2rIII.h */

/*
 * This function contains 174 FP additions, 84 FP multiplications,
 * (or, 138 additions, 48 multiplications, 36 fused multiply/add),
 * 66 stack variables, and 64 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2rIII.h"

static void hc2rIII_32(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP1_913880671, +1.913880671464417729871595773960539938965698411);
     DK(KP580569354, +0.580569354508924735272384751634790549382952557);
     DK(KP942793473, +0.942793473651995297112775251810508755314920638);
     DK(KP1_763842528, +1.763842528696710059425513727320776699016885241);
     DK(KP1_546020906, +1.546020906725473921621813219516939601942082586);
     DK(KP1_268786568, +1.268786568327290996430343226450986741351374190);
     DK(KP196034280, +0.196034280659121203988391127777283691722273346);
     DK(KP1_990369453, +1.990369453344393772489673906218959843150949737);
     DK(KP765366864, +0.765366864730179543456919968060797733522689125);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP1_961570560, +1.961570560806460898252364472268478073947867462);
     DK(KP390180644, +0.390180644032256535696569736954044481855383236);
     DK(KP1_111140466, +1.111140466039204449485661627897065748749874382);
     DK(KP1_662939224, +1.662939224605090474157576755235811513477121624);
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     DK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E T7, T2i, T2F, Tz, T1k, T1I, T1Z, T1x, Te, T22, T2E, T2j, T1f, T1y, TK;
	  E T1J, Tm, T2B, TW, T1a, T1C, T1L, T28, T2l, Tt, T2A, T17, T1b, T1F, T1M;
	  E T2d, T2m;
	  {
	       E T3, Tv, T1j, T2h, T6, T1g, Ty, T2g;
	       {
		    E T1, T2, T1h, T1i;
		    T1 = ri[0];
		    T2 = ri[WS(ris, 15)];
		    T3 = T1 + T2;
		    Tv = T1 - T2;
		    T1h = ii[0];
		    T1i = ii[WS(iis, 15)];
		    T1j = T1h + T1i;
		    T2h = T1i - T1h;
	       }
	       {
		    E T4, T5, Tw, Tx;
		    T4 = ri[WS(ris, 8)];
		    T5 = ri[WS(ris, 7)];
		    T6 = T4 + T5;
		    T1g = T4 - T5;
		    Tw = ii[WS(iis, 8)];
		    Tx = ii[WS(iis, 7)];
		    Ty = Tw + Tx;
		    T2g = Tw - Tx;
	       }
	       T7 = T3 + T6;
	       T2i = T2g + T2h;
	       T2F = T2h - T2g;
	       Tz = Tv - Ty;
	       T1k = T1g + T1j;
	       T1I = T1g - T1j;
	       T1Z = T3 - T6;
	       T1x = Tv + Ty;
	  }
	  {
	       E Ta, TA, TD, T21, Td, TF, TI, T20;
	       {
		    E T8, T9, TB, TC;
		    T8 = ri[WS(ris, 4)];
		    T9 = ri[WS(ris, 11)];
		    Ta = T8 + T9;
		    TA = T8 - T9;
		    TB = ii[WS(iis, 4)];
		    TC = ii[WS(iis, 11)];
		    TD = TB + TC;
		    T21 = TB - TC;
	       }
	       {
		    E Tb, Tc, TG, TH;
		    Tb = ri[WS(ris, 3)];
		    Tc = ri[WS(ris, 12)];
		    Td = Tb + Tc;
		    TF = Tb - Tc;
		    TG = ii[WS(iis, 3)];
		    TH = ii[WS(iis, 12)];
		    TI = TG + TH;
		    T20 = TH - TG;
	       }
	       Te = Ta + Td;
	       T22 = T20 - T21;
	       T2E = T21 + T20;
	       T2j = Ta - Td;
	       {
		    E T1d, T1e, TE, TJ;
		    T1d = TA + TD;
		    T1e = TF + TI;
		    T1f = KP707106781 * (T1d - T1e);
		    T1y = KP707106781 * (T1d + T1e);
		    TE = TA - TD;
		    TJ = TF - TI;
		    TK = KP707106781 * (TE + TJ);
		    T1J = KP707106781 * (TE - TJ);
	       }
	  }
	  {
	       E Ti, TM, TU, T25, Tl, TR, TP, T26, TQ, TV;
	       {
		    E Tg, Th, TS, TT;
		    Tg = ri[WS(ris, 2)];
		    Th = ri[WS(ris, 13)];
		    Ti = Tg + Th;
		    TM = Tg - Th;
		    TS = ii[WS(iis, 2)];
		    TT = ii[WS(iis, 13)];
		    TU = TS + TT;
		    T25 = TS - TT;
	       }
	       {
		    E Tj, Tk, TN, TO;
		    Tj = ri[WS(ris, 10)];
		    Tk = ri[WS(ris, 5)];
		    Tl = Tj + Tk;
		    TR = Tj - Tk;
		    TN = ii[WS(iis, 10)];
		    TO = ii[WS(iis, 5)];
		    TP = TN + TO;
		    T26 = TN - TO;
	       }
	       Tm = Ti + Tl;
	       T2B = T26 + T25;
	       TQ = TM - TP;
	       TV = TR + TU;
	       TW = FNMS(KP382683432, TV, KP923879532 * TQ);
	       T1a = FMA(KP382683432, TQ, KP923879532 * TV);
	       {
		    E T1A, T1B, T24, T27;
		    T1A = TM + TP;
		    T1B = TU - TR;
		    T1C = FNMS(KP923879532, T1B, KP382683432 * T1A);
		    T1L = FMA(KP923879532, T1A, KP382683432 * T1B);
		    T24 = Ti - Tl;
		    T27 = T25 - T26;
		    T28 = T24 - T27;
		    T2l = T24 + T27;
	       }
	  }
	  {
	       E Tp, TX, T15, T2a, Ts, T12, T10, T2b, T11, T16;
	       {
		    E Tn, To, T13, T14;
		    Tn = ri[WS(ris, 1)];
		    To = ri[WS(ris, 14)];
		    Tp = Tn + To;
		    TX = Tn - To;
		    T13 = ii[WS(iis, 1)];
		    T14 = ii[WS(iis, 14)];
		    T15 = T13 + T14;
		    T2a = T14 - T13;
	       }
	       {
		    E Tq, Tr, TY, TZ;
		    Tq = ri[WS(ris, 6)];
		    Tr = ri[WS(ris, 9)];
		    Ts = Tq + Tr;
		    T12 = Tq - Tr;
		    TY = ii[WS(iis, 6)];
		    TZ = ii[WS(iis, 9)];
		    T10 = TY + TZ;
		    T2b = TY - TZ;
	       }
	       Tt = Tp + Ts;
	       T2A = T2b + T2a;
	       T11 = TX - T10;
	       T16 = T12 - T15;
	       T17 = FMA(KP923879532, T11, KP382683432 * T16);
	       T1b = FNMS(KP382683432, T11, KP923879532 * T16);
	       {
		    E T1D, T1E, T29, T2c;
		    T1D = TX + T10;
		    T1E = T12 + T15;
		    T1F = FNMS(KP923879532, T1E, KP382683432 * T1D);
		    T1M = FMA(KP923879532, T1D, KP382683432 * T1E);
		    T29 = Tp - Ts;
		    T2c = T2a - T2b;
		    T2d = T29 + T2c;
		    T2m = T2c - T29;
	       }
	  }
	  {
	       E Tf, Tu, T2L, T2M, T2N, T2O;
	       Tf = T7 + Te;
	       Tu = Tm + Tt;
	       T2L = Tf - Tu;
	       T2M = T2B + T2A;
	       T2N = T2F - T2E;
	       T2O = T2M + T2N;
	       O[0] = KP2_000000000 * (Tf + Tu);
	       O[WS(os, 16)] = KP2_000000000 * (T2N - T2M);
	       O[WS(os, 8)] = KP1_414213562 * (T2L + T2O);
	       O[WS(os, 24)] = KP1_414213562 * (T2O - T2L);
	  }
	  {
	       E T2t, T2x, T2w, T2y;
	       {
		    E T2r, T2s, T2u, T2v;
		    T2r = T1Z - T22;
		    T2s = KP707106781 * (T2m - T2l);
		    T2t = T2r + T2s;
		    T2x = T2r - T2s;
		    T2u = T2j + T2i;
		    T2v = KP707106781 * (T28 - T2d);
		    T2w = T2u - T2v;
		    T2y = T2v + T2u;
	       }
	       O[WS(os, 6)] = FMA(KP1_662939224, T2t, KP1_111140466 * T2w);
	       O[WS(os, 30)] = FNMS(KP1_961570560, T2x, KP390180644 * T2y);
	       O[WS(os, 22)] = FNMS(KP1_111140466, T2t, KP1_662939224 * T2w);
	       O[WS(os, 14)] = FMA(KP390180644, T2x, KP1_961570560 * T2y);
	  }
	  {
	       E T2D, T2J, T2I, T2K;
	       {
		    E T2z, T2C, T2G, T2H;
		    T2z = T7 - Te;
		    T2C = T2A - T2B;
		    T2D = T2z + T2C;
		    T2J = T2z - T2C;
		    T2G = T2E + T2F;
		    T2H = Tm - Tt;
		    T2I = T2G - T2H;
		    T2K = T2H + T2G;
	       }
	       O[WS(os, 4)] = FMA(KP1_847759065, T2D, KP765366864 * T2I);
	       O[WS(os, 28)] = FNMS(KP1_847759065, T2J, KP765366864 * T2K);
	       O[WS(os, 20)] = FNMS(KP765366864, T2D, KP1_847759065 * T2I);
	       O[WS(os, 12)] = FMA(KP765366864, T2J, KP1_847759065 * T2K);
	  }
	  {
	       E T19, T1n, T1m, T1o;
	       {
		    E TL, T18, T1c, T1l;
		    TL = Tz + TK;
		    T18 = TW + T17;
		    T19 = TL + T18;
		    T1n = TL - T18;
		    T1c = T1a + T1b;
		    T1l = T1f + T1k;
		    T1m = T1c + T1l;
		    T1o = T1c - T1l;
	       }
	       O[WS(os, 1)] = FNMS(KP196034280, T1m, KP1_990369453 * T19);
	       O[WS(os, 25)] = FNMS(KP1_546020906, T1n, KP1_268786568 * T1o);
	       O[WS(os, 17)] = -(FMA(KP196034280, T19, KP1_990369453 * T1m));
	       O[WS(os, 9)] = FMA(KP1_268786568, T1n, KP1_546020906 * T1o);
	  }
	  {
	       E T1r, T1v, T1u, T1w;
	       {
		    E T1p, T1q, T1s, T1t;
		    T1p = Tz - TK;
		    T1q = T1b - T1a;
		    T1r = T1p + T1q;
		    T1v = T1p - T1q;
		    T1s = T1f - T1k;
		    T1t = TW - T17;
		    T1u = T1s - T1t;
		    T1w = T1t + T1s;
	       }
	       O[WS(os, 5)] = FMA(KP1_763842528, T1r, KP942793473 * T1u);
	       O[WS(os, 29)] = FNMS(KP1_913880671, T1v, KP580569354 * T1w);
	       O[WS(os, 21)] = FNMS(KP942793473, T1r, KP1_763842528 * T1u);
	       O[WS(os, 13)] = FMA(KP580569354, T1v, KP1_913880671 * T1w);
	  }
	  {
	       E T1T, T1X, T1W, T1Y;
	       {
		    E T1R, T1S, T1U, T1V;
		    T1R = T1x + T1y;
		    T1S = T1L + T1M;
		    T1T = T1R - T1S;
		    T1X = T1R + T1S;
		    T1U = T1J + T1I;
		    T1V = T1C - T1F;
		    T1W = T1U - T1V;
		    T1Y = T1V + T1U;
	       }
	       O[WS(os, 7)] = FMA(KP1_546020906, T1T, KP1_268786568 * T1W);
	       O[WS(os, 31)] = FNMS(KP1_990369453, T1X, KP196034280 * T1Y);
	       O[WS(os, 23)] = FNMS(KP1_268786568, T1T, KP1_546020906 * T1W);
	       O[WS(os, 15)] = FMA(KP196034280, T1X, KP1_990369453 * T1Y);
	  }
	  {
	       E T2f, T2p, T2o, T2q;
	       {
		    E T23, T2e, T2k, T2n;
		    T23 = T1Z + T22;
		    T2e = KP707106781 * (T28 + T2d);
		    T2f = T23 + T2e;
		    T2p = T23 - T2e;
		    T2k = T2i - T2j;
		    T2n = KP707106781 * (T2l + T2m);
		    T2o = T2k - T2n;
		    T2q = T2n + T2k;
	       }
	       O[WS(os, 2)] = FMA(KP1_961570560, T2f, KP390180644 * T2o);
	       O[WS(os, 26)] = FNMS(KP1_662939224, T2p, KP1_111140466 * T2q);
	       O[WS(os, 18)] = FNMS(KP390180644, T2f, KP1_961570560 * T2o);
	       O[WS(os, 10)] = FMA(KP1_111140466, T2p, KP1_662939224 * T2q);
	  }
	  {
	       E T1H, T1P, T1O, T1Q;
	       {
		    E T1z, T1G, T1K, T1N;
		    T1z = T1x - T1y;
		    T1G = T1C + T1F;
		    T1H = T1z + T1G;
		    T1P = T1z - T1G;
		    T1K = T1I - T1J;
		    T1N = T1L - T1M;
		    T1O = T1K - T1N;
		    T1Q = T1N + T1K;
	       }
	       O[WS(os, 3)] = FMA(KP1_913880671, T1H, KP580569354 * T1O);
	       O[WS(os, 27)] = FNMS(KP1_763842528, T1P, KP942793473 * T1Q);
	       O[WS(os, 19)] = FNMS(KP580569354, T1H, KP1_913880671 * T1O);
	       O[WS(os, 11)] = FMA(KP942793473, T1P, KP1_763842528 * T1Q);
	  }
     }
}

static const khc2r_desc desc = { 32, "hc2rIII_32", {138, 48, 36, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2rIII_32) (planner *p) {
     X(khc2rIII_register) (p, hc2rIII_32, &desc);
}

#endif				/* HAVE_FMA */
