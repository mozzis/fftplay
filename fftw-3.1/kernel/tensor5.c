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

/* $Id: tensor5.c,v 1.4 2006-01-05 03:04:27 stevenj Exp $ */

#include "ifftw.h"

static void dimcpy(iodim *dst, const iodim *src, int rnk)
{
     int i;
     if (FINITE_RNK(rnk))
          for (i = 0; i < rnk; ++i)
               dst[i] = src[i];
}

tensor *X(tensor_copy)(const tensor *sz)
{
     tensor *x = X(mktensor)(sz->rnk);
     dimcpy(x->dims, sz->dims, sz->rnk);
     return x;
}

/* like X(tensor_copy), but makes strides in-place by
   setting os = is if k == INPLACE_IS or is = os if k == INPLACE_OS. */
tensor *X(tensor_copy_inplace)(const tensor *sz, inplace_kind k)
{
     tensor *x = X(tensor_copy)(sz);
     if (FINITE_RNK(x->rnk)) {
	  int i;
	  if (k == INPLACE_OS)
	       for (i = 0; i < x->rnk; ++i)
		    x->dims[i].is = x->dims[i].os;
	  else
	       for (i = 0; i < x->rnk; ++i)
		    x->dims[i].os = x->dims[i].is;
     }
     return x;
}

/* Like X(tensor_copy), but copy all of the dimensions *except*
   except_dim. */
tensor *X(tensor_copy_except)(const tensor *sz, int except_dim)
{
     tensor *x;

     A(FINITE_RNK(sz->rnk) && sz->rnk >= 1 && except_dim < sz->rnk);
     x = X(mktensor)(sz->rnk - 1);
     dimcpy(x->dims, sz->dims, except_dim);
     dimcpy(x->dims + except_dim, sz->dims + except_dim + 1,
            x->rnk - except_dim);
     return x;
}

/* Like X(tensor_copy), but copy only rnk dimensions starting
   with start_dim. */
tensor *X(tensor_copy_sub)(const tensor *sz, int start_dim, int rnk)
{
     tensor *x;

     A(FINITE_RNK(sz->rnk) && start_dim + rnk <= sz->rnk);
     x = X(mktensor)(rnk);
     dimcpy(x->dims, sz->dims + start_dim, rnk);
     return x;
}

tensor *X(tensor_append)(const tensor *a, const tensor *b)
{
     if (!FINITE_RNK(a->rnk) || !FINITE_RNK(b->rnk)) {
          return X(mktensor)(RNK_MINFTY);
     } else {
	  tensor *x = X(mktensor)(a->rnk + b->rnk);
          dimcpy(x->dims, a->dims, a->rnk);
          dimcpy(x->dims + a->rnk, b->dims, b->rnk);
	  return x;
     }
}
