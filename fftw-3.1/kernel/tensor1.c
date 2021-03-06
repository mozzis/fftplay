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

/* $Id: tensor1.c,v 1.6 2006-01-05 03:04:27 stevenj Exp $ */

#include "ifftw.h"

tensor *X(mktensor_0d)(void)
{
     return X(mktensor(0));
}

tensor *X(mktensor_1d)(INT n, INT is, INT os)
{
     tensor *x = X(mktensor)(1);
     x->dims[0].n = n;
     x->dims[0].is = is;
     x->dims[0].os = os;
     return x;
}
