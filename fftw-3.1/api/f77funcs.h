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

/* Functions in the FFTW Fortran API, mangled according to the
   F77(...) macro.  This file is designed to be #included by
   f77api.c, possibly multiple times in order to support multiple
   compiler manglings (via redefinition of F77). */

void F77(execute, EXECUTE)(X(plan) * const p)
WITH_ALIGNED_STACK({
     plan *pln = (*p)->pln;
     pln->adt->solve(pln, (*p)->prb);
})

void F77(destroy_plan, DESTROY_PLAN)(X(plan) *p)
{
     X(destroy_plan)(*p);
}

void F77(cleanup, CLEANUP)(void)
{
     X(cleanup)();
}

void F77(forget_wisdom, FORGET_WISDOM)(void)
{
     X(forget_wisdom)();
}

void F77(export_wisdom, EXPORT_WISDOM)(void (*f77_write_char)(char *, void *),
				       void *data)
{
     write_char_data ad;
     ad.f77_write_char = f77_write_char;
     ad.data = data;
     X(export_wisdom)(write_char, (void *) &ad);
}

void F77(import_wisdom, IMPORT_WISDOM)(int *isuccess,
				       void (*f77_read_char)(int *, void *),
				       void *data)
{
     read_char_data ed;
     ed.f77_read_char = f77_read_char;
     ed.data = data;
     *isuccess = X(import_wisdom)(read_char, (void *) &ed);
}

void F77(import_system_wisdom, IMPORT_SYSTEM_WISDOM)(int *isuccess)
{
     *isuccess = X(import_system_wisdom)();
}

void F77(print_plan, PRINT_PLAN)(X(plan) * const p)
{
     X(print_plan)(*p);
     fflush(stdout);
}

void F77(flops,FLOPS)(X(plan) *p, double *add, double *mul, double *fma)
{
     X(flops)(*p, add, mul, fma);
}

/******************************** DFT ***********************************/

void F77(plan_dft, PLAN_DFT)(X(plan) *p, int *rank, const int *n,
			     C *in, C *out, int *sign, int *flags)
{
     int *nrev = reverse_n(*rank, n);
     *p = X(plan_dft)(*rank, nrev, in, out, *sign, *flags);
     X(ifree0)(nrev);
}

void F77(plan_dft_1d, PLAN_DFT_1D)(X(plan) *p, int *n, C *in, C *out,
				   int *sign, int *flags)
{
     *p = X(plan_dft_1d)(*n, in, out, *sign, *flags);
}

void F77(plan_dft_2d, PLAN_DFT_2D)(X(plan) *p, int *nx, int *ny,
				   C *in, C *out, int *sign, int *flags)
{
     *p = X(plan_dft_2d)(*ny, *nx, in, out, *sign, *flags);
}

void F77(plan_dft_3d, PLAN_DFT_3D)(X(plan) *p, int *nx, int *ny, int *nz,
				   C *in, C *out,
				   int *sign, int *flags)
{
     *p = X(plan_dft_3d)(*nz, *ny, *nx, in, out, *sign, *flags);
}

void F77(plan_many_dft, PLAN_MANY_DFT)(X(plan) *p, int *rank, const int *n,
				       int *howmany,
				       C *in, const int *inembed,
				       int *istride, int *idist,
				       C *out, const int *onembed,
				       int *ostride, int *odist,
				       int *sign, int *flags)
{
     int *nrev = reverse_n(*rank, n);
     int *inembedrev = reverse_n(*rank, inembed);
     int *onembedrev = reverse_n(*rank, onembed);
     *p = X(plan_many_dft)(*rank, nrev, *howmany,
			   in, inembedrev, *istride, *idist,
			   out, onembedrev, *ostride, *odist,
			   *sign, *flags);
     X(ifree0)(onembedrev);
     X(ifree0)(inembedrev);
     X(ifree0)(nrev);
}

void F77(plan_guru_dft, PLAN_GURU_DFT)(X(plan) *p, int *rank, const int *n,
				       const int *is, const int *os,
				       int *howmany_rank, const int *h_n,
				       const int *h_is, const int *h_os,
				       C *in, C *out, int *sign, int *flags)
{
     X(iodim) *dims = make_dims(*rank, n, is, os);
     X(iodim) *howmany_dims = make_dims(*howmany_rank, h_n, h_is, h_os);
     *p = X(plan_guru_dft)(*rank, dims, *howmany_rank, howmany_dims,
			   in, out, *sign, *flags);
     X(ifree0)(howmany_dims);
     X(ifree0)(dims);
}

void F77(plan_guru_split_dft, PLAN_GURU_SPLIT_DFT)(X(plan) *p, int *rank, const int *n,
				       const int *is, const int *os,
				       int *howmany_rank, const int *h_n,
				       const int *h_is, const int *h_os,
				       R *ri, R *ii, R *ro, R *io, int *flags)
{
     X(iodim) *dims = make_dims(*rank, n, is, os);
     X(iodim) *howmany_dims = make_dims(*howmany_rank, h_n, h_is, h_os);
     *p = X(plan_guru_split_dft)(*rank, dims, *howmany_rank, howmany_dims,
			   ri, ii, ro, io, *flags);
     X(ifree0)(howmany_dims);
     X(ifree0)(dims);
}

void F77(execute_dft, EXECUTE_DFT)(X(plan) * const p, C *in, C *out)
WITH_ALIGNED_STACK({
     plan_dft *pln = (plan_dft *) (*p)->pln;
     if ((*p)->sign == FFT_SIGN)
          pln->apply((plan *) pln, in[0], in[0]+1, out[0], out[0]+1);
     else
          pln->apply((plan *) pln, in[0]+1, in[0], out[0]+1, out[0]);
})

void F77(execute_split_dft, EXECUTE_SPLIT_DFT)(X(plan) * const p,
					       R *ri, R *ii, R *ro, R *io)
WITH_ALIGNED_STACK({
     plan_dft *pln = (plan_dft *) (*p)->pln;
     pln->apply((plan *) pln, ri, ii, ro, io);
})

/****************************** DFT r2c *********************************/

void F77(plan_dft_r2c, PLAN_DFT_R2C)(X(plan) *p, int *rank, const int *n,
				     R *in, C *out, int *flags)
{
     int *nrev = reverse_n(*rank, n);
     *p = X(plan_dft_r2c)(*rank, nrev, in, out, *flags);
     X(ifree0)(nrev);
}

void F77(plan_dft_r2c_1d, PLAN_DFT_R2C_1D)(X(plan) *p, int *n, R *in, C *out,
					   int *flags)
{
     *p = X(plan_dft_r2c_1d)(*n, in, out, *flags);
}

void F77(plan_dft_r2c_2d, PLAN_DFT_R2C_2D)(X(plan) *p, int *nx, int *ny,
					   R *in, C *out, int *flags)
{
     *p = X(plan_dft_r2c_2d)(*ny, *nx, in, out, *flags);
}

void F77(plan_dft_r2c_3d, PLAN_DFT_R2C_3D)(X(plan) *p,
					   int *nx, int *ny, int *nz,
					   R *in, C *out,
					   int *flags)
{
     *p = X(plan_dft_r2c_3d)(*nz, *ny, *nx, in, out, *flags);
}

void F77(plan_many_dft_r2c, PLAN_MANY_DFT_R2C)(
     X(plan) *p, int *rank, const int *n,
     int *howmany,
     R *in, const int *inembed, int *istride, int *idist,
     C *out, const int *onembed, int *ostride, int *odist,
     int *flags)
{
     int *nrev = reverse_n(*rank, n);
     int *inembedrev = reverse_n(*rank, inembed);
     int *onembedrev = reverse_n(*rank, onembed);
     *p = X(plan_many_dft_r2c)(*rank, nrev, *howmany,
			       in, inembedrev, *istride, *idist,
			       out, onembedrev, *ostride, *odist,
			       *flags);
     X(ifree0)(onembedrev);
     X(ifree0)(inembedrev);
     X(ifree0)(nrev);
}

void F77(plan_guru_dft_r2c, PLAN_GURU_DFT_R2C)(
     X(plan) *p, int *rank, const int *n,
     const int *is, const int *os,
     int *howmany_rank, const int *h_n,
     const int *h_is, const int *h_os,
     R *in, C *out, int *flags)
{
     X(iodim) *dims = make_dims(*rank, n, is, os);
     X(iodim) *howmany_dims = make_dims(*howmany_rank, h_n, h_is, h_os);
     *p = X(plan_guru_dft_r2c)(*rank, dims, *howmany_rank, howmany_dims,
			       in, out, *flags);
     X(ifree0)(howmany_dims);
     X(ifree0)(dims);
}

void F77(plan_guru_split_dft_r2c, PLAN_GURU_SPLIT_DFT_R2C)(
     X(plan) *p, int *rank, const int *n,
     const int *is, const int *os,
     int *howmany_rank, const int *h_n,
     const int *h_is, const int *h_os,
     R *in, R *ro, R *io, int *flags)
{
     X(iodim) *dims = make_dims(*rank, n, is, os);
     X(iodim) *howmany_dims = make_dims(*howmany_rank, h_n, h_is, h_os);
     *p = X(plan_guru_split_dft_r2c)(*rank, dims, *howmany_rank, howmany_dims,
			       in, ro, io, *flags);
     X(ifree0)(howmany_dims);
     X(ifree0)(dims);
}

void F77(execute_dft_r2c, EXECUTE_DFT_R2C)(X(plan) * const p, R *in, C *out)
WITH_ALIGNED_STACK({
     plan_rdft2 *pln = (plan_rdft2 *) (*p)->pln;
     pln->apply((plan *) pln, in, out[0], out[0]+1);
})

void F77(execute_split_dft_r2c, EXECUTE_SPLIT_DFT_R2C)(X(plan) * const p,
						       R *in, R *ro, R *io)
WITH_ALIGNED_STACK({
     plan_rdft2 *pln = (plan_rdft2 *) (*p)->pln;
     pln->apply((plan *) pln, in, ro, io);
})

/****************************** DFT c2r *********************************/

void F77(plan_dft_c2r, PLAN_DFT_C2R)(X(plan) *p, int *rank, const int *n,
				     C *in, R *out, int *flags)
{
     int *nrev = reverse_n(*rank, n);
     *p = X(plan_dft_c2r)(*rank, nrev, in, out, *flags);
     X(ifree0)(nrev);
}

void F77(plan_dft_c2r_1d, PLAN_DFT_C2R_1D)(X(plan) *p, int *n, C *in, R *out,
					   int *flags)
{
     *p = X(plan_dft_c2r_1d)(*n, in, out, *flags);
}

void F77(plan_dft_c2r_2d, PLAN_DFT_C2R_2D)(X(plan) *p, int *nx, int *ny,
					   C *in, R *out, int *flags)
{
     *p = X(plan_dft_c2r_2d)(*ny, *nx, in, out, *flags);
}

void F77(plan_dft_c2r_3d, PLAN_DFT_C2R_3D)(X(plan) *p,
					   int *nx, int *ny, int *nz,
					   C *in, R *out,
					   int *flags)
{
     *p = X(plan_dft_c2r_3d)(*nz, *ny, *nx, in, out, *flags);
}

void F77(plan_many_dft_c2r, PLAN_MANY_DFT_C2R)(
     X(plan) *p, int *rank, const int *n,
     int *howmany,
     C *in, const int *inembed, int *istride, int *idist,
     R *out, const int *onembed, int *ostride, int *odist,
     int *flags)
{
     int *nrev = reverse_n(*rank, n);
     int *inembedrev = reverse_n(*rank, inembed);
     int *onembedrev = reverse_n(*rank, onembed);
     *p = X(plan_many_dft_c2r)(*rank, nrev, *howmany,
			       in, inembedrev, *istride, *idist,
			       out, onembedrev, *ostride, *odist,
			       *flags);
     X(ifree0)(onembedrev);
     X(ifree0)(inembedrev);
     X(ifree0)(nrev);
}

void F77(plan_guru_dft_c2r, PLAN_GURU_DFT_C2R)(
     X(plan) *p, int *rank, const int *n,
     const int *is, const int *os,
     int *howmany_rank, const int *h_n,
     const int *h_is, const int *h_os,
     C *in, R *out, int *flags)
{
     X(iodim) *dims = make_dims(*rank, n, is, os);
     X(iodim) *howmany_dims = make_dims(*howmany_rank, h_n, h_is, h_os);
     *p = X(plan_guru_dft_c2r)(*rank, dims, *howmany_rank, howmany_dims,
			       in, out, *flags);
     X(ifree0)(howmany_dims);
     X(ifree0)(dims);
}

void F77(plan_guru_split_dft_c2r, PLAN_GURU_SPLIT_DFT_C2R)(
     X(plan) *p, int *rank, const int *n,
     const int *is, const int *os,
     int *howmany_rank, const int *h_n,
     const int *h_is, const int *h_os,
     R *ri, R *ii, R *out, int *flags)
{
     X(iodim) *dims = make_dims(*rank, n, is, os);
     X(iodim) *howmany_dims = make_dims(*howmany_rank, h_n, h_is, h_os);
     *p = X(plan_guru_split_dft_c2r)(*rank, dims, *howmany_rank, howmany_dims,
			       ri, ii, out, *flags);
     X(ifree0)(howmany_dims);
     X(ifree0)(dims);
}

void F77(execute_dft_c2r, EXECUTE_DFT_C2R)(X(plan) * const p, C *in, R *out)
WITH_ALIGNED_STACK({
     plan_rdft2 *pln = (plan_rdft2 *) (*p)->pln;
     pln->apply((plan *) pln, out, in[0], in[0]+1);
})

void F77(execute_split_dft_c2r, EXECUTE_SPLIT_DFT_C2R)(X(plan) * const p,
					   R *ri, R *ii, R *out)
WITH_ALIGNED_STACK({
     plan_rdft2 *pln = (plan_rdft2 *) (*p)->pln;
     pln->apply((plan *) pln, out, ri, ii);
})

/****************************** r2r *********************************/

void F77(plan_r2r, PLAN_R2R)(X(plan) *p, int *rank, const int *n,
			     R *in, R *out,
			     int *kind, int *flags)
{
     int *nrev = reverse_n(*rank, n);
     X(r2r_kind) *k = ints2kinds(*rank, kind);
     *p = X(plan_r2r)(*rank, nrev, in, out, k, *flags);
     X(ifree0)(k);
     X(ifree0)(nrev);
}

void F77(plan_r2r_1d, PLAN_R2R_1D)(X(plan) *p, int *n, R *in, R *out,
				   int *kind, int *flags)
{
     *p = X(plan_r2r_1d)(*n, in, out, (X(r2r_kind)) *kind, *flags);
}

void F77(plan_r2r_2d, PLAN_R2R_2D)(X(plan) *p, int *nx, int *ny,
				   R *in, R *out, 
				   int *kindx, int *kindy, int *flags)
{
     *p = X(plan_r2r_2d)(*ny, *nx, in, out,
			 (X(r2r_kind)) *kindy, (X(r2r_kind)) *kindx, *flags);
}

void F77(plan_r2r_3d, PLAN_R2R_3D)(X(plan) *p,
				   int *nx, int *ny, int *nz,
				   R *in, R *out,
				   int *kindx, int *kindy, int *kindz,
				   int *flags)
{
     *p = X(plan_r2r_3d)(*nz, *ny, *nx, in, out,
			 (X(r2r_kind)) *kindz, (X(r2r_kind)) *kindy, 
			 (X(r2r_kind)) *kindx, *flags);
}

void F77(plan_many_r2r, PLAN_MANY_R2R)(
     X(plan) *p, int *rank, const int *n,
     int *howmany,
     R *in, const int *inembed, int *istride, int *idist,
     R *out, const int *onembed, int *ostride, int *odist,
     int *kind, int *flags)
{
     int *nrev = reverse_n(*rank, n);
     int *inembedrev = reverse_n(*rank, inembed);
     int *onembedrev = reverse_n(*rank, onembed);
     X(r2r_kind) *k = ints2kinds(*rank, kind);
     *p = X(plan_many_r2r)(*rank, nrev, *howmany,
			       in, inembedrev, *istride, *idist,
			       out, onembedrev, *ostride, *odist,
			       k, *flags);
     X(ifree0)(k);
     X(ifree0)(onembedrev);
     X(ifree0)(inembedrev);
     X(ifree0)(nrev);
}

void F77(plan_guru_r2r, PLAN_GURU_R2R)(
     X(plan) *p, int *rank, const int *n,
     const int *is, const int *os,
     int *howmany_rank, const int *h_n,
     const int *h_is, const int *h_os,
     R *in, R *out, int *kind, int *flags)
{
     X(iodim) *dims = make_dims(*rank, n, is, os);
     X(iodim) *howmany_dims = make_dims(*howmany_rank, h_n, h_is, h_os);
     X(r2r_kind) *k = ints2kinds(*rank, kind);
     *p = X(plan_guru_r2r)(*rank, dims, *howmany_rank, howmany_dims,
			       in, out, k, *flags);
     X(ifree0)(k);
     X(ifree0)(howmany_dims);
     X(ifree0)(dims);
}

void F77(execute_r2r, EXECUTE_R2R)(X(plan) * const p, R *in, R *out)
WITH_ALIGNED_STACK({
     plan_rdft *pln = (plan_rdft *) (*p)->pln;
     pln->apply((plan *) pln, in, out);
})
