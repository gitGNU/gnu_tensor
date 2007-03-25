/* tensor/tensor_short.h
 *
 * Copyright (C) 2003, 2004 Jordi Burguet-Castell
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * This file contains the basic declarations for a tensor.
 */
#ifndef __TENSOR_short_H__
#define __TENSOR_short_H__

#include <stdlib.h>


#include <gsl/gsl_types.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_check_range.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "tensor_utilities.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/*
 * A tensor is a struct with the rank of the tensor (number of indices
 * it has), the dimension of the vectorial space it is in (number of
 * different possible values for each index) and an array to store the
 * dimension^rank values.
 *
 * For the moment, there is no tda, as opossed to matrices, because it
 * would complicate quite a bit the algorithms and probably it is not
 * worth it.
 */
typedef struct
{
  unsigned int rank;
  size_t dimension;
  size_t size;
  short * data;
} tensor_short;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_short *
tensor_short_alloc(const unsigned int rank, const size_t dimension);

tensor_short *
tensor_short_calloc(const unsigned int rank, const size_t dimension);

tensor_short *
tensor_short_copy(tensor_short * t);

void tensor_short_free(tensor_short * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_short * tensor_short_2matrix(tensor_short * t);
gsl_vector_short * tensor_short_2vector(tensor_short * t);


/* Operations */

short tensor_short_get(const tensor_short * t, const size_t * indices);
void tensor_short_set(tensor_short * t, const size_t * indices, const short x);


short * tensor_short_ptr(tensor_short * t, const size_t * indices);
const short * tensor_short_const_ptr(const tensor_short * t, const size_t * indices);

void tensor_short_set_zero(tensor_short * t);
void tensor_short_set_all(tensor_short * t, short x);

int tensor_short_fread(FILE * stream, tensor_short * t);
int tensor_short_fwrite(FILE * stream, const tensor_short * t);
int tensor_short_fscanf(FILE * stream, tensor_short * t);
int tensor_short_fprintf(FILE * stream, const tensor_short * t,
                        const char * format);

int tensor_short_memcpy(tensor_short * dest, const tensor_short * src);
int tensor_short_swap(tensor_short * t1, tensor_short * t2);

tensor_short *
tensor_short_swap_indices(const tensor_short * t, size_t i, size_t j);

short tensor_short_max(const tensor_short * t);
short tensor_short_min(const tensor_short * t);
void tensor_short_minmax(const tensor_short * t,
                        short * min_out, short * max_out);

void tensor_short_max_index(const tensor_short * t, size_t * indices);
void tensor_short_min_index(const tensor_short * t, size_t * indices);
void tensor_short_minmax_index(const tensor_short * t,
                              size_t * imin, size_t * imax);

int tensor_short_isnull(const tensor_short * t);

int tensor_short_add(tensor_short * a, const tensor_short * b);
int tensor_short_sub(tensor_short * a, const tensor_short * b);
int tensor_short_mul_elements(tensor_short * a, const tensor_short * b);
int tensor_short_div_elements(tensor_short * a, const tensor_short * b);
int tensor_short_scale(tensor_short * a, const double x);
int tensor_short_add_constant(tensor_short * a, const double x);
int tensor_short_add_diagonal(tensor_short * a, const double x);
tensor_short * tensor_short_product(const tensor_short * a,
                                  const tensor_short * b);
tensor_short * tensor_short_contract(const tensor_short * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_short_position(const size_t * indices, const tensor_short * t)
{
  size_t shift, position;
  unsigned int i;

  shift = t->size/t->dimension;  /* == quick_pow(t->dimension, t->rank-1) */

  position = 0;
  for (i = 0; i < t->rank; i++)
    {
#if GSL_RANGE_CHECK
      if (indices[i] >= t->dimension)
        return t->size;
#endif

      position += shift * indices[i];
      shift /= t->dimension;
    }

  return position;
}


extern inline 
short
tensor_short_get(const tensor_short * t, const size_t * indices)
{
  size_t position;

  position = tensor_short_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_short_set(tensor_short * t, const size_t * indices, const short x)
{
  size_t position;
  
  position = tensor_short_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
short *
tensor_short_ptr(tensor_short * t, const size_t * indices)
{
  size_t position;

  position = tensor_short_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (short *) (t->data + position);
}


extern inline 
const short *
tensor_short_const_ptr(const tensor_short * t, const size_t * indices)
{
  size_t position;

  position = tensor_short_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const short *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_short_H__ */
