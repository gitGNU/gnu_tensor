/* tensor/tensor_long.h
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
#ifndef __TENSOR_long_H__
#define __TENSOR_long_H__

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
  long * data;
} tensor_long;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_long *
tensor_long_alloc(const unsigned int rank, const size_t dimension);

tensor_long *
tensor_long_calloc(const unsigned int rank, const size_t dimension);

tensor_long *
tensor_long_copy(tensor_long * t);

void tensor_long_free(tensor_long * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_long * tensor_long_2matrix(tensor_long * t);
gsl_vector_long * tensor_long_2vector(tensor_long * t);


/* Operations */

long tensor_long_get(const tensor_long * t, const size_t * indices);
void tensor_long_set(tensor_long * t, const size_t * indices, const long x);


long * tensor_long_ptr(tensor_long * t, const size_t * indices);
const long * tensor_long_const_ptr(const tensor_long * t, const size_t * indices);

void tensor_long_set_zero(tensor_long * t);
void tensor_long_set_all(tensor_long * t, long x);

int tensor_long_fread(FILE * stream, tensor_long * t);
int tensor_long_fwrite(FILE * stream, const tensor_long * t);
int tensor_long_fscanf(FILE * stream, tensor_long * t);
int tensor_long_fprintf(FILE * stream, const tensor_long * t,
                        const char * format);

int tensor_long_memcpy(tensor_long * dest, const tensor_long * src);
int tensor_long_swap(tensor_long * t1, tensor_long * t2);

tensor_long *
tensor_long_swap_indices(const tensor_long * t, size_t i, size_t j);

long tensor_long_max(const tensor_long * t);
long tensor_long_min(const tensor_long * t);
void tensor_long_minmax(const tensor_long * t,
                        long * min_out, long * max_out);

void tensor_long_max_index(const tensor_long * t, size_t * indices);
void tensor_long_min_index(const tensor_long * t, size_t * indices);
void tensor_long_minmax_index(const tensor_long * t,
                              size_t * imin, size_t * imax);

int tensor_long_isnull(const tensor_long * t);

int tensor_long_add(tensor_long * a, const tensor_long * b);
int tensor_long_sub(tensor_long * a, const tensor_long * b);
int tensor_long_mul_elements(tensor_long * a, const tensor_long * b);
int tensor_long_div_elements(tensor_long * a, const tensor_long * b);
int tensor_long_scale(tensor_long * a, const double x);
int tensor_long_add_constant(tensor_long * a, const double x);
int tensor_long_add_diagonal(tensor_long * a, const double x);
tensor_long * tensor_long_product(const tensor_long * a,
                                  const tensor_long * b);
tensor_long * tensor_long_contract(const tensor_long * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_long_position(const size_t * indices, const tensor_long * t)
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
long
tensor_long_get(const tensor_long * t, const size_t * indices)
{
  size_t position;

  position = tensor_long_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_long_set(tensor_long * t, const size_t * indices, const long x)
{
  size_t position;
  
  position = tensor_long_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
long *
tensor_long_ptr(tensor_long * t, const size_t * indices)
{
  size_t position;

  position = tensor_long_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (long *) (t->data + position);
}


extern inline 
const long *
tensor_long_const_ptr(const tensor_long * t, const size_t * indices)
{
  size_t position;

  position = tensor_long_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const long *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_long_H__ */
