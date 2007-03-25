/* tensor/tensor_long_double.h
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
#ifndef __TENSOR_long_double_H__
#define __TENSOR_long_double_H__

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
  long double * data;
} tensor_long_double;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_long_double *
tensor_long_double_alloc(const unsigned int rank, const size_t dimension);

tensor_long_double *
tensor_long_double_calloc(const unsigned int rank, const size_t dimension);

tensor_long_double *
tensor_long_double_copy(tensor_long_double * t);

void tensor_long_double_free(tensor_long_double * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_long_double * tensor_long_double_2matrix(tensor_long_double * t);
gsl_vector_long_double * tensor_long_double_2vector(tensor_long_double * t);


/* Operations */

long double tensor_long_double_get(const tensor_long_double * t, const size_t * indices);
void tensor_long_double_set(tensor_long_double * t, const size_t * indices, const long double x);


long double * tensor_long_double_ptr(tensor_long_double * t, const size_t * indices);
const long double * tensor_long_double_const_ptr(const tensor_long_double * t, const size_t * indices);

void tensor_long_double_set_zero(tensor_long_double * t);
void tensor_long_double_set_all(tensor_long_double * t, long double x);

int tensor_long_double_fread(FILE * stream, tensor_long_double * t);
int tensor_long_double_fwrite(FILE * stream, const tensor_long_double * t);
int tensor_long_double_fscanf(FILE * stream, tensor_long_double * t);
int tensor_long_double_fprintf(FILE * stream, const tensor_long_double * t,
                        const char * format);

int tensor_long_double_memcpy(tensor_long_double * dest, const tensor_long_double * src);
int tensor_long_double_swap(tensor_long_double * t1, tensor_long_double * t2);

tensor_long_double *
tensor_long_double_swap_indices(const tensor_long_double * t, size_t i, size_t j);

long double tensor_long_double_max(const tensor_long_double * t);
long double tensor_long_double_min(const tensor_long_double * t);
void tensor_long_double_minmax(const tensor_long_double * t,
                        long double * min_out, long double * max_out);

void tensor_long_double_max_index(const tensor_long_double * t, size_t * indices);
void tensor_long_double_min_index(const tensor_long_double * t, size_t * indices);
void tensor_long_double_minmax_index(const tensor_long_double * t,
                              size_t * imin, size_t * imax);

int tensor_long_double_isnull(const tensor_long_double * t);

int tensor_long_double_add(tensor_long_double * a, const tensor_long_double * b);
int tensor_long_double_sub(tensor_long_double * a, const tensor_long_double * b);
int tensor_long_double_mul_elements(tensor_long_double * a, const tensor_long_double * b);
int tensor_long_double_div_elements(tensor_long_double * a, const tensor_long_double * b);
int tensor_long_double_scale(tensor_long_double * a, const double x);
int tensor_long_double_add_constant(tensor_long_double * a, const double x);
int tensor_long_double_add_diagonal(tensor_long_double * a, const double x);
tensor_long_double * tensor_long_double_product(const tensor_long_double * a,
                                  const tensor_long_double * b);
tensor_long_double * tensor_long_double_contract(const tensor_long_double * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_long_double_position(const size_t * indices, const tensor_long_double * t)
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
long double
tensor_long_double_get(const tensor_long_double * t, const size_t * indices)
{
  size_t position;

  position = tensor_long_double_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_long_double_set(tensor_long_double * t, const size_t * indices, const long double x)
{
  size_t position;
  
  position = tensor_long_double_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
long double *
tensor_long_double_ptr(tensor_long_double * t, const size_t * indices)
{
  size_t position;

  position = tensor_long_double_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (long double *) (t->data + position);
}


extern inline 
const long double *
tensor_long_double_const_ptr(const tensor_long_double * t, const size_t * indices)
{
  size_t position;

  position = tensor_long_double_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const long double *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_long_double_H__ */
