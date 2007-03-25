/* tensor/tensor_float.h
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
#ifndef __TENSOR_float_H__
#define __TENSOR_float_H__

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
  float * data;
} tensor_float;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_float *
tensor_float_alloc(const unsigned int rank, const size_t dimension);

tensor_float *
tensor_float_calloc(const unsigned int rank, const size_t dimension);

tensor_float *
tensor_float_copy(tensor_float * t);

void tensor_float_free(tensor_float * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_float * tensor_float_2matrix(tensor_float * t);
gsl_vector_float * tensor_float_2vector(tensor_float * t);


/* Operations */

float tensor_float_get(const tensor_float * t, const size_t * indices);
void tensor_float_set(tensor_float * t, const size_t * indices, const float x);


float * tensor_float_ptr(tensor_float * t, const size_t * indices);
const float * tensor_float_const_ptr(const tensor_float * t, const size_t * indices);

void tensor_float_set_zero(tensor_float * t);
void tensor_float_set_all(tensor_float * t, float x);

int tensor_float_fread(FILE * stream, tensor_float * t);
int tensor_float_fwrite(FILE * stream, const tensor_float * t);
int tensor_float_fscanf(FILE * stream, tensor_float * t);
int tensor_float_fprintf(FILE * stream, const tensor_float * t,
                        const char * format);

int tensor_float_memcpy(tensor_float * dest, const tensor_float * src);
int tensor_float_swap(tensor_float * t1, tensor_float * t2);

tensor_float *
tensor_float_swap_indices(const tensor_float * t, size_t i, size_t j);

float tensor_float_max(const tensor_float * t);
float tensor_float_min(const tensor_float * t);
void tensor_float_minmax(const tensor_float * t,
                        float * min_out, float * max_out);

void tensor_float_max_index(const tensor_float * t, size_t * indices);
void tensor_float_min_index(const tensor_float * t, size_t * indices);
void tensor_float_minmax_index(const tensor_float * t,
                              size_t * imin, size_t * imax);

int tensor_float_isnull(const tensor_float * t);

int tensor_float_add(tensor_float * a, const tensor_float * b);
int tensor_float_sub(tensor_float * a, const tensor_float * b);
int tensor_float_mul_elements(tensor_float * a, const tensor_float * b);
int tensor_float_div_elements(tensor_float * a, const tensor_float * b);
int tensor_float_scale(tensor_float * a, const double x);
int tensor_float_add_constant(tensor_float * a, const double x);
int tensor_float_add_diagonal(tensor_float * a, const double x);
tensor_float * tensor_float_product(const tensor_float * a,
                                  const tensor_float * b);
tensor_float * tensor_float_contract(const tensor_float * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_float_position(const size_t * indices, const tensor_float * t)
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
float
tensor_float_get(const tensor_float * t, const size_t * indices)
{
  size_t position;

  position = tensor_float_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_float_set(tensor_float * t, const size_t * indices, const float x)
{
  size_t position;
  
  position = tensor_float_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
float *
tensor_float_ptr(tensor_float * t, const size_t * indices)
{
  size_t position;

  position = tensor_float_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (float *) (t->data + position);
}


extern inline 
const float *
tensor_float_const_ptr(const tensor_float * t, const size_t * indices)
{
  size_t position;

  position = tensor_float_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const float *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_float_H__ */
