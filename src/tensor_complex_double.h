/* tensor/tensor_complex_double.h
 *
 * Copyright (C) 2010 Jordi Burguet-Castell
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
 * along with this program; if not, write to:
 *   Free Software Foundation, Inc.
 *   51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 */

/*
 * This file contains the basic declarations for a tensor.
 */
#ifndef __TENSOR_COMPLEX_DOUBLE_H__
#define __TENSOR_COMPLEX_DOUBLE_H__

/* I am not sure they are neede here, but anyway */
#include <stdlib.h>
#include <complex.h>

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
  complex double * data;
} tensor_complex;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_complex *
tensor_complex_alloc(const unsigned int rank, const size_t dimension);

tensor_complex *
tensor_complex_calloc(const unsigned int rank, const size_t dimension);

tensor_complex *
tensor_complex_copy(tensor_complex * t);

void tensor_complex_free(tensor_complex * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_complex * tensor_complex_2matrix(tensor_complex * t);
gsl_vector_complex * tensor_complex_2vector(tensor_complex * t);


/* Operations */

complex double tensor_complex_get(const tensor_complex * t, const size_t * indices);
void tensor_complex_set(tensor_complex * t, const size_t * indices, const complex double x);


complex double * tensor_complex_ptr(tensor_complex * t, const size_t * indices);
const complex double * tensor_complex_const_ptr(const tensor_complex * t, const size_t * indices);

void tensor_complex_set_zero(tensor_complex * t);
void tensor_complex_set_all(tensor_complex * t, complex double x);

int tensor_complex_fread(FILE * stream, tensor_complex * t);
int tensor_complex_fwrite(FILE * stream, const tensor_complex * t);
int tensor_complex_fscanf(FILE * stream, tensor_complex * t);
int tensor_complex_fprintf(FILE * stream, const tensor_complex * t, const char * format);

int tensor_complex_memcpy(tensor_complex * dest, const tensor_complex * src);
int tensor_complex_swap(tensor_complex * t1, tensor_complex * t2);

tensor_complex *
tensor_complex_swap_indices(const tensor_complex * t_ij, size_t i, size_t j);

int tensor_complex_isnull(const tensor_complex * t);

int tensor_complex_add(tensor_complex * a, const tensor_complex * b);
int tensor_complex_sub(tensor_complex * a, const tensor_complex * b);
int tensor_complex_mul_elements(tensor_complex * a, const tensor_complex * b);
int tensor_complex_div_elements(tensor_complex * a, const tensor_complex * b);
int tensor_complex_scale(tensor_complex * a, const double x);
int tensor_complex_add_constant(tensor_complex * a, const double x);
int tensor_complex_add_diagonal(tensor_complex * a, const double x);
tensor_complex * tensor_complex_product(const tensor_complex * a, const tensor_complex * b);
tensor_complex * tensor_complex_contract(const tensor_complex * t_ij, size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_complex_position(const size_t * indices, const tensor_complex * t)
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
complex double
tensor_complex_get(const tensor_complex * t, const size_t * indices)
{
  size_t position;

  position = tensor_complex_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}

extern inline
void
tensor_complex_set(tensor_complex * t, const size_t * indices, const complex double x)
{
  size_t position;

  position = tensor_complex_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
complex double *
tensor_complex_ptr(tensor_complex * t, const size_t * indices)
{
  size_t position;

  position = tensor_complex_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (complex double *) (t->data + position);
} 

extern inline 
const complex double *
tensor_complex_const_ptr(const tensor_complex * t, const size_t * indices)
{
  size_t position;

  position = tensor_complex_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const complex double *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_COMPLEX_DOUBLE_H__ */
