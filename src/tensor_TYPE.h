/* tensor/tensor_NAME.h
 *
 * Copyright (C) 2003, 2004, 2007 Jordi Burguet-Castell
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
#ifndef __TENSOR_NAME_H__
#define __TENSOR_NAME_H__

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
  TYPE * data;
} tensor_NAME;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_NAME *
tensor_NAME_alloc(const unsigned int rank, const size_t dimension);

tensor_NAME *
tensor_NAME_calloc(const unsigned int rank, const size_t dimension);

tensor_NAME *
tensor_NAME_copy(tensor_NAME * t);

void tensor_NAME_free(tensor_NAME * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_NAME * tensor_NAME_2matrix(tensor_NAME * t);
gsl_vector_NAME * tensor_NAME_2vector(tensor_NAME * t);


/* Operations */

TYPE tensor_NAME_get(const tensor_NAME * t, const size_t * indices);
void tensor_NAME_set(tensor_NAME * t, const size_t * indices, const TYPE x);


TYPE * tensor_NAME_ptr(tensor_NAME * t, const size_t * indices);
const TYPE * tensor_NAME_const_ptr(const tensor_NAME * t, const size_t * indices);

void tensor_NAME_set_zero(tensor_NAME * t);
void tensor_NAME_set_all(tensor_NAME * t, TYPE x);

int tensor_NAME_fread(FILE * stream, tensor_NAME * t);
int tensor_NAME_fwrite(FILE * stream, const tensor_NAME * t);
int tensor_NAME_fscanf(FILE * stream, tensor_NAME * t);
int tensor_NAME_fprintf(FILE * stream, const tensor_NAME * t,
                        const char * format);

int tensor_NAME_memcpy(tensor_NAME * dest, const tensor_NAME * src);
int tensor_NAME_swap(tensor_NAME * t1, tensor_NAME * t2);

tensor_NAME *
tensor_NAME_swap_indices(const tensor_NAME * t, size_t i, size_t j);

TYPE tensor_NAME_max(const tensor_NAME * t);
TYPE tensor_NAME_min(const tensor_NAME * t);
void tensor_NAME_minmax(const tensor_NAME * t,
                        TYPE * min_out, TYPE * max_out);

void tensor_NAME_max_index(const tensor_NAME * t, size_t * indices);
void tensor_NAME_min_index(const tensor_NAME * t, size_t * indices);
void tensor_NAME_minmax_index(const tensor_NAME * t,
                              size_t * imin, size_t * imax);

int tensor_NAME_isnull(const tensor_NAME * t);

int tensor_NAME_add(tensor_NAME * a, const tensor_NAME * b);
int tensor_NAME_sub(tensor_NAME * a, const tensor_NAME * b);
int tensor_NAME_mul_elements(tensor_NAME * a, const tensor_NAME * b);
int tensor_NAME_div_elements(tensor_NAME * a, const tensor_NAME * b);
int tensor_NAME_scale(tensor_NAME * a, const double x);
int tensor_NAME_add_constant(tensor_NAME * a, const double x);
int tensor_NAME_add_diagonal(tensor_NAME * a, const double x);
tensor_NAME * tensor_NAME_product(const tensor_NAME * a,
                                  const tensor_NAME * b);
tensor_NAME * tensor_NAME_contract(const tensor_NAME * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_NAME_position(const size_t * indices, const tensor_NAME * t)
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
TYPE
tensor_NAME_get(const tensor_NAME * t, const size_t * indices)
{
  size_t position;

  position = tensor_NAME_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_NAME_set(tensor_NAME * t, const size_t * indices, const TYPE x)
{
  size_t position;
  
  position = tensor_NAME_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
TYPE *
tensor_NAME_ptr(tensor_NAME * t, const size_t * indices)
{
  size_t position;

  position = tensor_NAME_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (TYPE *) (t->data + position);
}


extern inline 
const TYPE *
tensor_NAME_const_ptr(const tensor_NAME * t, const size_t * indices)
{
  size_t position;

  position = tensor_NAME_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const TYPE *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_NAME_H__ */
