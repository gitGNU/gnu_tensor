/* tensor/tensor_uint.h
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
#ifndef __TENSOR_uint_H__
#define __TENSOR_uint_H__

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
  unsigned int * data;
} tensor_uint;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_uint *
tensor_uint_alloc(const unsigned int rank, const size_t dimension);

tensor_uint *
tensor_uint_calloc(const unsigned int rank, const size_t dimension);

tensor_uint *
tensor_uint_copy(tensor_uint * t);

void tensor_uint_free(tensor_uint * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_uint * tensor_uint_2matrix(tensor_uint * t);
gsl_vector_uint * tensor_uint_2vector(tensor_uint * t);


/* Operations */

unsigned int tensor_uint_get(const tensor_uint * t, const size_t * indices);
void tensor_uint_set(tensor_uint * t, const size_t * indices, const unsigned int x);


unsigned int * tensor_uint_ptr(tensor_uint * t, const size_t * indices);
const unsigned int * tensor_uint_const_ptr(const tensor_uint * t, const size_t * indices);

void tensor_uint_set_zero(tensor_uint * t);
void tensor_uint_set_all(tensor_uint * t, unsigned int x);

int tensor_uint_fread(FILE * stream, tensor_uint * t);
int tensor_uint_fwrite(FILE * stream, const tensor_uint * t);
int tensor_uint_fscanf(FILE * stream, tensor_uint * t);
int tensor_uint_fprintf(FILE * stream, const tensor_uint * t,
                        const char * format);

int tensor_uint_memcpy(tensor_uint * dest, const tensor_uint * src);
int tensor_uint_swap(tensor_uint * t1, tensor_uint * t2);

tensor_uint *
tensor_uint_swap_indices(const tensor_uint * t, size_t i, size_t j);

unsigned int tensor_uint_max(const tensor_uint * t);
unsigned int tensor_uint_min(const tensor_uint * t);
void tensor_uint_minmax(const tensor_uint * t,
                        unsigned int * min_out, unsigned int * max_out);

void tensor_uint_max_index(const tensor_uint * t, size_t * indices);
void tensor_uint_min_index(const tensor_uint * t, size_t * indices);
void tensor_uint_minmax_index(const tensor_uint * t,
                              size_t * imin, size_t * imax);

int tensor_uint_isnull(const tensor_uint * t);

int tensor_uint_add(tensor_uint * a, const tensor_uint * b);
int tensor_uint_sub(tensor_uint * a, const tensor_uint * b);
int tensor_uint_mul_elements(tensor_uint * a, const tensor_uint * b);
int tensor_uint_div_elements(tensor_uint * a, const tensor_uint * b);
int tensor_uint_scale(tensor_uint * a, const double x);
int tensor_uint_add_constant(tensor_uint * a, const double x);
int tensor_uint_add_diagonal(tensor_uint * a, const double x);
tensor_uint * tensor_uint_product(const tensor_uint * a,
                                  const tensor_uint * b);
tensor_uint * tensor_uint_contract(const tensor_uint * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_uint_position(const size_t * indices, const tensor_uint * t)
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
unsigned int
tensor_uint_get(const tensor_uint * t, const size_t * indices)
{
  size_t position;

  position = tensor_uint_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_uint_set(tensor_uint * t, const size_t * indices, const unsigned int x)
{
  size_t position;
  
  position = tensor_uint_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
unsigned int *
tensor_uint_ptr(tensor_uint * t, const size_t * indices)
{
  size_t position;

  position = tensor_uint_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (unsigned int *) (t->data + position);
}


extern inline 
const unsigned int *
tensor_uint_const_ptr(const tensor_uint * t, const size_t * indices)
{
  size_t position;

  position = tensor_uint_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const unsigned int *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_uint_H__ */
