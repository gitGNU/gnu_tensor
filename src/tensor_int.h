/* tensor/tensor_int.h
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
#ifndef __TENSOR_int_H__
#define __TENSOR_int_H__

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
  int * data;
} tensor_int;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_int *
tensor_int_alloc(const unsigned int rank, const size_t dimension);

tensor_int *
tensor_int_calloc(const unsigned int rank, const size_t dimension);

tensor_int *
tensor_int_copy(tensor_int * t);

void tensor_int_free(tensor_int * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_int * tensor_int_2matrix(tensor_int * t);
gsl_vector_int * tensor_int_2vector(tensor_int * t);


/* Operations */

int tensor_int_get(const tensor_int * t, const size_t * indices);
void tensor_int_set(tensor_int * t, const size_t * indices, const int x);


int * tensor_int_ptr(tensor_int * t, const size_t * indices);
const int * tensor_int_const_ptr(const tensor_int * t, const size_t * indices);

void tensor_int_set_zero(tensor_int * t);
void tensor_int_set_all(tensor_int * t, int x);

int tensor_int_fread(FILE * stream, tensor_int * t);
int tensor_int_fwrite(FILE * stream, const tensor_int * t);
int tensor_int_fscanf(FILE * stream, tensor_int * t);
int tensor_int_fprintf(FILE * stream, const tensor_int * t,
                        const char * format);

int tensor_int_memcpy(tensor_int * dest, const tensor_int * src);
int tensor_int_swap(tensor_int * t1, tensor_int * t2);

tensor_int *
tensor_int_swap_indices(const tensor_int * t, size_t i, size_t j);

int tensor_int_max(const tensor_int * t);
int tensor_int_min(const tensor_int * t);
void tensor_int_minmax(const tensor_int * t,
                        int * min_out, int * max_out);

void tensor_int_max_index(const tensor_int * t, size_t * indices);
void tensor_int_min_index(const tensor_int * t, size_t * indices);
void tensor_int_minmax_index(const tensor_int * t,
                              size_t * imin, size_t * imax);

int tensor_int_isnull(const tensor_int * t);

int tensor_int_add(tensor_int * a, const tensor_int * b);
int tensor_int_sub(tensor_int * a, const tensor_int * b);
int tensor_int_mul_elements(tensor_int * a, const tensor_int * b);
int tensor_int_div_elements(tensor_int * a, const tensor_int * b);
int tensor_int_scale(tensor_int * a, const double x);
int tensor_int_add_constant(tensor_int * a, const double x);
int tensor_int_add_diagonal(tensor_int * a, const double x);
tensor_int * tensor_int_product(const tensor_int * a,
                                  const tensor_int * b);
tensor_int * tensor_int_contract(const tensor_int * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_int_position(const size_t * indices, const tensor_int * t)
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
int
tensor_int_get(const tensor_int * t, const size_t * indices)
{
  size_t position;

  position = tensor_int_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_int_set(tensor_int * t, const size_t * indices, const int x)
{
  size_t position;
  
  position = tensor_int_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
int *
tensor_int_ptr(tensor_int * t, const size_t * indices)
{
  size_t position;

  position = tensor_int_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (int *) (t->data + position);
}


extern inline 
const int *
tensor_int_const_ptr(const tensor_int * t, const size_t * indices)
{
  size_t position;

  position = tensor_int_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const int *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_int_H__ */
