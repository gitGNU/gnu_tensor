/* tensor/tensor_ulong.h
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
#ifndef __TENSOR_ulong_H__
#define __TENSOR_ulong_H__

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
  unsigned long * data;
} tensor_ulong;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_ulong *
tensor_ulong_alloc(const unsigned int rank, const size_t dimension);

tensor_ulong *
tensor_ulong_calloc(const unsigned int rank, const size_t dimension);

tensor_ulong *
tensor_ulong_copy(tensor_ulong * t);

void tensor_ulong_free(tensor_ulong * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_ulong * tensor_ulong_2matrix(tensor_ulong * t);
gsl_vector_ulong * tensor_ulong_2vector(tensor_ulong * t);


/* Operations */

unsigned long tensor_ulong_get(const tensor_ulong * t, const size_t * indices);
void tensor_ulong_set(tensor_ulong * t, const size_t * indices, const unsigned long x);


unsigned long * tensor_ulong_ptr(tensor_ulong * t, const size_t * indices);
const unsigned long * tensor_ulong_const_ptr(const tensor_ulong * t, const size_t * indices);

void tensor_ulong_set_zero(tensor_ulong * t);
void tensor_ulong_set_all(tensor_ulong * t, unsigned long x);

int tensor_ulong_fread(FILE * stream, tensor_ulong * t);
int tensor_ulong_fwrite(FILE * stream, const tensor_ulong * t);
int tensor_ulong_fscanf(FILE * stream, tensor_ulong * t);
int tensor_ulong_fprintf(FILE * stream, const tensor_ulong * t,
                        const char * format);

int tensor_ulong_memcpy(tensor_ulong * dest, const tensor_ulong * src);
int tensor_ulong_swap(tensor_ulong * t1, tensor_ulong * t2);

tensor_ulong *
tensor_ulong_swap_indices(const tensor_ulong * t, size_t i, size_t j);

unsigned long tensor_ulong_max(const tensor_ulong * t);
unsigned long tensor_ulong_min(const tensor_ulong * t);
void tensor_ulong_minmax(const tensor_ulong * t,
                        unsigned long * min_out, unsigned long * max_out);

void tensor_ulong_max_index(const tensor_ulong * t, size_t * indices);
void tensor_ulong_min_index(const tensor_ulong * t, size_t * indices);
void tensor_ulong_minmax_index(const tensor_ulong * t,
                              size_t * imin, size_t * imax);

int tensor_ulong_isnull(const tensor_ulong * t);

int tensor_ulong_add(tensor_ulong * a, const tensor_ulong * b);
int tensor_ulong_sub(tensor_ulong * a, const tensor_ulong * b);
int tensor_ulong_mul_elements(tensor_ulong * a, const tensor_ulong * b);
int tensor_ulong_div_elements(tensor_ulong * a, const tensor_ulong * b);
int tensor_ulong_scale(tensor_ulong * a, const double x);
int tensor_ulong_add_constant(tensor_ulong * a, const double x);
int tensor_ulong_add_diagonal(tensor_ulong * a, const double x);
tensor_ulong * tensor_ulong_product(const tensor_ulong * a,
                                  const tensor_ulong * b);
tensor_ulong * tensor_ulong_contract(const tensor_ulong * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_ulong_position(const size_t * indices, const tensor_ulong * t)
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
unsigned long
tensor_ulong_get(const tensor_ulong * t, const size_t * indices)
{
  size_t position;

  position = tensor_ulong_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_ulong_set(tensor_ulong * t, const size_t * indices, const unsigned long x)
{
  size_t position;
  
  position = tensor_ulong_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
unsigned long *
tensor_ulong_ptr(tensor_ulong * t, const size_t * indices)
{
  size_t position;

  position = tensor_ulong_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (unsigned long *) (t->data + position);
}


extern inline 
const unsigned long *
tensor_ulong_const_ptr(const tensor_ulong * t, const size_t * indices)
{
  size_t position;

  position = tensor_ulong_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const unsigned long *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_ulong_H__ */
