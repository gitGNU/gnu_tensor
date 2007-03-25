/* tensor/tensor_ushort.h
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
#ifndef __TENSOR_ushort_H__
#define __TENSOR_ushort_H__

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
  unsigned short * data;
} tensor_ushort;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_ushort *
tensor_ushort_alloc(const unsigned int rank, const size_t dimension);

tensor_ushort *
tensor_ushort_calloc(const unsigned int rank, const size_t dimension);

tensor_ushort *
tensor_ushort_copy(tensor_ushort * t);

void tensor_ushort_free(tensor_ushort * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_ushort * tensor_ushort_2matrix(tensor_ushort * t);
gsl_vector_ushort * tensor_ushort_2vector(tensor_ushort * t);


/* Operations */

unsigned short tensor_ushort_get(const tensor_ushort * t, const size_t * indices);
void tensor_ushort_set(tensor_ushort * t, const size_t * indices, const unsigned short x);


unsigned short * tensor_ushort_ptr(tensor_ushort * t, const size_t * indices);
const unsigned short * tensor_ushort_const_ptr(const tensor_ushort * t, const size_t * indices);

void tensor_ushort_set_zero(tensor_ushort * t);
void tensor_ushort_set_all(tensor_ushort * t, unsigned short x);

int tensor_ushort_fread(FILE * stream, tensor_ushort * t);
int tensor_ushort_fwrite(FILE * stream, const tensor_ushort * t);
int tensor_ushort_fscanf(FILE * stream, tensor_ushort * t);
int tensor_ushort_fprintf(FILE * stream, const tensor_ushort * t,
                        const char * format);

int tensor_ushort_memcpy(tensor_ushort * dest, const tensor_ushort * src);
int tensor_ushort_swap(tensor_ushort * t1, tensor_ushort * t2);

tensor_ushort *
tensor_ushort_swap_indices(const tensor_ushort * t, size_t i, size_t j);

unsigned short tensor_ushort_max(const tensor_ushort * t);
unsigned short tensor_ushort_min(const tensor_ushort * t);
void tensor_ushort_minmax(const tensor_ushort * t,
                        unsigned short * min_out, unsigned short * max_out);

void tensor_ushort_max_index(const tensor_ushort * t, size_t * indices);
void tensor_ushort_min_index(const tensor_ushort * t, size_t * indices);
void tensor_ushort_minmax_index(const tensor_ushort * t,
                              size_t * imin, size_t * imax);

int tensor_ushort_isnull(const tensor_ushort * t);

int tensor_ushort_add(tensor_ushort * a, const tensor_ushort * b);
int tensor_ushort_sub(tensor_ushort * a, const tensor_ushort * b);
int tensor_ushort_mul_elements(tensor_ushort * a, const tensor_ushort * b);
int tensor_ushort_div_elements(tensor_ushort * a, const tensor_ushort * b);
int tensor_ushort_scale(tensor_ushort * a, const double x);
int tensor_ushort_add_constant(tensor_ushort * a, const double x);
int tensor_ushort_add_diagonal(tensor_ushort * a, const double x);
tensor_ushort * tensor_ushort_product(const tensor_ushort * a,
                                  const tensor_ushort * b);
tensor_ushort * tensor_ushort_contract(const tensor_ushort * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_ushort_position(const size_t * indices, const tensor_ushort * t)
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
unsigned short
tensor_ushort_get(const tensor_ushort * t, const size_t * indices)
{
  size_t position;

  position = tensor_ushort_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_ushort_set(tensor_ushort * t, const size_t * indices, const unsigned short x)
{
  size_t position;
  
  position = tensor_ushort_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
unsigned short *
tensor_ushort_ptr(tensor_ushort * t, const size_t * indices)
{
  size_t position;

  position = tensor_ushort_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (unsigned short *) (t->data + position);
}


extern inline 
const unsigned short *
tensor_ushort_const_ptr(const tensor_ushort * t, const size_t * indices)
{
  size_t position;

  position = tensor_ushort_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const unsigned short *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_ushort_H__ */
