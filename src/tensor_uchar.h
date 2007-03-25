/* tensor/tensor_uchar.h
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
#ifndef __TENSOR_uchar_H__
#define __TENSOR_uchar_H__

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
  unsigned char * data;
} tensor_uchar;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor_uchar *
tensor_uchar_alloc(const unsigned int rank, const size_t dimension);

tensor_uchar *
tensor_uchar_calloc(const unsigned int rank, const size_t dimension);

tensor_uchar *
tensor_uchar_copy(tensor_uchar * t);

void tensor_uchar_free(tensor_uchar * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix_uchar * tensor_uchar_2matrix(tensor_uchar * t);
gsl_vector_uchar * tensor_uchar_2vector(tensor_uchar * t);


/* Operations */

unsigned char tensor_uchar_get(const tensor_uchar * t, const size_t * indices);
void tensor_uchar_set(tensor_uchar * t, const size_t * indices, const unsigned char x);


unsigned char * tensor_uchar_ptr(tensor_uchar * t, const size_t * indices);
const unsigned char * tensor_uchar_const_ptr(const tensor_uchar * t, const size_t * indices);

void tensor_uchar_set_zero(tensor_uchar * t);
void tensor_uchar_set_all(tensor_uchar * t, unsigned char x);

int tensor_uchar_fread(FILE * stream, tensor_uchar * t);
int tensor_uchar_fwrite(FILE * stream, const tensor_uchar * t);
int tensor_uchar_fscanf(FILE * stream, tensor_uchar * t);
int tensor_uchar_fprintf(FILE * stream, const tensor_uchar * t,
                        const char * format);

int tensor_uchar_memcpy(tensor_uchar * dest, const tensor_uchar * src);
int tensor_uchar_swap(tensor_uchar * t1, tensor_uchar * t2);

tensor_uchar *
tensor_uchar_swap_indices(const tensor_uchar * t, size_t i, size_t j);

unsigned char tensor_uchar_max(const tensor_uchar * t);
unsigned char tensor_uchar_min(const tensor_uchar * t);
void tensor_uchar_minmax(const tensor_uchar * t,
                        unsigned char * min_out, unsigned char * max_out);

void tensor_uchar_max_index(const tensor_uchar * t, size_t * indices);
void tensor_uchar_min_index(const tensor_uchar * t, size_t * indices);
void tensor_uchar_minmax_index(const tensor_uchar * t,
                              size_t * imin, size_t * imax);

int tensor_uchar_isnull(const tensor_uchar * t);

int tensor_uchar_add(tensor_uchar * a, const tensor_uchar * b);
int tensor_uchar_sub(tensor_uchar * a, const tensor_uchar * b);
int tensor_uchar_mul_elements(tensor_uchar * a, const tensor_uchar * b);
int tensor_uchar_div_elements(tensor_uchar * a, const tensor_uchar * b);
int tensor_uchar_scale(tensor_uchar * a, const double x);
int tensor_uchar_add_constant(tensor_uchar * a, const double x);
int tensor_uchar_add_diagonal(tensor_uchar * a, const double x);
tensor_uchar * tensor_uchar_product(const tensor_uchar * a,
                                  const tensor_uchar * b);
tensor_uchar * tensor_uchar_contract(const tensor_uchar * t_ij,
                                   size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_uchar_position(const size_t * indices, const tensor_uchar * t)
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
unsigned char
tensor_uchar_get(const tensor_uchar * t, const size_t * indices)
{
  size_t position;

  position = tensor_uchar_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}


extern inline
void
tensor_uchar_set(tensor_uchar * t, const size_t * indices, const unsigned char x)
{
  size_t position;
  
  position = tensor_uchar_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
unsigned char *
tensor_uchar_ptr(tensor_uchar * t, const size_t * indices)
{
  size_t position;

  position = tensor_uchar_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (unsigned char *) (t->data + position);
}


extern inline 
const unsigned char *
tensor_uchar_const_ptr(const tensor_uchar * t, const size_t * indices)
{
  size_t position;

  position = tensor_uchar_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const unsigned char *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_uchar_H__ */
