/* tensor/tensor_double.h
 *
 * Copyright (C) 2002, 2003, 2004, 2007 Jordi Burguet-Castell
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
#ifndef __TENSOR_DOUBLE_H__
#define __TENSOR_DOUBLE_H__

/* I am not sure they are neede here, but anyway */
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
  double * data;
} tensor;


/*
 * There is not such a thing as "tensor views", in contrast with the
 * case for gsl_matrix.
 */

/* Allocation */

tensor *
tensor_alloc(const unsigned int rank, const size_t dimension);

tensor *
tensor_calloc(const unsigned int rank, const size_t dimension);

tensor *
tensor_copy(tensor * t);

void tensor_free(tensor * t);


/* Views */

/*
 * There are no views.
 */


/* Conversions */

gsl_matrix * tensor_2matrix(tensor * t);
gsl_vector * tensor_2vector(tensor * t);


/* Operations */

double tensor_get(const tensor * t, const size_t * indices);
void tensor_set(tensor * t, const size_t * indices, const double x);


double * tensor_ptr(tensor * t, const size_t * indices);
const double * tensor_const_ptr(const tensor * t, const size_t * indices);

void tensor_set_zero(tensor * t);
void tensor_set_all(tensor * t, double x);

int tensor_fread(FILE * stream, tensor * t);
int tensor_fwrite(FILE * stream, const tensor * t);
int tensor_fscanf(FILE * stream, tensor * t);
int tensor_fprintf(FILE * stream, const tensor * t, const char * format);

int tensor_memcpy(tensor * dest, const tensor * src);
int tensor_swap(tensor * t1, tensor * t2);

tensor *
tensor_swap_indices(const tensor * t_ij, size_t i, size_t j);

double tensor_max(const tensor * t);
double tensor_min(const tensor * t);
void tensor_minmax(const tensor * t, double * min_out, double * max_out);

void tensor_max_index(const tensor * t, size_t * indices);
void tensor_min_index(const tensor * t, size_t * indices);
void tensor_minmax_index(const tensor * t, size_t * imin, size_t * imax);

int tensor_isnull(const tensor * t);

int tensor_add(tensor * a, const tensor * b);
int tensor_sub(tensor * a, const tensor * b);
int tensor_mul_elements(tensor * a, const tensor * b);
int tensor_div_elements(tensor * a, const tensor * b);
int tensor_scale(tensor * a, const double x);
int tensor_add_constant(tensor * a, const double x);
int tensor_add_diagonal(tensor * a, const double x);
tensor * tensor_product(const tensor * a, const tensor * b);
tensor * tensor_contract(const tensor * t_ij, size_t i, size_t j);


/* inline functions if you are using GCC */

#ifdef HAVE_INLINE
extern inline
size_t
tensor_position(const size_t * indices, const tensor * t)
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
double
tensor_get(const tensor * t, const size_t * indices)
{
  size_t position;

  position = tensor_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VAL("index out of range", GSL_EINVAL, 0);
#endif

  return t->data[position];
}

extern inline
void
tensor_set(tensor * t, const size_t * indices, const double x)
{
  size_t position;

  position = tensor_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_VOID("index out of range", GSL_EINVAL);
#endif

  t->data[position] = x;
}


extern inline 
double *
tensor_ptr(tensor * t, const size_t * indices)
{
  size_t position;

  position = tensor_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (double *) (t->data + position);
} 

extern inline 
const double *
tensor_const_ptr(const tensor * t, const size_t * indices)
{
  size_t position;

  position = tensor_position(indices, t);
#if GSL_RANGE_CHECK
  if (position >= t->size)
    GSL_ERROR_NULL("index out of range", GSL_EINVAL);
#endif

  return (const double *) (t->data + position);
} 

#endif

__END_DECLS

#endif /* __TENSOR_DOUBLE_H__ */
