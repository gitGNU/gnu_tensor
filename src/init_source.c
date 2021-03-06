/* tensor/init_source.c
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
 * This code follows as close as possible that of init_source.c in
 * the gsl/matrix directory
 */

#include "tensor_utilities.h"

#include <string.h>  /* to use memcpy() */


/* ------ Allocation ------ */

/*
 * Allocate memory for a tensor and return a pointer to it.
 */
TYPE(tensor) *
FUNCTION(tensor, alloc) (const unsigned int rank, const size_t dimension)
{
  size_t n;
  TYPE(tensor) * t;

  if (dimension == 0)
    {
      GSL_ERROR_VAL ("tensor dimension must be positive integer",
		     GSL_EINVAL, 0);
    }
  
  t = (TYPE(tensor) *) malloc (sizeof (TYPE(tensor)));

  if (t == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for tensor struct",
		     GSL_ENOMEM, 0);
    }

  n = quick_pow(dimension, rank);
  t->data = (ATOMIC *) malloc (n * sizeof (ATOMIC));

  if (t->data == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for data",
		     GSL_ENOMEM, 0);
    }

  t->rank = rank;
  t->dimension = dimension;
  t->size = n;

  return t;
}


/*
 * Same as tensor_alloc, but put all elements to 0.
 */
TYPE(tensor) *
FUNCTION(tensor, calloc) (const unsigned int rank, const size_t dimension)
{
  size_t i;
  size_t n;

  TYPE(tensor) * t = FUNCTION(tensor, alloc) (rank, dimension);

  if (t == 0)
    return NULL;
  
  /* initialize tensor to zero */
  
  n = t->size;
  for (i = 0; i < n; i++)
    t->data[i] = 0;

  return t;
}


/*
 * Copy from an existing tensor.
 */
TYPE(tensor) *
FUNCTION(tensor, copy) (TYPE(tensor) * tt)
{
  TYPE(tensor) * t = FUNCTION(tensor, alloc) (tt->rank, tt->dimension);

  if (t == 0)
    return NULL;
  
  t->rank = tt->rank;
  t->dimension = tt->dimension;
  t->size = tt->size;
  memcpy(t->data, tt->data, sizeof(BASE) * tt->size);

  return t;
}


/*
 * Free memory.
 */
void
FUNCTION(tensor, free) (TYPE(tensor) * t)
{
  free(t->data);
  free(t);
}



/* ------ Conversions ------ */


/*
 * Convert a rank 2 tensor to a NxN (square) matrix.
 */
TYPE (gsl_matrix) *
FUNCTION (tensor, 2matrix) (TYPE (tensor) * t)
{
  size_t n = t->dimension;
  TYPE(gsl_matrix) * m;

  if (t->rank != 2)
    GSL_ERROR_NULL("tensor of rank != 2", GSL_EINVAL);


  m = (TYPE (gsl_matrix) *) malloc (sizeof (TYPE (gsl_matrix)));
  if (m == 0)
    GSL_ERROR_VAL ("failed to allocate space for matrix struct",
                   GSL_ENOMEM, 0);

#if defined(BASE_COMPLEX_DOUBLE)
  m->data = (double *) t->data;
#else
  m->data = t->data;
#endif
  m->size1 = n;
  m->size2 = n;
  m->tda = n;
  m->block = NULL;  /* note that this is no problem because owner=0 */
  m->owner = 0;

  return m;
}


/*
 * Convert a rank 1 tensor to a vector.
 */
TYPE (gsl_vector) *
FUNCTION (tensor, 2vector) (TYPE (tensor) * t)
{
  size_t n = t->dimension;
  TYPE(gsl_vector) * v;

  if (t->rank != 1)
    GSL_ERROR_NULL("tensor of rank != 1", GSL_EINVAL);


  v = (TYPE (gsl_vector) *) malloc (sizeof (TYPE (gsl_vector)));
  if (v == 0)
    GSL_ERROR_VAL ("failed to allocate space for vector struct",
                   GSL_ENOMEM, 0);

#if defined(BASE_COMPLEX_DOUBLE)
  v->data = (double *) t->data;
#else
  v->data = t->data;
#endif
  v->size = n;
  v->stride = 1;
  v->block = NULL;  /* note that this is no problem because owner=0 */
  v->owner = 0;

  return v;
}




/* ------ Operations ------ */


/*
 * t = 0  (all elements = 0)
 */
void
FUNCTION(tensor, set_zero) (TYPE(tensor) * t)
{
  ATOMIC * const data = t->data;
  size_t i, n;

  n = t->size;

  for (i = 0; i < n; i++)
    {
      *(BASE *) (data + i) = 0;
    }
}


/*
 * t_ijk = x  (for i,j,k = 0,1,...,dimension-1)
 */
void
FUNCTION(tensor, set_all) (TYPE(tensor) * t, BASE x)
{
  size_t i, n;
  ATOMIC * const data = t->data;

  n = t->size;
  for (i = 0; i < n; i ++)
    *(BASE *) (data + i) = x;
}
