/* tensor/oper_complex_source.c
 * 
 * Copyright (C) 2002, 2003, 2004 Jordi Burguet-Castell
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
 * Note this is work in progress, and is not yet functional. Maybe in
 * version 2?
 */


int
FUNCTION(tensor, add) (TYPE (tensor) * a,
                           const TYPE (tensor) * b)
{
  return tensor_add(a, b);
}


int
FUNCTION(tensor, sub) (TYPE (tensor) * a,
                           const TYPE (tensor) * b)
{
  return tensor_sub(a, b);
}


int
FUNCTION(tensor, mul_elements) (TYPE (tensor) * a,
                                    const TYPE (tensor) * b)
{
  const unsigned int rank = a->rank;
  const size_t dimension  = a->dimension;
  size_t i, n;

  if (b->rank != rank || b->dimension != dimension)
    {
      GSL_ERROR ("tensors must have same dimensions", GSL_EBADLEN);
      return 1;
    }

  
  n = a->size / 2;

  for (i = 0; i < n; i++)
    {
      ATOMIC ar = a->data[2*i];
      ATOMIC ai = a->data[2*i+1];

      ATOMIC br = b->data[2*i];
      ATOMIC bi = b->data[2*i+1];

      a->data[2*i]   = ar * br - ai * bi;
      a->data[2*i+1] = ar * bi + ai * br;
    }

  return GSL_SUCCESS;
}


int
FUNCTION(tensor, div_elements) (TYPE (tensor) * a,
                                    const TYPE (tensor) * b)
{
  const unsigned int rank = a->rank;
  const size_t dimension  = a->dimension;
  size_t i, n;

  if (b->rank != rank || b->dimension != dimension)
    {
      GSL_ERROR ("tensors must have same dimensions", GSL_EBADLEN);
      return 1;
    }

  
  n = a->size / 2;

  for (i = 0; i < n; i++)
    {
      ATOMIC ar = a->data[2*i];
      ATOMIC ai = a->data[2*i+1];

      ATOMIC br = b->data[2*i];
      ATOMIC bi = b->data[2*i+1];

      ATOMIC s = 1.0 / hypot(br, bi);

      ATOMIC sbr = s * br;
      ATOMIC sbi = s * bi;
      
      a->data[2*i]     = (ar * sbr + ai * sbi) * s;
      a->data[2*i + 1] = (ai * sbr - ar * sbi) * s;
    }

  return GSL_SUCCESS;
}


int FUNCTION(tensor, scale) (TYPE (tensor) * a, const BASE x)
{
  size_t i, n;

  ATOMIC xr = GSL_REAL(x);
  ATOMIC xi = GSL_IMAG(x);

  n = a->size / 2;

  for (i = 0; i < n; i++)
    {
      ATOMIC ar = a->data[2*i];
      ATOMIC ai = a->data[2*i + 1];
      
      a->data[2*i]     = ar * xr - ai * xi;
      a->data[2*i + 1] = ar * xi + ai * xr;
    }
  
  return GSL_SUCCESS;
}


int FUNCTION(tensor, add_constant) (TYPE (tensor) * a, const BASE x)
{
  size_t i, n;

  ATOMIC xr = GSL_REAL(x);
  ATOMIC xi = GSL_IMAG(x);

  n = a->size / 2;

  for (i = 0; i < n; i++)
    {
      a->data[2*i]     += xr;
      a->data[2*i + 1] += xi;
    }

  return GSL_SUCCESS;
}


int FUNCTION(tensor, add_diagonal) (TYPE (tensor) * a, const BASE x)
{
  int i, j;
  size_t * index;
  size_t position;

  index = (size_t *) malloc(a->rank * sizeof(size_t));


  for (i = 0; i < a->rank; i++)
    {
      for (j = 0; j < a->rank; j++)
        index[j] = i;

      position = index2position(a->rank, a->dimension, index);

      a->data[position]     += GSL_REAL(x);
      a->data[position + 1] += GSL_IMAG(x);
    }


  free(index);

  return GSL_SUCCESS;
}
