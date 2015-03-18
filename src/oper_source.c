/* tensor/oper_source.c
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

int
FUNCTION(tensor, add) (TYPE(tensor) * a, const TYPE(tensor) * b)
{
  const unsigned int rank = a->rank;
  const size_t dimension  = a->dimension;
  size_t i, n;

  if (b->rank != rank || b->dimension != dimension)
    {
      GSL_ERROR ("tensors must have same dimensions", GSL_EBADLEN);
      return 1;
    }


  n = a->size;

  for (i = 0; i < n; i++)
    a->data[i] += b->data[i];

  return GSL_SUCCESS;
}


int
FUNCTION(tensor, sub) (TYPE(tensor) * a, const TYPE(tensor) * b)
{
  const unsigned int rank = a->rank;
  const size_t dimension  = a->dimension;
  size_t i, n;

  if (b->rank != rank || b->dimension != dimension)
    {
      GSL_ERROR ("tensors must have same dimensions", GSL_EBADLEN);
      return 1;
    }


  n = a->size;

  for (i = 0; i < n; i++)
    a->data[i] -= b->data[i];

  return GSL_SUCCESS;
}


int
FUNCTION(tensor, mul_elements) (TYPE(tensor) * a,
                                const TYPE(tensor) * b)
{
  const unsigned int rank = a->rank;
  const size_t dimension  = a->dimension;
  size_t i, n;

  if (b->rank != rank || b->dimension != dimension)
    {
      GSL_ERROR ("tensors must have same dimensions", GSL_EBADLEN);
      return 1;
    }


  n = a->size;

  for (i = 0; i < n; i++)
    a->data[i] *= b->data[i];

  return GSL_SUCCESS;
}


int
FUNCTION(tensor, div_elements) (TYPE(tensor) * a,
                                const TYPE(tensor) * b)
{
  const unsigned int rank = a->rank;
  const size_t dimension  = a->dimension;
  size_t i, n;

  if (b->rank != rank || b->dimension != dimension)
    {
      GSL_ERROR ("tensors must have same dimensions", GSL_EBADLEN);
      return 1;
    }


  n = a->size;

  for (i = 0; i < n; i++)
    a->data[i] /= b->data[i];

  return GSL_SUCCESS;
}


int
FUNCTION(tensor, scale) (TYPE(tensor) * a, const double x)
{
  size_t i, n;

  n = a->size;

  for (i = 0; i < n; i++)
    a->data[i] *= x;

  return GSL_SUCCESS;
}


int
FUNCTION(tensor, add_constant) (TYPE(tensor) * a, const double x)
{
  size_t i, n;

  n = a->size;

  for (i = 0; i < n; i++)
    a->data[i] += x;

  return GSL_SUCCESS;
}


int
FUNCTION(tensor, add_diagonal) (TYPE(tensor) * a, const double x)
{
  unsigned int i, j;
  size_t * index;
  size_t position;

  index = (size_t *) malloc(a->rank * sizeof(size_t));

  for (i = 0; i < a->rank; i++)
    {
      for (j = 0; j < a->rank; j++)
        index[j] = i;

      position = index2position(a->rank, a->dimension, index);

      a->data[position] += x;
    }


  free(index);

  return GSL_SUCCESS;
}


/*
 * Tensorial product of two tensors.
 *
 * This returns
 *   C_i..jk..l = A_i..j * B_k..l
 */
TYPE(tensor) *
FUNCTION(tensor, product) (const TYPE(tensor) * a,
                           const TYPE(tensor) * b)
{
  size_t i, j;
  size_t position;

  if (a->dimension != b->dimension)
    {
      GSL_ERROR_VAL("tensors must have same underlying dimension",
                    GSL_EBADLEN, 0);
      return NULL;
    }

  TYPE(tensor) * c =
    FUNCTION(tensor, alloc) (a->rank + b->rank, a->dimension);

  position = 0;
  for (i = 0; i < a->size; i++)
    for (j = 0; j < b->size; j++)
      c->data[position++] = a->data[i] * b->data[j];

  return c;
}


/*
 * Contracts the indices i and j of a tensor.
 *
 * This function, unlike the others above, is not trivial. The
 * algorithm is based on the correspondence of the indices in the
 * tensor and the serialized corresponding position.
 */
TYPE(tensor) *
FUNCTION(tensor, contract) (const TYPE(tensor) * t_ij,
                            size_t i, size_t j)
{
  size_t dimension;
  unsigned int rank;
  size_t n;

  size_t k;
  size_t step;
  size_t * indices;
  size_t pos;
  size_t pos_00;
  ATOMIC sum;

  dimension = t_ij->dimension;  /* to write less later */
  rank = t_ij->rank;

  if (i >= rank || j >= rank || i == j)
    {
      GSL_ERROR_VAL("bad indices to contract tensor", GSL_EINVAL, 0);
      return NULL;
    }

  /* Create a new tensor with the appropriate rank */
  TYPE(tensor) * t_ii = FUNCTION(tensor, alloc) (rank - 2, dimension);

  if (t_ii == NULL)
    {
      GSL_ERROR_VAL("no memory to allocate tensor", GSL_EINVAL, 0);
      return NULL;
    }

  /* Number of elements of the new tensor */
  n = t_ii->size;


  /* Start counting the indices in the opposite direction,
   * for clarity in the algorithm (but potentially confusing!)
   * Now the last two indices, for instance, will be 1,0.
   */
  i = (rank - 1) - i;
  j = (rank - 1) - j;

  if (i > j) {         /* swap indices if necessary, so j > i */
    k = i;
    i = j;
    j = k;
  }

  step = quick_pow(dimension, i) + quick_pow(dimension, j);

  indices = (size_t *) malloc(rank * sizeof(size_t));

  for (pos = 0; pos < n; pos++)
    {
      sum = 0;

      position2index(rank - 2, dimension, pos, indices);

      vec_insert(indices, rank - 2, i, 0);
      vec_insert(indices, rank - 1, j, 0);

      pos_00 = index2position(rank, dimension, indices);

      for (k = 0; k < dimension; k++)
        sum += t_ij->data[pos_00 + step * k];

      t_ii->data[pos] = sum;
    }

  free(indices);

  return t_ii;
}
