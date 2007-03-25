/* tensor/swap_source.c
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
 * t_.i.j. -> t_.j.i.
 */
TYPE(tensor) *
FUNCTION (tensor, swap_indices) (const TYPE (tensor) * t_ij,
                                     size_t i, size_t j)
{
  size_t pos;
  size_t n = t_ij->size;
  unsigned int rank = t_ij->rank;
  size_t dimension = t_ij->dimension;
  size_t * pos_in_base;

  if (i >= rank || j >= rank || i == j)
    {
      GSL_ERROR_VAL("bad indices in swap_indices request", GSL_EINVAL, 0);
      return NULL;
    }

  /* Create a new tensor with the appropiate rank */
  TYPE(tensor) * t_ji = FUNCTION(tensor, alloc) (rank, dimension);
    
  
  /* Start counting the indices in the opposite direction,
   * for clarity in the algorithm (but potentially confusing!)
   * Now the last two indices, for instance, will be 1,0.
   */
  i = (rank - 1) - i;
  j = (rank - 1) - j;
    
  pos_in_base = (size_t *) malloc(rank * sizeof(size_t));

  for (pos = 0; pos < n; pos++)
    {
      size_t newpos;

      position2index(rank, dimension, pos, pos_in_base);
      
      vec_swap(pos_in_base, i, j);
      
      newpos = index2position(rank, dimension, pos_in_base);
      
      t_ji->data[newpos] = t_ij->data[pos];
    }

  free(pos_in_base);

  return t_ji;
}  
