/* tensor/copy_source.c
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

#include <string.h>  /* to use memcpy() */


/*
 * Overwrites dest with the contents of tensor src.
 */
int
FUNCTION (tensor, memcpy) (TYPE (tensor) * dest, const TYPE (tensor) * src)
{
  if (dest->rank != src->rank || dest->dimension != src->dimension)
    {
      GSL_ERROR ("tensor sizes are different", GSL_EBADLEN);
    }

  memcpy(dest->data, src->data, sizeof(BASE) * src->size);

  return GSL_SUCCESS;
}


/*
 * Interchanges the values of tensors t1 and t2
 */
int
FUNCTION (tensor, swap) (TYPE (tensor) * t1, TYPE (tensor) * t2)
{
  if (t1->rank != t2->rank || t1->dimension != t2->dimension)
    {
      GSL_ERROR ("tensor sizes are different", GSL_EBADLEN);
    }

  {
    size_t i, j;
    size_t n = t1->size;

    for (i = 0; i < n; i++)
      {
        ATOMIC tmp = t1->data[i];
	  
        t1->data[i] = t2->data[i];

        t2->data[i] = tmp;
      }
  }

  return GSL_SUCCESS;
}
