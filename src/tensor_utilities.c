/* tensor/tensor_utilities.c
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

#include <config.h>
#include <gsl/gsl_errno.h>

#include "tensor_utilities.h"

/*
 * Auxiliary functions.
 */


/*
 * We should probably call a modified version of sys/pow_int.c, but
 * let's put a hack to show how it should work.
 *
 * This obviously returns x^n
 */
size_t quick_pow(size_t x, unsigned int n)
{
  double value = 1.0;

  /*
   * repeated squaring method
   */
  do {
    if (n & 1)
      value *= x;  /* for n odd */
    n >>= 1;  /* n = n/2 */
    x *= x;
  } while (n);

  return value;
}



/*
 * Changes the absolute position ("n") in the tensor viewed as an array
 * to its set of indices ("digits"). The naming of the variables comes
 * from the fact that it is equivalent to a base change of number "n"
 * to "base".
 *
 * It returns in the argument "digits" an array with digits[i] being
 * the ith digit of number n in base "base".
 *
 * The "digits" array must be allocated with enough space *before*
 * calling this function.
 */
void position2index(unsigned int n_digits, size_t base, size_t n,
                    size_t * digits)
{
  unsigned int i;

  for (i = 0; i < n_digits; i++)
    {
      digits[i] = n - base * (n / base);
      n /= base;
    }
}



/*
 * Changes the set of indices ("digits") of a tensor to its absolute
 * position when the tensor is viewed as an array. That is equivalent
 * to changing a number in base "base" to a decimal.
 *
 * This is the complementary of function position2index().
 */
size_t index2position(unsigned int n_digits, size_t base,
                      const size_t * digits)
{
  unsigned int i;
  size_t scale = 1;
  size_t n;

  n = 0;
  for (i = 0; i < n_digits; i++)
    {
      n += digits[i] * scale;
      scale *= base;
    }
  
  return n;
}


/*
 * Inserts an element (value) into an array (v) at a given position.
 *
 * The array must have enough space, that is, have more than n elements.
 */
void vec_insert(size_t * v, unsigned int n, unsigned int position,
                size_t value)
{
  unsigned int i;

  /* Shift all the elements after "position" to leave space */
  for (i = n; i > position; i--)
    v[i] = v[i-1];
  
  v[position] = value;
}


/*
 * Swaps indices i and j of vector.
 */
void vec_swap(size_t * v, unsigned int i, unsigned int j)
{
  size_t temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}
