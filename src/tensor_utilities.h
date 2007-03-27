/* tensor/tensor_utilities.h
 *
 * Copyright (C) 2002, 2003, 2004, 2007 Jordi Burguet Castell
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

size_t quick_pow(size_t x, unsigned int n);

void position2index(unsigned int n_digits, size_t base, size_t n,
                    size_t * digits);

size_t index2position(unsigned int n_digits, size_t base,
                      const size_t * digits);

void vec_insert(size_t * v, unsigned int n, unsigned int position,
                size_t value);

void vec_swap(size_t * v, unsigned int i, unsigned int j);
