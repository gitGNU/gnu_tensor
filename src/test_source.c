/* tensor/test_source.c
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

void FUNCTION(test, func) (void);
void FUNCTION(test, trap) (void);
void FUNCTION(test, text) (void);
void FUNCTION(test, binary) (void);


void
FUNCTION(test, func) (void)
{
  size_t i, j, k;    /* tensor indices */
  size_t counter;    /* to fill a tensor */
  size_t indices[RANK];  /* to pass tensor indices */


  /*
   * Test allocation.
   */
  TYPE(tensor) * t = FUNCTION(tensor, alloc) (RANK, DIMENSION);

  gsl_test (t->data == 0,
            NAME (tensor) "_alloc returns valid pointer");
  gsl_test (t->rank != RANK,
            NAME (tensor) "_alloc returns valid rank");
  gsl_test (t->dimension != DIMENSION,
            NAME (tensor) "_alloc returns valid dimension");


  /*
   * Test set.
   */
  counter = 0;
  for (i = 0; i < DIMENSION; i++)
    {
      for (j = 0; j < DIMENSION; j++)
        {
          for (k = 0; k < DIMENSION; k++)
            {
              counter++;
              indices[0] = i;  indices[1] = j;  indices[2] = k;
              FUNCTION(tensor, set) (t, indices, (BASE) counter);
            }
        }
    }


  {
    status = 0;
    counter = 0;
    for (i = 0; i < DIMENSION; i++)
      {
        for (j = 0; j < DIMENSION; j++)
          {
            for (k = 0; k < DIMENSION; k++)
              {
                counter++;
                if (t->data[DIMENSION*DIMENSION*i + DIMENSION*j + k] !=
                    (BASE) counter)
                  status = 1;
              }
          }
      }
    
    gsl_test (status, NAME (tensor) "_set writes into array");
  }
  



  /* Test get */
  {
    status = 0;
    counter = 0;
    for (i = 0; i < DIMENSION; i++)
      {
        for (j = 0; j < DIMENSION; j++)
          {
            for (k = 0; k < DIMENSION; k++)
              {
                counter++;
                indices[0] = i;  indices[1] = j;  indices[2] = k;
                if (FUNCTION(tensor, get) (t, indices) != (BASE) counter)
                  status = 1;
              }
          }
      }

    gsl_test (status, NAME (tensor) "_get reads from array");
  }


  FUNCTION(tensor, free) (t);   /* free whatever is in t */

  /* New round of tests */

  t = FUNCTION(tensor, calloc) (RANK, DIMENSION);

  /* Fill the tensor again */
  counter = 0;
  for (i = 0; i < DIMENSION; i++)
    {
      for (j = 0; j < DIMENSION; j++)
        {
          for (k = 0; k < DIMENSION; k++)
            {
              counter++;
              indices[0] = i;  indices[1] = j;  indices[2] = k;
              FUNCTION(tensor, set) (t, indices, (BASE) counter);
            }
        }
    }


  /*
   * Test maximum and minimum.
   */
  {
    indices[0] = 0;  indices[1] = 0;  indices[2] = 0;
    BASE exp_max = FUNCTION(tensor, get) (t, indices);
    BASE exp_min = FUNCTION(tensor, get) (t, indices);
    size_t exp_imax = 0, exp_jmax = 0, exp_kmax = 0;
    size_t exp_imin = 0, exp_jmin = 0, exp_kmin = 0;

    for (i = 0; i < DIMENSION; i++)
      {
        for (j = 0; j < DIMENSION; j++)
          {
            for (k = 0; k < DIMENSION; k++)
              {
                indices[0] = i;  indices[1] = j;  indices[2] = k;
                BASE value = FUNCTION(tensor, get) (t, indices);

                if (value > exp_max) {
                  exp_max = value;
                  exp_imax = i;  exp_jmax = j;  exp_kmax = k;
                }
                if (value < exp_min) {
                  exp_min = FUNCTION(tensor, get) (t, indices);
                  exp_imin = i;  exp_jmin = j;  exp_kmin = k;
                }
              }
          }
      }

    /* Test maximum */
    {
      BASE max = FUNCTION(tensor, max) (t);

      gsl_test (max != exp_max,
                NAME(tensor) "_max returns correct maximum value");
    }

    /* Test minimum */
    {
      BASE min = FUNCTION(tensor, min) (t);
      
      gsl_test (min != exp_min,
                NAME(tensor) "_min returns correct minimum value");
    }

    /* Test minmax */
    {
      BASE min, max;
      FUNCTION(tensor, minmax) (t, &min, &max);

      gsl_test (max != exp_max,
                NAME(tensor) "_minmax returns correct maximum value");
      gsl_test (min != exp_min,
                NAME(tensor) "_minmax returns correct minimum value");
    }

    /* Test min/max index */
    {
      size_t imax[RANK];

      FUNCTION(tensor, max_index) (t, imax);

      status = 0;

      if (imax[0] != exp_imax)
        status = 1;

      if (imax[1] != exp_jmax)
        status = 1;

      if (imax[2] != exp_kmax)
        status = 1;

      gsl_test (status,
                NAME(tensor) "_max_index returns correct maximum indices");
    }

    {
      size_t imin[RANK];

      FUNCTION(tensor, min_index) (t, imin);

      status = 0;

      if (imin[0] != exp_imin)
        status = 1;

      if (imin[1] != exp_jmin)
        status = 1;

      if (imin[2] != exp_kmin)
        status = 1;

      gsl_test (status,
                NAME(tensor) "_min_index returns correct minimum indices");
    }

    {
      size_t imin[RANK];
      size_t imax[RANK];

      FUNCTION(tensor, minmax_index) (t, imin, imax);

      status = 0;

      if (imax[0] != exp_imax)
        status = 1;

      if (imax[1] != exp_jmax)
        status = 1;

      if (imax[2] != exp_kmax)
        status = 1;

      if (imin[0] != exp_imin)
        status = 1;

      if (imin[1] != exp_jmin)
        status = 1;

      if (imin[2] != exp_kmin)
        status = 1;
      
      gsl_test (status,
                NAME(tensor) "_minmax_index returns correct indices");
    }
  }


  /*
   * Operations.
   */
  {
    TYPE (tensor) * a = FUNCTION(tensor, calloc) (RANK, DIMENSION);
    TYPE (tensor) * b = FUNCTION(tensor, calloc) (RANK, DIMENSION);

    for (i = 0; i < DIMENSION; i++)
      {
        for (j = 0; j < DIMENSION; j++)
          {
            for (k = 0; k < DIMENSION; k++)
              {
                indices[0] = i;  indices[1] = j;  indices[2] = k;
                FUNCTION(tensor, set) (a, indices,
                                           (BASE) (3 + i + 5 * j + 2 * k));
                FUNCTION(tensor, set) (b, indices,
                                           (BASE) (3 + 2 * i + 4 * j + k));
              }
          }
      }
    
    /* Addition */
    FUNCTION(tensor, memcpy) (t, a);
    FUNCTION(tensor, add) (t, b);

    {
      status = 0;

      for (i = 0; i < DIMENSION; i++)
        {
          for (j = 0; j < DIMENSION; j++)
            {
              for (k = 0; k < DIMENSION; k++)
                {
                  indices[0] = i;  indices[1] = j;  indices[2] = k;
                  BASE r = FUNCTION(tensor, get) (t, indices);
                  BASE x = FUNCTION(tensor, get) (a, indices);
                  BASE y = FUNCTION(tensor, get) (b, indices);
                  BASE z = x + y;
                  if (r != z)
                    status = 1;
                }
            }
        }
      
      gsl_test(status, NAME(tensor) "_add tensor addition");
    }

    /* Subtraction */
    FUNCTION(tensor, memcpy) (t, a);
    FUNCTION(tensor, sub) (t, b);
    
    {
      status = 0;

      for (i = 0; i < DIMENSION; i++)
        {
          for (j = 0; j < DIMENSION; j++)
            {
              for (k = 0; k < DIMENSION; k++)
                {
                  indices[0] = i;  indices[1] = j;  indices[2] = k;
                  BASE r = FUNCTION(tensor, get) (t, indices);
                  BASE x = FUNCTION(tensor, get) (a, indices);
                  BASE y = FUNCTION(tensor, get) (b, indices);
                  BASE z = x - y;
                  if (r != z)
                    status = 1;
                }
            }
        }
      
      gsl_test(status, NAME(tensor) "_sub tensor subtraction");
    }

    /* Element multiplication */
    FUNCTION(tensor, memcpy) (t, a);
    FUNCTION(tensor, mul_elements) (t, b);

    {
      status = 0;

      for (i = 0; i < DIMENSION; i++)
        {
          for (j = 0; j < DIMENSION; j++)
            {
              for (k = 0; k < DIMENSION; k++)
                {
                  indices[0] = i;  indices[1] = j;  indices[2] = k;
                  BASE r = FUNCTION(tensor, get) (t, indices);
                  BASE x = FUNCTION(tensor, get) (a, indices);
                  BASE y = FUNCTION(tensor, get) (b, indices);
                  BASE z = x * y;
                  if (r != z)
                    status = 1;
                }
            }
        }
      
      gsl_test(status,
               NAME(tensor) "_mul_elements elements multiplication");
    }

    /* Element division */
    FUNCTION(tensor, memcpy) (t, a);
    FUNCTION(tensor, div_elements) (t, b);

    {
      status = 0;

      for (i = 0; i < DIMENSION; i++)
        {
          for (j = 0; j < DIMENSION; j++)
            {
              for (k = 0; k < DIMENSION; k++)
                {
                  indices[0] = i;  indices[1] = j;  indices[2] = k;
                  BASE r = FUNCTION(tensor, get) (t, indices);
                  BASE x = FUNCTION(tensor, get) (a, indices);
                  BASE y = FUNCTION(tensor, get) (b, indices);
                  BASE z = x / y;
                  if (fabs(r - z) > 2 * GSL_FLT_EPSILON * fabs(z))
                      status = 1;
                }
            }
        }
      
      gsl_test(status, NAME(tensor) "_div_elements elements division");
    }

    /* Tensor product */
    {
      size_t l, m, n;
      size_t indices_c[6];

      TYPE(tensor) * c = FUNCTION(tensor, product) (a, b);

      status = 0;

      for (i = 0; i < DIMENSION; i++)
        for (j = 0; j < DIMENSION; j++)
          for (k = 0; k < DIMENSION; k++)
            for (l = 0; l < DIMENSION; l++)
              for (m = 0; m < DIMENSION; m++)
                for (n = 0; n < DIMENSION; n++)
                  {
                    indices_c[0] = i;  indices_c[1] = j;
                    indices_c[2] = k;  indices_c[3] = l;
                    indices_c[4] = m;  indices_c[5] = n;
                    BASE r = FUNCTION(tensor, get) (c, indices_c);
                    indices[0] = i;  indices[1] = j;  indices[2] = k;
                    BASE x = FUNCTION(tensor, get) (a, indices);
                    indices[0] = l;  indices[1] = m;  indices[2] = n;
                    BASE y = FUNCTION(tensor, get) (b, indices);
                    BASE z = x * y;
                    if (r != z)
                      status = 1;
                  }
      
      gsl_test(status, NAME(tensor) "_product tensorial product");

      FUNCTION(tensor, free) (c);
    }

    /* Index contraction */
    {
      TYPE(tensor) * tt = FUNCTION(tensor, contract) (a, 0, 1);

      gsl_test(tt->data == 0,
               NAME(tensor) "_contract returns valid pointer");
      gsl_test(tt->rank != RANK-2,
               NAME(tensor) "_contract returns valid rank");
      gsl_test(tt->dimension != DIMENSION,
               NAME (tensor) "_contract returns valid dimension");


      FUNCTION(tensor, free) (tt);
    }

    /* Swap indices */
    {
      TYPE(tensor) * a_102 = FUNCTION(tensor, swap_indices) (a, 0, 1);
      TYPE(tensor) * a_210 = FUNCTION(tensor, swap_indices) (a, 0, 2);
      TYPE(tensor) * a_021 = FUNCTION(tensor, swap_indices) (a, 1, 2);

      status = 0;
      for (i = 0; i < DIMENSION; i++)
        for (j = 0; j < DIMENSION; j++)
          for (k = 0; k < DIMENSION; k++)
            {
              indices[0] = i;  indices[1] = j;  indices[2] = k;
              BASE x     = FUNCTION(tensor, get) (a,     indices);
              
              indices[0] = j;  indices[1] = i;  indices[2] = k;
              BASE x_102 = FUNCTION(tensor, get) (a_102, indices);
              
              indices[0] = k;  indices[1] = j;  indices[2] = i;
              BASE x_210 = FUNCTION(tensor, get) (a_210, indices);
              
              indices[0] = i;  indices[1] = k;  indices[2] = j;
              BASE x_021 = FUNCTION(tensor, get) (a_021, indices);
              
              if (x != x_102 || x != x_210 || x != x_021)
                status = 1;
            }
      
      gsl_test(status, NAME(tensor) "_swap_indices swaps indices");

      FUNCTION(tensor, free) (a_102);
      FUNCTION(tensor, free) (a_210);
      FUNCTION(tensor, free) (a_021);
    }


    FUNCTION(tensor, free) (a);
    FUNCTION(tensor, free) (b);
  }


  FUNCTION(tensor, free) (t);
}


#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
void
FUNCTION(test, text) (void)
{
  size_t i, j, k;    /* tensor indices */
  size_t counter;    /* to fill a tensor */
  size_t indices[RANK];  /* to pass tensor indices */

  TYPE(tensor) * t = FUNCTION(tensor, alloc) (RANK, DIMENSION);

  {
    FILE *f = fopen("test.txt", "w");

    counter = 0;
    for (i = 0; i < DIMENSION; i++)
      for (j = 0; j < DIMENSION; j++)
        for (k = 0; k < DIMENSION; k++)
          {
            counter++;
            indices[0] = i;  indices[1] = j;  indices[2] = k;
            FUNCTION(tensor, set) (t, indices, (BASE) counter);
          }
    
    FUNCTION(tensor, fprintf) (f, t, OUT_FORMAT);
    fclose(f);
  }

  {
    FILE *f = fopen("test.txt", "r");
    TYPE(tensor) * tt = FUNCTION(tensor, alloc) (RANK, DIMENSION);
    status = 0;

    FUNCTION(tensor, fscanf) (f, tt);
    
    counter = 0;
    for (i = 0; i < DIMENSION; i++)
      {
        for (j = 0; j < DIMENSION; j++)
          {
            for (k = 0; k < DIMENSION; k++)
              {
                counter++;
                if (tt->data[DIMENSION*DIMENSION*i + DIMENSION*j + k] !=
                    (BASE) counter)
                  status = 1;
              }
          }
      }

    gsl_test(status, NAME (tensor) "_fprintf and fscanf");

    fclose(f);
    FUNCTION(tensor, free) (tt);
  }

  FUNCTION(tensor, free) (t);
}
#endif


void
FUNCTION(test, binary) (void)
{
  size_t i, j, k;    /* tensor indices */
  size_t counter;    /* to fill a tensor */
  size_t indices[RANK];  /* to pass tensor indices */

  TYPE(tensor) * t = FUNCTION(tensor, calloc) (RANK, DIMENSION);

  {
    FILE *f = fopen("test.dat", "wb");

    counter = 0;
    for (i = 0; i < DIMENSION; i++)
      {
        for (j = 0; j < DIMENSION; j++)
          {
            for (k = 0; k < DIMENSION; k++)
              {
                counter++;
                indices[0] = i;  indices[1] = j;  indices[2] = k;
                FUNCTION(tensor, set) (t, indices, (BASE) counter);
              }
          }
      }
    
    FUNCTION(tensor, fwrite) (f, t);
    fclose (f);
  }


  {
    FILE *f = fopen("test.dat", "rb");
    TYPE(tensor) * tt = FUNCTION(tensor, alloc) (RANK, DIMENSION);
    status = 0;

    FUNCTION(tensor, fread) (f, tt);
    
    counter = 0;
    for (i = 0; i < DIMENSION; i++)
      {
        for (j = 0; j < DIMENSION; j++)
          {
            for (k = 0; k < DIMENSION; k++)
              {
                counter++;
                if (tt->data[DIMENSION*DIMENSION*i + DIMENSION*j + k] !=
                    (BASE) counter)
                  status = 1;
              }
          }
      }

    gsl_test(status, NAME (tensor) "_write and read");

    fclose(f);
    FUNCTION(tensor, free) (tt);
  }
  
  FUNCTION(tensor, free) (t);
}



void
FUNCTION(test, trap) (void)
{
  size_t indices[RANK];  /* to pass tensor indices */
  double x;

  TYPE(tensor) * t = FUNCTION(tensor, calloc) (RANK, DIMENSION);

  /*
   * tensor_set tests.
   */

  /* Overflow above upper bound */

  status = 0;
  indices[0] = DIMENSION+1;  indices[1] = 0;  indices[2] = 0;
  FUNCTION(tensor, set) (t, indices, (BASE) 1.2);

  gsl_test (!status,
            NAME (tensor) "_set traps 1st index above upper bound");
  

  status = 0;
  indices[0] = 0;  indices[1] = DIMENSION+1;  indices[2] = 0;
  FUNCTION(tensor, set) (t, indices, (BASE) 1.2);

  gsl_test (!status,
            NAME (tensor) "_set traps 2nd index above upper bound");

  status = 0;
  indices[0] = 0;  indices[1] = 0;  indices[2] = DIMENSION+1;
  FUNCTION(tensor, set) (t, indices, (BASE) 1.2);

  gsl_test (!status,
            NAME (tensor) "_set traps 3rd index above upper bound");
  

  /* Overflow at upper bound */

  status = 0;
  indices[0] = 0;  indices[1] = DIMENSION;  indices[2] = 0;
  FUNCTION(tensor, set) (t, indices, (BASE) 1.2);

  gsl_test (!status,
            NAME (tensor) "_set traps 2nd index at upper bound");


  /* Underflow */

  status = 0;
  indices[0] = 0;  indices[1] = -1;  indices[2] = 0;
  FUNCTION(tensor, set) (t, indices, (BASE) 1.2);

  gsl_test (!status,
            NAME (tensor) "_set traps 2nd index below lower bound");


  /*
   * tensor_get tests.
   */

  /* Overflow above upper bound */

  status = 0;
  indices[0] = DIMENSION+1;  indices[1] = 0;  indices[2] = 0;
  x = FUNCTION(tensor, get) (t, indices);

  gsl_test (!status,
            NAME (tensor) "_get traps 1st index above upper bound");
  

  status = 0;
  indices[0] = 0;  indices[1] = DIMENSION+1;  indices[2] = 0;
  x = FUNCTION(tensor, get) (t, indices);

  gsl_test (!status,
            NAME (tensor) "_get traps 2nd index above upper bound");

  status = 0;
  indices[0] = 0;  indices[1] = 0;  indices[2] = DIMENSION+1;
  x = FUNCTION(tensor, get) (t, indices);

  gsl_test (!status,
            NAME (tensor) "_get traps 3rd index above upper bound");
  

  /* Overflow at upper bound */

  status = 0;
  indices[0] = 0;  indices[1] = DIMENSION;  indices[2] = 0;
  x = FUNCTION(tensor, get) (t, indices);

  gsl_test (!status,
            NAME (tensor) "_get traps 2nd index at upper bound");


  /* Underflow */

  status = 0;
  indices[0] = 0;  indices[1] = -1;  indices[2] = 0;
  x = FUNCTION(tensor, get) (t, indices);

  gsl_test (!status,
            NAME (tensor) "_get traps 2nd index below lower bound");



  FUNCTION(tensor, free) (t);
}
