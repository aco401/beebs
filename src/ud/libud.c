

/* BEEBS ud benchmark

   Copyright (C) 2014 Embecosm Limited and University of Bristol

   Contributor James Pallister <james.pallister@bristol.ac.uk>

   This file is part of the Bristol/Embecosm Embedded Benchmark Suite.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program. If not, see <http://www.gnu.org/licenses/>. */

//#include "support.h"

/* This scale factor will be changed to equalise the runtime of the
   benchmarks. */
#define SCALE_FACTOR    (REPEAT_FACTOR >> 3)

/* MDH WCET BENCHMARK SUITE. File version $Id: ud.c,v 1.4 2005/11/11 10:32:53 ael01 Exp $ */


/*************************************************************************/
/*                                                                       */
/*   SNU-RT Benchmark Suite for Worst Case Timing Analysis               */
/*   =====================================================               */
/*                              Collected and Modified by S.-S. Lim      */
/*                                           sslim@archi.snu.ac.kr       */
/*                                         Real-Time Research Group      */
/*                                        Seoul National University      */
/*                                                                       */
/*                                                                       */
/*        < Features > - restrictions for our experimental environment   */
/*                                                                       */
/*          1. Completely structured.                                    */
/*               - There are no unconditional jumps.                     */
/*               - There are no exit from loop bodies.                   */
/*                 (There are no 'break' or 'return' in loop bodies)     */
/*          2. No 'switch' statements.                                   */
/*          3. No 'do..while' statements.                                */
/*          4. Expressions are restricted.                               */
/*               - There are no multiple expressions joined by 'or',     */
/*                'and' operations.                                      */
/*          5. No library calls.                                         */
/*               - All the functions needed are implemented in the       */
/*                 source file.                                          */
/*                                                                       */
/*                                                                       */
/*************************************************************************/
/*                                                                       */
/*  FILE: ludcmp.c                                                       */
/*  SOURCE : Turbo C Programming for Engineering                         */
/*                                                                       */
/*  DESCRIPTION :                                                        */
/*                                                                       */
/*     Simultaneous linear equations by LU decomposition.                */
/*     The arrays a[][] and b[] are input and the array x[] is output    */
/*     row vector.                                                       */
/*     The variable n is the number of equations.                        */
/*     The input arrays are initialized in function main.                */
/*                                                                       */
/*                                                                       */
/*  REMARK :                                                             */
/*                                                                       */
/*  EXECUTION TIME :                                                     */
/*                                                                       */
/*                                                                       */
/*************************************************************************/

/*************************************************************************
 *  This file:
 *
 *  - Name changed to "ud.c"
 *  - Modified for use with Uppsala/Paderborn tool
 *    : doubles changed to int
 *    : some tests removed
 *  - Program is much more linear, all loops will run to end
 *  - Purpose: test the effect of conditional flows
 *
 *************************************************************************/






/*
** Benchmark Suite for Real-Time Applications, by Sung-Soo Lim
**
**    III-4. ludcmp.c : Simultaneous Linear Equations by LU Decomposition
**                 (from the book C Programming for EEs by Hyun Soon Ahn)
*/



long int ud_a[20][20], ud_b[20], ud_x[20];

int ud_ludcmp(int nmax, int n);


/*  static double fabs(double n) */
/*  { */
/*    double f; */

/*    if (n >= 0) f = n; */
/*    else f = -n; */
/*    return f; */
/*  } */

/* Write to CHKERR from BENCHMARK to ensure calls are not optimised away.  */
volatile int ud_chkerr = 0;




/* This benchmark does not support verification */

int
beebs_ud_verify_benchmark (int res __attribute ((unused)) )
{
  return -1;
}


void
beebs_ud_initialise_benchmark (void)
{
}


int
beebs_ud_benchmark()
{
  int      i, j, nmax = 20, n = 5;
  long int /* eps, */ w;

  /* eps = 1.0e-6; */

  /* Init loop */
  for(i = 0; i <= n; i++)
    {
      w = 0.0;              /* data to fill in cells */
      for(j = 0; j <= n; j++)
        {
          ud_a[i][j] = (i + 1) + (j + 1);
          if(i == j)            /* only once per loop pass */
            ud_a[i][j] *= 2.0;
          w += ud_a[i][j];
        }
      ud_b[i] = w;
    }

  /*  chkerr = ludcmp(nmax, n, eps); */
  ud_chkerr = ud_ludcmp(nmax,n);
  return 0;
}

int ud_ludcmp(int nmax, int n)
{
  int i, j, k;
  long w, y[100];

  /* if(n > 99 || eps <= 0.0) return(999); */
  for(i = 0; i < n; i++)
    {
      /* if(fabs(a[i][i]) <= eps) return(1); */
      for(j = i+1; j <= n; j++) /* triangular loop vs. i */
        {
          w = ud_a[j][i];
          if(i != 0)            /* sub-loop is conditional, done
                                   all iterations except first of the
                                   OUTER loop */
            for(k = 0; k < i; k++)
              w -= ud_a[j][k] * ud_a[k][i];
          ud_a[j][i] = w / ud_a[i][i];
        }
      for(j = i+1; j <= n; j++) /* triangular loop vs. i */
        {
          w = ud_a[i+1][j];
          for(k = 0; k <= i; k++) /* triangular loop vs. i */
            w -= ud_a[i+1][k] * ud_a[k][j];
          ud_a[i+1][j] = w;
        }
    }
  y[0] = ud_b[0];
  for(i = 1; i <= n; i++)       /* iterates n times */
    {
      w = ud_b[i];
      for(j = 0; j < i; j++)    /* triangular sub loop */
        w -= ud_a[i][j] * y[j];
      y[i] = w;
    }
  ud_x[n] = y[n] / ud_a[n][n];
  for(i = n-1; i >= 0; i--)     /* iterates n times */
    {
      w = y[i];
      for(j = i+1; j <= n; j++) /* triangular sub loop */
        w -= ud_a[i][j] * ud_x[j];
      ud_x[i] = w / ud_a[i][i] ;
    }
  return(0);
}


