
/* BEEBS select benchmark

   SOURCE : Numerical Recipes in C - The Second Edition

   Copyright (C) 2014 Embecosm Limited and University of Bristol

   Contributor Pierre Langlois <pierre.langlois@embecosm.com>

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
#define SCALE_FACTOR    (REPEAT_FACTOR >> 0)

#define SWAP(a,b) \
do{               \
    temp=(a);     \
    (a)=(b);      \
    (b)=temp;     \
}while(0);        \

float select_arr[20] = {
  5, 4, 10.3, 1.1, 5.7, 100, 231, 111, 49.5, 99,
  10, 150, 222.22, 101, 77, 44, 35, 20.54, 99.99, 888.88
};


float select(unsigned long k, unsigned long n)
{
	unsigned long i,ir,j,l,mid;
	float a,temp;
	int flag, flag2;

	l=1;
	ir=n;
	flag = flag2 = 0;
	while (!flag) {
		if (ir <= l+1) {
			if (ir == l+1)
			  if (select_arr[ir] < select_arr[l]) {
			    SWAP(select_arr[l],select_arr[ir])
			      }
			flag = 1;
		} else if (!flag) {
			mid=(l+ir) >> 1;
			SWAP(select_arr[mid],select_arr[l+1])
			if (select_arr[l+1] > select_arr[ir]) {
				SWAP(select_arr[l+1],select_arr[ir])
			}
			if (select_arr[l] > select_arr[ir]) {
				SWAP(select_arr[l],select_arr[ir])
			}
			if (select_arr[l+1]> select_arr[l]) {
				SWAP(select_arr[l+1],select_arr[l])
			}
			i=l+1;
			j=ir;
			a=select_arr[l];
			while (!flag2) {
				i++;
				while (select_arr[i] < a) i++;
				j--;
				while (select_arr[j] > a) j--;
				if (j < i) flag2 = 1;
				if (!flag2) SWAP(select_arr[i],select_arr[j]);

			}
			select_arr[l]=select_arr[j];
			select_arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}

	}
	return select_arr[k];
}

static int x, y;

int
beebs_select_benchmark (void)
{
  select(x, y);
  return 0;
}

void beebs_select_initialise_benchmark() {
  x = 10;
  y = 20;
}


/* This benchmark does not support verification */

int
beebs_select_verify_benchmark (int res __attribute ((unused)) )
{
  return -1;
}
