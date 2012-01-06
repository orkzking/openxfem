/************************************************************************************ 

   Copyright (C) 2005
   Stephane BORDAS, Cyrille DUNANT, Vinh Phu NGUYEN, Quang Tri TRUONG, Ravindra DUDDU

   This file is part of the XFEM C++ Library (OpenXFEM++) written 
   and maintained by above authors.

   This program is free software; you can redistribute it and/or modify it.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
	Any feedback is welcome. Emails : nvinhphu@gmail.com, ...

***************************************************************************************/

#include "mathfem.h"
#include <math.h>

#ifndef M_PI 
#define M_PI  3.14159265358979323846
#endif

// measure dependent constant
#define CUBIC_ZERO 1.0e-100

void cubic (double a, double b, double c, double d, double*r1, double *r2, double *r3,int *num)
//
//   solves cubic equation for real roots
//
//   input:
//     a,b,c,d - coefficients of equation in form:
//     ax^3 + bx^2 + cx + d = 0
//
//   output:
//     r1,r2,r3 - roots (only first num roots is valid)
//     num      - number of roots resolved
//
{
 double  aa, p, q, D, pq, u, v, phi;
 double help;
 
 if(fabs(a) < CUBIC_ZERO)
 {
  if(fabs(b) < CUBIC_ZERO)
  {
   *r1 = -d / c;
   *num = 1;
  }
  else
  {
   if((D = c * c - 4.0 * b * d) < 0.0) 
	{
    *num = 0;
    return;
   }else
	{
    *r1 = (-c + sqrt(D)) / 2.0 / b;
    *r2 = (-c - sqrt(D)) / 2.0 / b;
    *num = 2;
   }
  }
 }
 else 
 {
  aa = a;
  a = b/aa;
  b = c/aa;
  c = d/aa;
  p = b - a * a / 3.0;
  q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
  pq = p * q;
  D = q * q / 4.0 + p * p * p / 27.0;
  if(fabs(D) < CUBIC_ZERO)
  {
   if(fabs(pq) < CUBIC_ZERO)
	{
    *r1 = 0.0 - a / 3.0;
    *r2 = *r1;
    *r3 = *r1;
    *num = 3;
   }
   else 
	{
    if(q < 0.0)
     *r2 = -exp(log(-q / 2.0) / 3.0);
    else
     *r2 = exp(log(q / 2.0) / 3.0);
    *r1 = -2.0 * *r2 - a / 3.0;
    *r2 -= a / 3.0;
    *num = 2;
   }
  }
  else
  {
   if(D > 0.0){
    u = -q / 2.0 + sqrt(D);
    v = -q / 2.0 - sqrt(D);
    if(u < 0.0)
     u = -exp(log(-u) / 3.0);
    else
     u = exp(log(u) / 3.0);
    if(v < 0.0)
     v = -exp(log(-v) / 3.0);
    else
     v = exp(log(v) / 3.0);
    *r1 = u + v - a / 3.0;
    *num = 1;
   }
   else 
	{
    p = sqrt(fabs(p) / 3.0);
    help = (-q / (2.0 * p * p * p));
    if (fabs (help) > 1.0) help = sgn(help); // prevent rounding errors
    
    phi = acos(help) / 3.0;
    *r1 = 2.0 * p * cos(phi) - a / 3.0;
    *r2 = -2.0 * p * cos(phi - M_PI / 3.0) - a / 3.0;
    *r3 = -2.0 * p * cos(phi + M_PI / 3.0) - a / 3.0;
    *num = 3;
   }
  }
 }
 return ;
} 
