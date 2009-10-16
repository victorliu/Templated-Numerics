//*************************************************************************
//              2D and 3D determinant sign computation                     
//*************************************************************************

//*************************************************************************
// ABDPY.h
//*************************************************************************

//*************************************************************************
// Author : Olivier Devillers                                              
// Olivier.Devillers@sophia.inria.fr                                       
// http:/www.inria.fr:/prisme/personnel/devillers/anglais/determinant.html
//*************************************************************************

//*************************************************************************
//              Copyright (c) 1995  by  INRIA Prisme Project               
//                  BP 93 06902 Sophia Antipolis Cedex, France.            
//                           All rights reserved                           
//*************************************************************************

//*************************************************************************
// This code implement ideas presented in references given below
// 
// The given procedures compute the exact sign of a 2x2 or 3x3 determinant
// whose entries are 52 bits integers in 2d case or 50 bits integers in
// 3x3 case. The entries must be stored in variables of type "double".
// 
// The code use some feature of IEEE standard 754 arithmetic
// 
// 
// 
// REFERENCES :
// 
// @techreport{abdpy-esdus-94
// , author =      "F. Avnaim and J.-D. Boissonnat and O. Devillers and F. Preparata and M. Yvinec"
// , title =       "Evaluating signs of determinants using single-precision arithmetic"
// , type =        "Research {Report}"
// , number =      2306
// , institution = "INRIA"
// , address =     "BP93, 06902 Sophia-Antipolis, France"
// , year =        1994
// , url = "wais:/zenon.inria.fr:210/ra-mime-zenon-inria-fr?RR2306"
// , abstract =    "We propose a method to evaluate signs of  $2\times 2$ and $3\times 3$ determinants with $b$-bit integer entries using only $b$ and $(b+1)$-bit arithmetic respectively. This algorithm has numerous applications in geometric computation and provides a general and practical approach to robustness.  The algorithm has been implemented and experimental results show that it slows down the computing time by only a small factor  with respect to floating-point calculation."
// }
// 
// @inproceedings{bddp-spsr-92
// , author =      "F. Avnaim and J.-D. Boissonnat and O. Devillers and F. Preparata and M. Yvinec"
// , title =       "Evaluating signs of determinants using single-precision arithmetic"
// , booktitle =   "Proc. 11th Annu. ACM Sympos. Comput. Geom."
// , year =        1995
// , pages =       ""
// }
// 
//*************************************************************************/






int ABDPY__not_lazy_det2x2         (double, double, double, double);
int ABDPY_det2x2                   (double, double, double, double);
int ABDPY_secure_det2x2            (double, double, double, double);
//   These three procedures compute the sign of a 2x2 determinant.
//
//  PRECONDITION :
//      the double must contain integer values between -2^53 and 2^53
//   
//  the return code is      -1      if the determinant is negative
//  the return code is       0      if the determinant is null
//  the return code is       1      if the determinant is positive
//
//  procedure ABDPY_det2x2 and ABDPY_secure_det2x2 use some filter
//  to avoid expensive computation when the rounded value of the determinant
//  is big enough to decide the sign.
//  ABDPY__not_lazy_det2x2 does not use that, and does not use specifically
//  IEEE standard
//
//  ABDPY_secure_det2x2 verify the precondition
//  the return code is       2      if one entry is too big
//  the return code is       3      if one entry is not an integer

int ABDPY__not_lazy_det3x3
	   (double, double, double, double, double, double, double, double, double);
int ABDPY_det3x3
	   (double, double, double, double, double, double, double, double, double);
int ABDPY_secure_det3x3
	   (double, double, double, double, double, double, double, double, double);
//  These three procedures compute the sign of a 3x3 determinant.
//
//  PRECONDITION :
//      the double must contain integer values between -2^51 and 2^51
//   
//  the return code is      -1      if the determinant is negative
//  the return code is       0      if the determinant is null
//  the return code is       1      if the determinant is positive
//
//  procedure ABDPY_det3x3 and ABDPY_secure_det3x3 use some filter
//  to avoid expensive computation when the rounded value of the determinant
//  is big enough to decide the sign.
//  ABDPY__not_lazy_det3x3 does not use that, and does not use specifically
//  IEEE standard
//
//  ABDPY_secure_det3x3 verify the precondition
//  the return code is       2      if one entry is too big
//  the return code is       3      if one entry is not an integer


