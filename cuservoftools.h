//---------------------------------------------------------------------//
//---------------------------------------------------------------------//
//                          cuservoftools.h                            //
//---------------------------------------------------------------------//
//            Copyright (C) 2018 J. Lopez and J. Hernandez             //
//---------------------------------------------------------------------//
//      Declaration of the subroutines of the VOFTools user tools to   //
//      be used in C programs                                          //
//---------------------------------------------------------------------//
// This file is part of VOFTools.                                      //
//                                                                     //
// VOFTools is free software: you can redistribute it and/or           //
// modify it under the terms of the GNU General Public License as      //
// published by the Free Software Foundation, either version 3 of      //
// the License, or (at your option) any later version.                 //
//                                                                     //
// VOFTools is distributed in the hope that it will be useful,         //
// but WITHOUT ANY WARRANTY; without even the implied warranty of      //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       //
// GNU General Public License for more details.                        //
//                                                                     //
// You should have received a copy of the GNU General Public License   //
// along with VOFTools. If not, see <http://www.gnu.org/licenses/>.    //
//---------------------------------------------------------------------//

// Declaration of subroutines

extern void vofvardef_(double *f,int *icelltype,int *ishape,int *nc,
		    double *tol,double *xnc,double *xp,double *ync,
		    double *yp,double *znc,double *zp);

extern double func2d1_(double *x,double *y);

extern double func2d2_(double *x,double *y);

extern double func3d1_(double *x,double *y,double *z);

extern double func3d2_(double *x,double *y,double *z);

