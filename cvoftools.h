//---------------------------------------------------------------------//
//---------------------------------------------------------------------//
//                           cvoftools.h                               //
//---------------------------------------------------------------------//
//            Copyright (C) 2020 J. Lopez and J. Hernandez             //
//---------------------------------------------------------------------//
//      Declaration of the subroutines of the VOFTools library         //
//      to be used in C programs                                       //
//---------------------------------------------------------------------//
// This file is part of VOFTools.                                      //
//                                                                     //
// VOFTools is free software: you can redistribute it and/or           //
// modify it under the terms of the GNU Lesser General Public License  //
// as published by the Free Software Foundation, either version 3 of   //
// the License, or (at your option) any later version.                 //
//                                                                     //
// VOFTools is distributed in the hope that it will be useful,         //
// but WITHOUT ANY WARRANTY; without even the implied warranty of      //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       //
// GNU Lesser General Public License for more details.                 //
//                                                                     //
// You should have received a copy of the GNU Lesser General Public    //
// License along with VOFTools.  If not, see                           //
// <http://www.gnu.org/licenses/>.                                     //
//---------------------------------------------------------------------//

// Declaration of 3D subroutines

extern void enforv3d_(double *c,int *ipv,int *nipv,int *ntp,int *nts,
		      int *ntv,double *v,double *vt,double *vertp,
		      double *xnc, double *xns,double *ync,double *yns,
		      double *znc,double *zns);

extern void enforv3dsz_(double *c,double *dx,double *dy,double *dz,
			double *v,double *vertp,double *xnc,double *ync,
			double *znc);

extern void newpol3d_(int *ia,int *ipia0,int *ipia1,int *ipv0,
		      int *iscut,int *nipv0,int *ntp0,int *nts0,
		      int *ntv0,double *xnc,double *xns0,double *ync,
		      double *yns0,double *znc,double *zns0);

extern void inte3d_(double *c,int *icontn,int *icontp,int *ipv0,
		    int *nipv0,int *ntp0,int *nts0,int *ntv0,
		    double *vertp0,double *xnc,double *xns0,
		    double *ync,double *yns0,double *znc,double *zns0);

extern void toolv3d_(int *ipv,int *nipv,int *nts,double *vertp,
		     double *vt,double *xns,double *yns,double *zns);

extern void cppol3d_(double *cs,double *cs0,int *ipv,int *ipv0,
		     int *nipv,int *nipv0,int *ntp,int *ntp0,
		     int *nts,int *nts0,int *ntv,int *ntv0,
		     double *verti,double *verti0,double *xns,
		     double *xns0,double *yns,double *yns0,double *zns,
		     double *zns0);

extern void restore3d_(double *cs,int *ipv,int *nipv,int *ntp,int *nts,
		       int *ntv,double *vertp,double *xns,double *yns,
		       double *zns);

extern void eqsol3d_(double *c0,double *c1,double *c2,double *c3,
		     double *cmin,double *cmax,double *csol);

extern void newton3d_(double *a,double *b,double *c,double *d,
		      double *cmin,double *cmax,double *csol,int *isol);

extern void dist3d_(double *d,int *n,double *x,double *y,double *z,
		    double *xp,double *yp,double *zp);

extern double func3d_(double *x,double *y,double *z);

extern void initf3d_(double (*func3d_)(double *,double *,double *z),int *ipv,
		    int *nc,int *nipv,int *ntp,int *nts,int *ntv,
		    double *tol,double *vertp,double *vf,double *xns,
		    double *yns,double *zns);

extern void polout3d_(int *ifile,int *ipv,int *nipv,int *ntp,int *nts,
		      double *vertp);

// Declaration of 2D subroutines

extern void enforv2d_(double *c,int *ipv,int *ntp,int *ntv,double *v,
		      double *vt,double *vertp,double *xnc,double *ync);

extern void enforv2dsz_(double *c,double *dx,double *dy,double *v,
			double *vertp,double *xnc,double *ync);

extern void newpol2d_(int *ia,int *ipia0,int *ipia1,int *ipv0,int *ntp0,
		      int *ntv0,double *vertp0,double *xncut,double *yncut);

extern void inte2d_(double *c,int *icontn,int *icontp,int *ipv0,int *ntp0,
		    int *ntv0,double *vertp0,double *xnc,double *ync);

extern void toolv2d_(int *ipv,int *ntv,double *vertp,double *vol);

extern void cppol2d_(int *ipv,int *ipv1,int *ntp,int *ntp1,int *ntv,
		     int *ntv1,double *vertp,double *vertp1);

extern void restore2d_(int *ipv,int *ntp,int *ntv,double *vertp);

extern void dist2d_(double *d,double *x,double *y,double *xp,double *yp);

extern double func2d_(double *x,double *y);

extern void initf2d_(double (*func2d_)(double *,double *),int *ipv,int *nc,
		     int *ntp,int *ntv,double *tol,double *vertp,double *vf);

extern void polout2d_(int *ifile,int *ipv,int *ntv,double *vertp);
