//---------------------------------------------------------------------//
//---------------------------------------------------------------------//
//                           TEST2D in C                               //
//---------------------------------------------------------------------//
//            Copyright (C) 2020 J. Lopez and J. Hernandez             //
//---------------------------------------------------------------------//
//      Test program in C to solve the local volume enforcement        //     
//      problem, calculate the distance from a given point to the      //
//      interfacial segment and initialize the material volume         //
//      fraction in a cell                                             //
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
#include <stdio.h>
#include <stdlib.h>
#include "dimc.h"
#include "cvoftools.h"
#include "cmesh.h"
#include "cuservoftools.h"

int main()
{
  //* Original polyhedron
  double vertp[nv*2];
  int ipv[nv];
  //* working polyhedron 0
  double vertp0[nv*2];
  int ipv0[nv];
  //* working polyhedron 1
  double vertp1[nv*2];
  int ipv1[nv];
  //* working polyhedron 2
  double vertp2[nv*2];
  int ipv2[nv];
  //* Interfacial polygon
  double x[2],y[2];
  //* other.......
  int i, icelltype, icontn, icontp, ip, ishape, n, nc, ntp, ntp0, ntv, ntv0;
  int ifile;
  double c,d,dx,dy,f,tol,v,vf,vt,xnc,xp,ync,yp,znc,zp;

  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|          TEST PROGRAM IN 2D OF VOFTools           |\n");
  printf("|                                                   |\n");
  printf("|            (Version 5, January 2020)              |\n");
  printf("|                                                   |\n");
  printf("|                       by                          |\n");
  printf("|                                                   |\n");
  printf("|            J. Lopez and J. Hernandez              |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  //* Variable definition for the test cases:
  //----------------------------------------
  //ICELLTYPE=  1, square mesh
  //         =  2, hexagonal mesh      
  //         =  3, triangular mesh
  //         =  4, quadrangular mesh
  //         =  5, pentagonal mesh
  //         =  6, irregular hexagonal mesh
  //         = 11, cubic mesh
  //         = 12, general hexahedrical mesh
  //         = 13, tetrahedrical mesh
  //         = 14, dodecahedral mesh
  //         = 15, icosahedral mesh
  //         = 16, complex mesh
  //         =101, non-convex quadrangle         
  //         =102, non-convex pentagon
  //         =103, non-convex hexagon
  //         =104, non-convex stellated hexagon        
  //         =105, hollowed square
  //         =106, non-convex multi-square cell
  //         =111, non-convex pentagonal pyramid
  //         =112, non-convex cell obtained by subtracting a pyramid to a unit cube
  //         =113, stellated cube        
  //         =114, non-convex hexahedron mesh
  //         =115, stellated dodecahedron
  //         =116, stellated icosahedron
  //         =117, hollowed cube
  //         =118, drilled cube
  //         =119, zig-zag prism
  //         =120, VOFTools logo
  // ISHAPE  =  1, circle with radious 0.25 centered at (0.5,0.5)
  //         =  2, ellipse with semi-major axis 0.5, semi-minor axis 0.2
  //               and centered at (0.5,0.5)
  //         = 11, sphere with radious 0.25 centered at (0.5,0.5,0.5)
  //         = 12, torus with major radius 2/3, minor radius 1/3 and
  //               centered at (0.5,0.5,0.5)
  vofvardef_(&f,&icelltype,&ishape,&nc,&tol,&xnc,&xp,&ync,&yp,&znc,&zp);
  if(icelltype<1||icelltype>106||(icelltype>6&&icelltype<101)){
    printf("-----------------------------------------------------\n");
    printf("|---------------------------------------------------|\n");
    printf("|*********** WARNING FOR CELL SELECTION ************|\n");
    printf("|---------------------------------------------------|\n");
    printf("-----------------------------------------------------\n");
    printf("1.- Edit the vofvardef file.\n");    
    printf("2.- Choose an appropriate ICELLTYPE value (between 1 and 6 for convex cells or between 101 and 106 for non-convex cells).\n");
    exit(-1);
  }
  if(ishape<1||ishape>2){
    printf("-----------------------------------------------------\n");
    printf("|---------------------------------------------------|\n");
    printf("|**** WARNING FOR MATERIAL BODY SHAPE SELECTION ****|\n");
    printf("|---------------------------------------------------|\n");
    printf("-----------------------------------------------------\n");
    printf("1.- Edit the vofvardef file.\n");    
    printf("2.- Choose an appropriate ISHAPE value (1 for a circle or 2 for an ellipse).\n");
    exit(-1);
  }
  
  //* Convex cells:
  //
  if(icelltype==1){
  //* Square mesh
    squaremesh_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==2){
  //* Regular hexagonal mesh
    hexagomesh_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==3){
  //* Triangular mesh
    trianglemesh_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==4){
  //* Quadrangular mesh
    quadranglemesh_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==5){
  //* Pentagonal mesh
    pentagonmesh_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==6){
  //* Irregular hexagonal mesh
    hexagonmesh_(ipv,&ntp,&ntv,vertp);
  //
  //* Non-convex cells:
  //
  }
  else if(icelltype==101){
  //* Non-convex quadrangle
    ncquadranglemesh_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==102){
  //* Non-convex pentagon
    ncpentagonmesh_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==103){
  //* Non-convex hexagon
    nchexagonmesh_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==104){
  //* Non-convex stellated hexagon
    ncshexagonmesh_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==105){
  //* Hollowed square
    nchollowedsquare_(ipv,&ntp,&ntv,vertp);
  }
  else if(icelltype==106){
  //* Non-convex multi-square cell
    ncmultisquare_(ipv,&ntp,&ntv,vertp);
  }
  //* Print the original polygonal cell to the file 'pol00000.out'
  //--------------------------------------------------------------
  ifile=0;
  polout2d_(&ifile,ipv,&ntv,vertp);
  //* Calculate the area VT of the cell:
  //------------------------------------
  toolv2d_(ipv,&ntv,vertp,&vt);
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE TOOLV2D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Area of the selected cell:%f\n",vt);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* calculate the volume of liquid in the cell
  v=f*vt;
  //* Solve the local volume enforcement problem:
  //--------------------------------------------
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|        OUTPUT OF THE ENFORV2D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  enforv2d_(&c,ipv,&ntp,&ntv,&v,&vt,vertp,&xnc,&ync);
  //* For rectangular cells, like that defined in the squaremesh_ subroutine,
  //* it is more efficient to use the analytical method of Scardovelli and 
  //* Zaleski [Journal of Computational Physics, 164 (2000) 228-237], which
  //* was proposed to be used specifically for this type of cells. To use 
  //* this method, uncomment the following three lines.
  //dx=vertp[0*nv+2-1]-vertp[0*nv+1-1];    // x-side length
  //dy=vertp[1*nv+3-1]-vertp[1*nv+2-1];    // y-side length
  //enforv2dsz_(&c,&dx,&dy,&v,vertp,&xnc,&ync);
  printf("Solution of the problem:%f\n",c);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* copy the original polygon that defines the cell to the working
  //* polygon 0
  cppol2d_(ipv,ipv0,&ntp,&ntp0,&ntv,&ntv0,vertp,vertp0);
  //* intersection between the cell and the line defined as X·NC+C=0
  inte2d_(&c,&icontn,&icontp,ipv0,&ntp0,&ntv0,vertp0,&xnc,&ync);
  //* Print the truncated polygonal cell to the file 'pol00001.out'
  //---------------------------------------------------------------
  ifile=1;
  polout2d_(&ifile,ipv0,&ntv0,vertp0);
  if(icelltype>=1&&icelltype<=6){
  //* Calculate the distance from a given point P to the interfacial segment:
  //------------------------------------------------------------------------
  //* interfacial segment defined as the last edge of the truncated polygon 0
    ip=ntp0-1;
    x[0]=vertp0[0*nv+ip-1];
    y[0]=vertp0[1*nv+ip-1];
    ip=ntp0;
    x[1]=vertp0[0*nv+ip-1];
    y[1]=vertp0[1*nv+ip-1];
  //* calculate the distance
    dist2d_(&d,x,y,&xp,&yp);
    printf("-----------------------------------------------------\n");
    printf("|---------------------------------------------------|\n");
    printf("|         OUTPUT OF THE DIST2D SUBROUTINE           |\n");
    printf("|---------------------------------------------------|\n");
    printf("-----------------------------------------------------\n");
    printf("\n");
    printf("Distance from P to the interfacial segment:%f\n",d);
    printf("\n");
    printf("-----------------------------------------------------\n");
  }
  //* Initialize the material volume fraction in the cell:
  //-----------------------------------------------------      
  if(ishape==1){
    initf2d_(func2d1_,ipv,&nc,&ntp,&ntv,&tol,vertp,&vf);
  }
  else if(ishape==2){
    initf2d_(func2d2_,ipv,&nc,&ntp,&ntv,&tol,vertp,&vf);
  }
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE INITF2D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Material volume fraction in the selected cell:%f\n",vf);
  printf("\n");
  printf("-----------------------------------------------------\n");
} 

