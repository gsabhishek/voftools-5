//---------------------------------------------------------------------//
//---------------------------------------------------------------------//
//                              TEST3D                                 //
//---------------------------------------------------------------------//
//            Copyright (C) 2020 J. Lopez and J. Hernandez             //
//---------------------------------------------------------------------//
//      Test program in C to solve the local volume enforcement        //
//      problem, calculate the distance from a given point to the      //
//      interfacial polygon and initialize the material volume         //
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
  double cs[ns],vertp[nv*3],xns[ns],yns[ns],zns[ns];
  int ipv[ns*nv],nipv[ns];
  //* working polyhedron 0
  double cs0[ns],vertp0[nv*3],xns0[ns],yns0[ns],zns0[ns];
  int ipv0[ns*nv],nipv0[ns];
  //* working polyhedron 1
  double cs1[ns],vertp1[nv*3],xns1[ns],yns1[ns],zns1[ns];
  int ipv1[ns*nv],nipv1[ns];
  //* working polyhedron 2
  double cs2[ns],vertp2[nv*3],xns2[ns],yns2[ns],zns2[ns];
  int ipv2[ns*nv],nipv2[ns];
  //* Interfacial polygon
  double x[nv],y[nv],z[nv];
  //* other.......
  int i, icelltype, icontn, icontp, ip, ishape, n, nc, ntp, ntp0, nts,nts0, ntv, ntv0;
  int ifile;
  double c,d,dx,dy,dz,f,tol,v,vf,vt,xnc,xp,ync,yp,znc,zp;  
  
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|          TEST PROGRAM IN 3D OF VOFTools           |\n");
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
  if(icelltype<11||icelltype>120||(icelltype>16&&icelltype<111)){
    printf("-----------------------------------------------------\n");
    printf("|---------------------------------------------------|\n");
    printf("|*********** WARNING FOR CELL SELECTION ************|\n");
    printf("|---------------------------------------------------|\n");
    printf("-----------------------------------------------------\n");
    printf("1.- Edit the vofvardef file.\n");    
    printf("2.- Choose an appropriate ICELLTYPE value (between 11 and 16 for convex cells or between 111 and 120 for non-convex cells).\n");
    exit(-1);
  }
  if(ishape<11||ishape>12){
    printf("-----------------------------------------------------\n");
    printf("|---------------------------------------------------|\n");
    printf("|**** WARNING FOR MATERIAL BODY SHAPE SELECTION ****|\n");
    printf("|---------------------------------------------------|\n");
    printf("-----------------------------------------------------\n");
    printf("1.- Edit the vofvardef file.\n");    
    printf("2.- Choose an appropriate ISHAPE value (11 for a sphere or 12 for a torus).\n");
    exit(-1);
  }

  //
  //* Convex cells:
  //
  if(icelltype==11){
  //* Cubic mesh
    cubicmesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==12){
  //* general hexahedrical mesh
    hexahemesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==13){
  //* Tetrahedrical mesh
    tetramesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==14){
  //* dodecahedral mesh
    dodecamesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==15){
  //* icosahedral mesh
    icosamesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==16){
  //* complex mesh
    complexmesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  //
  //* Non-convex cells:
  //
  }
  else if(icelltype==111){
  //* Non-convex pentagonal pyramid
  ncpentapyramid_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==112){
  //* Non-convex cell: obtained by subtracting a pyramid to a unit cube
    nccubicpyramid_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==113){
  //* Small stellated cube
    ncscubicmesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==114){
  //* Non-convex hexahedron mesh 
    nchexahemesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==115){
  //* small stellated dodecahedron
    ncdodecamesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==116){
  //* Small stellated icosahedron
    ncicosamesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==117){
  //* Hollowed cube
    nchollowedcube_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==118){
  //* Drilled cube
    drilledcube_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==119){
  //* Zig-zag prism
    zigzagmesh_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  else if(icelltype==120){
  //* VOFTools logo
    voftoolslogo_(ipv,nipv,&ntp,&nts,&ntv,vertp,xns,yns,zns);
  }
  //* Print the original polyhedral cell to the file 'pol00000.vtk'
  //--------------------------------------------------------------
  ifile=0;
  polout3d_(&ifile,ipv,nipv,&ntp,&nts,vertp);
  //* Calculate the volume VT of the cell:
  //-------------------------------------
  toolv3d_(ipv,nipv,&nts,vertp,&vt,xns,yns,zns);
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE TOOLV3D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Volume of the selected cell:%f\n",vt);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* calculate the volume of liquid in the cell
  v=f*vt;
  //* Solve the local volume enforcement problem:
  //--------------------------------------------
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|        OUTPUT OF THE ENFORV3D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  enforv3d_(&c,ipv,nipv,&ntp,&nts,&ntv,&v,&vt,vertp,&xnc,xns,&ync,yns,
	    &znc,zns);
  //* For rectangular parallelepiped cells, like that defined in the 
  //* CUBICMESH subroutine,it is more efficient to use the analytical method 
  //* of Scardovelli and Zaleski [Journal of Computational Physics, 164 
  //* (2000) 228-237], which was proposed to be used specifically for this 
  //* type of cells. To use this method, uncomment the following four lines.
  //dx=vertp[0*nv+1-1]-vertp[0*nv+5-1];     // x-side length
  //dy=vertp[1*nv+4-1]-vertp[1*nv+1-1];     // y-side length
  //dz=vertp[2*nv+1-1]-vertp[2*nv+2-1];     // z-side length
  //enforv3dsz_(&c,&dx,&dy,&dz,&v,vertp,&xnc,&ync,&znc);
  printf("Solution for c:%f\n",c);
  printf("\n");
  printf("-----------------------------------------------------\n");
  //* copy the original polyhedron that defines the cell to the working
  //* polyhedron 0
  cppol3d_(cs0,cs,ipv0,ipv,nipv0,nipv,&ntp0,&ntp,&nts0,&nts,&ntv0,
	   &ntv,vertp0,vertp,xns0,xns,yns0,yns,zns0,zns);
  //c* intersection between the cell and the plane defined as X·NC+C=0  
  inte3d_(&c,&icontn,&icontp,ipv0,nipv0,&ntp0,&nts0,&ntv0,vertp0,&xnc,
	  xns0,&ync,yns0,&znc,zns0);
  //* Print the truncated polyhedral cell to the file 'pol00001.vtk'
  //---------------------------------------------------------------
  ifile=1;
  polout3d_(&ifile,ipv0,nipv0,&ntp0,&nts0,vertp0);
  if(icelltype>=11&&icelltype<=16){
  //* Calculate the distance from a given point P to the interfacial polygon:
  //------------------------------------------------------------------------
  //* interfacial polygon defined as the last face of the truncated polyhedron 0
    n=nipv0[nts0-1];
    for(i=0; i<n; ++i){
      ip=ipv0[i*ns+nts0-1];
      x[i]=vertp0[0*nv+ip-1];
      y[i]=vertp0[1*nv+ip-1];
      z[i]=vertp0[2*nv+ip-1];
    }
  //* calculate the distance
    dist3d_(&d,&n,x,y,z,&xp,&yp,&zp);
    printf("-----------------------------------------------------\n");
    printf("|---------------------------------------------------|\n");
    printf("|         OUTPUT OF THE DIST3D SUBROUTINE           |\n");
    printf("|---------------------------------------------------|\n");
    printf("-----------------------------------------------------\n");
    printf("\n");
    printf("Distance from P to the interfacial polygon:%f\n",d);
    printf("\n");
    printf("-----------------------------------------------------\n");
  }
  //* Initialize the material volume fraction in the cell:
  //-----------------------------------------------------      
  if(ishape==11){
    initf3d_(func3d1_,ipv,&nc,nipv,&ntp,&nts,&ntv,&tol,vertp,&vf,xns,yns,zns);
  }
  else if(ishape==12){
    initf3d_(func3d2_,ipv,&nc,nipv,&ntp,&nts,&ntv,&tol,vertp,&vf,xns,yns,zns);
  }
  printf("-----------------------------------------------------\n");
  printf("|---------------------------------------------------|\n");
  printf("|         OUTPUT OF THE INITF3D SUBROUTINE          |\n");
  printf("|---------------------------------------------------|\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf("Material volume fraction in the selected cell:%f\n",vf);
  printf("\n");
  printf("-----------------------------------------------------\n");
} 
