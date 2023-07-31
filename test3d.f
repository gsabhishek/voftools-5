c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              TEST3D                                 c
c---------------------------------------------------------------------c
c            Copyright (C) 2020 J. Lopez and J. Hernandez             c
c---------------------------------------------------------------------c
c      Test program to solve the local volume enforcement problem,    c
c      calculate the distance from a given point to the interfacial   c
c      polygon and initialize the material volume fraction in a cell  c 
c---------------------------------------------------------------------c
c This file is part of VOFTools.                                      c
c                                                                     c
c VOFTools is free software: you can redistribute it and/or           c
c modify it under the terms of the GNU Lesser General Public License  c
c as published by the Free Software Foundation, either version 3 of   c
c the License, or (at your option) any later version.                 c
c                                                                     c
c VOFTools is distributed in the hope that it will be useful,         c
c but WITHOUT ANY WARRANTY; without even the implied warranty of      c
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       c
c GNU Lesser General Public License for more details.                 c
c                                                                     c
c You should have received a copy of the GNU Lesser General Public    c
c License along with VOFTools.  If not, see                           c
c <http://www.gnu.org/licenses/>.                                     c
c---------------------------------------------------------------------c 
      PROGRAM TEST3D
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polyhedron
      DIMENSION CS(NS),IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),YNS(NS),
     -     ZNS(NS)
c* Working polyhedron 0
      DIMENSION CS0(NS),IPV0(NS,NV),NIPV0(NS),VERTP0(NV,3),XNS0(NS),
     -     YNS0(NS),ZNS0(NS)
c* Interfacial polygon
      DIMENSION X(NV),Y(NV),Z(NV)
c* Level contour
      DIMENSION IPVLC(NS,NV),ISLC(NV),NIPVLC(NS),PHIV(NV),VERTLC(NV,3)
c* External function to define the interface shape      
      EXTERNAL FUNC3D1,FUNC3D2
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|          TEST PROGRAM IN 3D OF VOFTools           |'
      WRITE(6,*)'|                                                   |'
      WRITE(6,*)'|            (Version 5, January 2020)              |'
      WRITE(6,*)'|                                                   |'
      WRITE(6,*)'|                       by                          |'
      WRITE(6,*)'|                                                   |'
      WRITE(6,*)'|            J. Lopez and J. Hernandez              |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
c. Variable definition for the test cases:
c.ICELLTYPE=  1, square mesh
c.         =  2, hexagonal mesh      
c.         =  3, triangular mesh
c.         =  4, quadrangular mesh
c.         =  5, pentagonal mesh
c.         =  6, irregular hexagonal mesh
c.         = 11, cubic mesh
c.         = 12, general hexahedrical mesh
c.         = 13, tetrahedrical mesh
c.         = 14, dodecahedral mesh
c.         = 15, icosahedral mesh
c.         = 16, complex mesh
c.         =101, non-convex quadrangle         
c.         =102, non-convex pentagon
c.         =103, non-convex hexagon
c.         =104, non-convex stellated hexagon        
c.         =105, hollowed square
c.         =106, non-convex multi-square cell
c.         =111, non-convex pentagonal pyramid
c.         =112, non-convex cell obtained by subtracting a pyramid to a unit cube
c.         =113, stellated cube        
c.         =114, non-convex hexahedron mesh
c.         =115, stellated dodecahedron
c.         =116, stellated icosahedron
c.         =117, hollowed cube
c.         =118, drilled cube
c.         =119, zig-zag prism
c.         =120, VOFTools logo
c. ISHAPE  =  1, circle with radious 0.25 centered at (0.5,0.5)
c.         =  2, ellipse with semi-major axis 0.5, semi-minor axis 0.2
c.               and centered at (0.5,0.5)
c.         = 11, sphere with radious 0.25 centered at (0.5,0.5,0.5)
c.         = 12, torus with major radius 2/3, minor radius 1/3 and
c.               centered at (0.5,0.5,0.5)
      CALL VOFVARDEF(F,ICELLTYPE,ISHAPE,NC,TOL,XNC,XP,YNC,YP,ZNC,ZP)
      IF(ICELLTYPE.LT.11.OR.ICELLTYPE.GT.120.OR.(ICELLTYPE.GT.16.AND.
     -     ICELLTYPE.LT.111)) THEN
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|*********** WARNING FOR CELL SELECTION ************|'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,*) "1.- Edit the vofvardef file."
         WRITE(6,*) "2.- Choose an appropriate ICELLTYPE value (",
     -               "between 11 and 16 for convex cells or between ", 
     -               "111 and 120 for non-convex cells)."
         STOP
      END IF
      IF(ISHAPE.LT.11.OR.ISHAPE.GT.12) THEN
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|**** WARNING FOR MATERIAL BODY SHAPE SELECTION ****|'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,*) "1.- Edit the vofvardef file."
         WRITE(6,*) "2.- Choose an appropriate ISHAPE value (11 for a ",
     -               "sphere or 12 for a torus)."
         STOP
      END IF
      IF(ICELLTYPE.EQ.11) THEN
c* Cubic mesh
         CALL CUBICMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.12) THEN
c* General hexahedrical mesh
         CALL HEXAHEMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.13) THEN
c* Tetrahedrical mesh
         CALL TETRAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.14) THEN
c* Dodecahedral mesh
         CALL DODECAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.15) THEN
* Icosahedral mesh
         CALL ICOSAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.16) THEN
c* Complex mesh
         CALL COMPLEXMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)

c* Non-convex cells:
      ELSEIF(ICELLTYPE.EQ.111) THEN
c* Non-convex pentagonal pyramid
         CALL NCPENTAPYRAMID(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.112) THEN
c* Non-convex cell: obtained by subtracting a pyramid to a unit cube
         CALL NCCUBICPYRAMID(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.113) THEN
c* Small stellated cube
         CALL NCSCUBICMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.114) THEN
c* Non-convex hexahedron mesh 
         CALL NCHEXAHEMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.115) THEN
c* Small stellated dodecahedron
         CALL NCDODECAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.116) THEN
c* Small stellated icosahedron
         CALL NCICOSAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.117) THEN
c* Hollowed cube
         CALL NCHOLLOWEDCUBE(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.118) THEN
c* Drilled cube
         CALL DRILLEDCUBE(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.119) THEN
c* Zig-zag prism
         CALL ZIGZAGMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      ELSEIF(ICELLTYPE.EQ.120) THEN
c* VOFTools logo
         CALL VOFTOOLSLOGO(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      END IF
c* Print the original polyhedral cell to the file 'pol00000.vtk'
      IFILE=0
      CALL POLOUT3D(IFILE,IPV,NIPV,NTP,NTS,VERTP)
c* calculate the volume VT of the cell:
c-------------------------------------
      CALL TOOLV3D(IPV,NIPV,NTS,VERTP,VT,XNS,YNS,ZNS)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE TOOLV3D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Volume of the selected cell:',VT
      WRITE(6,*)''
c* calculate the volume of liquid in the cell
      V=F*VT
c* Solve the local volume enforcement problem:
c--------------------------------------------
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|        OUTPUT OF THE ENFORV3D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      CALL ENFORV3D(C,IPV,NIPV,NTP,NTS,NTV,V,VT,VERTP,XNC,XNS,YNC,YNS,
     -     ZNC,ZNS)
c.. For rectangular parallelepiped cells, like that defined in the 
c.. CUBICMESH subroutine,it is more efficient to use the analytical method 
c.. of Scardovelli and Zaleski [Journal of Computational Physics, 164 
c.. (2000) 228-237], which was proposed to be used specifically for this 
c.. type of cells. To use this method, uncomment the following four lines.
c      DX=DABS(VERTP(1,1)-VERTP(5,1))     !x-side length
c      DY=DABS(VERTP(1,2)-VERTP(4,2))     !y-side length
c      DZ=DABS(VERTP(2,3)-VERTP(1,3))     !z-side length
c      CALL ENFORV3DSZ(C,DX,DY,DZ,V,VERTP,XNC,YNC,ZNC)
      WRITE(6,*)'Solution for C:',C
      WRITE(6,*)''
      WRITE(6,*)'-----------------------------------------------------'
c* copy the original polyhedron that defines the cell to the working
c* polyhedron 0
      CALL CPPOL3D(CS0,CS,IPV0,IPV,NIPV0,NIPV,NTP0,NTP,NTS0,NTS,NTV0,
     -     NTV,VERTP0,VERTP,XNS0,XNS,YNS0,YNS,ZNS0,ZNS)
c* intersection between the cell and the plane defined as X·NC+C=0  
      CALL INTE3D(C,ICONTN,ICONTP,IPV0,NIPV0,NTP0,NTS0,NTV0,VERTP0,
     -     XNC,XNS0,YNC,YNS0,ZNC,ZNS0)
c* Print the truncated polyhedral cell to the file 'pol00001.vtk'
      IFILE=1
      CALL POLOUT3D(IFILE,IPV0,NIPV0,NTP0,NTS0,VERTP0)
      IF(ICELLTYPE.GE.11.AND.ICELLTYPE.LE.16) THEN
c* Calculate the distance from a given point P to the interfacial polygon:
c------------------------------------------------------------------------
c* interfacial polygon defined as the last face of the truncated polyhedron 0
         N=NIPV0(NTS0)
         DO I=1,N
            IP=IPV0(NTS0,I)
            X(I)=VERTP0(IP,1)
            Y(I)=VERTP0(IP,2)
            Z(I)=VERTP0(IP,3)
         END DO
c* calculate the distance
         CALL DIST3D(D,N,X,Y,Z,XP,YP,ZP)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE DIST3D SUBROUTINE           |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Distance from P to the interfacial polygon:',D
      WRITE(6,*)''
      WRITE(6,*)'-----------------------------------------------------'
      END IF
c* Initialize the material volume fraction in the cell:
c-----------------------------------------------------      
      IF(ISHAPE.EQ.11) THEN
         CALL INITF3D(FUNC3D1,IPV,NC,NIPV,NTP,NTS,NTV,TOL,VERTP,VF,
     -        XNS,YNS,ZNS)
      ELSE IF(ISHAPE.EQ.12) THEN
         CALL INITF3D(FUNC3D2,IPV,NC,NIPV,NTP,NTS,NTV,TOL,VERTP,VF,
     -        XNS,YNS,ZNS)
      END IF
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE INITF3D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Material volume fraction in the selected cell:',VF
      WRITE(6,*)''
      WRITE(6,*)'-----------------------------------------------------'
      END PROGRAM TEST3D
c-------------------------- END OF TEST3D ----------------------------c
c---------------------------------------------------------------------c
