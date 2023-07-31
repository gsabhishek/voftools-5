c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              TEST2D                                 c
c---------------------------------------------------------------------c
c            Copyright (C) 2020 J. Lopez and J. Hernandez             c
c---------------------------------------------------------------------c
c      Test program to solve the local volume enforcement problem,    c
c      calculate the distance from a given point to the interfacial   c
c      segment and initialize the material volume fraction in a cell  c
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
      PROGRAM TEST2D
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polygon
      DIMENSION IPV(NV),VERTP(NV,2)
c* Work polygon 0
      DIMENSION IPV0(NV),VERTP0(NV,2)
c* Segment
      DIMENSION X(2),Y(2)
c* External function to define the interface shape      
      EXTERNAL FUNC2D1, FUNC2D2
      
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|          TEST PROGRAM IN 2D OF VOFTools           |'
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
      IF(ICELLTYPE.LT.1.OR.ICELLTYPE.GT.106.OR.(ICELLTYPE.GT.6.AND.
     -     ICELLTYPE.LT.101)) THEN
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|*********** WARNING FOR CELL SELECTION ************|'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,*) "1.- Edit the vofvardef file."
         WRITE(6,*) "2.- Choose an appropriate ICELLTYPE value (",
     -               "between 1 and 6 for convex cells or between 101 ", 
     -               "and 106 for non-convex cells)."
         STOP
      END IF
      IF(ISHAPE.LT.1.OR.ISHAPE.GT.2) THEN
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|**** WARNING FOR MATERIAL BODY SHAPE SELECTION ****|'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
         WRITE(6,*) "1.- Edit the vofvardef file."
         WRITE(6,*) "2.- Choose an appropriate ISHAPE value (1 for a ",
     -               "circle or 2 for an ellipse)."
         STOP
      END IF
c* The ZNC and ZP components are ignored for this 2D test program      
      IF(ICELLTYPE.EQ.1) THEN
c* Square mesh
         CALL SQUAREMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.2) THEN
c* Hexagonal mesh
         CALL HEXAGOMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.3) THEN         
c* Triangular mesh
         CALL TRIANGLEMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.4) THEN
c* Quadrangular mesh
         CALL QUADRANGLEMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.5) THEN
c* Pentagonal mesh
         CALL PENTAGONMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.6) THEN
c* Irregular hexagonal mesh
         CALL HEXAGONMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.101) THEN
c* Non-convex quadrangle
         CALL NCQUADRANGLEMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.102) THEN
c* Non-convex pentagon
         CALL NCPENTAGONMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.103) THEN
c* Non-convex hexagon
         CALL NCHEXAGONMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.104) THEN
c* Non-convex stellated hexagon
         CALL NCSHEXAGONMESH(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.105) THEN
c* Hollowed square
         CALL NCHOLLOWEDSQUARE(IPV,NTP,NTV,VERTP)
      ELSEIF(ICELLTYPE.EQ.106) THEN
c* Non-convex multi-square cell
         CALL NCMULTISQUARE(IPV,NTP,NTV,VERTP)
      ENDIF
c* Print the original polygonal cell to the file 'pol00000.out'
      IFILE=0
      CALL POLOUT2D(IFILE,IPV,NTV,VERTP)
c* Calculate the area VT of the cell:
c------------------------------------
      CALL TOOLV2D(IPV,NTV,VERTP,VT)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE TOOLV2D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Area of the selected cell:',VT
      WRITE(6,*)''
c* Calculate the volume of liquid in the cell
      V=F*VT
c* Solve the local volume enforcement problem:
c--------------------------------------------
      CALL ENFORV2D(C,IPV,NTP,NTV,V,VT,VERTP,XNC,YNC)
c.. For rectangular cells, like that defined in the SQUAREMESH subroutine,
c.. it is more efficient to use the analytical method of Scardovelli and 
c.. Zaleski [Journal of Computational Physics, 164 (2000) 228-237], which
c.. was proposed to be used specifically for this type of cells. To use
c.. this method, uncomment the following three lines.
c      DX=DABS(VERTP(1,1)-VERTP(2,1))    !x-side length
c      DY=DABS(VERTP(2,2)-VERTP(3,2))    !y-side length
c      CALL ENFORV2DSZ(C,DX,DY,V,VERTP,XNC,YNC)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|        OUTPUT OF THE ENFORV2D SUBROUTINE          |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Solution of the problem:',C
      WRITE(6,*)''
c* copy the original polygon that defines the cell to the working
c* polygon 0
      CALL CPPOL2D(IPV,IPV0,NTP,NTP0,NTV,NTV0,VERTP,VERTP0)
c* intersection between the cell and the line defined as X·NC+C=0
      CALL INTE2D(C,ICONTN,ICONTP,IPV0,NTP0,NTV0,VERTP0,XNC,YNC)
c* Print the truncated polygonal cell to the file 'pol00001.out'
      IFILE=1
      CALL POLOUT2D(IFILE,IPV0,NTV0,VERTP0)
      IF(ICELLTYPE.GE.1.AND.ICELLTYPE.LE.6) THEN
c* Calculate the distance from a given point P to the interfacial segment:
c------------------------------------------------------------------------
c* interfacial segment defined as the last edge of the truncated polygon 0
         X(1)=VERTP0(NTP0-1,1)
         Y(1)=VERTP0(NTP0-1,2)
         X(2)=VERTP0(NTP0,1)
         Y(2)=VERTP0(NTP0,2)
c* calculate the distance
         CALL DIST2D(D,X,Y,XP,YP)
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'|         OUTPUT OF THE DIST2D SUBROUTINE           |'
      WRITE(6,*)'|---------------------------------------------------|'
      WRITE(6,*)'-----------------------------------------------------'
      WRITE(6,*)''
      WRITE(6,*)'Distance from P to the interfacial segment:',D
      WRITE(6,*)''
      WRITE(6,*)'-----------------------------------------------------'
      END IF
c* Initialize the material volume fraction in the cell:
c-----------------------------------------------------      
      IF(ISHAPE.EQ.1) THEN
         CALL INITF2D(FUNC2D1,IPV,NC,NTP,NTV,TOL,VERTP,VF)
      ELSEIF(ISHAPE.EQ.2) THEN
         CALL INITF2D(FUNC2D2,IPV,NC,NTP,NTV,TOL,VERTP,VF)
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
      END PROGRAM TEST2D
c-------------------------- END OF TEST2D ----------------------------c
c---------------------------------------------------------------------c

