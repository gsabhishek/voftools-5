c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                           V O F T o o l s                           c
c                                                                     c
c        A package of FORTRAN subroutines with analytical and         c
c         geometrical tools for VOF methods in general grids          c
c                       and Cartesian geometry                        c
c                                                                     c
c                      (Version 5, January 2020)                      c
c                                                                     c
c                                                                     c
c           Copyright (C) 2020 J. Lopez* and J. Hernandez**           c
c                                                                     c
c                                                                     c
c  *Dpto. Ingenieria Mecanica, Materiales y Fabricacion, UPCT. 30202, c
c   Cartagena, Spain.                                                 c
c **Dpto. Mecanica, UNED. 28040, Madrid, Spain.                       c
c                                                                     c
c   For more information, please contact: joaquin.lopez@upct.es       c
c---------------------------------------------------------------------c
c This library is free software: you can redistribute it and/or       c
c modify it under the terms of the GNU Lesser General Public License  c
c as published by the Free Software Foundation, either version 3 of   c
c the License, or (at your option) any later version.                 c
c                                                                     c
c This library is distributed in the hope that it will be useful,     c
c but WITHOUT ANY WARRANTY; without even the implied warranty of      c
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       c
c GNU Lesser General Public License for more details.                 c
c                                                                     c
c You should have received a copy of the GNU Lesser General Public    c
c License along with this library.  If not, see                       c
c <http://www.gnu.org/licenses/>.                                     c
c---------------------------------------------------------------------c
c List of subroutines:                                                c
c=====================                                                c
c                                                                     c
c 3D subroutines:                                                     c
c---------------                                                      c
c ENFORV3D ----> solve the local volume enforcement problem in 3D     c
c                using the CIBRAVE method                             c
c ENFORV3DSZ---> solve the local volume enforcement problem for       c
c                rectangular parallelepiped cells using the method    c 
c                of Scardovelli and Zaleski [Journal of Computational c
c                Physics, 164 (2000) 228-237]                         c
c NEWPOL3D ----> vertex indices arrangement of the truncated          c
c                polyhedron                                           c
c INTE3D   ----> obtain the polyhedron truncated by a plane           c
c TOOLV3D  ----> compute the volume of a polyhedron                   c
c CPPOL3D  ----> copy a polyhedron into a new one                     c
c RESTORE3D----> restore the structure of a polyhedron                c
c EQSOL3D  ----> solve analytically a cubic equation                  c
c NEWTON3D ----> solve a cubic eq. using Newton-Raphson iterations    c
c DIST3D   ----> compute the distance from a point to a polygon       c
c INITF3D  ----> initialize the material volume fraction in a cell    c
c POLOUT3D ----> write in an external VTK-format file the geometry of c
c                a polyhedron                                         c
c                                                                     c
c 2D subroutines:                                                     c
c---------------                                                      c
c ENFORV2D ----> solve the local volume enforcement problem in 2D     c
c                using the CIBRAVE method                             c
c ENFORV2DSZ---> solve the local volume enforcement problem for       c
c                rectangular cells using the method of Scardovelli    c 
c                and Zaleski [Journal of Computational Physics, 164   c
c                (2000) 228-237]                                      c
c NEWPOL2D ----> vertex indices arrangement of the truncated          c
c                polygon                                              c
c INTE2D   ----> obtain the polygon truncated by a line               c
c TOOLV2D  ----> compute the area of a polygon                        c
c CPPOL2D  ----> copy a polygon into a new one                        c
c RESTORE2D----> restore the structure of a polygon                   c
c DIST2D   ----> compute the distance from a point to a segment       c
c INITF2D  ----> initialize the material area fraction in a cell      c
c POLOUT2D ----> write in an external file the vertex coordinatex x   c
c                and y of a polygon in two columns format             c
c                                                                     c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            ENFORV3D                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c IPV      = array containing the global indices of the polyhedron    c
c            vertices                                                 c
c NIPV     = number of vertices of each face                          c
c NTP      = last global vertex index (note that if the polyhedron    c
c            is not previously truncated, then NTP=NTV)               c
c NTS      = total number of faces                                    c
c NTV      = total number of vertices                                 c
c V        = liquid volume                                            c
c VT       = total volume of the polyhedron                           c
c VERTP    = vertex coordinates of the polyhedron                     c
c XNS, ... = unit-lenght normals to the faces of the polyhedron       c
c XNC, ... = unit-lenght normal to the new face boundaries on         c
c            \Gamma_c                                                 c
c On return:                                                          c
c===========                                                          c
c C        = solution of the problem                                  c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE ENFORV3D(C,IPV,NIPV,NTP,NTS,NTV,V,VT,VERTP,XNC,
     -     XNS,YNC,YNS,ZNC,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polyhedron
      DIMENSION CS(NS),IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),YNS(NS),
     -     ZNS(NS)
c* Working polyhedron 0
      DIMENSION CS0(NS),IPV0(NS,NV),NIPV0(NS),VERTP0(NV,3),XNS0(NS),
     -     YNS0(NS),ZNS0(NS)
c* Working variables
      DIMENSION BETXE(NS,NV),BETYE(NS,NV),BETZE(NS,NV),IA(NV),IPIA0(NV),
     -     IPIA1(NV),ISCUT(NS),LISTV(NV),PHIV(NV),X0(NS,NV),Y0(NS,NV),
     -     Z0(NS,NV)
      DIMENSION SUMK(NS),SUML(NS),SUMM(NS)
      
      IF(VT.LE.0.0) THEN
         WRITE(6,*) 'THE POLYHEDRON HAS NULL OR NEGATIVE AREA.'
         RETURN
      END IF
C.. To execute in quad precision may be convinient to use a lower value 
c.. for TOLC such as, for example, 1.0D-20
      TOLC=1.0D-12
      DO IS=1,NTS
         IP=IPV(IS,1)
         CS(IS)=-XNS(IS)*VERTP(IP,1)-YNS(IS)*VERTP(IP,2)-ZNS(IS)*
     -        VERTP(IP,3)
      END DO
c* If the polyhedron was previously truncated, the polyhedron must   *c 
c* be restored before enforcing volume conservation. Note that       *c
c* the index NTP of the last vertex inserted after the truncation    *c
c* procudure may be higher than the total number of vertices NTV of  *c
c* the truncated polyhedron. It should be mentioned that in the      *c 
c* context of a PLIC-VOF method the volume enforcement operations are*c
c* generally made before any intersection operations, and therefore, *c 
c* the restoration procedure is not generally required. At any case, *c
c* the corresponding RESTORE3D routine is also supplied.             *c
      IF(NTP.GT.NTV) CALL RESTORE3D(CS,IPV,NIPV,NTP,NTS,NTV,VERTP,
     -     XNS,YNS,ZNS)
      VAUX=V
      LISTV(1)=1
      DO IV=1,NTV
         PHIV(IV)=XNC*VERTP(IV,1)+YNC*VERTP(IV,2)+ZNC*VERTP(IV,3)
c* Ordered list of global vertex indices
         DO I=1,IV-1
            IF(PHIV(IV).GT.PHIV(LISTV(I))) THEN
               DO II=IV,I+1,-1
                  LISTV(II)=LISTV(II-1)
               END DO
               LISTV(I)=IV
               GOTO 10
            END IF
         END DO
         LISTV(IV)=IV
 10      CONTINUE
      END DO

      INVERT=0
      XNCOR=XNC
      YNCOR=YNC
      ZNCOR=ZNC
      IMIN=1
      IMAX=NTP
      VMIN=DBLE(0)
      VMAX=VT

C* Obtain the tentative solution bracketing by interpolation
      IMAXLOLD=NTP+1
 22   CONTINUE
      PHIINT=PHIV(LISTV(IMIN))-(PHIV(LISTV(IMIN))-PHIV(LISTV(IMAX)))*
     -     (V-VMIN)/(VMAX-VMIN)
      IMAXL=0
      IMINL=0
      DO IP=IMIN+1,IMAX
         I=IP
         IF(PHIV(LISTV(IP)).LT.PHIINT) THEN
            IMAXL=IP
            IMINL=IP-1
            GOTO 11
         END IF
      END DO
      IF(IMAXL.EQ.0.AND.IMINL.EQ.0) THEN
         C=-PHIINT
         RETURN
      END IF      
 11   CONTINUE

      CMAX=PHIV(LISTV(IMINL))
      CMIN=PHIV(LISTV(IMAXL))
            
      IF((NTP-IMAXL).LT.(IMINL-1)) THEN
         INVERT=1
         CAUX=CMIN
         CMIN=-CMAX
         CMAX=-CAUX
         VAUX=VT-V
         XNC=-XNCOR
         YNC=-YNCOR
         ZNC=-ZNCOR
      ELSE
         INVERT=0
         VAUX=V
         XNC=XNCOR
         YNC=YNCOR
         ZNC=ZNCOR
      END IF
      DO I=1,NTP
         IF(I.LE.IMINL) THEN
            IA(LISTV(I))=1-INVERT
         ELSE
            IA(LISTV(I))=INVERT
         END IF
      END DO
c* End of procedure SETIA
      CALL CPPOL3D(CS0,CS,IPV0,IPV,NIPV0,NIPV,NTP0,NTP,NTS0,NTS,NTV0,
     -     NTV,VERTP0,VERTP,XNS0,XNS,YNS0,YNS,ZNS0,ZNS)
c* Construction of the new polyhedron
      NTS00=NTS0
      CALL NEWPOL3D(IA,IPIA0,IPIA1,IPV0,ISCUT,NIPV0,NTP0,NTS0,
     -     NTV0,XNC,XNS0,YNC,YNS0,ZNC,ZNS0)
      if(nts0.le.nts00) then ! disjoint regions may produce this situation
         if((imax-imaxl).gt.(iminl-imin)) then
            imaxl=imaxl+1
            iminl=imaxl-1
         else
            iminl=iminl-1
            imaxl=iminl+1
         end if
         goto 11
      end if      
c* Contributions of the new faces \Gamma_c
      DO IS=NTS00+1,NTS0
      DO IV=1,NIPV0(IS)
         IPI=IPV0(IS,IV)
         IPF=IPIA1(IPI)
         IP=IPIA0(IPI)
         XV=VERTP0(IPF,1)-VERTP0(IP,1)
         YV=VERTP0(IPF,2)-VERTP0(IP,2)
         ZV=VERTP0(IPF,3)-VERTP0(IP,3)
         CIS=XNS0(IS)*XV+YNS0(IS)*YV+ZNS0(IS)*ZV
         if(cis.eq.0.0) then !cut edge and \Gamma_c normal are perpendicular
            BETXE(IS,IV)=DBLE(0)
            BETYE(IS,IV)=DBLE(0)
            BETZE(IS,IV)=DBLE(0)
         else
            BETXE(IS,IV)=-XV/CIS
            BETYE(IS,IV)=-YV/CIS
            BETZE(IS,IV)=-ZV/CIS
         end if
         COEF=XNS0(IS)*VERTP0(IPF,1)+YNS0(IS)*VERTP0(
     -        IPF,2)+ZNS0(IS)*VERTP0(IPF,3)
         X0(IS,IV)=VERTP0(IPF,1)+BETXE(IS,IV)*COEF
         Y0(IS,IV)=VERTP0(IPF,2)+BETYE(IS,IV)*COEF
         Z0(IS,IV)=VERTP0(IPF,3)+BETZE(IS,IV)*COEF         
      END DO
      END DO
C* Contributions of the rest of faces
      DO IS=1,NTS00
         DO IV=1,NIPV0(IS)
            IP=IPV0(IS,IV)
            IF(IA(IP).EQ.1) THEN
               X0(IS,IV)=VERTP0(IP,1)
               Y0(IS,IV)=VERTP0(IP,2)
               Z0(IS,IV)=VERTP0(IP,3)  
               BETXE(IS,IV)=DBLE(0)
               BETYE(IS,IV)=DBLE(0)
               BETZE(IS,IV)=DBLE(0)
            ELSE
               DO IS1=NTS00+1,NTS0
               ISC=IS1
               DO IV1=1,NIPV0(IS1)
                  IP1=IPV0(IS1,IV1)
                  IF(IP1.EQ.IP) GOTO 100
               END DO
               END DO
 100           CONTINUE
               BETXE(IS,IV)=BETXE(ISC,IV1)
               BETYE(IS,IV)=BETYE(ISC,IV1)
               BETZE(IS,IV)=BETZE(ISC,IV1)
               X0(IS,IV)=X0(ISC,IV1)
               Y0(IS,IV)=Y0(ISC,IV1)
               Z0(IS,IV)=Z0(ISC,IV1)               
            END IF
         END DO
      END DO
C** Vectors K,L and M
      DO IS=NTS00+1,NTS0
         ISCUT(IS)=1
      END DO
      DO IS=1,NTS0
         SUMK(IS)=DBLE(0)
         SUML(IS)=DBLE(0)
         SUMM(IS)=DBLE(0)
         IF(DABS(YNS0(IS)).GE.DABS(XNS0(IS)).AND.DABS(YNS0(IS)).GE.
     -        DABS(ZNS0(IS))) THEN
            IPROJ=2
            DNMAX=YNS0(IS)               
         ELSEIF(DABS(ZNS0(IS)).GE.DABS(XNS0(IS)).AND.DABS(ZNS0(IS)).GE.
     -           DABS(YNS0(IS))) THEN
            IPROJ=3
            DNMAX=ZNS0(IS)
         ELSE
            IPROJ=1
            DNMAX=XNS0(IS)
         END IF
         IH=INT((NIPV0(IS)-2)/2)
         DO I=2,IH+1
            IP=2*I
            IP1=IP-1
            IP2=IP-2
            IF(IPROJ.EQ.1) THEN
               YV1=Y0(IS,IP1)-Y0(IS,1)
               ZV1=Z0(IS,IP1)-Z0(IS,1)
               YV2=Y0(IS,IP)-Y0(IS,IP2)
               ZV2=Z0(IS,IP)-Z0(IS,IP2)
               YE1=BETYE(IS,IP1)-BETYE(IS,1)
               ZE1=BETZE(IS,IP1)-BETZE(IS,1)
               YE2=BETYE(IS,IP)-BETYE(IS,IP2)
               ZE2=BETZE(IS,IP)-BETZE(IS,IP2)
               SUMK(IS)=SUMK(IS)+YV1*ZV2-ZV1*YV2
               SUML(IS)=SUML(IS)+YV1*ZE2-ZV1*YE2-(YV2*ZE1-ZV2*YE1)
               SUMM(IS)=SUMM(IS)+YE1*ZE2-ZE1*YE2
            ELSEIF(IPROJ.EQ.2) THEN
               XV1=X0(IS,IP1)-X0(IS,1)
               ZV1=Z0(IS,IP1)-Z0(IS,1)
               XV2=X0(IS,IP)-X0(IS,IP2)
               ZV2=Z0(IS,IP)-Z0(IS,IP2)
               XE1=BETXE(IS,IP1)-BETXE(IS,1)
               ZE1=BETZE(IS,IP1)-BETZE(IS,1)
               XE2=BETXE(IS,IP)-BETXE(IS,IP2)
               ZE2=BETZE(IS,IP)-BETZE(IS,IP2)
               SUMK(IS)=SUMK(IS)+ZV1*XV2-XV1*ZV2
               SUML(IS)=SUML(IS)+ZV1*XE2-XV1*ZE2-(ZV2*XE1-XV2*ZE1)
               SUMM(IS)=SUMM(IS)+ZE1*XE2-XE1*ZE2
            ELSE
               XV1=X0(IS,IP1)-X0(IS,1)
               YV1=Y0(IS,IP1)-Y0(IS,1)
               XV2=X0(IS,IP)-X0(IS,IP2)
               YV2=Y0(IS,IP)-Y0(IS,IP2)
               XE1=BETXE(IS,IP1)-BETXE(IS,1)
               YE1=BETYE(IS,IP1)-BETYE(IS,1)
               XE2=BETXE(IS,IP)-BETXE(IS,IP2)
               YE2=BETYE(IS,IP)-BETYE(IS,IP2)
               SUMK(IS)=SUMK(IS)+XV1*YV2-YV1*XV2
               SUML(IS)=SUML(IS)+XV1*YE2-YV1*XE2-(XV2*YE1-YV2*XE1)
               SUMM(IS)=SUMM(IS)+XE1*YE2-YE1*XE2
            END IF
         END DO
         IF(2*(IH+1).LT.NIPV0(IS)) THEN
            IF(IPROJ.EQ.1) THEN
               YV1=Y0(IS,NIPV0(IS))-Y0(IS,1)
               ZV1=Z0(IS,NIPV0(IS))-Z0(IS,1)
               YV2=Y0(IS,1)-Y0(IS,NIPV0(IS)-1)
               ZV2=Z0(IS,1)-Z0(IS,NIPV0(IS)-1)
               YE1=BETYE(IS,NIPV0(IS))-BETYE(IS,1)
               ZE1=BETZE(IS,NIPV0(IS))-BETZE(IS,1)
               YE2=BETYE(IS,1)-BETYE(IS,NIPV0(IS)-1)
               ZE2=BETZE(IS,1)-BETZE(IS,NIPV0(IS)-1)
               SUMK(IS)=SUMK(IS)+YV1*ZV2-ZV1*YV2
               SUML(IS)=SUML(IS)+YV1*ZE2-ZV1*YE2-(YV2*ZE1-ZV2*YE1)
               SUMM(IS)=SUMM(IS)+YE1*ZE2-ZE1*YE2
            ELSEIF(IPROJ.EQ.2) THEN
               XV1=X0(IS,NIPV0(IS))-X0(IS,1)
               ZV1=Z0(IS,NIPV0(IS))-Z0(IS,1)
               XV2=X0(IS,1)-X0(IS,NIPV0(IS)-1)
               ZV2=Z0(IS,1)-Z0(IS,NIPV0(IS)-1)
               XE1=BETXE(IS,NIPV0(IS))-BETXE(IS,1)
               ZE1=BETZE(IS,NIPV0(IS))-BETZE(IS,1)
               XE2=BETXE(IS,1)-BETXE(IS,NIPV0(IS)-1)
               ZE2=BETZE(IS,1)-BETZE(IS,NIPV0(IS)-1)
               SUMK(IS)=SUMK(IS)+ZV1*XV2-XV1*ZV2
               SUML(IS)=SUML(IS)+ZV1*XE2-XV1*ZE2-(ZV2*XE1-XV2*ZE1)
               SUMM(IS)=SUMM(IS)+ZE1*XE2-XE1*ZE2
            ELSE
               XV1=X0(IS,NIPV0(IS))-X0(IS,1)
               YV1=Y0(IS,NIPV0(IS))-Y0(IS,1)
               XV2=X0(IS,1)-X0(IS,NIPV0(IS)-1)
               YV2=Y0(IS,1)-Y0(IS,NIPV0(IS)-1)
               XE1=BETXE(IS,NIPV0(IS))-BETXE(IS,1)
               YE1=BETYE(IS,NIPV0(IS))-BETYE(IS,1)
               XE2=BETXE(IS,1)-BETXE(IS,NIPV0(IS)-1)
               YE2=BETYE(IS,1)-BETYE(IS,NIPV0(IS)-1)
               SUMK(IS)=SUMK(IS)+XV1*YV2-YV1*XV2
               SUML(IS)=SUML(IS)+XV1*YE2-YV1*XE2-(XV2*YE1-YV2*XE1)
               SUMM(IS)=SUMM(IS)+XE1*YE2-YE1*XE2
            END IF
         END IF
         SUMK(IS)=SUMK(IS)/DNMAX
         SUML(IS)=SUML(IS)/DNMAX
         SUMM(IS)=SUMM(IS)/DNMAX
      END DO
C** Coefficients of the analytical equation: C3·x^3+C2·x^2+C1·x+C0=0
      C3=DBLE(0)
      C2=DBLE(0)
      C1=DBLE(0)
      DO IS=NTS00+1,NTS0
         C3=C3+SUMM(IS)
         C2=C2+SUML(IS)
         C1=C1+SUMK(IS)
      END DO
      C0=DBLE(6)*VAUX
      DO IS=1,NTS00
         C2=C2+SUMM(IS)*CS0(IS)
         C1=C1+SUML(IS)*CS0(IS)
         C0=C0+SUMK(IS)*CS0(IS)
      END DO
      VMAXL=VAUX-(C3*CMIN*CMIN*CMIN+C2*CMIN*CMIN+C1*
     -     CMIN+C0)/DBLE(6)
      VMINL=VAUX-(C3*CMAX*CMAX*CMAX+C2*CMAX*CMAX+C1*
     -     CMAX+C0)/DBLE(6)
      
      IF(INVERT.EQ.1) THEN
         VMINLL=VT-VMAXL
         VMAXL=VT-VMINL
         VMINL=VMINLL
      END IF
      SV=(VMINL-V)*(VMAXL-V)
      IF(SV.GT.DBLE(0).AND.(IMAX-IMIN).GT.1.AND.IMAXLOLD.NE.IMAXL) THEN
         IF(VMAXL.GT.V) THEN
            VMAX=VMINL
            IMAX=IMINL
         ELSE
            VMIN=VMAXL
            IMIN=IMAXL
         END IF
         IMAXLOLD=IMAXL
         GOTO 22
      END IF
      CALL EQSOL3D(C0,C1,C2,C3,CMIN,CMAX,C)
      ERRV=C3*C*C*C+C2*C*C+C1*C+C0
      IF(DABS(ERRV).GT.TOLC)CALL NEWTON3D(C0,C1,C2,C3,CMIN,CMAX,C,
     -     ISOL)
      IF(INVERT.EQ.0) THEN 
         C=-C
         CMAX2=-CMIN
         CMIN2=-CMAX
      ELSE
         XNC=-XNC
         YNC=-YNC
         ZNC=-ZNC         
         CMAX2=CMAX
         CMIN2=CMIN
      END IF
c..   Forced ajustment 
      IF(C.LT.CMIN2.OR.C.GT.CMAX2) THEN
         CALL CPPOL3D(CS0,CS,IPV0,IPV,NIPV0,NIPV,NTP0,NTP,NTS0,NTS,NTV0,
     -        NTV,VERTP0,VERTP,XNS0,XNS,YNS0,YNS,ZNS0,ZNS)
         CALL INTE3D(Cmin2,ICONTN,ICONTP,IPV0,NIPV0,NTP0,NTS0,NTV0,
     -        VERTP0,XNC,XNS0,YNC,YNS0,ZNC,ZNS0)
         IF(ICONTP.EQ.0) THEN
            VOLMIN=DBLE(0)
         ELSE
            CALL TOOLV3D(IPV0,NIPV0,NTS0,VERTP0,VOLmin,XNS0,YNS0,ZNS0)
         END IF
         CALL CPPOL3D(CS0,CS,IPV0,IPV,NIPV0,NIPV,NTP0,NTP,NTS0,NTS,NTV0,
     -        NTV,VERTP0,VERTP,XNS0,XNS,YNS0,YNS,ZNS0,ZNS)
         CALL INTE3D(CMAX2,ICONTN,ICONTP,IPV0,NIPV0,NTP0,NTS0,NTV0,
     -        VERTP0,XNC,XNS0,YNC,YNS0,ZNC,ZNS0)
         IF(ICONTP.EQ.0) THEN
            VOLMAX=0.0
         ELSE
            CALL TOOLV3D(IPV0,NIPV0,NTS0,VERTP0,VOLmax,XNS0,YNS0,ZNS0)
         END IF
         IF(DABS(VOLMIN-V).LT.DABS(VOLMAX-V)) THEN
            C=CMIN2
         ELSE
            C=CMAX2
         END IF
      END IF
      RETURN
      END
c---------------------   END OF ENFORV3D   ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                          ENFORV3DSZ                                 c
c... Scardovelli and Zaleski version for rectangular parallelepiped   c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c DX, ...  = side lengths of the rectangular parallelepiped           c
c VERTP    = vertex coordinates of the polyhedron                     c
c XNC, ... = unit-lenght normal to the new face \Gamma_c              c
c V        = liquid volume                                            c
c On return:                                                          c
c===========                                                          c
c C        = solution of the problem                                  c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE ENFORV3DSZ(C,DX,DY,DZ,V,VERTP,XNC,YNC,ZNC)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c.. Original polyhedron
      DIMENSION VERTP(NV,3)
      DOUBLE PRECISION M(3),M1,M12,M2,M3,MAUX
      TOLE=1.0D-09
      CMIN=1.0D+14
      CMAX=-1.0D+14
      VT=DX*DY*DZ
      VBACK=V
      V=V/VT
c.. The vertex indices of the rectangular parallelepiped are supposed 
c.. to be listed from 1 to 8.
      DO 10 I=1,8
         CI=-(VERTP(I,1)*XNC+VERTP(I,2)*YNC+VERTP(I,3)*ZNC)
         IF(CI.LE.CMIN) THEN
            CMIN=CI
            IMIN=I
         END IF
         IF(CI.GE.CMAX) THEN
            CMAX=CI
            IMAX=I
         END IF
 10   CONTINUE
c.. If the liquid volume fraction is higher than 0.5, solve the inverse problem
      IF((VBACK/VT).LE.(DBLE(1)/DBLE(2))) THEN
         CI=CMIN
         I=IMIN
      ELSE
         CI=CMAX
         I=IMAX
         V=DBLE(1)-V
      END IF
c.. Normalize the plane equation
      SN=DABS(XNC)+DABS(YNC)+DABS(ZNC)
      XM=XNC/SN
      YM=YNC/SN
      ZM=ZNC/SN      
      XMI=XM*DX
      YMI=YM*DY
      ZMI=ZM*DZ
      SN=DABS(XMI)+DABS(YMI)+DABS(ZMI)
      XM=DABS(XMI)/SN
      YM=DABS(YMI)/SN
      ZM=DABS(ZMI)/SN
c.. Region limits
      M(1)=XM
      M(2)=YM
      M(3)=ZM
      DO I=1,2
         I2=I+1
         IF(M(I2).LT.M(I)) THEN
            MAUX=M(I)
            M(I)=M(I2)
            M(I2)=MAUX
         END IF
      END DO
      IF(M(2).LT.M(1)) THEN
         MAUX=M(1)
         M(1)=M(2)
         M(2)=MAUX         
      END IF
      M1=M(1)
      M2=M(2)
      M3=M(3)
      M12=M1+M2
      V1=(M1*M1)/DMAX1(DBLE(6)*M2*M3,TOLE)
      V2=V1+(M2-M1)/(DBLE(2)*M3)
      IF(M3.LT.M12) THEN
         V3=((DBLE(3)*M12-M3)*(M3*M3)+(M1-DBLE(3)*M3)*(M1*M1)+
     -        (M2-DBLE(3)*M3)*(M2*M2))/(DBLE(6)*M1*M2*M3)
      ELSE
         V3=M12/(DBLE(2)*M3)
      END IF
      IF(V.GE.V2.AND.V.LT.V3) THEN
         A3=DBLE(-1)
         A2=DBLE(3)*M12/A3
         A1=DBLE(-3)*(M1*M1+M2*M2)/A3
         A0=(M1*M1*M1+M2*M2*M2-(DBLE(6)*M1*M2*M3*V))/A3
         A3=DBLE(1)
      ELSE
         IF(V.GE.V3.AND.V.LE.(DBLE(1)/DBLE(2)).AND.M3.LT.M12) THEN
            A3=DBLE(-2)
            A2=DBLE(3)/A3
            A1=DBLE(-3)*(M1*M1+M2*M2+M3*M3)/A3
            A0=(M1*M1*M1+M2*M2*M2+M3*M3*M3-(DBLE(6)*M1*M2*M3*V))/A3
            A3=DBLE(1)
         END IF
      END IF
c.. Solution of the inverse problem
      IF(V.GE.DBLE(0).AND.V.LT.V1) THEN
         ALPHA=(DBLE(6)*M1*M2*M3*V)**(DBLE(1)/DBLE(3))
         GOTO 20
      END IF
      IF(V.GE.V1.AND.V.LT.V2) THEN
         ALPHA=(DBLE(1)/DBLE(2))*(M1+(M1*M1+DBLE(8)*M2*M3*(V-V1))**(
     -        DBLE(1)/DBLE(2)))
         GOTO 20
      END IF
      IF((V.GE.V2.AND.V.LT.V3).OR.(V.GE.V3.AND.V.LE.(DBLE(1)/DBLE(2)).
     -     AND.M3.LT.M12)) THEN
         P0=(A1/DBLE(3))-((A2*A2)/DBLE(9))
         Q0=((A1*A2-DBLE(3)*A0)/DBLE(6))-(A2*A2*A2)/DBLE(27)
         THETA=DACOS(Q0/((DBLE(-1)*P0*P0*P0)**(DBLE(1)/DBLE(2))))/
     -        DBLE(3)
         ALPHA=((DBLE(-1)*P0)**(DBLE(1)/DBLE(2)))*(DSIN(THETA)*(
     -        DBLE(3)**(DBLE(1)/DBLE(2)))-DCOS(THETA))-(A2/DBLE(3))
      ELSE
         ALPHA=M3*V+M12/DBLE(2)
      END IF
 20   CONTINUE
      IF((VBACK/VT).LE.(DBLE(1)/DBLE(2))) THEN
         C=CMIN+ALPHA*DABS(CMAX-CMIN)
      ELSE
         C=CMAX-ALPHA*DABS(CMAX-CMIN)
      END IF
      V=VBACK
      RETURN
      END
c---------------------   END OF ENFORV3DSZ   -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            NEWPOL3D                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c XNS0, ...= unit-lenght normals to the faces of the original polyh.  c
c XNC, ... = unit-lenght normal to the new face \Gamma_c              c
c IA       = 0 if \Gamma_c points out of the vertex, 1 otherwise      c
c IPV0     = array containing the global indices of the original pol. c
c            vertices                                                 c
c NIPV0    = number of vertices of each face                          c
c NTP0     = last global vertex index                                 c
c NTS0     = total number of faces                                    c
c NTV0     = total number of vertices                                 c
c On return:                                                          c
c===========                                                          c
c XNS0, ...= unit-lenght normals to the faces of the truncated pol.   c
c XNC, ... = unit-lenght normal to the new face \Gamma_c              c
c IPV0     = array containing the global indices of the truncat. pol. c
c            vertices                                                 c
c NIPV0    = number of vertices of each face                          c
c NTP0     = last global vertex index                                 c
c NTS0     = total number of faces                                    c
c NTV0     = total number of vertices                                 c
c IPIA0    = global vertex index of the original polihedron with IA=0 c
c            and which is in the edge containing the intersection     c
c            point                                                    c
c IPIA1    = global vertex index of the original polihedron with IA=1 c
c            and which is in the edge containing the intersection     c
c            point                                                    c
c ISCUT    = 1 if the faces is truncated. 0 if the faces is not       c
c            truncated                                                c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE NEWPOL3D(IA,IPIA0,IPIA1,IPV0,ISCUT,NIPV0,NTP0,NTS0,
     -     NTV0,XNC,XNS0,YNC,YNS0,ZNC,ZNS0)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original/truncated polyhedron
      DIMENSION IPV0(NS,NV),NIPV0(NS),XNS0(NS),YNS0(NS),ZNS0(NS)
c* Work polyhedron 1
      DIMENSION IPV1(NS,NV),NIPV1(NS)
c* Work variables 1
      DIMENSION IA(NV),IPIA0(NV),IPIA1(NV),ISCUT(NS)
c* Work variables 2
      DIMENSION NEDGE(NS)
c*  OTRAS VARIABLES:
      DIMENSION ISE(NS,NV),IVISE(NS,NV),IPISE(NV,2),IPMARK(NV)
c* Determination of the cut faces
      ISINI=0
      NISCUT=0
c* NEDGE(IS) = Number of cut edges of the face IS
      DO IS=1,NTS0
         NEDGE(IS)=0
         IF(NIPV0(IS).GT.0) THEN
            ISCUT(IS)=0
            DO IV=1,NIPV0(IS)
               IP=IPV0(IS,IV)
               IV1=IV+1
               IF(IV.EQ.NIPV0(IS)) IV1=1
               IP1=IPV0(IS,IV1)
               IF(IA(IP).NE.IA(IP1)) THEN
                  ISCUT(IS)=1
                  NISCUT=NISCUT+1
                  IF(ISINI.EQ.0) ISINI=IS
                  NEDGE(IS)=NEDGE(IS)+1
               END IF
            END DO
            IF(ISCUT(IS).EQ.0.AND.IA(IPV0(IS,1)).EQ.0) NIPV0(IS)=
     -           -NIPV0(IS)
         END IF
      END DO
c* Disjoint regions may produce NISCUT=0 and both ICONTP and ICONTN \NEQ 0
c.. Cuando el polihedro es no convexo, puede ocurrir que despues de varias
c.. intersecciones sucesivas, aparezcan varios polihedros individuales. Si
c.. un nuevo plano de corte pasa entre ellos, no aparecera ninguna cara 
c.. cortada (NISCUT=0) aunque ICONTP y ICONTN sean ambos distintos de 0.
      IF(NISCUT.EQ.0) THEN
         NIPNEW=0
         NIVNEW=0
         ISNEW=0
         GOTO 50
      END IF
c* Construction of the cut faces
      NIPNEW=NTP0
      DO IS=1,NTS0
         IF(ISCUT(IS).EQ.1) THEN
            NIV=0
            NINT=0
            DO IV=1,NIPV0(IS)
               IP=IPV0(IS,IV)
               IV1=IV+1
               IF(IV1.GT.NIPV0(IS))IV1=1
               IP1=IPV0(IS,IV1)
               IF(IA(IP).EQ.1) THEN
                  NIV=NIV+1
                  IPV1(IS,NIV)=IPV0(IS,IV)
               END IF
               IF(IA(IP).NE.IA(IP1)) THEN
                  NINT=NINT+1
                  NIV=NIV+1
                  IF(IA(IP).EQ.1) THEN
                     IP1I=IP
                     IP0I=IP1
                     ITYPE=2
                  ELSE
                     IP1I=IP1
                     IP0I=IP
                     ITYPE=1
                  END IF
                  DO IS1=1,IS-1
                     DO IE=1,NEDGE(IS1)
                        IPNEW=ISE(IS1,IE)
                        IP0N=IPIA0(IPNEW)
                        IP1N=IPIA1(IPNEW)
                        IF(IP0N.EQ.IP0I.AND.IP1N.EQ.IP1I) THEN
                           ISE(IS,NINT)=IPNEW
                           IPV1(IS,NIV)=IPNEW
                           IVISE(IS,IPNEW)=NIV
                           IPISE(IPNEW,ITYPE)=IS
                           GOTO 10
                        END IF
                     END DO     
                  END DO        
                  NIPNEW=NIPNEW+1
                  IPIA0(NIPNEW)=IP0I
                  IPIA1(NIPNEW)=IP1I
                  IPV1(IS,NIV)=NIPNEW
                  ISE(IS,NINT)=NIPNEW
                  IVISE(IS,NIPNEW)=NIV
                  IPISE(NIPNEW,ITYPE)=IS
               END IF           
 10            CONTINUE
            END DO              
            NIPV1(IS)=NIV
         END IF                 
      END DO                    
c* Construction of the new faces
      NIVNEW=NIPNEW-NTP0
      ISNEW=NTS0
      DO IP=NTP0+1,NIPNEW
         IPMARK(IP)=0
      END DO
      IVNEWT=0
      IPNEW=NTP0+1
c* First point
 40   CONTINUE
      IVNEW=1
      IVNEWT=IVNEWT+1
      ISNEW=ISNEW+1
      IPINI=IPNEW
      IPV0(ISNEW,IVNEW)=IPNEW
      IPMARK(IPNEW)=1
 20   CONTINUE
      IS=IPISE(IPNEW,1)
      IV=IVISE(IS,IPNEW)
      IV1=IV-1
      IF(IV1.EQ.0) IV1=NIPV1(IS)
      IPNEW=IPV1(IS,IV1)
      IF(IPNEW.NE.IPINI) THEN
         IVNEW=IVNEW+1
         IVNEWT=IVNEWT+1
         IPV0(ISNEW,IVNEW)=IPNEW
         IPMARK(IPNEW)=1
         IF(IVNEWT.EQ.NIVNEW) GOTO 30
         GOTO 20
      END IF
      NIPV0(ISNEW)=IVNEW
      DO IPNEW=NTP0+2,NIPNEW
         IF(IPMARK(IPNEW).EQ.0) GOTO 40
      END DO
 30   CONTINUE
      NIPV0(ISNEW)=IVNEW
c* Assign the vertices of the new truncated polyhedron
 50   CONTINUE
      NIV=NIVNEW
      NTPMAX=NIPNEW
      NTSMAX=ISNEW
      DO IS=1,NTS0
         IF(NIPV0(IS).GT.0) THEN
            IF(ISCUT(IS).EQ.1) THEN
               NIPV0(IS)=NIPV1(IS)
               DO IV=1,NIPV1(IS)
                  IPV0(IS,IV)=IPV1(IS,IV)
                  IF(IA(IPV1(IS,IV)).EQ.1) THEN
                     NIV=NIV+1
                     IA(IPV1(IS,IV))=-1
                  END IF
               END DO
            ELSE
               IF(ISCUT(IS).EQ.0.AND.NIPV0(IS).LT.0) NIPV0(IS)=0
               DO IV=1,NIPV0(IS)
                  NTSMAX=MAX0(NTSMAX,IS)
                  IF(IA(IPV0(IS,IV)).EQ.1) THEN
                     NTPMAX=MAX0(NTPMAX,IPV0(IS,IV))
                     NIV=NIV+1
                     IA(IPV0(IS,IV))=-1
                  END IF
               END DO
            END IF
         END IF
      END DO
      DO IP=1,NTP0
         IF(IA(IP).EQ.-1) IA(IP)=1
      END DO
      DO IS=NTS0+1,ISNEW
         XNS0(IS)=-1.0*XNC
         YNS0(IS)=-1.0*YNC
         ZNS0(IS)=-1.0*ZNC
      END DO
      NTV0=NIV
      NTP0=NTPMAX
      NTS0=NTSMAX
      RETURN
      END
c-------------------------- END OF NEWPOL3D --------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                               INTE3D                                c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c XNC, ... = unit-lenght normal to the new face \Gamma_c              c
c C        = constant of the plane containing the new face \Gamma_c   c
c XNS0, ...= unit-lenght normals to the faces of the original pol.    c
c IPV0     = array containing the global indices of the original pol. c
c            vertices                                                 c
c NIPV0    = number of vertices of each face                          c
c NTS0     = total number of faces                                    c
c NTP0     = last global vertex index                                 c
c NTV0     = total number of vertices                                 c
c VERTP0   = vertex coordinates of the original polyhedron            c
c On return:                                                          c
c===========                                                          c
c XNS0, ...= unit-lenght normals to the faces of the truncated pol.   c
c IPV0     = array containing the global indices of the truncat. pol. c
c            vertices                                                 c
c NIPV0    = number of vertices of each face                          c
c NTS0     = total number of faces                                    c
c NTP0     = last global vertex index                                 c
c NTV0     = total number of vertices                                 c
c ICONTN   = num. of vertices of the original region that are outside c
c            the truncated region                                     c
c ICONTP   = num. of vertices of the original region that remain in   c
c            the truncated region                                     c
c VERTP0   = vertex coordinates of the truncated polyhedron           c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE INTE3D(C,ICONTN,ICONTP,IPV0,NIPV0,NTP0,NTS0,NTV0,
     -     VERTP0,XNC,XNS0,YNC,YNS0,ZNC,ZNS0)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original/truncated polyhedron
      DIMENSION IPV0(NS,NV),NIPV0(NS),VERTP0(NV,3),XNS0(NS),
     -     YNS0(NS),ZNS0(NS)
c* Work variables
      DIMENSION IA(NV),IPIA0(NV),IPIA1(NV),ISCUT(NS),PHIV(NV)
      ICONTP=0
      ICONTN=0
      SUMP=0.0
      SUMN=0.0
      DO IP=1,NTP0
         IA(IP)=-1
      END DO
c* Distance function and values of IA
      DO IS=1,NTS0
         DO IV=1,NIPV0(IS)
            IP=IPV0(IS,IV)
            IF(IA(IP).EQ.(-1)) THEN
               PHIV(IP)=XNC*VERTP0(IP,1)+YNC*VERTP0(IP,2)+
     -              ZNC*VERTP0(IP,3)+C
               IF(PHIV(IP).GT.0.0) THEN
                  IA(IP)=1
                  ICONTP=ICONTP+1
                  SUMP=SUMP+PHIV(IP)
               ELSE
                  IA(IP)=0
                  ICONTN=ICONTN+1
                  SUMN=SUMN-PHIV(IP)
               END IF
            END IF
         END DO
      END DO
      IF(ICONTP.NE.0.AND.ICONTN.NE.0) THEN
c* Construction of the new polyhedron
         NTS00=NTS0
         CALL NEWPOL3D(IA,IPIA0,IPIA1,IPV0,ISCUT,NIPV0,NTP0,NTS0,
     -        NTV0,XNC,XNS0,YNC,YNS0,ZNC,ZNS0)
         IF(NTS0.LT.0) THEN
            NTS0=-NTS0
            ICONTN=ICONTN+ICONTP
            ICONTP=0
            RETURN
         END IF
c* Position of the new vertices
         DO IS=NTS00+1,NTS0
         DO IV=1,NIPV0(IS)
            IP=IPV0(IS,IV)
            IP0=IPIA0(Ip)
            IP1=IPIA1(Ip)
            VERTP0(IP,1)=VERTP0(IP0,1)-PHIV(IP0)*(VERTP0(IP1,1)-
     -           VERTP0(IP0,1))/(PHIV(IP1)-PHIV(IP0))
            VERTP0(IP,2)=VERTP0(IP0,2)-PHIV(IP0)*(VERTP0(IP1,2)-
     -           VERTP0(IP0,2))/(PHIV(IP1)-PHIV(IP0))
            VERTP0(IP,3)=VERTP0(IP0,3)-PHIV(IP0)*(VERTP0(IP1,3)-
     -           VERTP0(IP0,3))/(PHIV(IP1)-PHIV(IP0))
         END DO
c. Faces with less than 3 vertices are supressed.         
         IF(NIPV0(IS).LT.3) NIPV0(IS)=0
         END DO
      END IF      
      RETURN
      END
c--------------------------- END OF INTE3D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              TOOLV3D                                c
c---------------------------------------------------------------------c
c          This routine computes the volume of a polyhedron           c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c VERTI    = vertex coordinates of the polyhedron                     c
c XNS, ... = unit-lenght normals to the faces of the polyhedron       c
c IPV      = array containing the global indices of the polyhedron    c
c            vertices                                                 c
c NIPV     = number of vertices of each face                          c
c NTS      = total number of faces                                    c
c On return:                                                          c
c===========                                                          c
c VOL      = volume of the polyhedron                                 c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE TOOLV3D(IPV,NIPV,NTS,VERTI,VOL,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTI(NV,3),XNS(NS),YNS(NS),
     -     ZNS(NS)
      TOL=1.0D-16
      SUMS=0.0
      DO 20 IS=1,NTS
         IF(NIPV(IS).GT.0) THEN
            SUMP=0.0
            CNMAX=0.0
            IF(ABS(YNS(IS)).GE.ABS(XNS(IS)).AND.ABS(YNS(IS)).GE.
     -           ABS(ZNS(IS))) THEN
               IPROJ=2
               DNMAX=YNS(IS)
               IF(DABS(DNMAX).GT.TOL) CNMAX=VERTI(IPV(IS,1),2)+(XNS(IS)*
     -           VERTI(IPV(IS,1),1)+ZNS(IS)*VERTI(IPV(IS,1),3))/DNMAX
            ELSEIF(ABS(ZNS(IS)).GE.ABS(XNS(IS)).AND.ABS(ZNS(IS)).GE.
     -           ABS(YNS(IS))) THEN
               IPROJ=3
               DNMAX=ZNS(IS)
               IF(DABS(DNMAX).GT.TOL) CNMAX=VERTI(IPV(IS,1),3)+(XNS(IS)*
     -           VERTI(IPV(IS,1),1)+YNS(IS)*VERTI(IPV(IS,1),2))/DNMAX
            ELSE
               IPROJ=1
               DNMAX=XNS(IS)
               IF(DABS(DNMAX).GT.TOL) CNMAX=VERTI(IPV(IS,1),1)+(YNS(IS)*
     -           VERTI(IPV(IS,1),2)+ZNS(IS)*VERTI(IPV(IS,1),3))/DNMAX
            END IF
            IH=INT((NIPV(IS)-2)/2)
            DO I=2,IH+1
               IP=2*I
               IP1=IP-1
               IP2=IP-2
               IF(IPROJ.EQ.1) THEN
                  YV1=VERTI(IPV(IS,IP1),2)-VERTI(IPV(IS,1),2)
                  ZV1=VERTI(IPV(IS,IP1),3)-VERTI(IPV(IS,1),3)
                  YV2=VERTI(IPV(IS,IP),2)-VERTI(IPV(IS,IP2),2)
                  ZV2=VERTI(IPV(IS,IP),3)-VERTI(IPV(IS,IP2),3)
                  SUMP=SUMP+YV1*ZV2-ZV1*YV2                  
               ELSEIF(IPROJ.EQ.2) THEN
                  XV1=VERTI(IPV(IS,IP1),1)-VERTI(IPV(IS,1),1)
                  ZV1=VERTI(IPV(IS,IP1),3)-VERTI(IPV(IS,1),3)
                  XV2=VERTI(IPV(IS,IP),1)-VERTI(IPV(IS,IP2),1)
                  ZV2=VERTI(IPV(IS,IP),3)-VERTI(IPV(IS,IP2),3)
                  SUMP=SUMP+ZV1*XV2-XV1*ZV2
               ELSE
                  XV1=VERTI(IPV(IS,IP1),1)-VERTI(IPV(IS,1),1)
                  YV1=VERTI(IPV(IS,IP1),2)-VERTI(IPV(IS,1),2)
                  XV2=VERTI(IPV(IS,IP),1)-VERTI(IPV(IS,IP2),1)
                  YV2=VERTI(IPV(IS,IP),2)-VERTI(IPV(IS,IP2),2)
                  SUMP=SUMP+XV1*YV2-YV1*XV2
               END IF
            END DO
            IF(2*(IH+1).LT.NIPV(IS)) THEN
               IF(IPROJ.EQ.1) THEN
                  YV1=VERTI(IPV(IS,NIPV(IS)),2)-VERTI(IPV(IS,1),2)
                  ZV1=VERTI(IPV(IS,NIPV(IS)),3)-VERTI(IPV(IS,1),3)
                  YV2=VERTI(IPV(IS,1),2)-VERTI(IPV(IS,NIPV(IS)-1),2)
                  ZV2=VERTI(IPV(IS,1),3)-VERTI(IPV(IS,NIPV(IS)-1),3)
                  SUMP=SUMP+YV1*ZV2-ZV1*YV2
               ELSEIF(IPROJ.EQ.2) THEN
                  XV1=VERTI(IPV(IS,NIPV(IS)),1)-VERTI(IPV(IS,1),1)
                  ZV1=VERTI(IPV(IS,NIPV(IS)),3)-VERTI(IPV(IS,1),3)
                  XV2=VERTI(IPV(IS,1),1)-VERTI(IPV(IS,NIPV(IS)-1),1)
                  ZV2=VERTI(IPV(IS,1),3)-VERTI(IPV(IS,NIPV(IS)-1),3)
                  SUMP=SUMP+ZV1*XV2-XV1*ZV2
               ELSE
                  XV1=VERTI(IPV(IS,NIPV(IS)),1)-VERTI(IPV(IS,1),1)
                  YV1=VERTI(IPV(IS,NIPV(IS)),2)-VERTI(IPV(IS,1),2)
                  XV2=VERTI(IPV(IS,1),1)-VERTI(IPV(IS,NIPV(IS)-1),1)
                  YV2=VERTI(IPV(IS,1),2)-VERTI(IPV(IS,NIPV(IS)-1),2)
                  SUMP=SUMP+XV1*YV2-YV1*XV2
               END IF
            ENDIF
            IF(DABS(DNMAX).GT.TOL) SUMS=SUMS+CNMAX*SUMP
         END IF
 20   CONTINUE
      VOL=SUMS/6.0
      RETURN
      END
c-------------------------- END OF TOOLV3D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              CPPOL3D                                c
c---------------------------------------------------------------------c
c      This version dated: April 10, 2007 by Lopez and Hernandez      c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c CS0      = constants of the planes containing the faces of the      c
c            original polyhedron                                      c
c XNS0, ...= unit-lenght normals to the faces of the original pol.    c
c IPV0     = array containing the global indices of the original pol. c
c            vertices                                                 c
c NIPV0    = number of vertices of each face                          c
c NTS0     = total number of faces                                    c
c NTP0     = last global vertex index                                 c
c NTV0     = total number of vertices                                 c
c VERTI0   = vertex coordinates of the original polyhedron            c
c On return:                                                          c
c===========                                                          c
c CS       = constants of the planes containing the faces of the      c
c            copied polyhedron                                        c
c XNS,  ...= unit-lenght normals to the faces of the copied pol.      c
c IPV      = array containing the global indices of the copied pol.   c
c            vertices                                                 c
c NIPV     = number of vertices of each face                          c
c NTS      = total number of faces                                    c
c NTP      = last global vertex index                                 c
c NTV      = total number of vertices                                 c
c VERTI    = vertex coordinates of the copied polyhedron              c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE CPPOL3D(CS,CS0,IPV,IPV0,NIPV,NIPV0,NTP,NTP0,NTS,NTS0,
     -     NTV,NTV0,VERTI,VERTI0,XNS,XNS0,YNS,YNS0,ZNS,ZNS0)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CS(NS),CS0(NS),IPV(NS,NV),IPV0(NS,NV),NIPV(NS),
     -     NIPV0(NS),VERTI0(NV,3),VERTI(NV,3),XNS(NS),XNS0(NS),YNS(NS),
     -     YNS0(NS),ZNS(NS),ZNS0(NS)
      NTS=NTS0
      NTV=NTV0
      NTP=NTP0
      DO 10 IP=1,NTP0
         DO 20 J=1,3
            VERTI(IP,J)=VERTI0(IP,J)
 20      CONTINUE
 10   CONTINUE
      DO 30 I=1,NTS0
         XNS(I)=XNS0(I)
         YNS(I)=YNS0(I)
         ZNS(I)=ZNS0(I)
         NIPV(I)=NIPV0(I)
         CS(I)=CS0(I)
         DO 40 J=1,NIPV0(I)
            IPV(I,J)=IPV0(I,J)
 40      CONTINUE
 30   CONTINUE
      RETURN
      END
c-------------------------- END OF CPPOL3D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            RESTORE3D                                c
c---------------------------------------------------------------------c
c          This routine restores the structure a polyhedron           c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c CS       = constants of the planes containing the faces of the      c
c            original polyhedron                                      c
c XNS, ... = unit-lenght normals to the faces of the original pol.    c
c IPV      = array containing the global indices of the original pol. c
c            vertices                                                 c
c NIPV     = number of vertices of each face                          c
c NTS      = total number of faces                                    c
c NTP      = last global vertex index                                 c
c NTV      = total number of vertices                                 c
c VERTI    = vertex coordinates of the original polyhedron            c
c On return:                                                          c
c===========                                                          c
c CS       = constants of the planes containing the faces of the      c
c            restored polyhedron                                      c
c XNS,  ...= unit-lenght normals to the faces of the restored pol.    c
c IPV      = array containing the global indices of the restored pol. c
c            vertices                                                 c
c NIPV     = number of vertices of each face                          c
c NTS      = total number of faces                                    c
c NTP      = last global vertex index                                 c
c NTV      = total number of vertices                                 c
c VERTI    = vertex coordinates of the restored polyhedron            c
c---------------------------------------------------------------------c 
c---------------------------------------------------------------------c 
      SUBROUTINE RESTORE3D(CS,IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polyhedron
      DIMENSION CS(NS),IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),YNS(NS),
     -     ZNS(NS)
c* Work polyhedron 0
      DIMENSION CS0(NS),IPV0(NS,NV),NIPV0(NS),VERTP0(NV,3),XNS0(NS),
     -     YNS0(NS),ZNS0(NS)
      DIMENSION IPT0(NV)
c* Obtain the work polyhedron 
      CALL CPPOL3D(CS0,CS,IPV0,IPV,NIPV0,NIPV,NTP0,NTP,NTS0,NTS,NTV0,
     -     NTV,VERTP0,VERTP,XNS0,XNS,YNS0,YNS,ZNS0,ZNS)
c* In each face, consecutive vertices with the same vector position 
c* are eliminated. We use the tolerance TOLP
      TOLP=1.0D-16
      DO IS=1,NTS0
         IPV(IS,1)=IPV0(IS,1)
         IVT=0
         DO IV=1,NIPV0(IS)
            IP=IPV0(IS,IV)
            IV0=IV-1
            IF(IV0.EQ.0) IV0=NIPV0(IS)
            IP0=IPV0(IS,IV0)
            DMOD=((VERTP0(IP,1)-VERTP0(IP0,1))**2+(VERTP0(IP,2)-
     -           VERTP0(IP0,2))**2+(VERTP0(IP,3)-
     -           VERTP0(IP0,3))**2)**0.5
            IF(DMOD.GT.TOLP) THEN
               IVT=IVT+1
               IPV(IS,IVT)=IPV0(IS,IV)
            END IF
         END DO
         NIPV(IS)=IVT
      END DO
      CALL CPPOL3D(CS0,CS,IPV0,IPV,NIPV0,NIPV,NTP0,NTP,NTS0,NTS,NTV0,
     -     NTV,VERTP0,VERTP,XNS0,XNS,YNS0,YNS,ZNS0,ZNS)
c* Eliminate faces with zero or only one vertex
      NTS=0
      DO IS=1,NTS0
         IF(NIPV0(IS).GT.1) THEN
            NTS=NTS+1
            NIPV(NTS)=NIPV0(IS)
            DO IV=1,NIPV0(IS)
               IPV(NTS,IV)=IPV0(IS,IV)
               CS(NTS)=CS0(IS)
               XNS(NTS)=XNS0(IS)
               YNS(NTS)=YNS0(IS)
               ZNS(NTS)=ZNS0(IS)
            END DO
         END IF
      END DO
c* Link coincident vertices of different faces
      DO IP1=1,NTP-1
         DO IP2=IP1+1,NTP
            DMOD=((VERTP(IP1,1)-VERTP(IP2,1))**2+(VERTP(IP1,2)-
     -           VERTP(IP2,2))**2+(VERTP(IP1,3)-
     -           VERTP(IP2,3))**2)**0.5
            IF(DMOD.LE.TOLP) THEN
               DO IS=1,NTS
                  DO IV=1,NIPV(IS)
                     IF(IPV(IS,IV).EQ.IP2) IPV(IS,IV)=IP1
                  END DO
               END DO
            END IF
         END DO
      END DO
c* Renumber consecutively all the vertex indices making NTP=NTV
      IPT=0
      DO IP=1,NTP
         IC=0
         DO IS=1,NTS
            DO IV=1,NIPV(IS)
               IF(IPV(IS,IV).EQ.IP.AND.IC.EQ.0) THEN
                  IC=1
                  IPT=IPT+1
                  IPT0(IP)=IPT
               END IF
            END DO
         END DO         
      END DO
      CALL CPPOL3D(CS0,CS,IPV0,IPV,NIPV0,NIPV,NTP0,NTP,NTS0,NTS,NTV0,
     -     NTV,VERTP0,VERTP,XNS0,XNS,YNS0,YNS,ZNS0,ZNS)
      NTP=IPT
      NTV=IPT
      DO IS=1,NTS0
         DO IV=1,NIPV0(IS)
            IP=IPV0(IS,IV)
            IP1=IPT0(IP)
            IPV(IS,IV)=IP1
            VERTP(IP1,1)=VERTP0(IP,1)
            VERTP(IP1,2)=VERTP0(IP,2)
            VERTP(IP1,3)=VERTP0(IP,3)
         END DO
      END DO
      RETURN
      END
c-------------------------- END OF RESTORE3D -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              EQSOL3D                                c
c---------------------------------------------------------------------c
c        This routine solves analyticaly a cubic equation             c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c Coefficients of the equation C3·x^3+C2·x^2+C1·x+C0=0                c
c CMIN,CMAX = brackets of the solution                                c
c On return:                                                          c
c===========                                                          c
c CSOL      = solution of the cubic equation bracketed by CMIN and    c
c             CMAX                                                    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE EQSOL3D(C0,C1,C2,C3,CMIN,CMAX,CSOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TOLC=1.0D-12
      TOLC1=1.0D-12
      PI=.31415926535897932384626433832795D+01
      FSIGN=DBLE(1)
c* Eq. D*X^3+C*X^2+B*X+A=0
      D=C3
      C=C2
      B=C1
      A=C0
      IF(DABS(D).LE.TOLC.AND.DABS(C).LE.TOLC) THEN
         CSOL=-A/B
         GOTO 10
      END IF
      IF(DABS(D).LE.TOLC1) THEN
c* Eq. C*X^2+B*X+A=0
         E=B*B-DBLE(4)*C*A
         IF(E.LT.DBLE(0)) THEN
            Q=-(DBLE(1)/DBLE(2))*B
         ELSE
            Q=-(DBLE(1)/DBLE(2))*(B+DSIGN(FSIGN,B)*DSQRT(E))
         END IF
         CSOL_1=Q/C
         CSOL_2=A/Q
         DSOL_1=DABS(CSOL_1-CMIN)+DABS(CSOL_1-CMAX)
         DSOL_2=DABS(CSOL_2-CMIN)+DABS(CSOL_2-CMAX)
         IF(DSOL_1.LT.DSOL_2) THEN
            CSOL=CSOL_1
         ELSE
            CSOL=CSOL_2
         END IF
         GOTO 10
      END IF
      E=C/D
      F=B/D
      G=A/D
      E2=E*E
c* Eq. X^3+E*X^2+F*X+G=0
      Q=(E2-DBLE(3)*F)/DBLE(9)
      R=(DBLE(2)*E2*E-DBLE(9)*E*F+DBLE(27)*G)/DBLE(54)
      R2=R*R
      Q3=Q*Q*Q
      E3=E/DBLE(3)
      IF((R2).LT.(Q3)) THEN
         THETA=DACOS(R/DSQRT(Q3))
         CSOL_1=DBLE(-2)*(DSQRT(Q))*DCOS(THETA/DBLE(3))-
     -        (E3)
         CSOL_2=DBLE(-2)*(DSQRT(Q))*DCOS((THETA+DBLE(2)*
     -        PI)/DBLE(3))-(E3)
         CSOL_3=DBLE(-2)*(DSQRT(Q))*DCOS((THETA-DBLE(2)*
     -        PI)/DBLE(3))-(E3)
         DSOL_1=DABS(CSOL_1-CMIN)+DABS(CSOL_1-CMAX)
         DSOL_2=DABS(CSOL_2-CMIN)+DABS(CSOL_2-CMAX)
         DSOL_3=DABS(CSOL_3-CMIN)+DABS(CSOL_3-CMAX)
         IF(DSOL_1.LT.DSOL_2.AND.DSOL_1.LT.DSOL_3) THEN
            CSOL=CSOL_1
         ELSEIF(DSOL_2.LT.DSOL_1.AND.DSOL_2.LT.DSOL_3) THEN
            CSOL=CSOL_2
         ELSE
            CSOL=CSOL_3
         END IF
      ELSEIF((R2).EQ.(Q3)) THEN
         IF(R.GT.DBLE(0)) THEN
            CSOL_1=DBLE(-2)*DSQRT(Q)-E3
            CSOL_2=DSQRT(Q)-E3
            DSOL_1=DABS(CSOL_1-CMIN)+DABS(CSOL_1-CMAX)
            DSOL_2=DABS(CSOL_2-CMIN)+DABS(CSOL_2-CMAX)
            IF(DSOL_1.LT.DSOL_2) THEN
               CSOL=CSOL_1
            ELSE
               CSOL=CSOL_2
            END IF
         ELSEIF(R.LT.DBLE(0)) THEN
            CSOL_1=DBLE(2)*DSQRT(Q)-E3
            CSOL_2=-DSQRT(Q)-E3
            DSOL_1=DABS(CSOL_1-CMIN)+DABS(CSOL_1-CMAX)
            DSOL_2=DABS(CSOL_2-CMIN)+DABS(CSOL_2-CMAX)
            IF(DSOL_1.LT.DSOL_2) THEN
               CSOL=CSOL_1
            ELSE
               CSOL=CSOL_2
            END IF
         ELSE
            CSOL=-E3
         END IF
      ELSE
         P=-DSIGN(FSIGN,R)*(DABS(R)+DSQRT(R2-Q3))**(
     -        DBLE(1)/DBLE(3))
         IF(DABS(P).LE.TOLC) THEN
            T=DBLE(0)
         ELSE
            T=Q/P
         END IF
         CSOL=(P+T)-(E3)
      END IF
 10   CONTINUE
      RETURN
      END
c-------------------------- END OF EQSOL3D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              NEWTON3D                               c
c---------------------------------------------------------------------c
c   This routine solves a cubic eq. using the Newton-Raphson method   c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c Coefficients of the equation D·x^3+C·x^2+B·x+A=0                    c
c CMIN,CMAX = brackets of the solution                                c
c On return:                                                          c
c===========                                                          c
c CSOL      = solution of the cubic equation bracketed by CMIN and    c
c             CMAX                                                    c
c isol      = 0 si la solucion se encuentra entre CMIN y CMAX,        c
c            -1 si una posible sol se encuentra a la izda del bracket,c         c            +1 si una posible sol se encuentra a la dcha del bracket.c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE NEWTON3D(A,B,C,D,CMIN,CMAX,CSOL,ISOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ITEMAX=100
      ISOL=0
      TOLC=1.0D-14
      TOLD=1.0D-20
      CSOL=(CMIN+CMAX)/DBLE(2)
      CL=CMIN
      CH=CMAX
      COLD=CSOL
      CSOL2=CSOL*CSOL
      CMIN2=CMIN*CMIN
      CMAX2=CMAX*CMAX
      FUNC=D*CSOL2*CSOL+C*CSOL2+B*CSOL+A
      DFUNC=DBLE(3)*D*CSOL2+DBLE(2)*C*CSOL+B
      FUNCMIN=D*CMIN2*CMIN+C*CMIN2+B*CMIN+A
      FUNCMAX=D*CMAX2*CMAX+C*CMAX2+B*CMAX+A
      IF(FUNCMIN*FUNCMAX.GT.DBLE(0)) THEN
         IF(DFUNC*FUNCMIN.GT.DBLE(0)) THEN
            ISOL=-1
            CSOL=CMIN
         ELSE
            ISOL=+1
            CSOL=CMAX
         ENDIF
         RETURN
      ENDIF
      DO JITER=1,ITEMAX
         IF(DABS(DFUNC).GT.TOLD) DCSOL=FUNC/DFUNC
c* Use a simple bisection when the Newton-Raphson iteration sends the
c* solution out of bounds or when DFUNC=0.0
         IF(DABS(DFUNC).LE.TOLD.OR.(CL-CSOL+DCSOL)*(CSOL-DCSOL-CH).
     -        LT.DBLE(0)) THEN
            DCSOL=(CH-CL)/DBLE(2)
            CSOL=CL+DCSOL
            IF(DABS(CSOL-COLD).LT.TOLC.AND.JITER.GT.1) then
               RETURN
            END IF
         ELSE
            CSOL=CSOL-DCSOL
            IF(DABS(DCSOL).LT.TOLC) RETURN
         END IF
         CSOL2=CSOL*CSOL
         FUNC=D*CSOL2*CSOL+C*CSOL2+B*CSOL+A
         DFUNC=DBLE(3)*D*CSOL2+DBLE(2)*C*CSOL+B
         IF(FUNC*FUNCMIN.GT.DBLE(0)) THEN
            CL=CSOL
         ELSE
            CH=CSOL
         END IF
         COLD=CSOL
      END DO
      RETURN
      END
c-------------------------- END OF NEWTON3D --------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                               DIST3D                                c
c---------------------------------------------------------------------c
c   This routine computes the distance from a point to a polygon      c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c N         = number of vertices of the polygon                       c
c X,Y,Z     = vertex coordinates of the polygon                       c
c XP,YP,ZP  = coordinates of the point                                c
c On return:                                                          c
c===========                                                          c
c D         = exact distance from the point (XP,YP,ZP) to the polygon c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE DIST3D(D,N,X,Y,Z,XP,YP,ZP) 
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Polygon
      DIMENSION X(NV),Y(NV),Z(NV)
c* Working variables and arrays
      DIMENSION CI(NV),PHI(NV),XNI(NV),XNT(NV),YNI(NV),YNT(NV),ZNI(NV),
     -     ZNT(NV)
      LOGICAL INSIDE
c*
c*                  ^N                            
c*            o-----|---------o                                  
c*           /      |         |                     
c*          /      _|         |                     
c*         /      | |         |                     
c*        /                   o                    
c*   I-1 o                   /                     
c*        \        ^NI      /                     
c*         \       |       /                         
c*          o-------------o                       
c*         I      -->     I+1     
c*                NTI
c*                     
c* Obtain the eq. X·N+C=0 that defines the plane containing the polygon
      TOL=1.0D-20
      VMOD=0.0
      DO I=1,N-2
         I2=I+1
         I3=I2+1
         XV1=X(I2)-X(I)
         YV1=Y(I2)-Y(I)
         ZV1=Z(I2)-Z(I)
         XV2=X(I3)-X(I2)
         YV2=Y(I3)-Y(I2)
         ZV2=Z(I3)-Z(I2)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         VMOD=(XN**2+YN**2+ZN**2)**0.5
         IF(VMOD.GT.0.0) GOTO 10
      END DO
 10   CONTINUE
      IF(VMOD.LT.TOL) THEN
         D=((XP-X(1))**2+(YP-Y(1))**2+(ZP-Z(1))**2)**0.5
         RETURN
      ENDIF
      XN=XN/VMOD
      YN=YN/VMOD
      ZN=ZN/VMOD
      C=-1.0*(XN*X(1)+YN*Y(1)+ZN*Z(1))
c* Compute the edges normal NT, NI = N x NT and the distance PHI from each
c* vertex I of the polygon to the plane defined as X·NI+CI=0
      DO I=1,N
         I2=I+1
         IF(I.EQ.N) I2=1
         XT=X(I2)-X(I)
         YT=Y(I2)-Y(I)
         ZT=Z(I2)-Z(I)
         TMOD=(XT**2+YT**2+ZT**2)**0.5
         IF(TMOD.LT.TOL) THEN
            PHI(I)=0.0
         ELSE
            XNT(I)=XT/TMOD
            YNT(I)=YT/TMOD
            ZNT(I)=ZT/TMOD
            XNI(I)=YN*ZNT(I)-ZN*YNT(I)
            YNI(I)=ZN*XNT(I)-XN*ZNT(I)
            ZNI(I)=XN*YNT(I)-YN*XNT(I)
            CI(I)=-1.0*(XNI(I)*X(I)+YNI(I)*Y(I)+ZNI(I)*Z(I))
            PHI(I)=XNI(I)*XP+YNI(I)*YP+ZNI(I)*ZP+CI(I)
         END IF
      END DO
      INSIDE=.TRUE.
c* Init loop
      DO I=1,N
         I2=I+1
         IF(I.EQ.N) I2=1
         IF(PHI(I).LT.0.0) THEN
            INSIDE=.FALSE.
            C1=-1.0*(XNT(I)*X(I)+YNT(I)*Y(I)+ZNT(I)*Z(I))
            C2=1.0*(XNT(I)*X(I2)+YNT(I)*Y(I2)+ZNT(I)*Z(I2))
            PHI1=XNT(I)*XP+YNT(I)*YP+ZNT(I)*ZP+C1
            PHI2=-1.0*(XNT(I)*XP+YNT(I)*YP+ZNT(I)*ZP)+C2
            IF(PHI1.GE.0.0.AND.PHI2.GE.0.0) THEN
               T0=XNT(I)*(XP-X(I))+YNT(I)*(YP-Y(I))+ZNT(I)*(ZP-Z(I))
               XQ=X(I)+T0*XNT(I)
               YQ=Y(I)+T0*YNT(I)
               ZQ=Z(I)+T0*ZNT(I)
               D=((XP-XQ)**2+(YP-YQ)**2+(ZP-ZQ)**2)**0.5
               RETURN
            ELSEIF(PHI1.LE.0.0) THEN
               D=((XP-X(I))**2+(YP-Y(I))**2+(ZP-Z(I))**2)**0.5
               IF(I.NE.1) RETURN
            ELSEIF(PHI2.LE.0.0.AND.PHI(I2).GE.0.0) THEN
               D=((XP-X(I2))**2+(YP-Y(I2))**2+(ZP-Z(I2))**2)**0.5
               RETURN
            END IF
         END IF
      END DO
      IF(INSIDE) D=DABS(XN*XP+YN*YP+ZN*ZP+C)
      RETURN
      END
c--------------------------- END OF DIST3D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              INITF3D                                c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c FUNC3D   = external user-supplied function where the interface      c
c            shape is analytically defined                            c
c IPV      = array containing the global indices of the original pol. c
c            vertices                                                 c
c NC       = number of sub-cells along each coordinate axis of the    c
c            superimposed Cartesian grid                              c
c NIPV     = number of vertices of each face                          c
c NTP      = last global vertex index                                 c
c NTS      = total number of faces                                    c
c NTV      = total number of vertices                                 c
c TOL      = prescribed positive tolerance for the distance to the    c
c            interface                                                c
c VERTP    = vertex coordinates of the original polyhedron            c
c XNS, ... = unit-lenght normals to the faces of the original polyh.  c
c On return:                                                          c
c===========                                                          c
c VF       = material volume fraction                                 c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE INITF3D(FUNC3D,IPV,NC,NIPV,NTP,NTS,NTV,TOL,VERTP,VF,
     -     XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polyhedron
      DIMENSION CS(NS),IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),YNS(NS),
     -     ZNS(NS)
c* Work polyhedron 0
      DIMENSION CS0(NS),IPV0(NS,NV),NIPV0(NS),VERTP0(NV,3),XNS0(NS),
     -     YNS0(NS),ZNS0(NS)
c* Work polyhedron 1
      DIMENSION CS1(NS),IPV1(NS,NV),NIPV1(NS),VERTP1(NV,3),XNS1(NS),
     -     YNS1(NS),ZNS1(NS)
c* Work polyhedron 2
      DIMENSION CS2(NS),IPV2(NS,NV),NIPV2(NS),VERTP2(NV,3),XNS2(NS),
     -     YNS2(NS),ZNS2(NS)
c* Work variables
      DIMENSION IA(NV),ICHECK(NV),IPIA0(NV),IPIA1(NV),ISCUT(NS),PHIV(NV)
c* External function where the interface shape is analytically defined
      EXTERNAL FUNC3D
c.. Coordinate extremes of the cell and vertex tagging
      XMIN=1.0D+20
      XMAX=-1.0D+20
      YMIN=1.0D+20
      YMAX=-1.0D+20
      ZMIN=1.0D+20
      ZMAX=-1.0D+20
      ICONTP=0
      ICONTN=0
      DO IP=1,NTP
         ICHECK(IP)=0
      END DO
      DO IS=1,NTS
         DO IV=1,NIPV(IS)
            IP=IPV(IS,IV)
            IF(ICHECK(IP).EQ.0) THEN
               ICHECK(IP)=1
               XP=VERTP(IP,1)
               YP=VERTP(IP,2)
               ZP=VERTP(IP,3)
               XMIN=DMIN1(XMIN,XP)
               XMAX=DMAX1(XMAX,XP)
               YMIN=DMIN1(YMIN,YP)
               YMAX=DMAX1(YMAX,YP)
               ZMIN=DMIN1(ZMIN,ZP)
               ZMAX=DMAX1(ZMAX,ZP)
               PHIV(IP)=FUNC3D(XP,YP,ZP)
               IF(PHIV(IP).GE.0.0) THEN
                  IA(IP)=1
                  ICONTP=ICONTP+1
               ELSE
                  IA(IP)=0
                  ICONTN=ICONTN+1
               END IF
            END IF
         END DO
      END DO
      DX=XMAX-XMIN
      DY=YMAX-YMIN
      DZ=ZMAX-ZMIN
c.. initialization
      IPHI=0
      PHIMIN=10.0*DMAX1(DX,DY,DZ)
      DO I=1,NTV
         PHIMIN=DMIN1(PHIMIN,DABS(PHIV(I)))
      END DO
      IF(PHIMIN.LT.TOL*DX) IPHI=1
      IF(IPHI.EQ.0) THEN
         IF(ICONTP.EQ.NTV) THEN
            VF=1.0
            RETURN
         END IF
         IF(ICONTN.EQ.NTV) THEN
            VF=0.0
            RETURN
         END IF
      END IF         
c.. compute the volume VOLT of the original polyhedron
      CALL TOOLV3D(IPV,NIPV,NTS,VERTP,VOLT,XNS,YNS,ZNS)
      DDX=DX/NC
      DDY=DY/NC
      DDZ=DZ/NC
      VF=0.0
      DO IC=1,NC
         XC=XMIN+(IC-1)*DDX
         CALL CPPOL3D(CS2,CS,IPV2,IPV,NIPV2,NIPV,NTP2,NTP,NTS2,
     -        NTS,NTV2,NTV,VERTP2,VERTP,XNS2,XNS,YNS2,YNS,ZNS2,ZNS)
         CX1=-XC
         IF(IC.GT.1) CALL INTE3D(CX1,ICONTN,ICONTP,IPV2,NIPV2,NTP2,NTS2,
     -        NTV2,VERTP2,1.0D0,XNS2,0.0D0,YNS2,0.0D0,ZNS2)
         CX2=XC+DDX
         CALL INTE3D(CX2,ICONTN,ICONTP,IPV2,NIPV2,NTP2,NTS2,NTV2,
     -        VERTP2,-1.0D0,XNS2,0.0D0,YNS2,0.0D0,ZNS2)
         DO JC=1,NC
            YC=YMIN+(JC-1)*DDY
            CALL CPPOL3D(CS1,CS2,IPV1,IPV2,NIPV1,NIPV2,NTP1,NTP2,
     -           NTS1,NTS2,NTV1,NTV2,VERTP1,VERTP2,XNS1,XNS2,YNS1,
     -           YNS2,ZNS1,ZNS2)
            CY1=-YC
            IF(JC.GT.1) CALL INTE3D(CY1,ICONTN,ICONTP,IPV1,NIPV1,NTP1,
     -           NTS1,NTV1,VERTP1,0.0D0,XNS1,1.0D0,YNS1,0.0D0,ZNS1)
            IF(ICONTP.NE.0.OR.JC.EQ.1) THEN
               CY2=YC+DDY               
               CALL INTE3D(CY2,ICONTN,ICONTP,IPV1,NIPV1,NTP1,NTS1,
     -              NTV1,VERTP1,0.0D0,XNS1,-1.0D0,YNS1,0.0D0,ZNS1)
               IF(ICONTP.NE.0) THEN
                  DO KC=1,NC
                     ZC=ZMIN+(KC-1)*DDZ
                     CALL CPPOL3D(CS0,CS1,IPV0,IPV1,NIPV0,NIPV1,
     -                    NTP0,NTP1,NTS0,NTS1,NTV0,NTV1,VERTP0,
     -                    VERTP1,XNS0,XNS1,YNS0,YNS1,ZNS0,ZNS1)
                     CZ1=-ZC
                     IF(KC.GT.1) CALL INTE3D(CZ1,ICONTN,ICONTP,IPV0,
     -                    NIPV0,NTP0,NTS0,NTV0,VERTP0,0.0D0,XNS0,0.0D0,
     -                    YNS0,1.0D0,ZNS0)
                     IF(ICONTP.NE.0.OR.KC.EQ.1) THEN
                        CZ2=ZC+DDZ
                        CALL INTE3D(CZ2,ICONTN,ICONTP,IPV0,NIPV0,
     -                       NTP0,NTS0,NTV0,VERTP0,0.0D0,XNS0,
     -                       0.0D0,YNS0,-1.0D0,ZNS0)
                        IF(ICONTP.NE.0) THEN
c..   Subcell dedtermination by truncation
                           ICONTP=0
                           ICONTN=0
                           DO IP=1,NTP0
                              ICHECK(IP)=0
                           END DO
                           DO IS=1,NTS0
                              DO IV=1,NIPV0(IS)
                                 IP=IPV0(IS,IV)
                                 IF(ICHECK(IP).EQ.0) THEN
                                    ICHECK(IP)=1
                                    X=VERTP0(IP,1)
                                    Y=VERTP0(IP,2)
                                    Z=VERTP0(IP,3)
                                    PHIV(IP)=FUNC3D(X,Y,Z)
                                    IF(PHIV(IP).GE.0.0) THEN
                                       IA(IP)=1
                                       ICONTP=ICONTP+1
                                    ELSE
                                       IA(IP)=0
                                       ICONTN=ICONTN+1
                                    END IF
                                 END IF
                              END DO
                           END DO                        
                           IF(ICONTN.EQ.0) THEN
                              CALL TOOLV3D(IPV0,NIPV0,NTS0,VERTP0,
     -                             VOLF,XNS0,YNS0,ZNS0)
                              VF=VF+VOLF                           
                           ELSEIF(ICONTN.GT.0.AND.ICONTP.GT.0)THEN
                              NTSINI=NTS0
                              CALL NEWPOL3D(IA,IPIA0,IPIA1,IPV0,
     -                             ISCUT,NIPV0,NTP0,NTS0,NTV0,
     -                             1.0d0,XNS0,0.0d0,YNS0,0.0d0,
     -                             ZNS0)
c.. Location of the new intersection points
                              IF(NTS0.GT.NTSINI) THEN
                                 IS=NTS0
                                 IS2=NTS0
                                 DO IS=NTSINI+1,NTS0
                                    SUMX=0.0
                                    SUMY=0.0
                                    SUMZ=0.0
                                    DO IV=1,NIPV0(IS)
                                       IP=IPV0(IS,IV)
                                       IP0=IPIA0(Ip)
                                       IP1=IPIA1(Ip)
                                       VERTP0(IP,1)=VERTP0(IP0,1)-
     -                                      PHIV(IP0)*(VERTP0(IP1,
     -                                      1)-VERTP0(IP0,1))/(
     -                                      PHIV(IP1)-PHIV(IP0))
                                       VERTP0(IP,2)=VERTP0(IP0,2)-
     -                                      PHIV(IP0)*(VERTP0(IP1,
     -                                      2)-VERTP0(IP0,2))/(
     -                                      PHIV(IP1)-PHIV(IP0))
                                       VERTP0(IP,3)=VERTP0(IP0,3)-
     -                                      PHIV(IP0)*(VERTP0(IP1,
     -                                      3)-VERTP0(IP0,3))/(
     -                                      PHIV(IP1)-PHIV(IP0))
                                       SUMX=SUMX+VERTP0(IP,1)
                                       SUMY=SUMY+VERTP0(IP,2)
                                       SUMZ=SUMZ+VERTP0(IP,3)
                                    END DO
                                    NTP0=NTP0+1
                                    VERTP0(NTP0,1)=SUMX/NIPV0(IS)
                                    VERTP0(NTP0,2)=SUMY/NIPV0(IS)
                                    VERTP0(NTP0,3)=SUMZ/NIPV0(IS)
c     : The new face IS is replaced by NIPV(IS) triangular faces
                                    DO IV=1,NIPV0(IS)
                                       IS2=IS2+1
                                       IV2=IV+1
                                       IF(IV2.GT.
     -                                      NIPV0(IS)) IV2=1
                                       NIPV0(IS2)=3
                                       IPV0(IS2,1)=NTP0
                                       IPV0(IS2,2)=IPV0(IS,IV)
                                       IPV0(IS2,3)=IPV0(IS,IV2)
                                       XV1=VERTP0(IPV0(IS2,2),1)-
     -                                      VERTP0(IPV0(IS2,1),1)
                                       YV1=VERTP0(IPV0(IS2,2),2)-
     -                                      VERTP0(IPV0(IS2,1),2)
                                       ZV1=VERTP0(IPV0(IS2,2),3)-
     -                                      VERTP0(IPV0(IS2,1),3)
                                       XV2=VERTP0(IPV0(IS2,3),1)-
     -                                      VERTP0(IPV0(IS2,2),1)
                                       YV2=VERTP0(IPV0(IS2,3),2)-
     -                                      VERTP0(IPV0(IS2,2),2)
                                       ZV2=VERTP0(IPV0(IS2,3),3)-
     -                                      VERTP0(IPV0(IS2,2),3)
                                       XM=YV1*ZV2-ZV1*YV2
                                       YM=ZV1*XV2-XV1*ZV2
                                       ZM=XV1*YV2-YV1*XV2
                                       AMOD=(XM**2.0+YM**2.0+
     -                                      ZM**2.0)**0.5
                                       IF(AMOD.NE.0.0) THEN
                                          XNS0(IS2)=XM/AMOD
                                          YNS0(IS2)=YM/AMOD
                                          ZNS0(IS2)=ZM/AMOD
                                       ELSE
                                          NIPV0(IS2)=0
                                       END IF
                                    END DO
c     * Cancel the IS face
                                    IF(IS2.GT.IS) NIPV0(IS)=0
                                 END DO
                                 NTS0=IS2
                              end if
                              CALL TOOLV3D(IPV0,NIPV0,NTS0,VERTP0,
     -                             VOLF,XNS0,YNS0,ZNS0)
                              VF=VF+VOLF
                           END IF                                 
                        END IF
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END DO
      VF=VF/VOLT
      RETURN
      END
c--------------------------- END OF INITF3D --------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              POLOUT3D                               c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c IFILE    = number # used to name the external VTK file              c
c IPV      = array containing the global indices of the polyhedron    c
c            vertices                                                 c
c NIPV     = number of vertices of each face                          c
c NTS      = total number of faces                                    c
c NTP      = last global vertex index (note that if the polyhedron    c
c            is not previously truncated, then NTP=NTV)               c
c VERTP    = vertex coordinates of the polyhedron                     c
c On return:                                                          c
c===========                                                          c
c pol#.vtk = external VTK file                                        c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE POLOUT3D(IFILE,IPV,NIPV,NTP,NTS,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polyhedron
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3)
      CHARACTER(LEN=12)FILENAME
      WRITE(FILENAME,'("pol",I5.5,".vtk")')IFILE
      OPEN(11, FILE=FILENAME)
      WRITE(11,'(A26)')'# vtk DataFile Version 2.0'
      WRITE(11,'(A6,I5.5)')'File: ', IFILE
      WRITE(11,'(A5)')'ASCII'
      WRITE(11,*)' '
      WRITE(11,'(A16)')'DATASET POLYDATA'
      WRITE(11,'(A6,I7,A6)')'POINTS',NTP,' float'
      DO IP=1,NTP
         WRITE(11,'(3F12.6)')VERTP(IP,1),VERTP(IP,2),VERTP(IP,3)
      END DO
      NPOLY=0
      NDATA=0
      DO IS=1,NTS
         IF(NIPV(IS).GT.0) THEN
            NPOLY=NPOLY+1
            NDATA=NDATA+NIPV(IS)+1
         END IF
      END DO
      WRITE(11,'(A8,I7,I7)')'POLYGONS',NPOLY,NDATA
      DO IS=1,NTS
         IF(NIPV(IS).GT.0) THEN
            WRITE(11,'(I7)')NIPV(IS)
            DO IV=1,NIPV(IS)
               WRITE(11,'(I7)')IPV(IS,IV)-1
            END DO
         END IF
      END DO
      CLOSE(11)
      RETURN
      END
c------------------------- END OF POLOUT3D ---------------------------c
c---------------------------------------------------------------------c
      
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              ENFORV2D                               c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c IPV      = array containing the global indices of the polygon       c
c            vertices                                                 c
c NTP      = last global vertex index (note that if the polygon is    c
c            not previously truncated, then NTP=NTV)                  c
c NTV      = total number of vertices                                 c
c V        = liquid volume                                            c
c VT       = total volume of the polygon                              c
c VERTP    = vertex coordinates of the polygon                        c
c XNC, ... = unit-lenght normal to the new edges on \Gamma_c          c
c On return:                                                          c
c===========                                                          c
c C        = solution of the problem                                  c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE ENFORV2D(C,IPV,NTP,NTV,V,VT,VERTP,XNC,YNC)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polygon
      DIMENSION IPV(NV),VERTP(NV,2)
c* Working polygon 0
      DIMENSION IPV0(NV),VERTP0(NV,2)
c* Working variables
      DIMENSION IA(NV),IPIA0(NV),IPIA1(NV),
     -     LISTV(NV),PHIV(NV),XNCUT(NV),YNCUT(NV)

      IF(VT.LE.0.0) THEN
         WRITE(6,*) 'THE POLYGON HAS NULL OR NEGATIVE AREA.'
         RETURN
      END IF

      TOLC=1.0D-12
c* If the polygon was previously truncated, the polygon must         *c 
c* be restored before enforcing volume conservation. Note that       *c
c* the index NTP of the last vertex inserted after the truncation    *c
c* procudure may be higher than the total number of vertices NTV of  *c
c* the truncated polygon. It should be mentioned that in the         *c 
c* context of a PLIC-VOF method the volume enforcement operations are*c
c* generally made before any intersection operations, and therefore, *c 
c* the restoration procedure is not generally required. At any case, *c
c* the corresponding RESTORE2D routine is also supplied.             *c
      IF(NTP.GT.NTV) CALL RESTORE2D(IPV,NTP,NTV,VERTP)
      VAUX=V
      LISTV(1)=1
      DO IV=1,NTV
         IA(IV)=0
         PHIV(IV)=XNC*VERTP(IV,1)+YNC*VERTP(IV,2)
c* Ordered list of global vertex indices
         DO I=1,IV-1
            IF(PHIV(IV).GT.PHIV(LISTV(I))) THEN
               DO II=IV,I+1,-1
                  LISTV(II)=LISTV(II-1)
               END DO
               LISTV(I)=IV
               GOTO 10
            END IF
         END DO
         LISTV(IV)=IV
 10      CONTINUE
      END DO
c* Begin the procedure SETIA used to obtain the brackets CMAX and CMIN
c* of the soluction C
      INVERT=0
      XNCOR=XNC
      YNCOR=YNC
      CMAX=PHIV(LISTV(1))
      CMIN=PHIV(LISTV(NTP))
      IMIN=1
      IMAX=NTP
      VMIN=DBLE(0)
      VMAX=VT
      IMAXLOLD=NTP+1
 22   CONTINUE
c* Obtain the tentative solution bracketing by interpolation
      PHIINT=PHIV(LISTV(IMIN))-(PHIV(LISTV(IMIN))-PHIV(LISTV(IMAX)))*
     -     (V-VMIN)/(VMAX-VMIN)
      IMAXL=0
      IMINL=0
      DO IP=IMIN+1,IMAX
         I=IP
         IF(PHIV(LISTV(IP)).LT.PHIINT) THEN
            IMAXL=IP
            IMINL=IP-1
            GOTO 11
         END IF
      END DO
      IF(IMAXL.EQ.0.AND.IMINL.EQ.0) THEN
         C=-PHIINT
         RETURN
      END IF      
 11   CONTINUE
      CMAX=PHIV(LISTV(IMINL))
      CMIN=PHIV(LISTV(IMAXL))
      IF((NTP-IMAXL).LT.(IMINL-1)) THEN
         INVERT=1
         CAUX=CMIN
         CMIN=-CMAX
         CMAX=-CAUX
         VAUX=VT-V
         XNC=-XNCOR
         YNC=-YNCOR
      ELSE
         INVERT=0
         VAUX=V
         XNC=XNCOR
         YNC=YNCOR
      END IF
      DO I=1,NTP
         IF(I.LE.IMINL) THEN
            IA(LISTV(I))=1-INVERT
         ELSE
            IA(LISTV(I))=INVERT
         END IF
      END DO
c* End of procedure SETIA
      CALL CPPOL2D(IPV,IPV0,NTP,NTP0,NTV,NTV0,VERTP,VERTP0)
c* Construction of the new polygon
      NTPINI=NTP0
      CALL NEWPOL2D(IA,IPIA0,IPIA1,IPV0,NTP0,NTV0,VERTP0,XNCUT,YNCUT)
      NCUT=NTP0-NTPINI
      DO IP=NTPINI+1,NTP0
         VERTP0(IP,1)=DBLE(0)
         VERTP0(IP,2)=DBLE(0)
         IA(IP)=0
      END DO
      C0=-DBLE(2)*VAUX
      C1=DBLE(0)
      C2=DBLE(0)
      C3=DBLE(0)
      IF((XNCUT(1)*YNC-YNCUT(1)*XNC).GT.DBLE(0)) THEN
         ICUTINI=1
         ICUTEND=NCUT
         ICUTNEX=1
      ELSE
         ICUTINI=NCUT
         ICUTEND=1
         ICUTNEX=-1
      END IF
      XNCS=-XNC
      YNCS=-YNC      
      DO I1=ICUTINI,ICUTEND,2*ICUTNEX
         I2=I1+ICUTNEX
         IPF1=IPIA1(I1)
         IPF2=IPIA1(I2)
         CF1=-(VERTP0(IPF1,1)*XNCUT(I1)+VERTP0(IPF1,2)*YNCUT(I1))
         CF2=-(VERTP0(IPF2,1)*XNCUT(I2)+VERTP0(IPF2,2)*YNCUT(I2))
         CUT1=XNCUT(I1)*YNCS-YNCUT(I1)*XNCS
         IF(CUT1.EQ.0.0) THEN !cut edge 1 and \Gamma_c normal are perpendicular
            BET1=DBLE(0)
         ELSE
            BET1=DBLE(1)/CUT1
         END IF
         CUT2=XNCUT(I2)*YNCS-YNCUT(I2)*XNCS
         IF(CUT2.EQ.0.0) THEN !cut edge 2 and \Gamma_c normal are perpendicular
            BET2=DBLE(0)            
         ELSE
            BET2=-DBLE(1)/CUT2
         END IF
         C2=C2+(YNCUT(I2)*XNCUT(I1)-XNCUT(I2)*YNCUT(I1))*BET1*BET2
         C1=C1-DBLE(2)*(CF2*BET2+CF1*BET1)
         C0=C0-(XNCS*VERTP0(IPF1,1)+YNCS*VERTP0(IPF1,2))*CF1*BET1-
     -        (XNCS*VERTP0(IPF2,1)+YNCS*VERTP0(IPF2,2))*CF2*BET2
      END DO
      IH=INT((NTV0-2)/2)
      IP0=IPV0(1)
      X1=VERTP0(IP0,1)*IA(IP0)
      Y1=VERTP0(IP0,2)*IA(IP0)
      DO I=2,IH+1
         IV=2*I
         IP=IPV0(IV)
         IP1=IPV0(IV-1)
         IP2=IPV0(IV-2)
         XV1=VERTP0(IP1,1)*IA(IP1)-X1
         YV1=VERTP0(IP1,2)*IA(IP1)-Y1
         XV2=VERTP0(IP,1)*IA(IP)-VERTP0(IP2,1)*IA(IP2)
         YV2=VERTP0(IP,2)*IA(IP)-VERTP0(IP2,2)*IA(IP2)
         C0=C0+XV1*YV2-YV1*XV2
      END DO
      IF(2*(IH+1).LT.NTV0) THEN
         IPN0=IPV0(NTV0-1)
         IPN=IPV0(NTV0)
         XV1=VERTP0(IPN,1)*IA(IPN)-X1
         YV1=VERTP0(IPN,2)*IA(IPN)-Y1
         XV2=X1-VERTP0(IPN0,1)*IA(IPN0)
         YV2=Y1-VERTP0(IPN0,2)*IA(IPN0)
         C0=C0+XV1*YV2-YV1*XV2
      ENDIF
      VMAXL=(C2*CMIN*CMIN+C1*
     -     CMIN+C0+DBLE(2)*VAUX)/DBLE(2)
      VMINL=(C2*CMAX*CMAX+C1*
     -     CMAX+C0+DBLE(2)*VAUX)/DBLE(2)
      
      IF(INVERT.EQ.1) THEN
         VMINLL=VT-VMAXL
         VMAXL=VT-VMINL
         VMINL=VMINLL
      END IF
      SV=(VMINL-V)*(VMAXL-V)
      IF(SV.GT.DBLE(0).AND.(IMAX-IMIN).GT.1.AND.IMAXLOLD.NE.IMAXL) THEN
         IF(VMAXL.GT.V) THEN
            VMAX=VMINL
            IMAX=IMINL
         ELSE
            VMIN=VMAXL
            IMIN=IMAXL
         END IF
         IMAXLOLD=IMAXL
         GOTO 22
      END IF
      CALL EQSOL3D(C0,C1,C2,C3,CMIN,CMAX,C)
      DSOL=DABS(C-CMIN)+DABS(C-CMAX)
      DREF=DABS(CMAX-CMIN)
      IF((DSOL/DREF).GT.(DBLE(1)+TOLC)) THEN
         CALL NEWTON3D(C0,C1,C2,C3,CMIN,CMAX,CSOLN,ISOL)
         DSOLN=DABS(CSOLN-CMIN)+DABS(CSOLN-CMAX)
         IF((DSOLN/DREF).GT.(DBLE(1)+TOLC)) C=CSOLN
      END IF
      IF(INVERT.EQ.0) THEN 
         C=-C
      ELSE
         XNC=-XNC
         YNC=-YNC
      END IF
      RETURN
      END
c----------------------    END OF ENFORV2D   -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                          ENFORV2DSZ                                 c
c... Scardovelli and Zaleski version for rectangle                    c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c DX, ...  = side lengths of the rectangle                            c
c VERTP    = vertex coordinates of the rectangle                      c
c XNC, ... = unit-lenght normal to the new face \Gamma_c              c
c V        = liquid volume                                            c
c On return:                                                          c
c===========                                                          c
c C        = solution of the problem                                  c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE ENFORV2DSZ(C,DX,DY,V,VERTP,XNC,YNC)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c.. Original polyhedron
      DIMENSION VERTP(NV,2)
      DOUBLE PRECISION M,M1
      CMIN=1.0D+14
      CMAX=-1.0D+14
      VT=DX*DY
      VBACK=V
      V=V/VT
c.. The vertex indices of the rectangular parallelepiped are supposed 
c.. to be listed from 1 to 4.
      DO 10 I=1,4
         CI=-(VERTP(I,1)*XNC+VERTP(I,2)*YNC)
         IF(CI.LE.CMIN) THEN
            CMIN=CI
            IMIN=I
         END IF
         IF(CI.GE.CMAX) THEN
            CMAX=CI
            IMAX=I
         END IF
 10   CONTINUE
c.. If the liquid volume fraction is higher than 0.5, solve the 
c.. inverse problem
      IF((VBACK/VT).LE.(DBLE(1)/DBLE(2))) THEN
         CI=CMIN
         I=IMIN
      ELSE
         CI=CMAX
         I=IMAX
         V=DBLE(1)-V
      END IF
c.. Normalize the plane equation
      SN=DABS(XNC)+DABS(YNC)
      XM=XNC/SN
      YM=YNC/SN
      XMI=XM*DX
      YMI=YM*DY
      SN=DABS(XMI)+DABS(YMI)
      XM=DABS(XMI)/SN
      YM=DABS(YMI)/SN
c.. Region limits
      M1=DMIN1(XM,YM)
      M=M1
      V1=M/(DBLE(2)*(DBLE(1)-M))
c.. Solution of the inverse problem
      IF(V.GE.DBLE(0).AND.V.LT.V1) THEN
         ALPHA=DSQRT(DBLE(2)*M*(DBLE(1)-M)*V)
      ELSE
         ALPHA=V*(DBLE(1)-M)+M/DBLE(2)
      END IF
      IF((VBACK/VT).LE.(DBLE(1)/DBLE(2))) THEN
         C=CMIN+ALPHA*DABS(CMAX-CMIN)
      ELSE
         C=CMAX-ALPHA*DABS(CMAX-CMIN)
      END IF
      V=VBACK
      RETURN
      END
c---------------------   END OF ENFORV2DSZ   -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             NEWPOL2D                                c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c IA       = 0 if \Gamma_c points out of the vertex, 1 otherwise      c
c IPV0     = array containing the global indices of the original pol. c
c            vertices                                                 c
c NTP0     = last global vertex index                                 c
c NTV0     = total number of vertices                                 c
c VERTP0   = vertex coordinates of the polygon                        c
c On return:                                                          c
c===========                                                          c
c IPV0     = array containing the global indices of the truncat. pol. c
c            vertices                                                 c
c NTP0     = last global vertex index                                 c
c NTV0     = total number of vertices                                 c
c IPIA0    = global vertex index of the original poligon with IA=0    c
c            and which is in the edge containing the intersection     c
c            point                                                    c
c IPIA1    = global vertex index of the original poligon with IA=1    c
c            and which is in the edge containing the intersection     c
c            point                                                    c
c XNCUT,...= unit-lenght normals to the edges cut by \Gamma_c         c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE NEWPOL2D(IA,IPIA0,IPIA1,IPV0,NTP0,NTV0,VERTP0,
     -     XNCUT,YNCUT)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original/truncated polygon
      DIMENSION IPV0(NV),VERTP0(NV,2)
c* Work polygon 1
      DIMENSION IPV1(NV)
c* Work variables
      DIMENSION IA(NV),IPIA0(NV),IPIA1(NV),XNCUT(NV),YNCUT(NV)
c* Determination of the cut edges
      NTV1=0
      ICUT=0
      DO IV=1,NTV0
         IP=IPV0(IV)
         IV2=IV+1
         IF(IV.EQ.NTV0) IV2=1
         IP2=IPV0(IV2)
         IF(IA(IP).EQ.1) THEN
            NTV1=NTV1+1
            IPV1(NTV1)=IPV0(IV)
         END IF
         IF(IA(IP).NE.IA(IP2)) THEN
            ICUT=ICUT+1
            NTP0=NTP0+1
            NTV1=NTV1+1
            IPV1(NTV1)=NTP0
            IA(NTP0)=0
            IF(IA(IP2).EQ.0) THEN
               IPIA0(ICUT)=IP2
               IPIA1(ICUT)=IP
            ELSE
               IPIA0(ICUT)=IP
               IPIA1(ICUT)=IP2
            END IF
            XV=VERTP0(IP2,1)-VERTP0(IP,1)
            YV=VERTP0(IP2,2)-VERTP0(IP,2)
            RMOD=(XV**2.0+YV**2.0)**0.5
            XNCUT(ICUT)=YV/RMOD
            YNCUT(ICUT)=-XV/RMOD
         END IF
      END DO
      NTV0=NTV1
      DO IV=1,NTV1
         IPV0(IV)=IPV1(IV)
      END DO
      RETURN
      END
c------------------------- END OF NEWPOL2D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                               INTE2D                                c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c C        = constant of the line containing the new edges on \Gamma_cc
c IPV0     = array containing the global indices of the original pol. c
c            vertices                                                 c
c NTP0     = last global vertex index                                 c
c NTV0     = total number of vertices                                 c
c VERTP0   = vertex coordinates of the original polygon               c
c XNC, ... = unit-lenght normal to the new edges on \Gamma_c          c
c On return:                                                          c
c===========                                                          c
c ICONTN   = num. of vertices of the original region that are outside c
c            the truncated region                                     c
c ICONTP   = num. of vertices of the original region that remain in   c
c            the truncated region                                     c
c IPV0     = array containing the global indices of the truncat. pol. c
c            vertices                                                 c
c NTP0     = last global vertex index                                 c
c NTV0     = total number of vertices                                 c
c VERTP0   = vertex coordinates of the truncated polygon              c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE INTE2D(C,ICONTN,ICONTP,IPV0,NTP0,NTV0,VERTP0,XNC,YNC)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original/truncated polyhedron
      DIMENSION IPV0(NV),VERTP0(NV,2)
c* Work variables
      DIMENSION IA(NV),IPIA0(NV),IPIA1(NV),PHIV(NV),XNCUT(NV),YNCUT(NV)
      ICONTP=0
      ICONTN=0
c* Distance function and values of IA
      DO IV=1,NTV0
         IP=IPV0(IV)
         PHIV(IP)=XNC*VERTP0(IP,1)+YNC*VERTP0(IP,2)+C
         IF(PHIV(IP).GT.0.0) THEN
            IA(IP)=1
            ICONTP=ICONTP+1
         ELSE
            IA(IP)=0
            ICONTN=ICONTN+1
         END IF
      END DO
      IF(ICONTP.NE.0.AND.ICONTN.NE.0) THEN
c* Construction of the new polygon
         NTPINI=NTP0
         CALL NEWPOL2D(IA,IPIA0,IPIA1,IPV0,NTP0,NTV0,VERTP0,XNCUT,
     -        YNCUT)
c* Position of the new vertices
         DO IP=NTPINI+1,NTP0
            IP0=IPIA0(IP-NTPINI)
            IP1=IPIA1(IP-NTPINI)            
               VERTP0(IP,1)=VERTP0(IP0,1)-PHIV(IP0)*(VERTP0(IP1,1)-
     -              VERTP0(IP0,1))/(PHIV(IP1)-PHIV(IP0))
               VERTP0(IP,2)=VERTP0(IP0,2)-PHIV(IP0)*(VERTP0(IP1,2)-
     -              VERTP0(IP0,2))/(PHIV(IP1)-PHIV(IP0))
         END DO
      END IF
      RETURN
      END
c--------------------------- END OF INTE2D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              TOOLV2D                                c
c---------------------------------------------------------------------c
c          This routine computes the volume of a polygon              c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c VERTP    = vertex coordinates of the polygon                        c
c IPV      = array containing the global indices of the polygon       c
c            vertices                                                 c
c NTV      = total number of vertices                                 c
c On return:                                                          c
c===========                                                          c
c VOL      = volume of the polygon                                    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE TOOLV2D(IPV,NTV,VERTP,VOL)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)

      SUMS=0.0
      IH=INT((ntv-2)/2)
      DO I=2,IH+1
         IP=2*I
         IP1=IP-1
         IP2=IP-2
         XV1=VERTP(IPV(IP1),1)-VERTP(IPV(1),1)
         YV1=VERTP(IPV(IP1),2)-VERTP(IPV(1),2)
         XV2=VERTP(IPV(IP),1)-VERTP(IPV(IP2),1)
         YV2=VERTP(IPV(IP),2)-VERTP(IPV(IP2),2)
         SUMS=SUMS+XV1*YV2-YV1*XV2
      END DO
      IF(2*(IH+1).LT.ntv) THEN
         XV1=VERTP(IPV(ntv),1)-VERTP(IPV(1),1)
         YV1=VERTP(IPV(ntv),2)-VERTP(IPV(1),2)
         XV2=VERTP(IPV(1),1)-VERTP(IPV(ntv-1),1)
         YV2=VERTP(IPV(1),2)-VERTP(IPV(ntv-1),2)
         SUMS=SUMS+XV1*YV2-YV1*XV2
      ENDIF
      VOL=SUMS/2.0
      
      RETURN
      END
c-------------------------- END OF TOOLV2D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              CPPOL2D                                c
c---------------------------------------------------------------------c
c         This routine copies a polygon into a new one                c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c VERTP    = vertex coordinates of the original polygon               c
c IPV      = array containing the global indices of the original      c
c            polygon vertices                                         c
c NTP      = last global vertex index (note that if the polygon is    c
c            not previously truncated, then NTP=NTV)                  c
c NTV      = total number of vertices                                 c
c On return:                                                          c
c===========                                                          c
c VERTP1   = vertex coordinates of the copied polygon                 c
c IPV1     = array containing the global indices of the copied        c
c            polygon vertices                                         c
c NTP1     = last global vertex index (note that if the polygon is    c
c            not previously truncated, then NTP=NTV)                  c
c NTV1     = total number of vertices                                 c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE CPPOL2D(IPV,IPV1,NTP,NTP1,NTV,NTV1,VERTP,VERTP1)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      DIMENSION IPV1(NV),VERTP1(NV,2)
      NTP1=NTP
      NTV1=NTV
      DO IV=1,NTV
         IP=IPV(IV)
         IPV1(IV)=IP
         VERTP1(IP,1)=VERTP(IP,1)
         VERTP1(IP,2)=VERTP(IP,2)
      END DO

      RETURN
      END
c-------------------------- END OF CPPOL2D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            RESTORE2D                                c
c---------------------------------------------------------------------c
c                 This routine restores a polygon                     c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c VERTP    = vertex coordinates of the original polygon               c
c IPV      = array containing the global indices of the original      c
c            polygon vertices                                         c
c NTP      = last global vertex index (note that if the polygon is    c
c            not previously truncated, then NTP=NTV)                  c
c NTV      = total number of vertices                                 c
c On return:                                                          c
c===========                                                          c
c VERTP    = vertex coordinates of the restored polygon               c
c IPV      = array containing the global indices of the restored      c
c            polygon vertices                                         c
c NTP      = last global vertex index (note that if the polygon is    c
c            not previously truncated, then NTP=NTV)                  c
c NTV      = total number of vertices                                 c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE RESTORE2D(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polygon
      DIMENSION IPV(NV),VERTP(NV,2)
c* Work polygon 0
      DIMENSION IPV0(NV),VERTP0(NV,2)
c* Obtain the work polygon
      CALL CPPOL2D(IPV,IPV0,NTP,NTP0,NTV,NTV0,VERTP,VERTP0)
c* Consecutive vertices with the same vector position are 
c* eliminated. We use the tolerance TOLP
      TOLP=1.0D-16
      IVT=0
      DO IV=1,NTV0
         IP=IPV0(IV)
         IV0=IV-1
         IF(IV0.EQ.0) IV0=NTV0
         IP0=IPV0(IV0)
         DMOD=((VERTP0(IP,1)-VERTP0(IP0,1))**2.0+(VERTP0(IP,2)-
     -        VERTP0(IP0,2))**2.0)**0.5
         IF(DMOD.GT.TOLP) THEN
            IVT=IVT+1
            IPV(IVT)=IVT
            VERTP(IVT,1)=VERTP0(IP,1)
            VERTP(IVT,2)=VERTP0(IP,2)
         END IF
      END DO
      NTV=IVT
      NTP=IVT
      RETURN
      END
c-------------------------- END OF RESTORE2D -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                               DIST2D                                c
c---------------------------------------------------------------------c
c   This routine computes the distance from a point to a segment      c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y       = vertex coordinates of the segment                       c
c XP,YP     = coordinates of the point                                c
c On return:                                                          c
c===========                                                          c
c D         = exact distance from the point (XP,YP) to the segment    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE DIST2D(D,X,Y,XP,YP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Segment
      DIMENSION X(2),Y(2)
      XNT=X(2)-X(1)
      YNT=Y(2)-Y(1)
      VMOD=(XNT**2.0+YNT**2.0)**0.5
      XNT=XNT/VMOD
      YNT=YNT/VMOD
      XN=-YNT
      YN=XNT
      C1=-1.0*(XNT*X(1)+YNT*Y(1))
      C2=1.0*(XNT*X(2)+YNT*Y(2))
      PHI1=XNT*XP+YNT*YP+C1
      PHI2=-XNT*XP-YNT*YP+C2
      IF(PHI1.GE.0.0.AND.PHI2.GE.0.0) THEN
         D=ABS(XN*X(1)+YN*Y(1)-(XN*XP+YN*YP))
      ELSEIF(PHI1.LE.0.0) THEN
         D=((XP-X(1))**2.0+(YP-Y(1))**2.0)**0.5
      ELSE
         D=((XP-X(2))**2.0+(YP-Y(2))**2.0)**0.5
      END IF
      RETURN
      END
c--------------------------- END OF DIST2D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              INITF2D                                c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c FUNC2D   = external user-supplied function where the interface      c
c            shape is analytically defined                            c
c IPV      = array containing the global indices of the original pol. c
c            vertices                                                 c
c NC       = number of sub-cells along each coordinate axis of the    c
c            superimposed Cartesian grid                              c
c NTP      = last global vertex index                                 c
c NTV      = total number of vertices                                 c
c TOL      = prescribed positive tolerance for the distance to the    c
c            interface                                                c
c VERTP    = vertex coordinates of the original polyhedron            c
c On return:                                                          c
c===========                                                          c
c VF       = material area fraction                                   c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE INITF2D(FUNC2D,IPV,NC,NTP,NTV,TOL,VERTP,VF)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polygon
      DIMENSION IPV(NV),VERTP(NV,2)
c* Work polygon 1
      DIMENSION IPV1(NV),VERTP1(NV,2)
c* Work polygon 2
      DIMENSION IPV2(NV),VERTP2(NV,2)
c* Work variables
      DIMENSION IA(NV),IPIA0(NV),IPIA1(NV),PHIV(NV)
      DIMENSION XNCUT(NV),YNCUT(NV)
c* External function where the interface shape is analytically defined
      EXTERNAL FUNC2D
c.. Coordinate extremes of the cell and vertex tagging
      XMIN=1.0D+20
      XMAX=-1.0D+20
      YMIN=1.0D+20
      YMAX=-1.0D+20
      ICONTP=0
      ICONTN=0
      DO IV=1,NTV
         IP=IPV(IV)
         XP=VERTP(IP,1)
         YP=VERTP(IP,2)
         XMIN=DMIN1(XMIN,XP)
         XMAX=DMAX1(XMAX,XP)
         YMIN=DMIN1(YMIN,YP)
         YMAX=DMAX1(YMAX,YP)
         PHIV(IP)=FUNC2D(XP,YP)
         IF(PHIV(IP).GE.0.0) THEN
            IA(IP)=1
            ICONTP=ICONTP+1
         ELSE
            IA(IP)=0
            ICONTN=ICONTN+1
         END IF         
      END DO
      DX=XMAX-XMIN
      DY=YMAX-YMIN
      IPHI=0
      PHIMIN=10.0*DMAX1(DX,DY)
      DO I=1,NTV
         PHIMIN=DMIN1(PHIMIN,DABS(PHIV(I)))
      END DO
      IF(PHIMIN.LT.TOL*DX) IPHI=1 
      IF(ICONTP.EQ.NTV.AND.IPHI.EQ.0) THEN
         VF=1.0
      ELSEIF(ICONTN.EQ.NTV.AND.IPHI.EQ.0) THEN
         VF=0.0
      ELSE
c.. Total volume VOLT of the original polygon
         CALL TOOLV2D(IPV,NTV,VERTP,VOLT)
         DDX=DX/NC
         DDY=DY/NC
         VF=0.0
         DO IC=1,NC
            XC=XMIN+(IC-1)*DDX
            CALL CPPOL2D(IPV,IPV2,NTP,NTP2,NTV,NTV2,VERTP,VERTP2)
            CX1=-XC
            IF(IC.GT.1) CALL INTE2D(CX1,ICONTN,ICONTP,IPV2,NTP2,NTV2,
     -           VERTP2,1.0D0,0.0D0)
            CX2=XC+DDX
            CALL INTE2D(CX2,ICONTN,ICONTP,IPV2,NTP2,NTV2,VERTP2,
     -           -1.0D0,0.0D0)
            DO JC=1,NC
               YC=YMIN+(JC-1)*DDY
               CALL CPPOL2D(IPV2,IPV1,NTP2,NTP1,NTV2,NTV1,VERTP2,
     -              VERTP1)
               CY1=-YC
               IF(JC.GT.1)CALL INTE2D(CY1,ICONTN,ICONTP,IPV1,NTP1,
     -              NTV1,VERTP1,0.0D0,1.0D0)
               IF(ICONTP.NE.0.OR.JC.EQ.1) THEN
                  CY2=YC+DDY               
                  CALL INTE2D(CY2,ICONTN,ICONTP,IPV1,NTP1,NTV1,
     -                 VERTP1,0.0D0,-1.0D0)
                  IF(ICONTP.NE.0) THEN
                     ICONTP=0
                     ICONTN=0
                     DO IV=1,NTV1
                        IP=IPV1(IV)
                        XP=VERTP1(IP,1)
                        YP=VERTP1(IP,2)
                        PHIV(IP)=FUNC2D(XP,YP)
                        IF(PHIV(IP).GE.0.0) THEN
                           IA(IP)=1
                           ICONTP=ICONTP+1
                        ELSE
                           IA(IP)=0
                           ICONTN=ICONTN+1
                        END IF         
                     END DO
                     
                     IF(ICONTN.EQ.0) THEN
                        CALL TOOLV2D(IPV1,NTV1,VERTP1,VOLF)
                        VF=VF+VOLF
                        
                     ELSEIF(ICONTN.GT.0.AND.ICONTP.GT.0) THEN
                        NTPINI=NTP1                        
                        CALL NEWPOL2D(IA,IPIA0,IPIA1,IPV1,NTP1,
     -                       NTV1,VERTP1,XNCUT,YNCUT)
c.. Location of the new intersection points
                        DO IP=NTPINI+1,NTP1
                           IP0=IPIA0(IP-NTPINI)
                           IP1=IPIA1(IP-NTPINI)            
                           VERTP1(IP,1)=VERTP1(IP0,1)-
     -                          PHIV(IP0)*(VERTP1(IP1,1)-
     -                          VERTP1(IP0,1))/(PHIV(IP1)-
     -                          PHIV(IP0))
                           VERTP1(IP,2)=VERTP1(IP0,2)-
     -                          PHIV(IP0)*(VERTP1(IP1,2)-
     -                          VERTP1(IP0,2))/(PHIV(IP1)-
     -                          PHIV(IP0))
                        END DO
                        CALL TOOLV2D(IPV1,NTV1,VERTP1,VOLF)
                        VF=VF+VOLF
                     END IF                                             
                  END IF
               END IF
            END DO
         END DO
         VF=VF/VOLT
      END IF
      RETURN
      END
c-------------------------- END OF INITF2D ---------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                              POLOUT2D                               c
c Write in an external file the vertex coordinates of the polygon in  c
c two columns format to be plotted, e.g., using GNUPLOT program       c      
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c IFILE    = number # used to name the external file                  c
c IPV      = array containing the global indices of the polygon       c
c            vertices                                                 c
c NTV      = total number of vertices                                 c
c VERTP    = vertex coordinates of the polygon                        c
c On return:                                                          c
c===========                                                          c
c pol#.out = external file                                            c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE POLOUT2D(IFILE,IPV,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c* Original polygon
      DIMENSION IPV(NV),VERTP(NV,2)
      CHARACTER(LEN=12)FILENAME
      WRITE(FILENAME,'("pol",I5.5,".out")')IFILE
      OPEN(11, FILE=FILENAME)
      WRITE(11,'(A7,I5.5)')'#File: ', IFILE
      WRITE(11,'(A9,I7)')'#VERTICES',NTV
      WRITE(11,'(A16)')'#X Y COORDINATES'
      DO IV=1,NTV
         IP=IPV(IV)
         WRITE(11,'(2F12.6)')VERTP(IP,1),VERTP(IP,2)
      END DO
      IP=IPV(1)
      WRITE(11,'(2F12.6)')VERTP(IP,1),VERTP(IP,2)
      CLOSE(11)
      RETURN
      END
c------------------------- END OF POLOUT2D ---------------------------c
c---------------------------------------------------------------------c
