c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                                MESH                                 c
c---------------------------------------------------------------------c
c            Copyright (C) 2020 J. Lopez and J. Hernandez             c
c---------------------------------------------------------------------c
c                           Mesh geometries                           c
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

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            CUBICMESH                                c
c---------------------------------------------------------------------c
c      This version dated: April 10, 2007 by Lopez and Hernandez      c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE CUBICMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
c          |               c          |                 
c          |______         c          |______          
c         /|      |        c         /      /|         
c        / |  ·3  |        c        /  ·4  / |         
c       /  |      |        c       /______/  |         
c       |·6------------    c       |      |·1|------   
c       | /      /         c       |  ·5  | /          
c       |/  ·2  /          c       |      |/           
c       /------/           c       /------/            
c      /                   c      /                    
c     /                    c     /                                
c    
c       7/----------/3
c       /|         /|
c      / |        / |
c    8/__|______4/  | 
c     |  |       |  |
c     |  /6------|--/2
c     | /        | /
c     |/_________|/
c     5           1
c
      D0=0.0D0
      D1=1.0D0

      XNS(1)=D1
      YNS(1)=D0
      ZNS(1)=D0
      XNS(2)=D0
      YNS(2)=-D1
      ZNS(2)=D0
      XNS(3)=D0
      YNS(3)=D0
      ZNS(3)=-D1
      XNS(4)=D0
      YNS(4)=D1
      ZNS(4)=D0
      XNS(5)=D0
      YNS(5)=D0
      ZNS(5)=D1
      XNS(6)=-D1
      YNS(6)=D0
      ZNS(6)=D0
      NTS=6
      NTV=8
      NTP=NTV
      NIPV(1)=4
      IPV(1,1)=1
      IPV(1,2)=2
      IPV(1,3)=3
      IPV(1,4)=4
      NIPV(2)=4
      IPV(2,1)=2
      IPV(2,2)=1
      IPV(2,3)=5
      IPV(2,4)=6
      NIPV(3)=4
      IPV(3,1)=3
      IPV(3,2)=2
      IPV(3,3)=6
      IPV(3,4)=7
      NIPV(4)=4
      IPV(4,1)=4
      IPV(4,2)=3
      IPV(4,3)=7
      IPV(4,4)=8
      NIPV(5)=4
      IPV(5,1)=1
      IPV(5,2)=4
      IPV(5,3)=8
      IPV(5,4)=5
      NIPV(6)=4
      IPV(6,1)=6
      IPV(6,2)=5
      IPV(6,3)=8
      IPV(6,4)=7
      VERTP(1,1)=D1
      VERTP(1,2)=D0
      VERTP(1,3)=D1
      VERTP(2,1)=D1
      VERTP(2,2)=D0
      VERTP(2,3)=D0
      VERTP(3,1)=D1
      VERTP(3,2)=D1
      VERTP(3,3)=D0
      VERTP(4,1)=D1
      VERTP(4,2)=D1
      VERTP(4,3)=D1
      VERTP(5,1)=D0
      VERTP(5,2)=D0
      VERTP(5,3)=D1
      VERTP(6,1)=D0
      VERTP(6,2)=D0
      VERTP(6,3)=D0
      VERTP(7,1)=D0
      VERTP(7,2)=D1
      VERTP(7,3)=D0
      VERTP(8,1)=D0
      VERTP(8,2)=D1
      VERTP(8,3)=D1
      RETURN
      END
c------------------------- END OF CUBICMESH --------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                        IRREGULAR HEXAHEDRON                         c
c---------------------------------------------------------------------c
c      This version dated: August 7, 2014 by Lopez and Hernandez      c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE HEXAHEMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CS(NS),IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
      CALL CUBICMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
c.. Cutting planes passing through the point (0,0,0)
      C=0.0D0
      XNC=1.0D0
      YNC=-0.1D0
      ZNC=-0.1D0
      DMOD=(XNC**2+YNC**2+ZNC**2)**0.5D0
      XNC=XNC/DMOD
      YNC=YNC/DMOD
      ZNC=ZNC/DMOD
      CALL INTE3D(C,ICONTN,ICONTP,IPV,NIPV,NTP,NTS,NTV,VERTP,XNC,
     -     XNS,YNC,YNS,ZNC,ZNS)
      XNC=-0.05D0
      YNC=1.0D0
      ZNC=-0.1D0
      DMOD=(XNC**2+YNC**2+ZNC**2)**0.5D0
      XNC=XNC/DMOD
      YNC=YNC/DMOD
      ZNC=ZNC/DMOD
      CALL INTE3D(C,ICONTN,ICONTP,IPV,NIPV,NTP,NTS,NTV,VERTP,XNC,
     -     XNS,YNC,YNS,ZNC,ZNS)
      XNC=-0.05D0
      YNC=-0.1D0
      ZNC=1.0D0
      DMOD=(XNC**2+YNC**2+ZNC**2)**0.5D0
      XNC=XNC/DMOD
      YNC=YNC/DMOD
      ZNC=ZNC/DMOD
      CALL INTE3D(C,ICONTN,ICONTP,IPV,NIPV,NTP,NTS,NTV,VERTP,XNC,
     -     XNS,YNC,YNS,ZNC,ZNS)
c.. Cutting planes passing through the point (1,1,1)
      XNC=-1.0D0
      YNC=0.05D0
      ZNC=0.1D0
      DMOD=(XNC**2+YNC**2+ZNC**2)**0.5D0
      XNC=XNC/DMOD
      YNC=YNC/DMOD
      ZNC=ZNC/DMOD
      C=-1.0D0*(XNC+YNC+ZNC)
      CALL INTE3D(C,ICONTN,ICONTP,IPV,NIPV,NTP,NTS,NTV,VERTP,XNC,
     -     XNS,YNC,YNS,ZNC,ZNS)
      XNC=0.05D0
      YNC=-1.0D0
      ZNC=0.025D0
      DMOD=(XNC**2+YNC**2+ZNC**2)**0.5D0
      XNC=XNC/DMOD
      YNC=YNC/DMOD
      ZNC=ZNC/DMOD
      C=-1.0D0*(XNC+YNC+ZNC)
      CALL INTE3D(C,ICONTN,ICONTP,IPV,NIPV,NTP,NTS,NTV,VERTP,XNC,
     -     XNS,YNC,YNS,ZNC,ZNS)
      XNC=0.05D0
      YNC=0.05D0
      ZNC=-1.0D0
      DMOD=(XNC**2+YNC**2+ZNC**2)**0.5D0
      XNC=XNC/DMOD
      YNC=YNC/DMOD
      ZNC=ZNC/DMOD
      C=-1.0D0*(XNC+YNC+ZNC)
      CALL INTE3D(C,ICONTN,ICONTP,IPV,NIPV,NTP,NTS,NTV,VERTP,XNC,
     -     XNS,YNC,YNS,ZNC,ZNS)
      CALL RESTORE3D(CS,IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      RETURN
      END
c------------------------ END OF HEXAHEMESH --------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            TETRAMESH                                c
c---------------------------------------------------------------------c
c      This version dated: August 9, 2014 by Lopez and Hernandez      c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE TETRAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
      NTS=4
      NTV=4
      NTP=NTV
      NIPV(1)=3
      IPV(1,1)=1
      IPV(1,2)=3
      IPV(1,3)=2
      NIPV(2)=3
      IPV(2,1)=3
      IPV(2,2)=1
      IPV(2,3)=4
      NIPV(3)=3
      IPV(3,1)=2
      IPV(3,2)=3
      IPV(3,3)=4
      NIPV(4)=3
      IPV(4,1)=1
      IPV(4,2)=2
      IPV(4,3)=4
      VERTP(1,1)=0.0D0
      VERTP(1,2)=0.0D0
      VERTP(1,3)=0.0D0
      VERTP(2,1)=0.91D0
      VERTP(2,2)=0.24D0
      VERTP(2,3)=1.0D0
      VERTP(3,1)=0.72D0
      VERTP(3,2)=0.16D0
      VERTP(3,3)=0.07D0
      VERTP(4,1)=1.0D0
      VERTP(4,2)=1.0D0
      VERTP(4,3)=1.0D0
      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5D0
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO

      RETURN
      END
c------------------------- END OF TETRAMESH --------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            DODECAMESH                               c
c---------------------------------------------------------------------c
c      This version dated: April 10, 2007 by Lopez and Hernandez      c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE DODECAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
      A=1.0D0/(3.0D0)**0.5D0
      B=((3.0D0-(5.0D0)**0.5D0)/6.0D0)**0.5D0
      C=((3.0D0+(5.0D0)**0.5D0)/6.0D0)**0.5D0
      NTS=12
      NTV=20
      NTP=NTV
      I=1
      VERTP(I,1)=A
      VERTP(I,2)=A
      VERTP(I,3)=A
      I=I+1
      VERTP(I,1)=A
      VERTP(I,2)=A
      VERTP(I,3)=-A
      I=I+1
      VERTP(I,1)=A
      VERTP(I,2)=-A
      VERTP(I,3)=A
      I=I+1
      VERTP(I,1)=A
      VERTP(I,2)=-A
      VERTP(I,3)=-A
      I=I+1
      VERTP(I,1)=-A
      VERTP(I,2)=A
      VERTP(I,3)=A
      I=I+1
      VERTP(I,1)=-A
      VERTP(I,2)=A
      VERTP(I,3)=-A
      I=I+1
      VERTP(I,1)=-A
      VERTP(I,2)=-A
      VERTP(I,3)=A
      I=I+1
      VERTP(I,1)=-A
      VERTP(I,2)=-A
      VERTP(I,3)=-A
      I=I+1
      VERTP(I,1)=B
      VERTP(I,2)=C
      VERTP(I,3)=0.0D0
      I=I+1
      VERTP(I,1)=-B
      VERTP(I,2)=C
      VERTP(I,3)=0.0D0
      I=I+1
      VERTP(I,1)=B
      VERTP(I,2)=-C
      VERTP(I,3)=0.0D0
      I=I+1
      VERTP(I,1)=-B
      VERTP(I,2)=-C
      VERTP(I,3)=0.0D0
      I=I+1
      VERTP(I,1)=C
      VERTP(I,2)=0.0D0
      VERTP(I,3)=B
      I=I+1
      VERTP(I,1)=C
      VERTP(I,2)=0.0D0
      VERTP(I,3)=-B
      I=I+1
      VERTP(I,1)=-C
      VERTP(I,2)=0.0D0
      VERTP(I,3)=B
      I=I+1
      VERTP(I,1)=-C
      VERTP(I,2)=0.0D0
      VERTP(I,3)=-B
      I=I+1
      VERTP(I,1)=0.0D0
      VERTP(I,2)=B
      VERTP(I,3)=C
      I=I+1
      VERTP(I,1)=0.0D0
      VERTP(I,2)=-B
      VERTP(I,3)=C
      I=I+1
      VERTP(I,1)=0.0D0
      VERTP(I,2)=B
      VERTP(I,3)=-C
      I=I+1
      VERTP(I,1)=0.0D0
      VERTP(I,2)=-B
      VERTP(I,3)=-C
      DO IS=1,NTS
         NIPV(IS)=5
      END DO
      IS=1
      I=1
      IPV(IS,I)=1
      I=I+1
      IPV(IS,I)=9
      I=I+1
      IPV(IS,I)=10
      I=I+1
      IPV(IS,I)=5
      I=I+1
      IPV(IS,I)=17
      IS=IS+1
      I=1
      IPV(IS,I)=1
      I=I+1
      IPV(IS,I)=17
      I=I+1
      IPV(IS,I)=18
      I=I+1
      IPV(IS,I)=3
      I=I+1
      IPV(IS,I)=13
      IS=IS+1
      I=1
      IPV(IS,I)=13
      I=I+1
      IPV(IS,I)=3
      I=I+1
      IPV(IS,I)=11
      I=I+1
      IPV(IS,I)=4
      I=I+1
      IPV(IS,I)=14
      IS=IS+1
      I=1
      IPV(IS,I)=10
      I=I+1
      IPV(IS,I)=6
      I=I+1
      IPV(IS,I)=16
      I=I+1
      IPV(IS,I)=15
      I=I+1
      IPV(IS,I)=5
      IS=IS+1
      I=1
      IPV(IS,I)=4
      I=I+1
      IPV(IS,I)=20
      I=I+1
      IPV(IS,I)=19
      I=I+1
      IPV(IS,I)=2
      I=I+1
      IPV(IS,I)=14
      IS=IS+1
      I=1
      IPV(IS,I)=8
      I=I+1
      IPV(IS,I)=12
      I=I+1
      IPV(IS,I)=7
      I=I+1
      IPV(IS,I)=15
      I=I+1
      IPV(IS,I)=16
      IS=IS+1
      I=1
      IPV(IS,I)=1
      I=I+1
      IPV(IS,I)=13
      I=I+1
      IPV(IS,I)=14
      I=I+1
      IPV(IS,I)=2
      I=I+1
      IPV(IS,I)=9
      IS=IS+1
      I=1
      IPV(IS,I)=9
      I=I+1
      IPV(IS,I)=2
      I=I+1
      IPV(IS,I)=19
      I=I+1
      IPV(IS,I)=6
      I=I+1
      IPV(IS,I)=10
      IS=IS+1
      I=1
      IPV(IS,I)=17
      I=I+1
      IPV(IS,I)=5
      I=I+1
      IPV(IS,I)=15
      I=I+1
      IPV(IS,I)=7
      I=I+1
      IPV(IS,I)=18
      IS=IS+1
      I=1
      IPV(IS,I)=7
      I=I+1
      IPV(IS,I)=12
      I=I+1
      IPV(IS,I)=11
      I=I+1
      IPV(IS,I)=3
      I=I+1
      IPV(IS,I)=18
      IS=IS+1
      I=1
      IPV(IS,I)=8
      I=I+1
      IPV(IS,I)=16
      I=I+1
      IPV(IS,I)=6
      I=I+1
      IPV(IS,I)=19
      I=I+1
      IPV(IS,I)=20
      IS=IS+1
      I=1
      IPV(IS,I)=8
      I=I+1
      IPV(IS,I)=20
      I=I+1
      IPV(IS,I)=4
      I=I+1
      IPV(IS,I)=11
      I=I+1
      IPV(IS,I)=12
      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5D0
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO
      RETURN
      END
c------------------------- END OF DODECAMESH -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            ICOSAMESH                                c
c---------------------------------------------------------------------c
c      This version dated: April 10, 2007 by Lopez and Hernandez      c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE ICOSAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
      T=(1.0D0+(5.0D0)**0.5D0)/2.0D0
      A=(1.0D0+T**2.0D0)**0.5D0
      NTS=20
      NTV=12
      NTP=NTV
      I=1
      VERTP(I,1)=T/A
      VERTP(I,2)=1.0D0/A
      VERTP(I,3)=0.0D0/A
      I=I+1
      VERTP(I,1)=-T/A
      VERTP(I,2)=1.0D0/A
      VERTP(I,3)=0.0D0/A
      I=I+1
      VERTP(I,1)=T/A
      VERTP(I,2)=-1.0D0/A
      VERTP(I,3)=0.0D0/A
      I=I+1
      VERTP(I,1)=-T/A
      VERTP(I,2)=-1.0D0/A
      VERTP(I,3)=0.0D0/A
      I=I+1
      VERTP(I,1)=1.0D0/A
      VERTP(I,2)=0.0D0/A
      VERTP(I,3)=T/A
      I=I+1
      VERTP(I,1)=1.0D0/A
      VERTP(I,2)=0.0D0/A
      VERTP(I,3)=-T/A
      I=I+1
      VERTP(I,1)=-1.0D0/A
      VERTP(I,2)=0.0D0/A
      VERTP(I,3)=T/A
      I=I+1
      VERTP(I,1)=-1.0D0/A
      VERTP(I,2)=0.0D0/A
      VERTP(I,3)=-T/A
      I=I+1
      VERTP(I,1)=0.0D0/A
      VERTP(I,2)=T/A
      VERTP(I,3)=1.0D0/A
      I=I+1
      VERTP(I,1)=0.0D0/A
      VERTP(I,2)=-T/A
      VERTP(I,3)=1.0D0/A
      I=I+1
      VERTP(I,1)=0.0D0/A
      VERTP(I,2)=T/A
      VERTP(I,3)=-1.0D0/A
      I=I+1
      VERTP(I,1)=0.0D0/A
      VERTP(I,2)=-T/A
      VERTP(I,3)=-1.0D0/A
      DO IS=1,NTS
         NIPV(IS)=3
      END DO
      IS=1
      I=1
      IPV(IS,I)=1
      I=I+1
      IPV(IS,I)=9
      I=I+1
      IPV(IS,I)=5
      IS=IS+1
      I=1
      IPV(IS,I)=1
      I=I+1
      IPV(IS,I)=6
      I=I+1
      IPV(IS,I)=11
      IS=IS+1
      I=1
      IPV(IS,I)=3
      I=I+1
      IPV(IS,I)=5
      I=I+1
      IPV(IS,I)=10
      IS=IS+1
      I=1
      IPV(IS,I)=3
      I=I+1
      IPV(IS,I)=12
      I=I+1
      IPV(IS,I)=6
      IS=IS+1
      I=1
      IPV(IS,I)=2
      I=I+1
      IPV(IS,I)=7
      I=I+1
      IPV(IS,I)=9
      IS=IS+1
      I=1
      IPV(IS,I)=2
      I=I+1
      IPV(IS,I)=11
      I=I+1
      IPV(IS,I)=8
      IS=IS+1
      I=1
      IPV(IS,I)=4
      I=I+1
      IPV(IS,I)=10
      I=I+1
      IPV(IS,I)=7
      IS=IS+1
      I=1
      IPV(IS,I)=4
      I=I+1
      IPV(IS,I)=8
      I=I+1
      IPV(IS,I)=12
      IS=IS+1
      I=1
      IPV(IS,I)=1
      I=I+1
      IPV(IS,I)=11
      I=I+1
      IPV(IS,I)=9
      IS=IS+1
      I=1
      IPV(IS,I)=2
      I=I+1
      IPV(IS,I)=9
      I=I+1
      IPV(IS,I)=11
      IS=IS+1
      I=1
      IPV(IS,I)=3
      I=I+1
      IPV(IS,I)=10
      I=I+1
      IPV(IS,I)=12
      IS=IS+1
      I=1
      IPV(IS,I)=4
      I=I+1
      IPV(IS,I)=12
      I=I+1
      IPV(IS,I)=10
      IS=IS+1
      I=1
      IPV(IS,I)=5
      I=I+1
      IPV(IS,I)=3
      I=I+1
      IPV(IS,I)=1
      IS=IS+1
      I=1
      IPV(IS,I)=6
      I=I+1
      IPV(IS,I)=1
      I=I+1
      IPV(IS,I)=3
      IS=IS+1
      I=1
      IPV(IS,I)=7
      I=I+1
      IPV(IS,I)=2
      I=I+1
      IPV(IS,I)=4
      IS=IS+1
      I=1
      IPV(IS,I)=8
      I=I+1
      IPV(IS,I)=4
      I=I+1
      IPV(IS,I)=2
      IS=IS+1
      I=1
      IPV(IS,I)=9
      I=I+1
      IPV(IS,I)=7
      I=I+1
      IPV(IS,I)=5
      IS=IS+1
      I=1
      IPV(IS,I)=10
      I=I+1
      IPV(IS,I)=5
      I=I+1
      IPV(IS,I)=7
      IS=IS+1
      I=1
      IPV(IS,I)=11
      I=I+1
      IPV(IS,I)=6
      I=I+1
      IPV(IS,I)=8
      IS=IS+1
      I=1
      IPV(IS,I)=12
      I=I+1
      IPV(IS,I)=8
      I=I+1
      IPV(IS,I)=6
      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO
      RETURN
      END
c------------------------- END OF ICOSAMESH --------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                           COMPLEXMESH                               c
c---------------------------------------------------------------------c
c      This version dated: FEB 25, 2015 by Lopez and Hernandez        c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE COMPLEXMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)

      NTP=32
      NTV=NTP
      NTS=18
      VERTP(1,1)=0.542827704131611521D0
      VERTP(1,2)=0.103161810115955432D0
      VERTP(1,3)=0.036812175979487702D0
      VERTP(2,1)=0.600855902479279003D0
      VERTP(2,2)=0.086471258131132003D0
      VERTP(2,3)=0.056436950730949877D0
      VERTP(3,1)=0.580277964506970667D0
      VERTP(3,2)=0.163039123851804080D0
      VERTP(3,3)=0.063539333285191152D0
      VERTP(4,1)=0.595353203042033874D0
      VERTP(4,2)=0.109843129311963592D0
      VERTP(4,3)=0.046191353733568030D0
      VERTP(5,1)=0.543760908872757964D0
      VERTP(5,2)=0.094099097938159251D0
      VERTP(5,3)=0.040511888762036471D0
      VERTP(6,1)=0.576819361018633514D0
      VERTP(6,2)=0.166401017225910497D0
      VERTP(6,3)=0.070695345095658793D0
      VERTP(7,1)=0.563208132688422847D0
      VERTP(7,2)=0.155406657743345611D0
      VERTP(7,3)=0.108888651199260125D0
      VERTP(8,1)=0.539375361775393247D0
      VERTP(8,2)=0.123540417024526422D0
      VERTP(8,3)=0.108841211139604044D0
      VERTP(9,1)=0.567496472062291590D0
      VERTP(9,2)=0.128649916760601390D0
      VERTP(9,3)=0.113035361390195654D0
      VERTP(10,1)=0.600028673063590978D0
      VERTP(10,2)=0.085143233872209040D0
      VERTP(10,3)=0.059111439263878407D0
      VERTP(11,1)=0.526213727798502173D0
      VERTP(11,2)=0.075636612108421930D0
      VERTP(11,3)=0.077638126137438007D0
      VERTP(12,1)=0.585556321506447985D0
      VERTP(12,2)=0.087406427365641859D0
      VERTP(12,3)=0.092250691458415329D0
      VERTP(13,1)=0.559954651117444135D0
      VERTP(13,2)=0.133408124832590513D0
      VERTP(13,3)=0.115142136143722096D0
      VERTP(14,1)=0.524564746594361031D0
      VERTP(14,2)=0.076122995267018601D0
      VERTP(14,3)=0.083047795609934638D0
      VERTP(15,1)=0.564747090248686634D0
      VERTP(15,2)=0.145410923059303421D0
      VERTP(15,3)=0.110587477743786994D0
      VERTP(16,1)=0.514395820566891815D0
      VERTP(16,2)=0.112406169868730282D0
      VERTP(16,3)=0.074847147170632553D0
      VERTP(17,1)=0.561075920034260656D0
      VERTP(17,2)=0.154023426807431363D0
      VERTP(17,3)=0.110265478973935432D0
      VERTP(18,1)=0.600086961628965465D0
      VERTP(18,2)=0.099753327716150253D0
      VERTP(18,3)=0.051150970749093465D0
      VERTP(19,1)=0.523529840650993616D0
      VERTP(19,2)=0.075758256340824406D0
      VERTP(19,3)=0.081652754048494577D0
      VERTP(20,1)=0.580306952541070453D0
      VERTP(20,2)=0.163300101039615730D0
      VERTP(20,3)=0.064077552617180802D0
      VERTP(21,1)=0.560784564570740995D0
      VERTP(21,2)=0.154556861208633767D0
      VERTP(21,3)=0.110001263915814051D0
      VERTP(22,1)=0.562766593882476296D0
      VERTP(22,2)=0.155209023571326599D0
      VERTP(22,3)=0.109311173701962666D0
      VERTP(23,1)=0.535792240519188834D0
      VERTP(23,2)=0.158071087967216972D0
      VERTP(23,3)=0.055846886029540778D0
      VERTP(24,1)=0.562461680790533380D0
      VERTP(24,2)=0.155100229495437003D0
      VERTP(24,3)=0.109460804789442229D0
      VERTP(25,1)=0.528789432867762255D0
      VERTP(25,2)=0.128164972986322317D0
      VERTP(25,3)=0.044046149723425451D0
      VERTP(26,1)=0.535717736632420283D0
      VERTP(26,2)=0.158416280041318303D0
      VERTP(26,3)=0.056565341810804984D0
      VERTP(27,1)=0.527966627100669772D0
      VERTP(27,2)=0.148427477131785279D0
      VERTP(27,3)=0.097455495587370988D0
      VERTP(28,1)=0.531197945823301376D0
      VERTP(28,2)=0.154922874892573614D0
      VERTP(28,3)=0.054094280101616717D0
      VERTP(29,1)=0.519703120392784435D0
      VERTP(29,2)=0.149093981287893751D0
      VERTP(29,3)=0.082950196610910618D0
      VERTP(30,1)=0.530093766604715744D0
      VERTP(30,2)=0.156612444958861563D0
      VERTP(30,3)=0.058283874375187901D0
      VERTP(31,1)=0.518229286434464531D0
      VERTP(31,2)=0.144572125378080146D0
      VERTP(31,3)=0.084271163208895244D0
      VERTP(32,1)=0.520436917533184662D0
      VERTP(32,2)=0.148367611160416940D0
      VERTP(32,3)=0.087663831740394063D0

      NIPV(1)= 5
      NIPV(2)=9
      NIPV(3)=4
      NIPV(4)=5
      NIPV(5)=4
      NIPV(6)=6
      NIPV(7)=6
      NIPV(8)=10
      NIPV(9)=3
      NIPV(10)=7
      NIPV(11)=5
      NIPV(12)=5
      NIPV(13)=6
      NIPV(14)=5
      NIPV(15)=6
      NIPV(16)=3
      NIPV(17)=4
      NIPV(18)=3

c      IPV(1,1:NIPV(1))=(/5,1,4,18,2/)
      IPV(1,1)=5
      IPV(1,2)=1
      IPV(1,3)=4
      IPV(1,4)=18
      IPV(1,5)=2

c      IPV(2,1:NIPV(2))=(/18,20,6,7,15,9,12,10,2/)
      IPV(2,1)=18
      IPV(2,2)=20
      IPV(2,3)=6
      IPV(2,4)=7
      IPV(2,5)=15
      IPV(2,6)=9
      IPV(2,7)=12
      IPV(2,8)=10
      IPV(2,9)=2

c      IPV(3,::)=(/10,11,5,2/)
      IPV(3,1)=10
      IPV(3,2)=11
      IPV(3,3)=5
      IPV(3,4)=2

c      IPV(4,1:NIPV(4))=(/23,26,6,20,3/)
      IPV(4,1)=23
      IPV(4,2)=26
      IPV(4,3)=6
      IPV(4,4)=20
      IPV(4,5)=3

c      IPV(5,1:NIPV(5))=(/20,18,4,3/)
      IPV(5,1)=20
      IPV(5,2)=18
      IPV(5,3)=4
      IPV(5,4)=3

c      IPV(6,1:NIPV(6))=(/4,1,25,28,23,3/)
      IPV(6,1)=4
      IPV(6,2)=1
      IPV(6,3)=25
      IPV(6,4)=28
      IPV(6,5)=23
      IPV(6,6)=3

c      IPV(7,1:NIPV(7))=(/11,19,16,25,1,5/)
      IPV(7,1)=11
      IPV(7,2)=19
      IPV(7,3)=16
      IPV(7,4)=25
      IPV(7,5)=1
      IPV(7,6)=5

c      IPV(8,1:NIPV(8))=(/26,30,29,32,27,21,24,22,7,6/)
      IPV(8,1)=26
      IPV(8,2)=30
      IPV(8,3)=29
      IPV(8,4)=32
      IPV(8,5)=27
      IPV(8,6)=21
      IPV(8,7)=24
      IPV(8,8)=22
      IPV(8,9)=7
      IPV(8,10)=6

c      IPV(9,1:NIPV(9))=(/22,15,7/)
      IPV(9,1)=22
      IPV(9,2)=15
      IPV(9,3)=7

c      IPV(10,1:NIPV(10))=(/27,32,31,16,19,14,8/) 
      IPV(10,1)=27 
      IPV(10,2)=32 
      IPV(10,3)=31 
      IPV(10,4)=16 
      IPV(10,5)=19 
      IPV(10,6)=14 
      IPV(10,7)=8 

c      IPV(11,1:NIPV(11))=(/14,12,9,13,8/)
      IPV(11,1)=14
      IPV(11,2)=12
      IPV(11,3)=9
      IPV(11,4)=13
      IPV(11,5)=8

c      IPV(12,1:NIPV(12))=(/13,17,21,27,8/)
      IPV(12,1)=13
      IPV(12,2)=17
      IPV(12,3)=21
      IPV(12,4)=27
      IPV(12,5)=8

c      IPV(13,1:NIPV(13))=(/15,22,24,17,13,9/)
      IPV(13,1)=15
      IPV(13,2)=22
      IPV(13,3)=24
      IPV(13,4)=17
      IPV(13,5)=13
      IPV(13,6)=9

c      IPV(14,1:NIPV(14))=(/12,14,19,11,10/)
      IPV(14,1)=12
      IPV(14,2)=14
      IPV(14,3)=19
      IPV(14,4)=11
      IPV(14,5)=10

c      IPV(15,1:NIPV(15))=(/31,29,30,28,25,16/)
      IPV(15,1)=31
      IPV(15,2)=29
      IPV(15,3)=30
      IPV(15,4)=28
      IPV(15,5)=25
      IPV(15,6)=16

c      IPV(16,1:NIPV(16))=(/24,21,17/)
      IPV(16,1)=24
      IPV(16,2)=21
      IPV(16,3)=17

c      IPV(17,1:NIPV(17))=(/28,30,26,23/)
      IPV(17,1)=28
      IPV(17,2)=30
      IPV(17,3)=26
      IPV(17,4)=23

c      IPV(18,1:NIPV(18))=(/31,32,29/)
      IPV(18,1)=31
      IPV(18,2)=32
      IPV(18,3)=29

      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO

      RETURN
      END
c------------------------ END OF COMPLEXMESH -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                         DISTORTCUBICMESH                            c
c---------------------------------------------------------------------c
c    CELL OBTAINED FROM THE ORIGNAL CUBICMESH BY MOVING VERTICES      c 
c    4 AND 7 ALONG THE Y AXIS AND DECOMPOSING THE ORIGINAL FACE 4     c
c               INTO THE NEW TRIANGULAR FACES 4 AND 7                 c     
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c      This version dated: March 2, 2017 by Lopez and Hernandez       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE DISTORTCUBICMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,
     -     ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
      D0=0.0D0
      D1=1.0D0

      XNS(1)=D1
      YNS(1)=D0
      ZNS(1)=D0
      XNS(2)=D0
      YNS(2)=-D1
      ZNS(2)=D0
      XNS(3)=D0
      YNS(3)=D0
      ZNS(3)=-D1
      XNS(5)=D0
      YNS(5)=D0
      ZNS(5)=D1
      XNS(6)=-D1
      YNS(6)=D0
      ZNS(6)=D0
      NTS=7
      NTV=8
      NTP=NTV
      NIPV(1)=4
      IPV(1,1)=1
      IPV(1,2)=2
      IPV(1,3)=3
      IPV(1,4)=4
      NIPV(2)=4
      IPV(2,1)=2
      IPV(2,2)=1
      IPV(2,3)=5
      IPV(2,4)=6
      NIPV(3)=4
      IPV(3,1)=3
      IPV(3,2)=2
      IPV(3,3)=6
      IPV(3,4)=7
      NIPV(4)=3
      IPV(4,1)=4
      IPV(4,2)=3
      IPV(4,3)=7
      NIPV(5)=4
      IPV(5,1)=1
      IPV(5,2)=4
      IPV(5,3)=8
      IPV(5,4)=5
      NIPV(6)=4
      IPV(6,1)=6
      IPV(6,2)=5
      IPV(6,3)=8
      IPV(6,4)=7
      NIPV(7)=3
      IPV(7,1)=4
      IPV(7,2)=7
      IPV(7,3)=8
      VERTP(1,1)=D1
      VERTP(1,2)=D0
      VERTP(1,3)=D1
      VERTP(2,1)=D1
      VERTP(2,2)=D0
      VERTP(2,3)=D0
      VERTP(3,1)=D1
      VERTP(3,2)=D1
      VERTP(3,3)=D0
      VERTP(4,1)=D1
      VERTP(4,2)=0.25D0
      VERTP(4,3)=D1
      VERTP(5,1)=D0
      VERTP(5,2)=D0
      VERTP(5,3)=D1
      VERTP(6,1)=D0
      VERTP(6,2)=D0
      VERTP(6,3)=D0
      VERTP(7,1)=D0
      VERTP(7,2)=0.25D0
      VERTP(7,3)=D0
      VERTP(8,1)=D0
      VERTP(8,2)=D1
      VERTP(8,3)=D1
c Orientation of the triangular face 4:
      IS=4
      XV1=VERTP(IPV(IS,2),1)-VERTP(IPv(IS,1),1)
      YV1=VERTP(IPv(IS,2),2)-VERTP(IPv(IS,1),2)
      ZV1=VERTP(IPv(IS,2),3)-VERTP(IPv(IS,1),3)
      XV2=VERTP(IPv(IS,3),1)-VERTP(IPv(IS,2),1)
      YV2=VERTP(IPv(IS,3),2)-VERTP(IPv(IS,2),2)
      ZV2=VERTP(IPv(IS,3),3)-VERTP(IPv(IS,2),3)
      XM=YV1*ZV2-ZV1*YV2
      YM=ZV1*XV2-XV1*ZV2
      ZM=XV1*YV2-YV1*XV2
      AMOD=(XM**2.0+YM**2.0+ZM**2.0)**0.5
      XNS(IS)=XM/AMOD
      YNS(IS)=YM/AMOD
      ZNS(IS)=ZM/AMOD
c Orientation of the triangular face 7:
      IS=7
      XV1=VERTP(IPV(IS,2),1)-VERTP(IPv(IS,1),1)
      YV1=VERTP(IPv(IS,2),2)-VERTP(IPv(IS,1),2)
      ZV1=VERTP(IPv(IS,2),3)-VERTP(IPv(IS,1),3)
      XV2=VERTP(IPv(IS,3),1)-VERTP(IPv(IS,2),1)
      YV2=VERTP(IPv(IS,3),2)-VERTP(IPv(IS,2),2)
      ZV2=VERTP(IPv(IS,3),3)-VERTP(IPv(IS,2),3)
      XM=YV1*ZV2-ZV1*YV2
      YM=ZV1*XV2-XV1*ZV2
      ZM=XV1*YV2-YV1*XV2
      AMOD=(XM**2.0+YM**2.0+ZM**2.0)**0.5
      XNS(IS)=XM/AMOD
      YNS(IS)=YM/AMOD
      ZNS(IS)=ZM/AMOD

      
      RETURN
      END
c--------------------- END OF DISTORTCUBICMESH -----------------------c
c---------------------------------------------------------------------c
      
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                   NON-CONVEX PENTAGONAL PYRAMID                     c
c---------------------------------------------------------------------c
c      This version dated: FEB 25, 2015 by Lopez and Hernandez        c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCPENTAPYRAMID(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),YNS(NS),ZNS(NS)

      NTS=6
      NTV=6
      NTP=NTV

      VERTP(1,1)=0.04
      VERTP(1,2)=0.77
      VERTP(1,3)=0.0
      VERTP(2,1)=0.0
      VERTP(2,2)=0.0
      VERTP(2,3)=0.0
      VERTP(3,1)=0.49
      VERTP(3,2)=0.22
      VERTP(3,3)=0.0
      VERTP(4,1)=1.0
      VERTP(4,2)=0.13
      VERTP(4,3)=0.0
      VERTP(5,1)=0.16
      VERTP(5,2)=1.0
      VERTP(5,3)=0.0
      VERTP(6,1)=0.1
      VERTP(6,2)=0.5
      VERTP(6,3)=1.0

c. IS=1
      NIPV(1)=5
      IPV(1,1)=1
      IPV(1,2)=5
      IPV(1,3)=4
      IPV(1,4)=3
      IPV(1,5)=2
c. IS=2
      NIPV(2)=3
      IPV(2,1)=1
      IPV(2,2)=2
      IPV(2,3)=6
c. IS=3
      NIPV(3)=3
      IPV(3,1)=2
      IPV(3,2)=3
      IPV(3,3)=6
c. IS=4
      NIPV(4)=3
      IPV(4,1)=3
      IPV(4,2)=4
      IPV(4,3)=6
c. IS=5
      NIPV(5)=3
      IPV(5,1)=4
      IPV(5,2)=5
      IPV(5,3)=6
c. IS=6
      NIPV(6)=3
      IPV(6,1)=5
      IPV(6,2)=1
      IPV(6,3)=6

      DO IS=1,NTS
            XV1=VERTP(IPv(IS,2),1)-VERTP(IPv(IS,1),1)
            YV1=VERTP(IPv(IS,2),2)-VERTP(IPv(IS,1),2)
            ZV1=VERTP(IPv(IS,2),3)-VERTP(IPv(IS,1),3)
            XV2=VERTP(IPv(IS,3),1)-VERTP(IPv(IS,2),1)
            YV2=VERTP(IPv(IS,3),2)-VERTP(IPv(IS,2),2)
            ZV2=VERTP(IPv(IS,3),3)-VERTP(IPv(IS,2),3)
            XM=YV1*ZV2-ZV1*YV2
            YM=ZV1*XV2-XV1*ZV2
            ZM=XV1*YV2-YV1*XV2
            AMOD=(XM**2.0+YM**2.0+ZM**2.0)**0.5
            XNS(IS)=XM/AMOD
            YNS(IS)=YM/AMOD
            ZNS(IS)=ZM/AMOD
      END DO

      RETURN
      END
c----------------------- END OF NCPENTAPYRAMID -----------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c   NON-CONVEX CELL: CELL OBTAINED BY SUBTRACTING A PYRAMID TO        c 
c                         A UNIT CUBIC CELL                           c
c---------------------------------------------------------------------c
c      This version dated: FEB 25, 2015 by Lopez and Hernandez        c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCCUBICPYRAMID(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
      D0=0.0D0
      D1=1.0D0

      NTS=9
      NTV=9
      NTP=NTV
      NIPV(1)=4
      IPV(1,1)=1
      IPV(1,2)=2
      IPV(1,3)=3
      IPV(1,4)=4
      NIPV(2)=4
      IPV(2,1)=2
      IPV(2,2)=1
      IPV(2,3)=5
      IPV(2,4)=6
      NIPV(3)=4
      IPV(3,1)=3
      IPV(3,2)=2
      IPV(3,3)=6
      IPV(3,4)=7
      NIPV(4)=3
      IPV(4,1)=4
      IPV(4,2)=9
      IPV(4,3)=8
      NIPV(5)=4
      IPV(5,1)=1
      IPV(5,2)=4
      IPV(5,3)=8
      IPV(5,4)=5
      NIPV(6)=4
      IPV(6,1)=6
      IPV(6,2)=5
      IPV(6,3)=8
      IPV(6,4)=7

      NIPV(7)=3
      IPV(7,1)=7
      IPV(7,2)=8
      IPV(7,3)=9

      NIPV(8)=3
      IPV(8,1)=3
      IPV(8,2)=7
      IPV(8,3)=9

      NIPV(9)=3
      IPV(9,1)=4
      IPV(9,2)=3
      IPV(9,3)=9

      VERTP(1,1)=D1
      VERTP(1,2)=D0
      VERTP(1,3)=D1
      VERTP(2,1)=D1
      VERTP(2,2)=D0
      VERTP(2,3)=D0
      VERTP(3,1)=D1
      VERTP(3,2)=D1
      VERTP(3,3)=D0
      VERTP(4,1)=D1
      VERTP(4,2)=D1
      VERTP(4,3)=D1
      VERTP(5,1)=D0
      VERTP(5,2)=D0
      VERTP(5,3)=D1
      VERTP(6,1)=D0
      VERTP(6,2)=D0
      VERTP(6,3)=D0
      VERTP(7,1)=D0
      VERTP(7,2)=D1
      VERTP(7,3)=D0
      VERTP(8,1)=D0
      VERTP(8,2)=D1
      VERTP(8,3)=D1

      VERTP(9,1)=0.5D0
      VERTP(9,2)=0.5D0
      VERTP(9,3)=0.5D0

      DO IS=1,NTS
            XV1=VERTP(IPv(IS,2),1)-VERTP(IPv(IS,1),1)
            YV1=VERTP(IPv(IS,2),2)-VERTP(IPv(IS,1),2)
            ZV1=VERTP(IPv(IS,2),3)-VERTP(IPv(IS,1),3)
            XV2=VERTP(IPv(IS,3),1)-VERTP(IPv(IS,2),1)
            YV2=VERTP(IPv(IS,3),2)-VERTP(IPv(IS,2),2)
            ZV2=VERTP(IPv(IS,3),3)-VERTP(IPv(IS,2),3)
            XM=YV1*ZV2-ZV1*YV2
            YM=ZV1*XV2-XV1*ZV2
            ZM=XV1*YV2-YV1*XV2
            AMOD=(XM**2.0+YM**2.0+ZM**2.0)**0.5
            XNS(IS)=XM/AMOD
            YNS(IS)=YM/AMOD
            ZNS(IS)=ZM/AMOD
      END DO

      RETURN
      END
c----------------------- END OF NCCUBICPYRAMID -----------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                      SMALL STELLATED CUBICMESH                      c
c---------------------------------------------------------------------c
c      This version dated: November 5, 2015 by Lopez and Hernandez    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCSCUBICMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)

      CALL CUBICMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      NTS0=NTS
      NTSI=NTS
      NTP0=NTP
      A=1D0
c.. Face centroid      
      DO IS=1,NTS
         XC=0.0
         YC=0.0
         ZC=0.0
         DO IV=1,NIPV(IS)
            IP=IPV(IS,IV)
            XC=XC+VERTP(IP,1)
            YC=YC+VERTP(IP,2)
            ZC=ZC+VERTP(IP,3)
         END DO
         NTP0=NTP0+1
         VERTP(NTP0,1)=XC/NIPV(IS)+A*XNS(IS)
         VERTP(NTP0,2)=YC/NIPV(IS)+A*YNS(IS)
         VERTP(NTP0,3)=ZC/NIPV(IS)+A*ZNS(IS)
         
         DO IV=1,NIPV(IS)
            IV2=IV+1
            IF(IV.EQ.NIPV(IS)) IV2=1
            NTS0=NTS0+1
            NIPV(NTS0)=3
            IPV(NTS0,1)=IPV(IS,IV)
            IPV(NTS0,2)=IPV(IS,IV2)
            IPV(NTS0,3)=NTP0
         END DO
      END DO

      NTS=NTS0-NTSI
      DO IS=1,NTS
         NIPV(IS)=NIPV(NTSI+IS)
         DO IV=1,NIPV(IS)
            IPV(IS,IV)=IPV(NTSI+IS,IV)
         END DO
      END DO
      NTP=NTP0
      NTV=NTP
      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5D0
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO
      RETURN
      END
c------------------------ END OF NCSCUBICMESH ------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                     NON-CONVEX HEXAHEDRON MESH                      c
c---------------------------------------------------------------------c
c      This version dated: Novemeber 4, 2015 by Lopez and Hernandez   c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCHEXAHEMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
      D0=0.0D0
      D1=1.0D0

      NTS=6
      NTV=8
      NTP=NTV
      NIPV(1)=4
      IPV(1,1)=1
      IPV(1,2)=2
      IPV(1,3)=3
      IPV(1,4)=4
      NIPV(2)=4
      IPV(2,1)=2
      IPV(2,2)=1
      IPV(2,3)=5
      IPV(2,4)=6
      NIPV(3)=4
      IPV(3,1)=3
      IPV(3,2)=2
      IPV(3,3)=6
      IPV(3,4)=7
      NIPV(4)=4
      IPV(4,1)=4
      IPV(4,2)=3
      IPV(4,3)=7
      IPV(4,4)=8
      NIPV(5)=4
      IPV(5,1)=1
      IPV(5,2)=4
      IPV(5,3)=8
      IPV(5,4)=5
      NIPV(6)=4
      IPV(6,1)=6
      IPV(6,2)=5
      IPV(6,3)=8
      IPV(6,4)=7
      VERTP(1,1)=0.5D0
      VERTP(1,2)=0.75D0
      VERTP(1,3)=D1
      VERTP(2,1)=0.5D0
      VERTP(2,2)=0.75D0
      VERTP(2,3)=D0
      VERTP(3,1)=D1
      VERTP(3,2)=D1
      VERTP(3,3)=D0
      VERTP(4,1)=D1
      VERTP(4,2)=D1
      VERTP(4,3)=D1
      VERTP(5,1)=D0
      VERTP(5,2)=D0
      VERTP(5,3)=D1
      VERTP(6,1)=D0
      VERTP(6,2)=D0
      VERTP(6,3)=D0
      VERTP(7,1)=D0
      VERTP(7,2)=D1
      VERTP(7,3)=D0
      VERTP(8,1)=D0
      VERTP(8,2)=D1
      VERTP(8,3)=D1

      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO

      RETURN
      END
c------------------------ END OF NCHEXAHEMESH ------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                     SMALL STELLATED DODECAMESH                      c
c---------------------------------------------------------------------c
c      This version dated: November 5, 2015 by Lopez and Hernandez    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCDODECAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)

      CALL DODECAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      NTS0=NTS
      NTSI=NTS
      NTP0=NTP
      A=0.5D0
c.. Face centroid      
      DO IS=1,NTS
         XC=0.0
         YC=0.0
         ZC=0.0
         DO IV=1,NIPV(IS)
            IP=IPV(IS,IV)
            XC=XC+VERTP(IP,1)
            YC=YC+VERTP(IP,2)
            ZC=ZC+VERTP(IP,3)
         END DO
         NTP0=NTP0+1
         VERTP(NTP0,1)=XC/NIPV(IS)+A*XNS(IS)
         VERTP(NTP0,2)=YC/NIPV(IS)+A*YNS(IS)
         VERTP(NTP0,3)=ZC/NIPV(IS)+A*ZNS(IS)
         
         DO IV=1,NIPV(IS)
            IV2=IV+1
            IF(IV.EQ.NIPV(IS)) IV2=1
            NTS0=NTS0+1
            NIPV(NTS0)=3
            IPV(NTS0,1)=IPV(IS,IV)
            IPV(NTS0,2)=IPV(IS,IV2)
            IPV(NTS0,3)=NTP0
         END DO
      END DO

      NTS=NTS0-NTSI
      DO IS=1,NTS
         NIPV(IS)=NIPV(NTSI+IS)
         DO IV=1,NIPV(IS)
            IPV(IS,IV)=IPV(NTSI+IS,IV)
         END DO
      END DO
      NTP=NTP0
      NTV=NTP
      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5D0
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO
      RETURN
      END
c------------------------ END OF NCDODECAMESH ------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                      SMALL STELLATED ICOSAMESH                      c
c---------------------------------------------------------------------c
c      This version dated: November 5, 2015 by Lopez and Hernandez    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCICOSAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)

      CALL ICOSAMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      NTS0=NTS
      NTSI=NTS
      NTP0=NTP
      A=1D0
c.. Face centroid      
      DO IS=1,NTS
         XC=0.0
         YC=0.0
         ZC=0.0
         DO IV=1,NIPV(IS)
            IP=IPV(IS,IV)
            XC=XC+VERTP(IP,1)
            YC=YC+VERTP(IP,2)
            ZC=ZC+VERTP(IP,3)
         END DO
         NTP0=NTP0+1
         VERTP(NTP0,1)=XC/NIPV(IS)+A*XNS(IS)
         VERTP(NTP0,2)=YC/NIPV(IS)+A*YNS(IS)
         VERTP(NTP0,3)=ZC/NIPV(IS)+A*ZNS(IS)
         
         DO IV=1,NIPV(IS)
            IV2=IV+1
            IF(IV.EQ.NIPV(IS)) IV2=1
            NTS0=NTS0+1
            NIPV(NTS0)=3
            IPV(NTS0,1)=IPV(IS,IV)
            IPV(NTS0,2)=IPV(IS,IV2)
            IPV(NTS0,3)=NTP0
         END DO
      END DO

      NTS=NTS0-NTSI
      DO IS=1,NTS
         NIPV(IS)=NIPV(NTSI+IS)
         DO IV=1,NIPV(IS)
            IPV(IS,IV)=IPV(NTSI+IS,IV)
         END DO
      END DO
      NTP=NTP0
      NTV=NTP
      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5D0
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO
      RETURN
      END
c------------------------- END OF NCICOSAMESH ------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                           NCHOLLOWEDCUBE                            c
c---------------------------------------------------------------------c
c      This version dated: October 16, 2018 by Lopez and Hernandez    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCHOLLOWEDCUBE(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
c Unit-length cube with a half-length cubic hollow in its center:
      D0=0.0D0      
      D1=1.0D0
      D14=0.25D0

      XNS(1)=D1
      YNS(1)=D0
      ZNS(1)=D0
      XNS(2)=D0
      YNS(2)=-D1
      ZNS(2)=D0
      XNS(3)=D0
      YNS(3)=D0
      ZNS(3)=-D1
      XNS(4)=D0
      YNS(4)=D1
      ZNS(4)=D0
      XNS(5)=D0
      YNS(5)=D0
      ZNS(5)=D1
      XNS(6)=-D1
      YNS(6)=D0
      ZNS(6)=D0
      NTS=12
      NTV=16
      NTP=NTV
      NIPV(1)=4
      IPV(1,1)=1
      IPV(1,2)=2
      IPV(1,3)=3
      IPV(1,4)=4
      NIPV(2)=4
      IPV(2,1)=2
      IPV(2,2)=1
      IPV(2,3)=5
      IPV(2,4)=6
      NIPV(3)=4
      IPV(3,1)=3
      IPV(3,2)=2
      IPV(3,3)=6
      IPV(3,4)=7
      NIPV(4)=4
      IPV(4,1)=4
      IPV(4,2)=3
      IPV(4,3)=7
      IPV(4,4)=8
      NIPV(5)=4
      IPV(5,1)=1
      IPV(5,2)=4
      IPV(5,3)=8
      IPV(5,4)=5
      NIPV(6)=4
      IPV(6,1)=6
      IPV(6,2)=5
      IPV(6,3)=8
      IPV(6,4)=7
      DO IS=1,6
         XNS(IS+6)=-XNS(IS)
         YNS(IS+6)=-YNS(IS)
         ZNS(IS+6)=-ZNS(IS)
         NIPV(IS+6)=NIPV(IS)
         DO IV=1,4
            IV2=4-IV+1
            IPV(IS+6,IV)=IPV(IS,IV2)+8
         END DO
      END DO
c       7/----------/3
c       /|         /|
c      / |        / |
c    8/__|______4/  | 
c     |  |       |  |
c     |  /6------|--/2
c     | /        | /
c     |/_________|/
c     5           1
      VERTP(1,1)=D1
      VERTP(1,2)=D0
      VERTP(1,3)=D1
      VERTP(2,1)=D1
      VERTP(2,2)=D0
      VERTP(2,3)=D0
      VERTP(3,1)=D1
      VERTP(3,2)=D1
      VERTP(3,3)=D0
      VERTP(4,1)=D1
      VERTP(4,2)=D1
      VERTP(4,3)=D1
      VERTP(5,1)=D0
      VERTP(5,2)=D0
      VERTP(5,3)=D1
      VERTP(6,1)=D0
      VERTP(6,2)=D0
      VERTP(6,3)=D0
      VERTP(7,1)=D0
      VERTP(7,2)=D1
      VERTP(7,3)=D0
      VERTP(8,1)=D0
      VERTP(8,2)=D1
      VERTP(8,3)=D1

      VERTP(9,1)=D1-D14
      VERTP(9,2)=D0+D14
      VERTP(9,3)=D1-D14
      VERTP(10,1)=D1-D14
      VERTP(10,2)=D0+D14
      VERTP(10,3)=D0+D14
      VERTP(11,1)=D1-D14
      VERTP(11,2)=D1-D14
      VERTP(11,3)=D0+D14
      VERTP(12,1)=D1-D14
      VERTP(12,2)=D1-D14
      VERTP(12,3)=D1-D14
      VERTP(13,1)=D0+D14
      VERTP(13,2)=D0+D14
      VERTP(13,3)=D1-D14
      VERTP(14,1)=D0+D14
      VERTP(14,2)=D0+D14
      VERTP(14,3)=D0+D14
      VERTP(15,1)=D0+D14
      VERTP(15,2)=D1-D14
      VERTP(15,3)=D0+D14
      VERTP(16,1)=D0+D14
      VERTP(16,2)=D1-D14
      VERTP(16,3)=D1-D14
      RETURN
      END
c---------------------- END OF NCHOLLOWEDCUBE ------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             DRILLEDCUBE                             c
c                  Non-simply connected polyhedron                    c
c---------------------------------------------------------------------c
c      This version dated: January 29, 2019 by Lopez and Hernandez    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE DRILLEDCUBE(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
      CALL CUBICMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      D0=0.0D0
      D1=1.0D0

      XNS(NTS+1)=-D1
      YNS(NTS+1)=-D0
      ZNS(NTS+1)=-D0
      XNS(NTS+2)=-D0
      YNS(NTS+2)=D1
      ZNS(NTS+2)=-D0
      XNS(NTS+3)=-D0
      YNS(NTS+3)=-D0
      ZNS(NTS+3)=D1
      XNS(NTS+4)=-D0
      YNS(NTS+4)=-D1
      ZNS(NTS+4)=-D0
      XNS(NTS+5)=-D0
      YNS(NTS+5)=-D0
      ZNS(NTS+5)=-D1
      XNS(NTS+6)=D1
      YNS(NTS+6)=-D0
      ZNS(NTS+6)=-D0
      
      NIPV(NTS+1)=4
      IPV(NTS+1,4)=NTP+1
      IPV(NTS+1,3)=NTP+2
      IPV(NTS+1,2)=NTP+3
      IPV(NTS+1,1)=NTP+4
      NIPV(NTS+2)=4
      IPV(NTS+2,4)=NTP+2
      IPV(NTS+2,3)=NTP+1
      IPV(NTS+2,2)=NTP+5
      IPV(NTS+2,1)=NTP+6
      NIPV(NTS+3)=4
      IPV(NTS+3,4)=NTP+3
      IPV(NTS+3,3)=NTP+2
      IPV(NTS+3,2)=NTP+6
      IPV(NTS+3,1)=NTP+7
      NIPV(NTS+4)=4
      IPV(NTS+4,4)=NTP+4
      IPV(NTS+4,3)=NTP+3
      IPV(NTS+4,2)=NTP+7
      IPV(NTS+4,1)=NTP+8
      NIPV(NTS+5)=4
      IPV(NTS+5,4)=NTP+1
      IPV(NTS+5,3)=NTP+4
      IPV(NTS+5,2)=NTP+8
      IPV(NTS+5,1)=NTP+5
      NIPV(NTS+6)=4
      IPV(NTS+6,4)=NTP+6
      IPV(NTS+6,3)=NTP+5
      IPV(NTS+6,2)=NTP+8
      IPV(NTS+6,1)=NTP+7

      D02=0.25D0
      D12=0.75D0

      VERTP(NTP+1,1)=D12
      VERTP(NTP+1,2)=D0
      VERTP(NTP+1,3)=D12
      VERTP(NTP+2,1)=D12
      VERTP(NTP+2,2)=D0
      VERTP(NTP+2,3)=D02
      VERTP(NTP+3,1)=D12
      VERTP(NTP+3,2)=D1
      VERTP(NTP+3,3)=D02
      VERTP(NTP+4,1)=D12
      VERTP(NTP+4,2)=D1
      VERTP(NTP+4,3)=D12
      VERTP(NTP+5,1)=D02
      VERTP(NTP+5,2)=D0
      VERTP(NTP+5,3)=D12
      VERTP(NTP+6,1)=D02
      VERTP(NTP+6,2)=D0
      VERTP(NTP+6,3)=D02
      VERTP(NTP+7,1)=D02
      VERTP(NTP+7,2)=D1
      VERTP(NTP+7,3)=D02
      VERTP(NTP+8,1)=D02
      VERTP(NTP+8,2)=D1
      VERTP(NTP+8,3)=D12

      NTS=NTS+6
      NTV=NTV+8
      NTP=NTV

      RETURN
      END
c------------------------ END OF DRILLEDCUBE -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             ZIGZAGMESH                              c
c---------------------------------------------------------------------c
c      This version dated: January 29, 2019 by Lopez and Hernandez    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE ZIGZAGMESH(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
c. Change NZIGS to modify the number of zig-zag sections
      NZIGS=5
      DOFF=0.1D0
      D0=0.0D0
      D1=1.0D0
      NTP=4*NZIGS
      NTV=NTP
      NTS=2*NZIGS+2
      DO IV=1,NZIGS
         VERTP(IV,1)=D1*(IV-1)
         VERTP(IV,2)=DOFF+D1*MOD(IV-1,2)
         VERTP(IV,3)=D0
         VERTP(IV+NZIGS,1)=D1*(NZIGS-(IV-1)-1)
         VERTP(IV+NZIGS,2)=D1*MOD(NZIGS-(IV-1)-1,2)
         VERTP(IV+NZIGS,3)=D0
         VERTP(IV+2*NZIGS,1)=D1*(IV-1)
         VERTP(IV+2*NZIGS,2)=DOFF+D1*MOD(IV-1,2)
         VERTP(IV+2*NZIGS,3)=D1
         VERTP(IV+3*NZIGS,1)=D1*(NZIGS-(IV-1)-1)
         VERTP(IV+3*NZIGS,2)=D1*MOD(NZIGS-(IV-1)-1,2)
         VERTP(IV+3*NZIGS,3)=D1
      END DO
      NIPV(1)=2*NZIGS
      NIPV(2)=2*NZIGS
      DO IV=1,2*NZIGS
         IPV(1,IV)=IV
         IPV(2,IV)=4*NZIGS-(IV-1)         
      END DO
      IV=1
      NIPV(3)=4
      IPV(3,1)=IV
      IPV(3,2)=2*NZIGS-(IV-1)
      IPV(3,3)=4*NZIGS-(IV-1)
      IPV(3,4)=IV+2*NZIGS
      IV=NZIGS
      NIPV(4)=4
      IPV(4,4)=IV
      IPV(4,3)=2*NZIGS-(IV-1)
      IPV(4,2)=4*NZIGS-(IV-1)
      IPV(4,1)=IV+2*NZIGS
      DO IV=1,NZIGS-1
         NIPV(IV+4)=4
         IPV(IV+4,1)=IV
         IPV(IV+4,2)=IV+2*NZIGS
         IPV(IV+4,3)=IV+2*NZIGS+1
         IPV(IV+4,4)=IV+1
         NIPV(IV+4+(NZIGS-1))=4
         IPV(IV+4+(NZIGS-1),1)=2*NZIGS-(IV-1)
         IPV(IV+4+(NZIGS-1),2)=2*NZIGS-(IV-1)-1
         IPV(IV+4+(NZIGS-1),3)=2*NZIGS-(IV-1)-1+2*NZIGS
         IPV(IV+4+(NZIGS-1),4)=2*NZIGS-(IV-1)+2*NZIGS
      END DO
      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5D0
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO
      RETURN
      END
c------------------------ END OF ZIGZAGMESH --------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            VOFTOOLSLOGO                             c
c---------------------------------------------------------------------c
c      This version dated: April 15, 2019 by Lopez and Hernandez      c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE VOFTOOLSLOGO(IPV,NIPV,NTP,NTS,NTV,VERTP,XNS,YNS,ZNS)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NS,NV),NIPV(NS),VERTP(NV,3),XNS(NS),
     -     YNS(NS),ZNS(NS)
      NTS=0
      NTP=0
      DEPTH=1.0D0
c.. V letter:      
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=1.0D0
      NTP=NTP+1
      VERTP(NTP,1)=0.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=2.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=4.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=6.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=5.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=3.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=1.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO

c..O letter:
      DISP=6.2D0
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=1.0D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO

      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO
c..
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=1.0D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO
c..F letter:
      DISP=DISP+4.6D0
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=1.0D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=2.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=2.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=3.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=3.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.4*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.4*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO

c..T letter:
      DISP=DISP+3.45D0
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=1.0D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+2.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.65*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.65*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.35*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.35*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+2.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO
c..o letter:
      DISP=DISP+4.2D0
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=0.7D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO
c..
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=0.7D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO

c..o letter:
      DISP=DISP+3.3D0
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=0.7D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO
c..
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=0.7D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO
c..l letter:
      DISP=DISP+2.7D0
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=1.0D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+2.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+2.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO
c..s letter:
      DISP=DISP+2.5D0
      NTSI=NTS+1
      NTPI=NTP+1
      FACTOR=0.7D0
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=0.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=3.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=3.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+1.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=4.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+4.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=5.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=2.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=2.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+3.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTP=NTP+1
      VERTP(NTP,1)=DISP+0.0*FACTOR
      VERTP(NTP,2)=1.0*FACTOR
      VERTP(NTP,3)=DEPTH
      NTS=NTS+1
      NIPV(NTS)=NTP-NTPI+1
      NTPL=NTP-NTPI+1
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+IV-1
      END DO
      DO I=1,NTPL
         NTP=NTP+1
         VERTP(NTP,1)=VERTP(NTPI+I-1,1)
         VERTP(NTP,2)=VERTP(NTPI+I-1,2)
         VERTP(NTP,3)=0.0
      END DO
      DO IS=2,NTPL+1
         NTS=NTS+1
         NIPV(NTS)=4
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,1)=(NTPI-1)+1
         ELSE
            IPV(NTS,1)=(NTPI-1)+IS
         ENDIF
         IPV(NTS,2)=(NTPI-1)+IS-1
         IPV(NTS,3)=(NTPI-1)+IS+NTPL-1
         IF(IS.EQ.(NTPL+1)) THEN
            IPV(NTS,4)=(NTPI-1)+1+NTPL
         ELSE
            IPV(NTS,4)=(NTPI-1)+IS+NTPL
         END IF
      END DO
      NTS=NTS+1
      NIPV(NTS)=NIPV(NTSI)
      DO IV=1,NIPV(NTS)
         IPV(NTS,IV)=NTPI+NTPL*2-IV
      END DO      
c............      
      NTV=NTP
      DO IS=1,NTS
         IP1=IPV(IS,1)
         IP2=IPV(IS,2)
         IP3=IPV(IS,3)         
         XV1=VERTP(IP2,1)-VERTP(IP1,1)
         YV1=VERTP(IP2,2)-VERTP(IP1,2)
         ZV1=VERTP(IP2,3)-VERTP(IP1,3)
         XV2=VERTP(IP3,1)-VERTP(IP2,1)
         YV2=VERTP(IP3,2)-VERTP(IP2,2)
         ZV2=VERTP(IP3,3)-VERTP(IP2,3)
         XN=YV1*ZV2-ZV1*YV2
         YN=ZV1*XV2-XV1*ZV2
         ZN=XV1*YV2-YV1*XV2
         DMOD=(XN**2.0D0+YN**2.0D0+ZN**2.0D0)**0.5D0
         XNS(IS)=XN/DMOD
         YNS(IS)=YN/DMOD
         ZNS(IS)=ZN/DMOD
      END DO
      RETURN
      END
c----------------------- END OF VOFTOOLSLOGO -------------------------c
c---------------------------------------------------------------------c
      
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            SQUAREMESH                               c
c---------------------------------------------------------------------c
c      This version dated: April 10, 2007 by Lopez and Hernandez      c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE SQUAREMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=4
      NTP=4
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=1.0
      VERTP(2,2)=0.0
      VERTP(3,1)=1.0
      VERTP(3,2)=1.0
      VERTP(4,1)=0.0
      VERTP(4,2)=1.0
      RETURN
      END
c------------------------- END OF SQUAREMESH -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                       HEXAGOMESH (Regular hexagon)                  c
c---------------------------------------------------------------------c
c      This version dated: FEB 25, 2015 by Lopez and Hernandez        c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE HEXAGOMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=6
      NTP=6
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.5
      VERTP(1,2)=0.0
      VERTP(2,1)=0.9330127
      VERTP(2,2)=0.25
      VERTP(3,1)=0.9330127
      VERTP(3,2)=0.75
      VERTP(4,1)=0.5
      VERTP(4,2)=1.0
      VERTP(5,1)=0.066987298
      VERTP(5,2)=0.75
      VERTP(6,1)=0.066987298
      VERTP(6,2)=0.25
      RETURN
      END
c------------------------- END OF HEXAGOMESH -------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            TRIANGLEMESH                             c
c---------------------------------------------------------------------c
c      This version dated: January, 2015 by Lopez and Hernandez       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE TRIANGLEMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=3
      NTP=3
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=0.72
      VERTP(2,2)=0.13
      VERTP(3,1)=1.0
      VERTP(3,2)=1.0
      RETURN
      END
c------------------------ END OF TRIANGLEMESH ------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                          QUADRANGLEMESH                             c
c---------------------------------------------------------------------c
c      This version dated: January, 2015 by Lopez and Hernandez       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE QUADRANGLEMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=4
      NTP=4
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=1.0
      VERTP(2,2)=0.13
      VERTP(3,1)=0.72
      VERTP(3,2)=1.0
      VERTP(4,1)=0.13
      VERTP(4,2)=0.56
      RETURN
      END
c----------------------- END OF QUADRANGLEMESH -----------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                           PENTAGONMESH                              c
c---------------------------------------------------------------------c
c      This version dated: January, 2015 by Lopez and Hernandez       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE PENTAGONMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=5
      NTP=5
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.15
      VERTP(2,1)=0.09
      VERTP(2,2)=0.0
      VERTP(3,1)=1.0
      VERTP(3,2)=0.13
      VERTP(4,1)=0.16
      VERTP(4,2)=1.0
      VERTP(5,1)=0.04
      VERTP(5,2)=0.77
      RETURN
      END
c------------------------ END OF PENTAGONMESH ------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                            HEXAGONMESH                              c
c---------------------------------------------------------------------c
c      This version dated: January, 2015 by Lopez and Hernandez       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE HEXAGONMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=6
      NTP=6
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=0.66
      VERTP(2,2)=0.03
      VERTP(3,1)=1.0
      VERTP(3,2)=0.22
      VERTP(4,1)=0.9
      VERTP(4,2)=0.77
      VERTP(5,1)=0.72
      VERTP(5,2)=1.0
      VERTP(6,1)=0.33
      VERTP(6,2)=0.86
      RETURN
      END
c------------------------- END OF HEXAGONMESH ------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                      NON-CONVEX QUADRANGLEMESH                      c
c---------------------------------------------------------------------c
c      This version dated: January, 2015 by Lopez and Hernandez       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCQUADRANGLEMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=4
      NTP=4
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=1.0
      VERTP(2,2)=0.13
      VERTP(3,1)=0.72
      VERTP(3,2)=1.0
      VERTP(4,1)=0.75
      VERTP(4,2)=0.56
      RETURN
      END
c--------------------- END OF NCQUADRANGLEMESH -----------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                        NON-CONVEX PENTAGONMESH                      c
c---------------------------------------------------------------------c
c      This version dated: January, 2015 by Lopez and Hernandez       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCPENTAGONMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=5
      NTP=5
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=0.49
      VERTP(2,2)=0.22
      VERTP(3,1)=1.0
      VERTP(3,2)=0.13
      VERTP(4,1)=0.16
      VERTP(4,2)=1.0
      VERTP(5,1)=0.04
      VERTP(5,2)=0.77
      RETURN
      END
c----------------------- END OF NCPENTAGONMESH -----------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                         NON-CONVEX HEXAGONMESH                      c
c---------------------------------------------------------------------c
c      This version dated: January, 2015 by Lopez and Hernandez       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCHEXAGONMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=6
      NTP=6
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=0.33
      VERTP(2,2)=0.43
      VERTP(3,1)=1.0
      VERTP(3,2)=0.22
      VERTP(4,1)=0.9
      VERTP(4,2)=0.77
      VERTP(5,1)=0.72
      VERTP(5,2)=1.0
      VERTP(6,1)=0.33
      VERTP(6,2)=0.86
      RETURN
      END
c------------------------ END OF NCHEXAGONMESH -----------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                    NON-CONVEX STELLATED HEXAGON                     c
c---------------------------------------------------------------------c
c      This version dated: January, 2015 by Lopez and Hernandez       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCSHEXAGONMESH(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=12
      NTP=12
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.2
      VERTP(2,1)=0.25
      VERTP(2,2)=0.25
      VERTP(3,1)=0.4
      VERTP(3,2)=0.0
      VERTP(4,1)=0.55
      VERTP(4,2)=0.3
      VERTP(5,1)=0.8
      VERTP(5,2)=0.35
      VERTP(6,1)=0.65
      VERTP(6,2)=0.55
      VERTP(7,1)=0.8
      VERTP(7,2)=0.75
      VERTP(8,1)=0.5
      VERTP(8,2)=0.7
      VERTP(9,1)=0.35
      VERTP(9,2)=0.9
      VERTP(10,1)=0.25
      VERTP(10,2)=0.7
      VERTP(11,1)=-0.1
      VERTP(11,2)=0.55
      VERTP(12,1)=0.2
      VERTP(12,2)=0.45
      RETURN
      END
c----------------------- END OF NCSHEXAGONMESH -----------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                           HOLLOWED SQUARE                           c
c---------------------------------------------------------------------c
c      This version dated: November 4, 2015 by Lopez and Hernandez    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCHOLLOWEDSQUARE(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=10
      NTP=10
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=1.0
      VERTP(2,2)=0.0
      VERTP(3,1)=1.0
      VERTP(3,2)=1.0
      VERTP(4,1)=0.0
      VERTP(4,2)=1.0
      VERTP(5,1)=VERTP(1,1)
      VERTP(5,2)=VERTP(1,2)
c.. hollow
      VERTP(6,1)=0.25
      VERTP(6,2)=0.25
      VERTP(7,1)=0.25
      VERTP(7,2)=0.75
      VERTP(8,1)=0.75
      VERTP(8,2)=0.75
      VERTP(9,1)=0.75
      VERTP(9,2)=0.25
      VERTP(10,1)=VERTP(6,1)
      VERTP(10,2)=VERTP(6,2)

      RETURN
      END
c----------------------- END OF NCHOLLOWEDSQUARE ---------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                   NON-CONVEX MULTI-SQUARE CELL                      c
c---------------------------------------------------------------------c
c      This version dated: November 4, 2015 by Lopez and Hernandez    c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c 
      SUBROUTINE NCMULTISQUARE(IPV,NTP,NTV,VERTP)
      INCLUDE "dim.h"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPV(NV),VERTP(NV,2)
      NTV=15
      NTP=15
      DO I=1,NTV
         IPV(I)=I
      END DO
      VERTP(1,1)=0.0
      VERTP(1,2)=0.0
      VERTP(2,1)=1.0
      VERTP(2,2)=0.0
      VERTP(3,1)=1.0
      VERTP(3,2)=1.0
      VERTP(4,1)=0.0
      VERTP(4,2)=1.0
      VERTP(5,1)=VERTP(1,1)
      VERTP(5,2)=VERTP(1,2)
c.. hollow
      VERTP(6,1)=0.25
      VERTP(6,2)=0.25
      VERTP(7,1)=0.25
      VERTP(7,2)=0.75
      VERTP(8,1)=0.75
      VERTP(8,2)=0.75
      VERTP(9,1)=0.75
      VERTP(9,2)=0.25
      VERTP(10,1)=VERTP(6,1)
      VERTP(10,2)=VERTP(6,2)
c.. liquid inside hollow
      VERTP(11,1)=0.4
      VERTP(11,2)=0.4
      VERTP(12,1)=0.6
      VERTP(12,2)=0.4
      VERTP(13,1)=0.6
      VERTP(13,2)=0.6
      VERTP(14,1)=0.4
      VERTP(14,2)=0.6
      VERTP(15,1)=VERTP(11,1)
      VERTP(15,2)=VERTP(11,2)
      RETURN
      END
c------------------------ END OF NCMULTISQUARE -----------------------c
c---------------------------------------------------------------------c

