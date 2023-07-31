c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                           VOFVARDEF                                 c
c... Define the main parameters of the test cases                     c
c---------------------------------------------------------------------c
c On return:                                                          c
c===========                                                          c
c F        = material volume/area fraction                            c
c ICELLTYPE= cell geometry                                            c
c ISHAPE   = material body shape                                      c
c NC       = subdivision number in the volume fraction cell           c
c            initialization                                           c
c TOL      = prescribed positive tolerance for the distance to the    c
c            interface for the initialization procedure               c
c XNC,..   = unit-length normal vector of the interface plane         c
c XP,..    = coordinates of the point P from which the distance to    c
c            the reconstructed interface is calculated                c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      SUBROUTINE VOFVARDEF(F,ICELLTYPE,ISHAPE,NC,TOL,XNC,XP,YNC,YP,ZNC,
     -     ZP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LEN=100)
      CHARACTER (LEN) :: LINE
      OPEN(1,FILE='vofvardef',STATUS='OLD',ACTION='READ')
      READ(1,'(A)') LINE
      READ(1,*) ICELLTYPE
      WRITE(6,*) 'ICELLTYPE:',ICELLTYPE
      READ(1,'(A)') LINE
      READ(1,*) ISHAPE
      WRITE(6,*) 'ISHAPE:',ISHAPE
      READ(1,'(A)') LINE
      READ(1,*) F
      WRITE(6,*) 'F:',F
      READ(1,'(A)') LINE
      READ(1,*) XNC
      WRITE(6,*) 'XNC:',XNC
      READ(1,'(A)') LINE
      READ(1,*) YNC
      WRITE(6,*) 'YNC:',YNC
      READ(1,'(A)') LINE
      READ(1,*) ZNC
      WRITE(6,*) 'ZNC:',ZNC
      READ(1,'(A)') LINE
      READ(1,*) XP
      WRITE(6,*) 'XP:',XP
      READ(1,'(A)') LINE
      READ(1,*) YP
      WRITE(6,*) 'YP:',YP
      READ(1,'(A)') LINE
      READ(1,*) ZP
      WRITE(6,*) 'ZP:',ZP
      READ(1,'(A)') LINE
      READ(1,*) NC
      WRITE(6,*) 'NC:',NC
      READ(1,'(A)') LINE
      READ(1,*) TOL
      WRITE(6,*) 'TOL:',TOL
      CLOSE(1)      
      RETURN
      END
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             FUNC2D1                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y      = coordinates of the point where VALUE is computed         c
c On return:                                                          c
c===========                                                          c
c FUNC2D1  = VALUE of the implicit interface shape function:          c
c            > 0 (inside the interface), < 0 (outside the             c
c            interface); = 0 (on the interface)                       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      FUNCTION FUNC2D1(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION FUNC2D1,X,Y
c* Circle with radious 0.325 centered at (0.5,0.5):      
c* Exact area encolsed by the interface=PI*0.325**2      
      FUNC2D1=-1.0*((X-0.5)**2+(Y-0.5)**2-0.325**2)
      RETURN
      END 
c------------------------- END OF FUNC2D1 ----------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             FUNC2D2                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y      = coordinates of the point where VALUE is computed         c
c On return:                                                          c
c===========                                                          c
c FUNC2D2  = VALUE of the implicit interface shape function:          c
c            > 0 (inside the interface), < 0 (outside the             c
c            interface); = 0 (on the interface)                       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      FUNCTION FUNC2D2(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION FUNC2D2,X,Y
c* Ellipse with semi-major axis 0.5, semi-minor axis 0.2 and
c* centered at (0.5,0.5):      
c* Exact area encolsed by the interface=PI*0.35*0.25      
      FUNC2D2=1.0-(((X-0.5)/0.5)**2+((Y-0.5)/0.2)**2)
      RETURN
      END 
c------------------------- END OF FUNC2D2 ----------------------------c
c---------------------------------------------------------------------c
      
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             FUNC3D1                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y,Z    = coordinates of the point where VALUE is computed         c
c On return:                                                          c
c===========                                                          c
c FUNC3D1  = VALUE of the implicit interface shape function:          c
c            > 0 (inside the interface), < 0 (outside the             c
c            interface); = 0 (on the interface)                       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      FUNCTION FUNC3D1(X,Y,Z)
      IMPLICIT NONE
      DOUBLE PRECISION FUNC3D1,X,Y,Z
c* Sphere with radious 0.325 centered at (0.5,0.5,0.5):      
c* Exact volume encolsed by the interface=(4/3)*PI*0.325**3      
      FUNC3D1=-1.0*((X-0.5)**2+(Y-0.5)**2+(Z-0.5)**2-0.325**2)
      RETURN
      END 
c------------------------- END OF FUNC3D1 ----------------------------c
c---------------------------------------------------------------------c

c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
c                             FUNC3D2                                 c
c---------------------------------------------------------------------c
c On entry:                                                           c
c==========                                                           c
c X,Y,Z    = coordinates of the point where VALUE is computed         c
c On return:                                                          c
c===========                                                          c
c FUNC3D2  = VALUE of the implicit interface shape function:          c
c            > 0 (inside the interface), < 0 (outside the             c
c            interface); = 0 (on the interface)                       c
c---------------------------------------------------------------------c
c---------------------------------------------------------------------c
      FUNCTION FUNC3D2(X,Y,Z)
      IMPLICIT NONE
      DOUBLE PRECISION FUNC3D2,X,Y,Z
c* Torus with major radius 2/3, minor radius 1/3 and centered at (0.5,0.5,0.5):
c* Exact volume encolsed by the interface=2*PI**2*(2/3)*(1/3)**2      
      FUNC3D2=(1.0/3.0)**2-((2.0/3.0)-((X-0.5)**2+(Z-0.5)**2)**0.5)**2-
     -     (Y-0.5)**2
      RETURN
      END 
c------------------------- END OF FUNC3D2 ----------------------------c
c---------------------------------------------------------------------c
