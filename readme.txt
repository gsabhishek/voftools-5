-----------------------------------------------------------------------
|                           V O F T o o l s                           |
|                                                                     |
|	 A package of FORTRAN subroutines with analytical and         |
|                                                                     |  
|         geometrical tools for 2D/3D VOF methods in general          |
|                                                                     |
|                     grids and Cartesian geometry                    |
|                                                                     |
|                                                                     |
|                       (Version 5, January 2020)                     |
|                                                                     |
|				                                      |
|             Copyright (C) 2020 J. Lopez and J. Hernandez            |
|                                                                     |
|                                                                     |
-----------------------------------------------------------------------

-----------------------------------------------------------------------

VOFTools is a collection of FORTRAN subroutines with analytical and 
geometrical tools for 2D/3D Volume of Fluid (VOF) methods in
general grids and Cartesian geometry.

This library is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
                                                                     
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
                                                                     
You should have received a copy of the GNU Lesser General Public
License along with this library.  If not, see
<http://www.gnu.org/licenses/>.

References about VOFTools are the following:

[1] J. Lopez and J. Hernandez, Analytical and geometrical tools for 3D
volume of fluid methods in general grids, Journal of Computational
Physics 227 (2008) 5939-5948. 

[2] J. López, P. Gómez, J. Hernández, F. Faura, A two-grid adaptive
volume of fluid approach for dendritic solidification, Computers &
Fluids, 86 (2013) 326-342. 

[3] J. Lopez, J. Hernandez, P. Gómez, F. Faura, A new volume 
conservation enforcement method for PLIC reconstruction in general 
convex grids, Journal of Computational Physics 316 (2016) 338-359.

[4] J. Lopez, P. Gómez, C. Zanzi, J. Hernandez, F. Faura, Application of 
non-convex analytic and geometric tools to a PLIC-VOF method, Proceedings 
of the ASME 2016 International Mechanical Engineering Congress and 
Exposition, November 11-17, 2016, Phoenix, Arizona, IMECE2016-67409.

[5] J. López, J. Hernández, P. Gómez and F. Faura, VOFTools - A software
package of calculation tools for volume of fluid methods using general
convex grids, Computer Physics Communications 223 (2018) 45-54.

[6] J. López, J. Hernández, P. Gómez and F. Faura, Non-convex analytical
and geometrical tools for volume truncation, initialization and
conservation enforcement in VOF methods, Journal of Computational Physics 
392 (2019) 666-693.

[7] J. López, J. Hernández, P. Gómez, C. Zanzi, R. Zamora, VOFTools 3.2: 
Added VOF functionality to initialize the liquid volume fraction in 
general convex cells, Computer Physics Communications 245 (2019) 106859.

For more information contact joaquin.lopez@upct.es 

-----------------------------------------------------------------------
                    V O F T o o l s  VERSION 5
-----------------------------------------------------------------------

In this directory you will find the following files: 

voftools.f      ---> Contains the source code of the analytical and 
                     geometrical tools

uservoftools.f  ---> Contains the source code of user-defined routines
                     and functions.

mesh.f          ---> Contains the definitions of different cell types

dim.h           ---> Contains the array dimensions for FORTRAN codes

dimc.h          ---> Contains the array dimensions for C codes

cvoftools.h     ---> Contains the declaration of the subroutines of the
                     VOFTools library to be used in C programs

cuservoftools.h ---> Contains the declaration of the user-defined routines
                     and functions included in the uservoftools.f file to
		     be used in C programs

cmesh.h         ---> Contains the declaration of the subroutines used to
                     define the cell geometries considered in the test
                     programs to be used in C programs

test2d.f        ---> Contains the 2D test program in FORTRAN

test2d.c        ---> Contains the 2D test program in C

test3d.f        ---> Contains the 3D test program in FORTRAN

test3d.c        ---> Contains the 3D test program in C

vofvardef       ---> Contains the input data for the test programs

Makefile.#      ---> Constructs in different platforms the VOFTools library 
                     and makes executable files for test programs in C and 
                     Fortran

change.log      ---> Contains a record of all notable changes made to the 
                     routines

COPYING         ---> Copy of the GNU General Pulic License, Version 3

