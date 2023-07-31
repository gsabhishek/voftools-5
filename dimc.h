//---------------------------------------------------------------------//
//---------------------------------------------------------------------//
//                             dimc.h                                  //
//---------------------------------------------------------------------//
//            Copyright (C) 2020 J. Lopez and J. Hernandez             //
//---------------------------------------------------------------------//
//      Parameters used to dimension the arrays that are passed to     //
//      the VOFTools routines from the C test program.                 //
//---------------------------------------------------------------------//
//      IMPORTANT ADVISE: the values of the parameters 'ns' and 'nv'   //
//      ----------------  must coincide with the values of the same    //
//                        parametes given in the file "dim.h" which is //
//                        used to construct the VOFTools library. If   //
//                        you need to modify these values (to use      //
//                        cells with a higher number of faces or       //
//                        vertices, for example), the VOFTools library //
//                        must be recompiled by changing accordingly   //
//                        the corresponding parameters in the file     //
//                        "dim.h".                                     //
//---------------------------------------------------------------------//
// This file is part of VOFTools.                                      //
//                                                                     //
// VOFTools is free software: you can redistribute it and/or           //
// modify it under the terms of the GNU General Public License         //
// as published by the Free Software Foundation, either version 3 of   //
// the License, or (at your option) any later version.                 //
//                                                                     //
// VOFTools is distributed in the hope that it will be useful,         //
// but WITHOUT ANY WARRANTY; without even the implied warranty of      //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       //
// GNU Lesser General Public License for more details.                 //
//                                                                     //
// You should have received a copy of the GNU General Public License   //
// along with VOFTools. If not, see <http://www.gnu.org/licenses/>.    //
//---------------------------------------------------------------------//
#define ns 200
#define nv 240
