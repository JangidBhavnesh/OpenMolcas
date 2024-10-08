!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine Get_Info_Static()
!***********************************************************************
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January 1992                                             *
!***********************************************************************

use Symmetry_Info, only: Symmetry_Info_Get
use Sizes_of_Seward, only: Size_Get
use DKH_Info, only: DKH_Info_Get
use Gateway_Info, only: Gateway_Info_Get
use RICD_Info, only: RICD_Info_Get
use NQ_Info, only: NQ_Info_Get
use External_Centers, only: External_Centers_Get

implicit none

call Symmetry_Info_Get()
call Size_Get()
call DKH_Info_Get()
call Gateway_Info_Get()
call RICD_Info_Get()
call NQ_Info_Get()
call External_Centers_Get()

return

end subroutine Get_Info_Static
