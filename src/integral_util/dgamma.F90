!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

function DGamma_Molcas(arg)

use Definitions, only: wp

implicit none
real(kind=wp) :: DGamma_Molcas
real(kind=wp), intent(in) :: arg
real(kind=wp), external :: gammln

DGamma_Molcas = exp(gammln(arg))

return

end function DGamma_Molcas
