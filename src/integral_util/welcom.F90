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

module welcom

use define_af, only: iTabMx
use Constants, only: Half, Quart
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: kMax = iTabMx+6, k2 = int(kmax*Half)+1, k4 = int(kmax*Quart)+1
integer(kind=iwp) :: ipot3(0:kmax+1)
real(kind=wp) :: anorm(0:kmax,0:k2,0:k4), binom(-1:kmax,-1:kmax), fac(0:kmax), fiint(0:kmax,0:kmax), tetint(0:kmax,0:k2)

public :: anorm, binom, fac, fiint, ipot3, kMax, tetint

end module welcom
