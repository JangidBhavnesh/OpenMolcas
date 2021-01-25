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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************

subroutine RotateGeoms(Q)

use Quater_globals, only: debug, ngeoms, list
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Q(0:3)
integer(kind=iwp) :: igeom

do igeom=3,ngeoms+2
  if (debug) then
    write(u6,*) 'Before rotation'
    call PrintGeom(u6,list(igeom)%nat,list(igeom)%title,list(igeom)%geo,list(igeom)%geolbl)
  end if
  call RotateGeom(Q,list(igeom)%nat,list(igeom)%geo)
  if (debug) then
    write(u6,*) 'After rotation'
    call PrintGeom(u6,list(igeom)%nat,list(igeom)%title,list(igeom)%geo,list(igeom)%geolbl)
  end if
end do

return

end subroutine RotateGeoms
