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
! Copyright (C) 2023, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine get_charge(charge)
  use definitions, only: wp, iwp
  use mcpdft_output, only: lf
  use OneDat, only: sNoOri
  implicit none

  integer(kind=iwp), intent(out) :: charge

#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
  integer(kind=iwp) :: RC, Opt, Comp, SyLbl, isym
  character(len=8) :: Label
  real(kind=wp), dimension(nTot1+4) :: overlap

  RC = -1
  Opt = ibset(0, sNoOri)
  Comp = 1
  SyLbl = 1
  Label = 'Mltpl  0'
  call rdone(rc, opt, label, comp, overlap, sylbl)
  if(RC /= 0) then
    write(lf, *) 'RC from call RdOne not 0'
    write(lf, *) 'Label = ', label
    write(lf, *) 'RC = ', RC
    call Abend
  end if
  tot_nuc_charge = overlap(nTot1+4)

  tot_el_charge = 0.0D0
  do iSym = 1, nSym
    tot_el_charge = tot_el_charge - 2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
  end do
  tot_el_charge = tot_el_charge - DBLE(nActEl)

  charge = INT(tot_nuc_charge + tot_el_charge)
end subroutine get_charge
