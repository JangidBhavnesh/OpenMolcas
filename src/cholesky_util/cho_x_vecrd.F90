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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  Cho_X_VecRd
!
!> @brief
!>   Read as many Cholesky vectors as possible in the range \p jVec1, \p jVec1+1, \p jVec1+2, ..., \p jVec2.
!> @author Thomas Bondo Pedersen
!>
!> @details
!> The vectors are returned in their native storage (reduced set).
!> Starting with vector \p jVec1, this routine will read as many
!> vectors as possible, although at most (\p jVec2-\p jVec1+1) vectors
!> are read.
!>
!> On entry as well as exit, \p iRedC identifies the
!> reduced set stored at location ``3`` (``-1`` if none or unknown).
!> On exit, \p jNum is the number
!> of vectors actually read and \p mUsed is the memory
!> (in ``real`` words) actually used.
!>
!> @note
!> The Cholesky module must be initialized.
!>
!> @param[out]    Scr   contains the vectors on exit
!> @param[in]     lScr  dimension of \p Scr
!> @param[in]     jVec1 first vector to read
!> @param[in]     jVec2 last vector allowed to read
!> @param[in]     iSym  vector symmetry
!> @param[out]    jNum  number of vectors actually read
!> @param[in,out] iRedC reduced set stored at location ``3`` on entry as well as exit
!> @param[out]    mUsed amount of memory actually used (in ``real`` words)
!***********************************************************************

subroutine Cho_X_VecRd(Scr,lScr,jVec1,jVec2,iSym,jNum,iRedC,mUsed)

use Cholesky, only: NumCho
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lScr, jVec1, jVec2, iSym
real(kind=wp), intent(out) :: Scr(lScr)
integer(kind=iwp), intent(out) :: jNum, mUsed
integer(kind=iwp), intent(inout) :: iRedC
integer(kind=iwp) :: l_jVec2

if ((iSym < 1) .or. (iSym > 8)) then
  jNum = 0
  mUsed = 0
else
  l_jVec2 = min(NumCho(iSym),jVec2)
  call Cho_VecRd(Scr,lScr,jVec1,l_jVec2,iSym,jNum,iRedC,mUsed)
end if

end subroutine Cho_X_VecRd
