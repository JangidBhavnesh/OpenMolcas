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

subroutine CHO_MCA_DIAGINT(ISHLA,ISHLB,SCR,LSCR)
!
! Purpose: call Seward to calculate diagonal shell (AB|AB).

use Integral_interfaces, only: Int_PostProcess, int_wrout
#ifdef _DEBUGPRINT_
use Cholesky, only: LuPri
use Gateway_Info, only: CutInt, ThrInt
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ISHLA, ISHLB, LSCR
real(kind=wp), intent(out) :: SCR(LSCR)
procedure(int_wrout) :: Integral_Wrout_Cho_Diag
#ifdef _DEBUGPRINT_
real(kind=wp) :: CUTINT1, CUTINT2, THRINT1, THRINT2
character(len=*), parameter :: SECNAM = 'CHO_MCA_DIAGINT'
#endif

Int_PostProcess => Integral_Wrout_Cho_Diag

#ifdef _DEBUGPRINT_
CUTINT1 = CutInt
THRINT1 = ThrInt
#endif

call EVAL_IJKL(ISHLA,ISHLB,ISHLA,ISHLB,SCR,LSCR)

nullify(Int_PostProcess)
#ifdef _DEBUGPRINT_
CUTINT2 = CutInt
THRINT2 = ThrInt
if ((CUTINT2 /= CUTINT1) .or. (THRINT2 /= THRINT1)) then
  write(LUPRI,*) SECNAM,': CutInt before Eval_Ints_: ',CUTINT1
  write(LUPRI,*) SECNAM,': CutInt after  Eval_Ints_: ',CUTINT2
  write(LUPRI,*) SECNAM,': ThrInt before Eval_Ints_: ',THRINT1
  write(LUPRI,*) SECNAM,': ThrInt after  Eval_Ints_: ',THRINT2
  call CHO_QUIT('Integral prescreening error detected in '//SECNAM,102)
end if
#endif

end subroutine CHO_MCA_DIAGINT
