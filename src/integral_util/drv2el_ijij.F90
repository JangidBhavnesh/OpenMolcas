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
! Copyright (C) 1990,1991,1993,1998,2025, Roland Lindh                 *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drv2El_ijij(Pair_Index,nPairs,TMax,nSkal)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals.                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified for k2 loop. August '91                         *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified driver. Jan. '98                                *
!             Modified driver for (ij|ij) integrals. Jan '25           *
!***********************************************************************

use iSD_data, only: iSD
use setup, only: mSkal
use Basis_Info, only: dbsc
use Gateway_Info, only: CutInt
use stdalloc, only: mma_allocate, mma_deallocate
use Integral_interfaces, only: Int_PostProcess, int_wrout
use Int_Options, only: DoIntegrals, DoFock, Disc_Mx
use Basis_Info, only: Shells
use Definitions, only: wp, iwp
use k2_arrays, only: Sew_Scr
use Constants, only: Zero

implicit none
integer(kind=iwp), intent(in):: nPairs, nSkal
integer(kind=iwp), intent(in):: Pair_Index(2,nPairs)
real(kind=wp), intent(inout) :: TMax(nSkal,nSkal)

integer(kind=iwp) :: iCnttp, ijS, iS, jCnttp, jS, iShll, jShll
real(kind=wp) :: A_int, Save_Disc_Mx
character(len=72) :: SLine
real(kind=wp), allocatable :: TInt(:)
integer(kind=iwp), parameter :: nTInt = 1
logical(kind=iwp) :: Save(2)
logical(kind=iwp) :: Deallocate_Sew_Scr
procedure(int_wrout) :: Integral_ijij
procedure(int_wrout), pointer :: Int_postprocess_Save => null()

!                                                                      *
!***********************************************************************
!                                                                      *
! Store the current setting of the integral environment
Save(1)=DoIntegrals
Save(2)=DoFock
Save_Disc_Mx=Disc_Mx
Int_PostProcess_Save => Int_PostProcess

! Setup the integral environment consistent with Drv2el_ijij

DoIntegrals=.True.
DoFock=.False.
Disc_Mx=Zero
Int_PostProcess =>  Integral_ijij

! check if the Seward scratch array is already allocated externally.
! If so do not deallocate on exit
Deallocate_Sew_Scr=.Not.Allocated(Sew_Scr)

SLine = 'Computing 2-electron integrals'
call StatusLine('Seward: ',SLine)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(TInt,nTint,Label='TInt')
!                                                                      *
!***********************************************************************
!                                                                      *

do ijS = 1, nPairs
   iS = Pair_Index(1,ijS)
   jS = Pair_Index(2,ijS)

    iCnttp = iSD(13,iS)
    jCnttp = iSD(13,jS)
    if (dbsc(iCnttp)%fMass /= dbsc(jCnttp)%fMass) Cycle

    iShll = iSD(0,iS)

    ! In case of auxiliary basis sets we want the iS index to point at the dummay shell.
    if (Shells(iShll)%Aux .and. (iS /= mSkal)) cycle

    jShll = iSD(0,jS)

    ! In case the first shell is the dummy auxiliary basis shell make sure that
    ! the second shell, jS, also is a auxiliary basis shell.
    if (Shells(iShll)%Aux .and. (.not. Shells(jShll)%Aux)) cycle

    ! Make sure that the second shell never is the dummy auxiliary basis shell.
    if (Shells(jShll)%Aux .and. (jS == mSkal)) cycle

    A_int = TMax(iS,jS)*TMax(iS,jS)
    if (A_Int < CutInt) Cycle

    call Eval_IJKL(iS,jS,iS,jS,TInt,nTInt)
!   Write (6,*) iS, jS, Tmax(iS,jS), Sqrt(Abs(TInt(1)))
    TMax(iS,jS)=Sqrt(Abs(TInt(1)))
    TMax(jS,iS)=TMax(iS,jS)

end do

call mma_deallocate(TInt)
nullify(Int_PostProcess)

! Restore the status the integral environment as before the call.

DoIntegrals=Save(1)
DoFock=Save(2)
Disc_Mx=Save_Disc_Mx
Int_PostProcess => Int_PostProcess_Save

If (Deallocate_Sew_Scr) Call mma_deallocate(Sew_Scr)

end subroutine Drv2El_ijij
