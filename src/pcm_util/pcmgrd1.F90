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
! Copyright (C) 1995,2001, Roland Lindh                                *
!***********************************************************************

subroutine PCMgrd1( &
#                  define _CALLING_
#                  include "grd_interface.fh"
                  )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, May '95                                 *
!                                                                      *
!             Modified to PCM gradients September 2001, Lund, by       *
!             R. Lindh.                                                *
!***********************************************************************

use PCM_arrays, only: PCMTess
use Center_Info, only: dc
use Index_Functions, only: nTri_Elem1
use Symmetry_Info, only: ChOper
use Rys_interfaces, only: cff2d_kernel, modu2_kernel, tval1_kernel
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
integer(kind=iwp) :: i, iAlpha, iAnga(4), iBeta, iCar, iDAO, iDCRT(0:7), ii, ipA, ipAOff, ipB, ipBOff, ipDAO, iPrint, iRout, &
                     iStb(0:7), iTs, iuvwx(4), iZeta, JndGrd(3,4), lDCRT, LmbdT, lOp(4), mGrad, mRys, nArray, nDAO, nDCRT, nDiff, &
                     nip, nStb, nT
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), Fact, Q, TC(3)
logical(kind=iwp) :: NoLoop, JfGrad(3,4)
procedure(cff2d_kernel) :: XCff2D
procedure(modu2_kernel) :: Fake
procedure(tval1_kernel) :: TNAI1
integer(kind=iwp), external :: NrOpr
#include "print.fh"

#include "macros.fh"
unused_var(rFinal)
unused_var(nHer)
unused_var(Ccoor(1))
unused_var(nComp)

iRout = 151
iPrint = nPrint(iRout)

! Modify the density matrix with the prefactor

nDAO = nTri_Elem1(la)*nTri_Elem1(lb)
do iDAO=1,nDAO
  do iZeta=1,nZeta
    Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
    DAO(iZeta,iDAO) = Fact*DAO(iZeta,iDAO)
  end do
end do
if (iPrint >= 99) call RecPrt('DAO',' ',DAO,nZeta,nDAO)

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
ipDAO = nip
nip = nip+nAlpha*nBeta*nTri_Elem1(la)*nTri_Elem1(lb)*nTri_Elem1(nOrdOp)
if (nip-1 > nZeta*nArr) then
  write(u6,*) 'nip-1 > nZeta*nArr'
  call AbEnd()
end if
nArray = nZeta*nArr-nip+1

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = nOrdOp
iAnga(4) = 0
Coori(:,1) = A
Coori(:,2) = RB

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  CoorAC(:,1) = A
else
  CoorAC(:,1) = RB
end if
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
lOp(1) = kOp(1)
lOp(2) = kOp(2)

ipAOff = ipA
do iBeta=1,nBeta
  call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
  ipAOff = ipAOff+nAlpha
end do

ipBOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
  ipBOff = ipBOff+1
end do

! Loop over the tiles

! pcm_solvent remove the loop
!do iTs=1,nTs
do iTs=1,1
! pcm_solvent end
  ! pcm_solvent put "charge" to 1
  !Q = PCM_SQ(1,iTs)+PCM_SQ(2,iTs)
  Q = One
  ! pcm_solvent end
  NoLoop = Q == Zero
  if (NoLoop) cycle
  ! Pick up the tile coordinates
  C(1:3) = PCMTess(1:3,iTs)

  if (iPrint >= 99) call RecPrt('C',' ',C,1,3)

  ! Generate stabilizer of C

  nStb = 1
  iStb(0) = 0

  ! Find the DCR for M and S

  call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
  Fact = -real(nStabM,kind=wp)/real(LmbdT,kind=wp)

  if (iPrint >= 99) then
    write(u6,*) ' Q=',Q
    write(u6,*) ' Fact=',Fact
    call RecPrt('DAO*Fact*Q',' ',Array(ipDAO),nZeta*nDAO,nTri_Elem1(nOrdOp))
    write(u6,*) ' m      =',nStabM
    write(u6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
    write(u6,*) ' s      =',nStb
    write(u6,'(9A)') '(S)=',(ChOper(iStb(ii)),ii=0,nStb-1)
    write(u6,*) ' LambdaT=',LmbdT
    write(u6,*) ' t      =',nDCRT
    write(u6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
  end if
  iuvwx(3) = nStb
  iuvwx(4) = nStb
  JndGrd(:,1:2) = IndGrd
  JfGrad(:,1:2) = IfGrad

  ! No derivatives with respect to the third or fourth center.
  ! The positions of the points in the external field are frozen.

  JndGrd(:,3:4) = 0
  JfGrad(:,3:4) = .false.
  mGrad = 0
  do iCar=1,3
    do i=1,2
      if (JfGrad(iCar,i)) mGrad = mGrad+1
    end do
  end do
  if (iPrint >= 99) write(u6,*) ' mGrad=',mGrad
  if (mGrad == 0) cycle

  do lDCRT=0,nDCRT-1
    lOp(3) = NrOpr(iDCRT(lDCRT))
    lOp(4) = lOp(3)
    call OA(iDCRT(lDCRT),C,TC)
    CoorAC(:,2) = TC
    Coori(:,3) = TC
    Coori(:,4) = TC

    Array(ipDAO:ipDAO+nZeta*nDAO-1) = Fact*Q*pack(DAO,.true.)

    ! Compute integrals with the Rys quadrature.

    nT = nZeta
    nDiff = 1
    mRys = (la+lb+2+nDiff+nOrdOp)/2
    call Rysg1(iAnga,mRys,nT,Array(ipA),Array(ipB),[One],[One], &
               Zeta,ZInv,nZeta,[One],[One],1, &
               P,nZeta,TC,1,Coori,Coori,CoorAC, &
               Array(nip),nArray,TNAI1,Fake,XCff2D,Array(ipDAO),nDAO*nTri_Elem1(nOrdOp),Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)

    !call RecPrt(' In PCMgrd1:Grad',' ',Grad,nGrad,1)
  end do  ! End loop over DCRs

end do     ! End loop over centers in the external field

return

end subroutine PCMgrd1
