************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2016,2017, Roland Lindh                                *
************************************************************************
      SubRoutine GrdClc(FstItr,iOpt)
      use SCF_Arrays
      use InfSCF
      Implicit Real*8 (a-h,o-z)
      Logical FstItr
      Integer   iOpt
#include "real.fh"
*
      nD = iUHF + 1
*                                                                      *
************************************************************************
*                                                                      *
*     iOpt    : SCF optimization scheme
*
*     Select case(iOpt)

*     Case (2,3,4)
         Call GrdClc_(FstItr,Dens,TwoHam,Vxc,nBT,nDens,nD,OneHam,
     &                CMO   ,nBB,Ovrlp,CMO)
*     Case (1)
*        Call GrdClc_(FstItr,Dens,TwoHam,Vxc,nBT,nDens,nD,OneHam,
*    &                Lowdin,nBB,Ovrlp,CMO)
*     Case Default
*        Write (6,*) 'GrdClc: Illegal iOpt Value:',iOpt
*        Call abend()
*     End Select
*
      Return
      End
      SubRoutine GrdClc_(FstItr,Dens,TwoHam,Vxc,mBT,mDens,nD,OneHam,
     &                   OCMO,mBB,Ovrlp,CMO)
************************************************************************
*                                                                      *
*     purpose: Compute gradients and write on disk.                    *
*                                                                      *
*                                                                      *
*     input:                                                           *
*       FstItr  : variable telling what gradients compute: .true. -    *
*                 all gradients, .False. - last gradient               *
*                                                                      *
*     called from: Wfctl_scf                                           *
*                                                                      *
*     calls to:         EGrad                                          *
*               uses SubRoutines and Functions from Module lnklst.f    *
*               -linked list implementation to store series of vectors *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Use Interfaces_SCF, Only: vOO2OV
      Use InfSO
      Use InfSCF
      use LnkLst, only: LLGrad
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "file.fh"
*
      Real*8 Dens(mBT,nD,mDens), TwoHam(mBT,nD,mDens), CMO(mBB,nD),
     &       OneHam(mBT), OCMO(mBB,nD), Ovrlp(mBT), Vxc(mBT,nD,mDens)
      Real*8, Dimension(:,:), Allocatable:: GrdOO,AuxD,AuxT,AuxV
      Real*8, Allocatable:: GrdOV(:)
      Logical FstItr
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*#define _DEBUGPRINT_
*
*--- Allocate memory for gradients and gradient contributions
      Call mma_allocate(GrdOO,nOO,nD,Label='GrdOO')
      Call mma_allocate(GrdOV,mOV,Label='GrdOV')

*--- Allocate memory for auxiliary matrices
      Call mma_allocate(AuxD,nBT,nD,Label='AuxD')
      Call mma_allocate(AuxT,nBT,nD,Label='AuxT')
      Call mma_allocate(AuxV,nBT,nD,Label='AuxV')

*--- Find the beginning of the loop
      If (FstItr) Then
         LpStrt = 1
         kOptim_=iter
         FstItr=.False.
      Else
         LpStrt = kOptim
         kOptim_= kOptim
      End If

*--- Compute all gradients / last gradient
*
      iter_d=iter-iter0
      Do iOpt = LpStrt, kOptim_
         iDT = iter_d - kOptim_ + iOpt
*
         GrdOV(:)=Zero
*
         jDT=MapDns(iDT)
         If (jDT.lt.0) Then
           Call RWDTG(-jDT,AuxD,nBT*nD,'R','DENS  ',iDisk,SIZE(iDisk,1))
           Call RWDTG(-jDT,AuxT,nBT*nD,'R','TWOHAM',iDisk,SIZE(iDisk,1))
           Call RWDTG(-jDT,AuxV,nBT*nD,'R','dVxcdR',iDisk,SIZE(iDisk,1))
*
            Call EGrad(OneHam,AuxT,AuxV,Ovrlp,AuxD,nBT,OCMO,nBO,
     &                 GrdOO,nOO,nD,CMO)
*
         Else
*
            Call EGrad(OneHam,TwoHam(1,1,jDT),Vxc(1,1,jDT),Ovrlp,
     &                 Dens(1,1,jDT),nBT,OCMO,nBO,GrdOO,nOO,nD,CMO)
*
         End If
*
         Call vOO2OV(GrdOO,nOO,GrdOV,mOV,nD,nOV,kOV)
*
*------- Write Gradient to linked list
*
         Call PutVec(GrdOV,mOV,iDT+iter0,'OVWR',LLGrad)
*
#ifdef _DEBUGPRINT_
         Write (6,*) 'GrdClc: Put Gradient iteration:',iDT+iter0
         Write (6,*) 'iOpt=',iOpt
         Call NrmClc(GrdOO,nOO*nD,'GrdClc','GrdOO')
         Call NrmClc(GrdOV,mOV,'GrdClc','GrdOV')
         Call RecPrt('GrdClc: g(i)',' ',GrdOV,1,mOV)
#endif
      End Do
*
*     Deallocate memory
*
      Call mma_deallocate(AuxD)
      Call mma_deallocate(AuxT)
      Call mma_deallocate(AuxV)
      Call mma_deallocate(GrdOV)
      Call mma_deallocate(GrdOO)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
