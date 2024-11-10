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
* Copyright (C) 1989, Bjorn O. Roos                                    *
*               1989, Per Ake Malmqvist                                *
*               1991,1993,1996, Markus P. Fuelscher                    *
************************************************************************
      Subroutine Fmat_m(CMO,D1A,FI,FA)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Update the Fock matrix for the active orbitals and transform     *
*     it to MO basis as well as the matrix FI (Fock matrix) for        *
*     frozen and inactive orbitals).                                   *
*                                                                      *
*     calling arguments:                                               *
*     CMO     : array of real*8                                        *
*               MO-coefficients                                        *
*     D1A     : array of real*8                                        *
*               active one body density matrix in AO-basis             *
*     FI      : array of real*8                                        *
*               inactive Fock matrix                                   *
*     FA      : array of real*8                                        *
*               active Fock matrix                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     B.O. Roos and P.Aa. Malmqvist                                    *
*     University of Lund, Sweden, 1989                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     - updated for MOLCAS version 2                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1991               *
*     - updated for MOLCAS version 3                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1993               *
*     - updated for integral direct and reaction field calculations    *
*       M.P. Fuelscher, University of Lund, Sweden, 1996               *
*                                                                      *
************************************************************************

      use printlevel, only: debug
      use mcpdft_output, only: lf, iPrLoc
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_global, only: ECAS, EMY, VIA

      Implicit None

      Real*8 CMO(*) , D1A(*) , FI(*) , FA(*)

#include "rasdim.fh"
#include "general.fh"
      Character(LEN=16), Parameter:: ROUTINE='FMAT    '


      Real*8, allocatable:: Tmp1(:), Tmp2(:)
      Integer iBas, iFro, iOff, iOff1, iOff2, iOff3, iOrb, iPrLev, iSym
      Real*8, External:: DDot_
      real*8 vaa

C Local print level (if any)
      IPRLEV=IPRLOC(4)

!************************************************************
! Here we should start the real work!
!************************************************************
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
!     create FA in AO basis
      Call mma_allocate(Tmp1,nTot1,Label='Tmp1')
      Call Fold(nSym,nBas,D1A,Tmp1)

      ! Active-Active contribution to ECAS
      VAA = 0.5D0*ddot_(nTot1,FA,1,Tmp1,1)

!     Inactive-active contribution to ECAS
      VIA=dDot_(nTot1,FI,1,Tmp1,1)
      ECAS = EMY + VIA + VAA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,'(A,ES20.10)') ' Total core energy:            ',EMY
        Write(LF,'(A,ES20.10)') ' inactive-active interaction:  ',VIA
        Write(LF,'(A,ES20.10)') ' active-active interaction:  ',VAA
        Write(LF,'(A,ES20.10)') ' CAS energy (core+interaction):',ECAS
      End If
      Call mma_deallocate(Tmp1)

!     transform FI from AO to MO basis
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        If (iBas==0) Cycle
        iOrb = nOrb(iSym)
        If (iOrb==0) Cycle
        iFro = nFro(iSym)
        Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
        Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
        Call Square(FI(iOff1),Tmp1,1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Tmp1,iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Tmp2,iBas)
        Call DGEMM_Tri('T','N',iOrb,iOrb,iBas,
     &                 1.0D0,Tmp2,iBas,
     &                       CMO(iOff2+(iFro*iBas)),iBas,
     &                 0.0D0,FI(iOff3),iOrb)
        Call mma_deallocate(Tmp2)
        Call mma_deallocate(Tmp1)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

! transform FA from AO to MO basis
      call ao2mo(cmo,fa,fa)

!     print FI and FA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI in MO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FI(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
        Write(LF,*)
        Write(LF,*) ' FA in MO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If

      End Subroutine Fmat_m
