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
* Copyright (C) 1997, Luis Serrano-Andres                              *
************************************************************************
      SUBROUTINE SUPSCH(SMAT,CMOO,CMON)
C
C     Program RASSCF
C
C     Objective: To check the order of the input of natural orbitals
C                to obtain the right labels for the Supersymmetry matrix.
C
C     Called from ortho, neworb, fckpt2, and natorb.
C
C     Luis Serrano-Andres
C     University of Lund, Sweden, 1997
C     **** Molcas-4 *** Release 97 04 01 **********
C
      use stdalloc, only: mma_allocate, mma_deallocate

      IMPLICIT None
#include "rasdim.fh"
#include "general.fh"
      Real*8 CMOO(*),CMON(*),SMAT(*)

      Real*8, Allocatable:: Temp1(:), Temp2(:)
      Integer, Allocatable:: IxSym2(:)
      Integer :: iSym, nOrb_Tot, nOrbMx
*
*
      nOrbMX=0
      nOrb_tot=0
      Do iSym=1,nSym
         nOrbMX=Max(nOrbMX,nBas(iSym))
         nOrb_tot=nOrb_tot+nBas(iSym)
      End Do
*
      Call mma_allocate(Temp1,nOrbMX*nOrbMX,Label='Temp1')
      Call mma_allocate(Temp2,nOrbMX*nOrbMX,Label='Temp2')
      Call mma_allocate(IxSym2,nOrb_tot,Label='IxSym2')
*
      Call SUPSCH_(SMAT,CMOO,CMON,Temp1,Temp2,nOrbMX,IxSym2,nOrb_tot)
*
      Call mma_deallocate(IxSym2)
      Call mma_deallocate(Temp2)
      Call mma_deallocate(Temp1)
*
*
      End SUBROUTINE SUPSCH
      SUBROUTINE SUPSCH_(SMAT,CMOO,CMON,Temp1,Temp2,nOrbMX,IxSym2,
     &                   nOrb_tot)
C
C     Program RASSCF
C
C     Objective: To check the order of the input of natural orbitals
C                to obtain the right labels for the Supersymmetry matrix.
C
C     Called from ortho, neworb, fckpt2, and natorb.
C
C     Luis Serrano-Andres
C     University of Lund, Sweden, 1997
C     **** Molcas-4 *** Release 97 04 01 **********
C
      use OneDat, only: sNoNuc, sNoOri
      use rasscf_global, only: FDIAG, iSupSM, Iter, ixsym
      use output_ras, only: LF,DEBUG,IPRLOC

      IMPLICIT None
#include "rasdim.fh"
#include "warnings.h"
#include "general.fh"
      Integer nOrbMX, nOrb_Tot
      Real*8 CMOO(*),CMON(*),SMAT(*)
      Real*8 Temp1(nOrbMX*nOrbMX),Temp2(nOrbMX*nOrbMX)
      INTEGER IxSym2(nOrb_tot)

      Character(LEN=16), Parameter ::  ROUTINE='SUPSCH_ '
      Integer pSij
      Real*8 DUM(1)
      character(len=8) :: Label
      Integer :: i_Component, i_Opt, i_RC, i_SymLbl, iGroup, iLabel,
     &           iOrb, iOrder, iPrLev, iSafe, jOrb, kCof, kGroup, kOrb,
     &           nBs, nnOrb, nOGr1, nOGr2, iSym
      Real*8 :: OldOvlp, Ovlp1, Ovlp2, xOvlp
C
C Local print level (if any)
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
*
*     Read overlap matrix SMAT:
*
      i_Rc=0
      i_Opt=ibset(ibset(0,sNoOri),sNoNuc)
      i_Component=1
      i_SymLbl=1
      Label='Mltpl  0'
      Call RdOne(i_Rc,i_Opt,Label,i_Component,Smat,i_SymLbl)
      If ( i_Rc.ne.0 ) Then
        Write(LF,*)
        Write(LF,*)' ********************* ERROR **********************'
        Write(LF,*)' SUPSCH: Failed to read overlap from ONEINT.       '
        Write(LF,*)' RASSCF is using overlaps to compare old and new   '
        Write(LF,*)' orbitals, but could not read overlaps from ONEINT.'
        Write(LF,*)' Something is wrong with the file, or possibly with'
        Write(LF,*)' the program. Please check.                        '
        Write(LF,*)' **************************************************'
        Call Quit(_RC_IO_ERROR_READ_)
      End If
C
C
C     Check order of the orbitals for supersymmetry option
C
      If( ( ISupsm.eq.1 ).and.( Iter.ge.1 ) ) then
C
       If(IPRLEV.GE.DEBUG) then
         CALL PRIMO_RASSCF('Testing old orb for supersymmetry',
     &               FDIAG,DUM,CMOO)
         CALL PRIMO_RASSCF('Testing new orb for supersymmetry',
     &               FDIAG,DUM,CMON)
       End if
C
       kOrb=0
       kCof=0
       pSij=0
       Do iOrb=1,nOrb_tot
         IxSym2(iOrb)=0
       End do
       Do iSym=1,nSym
         iSafe=0
         nBs=nBas(iSym)
         If(nBs.le.0) Goto 1966
C
C        Computing orbital overlaping Sum(p,q) C1kp* C2lq Spq
C
         Call Square(SMat(pSij+1),Temp1,1,nBs,nBs)
         Call DGEMM_('N','N',
     &               nBs,nBs,nBs,
     &               1.0d0,Temp1,nBs,
     &               CMON(kCof+1),nBs,
     &               0.0d0,Temp2,nBs)
         Call DGEMM_('T','N',
     &               nBs,nBs,nBs,
     &               1.0d0,CMOO(kCof+1),nBs,
     &               Temp2,nBs,
     &               0.0d0,Temp1,nBs)
C
C        Checking the maximum overlap between the orbitals
C        and building the new supersymmetry matrix
C
         nnOrb=1
         Do iOrb=1,nBs
           iLabel=IxSym(kOrb+iOrb)
           iOrder=1
           xOvlp=0.0d0
           Do jOrb=0,nBs-2
             Ovlp1=Abs(Temp1(nnOrb+(jOrb*nBs)))
             Ovlp2=Abs(Temp1(nnOrb+((jOrb+1)*nBs)))
             OldOvlp=xOvlp
             xOvlp=Max(Ovlp1,Ovlp2,xOvlp)
             If (jOrb.eq.0) then
               If(Ovlp1.lt.Ovlp2) iOrder=2
             Else
               If(xOvlp.ne.OldOvlp) iOrder=jOrb+2
             End if
           End do
           IxSym2(kOrb+iOrder)=iLabel
           nnOrb=nnOrb+1
         End do
C
C        Number of orbital groups on the symmetry
C
         kGroup=0
         Do iOrb=1,nBs-1
           kGroup=Max(IxSym(kOrb+iOrb),IxSym(kOrb+iOrb+1),kGroup)
         End do
         Do iGroup=0,kGroup
           nOGr1=0
           nOGr2=0
           Do iOrb=1,nBs
             if(IxSym(kOrb+iOrb).eq.iGroup) nOGr1=nOGr1+1
             if(IxSym2(kOrb+iOrb).eq.iGroup) nOGr2=nOGr2+1
           End do
C
C        Checking if we have the same number of orbitals
C        per group
C
           If((nOGr2.ne.nOGr1).and.(iGroup.ne.0)) then
             iSafe=1
             Call WarningMessage(1,'Supersymmetry may have failed.')
             Write(LF,*)' Check orbital order or try cleaning orbitals.'
           End if
         End do
C
C        New matrix replaces the old one
C
         If(iSafe.eq.0) then
           Do iOrb=1,nBs
             IxSym(kOrb+iOrb)=IxSym2(kOrb+iOrb)
           End do
         End if
         kCof=kCof+(nBs*nBs)
         pSij=pSij+(nBs*nBs+nBs)/2
         kOrb=kOrb+nBs
1966     Continue
       End do
      End if
C
      END SUBROUTINE SUPSCH_
