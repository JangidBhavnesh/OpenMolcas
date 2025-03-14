************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine DetCtl()
      use Arrays, only: pINT1, pINT2
      use stdalloc, only: mma_allocate
      use MCLR_Data, only: iST,i12
      use MCLR_Data, only: IPRSTR,IPRORB,IPRCIX
      use MCLR_Data, only: MS2,idc,PSSIGN
      use MCLR_Data, only: FnCSF2SD, LuCSF2SD
      use MCLR_Data, only: NOCSF,IDENMT,NOPART,IDIAG,INTIMP,INCORE,
     &                     ICISTR
      use input_mclr, only: nSym,PntGrp,nIrrep,nsMOB,iSpin,
     &                      nHole1,nActEl,nElec3,
     &                      nRs1,nRs2,nRs3,State_Sym
*
      Implicit None

      Integer iTmp, nTRas1,nTRas2,nTRas3,iSym,iDum,MNRS10,MXR4TP,MXRS30

      Call mma_Allocate(pINT1,nSym,Label='pInt1')
      pInt1(:)=0
      Call mma_Allocate(pINT2,nSym**3,Label='pInt2')
      pInt2(:)=0

      Pntgrp=1
      NOCSF  = 0
      idenmt=0
      nopart=0
      nIrrep=nSym
      nsmob=nSym
      mxr4tp=0
      idiag=1
      intimp=5
      incore=1
      icistr=1
      ist=1
      i12=2
      MS2=iSpin-1
      If (ms2.ne.0) Then
        idc=1
        pssign=0.0d0
      else
        itmp=(ispin-1)/2
        pssign=(-1.0d0)**itmp
        idc=2
      end if

      ntRas1=0
      ntRas2=0
      ntRas3=0
      Do iSym=1,nSym
         ntRas1=ntRas1+nRs1(iSym)
         ntRas2=ntRas2+nRs2(iSym)
         ntRas3=ntRas3+nRs3(iSym)
      End Do
      MNRS10 = Max(0,2*ntRas1-nHole1)
      MXRS30 = Max(0,Min(2*ntRas3,nElec3))
*. Initialize print flags
      IPRSTR =  0
      IPRORB =  0
      IPRCIX =  0
*. From shells to orbitals
      CALL ORBINF_MCLR(nSym,nSym,nRs1,nRs2,nRs3,mxr4tp,IPRORB) ! OK
*. Number of string types
      CALL STRTYP(ms2,nActEl,MNRS10,MXRS30,IPRSTR)   ! looks allright
*. Symmetry information
      CALL SYMINF_MCLR(nSym,IPRORB) ! looks allright
*. Internal string information
      CALL STRINF(IPRSTR)     ! looks allright, no!
*. Internal subspaces
      CALL ICISPC(MNRS10,MXRS30,IPRCIX)    ! looks allright
      CALL ICISPS(IPRCIX)        ! looks allright
*. CSF information
      CALL DANAME(LUCSF2SD,FNCSF2SD)
      CALL CSFINF(State_sym,iSpin,idum,1,IPRCIX,nsym)
*
      End Subroutine DetCtl
