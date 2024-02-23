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
      Subroutine CXInit(SGS,CIS,EXS)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Struct, only:  SGStruct, CIStruct, EXStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Type (EXStruct) EXS

      Integer, Allocatable:: IVR(:), ISgm(:), NRL(:), iLndw(:), Scr(:)
      Real*8,  Allocatable::         VSgm(:)
      Real*8,  Allocatable:: VTabTmp(:), Val(:)

      nSym   =SGS%nSym
      nLev   =SGS%nLev
      nVert  =SGS%nVert
      MidLev =SGS%MidLev
      MVSta  =SGS%MVSta
      MVEnd  =SGS%MVEnd

C Calculate segment values, and MVL and MVR tables:
      nMidV=1+MVEnd-MVSta
C nIpWlk: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.
      nIpWlk=1+(MidLev-1)/15
      nIpWlk=MAX(nIpWlk,1+(nLev-MidLev-1)/15)
      Call mma_allocate(IVR,2*nVert,Label='IVR')
      Call mma_allocate(EXS%MVR,2*nMidV)
      Call mma_allocate(EXS%MVL,2*nMidV)
      nSgmnt=26*nVert
      Call mma_allocate(ISgm,nSgmnt,Label='ISgm')
      Call mma_allocate(VSgm,nSgmnt,Label='VSgm')
      Call MkSeg(SGS,nLev,nVert,nMidv,SGS%DRT,SGS%Down,SGS%LTV,
     &           IVR,EXS%MVL,EXS%MVR,ISgm,VSgm)
      CIS%nMidV   =nMidV
      CIS%nIpWlk  = nIpWlk

C Various offset tables:
      nNOW=2*nMidV*nSym
      Call mma_allocate(CIS%NOW,nNOW,Label='CIS%NOW')
      Call mma_allocate(CIS%IOW,nNOW,Label='CIS%IOW')
      MxEO=(nLev*(nLev+5))/2
      nNOCP=MxEO*nMidV*nSym
      nIOCP=nNOCP
      nNRL=(1+MxEO)*nVert*nSym
      Call mma_allocate(EXS%NOCP,nNOCP,Label='EXS%NOCP')
      Call mma_allocate(EXS%IOCP,nIOCP,Label='EXS%IOCP')
      Call mma_allocate(CIS%NCSF,nSym,Label='CIS%NCSF')
      Call mma_allocate(NRL,nNRL,Label='NRL')
      nNOCSF=nMidV*(nSym**2)
      nIOCSF=nNOCSF

      Call mma_allocate(CIS%NOCSF,nNOCSF,Label='CIS%NOCSF')
      Call mma_allocate(CIS%IOCSF,nIOCSF,Label='CIS%IOCSF')
      EXS%MxEO =MxEO
      Call NrCoup(SGS,CIS,EXS,
     &         nVert,nMidV,MxEO,SGS%ISm,SGS%DRT,
     &         ISgm,CIS%NOW,CIS%IOW,EXS%NOCP,
     &         EXS%IOCP,CIS%NOCSF,CIS%IOCSF,
     &         CIS%NCSF,NRL,EXS%MVL,EXS%MVR)
      Call mma_deallocate(NRL)
C Computed in NrCoup:
      nWalk=CIS%nWalk
      nICoup=EXS%nICoup

      nICase=nWalk*nIpWlk
      Call mma_allocate(CIS%ICase,nICase,Label='CIS%ICase')
      nnICoup=3*nICoup
      Call mma_allocate(EXS%ICoup,nnICoup,Label='EXS%ICoup')
      nVMax=5000
      Call mma_allocate(VTabTmp,nVMax,Label='VTabTmp')
      nILNDW=nWalk
      Call mma_allocate(iLndw,niLndw,Label='iLndw')
      nScr=7*(nLev+1)
      Call mma_allocate(Scr,nScr,Label='Scr')
      Call mma_allocate(Val,nLev+1,Label='Val')
      EXS%nVTab =nVMax
      nVTab=nVMax
      Call MkCoup(nLev,SGS%Ism,nVert,MidLev,nMidV,MVSta,MVEnd,
     &            MxEO,nICoup,nWalk,nICase,nVTab,
     &            IVR,SGS%MAW,ISGM,VSGM,CIS%NOW,CIS%IOW,EXS%NOCP,
     &            EXS%IOCP,ILNDW,CIS%ICase,EXS%ICOUP,VTabTmp,SCR,VAL)

C nVTab has now been updated to the true size. Allocate final array:
      Call mma_allocate(EXS%Vtab,nVTab,Label='EXS%VTab')
      EXS%nVTab=nVTab
      EXS%VTab(1:nVTab)=VTabTmp(1:nVTab)
      Call mma_deallocate(VTabTmp)

      Call mma_deallocate(iLndw)
      Call mma_deallocate(Scr)
      Call mma_deallocate(Val)
      Call mma_deallocate(ISgm)
      Call mma_deallocate(VSgm)
      Call mma_deallocate(IVR)

      end Subroutine CXInit
