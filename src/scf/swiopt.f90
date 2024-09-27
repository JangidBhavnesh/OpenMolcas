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
! Copyright (C) 1996, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************
      SubRoutine SwiOpt(AllCnt,OneHam,Ovrlp,mBT,CMO,mBB,nD)
!***********************************************************************
!                                                                      *
!     purpose: Switch from HF AO to MO optimization                    *
!                                                                      *
!***********************************************************************
      use OneDat, only: sNoNuc, sNoOri
      Use InfSO, only: DltNth
      use InfSCF, only: DThr, EThr, FThr, nBO, nBT, nIterP, PotNuc
      use Files, only: FnDel, FnDGd, FnDSt, FnGrd, FnOSt, FnTSt, Fnx, Fny, &
                       LuDel, LuDGd, LuDSt, LuGrd, LuOSt, LuTSt, Lux, Luy
      use NDDO, only: twoel_NDDO
      use Constants, only: One
      Implicit None
!
!     declaration of subroutine parameter...
      Logical AllCnt
      Integer mBT, mBB, nD
      Real*8 OneHam(mBT), Ovrlp(mBT), CMO(mBB,nD)

!     declaration of some local variables...
      Real*8 EThr_o,FThr_o,DThr_o,DNTh_o, ThrInt_o
      Save EThr_o,FThr_o,DThr_o,DNTh_o, ThrInt_o
      Character(LEN=8) Label
      Integer iComp, iD, iOpt, iRC, lOper
      Real*8, External:: Get_ThrInt
!
      If (AllCnt .AND. twoel_NDDO) Then
        nIterP=1
!       read full overlap matrix from ONEINT file
        Label='Mltpl  0'
        iOpt = ibset(ibset(0,sNoOri),sNoNuc)
        iRC = -1
        iComp = 1
        Call RdOne(iRC,iOpt,Label,iComp,Ovrlp,lOper)
        If (iRC.ne.0) GoTo 9999
!       read full one-electron Hamiltonian from ONEINT file
        Label='OneHam  '
        iOpt=ibset(ibset(0,sNoOri),sNoNuc)
        iRC = -1
        Call RdOne(iRC,iOpt,Label,iComp,OneHam,lOper)
        If (iRc.ne.0) GoTo 9999
!       Call Get_PotNuc(PotNuc)
!       Call Get_dScalar('PotNuc',PotNuc)
        Call Peek_dScalar('PotNuc',PotNuc)
!       orthonormalize CMO
        Do iD = 1, nD
           Call Ortho(CMO(1,iD),nBO,Ovrlp,nBT)
        End Do
!       restore threshold values in WfCtl
        EThr=EThr_o
        Call Put_dScalar('EThr',EThr)
        FThr=FThr_o
        DThr=DThr_o
        DltNTh=DNTh_o
        Call xSet_ThrInt(ThrInt_o)
!       set twoel to AllCnt...
        twoel_NDDO=.FALSE.
!       close and reopen some DA files...
        Call DaClos(LuDSt)
        Call DaClos(LuOSt)
        Call DaClos(LuTSt)
        Call DaClos(LuGrd)
        Call DaClos(LuDGd)
        Call DaClos(Lux)
        Call DaClos(LuDel)
        Call DaClos(Luy)
        Call DAName(LuDSt,FnDSt)
        Call DAName(LuOSt,FnOSt)
        Call DAName(LuTSt,FnTSt)
        Call DAName(LuGrd,FnGrd)
        Call DAName(LuDGd,FnDGd)
        Call DAName(Lux,Fnx)
        Call DAName(LuDel,FnDel)
        Call DAName(Luy,Fny)
      Else
        nIterP=0
!       read kinetic energy matrix, use space for overlap matrix,
!       since that one is reread afterwards anyway...
        Label='Kinetic '
        iOpt = ibset(ibset(0,sNoOri),sNoNuc)
        iRC = -1
        iComp = 1
        Call RdOne(iRC,iOpt,Label,iComp,Ovrlp,lOper)
        If (iRC.ne.0) GoTo 9999
!       read NDDO NA matrix from ONEINT file...
        Label='AttractS'
        iOpt = ibset(ibset(0,sNoOri),sNoNuc)
        iRC = -1
        Call RdOne(iRC,iOpt,Label,iComp,OneHam,lOper)
        If (iRC.ne.0) GoTo 9999
!       and form NDDO one-electron Hamiltonian...
        Call DaXpY_(nBT,One,Ovrlp,1,OneHam,1)
!       read NDDO overlap matrix from ONEINT file...
        Label='MltplS 0'
        iOpt = ibset(ibset(0,sNoOri),sNoNuc)
        iRC = -1
        Call RdOne(iRC,iOpt,Label,iComp,Ovrlp,lOper)
        If (iRC.ne.0) GoTo 9999
!       save threshold values in WfCtl
        EThr_o=EThr
        FThr_o=FThr
        DThr_o=DThr
        DNTh_o=DltNTh
        ThrInt_o=Get_ThrInt()
!       and set new thresholds
!       EThr=EThr*1.0D+04
!       Call Put_dScalar('EThr',EThr)
!       FThr=FThr*1.0D+04
!       DThr=DThr*1.0D+04
!       DltNTh=DltNTh*1.0D+04
!       Call xSet_ThrInt(ThrInt_o*1.0D+04)
!       set twoel to OneCnt...
        twoel_NDDO=.TRUE.
      End If
!
      Return
!
!---- Error exit
 9999 Continue
      Write (6,*) 'SwiOpt: Error reading ONEINT'
      Write (6,'(A,A)') 'Label=',Label
      Call Abend()
!
      End SubRoutine SwiOpt
