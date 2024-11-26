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
      Subroutine InpCtl_MCLR(iPL)
************************************************************************
*                                                                      *
*     Read all relevant input data and display them                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use Str_Info, only: DTOC
      use negpre, only: nGP
      use ipPage, only: W
      use gugx, only: SGS, CIS, EXS
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: ipCI
      Implicit None
      Integer iPL
#include "Input.fh"
#include "Files_mclr.fh"
#include "sa.fh"
#include "detdim.fh"
#include "spinfo_mclr.fh"
#include "dmrginfo_mclr.fh"
      logical ldisk,ipopen
      Character(LEN=8) Method
      Real*8, Allocatable:: CIVec(:,:), CITmp(:)
      Integer i, ii, ipCII, iRC, iprDia, iSSM
      Integer, external:: ipGet, ipIn, ipOut
      Integer, External:: IsFreeUnit

! ==========================================================
      integer,allocatable::index_SD(:) ! not final version
      real*8,allocatable::vector_cidmrg(:)
! ==========================================================

*                                                                      *
************************************************************************
*                                                                      *
      Interface
        Subroutine RdJobIph_td(CIVec)
        Real*8, Allocatable:: CIVec(:,:)
        End Subroutine RdJobIph_td
        Subroutine RdJobIph(CIVec)
        Real*8, Allocatable:: CIVec(:,:)
        End Subroutine RdJobIph
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      !Read in interesting info from RUNFILE and ONEINT
      Call Rd1Int_MCLR()
      Call RdAB()   ! Read in orbitals, perturbation type, etc.
*                                                                      *
************************************************************************
*                                                                      *
      Call Rd2Int(iPL) ! Read in 2el header
*                                                                      *
************************************************************************
*                                                                      *
      Call RdInp_MCLR()  ! Read in input
*                                                                      *
************************************************************************
*                                                                      *
*     Default activate ippage utility
*
      ldisk  =ipopen(0,.True.)
*
      PT2 = .FALSE.
      Call Get_cArray('Relax Method',Method,8)
      If (Method.eq.'CASPT2  ') Then
        PT2 = .TRUE.
        !! Read the states requested by CASPT2
        !! This means that the root(s) specified in &MCLR is usually
        !! ignored for CASPT2 gradient/NAC.
        Call check_caspt2(1)
      End If
*
C     write(6,*) "iMethod:",iMethod,iCASSCF
      If (iMethod.eq.iCASSCF) Then
         If (TimeDep) Then
            Call RdJobIph_td(CIVec)
         Else
            Call RdJobIph(CIVec)
         End If

*        Write(6,*) 'Setup of Determinant tables'
         Call DetCtl   ! set up determinant tables
*....... Read in tables from disk
         Call InCsfSD(State_sym,State_sym,.true.)
*                                                                      *
************************************************************************
*                                                                      *
*        Write(6,*) 'Transformation of CI vector to symmetric '
*    &             ,'group from GUGA pepresentation'

         !> scratch  ! yma testing
!         if(doDMRG.and.doMCLR)then
!           call xflush(117)
!           close(117)
!         end if

         Do i=1,nroots
           ! yma
!          No need to copy,since there are no CI-vectors
           if(doDMRG.and.doMCLR)then
             Call mma_allocate(CITmp,ndets_RGLR,Label='CITmp')
           else
             Call mma_allocate(CITmp,nconf,Label='CITmp')
             call dcopy_(nconf,CIVec(:,i),1,CITmp,1)
           end if

           !> If doDMRG
           if(doDMRG.and.doMCLR)then ! yma
           else
             ! transform to sym. group
             Call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,
     &                    nRs1,nRs2,nRs3,
     &                    SGS,CIS,EXS,CITmp,1,State_Sym,State_Sym)
             NCSF(1:nSym)=CIS%NCSF(1:nSym)
             NCONF=CIS%NCSF(State_Sym)
             Call mkGuga_Free(SGS,CIS,EXS)

           end if

! Here should be the position for introducing the CI(SR) coefficients
!           iSSM=1     ! yma
!           write(6,*)"Set ISSM eq 1 ",ISSM

           if(doDMRG)then !yma
             call mma_allocate(index_SD,ndets_RGLR,label='index_SD')
             call mma_allocate(vector_cidmrg,ndets_RGLR,
     &                         label='vector_cidmrg')
             call ci_reconstruct(i,ndets_RGLR,vector_cidmrg,index_SD)
             do ii=1,ndets_RGLR
               if(abs(vector_cidmrg(ii)).lt.0.0d0)then
                 vector_cidmrg(ii)=0.0d0
               end if
             end do
             call CSDTVC_dmrg(CITmp,vector_cidmrg,2,DTOC,
     &                     index_SD,ISSM,1,IPRDIA)
             call mma_deallocate(index_SD)
             call mma_deallocate(vector_cidmrg)
           end if

           call dcopy_(nconf,CITmp,1,CIVec(:,i),1)
           Call mma_deallocate(CITmp)
        End Do
*                                                                      *
************************************************************************
*                                                                      *
         ldisk  =ipopen(nconf,page)
*
*        If we are computing Lagrangian multipliers we pick up all CI
*        vectors. For Hessian calculations we pick up just one vector.
*
C        Write (*,*) 'iState,SA,nroots=',iState,SA,nroots
         If (SA.or.iMCPD.or.PT2) Then
            ipcii=ipget(nconf*nroots)
            irc=ipin(ipcii)
            call dcopy_(nconf*nroots,CIVec,1,W(ipcii)%Vec,1)
            nDisp=1
         Else
            ipcii=ipget(nconf)
            irc=ipin(ipcii)
            call dcopy_(nConf,CIVec(:,iState),1,W(ipcii)%Vec,1)
            If (iRoot(iState).ne.1) Then
               Write (6,*) 'McKinley does not support computation of'
     &                   //' harmonic frequencies of excited states'
               Call Abend()
            End If
         End If
C        irc=ipin(ipcii)
C        Call RecPrt('CI vector',' ',W(ipcii)%Vec,1,nConf)
         Call mma_deallocate(CIVec)
*
*        At this point we change to ipci being the index of the CI
*        vector in the ipage utility.
*
         ipci=ipcii
         irc=ipout(ipci)
*                                                                      *
************************************************************************
*                                                                      *
         If (ngp) Call rdciv()
         If (PT2) Then
           LuPT2 = isFreeUnit(LuPT2)
           Call Molcas_Open(LuPT2,'PT2_Lag')
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call InpOne()         ! read in oneham
      Call PrInp_MCLR(iPL)  ! Print all info
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _WARNING_WORKAROUND_
      If (.False.) Then
         Call Unused_integer(irc)
         Call Unused_logical(ldisk)
      End If
#endif
      End Subroutine InpCtl_MCLR
