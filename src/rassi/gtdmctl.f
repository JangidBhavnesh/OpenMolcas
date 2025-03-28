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
      SUBROUTINE GTDMCTL(PROP,JOB1,JOB2,OVLP,DYSAMPS,NZ,IDISK)

      use rasdef, only: NRAS, NRASEL, NRS1, NRS1T, NRS2, NRS2T, NRS3,
     &                  NRS3T, NRSPRT
      use rassi_global_arrays, only: PART, OrbTab, HAM, SFDYS, SSTAB,
     &                               REST1, REST2,
     &                               CnfTab1, CnfTab2,
     &                               FSBTAB1, FSBTAB2,
     &                               SPNTAB1, SPNTAB2,
     &                               TRANS1, TRANS2
#ifdef _DMRG_
      use rassi_global_arrays, only: LROOT
#endif
      !> module dependencies
#ifdef _DMRG_
      use qcmaquis_interface_cfg
      use qcmaquis_interface_utility_routines, only:
     &    pretty_print_util
      use qcmaquis_info
      use qcmaquis_interface_mpssi
#endif
      use mspt2_eigenvectors
      use rasscf_global, only: DoDMRG
      use rassi_aux, only : AO_Mode, ipglob, iDisk_TDM, jDisk_TDM
      use gugx, only: SGStruct, CIStruct, EXStruct
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
C      use para_info, only: nProcs, is_real_par, king
#ifdef _HDF5_
      use mh5, only: mh5_put_dset
      use Cntrl, only: CIH5
      use RASSIWfn, only: wfn_CMO, wfn_CMO_OR, wfn_DetCoeff,
     &                    wfn_DetCoeff_OR, wfn_DetOcc,
     &                    wfn_DetOcc_OR
#endif
      use frenkel_global_vars, only: DoCoul
      use Constants, only: auToEV, Half, One, Zero
      use cntrl, only: sonatnstate
      use Cntrl, only: NSTATE, NPROP, LSYM1, LSYM2, IFHEXT,
     &                 IFEJOB, TDYS, DYSO, DCHS, NATO, DOGSOR,
     &                 QDPT2EV, ERFNUC, NCONF1, NCONF2, PRCI, DCHO,
     &                 IFNTO, SAVEDENS, QDPT2SC, IFTRD1, IFTRD2, NDET,
     &                 IFHAM, IFHEFF, SECOND_TIME, CITHR, IRREP, ISTAT,
     &                 JBNAME, MLTPLT, NACTE, NELE3, NHOLE1, NSTAT,
     &                 RASTYP
      use cntrl, only: iToc15, LuIph, LuTDM
      use Symmetry_Info, only: nSym=>nIrrep, MUL
      use rassi_data, only: NASHT,NISHT,NCMO,ENUC,NASH,NDEL,NFRO,NISH,
     &                      NOSH,NSSH,NTDMAB,NTDMZZ,NTRA

      IMPLICIT NONE
      Type (SGStruct), Target :: SGS(2)
      Type (CIStruct) :: CIS(2)
      Type (EXStruct) :: EXS(2)
      Real*8 PROP(NSTATE,NSTATE,NPROP)
      Integer NGASORB(100),NGASLIM(2,10)
      Integer NASHES(8)
      Real*8 OVLP(NSTATE,NSTATE)
      Real*8 DYSAMPS(NSTATE,NSTATE)
      LOGICAL IF00, IF10,IF01,IF20,IF11,IF02,IF21,IF12,IF22
      LOGICAL IFTWO,TRORB
      CHARACTER(LEN=8) WFTP1,WFTP2
      CHARACTER(LEN=48) STLNE2
      Real*8 Energies(1:20)
      Integer IAD,LUIPHn,LUCITH
      Real*8 Norm_fac
CC prepare the parallel infrastructure for (istate,jstate loop)
C      integer :: itask, ltask, ltaski, ltaskj, ntasks
C#ifdef _MOLCAS_MPP_
C      logical :: Rsv_Tsk
C      integer :: ID
C#endif
CC
C     transformed CI expansion in determinant basis
      real(kind=wp), allocatable :: detcoeff1(:), detcoeff2(:)
      character(len=(NASHT+1)), allocatable :: detocc(:)
CC CC
CC    NTO section
      Logical DoNTO
CC    NTO section

      type mixed_1pdensities
        real*8              :: overlap
        real*8, allocatable :: rtdm(:)
        real*8, allocatable :: stdm(:)
        real*8, allocatable :: wtdm(:)
      end type

      type(mixed_1pdensities), allocatable :: mstate_1pdens(:,:)

      logical               :: mstate_dens
      real*8                :: fac1, fac2
      real*8, Allocatable:: CMO1(:), CMO2(:)
      real*8, Allocatable:: TRAD(:), TRASD(:), WERD(:)
      real*8, Allocatable:: TDMAB(:), TSDMAB(:), WDMAB(:)
      real*8, Allocatable:: TDMZZ(:), TSDMZZ(:), WDMZZ(:)
      real*8, Allocatable:: TDM2(:), TRA1(:), TRA2(:), FMO(:), TUVX(:)
      real*8, Allocatable:: DYSCOF(:), DYSAB(:), DYSZZ(:)
      Integer NRT2M,NRT2MAB,AUGSPIN,NDCHSM,DCHIJ
      Integer ISY,JSY,LSY,NI,NJ,NL
      real*8, Allocatable:: RT2M(:),RT2MAB(:)
      real*8, Allocatable:: DCHSM(:)
      real*8 BEi,BEj,BEij
      Integer, Allocatable:: OMAP(:)
      real*8, allocatable:: CI1(:), CI2(:), CI2_o(:)
      real*8, pointer:: DET1(:), DET2(:)
      real*8, allocatable, target:: DETTOT1(:,:), DETTOT2(:,:)
      real*8, allocatable:: Theta1(:), ThetaN(:), ThetaM(:)

      Integer nActE1, JOB1, MPLET1, NHOL11, NELE31, NACTE2, JOB2,
     &        MPLET2, NHOL12, NASORB, NTDM1, NTDM2, ISY12, NDYSAB,
     &        NDYSZZ, NZ, NTRAD, I, J, MSPROJ1, MSPROJ2, NASPRT, KSPART,
     &        KOINFO, ISUM, ISYM, ISORB, NO, IO, ISOIND, IT, ITABS,
     &        NPART, NGAS, IE1MN, IE1MX, IE3MN, IE3MX, IE2MN, IE2MX,
     &        NGL11, NGL21, NGL12, NGL22, NGL13, NGL23, MAXOP, IE1,
     &        IOP1, IE3, IOP3, IE2, IOP2, IFORM, MINOP, NDET1, NDET2,
     &        IST, ISTATE, JST, JSTATE, IRC, IJ, IDISK, IOPT, IGO, IERR,
     &        NELE32, ISPART, IEMPTY
      REAL*8 ECORE, Dot_Prod, HZERO, HONE, HTWO, SIJ, DYSAMP, DYNORM,
     &       HII, HJJ, HIJ, OVERLAP_RASSI
      REAL*8, External:: DDot_
      Integer, External:: IsFreeUnit

      Interface
      Subroutine SGInit(nSym,nActEl,iSpin,SGS,CIS)
      use gugx, only: SGStruct, CIStruct
      IMPLICIT None
      Integer nSym, nActEl, iSpin
      Type (SGStruct), Target :: SGS
      Type (CIStruct) :: CIS
      End Subroutine SGInit
      End Interface

#define _TIME_GTDM
#ifdef _TIME_GTDM_
      Call CWTime(TCpu1,TWall1)
#endif
#ifdef _WARNING_WORKAROUND_
* Avoid compiler warnings about possibly unitialised mstate_1pdens
* The below can be removed if the file is compiled with
* -Wno-error=maybe-uninitialized
      allocate(mstate_1pdens(0,0))
      deallocate(mstate_1pdens)
#endif
C WF parameters for ISTATE and JSTATE

      NACTE1=NACTE(JOB1)
      MPLET1=MLTPLT(JOB1)
      LSYM1=IRREP(JOB1)
      NHOL11=NHOLE1(JOB1)
      NELE31=NELE3(JOB1)
      WFTP1=RASTYP(JOB1)
      NACTE2=NACTE(JOB2)
      MPLET2=MLTPLT(JOB2)
      LSYM2=IRREP(JOB2)
      NHOL12=NHOLE1(JOB2)
      NELE32=NELE3(JOB2)
      WFTP2=RASTYP(JOB2)
      SGS(1)%IFRAS=1
      SGS(2)%IFRAS=1
      IF(IPGLOB.GE.4) THEN
        WRITE(6,*)' Entered GTDMCTL.'
        WRITE(6,'(1X,A,I3,A,I3)')'  JOB1:',  JOB1,'        JOB2:',  JOB2
        WRITE(6,'(1X,A,I3,A,I3)')'NACTE1:',NACTE1,'      NACTE2:',NACTE2
        WRITE(6,'(1X,A,I3,A,I3)')'MPLET1:',MPLET1,'      MPLET2:',MPLET2
        WRITE(6,'(1X,A,I3,A,I3)')' LSYM1:', LSYM1,'       LSYM2:', LSYM2
        WRITE(6,'(1X,A,A8,A,A8)')' WFTP1:',WFTP1, '       WFTP2:',WFTP2
        WRITE(6,'(1X,A,I3,A,I3)')' NPROP:', NPROP
      END IF

      if(doDMRG .and. (nacte1 /= nacte2 ))then
        call WarningMessage(2,'Problem in gtdmctl for MPS-SI: no '//
     &   ' match for #e- in bra/ket')
        call abend()
      end if

      !> Logical variables, controlling which GTDM''s to compute

      !>  Overlap
      IF00 = NACTE1.EQ.NACTE2.AND. MPLET1.EQ.MPLET2
      IF00 = IF00.AND.LSYM1.EQ.LSYM2

      !> Dyson amplitudes
      IF10 = (NACTE1-NACTE2).EQ. 1
      IF10 = IF10.AND.( ABS(MPLET1-MPLET2).EQ.1)
      IF01 = (NACTE1-NACTE2).EQ.-1
      IF01 = IF01.AND.( ABS(MPLET1-MPLET2).EQ.1)

      !> Pair amplitudes:
      IF20 = (NACTE1-NACTE2).EQ. 2
      IF02 = (NACTE1-NACTE2).EQ.-2

      !> 1-TDMs and transition spin densities
      IF11 = ( NACTE1.EQ.NACTE2.AND.NACTE1.GE.1)
      IF11 = IF11.AND.( ABS(MPLET1-MPLET2).LE.2)

      !> 2h1p and 1h2p amplitudes:
      IF21 = IF10.AND.NACTE2.GE.1
      IF21 = IF21.AND.( ABS(MPLET1-MPLET2).LE.3)
      IF12 = IF01.AND.NACTE1.GE.1
      IF12 = IF12.AND.( ABS(MPLET1-MPLET2).LE.3)

      !> 2-TDMs and transition spin densities
      IF22 = (NACTE1.EQ.NACTE2.AND.NACTE1.GE.2)
      IF22 = IF22.AND.(ABS(MPLET1-MPLET2).LE.4)

      !> check if they are needed at all:
        !> It may be that the Hamiltonian matrix should be used in
        !> diagonalization (IFHAM is .TRUE.), but it does not have
        !> to be computed (because IFHEXT or IFHEFF or IFEJOB are true).
      IFTWO = IFHAM .AND..NOT.(IFHEXT.OR.IFHEFF.OR.IFEJOB)

      !> For the moment, we have no use for the two-electron density
      !> except when used for the scalar two-body Hamiltonian matrix:
      IF22 = IF22.AND.IFTWO.AND.(MPLET1.EQ.MPLET2).AND.(LSYM1.EQ.LSYM2)

      IF(IPGLOB.GE.4) THEN
        IF(IF00) WRITE(6,*)' Overlap will be computed.'
        IF(IF10.or.IF01) WRITE(6,*)' Dyson orbital will be computed.'
        IF(IF20.or.IF02) WRITE(6,*)' Pair amplitudes will be computed.'
        IF(IF11) WRITE(6,*)' Density 1-matrix will be computed.'
        IF(IF21.or.IF12) WRITE(6,*)' 2h1p amplitudes will be computed.'
        IF(IF22) WRITE(6,*)' Density 2-matrix will be computed.'
      END IF

C Pick up orbitals of ket and bra states.
      Call mma_allocate(CMO1,nCMO,Label='CMO1')
      Call mma_allocate(CMO2,nCMO,Label='CMO2')
      CALL RDCMO_RASSI(JOB1,CMO1)
      CALL RDCMO_RASSI(JOB2,CMO2)

C Nr of active spin-orbitals
      NASORB=2*NASHT
      NTDM1=NASHT**2
      NTDM2=(NTDM1*(NTDM1+1))/2

C Size of some data sets of reduced-2TDM in terms of active
C orbitals NASHT (For Auger matrix elements):
      NRT2M=NASHT**3
C Size of Symmetry blocks
      ISY12=MUL(LSYM1,LSYM2)
      NRT2MAB=0
      DO ISY=1,NSYM
       NI=NOSH(ISY)
       IF(NI.EQ.0) GOTO 200
       DO JSY=1,NSYM
        NJ=NOSH(JSY)
        IF(NJ.EQ.0) GOTO 300
        DO LSY=1,NSYM
         NL=NOSH(LSY)
         IF(NL.EQ.0) GOTO 400
          IF(MUL(ISY,MUL(JSY,LSY)).EQ.ISY12) THEN
           NRT2MAB=NRT2MAB+NI*NJ*NL
          END IF
400     CONTINUE
        END DO
300    CONTINUE
       END DO
200   CONTINUE
      END DO
      IF (TDYS.and..not.DYSO) THEN
       Write(6,*) ' '
       Write(6,*) 'Auger (TDYS) requires Dyson calculation.'
       Write(6,*) 'Make sure to activate Dyson in your RASSI input.'
       Write(6,*) 'For now, Auger computation will be skipped.'
      END IF
C evaluation of DCH
      IF(DCHS) THEN
       NDCHSM=NASHT**2
      END IF

! +++ J. Norell 13/7 - 2018
C 1D arrays for Dyson orbital coefficients
C COF = active biorthonormal orbital base
C AB  = inactive+active biorthonormal orbital base
C ZZ  = atomic (basis function) base
      IF ((IF10.or.IF01).and.DYSO) THEN
        Call mma_allocate(DYSCOF,NASORB,Label='DYSCOF')
        ! Number of inactive+active orbitals
        NDYSAB = NASHT+NISHT
        Call mma_allocate(DYSAB,nDYSAB,Label='DYSAB')
        ! Number of atomic / basis functions
        NDYSZZ = NZ
        Call mma_allocate(DYSZZ,nDYSZZ,Label='DYSZZ')
        DYSZZ(:) = Zero
      END IF
! +++

C Transition density matrices, TDMAB is for active biorthonormal
C orbitals only, while TDMZZ is in the fixed AO basis.
C WDMAB, WDMZZ similar, but WE-reduced 'triplet' densities.
      IF(IF11.AND.(NATO.OR.NPROP.GT.0)) THEN
        Call mma_allocate(TDMAB,nTDMAB,Label='TDMAB')
        Call mma_allocate(TSDMAB,nTDMAB,Label='TSDMAB')
        Call mma_allocate(WDMAB,nTDMAB,Label='WDMAB')
        Call mma_allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
        Call mma_allocate(TSDMZZ,nTDMZZ,Label='TSDMZZ')
        Call mma_allocate(WDMZZ,nTDMZZ,Label='WDMZZ')
      END IF

      IF (IF11) THEN
        NTRAD=NASHT**2 ! NTRAD == NWERD == NTRASD
        Call mma_allocate(TRAD,nTRAD+1,Label='TRAD')
        Call mma_allocate(TRASD,nTRAD+1,Label='TRASD')
        Call mma_allocate(WERD,nTRAD+1,Label='WERD')
      END IF
      IF (IF22) THEN
        Call mma_allocate(TDM2,nTDM2,Label='TDM2')
      ELSE
        ! To avoid passing an unallocated argument
        Call mma_allocate(TDM2,0,Label='TDM2')
      END IF

      IF(JOB1.NE.JOB2) THEN
C Transform to biorthonormal orbital system
        IF (DoGSOR) Then
          Call FCopy(Trim(JBNAME(JOB2)),'JOBGS',ierr)
          Call DANAME(LUIPH,'JOBGS')
          IAD = 0
          Call IDAFile(LUIPH,2,ITOC15,30,IAD)
          IAD=ITOC15(2)
          Call DDAFile(LUIPH,1,CMO2,nCMO,IAD)
          Call DACLOS(LUIPH)
        End if !DoGSOR


        Call mma_allocate(TRA1,nTRA,Label='TRA1')
        Call mma_allocate(TRA2,nTRA,Label='TRA2')
        CALL FINDT(CMO1,CMO2,TRA1,TRA2)
#ifdef _HDF5_
C       put the pair of transformed orbitals to h5
        if (CIH5) then
          call mh5_put_dset(wfn_cmo,CMO1,[nCMO,1],[0,JOB1-1])
          call mh5_put_dset(wfn_cmo,CMO2,[nCMO,1],[0,JOB2-1])
        endif
#endif
        TRORB = .true.
      else
#ifdef _HDF5_
C       put original orbitals to hdf5 file
        if (CIH5) call mh5_put_dset(wfn_cmo_or,CMO1,
     &                               [nCMO,1],[0,JOB1-1])
#endif
        TRORB = .false.
      end if

!     > check whether we do RASSI with an effective multi-state PT2 Hamiltonian
!     > whose eigenvectors are stored in Heff_evc
!     > i.e., we do not use mixed CI coefficients / MPS wave functions but rather mix the TDMs

      mstate_dens = job1.eq.job2.and.
     &              (allocated(Heff_evc(job1)%pc).or.
     &               allocated(Heff_evc(job1)%sc))
      mstate_dens = mstate_dens.and.if11
      mstate_dens = mstate_dens.and.qdpt2ev

      if(mstate_dens)then
        allocate(mstate_1pdens(nstat(job1),nstat(job1)))
        do i = 1, nstat(job1)
          do j = 1, nstat(job1)
            call mma_allocate(mstate_1pdens(i,j)%rtdm, NTDMZZ)
            call mma_allocate(mstate_1pdens(i,j)%stdm, NTDMZZ)
            call mma_allocate(mstate_1pdens(i,j)%wtdm,NTDMZZ)
            mstate_1pdens(i,j)%rtdm = 0; mstate_1pdens(i,j)%stdm    = 0
            mstate_1pdens(i,j)%wtdm = 0; mstate_1pdens(i,j)%overlap = 0
          end do
        end do
      end if

#ifdef _DMRG_
      dmrg_external%MPSrotated = trorb
#endif

C OBTAIN CORE ENERGY, FOCK MATRIX, AND TWO-ELECTRON INTEGRALS
C IN THE MIXED ACTIVE MO BASIS:
      ECORE = Zero
      IF (IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
       Call mma_allocate(FMO,nTDM1,Label='FMO')
       Call mma_allocate(TUVX,nTDM2,Label='TUVX')
       TUVX(:) = Zero
CTEST       write(*,*)'GTDMCTL calling TRINT.'
       CALL TRINT(CMO1,CMO2,ECORE,nTDM1,FMO,nTDM2,TUVX)
       ECORE=ENUC+ERFNUC+ECORE
CTEST       write(*,*)'GTDMCTL back from TRINT.'
CTEST       write(*,*)'ENUC  =',ENUC
CTEST       write(*,*)'ERFNUC=',ERFNUC
CTEST       write(*,*)'ECORE =',ECORE
      END IF

C In the calculation of matrix elements ( S1, S2 ), we will use
C the same Ms quantum numbers, if S1 and S2 differ by 0 or an int,
C else we will use Ms quantum numbers that differ by 1/2:
      IF( MOD(ABS(MPLET1-MPLET2),2).EQ.0 ) THEN
        MSPROJ1=MIN(MPLET1,MPLET2)-1
        MSPROJ2=MSPROJ1
      ELSE
        IF(MPLET1.GT.MPLET2) THEN
          MSPROJ2=MPLET2-1
          MSPROJ1=MSPROJ2+1
        ELSE
          MSPROJ1=MPLET1-1
          MSPROJ2=MSPROJ1+1
        END IF
      END IF

#ifdef _DMRG_
      !> set spin-up/spin-down # of electrons for target state(s)
      if(doDMRG)then
        dmrg_external%nalpha = (nacte1 + msproj1) / 2
        dmrg_external%nbeta  = (nacte1 - msproj1) / 2
      end if
#endif

C---------------  For all wave functions: ---------------------
C Define structures ('tables') pertinent all jobs.
C (Later, move this up before the GTDMCTL calls).
C These are at:
C PART
C ORBTAB
C SSTAB
      Call NEWPRTTAB(NSYM,NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NDEL)
      IF(IPGLOB.GE.4) CALL PRPRTTAB(PART)

      Call NEWORBTAB(PART)
      IF(IPGLOB.GE.4) CALL PRORBTAB(ORBTAB)

      Call NEWSSTAB(ORBTAB)
      IF(IPGLOB.GE.4) CALL PRSSTAB(SSTAB)

C Mapping from active spin-orbital to active orbital in external order.
C Note that these differ, not just because of the existence of two
C spin-orbitals for each orbital, but also because the active orbitals
C (external order) are grouped by symmetry and then RAS space, but the
C spin orbitals are grouped by subpartition.
      CALL mma_allocate(OMAP,NASORB,Label='OMAP')
      NASPRT=ORBTAB(9)
      KSPART=ORBTAB(10)
      KOINFO=19
      ISUM=0
      DO ISYM=1,NSYM
        NASHES(ISYM)=ISUM
        ISUM=ISUM+NASH(ISYM)
      END DO
      ISORB=0
      DO ISPART=1,NASPRT
        NO=OrbTab(KSPART-1+ISPART)
        DO IO=1,NO
          ISORB=ISORB+1
C Orbital symmetry:
          ISYM=OrbTab(KOINFO+1+(ISORB-1)*8)
C In-Symmetry orbital index:
          ISOIND=OrbTab(KOINFO+2+(ISORB-1)*8)
C Subtract nr of inactive orbitals in that symmetry:
          IT=ISOIND-NISH(ISYM)
C Add nr of actives in earlier symmetries:
          ITABS=NASHES(ISYM)+IT
          OMAP(ISORB)=ITABS
        END DO
      END DO

C---------------    JOB1 wave functions: ---------------------
C Initialize SGUGA tables for JOB1 functions.
C These are structures stored in user defined types:
C SGS(1),CIS(1) and EXS(1).

C Set variables in /RASDEF/, used by SGUGA codes, which define
C the SGUGA space of JOB1. General RAS:
      IF(WFTP1.EQ.'GENERAL ') THEN
        NRSPRT=3
        DO I=1,8
          NRAS(I,1)=NRS1(I)
          NRAS(I,2)=NRS2(I)
          NRAS(I,3)=NRS3(I)
        END DO
        NRASEL(1)=2*NRS1T-NHOL11
        NRASEL(2)=NACTE1-NELE31
        NRASEL(3)=NACTE1

        if(.not.doDMRG)then
          CALL SGINIT(NSYM,NACTE1,MPLET1,SGS(1),CIS(1))
          IF(IPGLOB.GT.4) THEN
            WRITE(6,*)'Split-graph structure for JOB1=',JOB1
            CALL SGPRINT(SGS(1))
          END IF
          CALL CXINIT(SGS(1),CIS(1),EXS(1))
C CI sizes, as function of symmetry, are now known.
          NCONF1=CIS(1)%NCSF(LSYM1)
        else
          NCONF1=1
        end if
      ELSE
C Presently, the only other cases are HISPIN, CLOSED or EMPTY.
* Note: the HISPIN case may be buggy and is not used presently.
        NCONF1=1
      END IF
      CALL mma_allocate(CI1,NCONF1,Label='CI1')

C Still JOB1, define structures ('tables') pertinent to JOB1
C These are at:
C REST1
C CNFTAB1
C FSBTAB1
C SPNTAB1

      NPART=3
      NGAS=NPART
      DO ISYM=1,NSYM
        NGASORB(ISYM)=NRS1(ISYM)
        NGASORB(ISYM+NSYM)=NRS2(ISYM)
        NGASORB(ISYM+2*NSYM)=NRS3(ISYM)
      END DO

*PAM2008: The old MAXOP was far too generous:
*      MAXOP=NASHT
*PAM2008: MAXOP is determined by RAS restrictions:
* Preliminary ranges of nr of electrons:
      IE1MN=MAX(0,2*NRS1T-NHOL11)
      IE1MX=2*NRS1T
      IE3MN=MAX(0,NACTE1-IE1MX-2*NRS2T)
      IE3MX=MIN(2*NRS3T,NELE31)
      IE2MN=MAX(0,NACTE1-IE1MX-IE3MX)
      IE2MX=MIN(2*NRS2T,NACTE1-IE1MN-IE3MN)
* Preliminary NGASLIM:
      NGL11=IE1MX
      NGL21=IE1MN
      NGL12=IE2MX
      NGL22=IE2MN
      NGL13=IE3MX
      NGL23=IE3MN
* Start with MAXOP=0, then increase:
      MAXOP=0
* Loop over possible ranges:
      DO IE1=IE1MN,IE1MX
       IOP1=MIN(IE1,(2*NRS1T-IE1))
       IF(IOP1.GE.0) THEN
        DO IE3=IE3MN,IE3MX
         IOP3=MIN(IE3,(2*NRS3T-IE3))
         IF(IOP3.GE.0) THEN
          IE2=NACTE1-IE1-IE3
          IOP2=MIN(IE2,(2*NRS2T-IE2))
          IF(IOP2.GE.0) THEN
* Actually possible combination:
           MAXOP=MAX(MAXOP,IOP1+IOP2+IOP3)
           NGL11=MIN(NGL11,IE1)
           NGL21=MAX(NGL21,IE1)
           NGL12=MIN(NGL12,IE2)
           NGL22=MAX(NGL22,IE2)
           NGL13=MIN(NGL13,IE3)
           NGL23=MAX(NGL23,IE3)
          END IF
         END IF
        END DO
       END IF
      END DO
      NGASLIM(1,1)=NGL11
      NGASLIM(2,1)=NGL21
      NGASLIM(1,2)=NGL12
      NGASLIM(2,2)=NGL22
      NGASLIM(1,3)=NGL13
      NGASLIM(2,3)=NGL23

      if(.not.doDMRG)then
        IFORM=1
        MINOP=0
        Call NEWGASTAB(NSYM,NGAS,NGASORB,NGASLIM,1)
        IF(IPGLOB.GE.4) CALL PRGASTAB(REST1)

C At present, we will only annihilate, at most 2 electrons will
C be removed. This limits the possible MAXOP:
        MAXOP=MIN(MAXOP+1,NACTE1,NASHT)
        Call NEWCNFTAB(NACTE1,NASHT,MINOP,MAXOP,LSYM1,NGAS,
     &                     NGASORB,NGASLIM,IFORM,1)
        IF(IPGLOB.GE.4) CALL PRCNFTAB(CNFTAB1,100)

        Call NEWFSBTAB(NACTE1,MSPROJ1,LSYM1,REST1,SSTAB,1)
        IF(IPGLOB.GE.4) CALL PRFSBTAB(FSBTAB1)
        NDET1=FSBTAB1(5)
        if (ndet1 /= ndet(job1)) ndet(job1) = ndet1
        Call NEWSCTAB(MINOP,MAXOP,MPLET1,MSPROJ1,1)
        IF (IPGLOB.GT.4) THEN
*PAM2009: Put in impossible call to PRSCTAB, just so code analyzers
* do not get their knickers into a twist.
          CALL PRSCTAB(SPNTAB1,TRANS1)
        END IF
      else
        NDET1 = 1 ! minimum to avoid runtime error
      end if
C---------------    JOB2 wave functions: ---------------------
C Initialize SGUGA tables for JOB2 functions.
C These are structures stored in arrays:
C SGS(2),CIS(2) and EXS(2).

C Set variables in /RASDEF/, used by SGUGA codes, which define
C the SGUGA space of JOB1. General RAS:
      IF(WFTP2.EQ.'GENERAL ') THEN
        NRSPRT=3
        DO I=1,8
          NRAS(I,1)=NRS1(I)
          NRAS(I,2)=NRS2(I)
          NRAS(I,3)=NRS3(I)
        END DO
        NRASEL(1)=2*NRS1T-NHOL12
        NRASEL(2)=NACTE2-NELE32
        NRASEL(3)=NACTE2

        IF(.not.doDMRG)then
          CALL SGINIT(NSYM,NACTE2,MPLET2,SGS(2),CIS(2))
          IF(IPGLOB.GT.4) THEN
            WRITE(6,*)'Split-graph structure for JOB2=',JOB2
            CALL SGPRINT(SGS(2))
          END IF
          CALL CXINIT(SGS(2),CIS(2),EXS(2))
C CI sizes, as function of symmetry, are now known.
          NCONF2=CIS(2)%NCSF(LSYM2)
        else
          NCONF2=1
        end if
      ELSE
C Presently, the only other cases are HISPIN, CLOSED or EMPTY.
* Note: the HISPIN case may be buggy and is not used presently.
        NCONF2=1
      END IF
      CALL mma_allocate(CI2,NCONF2,Label='CI2')
      If (DoGSOR) Then
        CALL mma_allocate(CI2_o,NCONF2,Label='CI2_o')
      end if!DoGSOR

      NPART=3
      NGAS=NPART
      DO ISYM=1,NSYM
        NGASORB(ISYM)=NRS1(ISYM)
        NGASORB(ISYM+NSYM)=NRS2(ISYM)
        NGASORB(ISYM+2*NSYM)=NRS3(ISYM)
      END DO
*PAM2008: The old MAXOP was far too generous:
*      MAXOP=NASHT
*PAM2008: MAXOP is determined by RAS restrictions:
* Preliminary ranges of nr of electrons:
      IE1MN=MAX(0,2*NRS1T-NHOL12)
      IE1MX=2*NRS1T
      IE3MN=MAX(0,NACTE2-IE1MX-2*NRS2T)
      IE3MX=MIN(2*NRS3T,NELE32)
      IE2MN=MAX(0,NACTE2-IE1MX-IE3MX)
      IE2MX=MIN(2*NRS2T,NACTE2-IE1MN-IE3MN)
* Preliminary NGASLIM:
      NGL11=IE1MX
      NGL21=IE1MN
      NGL12=IE2MX
      NGL22=IE2MN
      NGL13=IE3MX
      NGL23=IE3MN
* Start with MAXOP=0, then increase:
      MAXOP=0
* Loop over possible ranges:
      DO IE1=IE1MN,IE1MX
       IOP1=MIN(IE1,(2*NRS1T-IE1))
       IF(IOP1.GE.0) THEN
        DO IE3=IE3MN,IE3MX
         IOP3=MIN(IE3,(2*NRS3T-IE3))
         IF(IOP3.GE.0) THEN
          IE2=NACTE2-IE1-IE3
          IOP2=MIN(IE2,(2*NRS2T-IE2))
          IF(IOP2.GE.0) THEN
* Actually possible combination:
           MAXOP=MAX(MAXOP,IOP1+IOP2+IOP3)
           NGL11=MIN(NGL11,IE1)
           NGL21=MAX(NGL21,IE1)
           NGL12=MIN(NGL12,IE2)
           NGL22=MAX(NGL22,IE2)
           NGL13=MIN(NGL13,IE3)
           NGL23=MAX(NGL23,IE3)
          END IF
         END IF
        END DO
       END IF
      END DO
      NGASLIM(1,1)=NGL11
      NGASLIM(2,1)=NGL21
      NGASLIM(1,2)=NGL12
      NGASLIM(2,2)=NGL22
      NGASLIM(1,3)=NGL13
      NGASLIM(2,3)=NGL23

      if(.not.dodmrg)then
        CALL NEWGASTAB(NSYM,NGAS,NGASORB,NGASLIM,2)
        IF(IPGLOB.GE.4) CALL PRGASTAB(REST2)

        IFORM=1
        MINOP=0
C At present, we will only annihilate. This limits the possible MAXOP:
        MAXOP=MIN(MAXOP+1,NACTE2,NASHT)
        Call NEWCNFTAB(NACTE2,NASHT,MINOP,MAXOP,LSYM2,NGAS,
     &                     NGASORB,NGASLIM,IFORM,2)
        IF(IPGLOB.GE.4) CALL PRCNFTAB(CNFTAB2,100)

        Call NEWFSBTAB(NACTE2,MSPROJ2,LSYM2,REST2,SSTAB,2)
        IF(IPGLOB.GE.4) CALL PRFSBTAB(FSBTAB2)
        NDET2=FSBTAB2(5)
        if (ndet2 /= ndet(job2)) ndet(job2) = ndet2
        Call NEWSCTAB(MINOP,MAXOP,MPLET2,MSPROJ2,2)
        IF (IPGLOB.GT.4) THEN
*PAM2009: Put in impossible call to PRSCTAB, just so code analyzers
* do not get their knickers into a twist.
          CALL PRSCTAB(SPNTAB2,TRANS2)
        END IF
      else
        NDET2 = 1 ! minimum to avoid runtime error
      end if
C-------------------------------------------------------------
      call mma_allocate(DETTOT1,NDET1,NSTAT(JOB1),Label='DETTOT1')
      call mma_allocate(DETTOT2,NDET2,NSTAT(JOB2),Label='DETTOT2')
      call mma_allocate(detocc,max(nDet1,nDet2),label='detocc')

C Loop over the states of JOBIPH nr JOB1
      DO IST=1,NSTAT(JOB1)
        DET1=>DETTOT1(1:NDET1,IST)
        ISTATE=ISTAT(JOB1)-1+IST

        if(.not.doDMRG)then
C Read ISTATE wave function
          IF(WFTP1.EQ.'GENERAL ') THEN
            CALL READCI(ISTATE,SGS(1),CIS(1),NCONF1,CI1)
          ELSE
            CI1(1) = One
          END IF
          DET1(:)=0.0D0
C         Transform to bion basis, Split-Guga format
          If (TrOrb) CALL CITRA (WFTP1,SGS(1),CIS(1),EXS(1),
     &                           LSYM1,TRA1,NCONF1,CI1)
          call mma_allocate(detcoeff1,nDet1,label='detcoeff1')
          CALL PREPSD(WFTP1,SGS(1),CIS(1),LSYM1,
     &                CNFTAB1,SPNTAB1,
     &                SSTAB,FSBTAB1,NCONF1,CI1,
     &                DET1,detocc,detcoeff1,
     &                TRANS1)

C       print transformed ci expansion
        if (JOB1 /= JOB2) then
          if (PRCI) then
            call prwf_biorth(istate, job1, nconf1, ndet1, nasht,
     &                       detocc, detcoeff1, cithr)
          end if
#ifdef _HDF5_
C         put transformed ci coefficients for JOB1 to h5
          if (CIH5) then
            call mh5_put_dset(wfn_detcoeff,detcoeff1,
     &                        [nDet1,1],[0,istate-1])
            call mh5_put_dset(wfn_detocc, detocc(:nDet1),
     &                        [nDet1,1], [0,(JOB1-1)])
          end if
#endif
        else
#ifdef _HDF5_
C         JOB1=JOB2, put original ci coefficients for JOB1 to h5
          if (CIH5) then
            call mh5_put_dset(wfn_detcoeff_or,detcoeff1,
     &                        [nDet1,1],[0,istate-1])
            call mh5_put_dset(wfn_detocc_or, detocc(:nDet1),
     &                        [nDet1,1], [0,(JOB1-1)])
          end if
#endif
        end if

        call mma_deallocate(detcoeff1)

        else ! doDMRG
#ifdef _DMRG_
          call prepMPS(
     &                 TRORB,
     &                 LROOT(ISTATE),
     &                 LSYM1,
     &                 MPLET1,
     &                 MSPROJ1,
     &                 NACTE1,
     &                 TRA1,
     &                 NTRA,
     &                 NISH,
     &                 NASH,
     &                 NOSH,
     &                 NSYM,
     &                 6,
     &                 ISTATE,
     &                 job1,
     &                 ist
     &                )
#endif
        end if
      END DO

      If (DoGSOR) Then
        CALL mma_allocate(Theta1,NCONF2,Label='Theta1')
        Theta1(:)=0.0D0
      End If

C-------------------------------------------------------------

      DO JST=1,NSTAT(JOB2)
        DET2=>DETTOT2(1:NDET2,JST)
        JSTATE=ISTAT(JOB2)-1+JST
        if(.not.doDMRG)then
C Read JSTATE wave function
          IF(WFTP2.EQ.'GENERAL ') THEN
            CALL READCI(JSTATE,SGS(2),CIS(2),NCONF2,CI2)
          ELSE
            CI2(1) = One
          END IF
          If(DoGSOR) Then
            CALL DCOPY_(NCONF2,CI2,1,CI2_o,1)
          End If
          DET2(:)=0.0D0
C         Transform to bion basis, Split-Guga format
          If (TrOrb) CALL CITRA (WFTP2,SGS(2),CIS(2),EXS(2),
     &                           LSYM2,TRA2,NCONF2,CI2)
          call mma_allocate(detcoeff2,nDet2,label='detcoeff2')
          CALL PREPSD(WFTP2,SGS(2),CIS(2),LSYM2,
     &                CNFTAB2,SPNTAB2,
     &                SSTAB,FSBTAB2,NCONF2,CI2,
     &                DET2,detocc,detcoeff2,
     &                TRANS2)

C         print transformed ci expansion
          if (JOB1 /= JOB2) then
            if (PRCI) then
              call prwf_biorth(jstate, job2, nconf2, ndet2, nasht,
     &                         detocc, detcoeff2, cithr)
            end if
#ifdef _HDF5_
C           put ci coefficients for JOB2 to h5
            if (CIH5) then
              call mh5_put_dset(wfn_detcoeff,detcoeff2,
     &                          [nDet2,1],[0,jstate-1])
              call mh5_put_dset(wfn_detocc, detocc(:nDet2),
     &                          [nDet2,1], [0,(JOB2-1)])
            end if
#endif
          end if
          call mma_deallocate(detcoeff2)

        else
#ifdef _DMRG_
          call prepMPS(
     &                 TRORB,
     &                 lroot(JSTATE),
     &                 LSYM2,
     &                 MPLET2,
     &                 MSPROJ2,
     &                 NACTE2,
     &                 TRA2,
     &                 NTRA,
     &                 NISH,
     &                 NASH,
     &                 NOSH,
     &                 NSYM,
     &                 6,
     &                 JSTATE,
     &                 job2,
     &                 jst
     &                )
#endif
        end if
      end do

C Loop over the states of JOBIPH nr JOB2
      job2_loop: DO JST=1,NSTAT(JOB2)
        JSTATE=ISTAT(JOB2)-1+JST
C Loop over the states of JOBIPH nr JOB1
        job1_loop: DO IST=1,NSTAT(JOB1)
          ISTATE=ISTAT(JOB1)-1+IST
        IF(ISTATE.LT.JSTATE) cycle
CC-----------------------------------------------------------------------

C Entry into monitor: Status line
        WRITE(STLNE2,'(A33,I5,A5,I5)')
     &      'Trans. dens. matrices for states ',ISTATE,' and ',JSTATE
        Call StatusLine('RASSI: ',STLNE2)

C Read ISTATE WF from TOTDET1 and JSTATE WF from TOTDET2
#ifdef _DMRG_
      if(.not.doDMRG)then
#endif
      DET1=>DETTOT1(1:NDET1,IST)
      DET2=>DETTOT2(1:NDET2,JST)
#ifdef _DMRG_
      end if
#endif

       if(doGSOR) then
         if(JOB1.ne.JOB2) then
           Dot_prod = 0
           Dot_prod = DDOT_(NCONF2,CI1,1,CI2,1)
           Call DAXPY_(NCONF2,Dot_prod,CI2_o,1,THETA1,1)
         end if
       end if

C Calculate whatever type of GTDM that was requested, unless
C it is known to be zero.
      HZERO = Zero
      HONE = Zero
      HTWO = Zero

      SIJ = Zero
      DYSAMP = Zero
! +++ J. Norell 12/7 - 2018
C +++ Modified by Bruno Tenorio, 2020
C Dyson amplitudes:
C DYSAMP = D_ij for states i and j
C DYSCOF = Active orbital coefficents of the DO
      IF ((IF10.or.IF01).and.DYSO) THEN
        CALL DYSON(FSBTAB1,FSBTAB2,SSTAB,
     &            DET1,DET2,
     &            IF10,IF01,
     &            DYSAMP,DYSCOF,OrbTab)

C Write Dyson orbital coefficients in AO basis to disk.
C In full biorthonormal basis:
         CALL MKDYSAB(DYSCOF,DYSAB)
C Correct Dyson norms, for a biorth. basis. Add by Bruno
         DYNORM = Zero
         CALL DYSNORM(CMO2,DYSAB,DYNORM) !do not change CMO2
       IF (DYNORM.GT.1.0D-5) THEN
C In AO basis:
         CALL MKDYSZZ(CMO2,DYSAB,DYSZZ)  !do not change CMO2
        IF (DYSO) THEN
          SFDYS(:,JSTATE,ISTATE)=DYSZZ(:)
          SFDYS(:,ISTATE,JSTATE)=DYSZZ(:)
        END IF
        DYSZZ(:) = Zero
C DYSAMPS corresponds to the Dyson norms corrected
C for a MO biorth. basis
         DYSAMPS(ISTATE,JSTATE)=SQRT(DYNORM)
         DYSAMPS(JSTATE,ISTATE)=SQRT(DYNORM)
       END IF ! AMP THRS
      END IF ! IF01 IF10

C ------------------------------------------------------------
C This part computes the needed densities for Auger.
C (DOI:10.1021/acs.jctc.2c00252)
      IF ((IF21.or.IF12).and.TDYS.and.DYSO) THEN
       Call mma_allocate(RT2M,nRT2M,Label='RT2M')
       RT2M(:) = Zero
       Call mma_allocate(RT2MAB,nRT2MAB,Label='RT2MAB')
       RT2MAB(:) = Zero
C     Defining the Binding energy Ei-Ej
       BEi=HAM(ISTATE,ISTATE)
       BEj=HAM(JSTATE,JSTATE)
       BEij=ABS(BEi-BEj)*auToEV

       IF((MPLET1-MPLET2).eq.INT(1)) THEN
C evaluate K-2V spin+1 density
        AUGSPIN=1
        CALL MKRTDM2(FSBTAB1,FSBTAB2,
     &               SSTAB,OMAP,DET1,DET2,
     &               IF21,IF12,NRT2M,RT2M,AUGSPIN,OrbTab)
        CALL RTDM2_PRINT(ISTATE,JSTATE,BEij,NDYSAB,DYSAB,NRT2MAB,
     &                  RT2M,CMO1,CMO2,AUGSPIN)

        ELSE IF ((MPLET1-MPLET2).eq.INT(-1)) THEN
C evaluate K-2V spin-1 density
        AUGSPIN=-1
        CALL MKRTDM2(FSBTAB1,FSBTAB2,
     &               SSTAB,OMAP,DET1,DET2,
     &               IF21,IF12,NRT2M,RT2M,AUGSPIN,OrbTab)
        CALL RTDM2_PRINT(ISTATE,JSTATE,BEij,NDYSAB,DYSAB,NRT2MAB,
     &                  RT2M,CMO1,CMO2,AUGSPIN)
        ELSE ! write then both
        AUGSPIN=1
        CALL MKRTDM2(FSBTAB1,FSBTAB2,
     &               SSTAB,OMAP,DET1,DET2,
     &               IF21,IF12,NRT2M,RT2M,AUGSPIN,OrbTab)
        CALL RTDM2_PRINT(ISTATE,JSTATE,BEij,NDYSAB,DYSAB,NRT2MAB,
     &                  RT2M,CMO1,CMO2,AUGSPIN)

        AUGSPIN=-1
        CALL MKRTDM2(FSBTAB1,FSBTAB2,
     &               SSTAB,OMAP,DET1,DET2,
     &               IF21,IF12,NRT2M,RT2M,AUGSPIN,OrbTab)
        CALL RTDM2_PRINT(ISTATE,JSTATE,BEij,NDYSAB,DYSAB,NRT2MAB,
     &                  RT2M,CMO1,CMO2,AUGSPIN)
       END IF
       Call mma_deallocate(RT2M)
       Call mma_deallocate(RT2MAB)
      END IF
C ------------------------------------------------------------

C evaluation of DCH shake-up intensities (DOI:10.1063/5.0062130)
      IF ((IF20.or.IF02).and.DCHS) THEN
      DCHIJ=DCHO+NASHT*(DCHO-1)
C     Defining the Binding energy Ei-Ej
      BEi=HAM(ISTATE,ISTATE)
      BEj=HAM(JSTATE,JSTATE)
      BEij=ABS(BEi-BEj)*auToEV
      Call mma_allocate(DCHSM,nDCHSM,Label='DCHSM')
      DCHSM(:) = Zero
      CALL MKDCHS(FSBTAB1,FSBTAB2,
     &            SSTAB,OMAP,DET1,DET2,
     &            IF20,IF02,NDCHSM,DCHSM,OrbTab)
      Write(6,'(A,I5,I5,A,F14.5,ES23.14)') '  RASSI Pair States:',
     &      JSTATE,ISTATE,'  ssDCH BE(eV) and Norm:  ',BEij,
     &      DCHSM(DCHIJ)
      Call mma_deallocate(DCHSM)
      END IF
* ------------------------------------------------------------

C General 1-particle transition density matrix:
      IF (IF11) THEN
        CALL MKTDM1(LSYM1,MPLET1,MSPROJ1,FSBTAB1,
     &              LSYM2,MPLET2,MSPROJ2,FSBTAB2,SSTAB,
     &            OMAP,DET1,DET2,SIJ,NASHT,
     &            TRAD,TRASD,WERD,ISTATE,
     &            JSTATE,job1,job2,ist,jst,OrbTab)
C Calculate Natural Transition Orbital (NTO):
        IF (IFNTO) THEN
         IF (job1.ne.job2) THEN
           DoNTO=.true.
         Else
           DoNTO=.false.
         End If
         IF (DoNTO) Then
          Call NTOCalc(job1,job2,ISTATE,JSTATE,TRAD,TRASD,MPLET1)
          write(6,*) 'ntocalculation finished'
         End If
        End If
C End of Calculating NTO

        IF(IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
C Compute 1-electron contribution to Hamiltonian matrix element:
        HONE=DDOT_(NTRAD,TRAD,1,FMO,1)
        END IF

C BEGIN MODIFIED by Aquilante, Segatta and Kaiser (2022)
        if (DoCoul) then
          call EXCTDM(SIJ, TRAD, TDMAB, iRC, CMO1, CMO2, TDMZZ,
     &                TRASD, TSDMAB, TSDMZZ, ISTATE, JSTATE)
        end if
C END MODIFIED by Aquilante, Segatta and Kaiser(2022)



C             Write density 1-matrices in AO basis to disk.
            IF(NATO.OR.(NPROP.GT.0))THEN

              iEmpty=0
              !> regular-TDM
              CALL MKTDAB(SIJ,TRAD,TDMAB,iRC)
              !> transform to AO basis
              CALL MKTDZZ(CMO1,CMO2,TDMAB,TDMZZ,iRC)
              If (iRC.eq.1) iEmpty=1

              !> spin-TDM
              CALL MKTDAB(Zero,TRASD,TSDMAB,iRC)
              !> transform to AO basis
              CALL MKTDZZ(CMO1,CMO2,TSDMAB,TSDMZZ,iRC)
              If (iRC.eq.1) iEmpty=iEmpty+2

              !> WE-reduced TDM''s of triplet type:
              CALL MKTDAB(Zero,WERD,WDMAB,iRC)
              !> transform to AO basis
              CALL MKTDZZ(CMO1,CMO2,WDMAB,WDMZZ,iRC)
              If (iRC.eq.1) iEmpty=iEmpty+4

              if(.not.mstate_dens)then

                IF(SaveDens) THEN
*C Transition density matrices, TDMZZ, in AO or MO basis.
*C WDMZZ similar, but WE-reduced 'triplet' densities.
                  ij=ISTATE*(iSTATE-1)/2 + JSTATE
                  jDisk_TDM(1,ij)=IDISK
                  jDisk_TDM(2,ij)=iEmpty
                  iOpt=1
                  iGo=7
                  If (AO_Mode) Then
                     CALL dens2file(TDMZZ,TSDMZZ,WDMZZ,nTDMZZ,LUTDM,
     &                              IDISK,iEmpty,iOpt,iGo,IState,jState)
                  Else
                     iEmpty=0
                     TRAD(nTrad+1)=SIJ
                     iRC=0
                     If(DDot_(nTRAD+1,TRAD,1,TRAD,1) > Zero) iRC=1
                     If (iRC.eq.1) iEmpty=1
*
                     TRASD(nTrad+1) = Zero
                     iRC=0
                     If(DDot_(nTRAD+1,TRASD,1,TRASD,1) > Zero) iRC=1
                     If (iRC.eq.1) iEmpty=iEmpty+2
*
                     WERD(nTrad+1) = Zero
                     iRC=0
                     If(DDot_(nTRAD+1,WERD,1,WERD,1) > Zero) iRC=1
                     If (iRC.eq.1) iEmpty=iEmpty+4
*
                     CALL dens2file(TRAD,TRASD,WERD,nTRAD+1,LUTDM,
     &                              IDISK,iEmpty,iOpt,iGo,ISTATE,JSTATE)
                  End If
                END IF
                !> calculate property matrix elements
                CALL PROPER(PROP,ISTATE,JSTATE,TDMZZ,WDMZZ)
              else

!               > scale rdm elements with eigenvector coefficients of Heff of a multi-state (PT2) Hamiltonian
!               > accumulate data first and run PROPER and other utility routines later
                do i = 1, nstat(job1)
                  do j = 1, nstat(job2)

                    if(i < j) cycle

                    if(qdpt2sc)then
                      fac1 = Heff_evc(job1)%sc(ist,i)
                      fac2 = Heff_evc(job2)%sc(jst,j)
                    else
                      fac1 = Heff_evc(job1)%pc(ist,i)
                      fac2 = Heff_evc(job2)%pc(jst,j)
                    end if

                    !> regular-TDM
                    call daxpy_(ntdmzz,
     &                          fac1*fac2,
     &                          tdmzz,1,
     &                          mstate_1pdens(i,j)%rtdm,1
     &                         )
                    !> spin-TDM
                    call daxpy_(ntdmzz,
     &                          fac1*fac2,
     &                          tsdmzz,1,
     &                          mstate_1pdens(i,j)%stdm,1
     &                         )
                    !> WE-reduced TDM''s of triplet type:
                    call daxpy_(ntdmzz,
     &                          fac1*fac2,
     &                          wdmzz,1,
     &                          mstate_1pdens(i,j)%wtdm,1
     &                         )
                    !> overlap
                    mstate_1pdens(i,j)%overlap=
     &                mstate_1pdens(i,j)%overlap+fac1*fac2*sij
                  end do
                end do
              end if
            END IF

          ELSE ! IF11

            !> overlap
            IF (IF00) THEN
#ifdef _DMRG_
              if(.not.doDMRG)then
#endif
                SIJ=OVERLAP_RASSI(FSBTAB1,FSBTAB2,DET1,DET2)
#ifdef _DMRG_
              else
                sij = qcmaquis_mpssi_overlap(
     &            qcm_prefixes(job1),
     &            ist,
     &            qcm_prefixes(job2),
     &            jst,
     &            .true.)
              end if !doDMRG
#endif
            END IF ! IF00
          END IF ! IF11

          if(.not.mstate_dens)then
            OVLP(ISTATE,JSTATE)=SIJ
            OVLP(JSTATE,ISTATE)=SIJ
          end if

          !> General 2-particle transition density matrix:
          IF (IF22) THEN
            CALL MKTDM2(LSYM1,MPLET1,MSPROJ1,FSBTAB1,
     &                  LSYM2,MPLET2,MSPROJ2,FSBTAB2,
     &                  SSTAB,OMAP,
     &                  DET1,DET2,NTDM2,TDM2,
     &                  ISTATE,JSTATE,OrbTab)

!           > Compute 2-electron contribution to Hamiltonian matrix element:
            IF(IFTWO.AND.(MPLET1.EQ.MPLET2))
     &      HTWO=DDOT_(NTDM2,TDM2,1,TUVX,1)

          END IF ! IF22

          !> PAM 2011 Nov 3, writing transition matrices if requested
          IF ((IFTRD1.or.IFTRD2).and..not.mstate_dens) THEN
            call trd_print(ISTATE, JSTATE, IFTRD2.AND.IF22,
     &                    TDMAB,TDM2,CMO1,CMO2,SIJ)
          END IF

          !Store SIJ temporarily
          IF (IFEJOB.and.(ISTATE.ne.JSTATE)) THEN
            HAM(ISTATE,JSTATE) = SIJ
            HAM(JSTATE,ISTATE) = SIJ
          END IF
          IF(IFHAM.AND..NOT.(IFHEXT.or.IFHEFF.or.IFEJOB))THEN
            HZERO              = ECORE*SIJ
            HIJ                = HZERO+HONE+HTWO
            HAM(ISTATE,JSTATE) = HIJ
            HAM(JSTATE,ISTATE) = HIJ

            !SI-PDFT related code for "second_time" case
            if(second_time) then
              Energies(:) = Zero
              CALL DANAME(LUIPH,'JOBGS')
              IAD = 0
              Call IDAFILE(LUIPH,2,ITOC15,30,IAD)
              IAD=ITOC15(6)
              Call DDAFILE(LUIPH,2,Energies,NSTAT(JOB1),IAD)
              do i=1,NSTAT(JOB1)
                HAM(i,i) = Energies(i)
              end do
              Call DACLOS(LUIPH)
            end if

            IF(IPGLOB.GE.4) THEN
              WRITE(6,'(1x,a,2I5)')' ISTATE, JSTATE:',ISTATE,JSTATE
              WRITE(6,'(1x,a,f16.8)')' HZERO=',HZERO
              WRITE(6,'(1x,a,f16.8)')' HONE =',HONE
              WRITE(6,'(1x,a,f16.8)')' HTWO =',HTWO
              WRITE(6,'(1x,a,f16.8)')' HIJ  =',HIJ
            END IF
          END IF
        END DO job1_loop

      END DO job2_loop
*
** For ejob, create an approximate off-diagonal based on the overlap (temporarily stored in HIJ)
*
      IF (IFEJOB) THEN
        DO JST=1,NSTAT(JOB2)
          JSTATE=ISTAT(JOB2)-1+JST
          DO IST=1,NSTAT(JOB1)
            ISTATE=ISTAT(JOB1)-1+IST
            IF(ISTATE.LE.JSTATE) CYCLE
            SIJ=HAM(ISTATE,JSTATE)
            HII=HAM(ISTATE,ISTATE)
            HJJ=HAM(JSTATE,JSTATE)
            HAM(ISTATE,JSTATE) = SIJ*(HII+HJJ)*Half
            HAM(JSTATE,ISTATE) = SIJ*(HII+HJJ)*Half
          END DO
        END DO
      END IF

      IF(DoGSOR) then
        if(job1.ne.job2) then
        Norm_Fac = Zero
        dot_prod = DDOT_(NCONF2,THETA1,1,THETA1,1)
        Norm_Fac = One/sqrt(dot_prod)
        Call DSCAL_(NCONF2,Norm_Fac,THETA1,1)

      !Write theta1 to file.
        LUCITH=87
        LUCITH=IsFreeUnit(LUCITH)
        !Open(unit=87,file='CI_THETA', action='write',iostat=ios)
        Call Molcas_Open(LUCITH,'CI_THETA')
        do i=1,NCONF2
          write(LUCITH,*) Theta1(i)
        end do
        Close(LUCITH)

       !Now we need to build the other states.
        call mma_allocate(detcoeff2,nDet2,label='detcoeff2')
        DO JST=2,NSTAT(JOB2)
          JSTATE=ISTAT(JOB2)-1+JST
          CALL READCI(JSTATE,SGS(2),CIS(2),NCONF2,CI2)
          Call DCOPY_(NCONF2,CI2,1,CI2_o,1)
          CALL DCOPY_(NDET2,[Zero],0,DET2,1)
          If (TrOrb) CALL CITRA (WFTP2,SGS(2),CIS(2),EXS(2),
     &                           LSYM2,TRA2,NCONF2,CI2)
          CALL PREPSD(WFTP2,SGS(2),CIS(2),LSYM2,
     &                CNFTAB2,SPNTAB2,
     &                SSTAB,FSBTAB2,NCONF2,CI2,
     &                DET2,detocc,detcoeff2,
     &                TRANS2)

          CALL mma_allocate(ThetaN,NCONF2,Label='ThetaN')
          ThetaN(:)=0.0D0
          Norm_Fac = DDOT_(NCONF2,THETA1,1,CI2_o,1)
          Call DAXPY_(NCONF2,-Norm_Fac,THETA1,1,ThetaN,1)

          LUCITH=IsFreeUnit(LUCITH)
          Call Molcas_Open(LUCITH,'CI_THETA')
          !Open(unit=87,file='CI_THETA', action='read',iostat=ios)
          if(JST-1.ge.2) then
            do i=1,NCONF2
              Read(LUCITH,*) dot_prod ! dummy
            end do
          end if
          CALL mma_allocate(ThetaM,NCONF2,Label='ThetaM')
          DO IST=2,JST-1
            ThetaM(:)=0.0D0
            !Read in previous theta vectors
            do i=1,NCONF2
              Read(LUCITH,*) ThetaM(i)
            end do
            Dot_prod = DDOT_(NCONF2,ThetaM,1,CI2_o,1)
           Call DAXPY_(NCONF2,-Dot_prod,ThetaM,1,ThetaN,1)

          END DO
          call mma_deallocate(detcoeff2)
          Close(LUCITH)
          !Normalize
          dot_prod = DDOT_(NCONF2,ThetaN,1,ThetaN,1)
          Norm_Fac = One/sqrt(dot_prod)
          Call DSCAL_(NCONF2,Norm_Fac,ThetaN,1)

        !dot_prod = DDOT_(NCONF2,THETA1,1,THETA1,1)
        !dot_prod = DDOT_(NCONF2,ThetaN,1,THETA1,1)
        !dot_prod = DDOT_(NCONF2,ThetaN,1,ThetaN,1)

          !Write to file
          LUCITH=IsFreeUnit(LUCITH)
          Call Molcas_Open(LUCITH,'CI_THETA')
          Call Append_file(LUCITH)
          !Open(unit=87,file='CI_THETA', position='append',iostat=ios,
!    &    action='write')
          do i=1,nConf2
            write(LUCITH,*) ThetaN(i)
          end do
          close(LUCITH)
          !Deallocate
          CALL mma_deallocate(ThetaN)
        END DO
!Copy to new IPH file
        LUCITH=IsFreeUnit(LUCITH)
        Call Molcas_Open(LUCITH,'CI_THETA')
!       Open(unit=87,file='CI_THETA',iostat=ios,
!    &    action='read')
        CALL DANAME(LUIPHn,'JOBGS')
        IAD = 0
        Call IDAFILE(LUIPHn,2,ITOC15,30,IAD)
        IAD=ITOC15(4)
        do i=1,ISTAT(JOB1)-1
         ThetaM(:)=0.0D0
         do j=1,nCONF2
           read(LUCITH,*) ThetaM(i)
         end do
         Call DDafile(LUIPHn,1,ThetaM,nCONF2,IAD)
        end do

       IAD = ITOC15(4)
       ThetaM(:)=0.0D0
       Call DDAFILE(LUIPHn,2,ThetaM,nCONF2,IAD)
       Call DDAFILE(LUIPHn,2,ThetaM,nCONF2,IAD)

       Close(LUCITH)
       Call DACLOS(LUIPHn)
       CALL mma_deallocate(ThetaM)
       end if
       CALL mma_deallocate(Theta1)
      end if!DoGSOR


#ifdef _DMRG_
      IF(IPGLOB.GE.4) THEN
         write(6,*) 'full SF-HAMILTONIAN '
         write(6,*) 'dimension: ',nstate**2
         call pretty_print_util(HAM,1,nstate,1,nstate,
     &                          nstate,nstate,1,6)
      END IF
#endif

!     > create actual property data and put everything to file (if requested) in case of using eigenvectors of a multi-state (PT2) Hamiltonian
      if(mstate_dens)then
        DO JST=1,NSTAT(JOB2)
          JSTATE=ISTAT(JOB2)-1+JST
          DO IST=1,NSTAT(JOB1)
            ISTATE=ISTAT(JOB1)-1+IST
            if(istate < jstate) cycle

            ovlp(istate,jstate) = mstate_1pdens(ist,jst)%overlap
            ovlp(jstate,istate) = mstate_1pdens(ist,jst)%overlap

            call prpdata_mspt2_eigenvectors(
     &                                      mstate_1pdens(ist,jst)%rtdm,
     &                                      mstate_1pdens(ist,jst)%stdm,
     &                                      mstate_1pdens(ist,jst)%wtdm,
     &                                      prop,
     &                                      nprop,
     &                                      nstate,
     &                                      istate,
     &                                      jstate,
     &                                      ntdmzz,
     &                                      iDisk_TDM(JSTATE,ISTATE,1),
     &                                      iDisk_TDM(JSTATE,ISTATE,2),
     &                                      lutdm,
     &                                      (sonatnstate.gt.0),
     &                                      if11.and.(lsym1.eq.lsym2)
     &                                     )
          end do
        end do
      end if

      IF(WFTP1.EQ.'GENERAL ') THEN
        if(.not.doDMRG)then
          CALL MkGUGA_Free(SGS(1),CIS(1),EXS(1))
        end if
      END IF
      IF(WFTP2.EQ.'GENERAL ') THEN
        if(.not.doDMRG)then
          CALL MkGUGA_Free(SGS(2),CIS(2),EXS(2))
        end if
      END IF

      IF(JOB1.NE.JOB2) THEN
        Call mma_deallocate(TRA1)
        Call mma_deallocate(TRA2)
      END IF
      nullify(DET1,DET2)
      call mma_deallocate(DETTOT1)
      call mma_deallocate(DETTOT2)
      call mma_deallocate(detocc)
      CALL mma_deallocate(CI2)
      If (DoGSOR) CALL mma_deallocate(CI2_o)
      CALL mma_deallocate(CI1)
      if(.not.doDMRG)then
        Call mma_deallocate(TRANS2)
        Call mma_deallocate(TRANS1)
        Call mma_deallocate(SPNTAB1)
        Call mma_deallocate(SPNTAB2)
      end if
      IF ((IF10.or.IF01).and.DYSO) THEN
        Call mma_deallocate(DYSCOF)
        Call mma_deallocate(DYSAB)
        Call mma_deallocate(DYSZZ)
      END IF
      IF (IF11) THEN
        Call mma_deallocate(TRAD)
        Call mma_deallocate(TRASD)
        Call mma_deallocate(WERD)
        IF(NATO.OR.NPROP.GT.0) THEN
          Call mma_deallocate(TDMAB)
          Call mma_deallocate(TSDMAB)
          Call mma_deallocate(WDMAB)
          Call mma_deallocate(TDMZZ)
          Call mma_deallocate(TSDMZZ)
          Call mma_deallocate(WDMZZ)
        END IF
      END IF
      Call mma_deallocate(TDM2)

      IF(IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
        Call mma_deallocate(FMO)
        Call mma_deallocate(TUVX)
      END IF

      Call mma_deallocate(CMO2)
      Call mma_deallocate(CMO1)
      Call mma_deallocate(PART)
      Call mma_deallocate(OrbTab)
      Call mma_deallocate(SSTAB)
      if(.not.doDMRG)then
        Call mma_deallocate(REST2)
        Call mma_deallocate(REST1)
        Call mma_deallocate(CNFTAB2)
        Call mma_deallocate(CNFTAB1)
        Call mma_deallocate(FSBTAB2)
        Call mma_deallocate(FSBTAB1)
      end if
      CALL mma_deallocate(OMAP)

      !> release memory
      if(mstate_dens)then
        do i = 1, nstat(job1)
          do j = 1, nstat(job1)
            call mma_deallocate(mstate_1pdens(i,j)%rtdm,safe='*')
            call mma_deallocate(mstate_1pdens(i,j)%stdm,safe='*')
            call mma_deallocate(mstate_1pdens(i,j)%wtdm,safe='*')
          end do
        end do
        if(allocated(mstate_1pdens)) deallocate(mstate_1pdens)
      end if

#ifdef _TIME_GTDM_
      Call CWTime(TCpu2,TWall2)
      write(6,*) 'Time for GTDM : ',TCpu2-TCpu1,TWall2-TWall1
#endif

      END SUBROUTINE GTDMCTL
