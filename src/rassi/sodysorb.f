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
* Copyright (C) 1987, Per Ake Malmqvist                                *
*               2018, Jesper Norell                                    *
*               2018, Joel Creutzberg                                  *
*               2023, Ignacio Fdez. Galvan                             *
************************************************************************
      SUBROUTINE SODYSORB(NSS,USOR,USOI,DYSAMPS,NZ,SOENE)
      use rassi_global_arrays, only: SFDYS, SODYSAMPS,
     &                               SODYSAMPSR, SODYSAMPSI, JBNUM
      use OneDat, only: sNoNuc, sNoOri
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NSTATE, DYSEXPORT, DYSEXPSO, MLTPLT
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: NBASF,NOSH

      IMPLICIT None
      Integer NSS, NZ
      REAL*8 USOR(NSS,NSS), USOI(NSS,NSS)
      ! Array for calculation of amplitudes
      Real*8 DYSAMPS(NSTATE,NSTATE)
      ! Arrays for orbital export
      Real*8 SOENE(NSS)

      ! Arrays, bounds, and indices
      REAL*8    MSPROJI,MSPROJJ, CJR, CJI, CIR, CII, CREAL, CIMAG,
     &          AMPLITUDE, AMPR, AMPI
      INTEGER   SOTOT,SFTOT
      INTEGER   ORBNUM
      INTEGER   SODYSCIND
      INTEGER   SFI,SFJ,ZI,ZJ,NSZZ,NDUM
      INTEGER   ISTATE, JOB1, MPLET1, MSPROJ, JSTATE, NSSQ, NPROD, ISY,
     &          NO, NB, IRC, IOPT, ICMP, ISYLAB, NOFF, JEIG, IEIG, LUNIT
      INTEGER, EXTERNAL:: IsFreeUnit

      ! Arrays for calculation of amplitudes
      Real*8 SODYSCOFSR(NZ),SODYSCOFSI(NZ)
      Real*8 SZZFULL(NZ,NZ)

      ! Arrays for orbital export
      Real*8 SODYSCMOR(NZ*NSS)
      Real*8 SODYSCMOI(NZ*NSS)
      Real*8 DYSEN(NSS)
      Real*8 AMPS(NSS)
      Character(LEN=30) Filename
      Character(LEN=80) TITLE
      character(len=8) :: Label


      Integer, Allocatable:: SO2SF(:)
      Real*8, Allocatable:: MSPROJS(:), SZZ(:)

! +++ J.Norell 2018

C Calculates spin-orbit Dyson orbitals
C The routine was in some part adapted from DO_SONATORB

C Computes SO Dyson amplitudes by expanding the SF results with
C the SO eigenvectors (of the complex Hamiltonian)
C 1. (Fast): Compute the SO amplitudes directly from the SF amplitudes
C    (approximation) for all states
C 2. (Slower): Compute the full SO Dyson orbitals for the requested
C    initial states and export them to .molden format. SO amplitudes
C    are correctly calculated for these states.

! IFG: Added DysOrb export, not sure it's correct

! ****************************************************************

C Setup SO2SF list which contains the original SF state numbers
C as a function of the SO state number
C And MSPROJS which saves their ms projections for later use

      CALL mma_allocate(SO2SF,NSS,Label='SO2SF')
      CALL mma_allocate(MSPROJS,NSS,Label='MSPROJS')
      SOTOT=0
      SFTOT=0
      DO ISTATE=1,NSTATE
       JOB1=JBNUM(ISTATE)
       MPLET1=MLTPLT(JOB1)
       SFTOT=SFTOT+1

       DO MSPROJ=-MPLET1+1,MPLET1-1,2
        SOTOT=SOTOT+1
        SO2SF(SOTOT)=SFTOT
        MSPROJS(SOTOT)=MSPROJ

       END DO ! DO MSPROJ1=-MPLET1+1,MPLET1-1,2
      END DO ! DO ISTATE=1,NSTATE

      IF(.NOT.DYSEXPORT) THEN ! Approximative amplitude calculation

C Now construct the SF dysamps in the multiplicity expanded basis
C (initially all real, therefore put into SODYSAMPSR)
      SODYSAMPSR(:,:)=0.0D0
      SODYSAMPSI(:,:)=0.0D0
      DO JSTATE=1,NSS
         DO ISTATE=JSTATE+1,NSS
          SFJ=SO2SF(JSTATE)
          SFI=SO2SF(ISTATE)
          SODYSAMPSR(JSTATE,ISTATE)=DYSAMPS(SFJ,SFI)
          SODYSAMPSR(ISTATE,JSTATE)=DYSAMPS(SFJ,SFI)
         END DO
      END DO

C Now perform the transformation from SF dysamps to SO dysamps
C by combining the multiplicity expanded SF dysamps with the
C SO eigenvector in the ZTRNSF routine.
      CALL ZTRNSF(NSS,USOR,USOI,SODYSAMPSR,SODYSAMPSI)

C Compute the magnitude of the complex amplitudes as an approximation
      SODYSAMPSR(:,:)=SODYSAMPSR*SODYSAMPSR
      SODYSAMPSI(:,:)=SODYSAMPSI*SODYSAMPSI
      SODYSAMPS(:,:)=SQRT(SODYSAMPSR+SODYSAMPSI)

      END IF ! Approximative amplitude calculation

! ****************************************************************

      IF (.NOT.DYSEXPORT) THEN
       GOTO 100
      END IF

! Export part of the routine and exact calculation of amplitudes

! Read in the atomic overlap matrix, that will be needed below for
! for normalization of DOs
! (Code from mksxy)
      NSZZ=0
      NSSQ=0
      NPROD=0
      DO ISY=1,NSYM
        NO=NOSH(ISY)
        NB=NBASF(ISY)
        NSZZ=NSZZ+(NB*(NB+1))/2
        NSSQ=MAX(NSSQ,NB**2)
        NPROD=MAX(NPROD,NO*NB)
      END DO
      CALL mma_allocate(SZZ,NSZZ,Label='SZZ')
      IRC=-1
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      Label='MLTPL  0'
      CALL RDONE(IRC,IOPT,Label,ICMP,SZZ,ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE SODYSORB ***'
        WRITE(6,*)'     OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF

! SZZ is originally given in symmetry-blocked triangular form,
! lets make it a full matrix for convenience
      SZZFULL=0.0D0
      NDUM=1
      NOFF=0
      DO ISY=1,NSYM
       NB=NBASF(ISY)
       DO ZJ=1,NB
        DO ZI=1,ZJ
         SZZFULL(ZJ+NOFF,ZI+NOFF)=SZZ(NDUM)
         SZZFULL(ZI+NOFF,ZJ+NOFF)=SZZ(NDUM)
         NDUM=NDUM+1
        END DO
       END DO
       NOFF=NOFF+NB
      END DO
      CALL mma_deallocate(SZZ)

! ****************************************************************

C Multiply together with the SO eigenvector coefficients with the SF
C Dyson orbital coefficients in the atomic basis to obtain the full
C SO Dyson orbitals

C Multiply together with the SO eigenvector coefficients with the SF
C Dyson orbital coefficients in the atomic basis to obtain
C SO Dyson orbitals

      SODYSAMPS(:,:)=0.0D0
      ! For all requested initial states J and all final states I
      DO JSTATE=1,DYSEXPSO

         ! For each initial state JSTATE up to DYSEXPSFSO we will
         ! gather all the obtained Dysorbs
         ! and export to a shared .molden file
         SODYSCIND=0 ! Orbital coeff. index
         ORBNUM=0 ! Dysorb index for given JSTATE
         SODYSCMOR=0.0D0 ! Real orbital coefficients
         SODYSCMOI=0.0D0 ! Imaginary orbital coefficients
         DYSEN=0.0D0 ! Orbital energies
         AMPS=0.0D0 ! Transition amplitudes (shown as occupations)

         DO ISTATE=JSTATE+1,NSS

          ! Reset values for next state combination
          SODYSCOFSR=0.0D0
          SODYSCOFSI=0.0D0

          ! Iterate over the eigenvector components of both states
          DO JEIG=1,NSS

           ! Coefficient of first state
           CJR=USOR(JEIG,JSTATE)
           CJI=USOI(JEIG,JSTATE)
           ! Find the corresponding SF states
           SFJ=SO2SF(JEIG)

           DO IEIG=1,NSS

            ! Coefficient of second state
            CIR=USOR(IEIG,ISTATE)
            CII=USOI(IEIG,ISTATE)
            ! Find the corresponding SF states
            SFI=SO2SF(IEIG)

            ! Check change in ms projection
            MSPROJJ=MSPROJS(JEIG)
            MSPROJI=MSPROJS(IEIG)
            ! Check |delta ms|=0.5 selection rule
            IF(ABS(MSPROJJ-MSPROJI).NE.1) THEN
             CYCLE
            END IF

            IF (DYSAMPS(SFJ,SFI).GT.1.0D-5) THEN
             ! Multiply together coefficients
             CREAL=CJR*CIR+CJI*CII
             CIMAG=CJR*CII-CJI*CIR
             ! Multiply with the corresponding SF Dyson orbital
             SODYSCOFSR=SODYSCOFSR+CREAL*SFDYS(:,SFJ,SFI)
             SODYSCOFSI=SODYSCOFSI+CIMAG*SFDYS(:,SFJ,SFI)
            END IF

           END DO ! IEIG
          END DO ! JEIG

! Normalize the overlap of SODYSCOFS expanded orbitals with the
! atomic overlap matrix SZZ to obtain correct amplitudes

          AMPLITUDE=0.0D0

          DO ZJ=1,NZ
           DO ZI=1,NZ
            AMPR=SODYSCOFSR(ZJ)*SODYSCOFSR(ZI)
     &            +SODYSCOFSI(ZJ)*SODYSCOFSI(ZI)
            AMPI=SODYSCOFSI(ZJ)*SODYSCOFSR(ZI)
     &            -SODYSCOFSR(ZJ)*SODYSCOFSI(ZI)
            AMPLITUDE=AMPLITUDE+(AMPR+AMPI)*SZZFULL(ZJ,ZI)
           END DO ! ZI
          END DO ! ZJ

          AMPLITUDE=SQRT(AMPLITUDE)
          SODYSAMPS(JSTATE,ISTATE)=AMPLITUDE
          SODYSAMPS(ISTATE,JSTATE)=AMPLITUDE

          ! Export Re and Im part of the coefficients
          IF (AMPLITUDE.GT.1.0D-5) THEN
           DO NDUM=1,NZ
              SODYSCIND=SODYSCIND+1
              SODYSCMOR(SODYSCIND)=SODYSCOFSR(NDUM)
              SODYSCMOI(SODYSCIND)=SODYSCOFSI(NDUM)
           END DO
           ORBNUM=ORBNUM+1
           DYSEN(ORBNUM)=SOENE(ISTATE)-SOENE(JSTATE)
           AMPS(ORBNUM)=AMPLITUDE*AMPLITUDE
          END IF

        END DO ! ISTATE

! If at least one orbital was found, export it/them
        IF(ORBNUM.GT.0) THEN
         Write(filename,'(A,I0,A3)') 'MD_DYS.SO.',JSTATE,'.Re'
         Call Molden_DysOrb(filename,DYSEN,AMPS,SODYSCMOR,ORBNUM,NZ)
         Write(filename,'(A,I0,A3)') 'MD_DYS.SO.',JSTATE,'.Im'
         Call Molden_DysOrb(filename,DYSEN,AMPS,SODYSCMOI,ORBNUM,NZ)

         ! This does not work for SO-Dyson orbitals, because they may
         ! contain contributions from several irreps.
         ! Either that's a bug elsewhere in the code, or the InpOrb
         ! format is not adequate for these orbitals.
         Write(filename,'(A,I0,A3)') 'DYSORB.SO.',JSTATE,'.Re'
         LUNIT=IsFreeUnit(50)
         Write(TITLE,'(A,I0,A)') '* Spin-orbit Dyson orbitals for '//
     &                           'state ',JSTATE,' (real part)'
         Call WRVEC_DYSON(filename,LUNIT,NSYM,NBASF,ORBNUM,SODYSCMOR,
     &                    AMPS,DYSEN,Trim(TITLE),NZ)
         Write(filename,'(A,I0,A3)') 'DYSORB.SO.',JSTATE,'.Im'
         Write(TITLE,'(A,I0,A)') '* Spin-orbit Dyson orbitals for '//
     &                           'state ',JSTATE,' (imaginary part)'
         Call WRVEC_DYSON(filename,LUNIT,NSYM,NBASF,ORBNUM,SODYSCMOI,
     &                    AMPS,DYSEN,Trim(TITLE),NZ)
         Close(LUNIT)
        END IF

       END DO ! JSTATE

100    CONTINUE

! ****************************************************************

C Free all the allocated memory

      Call mma_deallocate(SO2SF)
      Call mma_deallocate(MSPROJS)

      END SUBROUTINE SODYSORB
