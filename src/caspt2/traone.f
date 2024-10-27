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
      SUBROUTINE TRAONE(CMO,NCMO)
      use OneDat, only: sNoNuc, sNoOri
      use caspt2_output, only:iPrGlb
      use caspt2_data, only: HONE
      use caspt2_data, only: LUONEM
      use PrintLevel, only: verbose
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
#include "warnings.h"
#include "caspt2.fh"
      INTEGER NCMO
      REAL*8 CMO(NCMO)

      INTEGER nBasXX(8),Keep(8)
      Logical iSquar, Found
      character(len=8) :: Label
      Real*8, allocatable:: WFLT(:), Temp(:), WDLT(:), WDSQ(:),
     &                      WFMO(:), WTMP(:)

c Objective: Transformation of one-electron integrals
c (effective one electron Hamiltonian) for CASPT2.

#ifdef _DEBUGPRINT_
      IFTEST=1
#else
      IFTEST=0
#endif

      Call GetOrd(IRC,iSquar,nSymXX,nBasXX,Keep)
      IF ( IPRGLB.GE.VERBOSE ) THEN
        If(iSquar)      WRITE(6,*) 'TRAONE OrdInt status: squared'
        If(.not.iSquar) WRITE(6,*) 'TRAONE OrdInt status: non-squared'
      ENDIF
      IERR=0
      DO ISYM=1,NSYM
        IF (NBAS(ISYM).NE.NBASXX(ISYM)) IERR=1
      END DO
      IF(IERR.NE.0) THEN
        WRITE(6,*)'     *** ERROR IN SUBROUTINE TRAONE ***'
        WRITE(6,*)'          INCOMPATIBLE BASIS DATA'
        WRITE(6,*)
        WRITE(6,*)' JOBIPH NR OF SYMM:', NSYM
        WRITE(6,*)' JOBIPH NR OF BASIS FUNCTIONS/SYMM:'
        WRITE(6,'(1x,8I5)')(NBAS(I),I=1,NSYM)
        WRITE(6,*)
        WRITE(6,*)' ORDINT NR OF SYMM:', NSYMXX
        WRITE(6,*)' ORDINT NR OF BASIS FUNCTIONS/SYMM:'
        WRITE(6,'(1x,8I5)')(NBASXX(I),I=1,NSYMXX)
        CALL ABEND()
      END IF
c Allocate FLT,DLT, and DSQ.
      CALL mma_allocate(WFLT,NBTRI,Label='WFLT')
c Read nuclear repulsion energy:
      IRC=-1
      IOPT=0
      ICOMP=0
      ISYLBL=1
      IF ( IFTEST.NE.0 ) WRITE(6,*)' GET POTNUC FROM RUNFILE'
      Call Get_dScalar('PotNuc',PotNuc)
      IF ( IFTEST.NE.0 ) WRITE(6,*)' POTNUC:',POTNUC
c Read one-electron hamiltonian matrix into FLT.
      IRC=-1
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICOMP=1
      ISYLBL=1
      Label='OneHam'
      IF ( IFTEST.NE.0 ) WRITE(6,*)' CALLING RDONE (ONEHAM)'
      CALL RDONE(IRC,IOPT,Label,ICOMP,WFLT,ISYLBL)
      IF ( IFTEST.NE.0 ) WRITE(6,*)' BACK FROM RDONE'
      IF(IRC.NE.0) THEN
        WRITE(6,*)'TRAONE Error: RDONE failed reading OneHam.'
        Call Quit(_RC_IO_ERROR_READ_)
      END IF

      IF ( IFTEST.NE.0 ) THEN
        WRITE(6,*)'     TEST PRINTS FROM TRAONE.'
        WRITE(6,*)'     NAKED 1-EL HAMILTONIAN IN AO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          IF ( NBAS(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',WFLT(ISTLT),NBAS(ISYM))
            ISTLT=ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
          END IF
        END DO
      END IF

c If this is a perturbative reaction field calculation then
c modifiy the one-electron Hamiltonian by the reaction field and
c the nuclear attraction by the cavity self-energy

      If ( RFpert ) then
         nTemp=0
         Do iSym=1,nSym
            nTemp=nTemp+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
         Call mma_allocate(Temp,nTemp,Label='Temp')
*
         Call f_Inquire('RUNOLD',Found)
         If (Found) Call NameRun('RUNOLD')
         Call Get_dScalar('RF Self Energy',ERFSelf)
         Call Get_dArray('Reaction field',Temp,nTemp)
         If (Found) Call NameRun('#Pop')
         PotNuc=PotNuc+ERFself
         Call Daxpy_(nTemp,1.0D0,Temp,1,WFLT,1)
*
         Call mma_deallocate(Temp)
         IF ( IFTEST.NE.0 ) THEN
           WRITE(6,*)' 1-EL HAMILTONIAN INCLUDING REACTION FIELD'
           ISTLT=1
           DO ISYM=1,NSYM
             IF ( NBAS(ISYM).GT.0 ) THEN
               WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
               CALL TRIPRT(' ',' ',WFLT(ISTLT),NBAS(ISYM))
               ISTLT=ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
             END IF
           END DO
         END IF
      End If

      EONE=0.0d0
      ETWO=0.0d0
c The following section is needed for frozen orbitals:
      IF(NFROT.EQ.0) GOTO 300
      CALL mma_allocate(WDLT,NBTRI,LABEL='WDLT')
      CALL mma_allocate(WDSQ,NBSQT,LABEL='WDSQ')
c Compute the density matrix of the frozen orbitals
c The DLT matrix contains the same data as DSQ, but
c with symmetry blocks in lower triangular format, and
c with non-diagonal elements doubled.
      WDLT(:)=0.0D0
      WDSQ(:)=0.0D0
      ISTMO=1
      ISTSQ=1
      ISTLT=1
      DO 100 ISYM=1,NSYM
        NF=NFRO(ISYM)
        NB=NBAS(ISYM)
        IF(NB.EQ.0) GOTO 100
        IF(NF.EQ.0) GOTO 110
        CALL DGEMM_('N','T',NB,NB,NF,2.0D0,CMO(ISTMO),NB,
     &             CMO(ISTMO),NB,0.0D0,WDSQ(ISTSQ),NB)
        IJ=ISTLT-1
        DO 130 IB=1,NB
          DO 140 JB=1,IB
            IJ=IJ+1
            WDLT(IJ)=2.0D0*WDSQ(ISTSQ+JB-1+(IB-1)*NB)
140       CONTINUE
          WDLT(IJ)=0.5D0*WDLT(IJ)
130     CONTINUE
110     CONTINUE
        ISTMO=ISTMO+NB*NB
        ISTSQ=ISTSQ+NB*NB
        ISTLT=ISTLT+NB*(NB+1)/2
100   CONTINUE

c One-electron contribution to the core energy.
c Note that FLT still contains only the naked
c  one-electron hamiltonian.
      EONE=DDOT_(NBTRI,WDLT,1,WFLT,1)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate Fock-matrix for frozen orbitals
*     and compute the total core energy
*     Look out-- we temporarily allocate all available memory.
*
      ExFac=1.0D0
         Call FTwo_Drv(nSym,nBas,nFro,KEEP,
     &                 WDLT,WDSQ,WFLT,NBTRI,
     &                 ExFac,nBMX,CMO)

*                                                                      *
************************************************************************
*                                                                      *
c Compute the two-electron contribution to the core energy
      ETWO=0.5D0*(DDOT_(NBTRI,WDLT,1,WFLT,1)-EONE)
      CALL mma_deallocate(WDSQ)
      CALL mma_deallocate(WDLT)
c Previous section was bypassed if NFROT.EQ.0.
 300  CONTINUE
      ECORE=POTNUC+EONE+ETWO
      IF ( IFTEST.NE.0 ) THEN
         WRITE(6,'(6X,A,ES20.10)') 'NUCLEAR REPULSION ENERGY:',POTNUC
         WRITE(6,'(6X,A,ES20.10)') 'ONE-ELECTRON CORE ENERGY:',EONE
         WRITE(6,'(6X,A,ES20.10)') 'TWO-ELECTRON CORE ENERGY:',ETWO
         WRITE(6,'(6X,A,ES20.10)') '       TOTAL CORE ENERGY:',ECORE
      ENDIF

c Allocate FMO, TMP:
      NWTMP=2*NBMX**2
      CALL mma_allocate(WFMO,notri,LABEL='WFMO')
      CALL mma_allocate(WTMP,NWTMP,LABEL='WTMP')

c Transform one-electron effective Hamiltonian:
      WFMO(:)=0.0D0
      WTMP(:)=0.0D0
      ICMO=1
      IAO =1
      IMO =1
      DO 200 ISYM=1,NSYM
         ICMO=ICMO+NBAS(ISYM)*NFRO(ISYM)
         IOFF=1+NBAS(ISYM)*NBAS(ISYM)
         IF(NORB(ISYM).GT.0) THEN
           CALL SQUARE(WFLT(IAO),WTMP,1,NBAS(ISYM),NBAS(ISYM))

           CALL DGEMM_('T','N',NORB(ISYM),NBAS(ISYM),NBAS(ISYM),
     &                  1.0d0,CMO(ICMO),NBAS(ISYM),WTMP,
     &                  NBAS(ISYM),0.0d0,WTMP(IOFF),NORB(ISYM))

           Call DGEMM_Tri('N','N',NORB(ISYM),NORB(ISYM),NBAS(ISYM),
     &                    1.0D0,WTMP(IOFF),NORB(ISYM),
     &                          CMO(ICMO),NBAS(ISYM),
     &                    0.0D0,WFMO(IMO),NORB(ISYM))
         END IF
         ICMO=ICMO+NBAS(ISYM)*(NORB(ISYM)+NDEL(ISYM))
         IAO =IAO +NBAS(ISYM)*(NBAS(ISYM)+1)/2
         IMO =IMO +NORB(ISYM)*(NORB(ISYM)+1)/2
200   CONTINUE

      IF ( IFTEST.NE.0 ) THEN
        WRITE(6,*)'      EFFECTIVE 1-EL HAMILTONIAN IN MO BASIS'
        ISTLT=1
        DO ISYM=1,NSYM
          IF ( NORB(ISYM).GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')' SYMMETRY SPECIES:',ISYM
            CALL TRIPRT(' ',' ',WFMO(ISTLT),NORB(ISYM))
            ISTLT=ISTLT+NORB(ISYM)*(NORB(ISYM)+1)/2
          END IF
        END DO
      END IF
      IDISK=IEOF1M
      IAD1M(3)=IDISK
      CALL DDAFILE(LUONEM,1,WFMO,notri,IDISK)
      IEOF1M=IDISK
      CALL DCOPY_(NOTRI,WFMO,1,HONE,1)
      CALL mma_deallocate(WTMP)
      CALL mma_deallocate(WFMO)
      CALL mma_deallocate(WFLT)

      End SUBROUTINE TRAONE
