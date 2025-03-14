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
* Copyright (C) 1992,2001, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE CSDTMT_GAS(IPRCSF)
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS, only: DFTP, CFTP, DTOC, Z_PTDT, REO_PTDT
*
* Construct in IDTFTP list of proto type combinations in IDFTP
* Construct in ICFTP list of proto type CSF's in ICFTP
* Construct in DTOC matrix expanding proto type CSF's in terms of
* prototype combinations in DTOC
*
* Construct also array for going from lexical address of prototype det
* to address in IDFTP
*
* where MXPDBL is size of largest prototype combination block .
*
* Jeppe Olsen
*
* Changed to Combination form, June 92
* Adapted to LUCIA December 2001
*
      use lucia_data, only:MAXOP,MINOP,NPCMCNF,NPCSCNF
      use lucia_data, only:MS2,MULTS,PSSIGN
      Implicit NONE
      Integer IPRCSF
*. Output
      Real*8,  Allocatable:: SCR1(:)
      Integer, Allocatable:: iSCR2(:)
      INTEGER NTEST,MAX_DC,IOPEN,L,LSCR,IDTBS,ICSBS,ITP,IFLAG,
     &        IALPHA,IDET,ICDCBS,NNCS,NNCM,NNDET
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRCSF)
*. Size of largest csf-sd block
      MAX_DC = 0
      DO IOPEN = 0, MAXOP
        L = NPCMCNF(IOPEN+1)
        MAX_DC = MAX(MAX_DC,L)
      END DO
      IF(NTEST.GE.100) WRITE(6,*) ' Size of largest D to C block ',
     &MAX_DC
      LSCR = MAX(MAX_DC,MAXOP)
                LSCR = MAX_DC*MAXOP+MAXOP
      Call mma_allocate(SCR1,LSCR,Label='SCR1')
*
*
* .. Set up combinations and upper determinants
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' **************************************'
        WRITE(6,*) ' Generation of proto type determinants '
        WRITE(6,*) ' **************************************'
        WRITE(6,*)
      END IF
*. Still tired of stupid compiler warnings
      IDTBS = 0
      ICSBS = 0
      DO IOPEN = 0, MAXOP
        ITP = IOPEN + 1
        IF( NTEST .GE. 5 ) THEN
          WRITE(6,*)
          WRITE(6,'(A,I3,A)')
     &    '       Type with ',IOPEN,' open orbitals '
          WRITE(6,'(A)')
     &    '       **********************************'
          WRITE(6,*)
        END IF
        IF( ITP .EQ. 1 ) THEN
          IDTBS = 1
          ICSBS = 1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NPCMCNF(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NPCSCNF(ITP-1)
        END IF
C
        IF( IOPEN .NE. 0 ) THEN
*. Proto type combinations and branching diagram for
*  proto type combinations
          IF( MS2+1 .EQ. MULTS ) THEN
            IFLAG = 2
            CALL SPNCOM_LUCIA(  IOPEN,
     &                            MS2,
     &                          NNDET,
     &                        DFTP(IDTBS),
     &                        CFTP(ICSBS),
*
     &                          IFLAG,
     &                         PSSIGN,
     &                         IPRCSF)
C                SPNCOM(NOPEN,MS2,NDET,IABDET,
C    &                  IABUPP,IFLAG,PSSIGN,IPRCSF)
          ELSE
            IFLAG = 1
            CALL SPNCOM_LUCIA(  IOPEN,
     &                            MS2,
     &                          NNDET,
     &                        DFTP(IDTBS),
     &                        CFTP(ICSBS),
*
     &                          IFLAG,
     &                         PSSIGN,
     &                         IPRCSF)
            IFLAG = 3
            CALL SPNCOM_LUCIA(  IOPEN,
     &                        MULTS-1,
     &                          NNDET,
     &                        DFTP(IDTBS),
     &                        CFTP(ICSBS),
*
     &                          IFLAG,
     &                         PSSIGN,
     &                         IPRCSF)
           END IF
         END IF
      END DO
*     ^ End of loop over number of open orbitals
*
* Set up z-matrices for addressing prototype determinants with
* a given number of open orbitals, and for readdressing to
* the order given in IDCNF
*. Scr : largest block of 2*NOPEN + (NALPHA+1)*(NOPEN+1)
      LSCR = 0
      DO IOPEN = MINOP, MAXOP
        IF(MOD(IOPEN-MS2,2).EQ.0) THEN
          IALPHA = (IOPEN+MS2)/2
          L = 2*IOPEN + (IALPHA+1)*(IOPEN+1)
          LSCR = MAX(L,LSCR)
        END IF
      END DO
      Call mma_allocate(iSCR2,LSCR,Label='iSCR2')
*
      IDTBS = 1
C-jwk      DO IOPEN = 0, MAXOP
      DO IOPEN = MINOP, MAXOP
        ITP = IOPEN + 1
        IF( ITP .EQ. 1 ) THEN
          IDTBS = 1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NPCMCNF(ITP-1)
        END IF

        IALPHA = (IOPEN+MS2)/2
        IDET = NPCMCNF(IOPEN+1)
C?      WRITE(6,*) ' IOPEN, IDET = ', IOPEN, IDET
        CALL REO_PTDET(   IOPEN,
     &                   IALPHA,
     &                 Z_PTDT(ITP)%I,
     &                 REO_PTDT(ITP)%I,
     &                 DFTP(IDTBS),
*
     &                     IDET,iSCR2)
      END DO

*
*. matrix expressing csf's in terms of determinants
*
*
*. Tired of compiler warnings
      IDTBS = 0
      ICSBS = 0
      ICDCBS = 0
      DO IOPEN = 0, MAXOP
        ITP = IOPEN + 1
        IF( ITP .EQ. 1 ) THEN
          IDTBS = 1
          ICSBS = 1
          ICDCBS =1
        ELSE
          IDTBS = IDTBS + (IOPEN-1)*NPCMCNF(ITP-1)
          ICSBS = ICSBS + (IOPEN-1)*NPCSCNF(ITP-1)
          ICDCBS = ICDCBS + NPCMCNF(ITP-1)*NPCSCNF(ITP-1)
        END IF
C       IF(NPCMCNF(ITP)*NPCSCNF(ITP).EQ.0) GOTO 30
        IF( NTEST .GE. 5 ) THEN
          WRITE(6,*)
          WRITE(6,*) ' ************************************'
          WRITE(6,*) ' CSF - SD/COMB transformation matrix '
          WRITE(6,*) ' ************************************'
          WRITE(6,'(A)')
          WRITE(6,'(A,I3,A)')
     &    '  Type with ',IOPEN,' open orbitals '
          WRITE(6,'(A)')
     &    '  ************************************'
          WRITE(6,*)
        END IF
        IF(IOPEN .EQ. 0 ) THEN
           DTOC(ICDCBS) = 1.0D0
        ELSE
          CALL CSFDET_LUCIA(  IOPEN,
     &                      DFTP(IDTBS),
     &                      NPCMCNF(ITP),
     &                      CFTP(ICSBS),
     &                      NPCSCNF(ITP),
*
     &                      DTOC(ICDCBS),SCR1,SIZE(SCR1),PSSIGN,IPRCSF)
C              CSFDET(NOPEN,IDET,NDET,ICSF,NCSF,CDC,WORK,PSSIGN,
C    &                IPRCSF)
        END IF
      END DO
*     ^ End of loop over number of open shells
*
      Call mma_deallocate(SCR1)
      Call mma_deallocate(iSCR2)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)  ' List of CSF-SD transformation matrices '
        WRITE(6,*)  ' ======================================='
        WRITE(6,*)
        IDTBS = 1
        ICSBS = 1
        ICDCBS = 1
        DO IOPEN = 0, MAXOP
          ITP = IOPEN + 1
          IF( ITP .EQ. 1 ) THEN
            IDTBS = 1
            ICSBS = 1
            ICDCBS =1
          ELSE
            IDTBS = IDTBS + (IOPEN-1)*NPCMCNF(ITP-1)
            ICSBS = ICSBS + (IOPEN-1)*NPCSCNF(ITP-1)
            ICDCBS = ICDCBS + NPCMCNF(ITP-1)*NPCSCNF(ITP-1)
          END IF
          NNCS = NPCSCNF(ITP)
          NNCM = NPCMCNF(ITP)
          IF(NNCS.GT.0.AND.NNCM.GT.0) THEN
            WRITE(6,*) ' Number of open shells : ', IOPEN
            WRITE(6,*) ' Number of combinations per conf ', NNCM
            WRITE(6,*) ' Number of CSFs per conf         ', NNCS
            CALL WRTMAT(DTOC(ICDCBS),NNCM,NNCS,NNCM,NNCS)
          END IF
        END DO
      END IF
*
      END SUBROUTINE CSDTMT_GAS
