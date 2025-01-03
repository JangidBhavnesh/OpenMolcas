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
* Copyright (C) 1995, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ADADST_GAS(    IOB,  IOBSM,  IOBTP,   NIOB,    JOB,
     &                        JOBSM,  JOBTP,   NJOB,  ISPGP,    ISM,
     &                          ITP,   KMIN,   KMAX,     I1,   XI1S,
     &                          LI1,     NK,   IEND,  IFRST,  KFRST,
     &                          I12,    K12, SCLFAC)
      use HIDSCR, only: ZSCR, ZOCSTR => OCSTR, REO, Z
      use lucia_data, only: NGAS
      use lucia_data, only: IBSPGPFTP,NELFSPGP,NELFTP
      use lucia_data, only: NELIS,NSTRKS
      use lucia_data, only: IOBPTS,NOBPT,NOCOB
*
*
*
* Obtain mappings
* a+IORB a+ JORB !KSTR> = +/-!ISTR>
* In the form
* I1(KSTR) =  ISTR if a+IORB a+ JORB !KSTR> = +/-!ISTR> , ISTR is in
* ISPGP,ISM,IGRP.
* (numbering relative to TS start)
*. Only excitations IOB. GE. JOB are included
* The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
* and IOB (JOB) is the first orbital to be used, and the number of orbitals
* to be checked is NIOB ( NJOB).
*
* Only orbital pairs IOB .gt. JOB are included
*
* The output is given in I1(KSTR,I,J) = I1 ((KSTR,(J-1)*NIOB + I)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* Only Kstrings with relative numbers from KMIN to KMAX are included
*
* If IEND .ne. 0 last string has been checked
*
* Jeppe Olsen , August of 95
*
* ======
*. Input
* ======
*
      IMPLICIT NONE
      INTEGER IOB, IOBSM, IOBTP, NIOB, JOB, JOBSM, JOBTP, NJOB, ISPGP,
     &        ISM, ITP, KMIN, KMAX, LI1, NK, IEND, IFRST, KFRST, I12,
     &        K12
      REAL*8 SCLFAC
*
      INTEGER I1(*)
      REAL*8 XI1S(*)
*
      INTEGER IDUM_ARR(1)
      INTEGER NTEST,ISPGPABS,K1SM,K1SPGPABS,KSM,KSPGPABS,NTEST2,NELI,
     &        NELK,NSTRK,IIOB,JJOB,NSTRI
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' ADADST_GAS in service '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        WRITE(6,*) ' IOB,IOBSM,IOBTP ', IOB,IOBSM,IOBTP
        WRITE(6,*) ' JOB,JOBSM,JOBTP ', JOB,JOBSM,JOBTP
      END IF
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) 'Problemo, ADADST '
C?      WRITE(6,*) ' SCLFAC = ',SCLFAC
C?    END IF

*
*. Internal affairs
*
      IF(I12.GT.SIZE(Z,2).OR.K12.GT.SIZE(ZOCSTR,2)) THEN
        WRITE(6,*) ' ADST_GAS : Illegal value of K12 = ', K12
        WRITE(6,*) ' ADST_GAS : Illegal value of I12 = ', I12
*        STOP' ADST_GAS : Illegal value of I12  '
        CALL SYSABENDMSG('lucia_util/adst_gas','Internal error',' ')
        RETURN
      END IF

*
*. Supergroup and symmetry of K strings
*
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      CALL NEWTYP(ISPGPABS,1,IOBTP,K1SPGPABS)
      CALL NEWTYP(K1SPGPABS,1,JOBTP,KSPGPABS)
      CALL SYMCOM(2,0,IOBSM,K1SM,ISM)
      CALL SYMCOM(2,0,JOBSM,KSM,K1SM)
      IF(NTEST.GE.100) WRITE(6,*)
     & ' K1SM,K1SPGPABS,KSM,KSPGPABS : ',
     &   K1SM,K1SPGPABS,KSM,KSPGPABS
* In ADADS1_GAS we need : Occupation of KSTRINGS
*                         lexical => Actual order for I strings
* Generate if required
*
      IF(IFRST.NE.0) THEN
*.. Generate information about I strings
*. Arc weights for ISPGP
        NTEST2 = NTEST
        CALL WEIGHT_SPGP(Z(:,I12),NGAS,NELFSPGP(1,ISPGPABS),
     &                      NOBPT,ZSCR,NTEST2)
        NELI = NELFTP(ITP)
        NELIS(I12) = NELI
*. Reorder array for I strings
        CALL GETSTR_TOTSM_SPGP(    ITP,  ISPGP,    ISM,   NELI,  NSTRI,
     &                         ZOCSTR(:,K12),NOCOB,1,
     &                         Z(:,I12),REO(:,I12))
      END IF
      NELK = NELIS(I12) - 2
      IF(KFRST.NE.0) THEN
*. Generate occupation of K STRINGS
       CALL GETSTR_TOTSM_SPGP(      1,KSPGPABS,   KSM,  NELK, NSTRK,
     &                        ZOCSTR(:,K12),NOCOB,  0,IDUM_ARR,IDUM_ARR)
       NSTRKS(K12) = NSTRK
      END IF
*
      NSTRK = NSTRKS(K12)
*
      IIOB = IOBPTS(IOBTP,IOBSM) + IOB - 1
      JJOB = IOBPTS(JOBTP,JOBSM) + JOB - 1
      CALL ADADS1_GAS(       NK,       I1,     XI1S,      LI1,     IIOB,
     &                     NIOB,     JJOB,    NJOB,ZOCSTR(:,K12),  NELK,
     &                    NSTRK,REO(:,I12),Z(:,I12),NOCOB, KMAX,
     &                     KMIN,     IEND,   SCLFAC)
*
      END SUBROUTINE ADADST_GAS
