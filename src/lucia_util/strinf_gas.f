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
* Copyright (C) 1994, Jeppe Olsen                                      *
*               2024, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE STRINF_GAS(IPRNT)
      use stdalloc, only: mma_allocate, mma_deallocate
      use strbas, only: ZMAT,NSTSGP,ISTSGP,STREO,OCSTR,STSTM,NSTSO,
     &                  ISTSO,IOCLS,SPGPAN,SPGPCR
      use lucia_data, only: NGAS,IGSOCC,IPHGAS,NMXOCCLS
* modification Jeppe + Giovanni + Dongxia.
* G. Li Manni, June 2024: Scale-up capability for single SD ROHF type calculations
      use distsym, only: INGRP_VAL,ISMDFGP,ISMSCR,NACTSYM
      use lucia_data, only: MS2
      use lucia_data, only: NGRP,NTSPGP,MXNSTR,MXSMCLS,MXSMCLSE,
     &                      MXSMCLSE1,MAX_STR_OC_BLK,MAX_STR_SPGP,
     &                      MINMAX_SM_GP,IBSPGPFTP,IGSFGP,ISPGPFTP,
     &                      ISTFSMGP,NELFGP,NELFSPGP,NELFTP,NHLFSPGP,
     &                      NSPGPFTP,NSTFGP,NSTFSMGP,NSTFSMSPGP
      use lucia_data, only: NACTEL
      use lucia_data, only: NACOB,NORB1,NORB2,NORB3,NOBPT
      use lucia_data, only: ISTAC
      use lucia_data, only: NSTTYP
      use lucia_data, only: MXPNSMST,MXPNGAS
      use csm_data, only: NSMST
*
* Obtain string information for GAS expansion
*
* =====
*.Input
* =====
*
* /LUCINP/,/ORBINP/,/CSM/, /CGAS/, /GASSTR/
*
* =====
*.Output
* =====
*
* /STRINP/,/STINF/,STRBAS and string information in STIN
*
      IMPLICIT NONE
      INTEGER IPRNT

      INTEGER ZERO_ARR(1), IDUM(1)
      Integer, Allocatable:: FREEL(:)
*. A bit of scratch
      INTEGER, External:: IELSUM
      INTEGER LAC,NTEST,IEL,IGRP,NACOB_EFFECTIVE,MAXSCR,IGAS,MNRS1X,
     &        MXRS1X,MNRS3X,MXRS3X,IOCTYPX,IGP,MX,ISM,MN,NGSOBP,IGSOB,
     &        NSTINI,IEC,LROW,IZERO,JGRP,ITP,IGRPABS,NSMCLS,NSMCLSE,
     &        NSMCLSE1,IISPGP,NHOLE,ISPGP,NSTR,NEL,ISTSM,ISTTYP,IIEL,
     &        ISTTYPC,JSTTYP,ISTTYPA,MXNSTRFSG
*
*
* Some dummy initializtions
      LAC = 0 ! jwk-cleanup
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
*
**.2 : Number of classes per string type and mappings between
**.    string types (/STINF/)
*
      If(NActEl.ne.MS2.or.NActEl.ne.NACOB) CALL ZSTINF_GAS(IPRNT)
*
**.3 : Static memory for string information
*
      CALL MEMSTR_GAS
*
** 4 : Info about group of strings
*
*.First free address
*
*     Find maximum needed length of scratch.
*       MAXSCR = 2*NACOB+(IEL+1)(NACOB+1)
*       with IEL = MAX(NELFGP(IGRP=1,NGRP)
*
      IEL = 0
      DO IGRP = 1, NGRP
         IEL = MAX(IEL, NELFGP(IGRP))
      ENDDO
      NACOB_EFFECTIVE = NACOB
      IF (NACOB .EQ. 0) NACOB_EFFECTIVE = 1
      MAXSCR = 2*NACOB_EFFECTIVE +(IEL+1)*(NACOB_EFFECTIVE+1) +NSMST
      Call mma_allocate(FREEL,MAXSCR,Label='FREEL')
      DO IGRP = 1, NGRP
*. A gas group can be considered as a RAS group with 0 electrons in
*  RAS1, RAS3 !
        IGAS = IGSFGP(IGRP)
        IF(IGAS.EQ.1) THEN
          NORB1 = 0
        ELSE
          NORB1 = IELSUM(NOBPT,IGAS-1)
        END IF
        NORB2 = NOBPT(IGAS)
        NORB3 = NACOB-NORB1-NORB2
        MNRS1X = 0
        MXRS1X = 0
        MNRS3X = 0
        MXRS3X = 0
        IEL = NELFGP(IGRP)
        IOCTYPX = 1
*. Reverse lexical adresing schemes for each group of string
        CALL WEIGHT_LUCIA(Zmat(IGRP)%I,IEL,  NORB1,  NORB2,  NORB3,
     &                      MNRS1X,MXRS1X, MNRS3X,
     &                      MXRS3X,FREEL,
     &                      IPRNT )
*. Number of strings per symmetry in a given group
        CALL NSTRSO_GAS(     IEL,   NORB1,   NORB2,   NORB3,  MNRS1X,
     &                    MXRS1X,  MNRS3X,  MXRS3X,FREEL,NACOB,
     &                  NSTSGP(1)%I,
     &                  ISTSGP(1)%I,IOCTYPX,NSMST,IGRP,IPRNT)
*. Construct the strings ordered by symmetry
        CALL GENSTR_GAS(     IEL,  MNRS1X,  MXRS1X,  MNRS3X,  MXRS3X,
     &                  ISTSGP(1)%I,
     &                  IGRP,IOCTYPX,NSMST,Zmat(IGRP)%I,
     &                  FREEL,STREO(IGRP)%I,OCSTR(IGRP)%I,
*
     &                  FREEL(1+IOCTYPX*NSMST),IGRP,IPRNT)
*
       CALL ICOPVE2(NSTSGP(1)%I,1+(IGRP-1)*NSMST,NSMST,
     &              NSTFSMGP(1,IGRP))
       CALL ICOPVE2(ISTSGP(1)%I,1+(IGRP-1)*NSMST,NSMST,
     &              ISTFSMGP(1,IGRP))
      END DO
      Call mma_deallocate(FREEL)
*
      INGRP_VAL = NGRP
      CALL mma_allocate(ISMDFGP, NSMST*NGRP,Label='ISMDFGP')
      CALL mma_allocate(NACTSYM, NGRP,Label='NACTSYM')
      CALL mma_allocate(ISMSCR, NGRP,Label='ISMSCR')
      call SMDFGP_GEN(NGRP,NSMST,MXPNSMST,NSTFSMGP,
     &                NACTSYM,ISMDFGP)
*
      IF(NTEST.GE.10) THEN
        write(6,*) 'NGRP', NGRP
        write(6,*) 'NSMST*NGRP', NSMST*NGRP
        WRITE(6,*) ' Number of strings per group and symmetry '
        CALL IWRTMA10(NSTSGP(1)%I,NSMST,NGRP,NSMST,NGRP)
        WRITE(6,*) ' Number of strings per group and symmetry(2) '
        CALL IWRTMA10(NSTFSMGP,NSMST,NGRP,MXPNSMST,NGRP)
      END IF
*
*. Min and max of sym for each group
*
      DO IGP = 1, NGRP
       MX = 1
       DO ISM = 1, NSMST
         IF(NSTFSMGP(ISM,IGP).GT.0) MX = ISM
       END DO
*
       MN = NSMST
       DO ISM = NSMST,1,-1
         IF(NSTFSMGP(ISM,IGP).GT.0) MN = ISM
       END DO
*
       MINMAX_SM_GP(1,IGP) = MN
       MINMAX_SM_GP(2,IGP) = MX
*
      END DO
      IF(NTEST.GT.5) THEN
        WRITE(6,*) ' MINMAX array for sym of groups '
        WRITE(6,*) ' =============================='
        CALL IWRTMA(MINMAX_SM_GP,2,NGRP,2,NGRP)
      END IF
*
*
* 4.5 : Creation/Annihilation mappings between different
*       types of strings
*
      DO IGRP = 1, NGRP
*
        IGAS = IGSFGP(IGRP)
        NGSOBP = NOBPT(IGAS)
*. First orbital in GAS spacce
        IGSOB = IELSUM(NOBPT,IGAS-1)+1
        IEL = NELFGP(IGRP)
        NSTINI = NSTFGP(IGRP)
*
*. Type of mapping : Only creation                  (LAC = 1)
*                    Only annihilation              (LAC = 2)
*                    Both annihilation and creation (LAC = 3)
* If only annihilation is present the string mapping arrays
* will only be over electronns
        IF(     ISTAC(IGRP,1).NE.0.AND.ISTAC(IGRP,2).NE.0) THEN
          LAC = 3
          IEC = 1
          LROW = NGSOBP
        ELSE IF(ISTAC(IGRP,1).NE.0.AND.ISTAC(IGRP,2).EQ.0) THEN
          LAC = 1
          IEC = 2
          LROW = IEL
        ELSE IF(ISTAC(IGRP,1).EQ.0.AND.ISTAC(IGRP,2).NE.0) THEN
          LAC = 2
          IEC = 0
          LROW = NGSOBP
        ELSE IF(ISTAC(IGRP,1).EQ.0.AND.ISTAC(IGRP,2).EQ.0) THEN
          LAC = 0
          IEC = 0
          LROW = 0
        END IF
*. Zero
        IF(LAC.NE.0) THEN
          IZERO = 0
          CALL ISETVC(STSTM(IGRP,1)%I,IZERO,LROW*NSTINI)
          CALL ISETVC(STSTM(IGRP,2)%I,IZERO,LROW*NSTINI)
        END IF
*
        IF(ISTAC(IGRP,2).NE.0) THEN
          JGRP = ISTAC(IGRP,2)
          CALL CRESTR_GAS(OCSTR(IGRP)%I,
     &                    NSTFGP(IGRP),NSTFGP(JGRP),IEL,NGSOBP,  IGSOB,
     &                    Zmat(JGRP)%I,STREO(JGRP)%I,
     &                    0,IDUM,IDUM,
     &                    STSTM(IGRP,1)%I,
     &                    STSTM(IGRP,2)%I,NACOB,IPRNT)
*
        END IF
        IF(ISTAC(IGRP,1).NE.0) THEN
          JGRP = ISTAC(IGRP,1)
          CALL ANNSTR_GAS(OCSTR(IGRP)%I,
     &                    NSTFGP(IGRP),NSTFGP(JGRP),IEL,NGSOBP,  IGSOB,
     &                    Zmat(JGRP)%I,STREO(JGRP)%I,
     &                    0,IDUM,IDUM,
     &                    STSTM(IGRP,1)%I,
     &                    STSTM(IGRP,2)%I,NACOB,IEC,LROW,IPRNT)
*
        END IF
      END DO
*
*
*. Now to supergroups , i.e. strings of with given number of elecs in
*  each GAspace
*
      CALL ISETVC(NSTFSMSPGP,0,MXPNSMST*NTSPGP)
      MXNSTR = -1
      DO ITP = 1, NSTTYP
*. Loop over supergroups of given type . i.e. strings
*  with given occupation in each GAS space
        DO IGRP = 1, NSPGPFTP(ITP)
          IGRPABS = IGRP-1 + IBSPGPFTP(ITP)
          CALL NSTPTP_GAS(    NGAS,
     &                    ISPGPFTP(1,IGRPABS),
     &                    NSTSGP(1)%I,
     &                   NSMST,NSTSO(ITP)%I,IGRP,MXNSTRFSG,NSMCLS,
     &                     NSMCLSE,NSMCLSE1)
*
*
          MXSMCLS   = MAX(MXSMCLS,NSMCLS)
          MXSMCLSE  = MAX(MXSMCLSE,NSMCLSE)
          MXSMCLSE1 = MAX(MXSMCLSE1,NSMCLSE1)
*
          MXNSTR = MAX(MXNSTR,MXNSTRFSG)
        END DO
*
        CALL ICOPMT(NSTSO(ITP)%I,
     &                  NSMST,
     &              NSPGPFTP(ITP),
     &              NSTFSMSPGP(1,IBSPGPFTP(ITP)),
     &               MXPNSMST,NSPGPFTP(ITP))
*. Corresponding offset array : Each supergroup is generated individually
*. so each supergroup starts with offset 1 !
        CALL ZSPGPIB(NSTSO(ITP)%I,ISTSO(ITP)%I,NSPGPFTP(ITP),NSMST)
*
        IF(NTEST.GE.5) THEN
          WRITE(6,*)
     &    ' Number of strings per sym (row) and supergroup(column)',
     &    ' for type = ', ITP
          CALL IWRTMA(NSTSO(ITP)%I,NSMST,NSPGPFTP(ITP),
     &                NSMST,NSPGPFTP(ITP))
          WRITE(6,'(A,3I6)') ' NSMCLS,NSMCLSE,NSMCLSE1=',
     &                         NSMCLS,NSMCLSE,NSMCLSE1
          WRITE(6,*)
        END IF
*
      END DO
*. Number of electron in each AS for each supergroup
      CALL ZNELFSPGP(IPRNT)
*
* Number of holes per supergroup
      DO IISPGP = 1, NTSPGP
        NHOLE = 0
        DO IGAS = 1, NGAS
          IF(IPHGAS(IGAS).EQ.2) NHOLE = NHOLE + NELFSPGP(IGAS,IISPGP)
        END DO
        NHLFSPGP(IISPGP) = NHOLE
      END DO
      IF(NTEST.GE.10) THEN
      WRITE(6,*) ' Number of electrons in hole spaces per supergroup '
      CALL IWRTMA(NHLFSPGP,1,NTSPGP,1,NTSPGP)
      END IF
*. Largest number of strings belonging to given supergroup
*. Largest Occupation block for given supergroup and sym
      MAX_STR_OC_BLK = -1
      MAX_STR_SPGP = 0
      DO ISPGP = 1, NTSPGP
        NSTR = IELSUM(NSTFSMSPGP(1,ISPGP),NSMST)
        MAX_STR_SPGP = MAX(MAX_STR_SPGP,NSTR)
        NEL = IELSUM(NELFSPGP(1,ISPGP),NGAS)
        DO ISTSM = 1, NSMST
          MAX_STR_OC_BLK
     &  = MAX(MAX_STR_OC_BLK,(NEL+4)*NSTFSMSPGP(ISTSM,ISPGP))
CMOD &  = MAX(MAX_STR_OC_BLK,NEL*NSTFSMSPGP(ISTSM,ISPGP))
        END DO
      END DO
*
      IF(NTEST.GE.2) THEN
      WRITE(6,*)
     & ' Largest number of strings of given supergroup        ',
     & MAX_STR_SPGP
      WRITE(6,*) ' Largest block of string occupations ',
     &              MAX_STR_OC_BLK
*
      WRITE(6,*)
     & ' Largest number of strings of given supergroup and sym', MXNSTR
      END IF
C?    WRITE(6,'(A,3I6)') ' MXSMCLS,MXSMCLSE,MXSMCLSE1 = ',
C?   &                     MXSMCLS,MXSMCLSE,MXSMCLSE1
*
*
* Possible occupation classes
*
      ZERO_ARR(1)=0
      CALL OCCLS(         2,  NMXOCCLS,IOCLS,  NACTEL,    NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),       0,ZERO_ARR,   NOBPT)
*
* Maps creation/annihilation of given gas orb from given supergroup
* gives new supergroup.
*
      IZERO = 0
      CALL ISETVC(SPGPCR,IZERO,NGAS*NTSPGP)
      CALL ISETVC(SPGPAN,IZERO,NGAS*NTSPGP)
*
      DO ISTTYP = 1,NSTTYP
*. Creation map from this type
        IIEL = NELFTP(ISTTYP)
*. Type of string with one elec more
        ISTTYPC = 0
        DO JSTTYP = 1, NSTTYP
          IF(MOD(ISTTYP,2).EQ.MOD(JSTTYP,2).AND.
     &       NELFTP(JSTTYP).EQ.IIEL+1           ) ISTTYPC = JSTTYP
        END DO
C?      WRITE(6,*) ' ISTTYP and ISTTYPC ',ISTTYP,ISTTYPC
        IF(NSPGPFTP(ISTTYP).GT.0) THEN
*
        IF(ISTTYPC.GE.1.AND.NSPGPFTP(ISTTYPC).GT.0) THEN
           CALL SPGP_AC(NELFSPGP(1,1),
     &                  NSPGPFTP(ISTTYP),
     &                  NELFSPGP(1,1),NSPGPFTP(ISTTYPC),NGAS,MXPNGAS,2,
     &                  SPGPCR,
     &                  IBSPGPFTP(ISTTYP),
     &                  IBSPGPFTP(ISTTYPC))
        ELSE
        END IF
*. Annihilation maps
        ISTTYPA = 0
        DO JSTTYP = 1, NSTTYP
          IF(MOD(ISTTYP,2).EQ.MOD(JSTTYP,2).AND.
     &       NELFTP(JSTTYP).EQ.IIEL-1           ) ISTTYPA = JSTTYP
        END DO
C?      WRITE(6,*) 'ISTTYP, ISTTYPA', ISTTYP,ISTTYPA
        IF(ISTTYPA.GE.1 .AND.NSPGPFTP(ISTTYPA).GT.0) THEN
           CALL SPGP_AC(NELFSPGP(1,1),
     &                  NSPGPFTP(ISTTYP),
     &                  NELFSPGP(1,1),NSPGPFTP(ISTTYPA),NGAS,MXPNGAS,1,
     &                  SPGPAN,
     &                  IBSPGPFTP(ISTTYP),
     &                  IBSPGPFTP(ISTTYPA))
        END IF
        END IF
      END DO
*
      END SUBROUTINE STRINF_GAS
