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
      SUBROUTINE GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTROC,  NSTR,   NEL,
     &                              NNSTSGP,IISTSGP)
      use strbas, only: OCSTR
      use lucia_data, only: NGAS
      use lucia_data, only: NELFGP
      use lucia_data, only: MXPNSMST,MXPNGAS
*
* Obtain all superstrings containing  strings of given sym and type
*
* ( Superstring :contains electrons belonging to all gasspaces
*        string :contains electrons belonging to a given GAS space
* A super string is thus a product of NGAS strings )
*
* Jeppe Olsen, Summer of 95
*              Optimized version, october 1995
*
*. In this subroutine the ordering of strings belonging to a given type
*  is defined !!
* Currently we are using the order
* Loop over GAS 1 strings
*  Loop over GAS 2 strings
*   Loop over GAS 3 strings --
*
*     Loop over gas N strings
*
      IMPLICIT NONE
      INTEGER NSTR,NEL
*. Specific input
      INTEGER ITPFGS(*), ISMFGS(*)
      INTEGER NNSTSGP(MXPNSMST,*), IISTSGP(MXPNSMST,*)
*. Local scratch
      INTEGER NSTFGS(MXPNGAS), IBSTFGS(MXPNGAS)
*. Output
      INTEGER ISTROC(NEL,*)

      INTEGER IGAS,IGASL,NSTRTOT,NELB,NELI,NSTA,JGAS,NSTB,NSTI
*. Number of strings per GAS space
C?    write(6,*) ' entering problem child '
      DO IGAS = 1, NGAS
        NSTFGS(IGAS)  = NNSTSGP(ISMFGS(IGAS),IGAS)
        IBSTFGS(IGAS) = IISTSGP(ISMFGS(IGAS),IGAS)
      END DO
*
#ifdef _DEBUGPRINT_
        WRITE(6,*) '  GETSTR_GASSM_SPGP speaking '
        WRITE(6,*) '  =========================== '
        WRITE(6,*) ' ISMFGS,ITPFGS (input) '
        CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' NSTFGS, IBSTFGS ( intermediate results ) '
        CALL IWRTMA(NSTFGS,1,NGAS,1,NGAS)
        CALL IWRTMA(IBSTFGS,1,NGAS,1,NGAS)
#endif
*. Last gasspace with a nonvanishing number of electrons
      IGASL = 0
      DO IGAS = 1, NGAS
        IF( NELFGP(ITPFGS(IGAS)) .NE. 0 ) IGASL = IGAS
      END DO
*
      NSTRTOT = 1
      DO IGAS = 1, NGAS
        NSTRTOT = NSTRTOT*NSTFGS(IGAS)
      END DO
      IF(IGASL.EQ.0) GOTO 2810
*
      IF(NSTRTOT.EQ.0) GOTO 1001
*. Loop over GAS spaces
      NELB = 0
      DO IGAS = 1, IGASL
*. Number of electrons in GAS = 1, IGAS - 1
        IF(IGAS.GT.1) THEN
          NELB = NELB +  NELFGP(ITPFGS(IGAS-1))
        END IF
*. Number of electron in IGAS
        NELI = NELFGP(ITPFGS(IGAS))
        IF(NELI.GT.0) THEN

*. The order of strings corresponds to a matrix A(I(after),Igas,I(before))
*. where I(after) loops over strings in IGAS+1 - IGASL and
*  I(before) loop over strings in 1 - IGAS -1
          NSTA = 1
          DO JGAS = IGAS+1, IGASL
            NSTA = NSTA * NSTFGS(JGAS)
          END DO
*
          NSTB =  1
          DO JGAS = 1, IGAS-1
            NSTB = NSTB * NSTFGS(JGAS)
          END DO
*
          NSTI = NSTFGS(IGAS)

#ifdef _DEBUGPRINT_
            WRITE(6,*) ' NSTI,NSTB,NSTA,NELB,NELI,NEL ',
     &                   NSTI,NSTB,NSTA,NELB,NELI,NEL
            WRITE(6,*) ' IBSTFGS(IGAS)',IBSTFGS(IGAS)
#endif
*
          CALL ADD_STR_GROUP(NSTI,IBSTFGS(IGAS),
     &                       OCSTR(ITPFGS(IGAS))%I,
     &                       NSTB,   NSTA, ISTROC, NELB+1,   NELI,
     &                           NEL)
*
*. Loop over strings in IGAS
        END IF
      END DO
 1001 CONTINUE
 2810 CONTINUE
      NSTR = NSTRTOT
*
#ifdef _DEBUGPRINT_
        WRITE(6,*) ' Info from  GETSTR_GASSM_SPGP '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry and type strings : '
        WRITE(6,*)
        WRITE(6,*) '   AS    Sym  Type '
        WRITE(6,*) ' =================='
        DO IGAS = 1, NGAS
          WRITE(6,'(3I6)') IGAS,ISMFGS(IGAS),ITPFGS(IGAS)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Number of strings generated : ', NSTR
        WRITE(6,*) ' Strings generated '
        CALL PRTSTR(ISTROC,NEL,NSTR)
#endif
*
      END SUBROUTINE GETSTRN_GASSM_SPGP
