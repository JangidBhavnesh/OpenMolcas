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
* Copyright (C) 1992,1994, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE POLY2(CI)
      use gugx, only: SGS
#ifdef _ENABLE_CHEMPS2_DMRG_
      use caspt2_output, only:iPrGlb
      use PrintLevel, only: debug
#endif
      IMPLICIT NONE
* PER-AAKE MALMQUIST, 92-12-07
* THIS PROGRAM CALCULATES 1-EL AND 2-EL
* DENSITY MATRICES FOR A CASSCF WAVE FUNCTION.
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"

      REAL*8, INTENT(IN) :: CI(NCONF)

      INTEGER LSGM1,LSGM2,LG1TMP,LG2TMP

      INTEGER I

#ifdef _ENABLE_CHEMPS2_DMRG_
      REAL*8, EXTERNAL :: DNRM2_
      INTEGER NAC4
#endif
      Integer :: nLev
      nLev = SGS%nLev

      IF(NLEV.GT.0) THEN
* NN.15 in case of DMRG-CASPT2, CI=1 and MXCI=1
        CALL GETMEM('LSGM1','ALLO','REAL',LSGM1 ,MXCI)
        CALL GETMEM('LSGM2','ALLO','REAL',LSGM2 ,MXCI)
        CALL GETMEM('LG1TMP','ALLO','REAL',LG1TMP,NG1)
        CALL GETMEM('LG2TMP','ALLO','REAL',LG2TMP,NG2)
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
        IF(.Not.DoCumulant) THEN
#endif
          CALL DENS2_RPT2(CI,WORK(LSGM1),WORK(LSGM2),
     &                    WORK(LG1TMP),WORK(LG2TMP),NLEV)
#ifdef _ENABLE_BLOCK_DMRG_
        ELSE
* this is provided by BLOCK code
* load 2PDM from scratch file generated by BLOCK
          CALL block_load2pdm(NASHT,WORK(LG2TMP),JSTATE,JSTATE)
* this is located under block_dmrg_util/
* compute 1PDM from 2PDM
          CALL TWO2ONERDM(NASHT,NACTEL,WORK(LG2TMP),WORK(LG1TMP))
        END IF
#elif _ENABLE_CHEMPS2_DMRG_
        ELSE
          NAC4 = NLEV * NLEV * NLEV * NLEV
          CALL chemps2_load2pdm( NASHT, WORK( LG2TMP ), MSTATE(JSTATE) )
          CALL TWO2ONERDM(NASHT,NACTEL,WORK(LG2TMP),WORK(LG1TMP))
          IF(iPrGlb.GE.DEBUG) THEN
            WRITE(6,'("DEBUG> ",A)')
     &        "CHEMPS2: norms of the density matrices:"
            WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G1:"
     &          , DNRM2_(NG1,WORK(LG1TMP),1)
            WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G2:"
     &          , DNRM2_(NG2,WORK(LG2TMP),1)
          ENDIF
        END IF
#endif
      END IF

* REINITIALIZE USE OF DMAT.
* The fields IADR10 and CLAB10 are kept in common included from pt2_guga.fh
* CLAB10 replaces older field called LABEL.
      DO I=1,64
        IADR10(I,1)=-1
        IADR10(I,2)=0
        CLAB10(I)='   EMPTY'
      END DO
      IADR10(1,1)=0
* HENCEFORTH, THE CALL PUT(NSIZE,LABEL,ARRAY) WILL ENTER AN
* ARRAY ON LUDMAT AND UPDATE THE TOC.
      IF(NLEV.GT.0) THEN
        CALL PT2_PUT(NG1,' GAMMA1',WORK(LG1TMP))
        CALL PT2_PUT(NG2,' GAMMA2',WORK(LG2TMP))

        CALL GETMEM('LSGM1','FREE','REAL',LSGM1 ,MXCI)
        CALL GETMEM('LSGM2','FREE','REAL',LSGM2 ,MXCI)
        CALL GETMEM('LG1TMP','FREE','REAL',LG1TMP,NG1)
        CALL GETMEM('LG2TMP','FREE','REAL',LG2TMP,NG2)
      END IF

      END SUBROUTINE POLY2
