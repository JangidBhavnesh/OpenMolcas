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
* Copyright (C) 2000, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2000  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE GETSGM2(ILEV,JLEV,ISYCI,CI,SGM)
      use Symmetry_Info, only: Mul
      use gugx, only:  SGS, CIS, EXS
      IMPLICIT None

#include "pt2_guga.fh"

      Integer :: ILEV, JLEV, ISYCI
      Real*8, Intent(In) ::  CI(MXCI)
      Real*8, Intent(Out)::  SGM(MXCI)
      Integer IS, JS, IJS, ISSG, NSGM, nMidV, nICoup, MxEO, nVTab
      nMidV = CIS%nMidV
      MxEO  = EXS%MxEO
      nICoup=Size(EXS%ICoup)/3
      nVTab =Size(EXS%VTab)

C GIVEN CI COUPLING LEVELS ILEV, JLEV, COMPUTE SGM=E(ILEV,JLEV)*CI
C ILEV,JLEV ARE IN PRINCIPLE ACTIVE ORBITAL NUMBERS, BUT POSSIBLY
C IN ANOTHER ORDER THAN THE USUAL ONE -- HERE WE USE THE ORDER
C FOLLOWED BY THE GUGA COUPLING SCHEME.
C
C THIS ROUTINE REPLACES EARLIER GETSGM, TO GET RID OF THE PACKING AND
C STORING USED EARLIER.

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C NOTE!! THE EARLIER CALL GETSGM(ILEV,JLEV,IDARR,SGM) IS REPLACED BY
C GETSGM2(ILEV,JLEV,CI,SGM)!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IS=SGS%ISM(ILEV)
      JS=SGS%ISM(JLEV)
      IJS=MUL(IS,JS)
      ISSG=MUL(IJS,ISYCI)
      NSGM=CIS%NCSF(ISSG)
      IF(NSGM.EQ.0) RETURN

      SGM(1:NSGM)=0.0D0
      CALL SIGMA1_CP2(SGS,CIS,EXS,
     &     ILEV,JLEV,1.0D00,ISYCI,CI,SGM,
     &      CIS%NOCSF,CIS%IOCSF,CIS%NOW,CIS%IOW,
     &      EXS%NOCP,EXS%IOCP,EXS%ICOUP,
     &      EXS%VTAB,EXS%MVL,EXS%MVR,nMidV,nICoup,MxEO,nVTab)

      END SUBROUTINE GETSGM2
