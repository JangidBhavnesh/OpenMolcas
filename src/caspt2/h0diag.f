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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE H0DIAG_CASPT2(ISYCI,DIAG,NOW,IOW)
      use pt2_guga_data
      IMPLICIT REAL*8 (A-H,O-Z)
C INPUT ARRAYS:

#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
      DIMENSION DIAG(MXCI),NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)

C PURPOSE: FORM AN ARRAY OF DIAGONAL HAMILTONIAN MATRIX ELEMENTS
C FOR THE SPECIFIED TOTAL SYMMETRY ISYCI

      CALL DCOPY_(MXCI,[0.0D0],0,DIAG,1)
      IEMU=1
      DO MV=1,NMIDV
        DO ISYUP=1,NSYM
          NUP=NOW(1,ISYUP,MV)
          IF(NUP.EQ.0) GOTO 30
          ISYDWN=MUL(ISYUP,ISYCI)
          NDWN=NOW(2,ISYDWN,MV)
          IF(NDWN.EQ.0) GOTO 30
          ICS=1+IOW(1,ISYUP,MV)
          JCS=1+IOW(2,ISYDWN,MV)
          NC=NUP*NDWN
          CALL DIELMV(ICASE(ICS),ICASE(JCS),NUP,NDWN,DIAG(IEMU))
          IEMU=IEMU+NC
  30      CONTINUE
        END DO
      END DO
      END SUBROUTINE H0DIAG_CASPT2
