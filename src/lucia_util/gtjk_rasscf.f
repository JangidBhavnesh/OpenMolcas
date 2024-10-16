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
C
C
      SUBROUTINE GTJK_RASSCF(RJ,RK,NAC,IREOST)
C
C     PURPOSE: GET ALL INTEGRALS COULOMB AND EXCHANGE INTEGRALS
C              WITH THE CHARGE DISTRIBUTION JK
C
*. Modified by addition of IREOST, August 2003.
      use wadr, only: TUVX
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input : Reorder array, symmetry => type (sic!)
      INTEGER IREOST(*)
      Real*8 RJ(*),RK(*)
C
C     FORM THE COULOMB (RJ) AND EXCHANGE (RK) INTEGRAL MATRICES FROM
C     THE TWO-ELECTRON INTEGRAL LIST
C
      NTUT=0
      DO 100 NT=1,NAC
      DO 101 NU=1,NT
       NT_REO = IREOST(NT)
       NU_REO = IREOST(NU)
*
       NTU_REO = NAC*(NT_REO-1) + NU_REO
       NUT_REO = NAC*(NU_REO-1) + NT_REO
       NTUT=NTUT+1
       NTUK=(NTUT**2+NTUT)/2
       RK(NTU_REO)=TUVX(NTUK)
       RK(NUT_REO)=TUVX(NTUK)
C
       NTT=(NT**2+NT)/2
       NTUJ=(NTT**2-NTT)/2+(NU**2+NU)/2
       RJ(NTU_REO)=TUVX(NTUJ)
       RJ(NUT_REO)=TUVX(NTUJ)
101   CONTINUE
100   CONTINUE
C
C     EXIT
C
      END
