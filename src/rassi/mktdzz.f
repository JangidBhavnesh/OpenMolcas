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
      SUBROUTINE MKTDZZ(CMOA,CMOB,TDMAB,TDMZZ,iRC)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: LSYM1, LSYM2
      use Symmetry_Info, only: nSym=>nIrrep, MUL
      use rassi_data, only: NCMO,NTDMAB,NTDMZZ,NBASF,NOSH
      IMPLICIT None
      Real*8 CMOA(NCMO),CMOB(NCMO)
      Real*8 TDMAB(NTDMAB),TDMZZ(NTDMZZ)
      Integer iRC

      Integer ISTCMO(8)
      Real*8, Allocatable:: SCR(:)
      Integer ISY12, NSCR, IST, ISY1, NO1, ISY2, ISTTA, ISTCA, ISTTZ,
     &        ISTCB, NO2, NB1, NB2

      If (iRC.eq.0) Then
         TDMZZ(:)=0.0D0
         Return
      End If
C ISTCMO()=START INDEX FOR CMO ARRAY SYMMETRY BLOCKS.
      ISY12=MUL(LSYM1,LSYM2)
C NSCR=SIZE NEEDED FOR TEMPORARY MATRIX PRODUCT.
      NSCR=0
      IST=1
      DO ISY1=1,NSYM
        ISTCMO(ISY1)=IST
        NO1=NOSH(ISY1)
        IST=IST+NO1*NBASF(ISY1)
        ISY2=MUL(ISY1,ISY12)
        NSCR=MAX(NSCR,NO1*NBASF(ISY2))
      end do
      CALL mma_allocate(SCR,NSCR,Label='SCR')
      ISTTA=1
      ISTCA=1
      ISTTZ=1
      DO ISY1=1,NSYM
        ISY2=MUL(ISY1,ISY12)
        ISTCB=ISTCMO(ISY2)
        NO1=NOSH(ISY1)
        NO2=NOSH(ISY2)
        NB1=NBASF(ISY1)
        NB2=NBASF(ISY2)
        if (NB1*NB2 /= 0) then
        IF(NO1*NO2.EQ.0) THEN
          CALL FZERO(TDMZZ(ISTTZ),NB1*NB2)
        ELSE
          CALL DGEMM_('N','T',NO1,NB2,NO2,1.0D0,
     &                TDMAB(ISTTA),NO1,CMOB(ISTCB),NB2,
     &         0.0D0,  SCR,NO1)
          CALL DGEMM_('N','N',NB1,NB2,NO1,1.0D0,
     &                 CMOA(ISTCA),NB1,SCR,NO1,
     &         0.0D0,  TDMZZ(ISTTZ),NB1)
          ISTTA=ISTTA+NO1*NO2
        END IF
        end if
        ISTCA=ISTCA+NB1*NO1
        ISTTZ=ISTTZ+NB1*NB2
      end do
      CALL mma_deallocate(SCR)

      END SUBROUTINE MKTDZZ
