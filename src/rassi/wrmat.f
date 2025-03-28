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
* Copyright (C) 1989, Per Ake Malmqvist                                *
************************************************************************
*****************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST
*  SUBROUTINE WRMAT     IBM-3090 RELEASE 89 01 31
*  PRINT OUT A SYMMETRY-BLOCKED MATRIX. THE COMBINED SYMMETRY
*  OF ROWS AND COLUMNS IS ISY12. NDIM1 GIVES THE NUMBER OF
*  OF ROWS WITHIN EACH SYMMETRY TYPE, SIMILAR NDIM2, COLUMNS.
*****************************************************************
      SUBROUTINE WRMAT(TEXT,ISY12,NDIM1,NDIM2,NMAT,XMAT)
      use Symmetry_Info, only: MUL, nSym=>nIrrep
      IMPLICIT NONE
      CHARACTER(Len=*) TEXT
      Integer ISY12, NMAT
      Integer NDIM1(NSYM),NDIM2(NSYM)
      Real*8 XMAT(NMAT)

      Integer ISTA, ISY1, ISY2, NN
      ISTA=1
      WRITE(6,'(/,1X,A,/)') TEXT
      DO ISY1=1,NSYM
        ISY2=MUL(ISY1,ISY12)
        NN=NDIM1(ISY1)*NDIM2(ISY2)
        IF (NN /= 0) THEN
          WRITE(6,*)
          WRITE(6,'(A,2I2)')' SYMMETRY LABELS OF ROWS/COLS:',ISY1,ISY2
          CALL WRMAT1(NDIM1(ISY1),NDIM2(ISY2),XMAT(ISTA))
        END IF
        ISTA=ISTA+NN
      END DO
      WRITE(6,*)
      WRITE(6,*) repeat('*',80)

      END SUBROUTINE WRMAT
