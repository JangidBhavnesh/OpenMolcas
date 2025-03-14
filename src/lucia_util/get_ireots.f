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
      SUBROUTINE GET_IREOTS(IARRAY,NACTOB)
      use lucia_data, only: IREOTS
      IMPLICIT NONE
      INTEGER NACTOB
      INTEGER IARRAY(NACTOB)
*
* Make the IREOTS reorder array available in ARRAY.
* Added in order to let MOLCAS use IREOTS.
*
      CALL ICOPY(NACTOB,IREOTS,1,IARRAY,1)

      END SUBROUTINE GET_IREOTS
