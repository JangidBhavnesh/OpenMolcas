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
      SUBROUTINE diag_master()
      use GLBBAS, only: INT1, INT1O
      use CandS, only: ISSM
      use lucia_data, only: IPRDIA
*
*  To do in this subroutine:
*
*     - Make sure all calling parameters are accounted for
*     - Make sure mxntts is set in previous subroutine
*     - Make sure nsmst is set in previous subroutine
*
*  Set up the diagonal for the CI calculation
*
      implicit None
#include "mxpdim.fh"
#include "orbinp.fh"
      REAL*8 EREF
*
      INT1(:)=INT1O(:)
      CALL GASCI(ISSM, 1, IPRDIA, EREF, 0, 0)
*
      END SUBROUTINE diag_master
*
