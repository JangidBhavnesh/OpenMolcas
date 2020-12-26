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
      Subroutine Freq2(nIter,dqInt,nInter,Delta,qInt)
      use Slapaf_Info, only: Shift
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 dqInt(nInter,nIter), qInt(nInter,nIter+1)
*
      iRout = 183
      iPrint = nPrint(iRout)
*
*-----Compute new parameters for numerical differentiation.
      Call NwShft(Shift,nInter,dqInt,nIter,Delta,qInt)

      If (iPrint.gt.6) Then
         Write (6,*) ' Accumulate the gradient for yet one',
     &         ' parameter set'
         Write (6,*)
      End If
*
      Return
      End
