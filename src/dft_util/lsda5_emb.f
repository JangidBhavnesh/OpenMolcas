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
      Subroutine LSDA5_emb(mGrid,Rho,nRho,
     &                     nDmat,F_xc)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
************************************************************************
      use OFembed, only: KEonly
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "hflda.fh"
      Real*8 Rho(nRho,mGrid),
     &       F_xc(mGrid)
*
************************************************************************
*
*---- Thomas-Fermi Kinetic energy functional
*
      Coeff=One
      Call TF_Ts(mGrid,nDmat,F_xc,Coeff)

      If (KEonly) Return
*
*---- LSDA5 for exchange-correlation energy functional
*
      Call LSDA5(mGrid,Rho,nRho,nDmat,F_xc)

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
