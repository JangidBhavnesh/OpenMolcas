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
* Copyright (C) 1997, Anders Bernhardsson                              *
************************************************************************
      SubRoutine AddGrad2(rMat,fact)
      use Arrays, only: F0SQMO
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: ipCM, ipMat
      use input_mclr, only: nSym,nOrb
*
*     Purpose:
*             Adds the contribution from the gradient to
*              [2]
*             E   Kappa. This is done to insure us about
*             a beautiful convergence of the PCG,
*             which is just the case if E is symmetric.
*
      Implicit None
      Real*8 fact
      Real*8 rMat(*)
      Integer iS
      Real*8, Allocatable:: Temp(:)

      Do iS=1,nSym
        If (nOrb(is)*nOrb(is)==0) Cycle
        Call mma_allocate(Temp,nOrb(is)**2,Label='Temp')
*
*    T=Brillouin matrix
*
        Call DGeSub(F0SQMO(ipCM(is)),nOrb(is),'N',
     &              F0SQMO(ipCM(is)),nOrb(is),'T',
     &              Temp,nOrb(is),
     &              nOrb(is),nOrb(is))
*
*               t           t
*   +1/2 { Kappa T - T kappa  }
*
*
        Call DaXpY_(nOrb(is)**2,-2.0d0*Fact,Temp,1,
     &              rMat(ipMat(is,is)),1)
        Call mma_deallocate(Temp)
      End Do

      End SubRoutine AddGrad2
