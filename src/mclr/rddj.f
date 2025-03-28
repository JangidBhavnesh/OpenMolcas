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
* Copyright (C) 2000, Jonna Stalring                                   *
************************************************************************
*
      Subroutine rddj(G1r,G1Q,G2r,iestate)
*
* Jonna 000411
*
* Reads the one and two electron densities for estate
* and returns them in rectangular and single triangular storage
*
*
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: LuJob
      use input_mclr, only: ntAsh,iTOC
      Implicit None
      Real*8 G1r(*), G1Q(*),G2r(*)
      integer iestate

#include "SysDef.fh"
      Real*8, Allocatable:: G2Q(:)
      Real*8 rdum(1)
      Integer nG1,nG2,iR,jDisk,i,j,iB,jB,iDij,iRij,kB,lB,iDkl,iRkl,
     &        iIJKL, iRijkl
      Real*8 Fact
*
      Integer itri
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

*
      ng1=itri(ntash,ntash)
      ng2=itri(ng1,ng1)
*
      Call mma_allocate(G2Q,ng2,Label='G2Q')
c
c     Read one and two el dens for state iestate
c
      iR=iestate
      jdisk=itoc(3)
      Do i=1,iR-1
         Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
         Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
      End Do
      Call dDaFile(LUJOB ,2,G1q,ng1,jDisk)
      Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
      Call dDaFile(LUJOB ,2,G2Q,Ng2,jDisk)
      Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
c
c Make one el rectangular and two el singel triang.
c
      Do iB=1,ntash
         Do jB=1,ntash
            iDij=iTri(ib,jB)
            iRij=jb+(ib-1)*ntash
            Do kB=1,ntash
               Do lB=1,ntash
                  iDkl=iTri(kB,lB)
                  iRkl=lb+(kb-1)*ntash
                  fact=1.0d00
                  If(iDij.ge.iDkl .and. kB.eq.lB) fact=2.0d00
                  If(iDij.lt.iDkl .and. iB.eq.jB) fact=2.0d00
                  iijkl=itri(iDij,iDkl)
                  iRijkl=itri(iRij,iRkl)
                  G2R(iRijkl)=Fact*G2Q(iijkl)
               End Do
            End Do
         End Do
      End Do
      Do iB=1,ntash
         Do jB=1,ntash
            G1R(+ib+(jb-1)*ntash)= g1q(+itri(ib,jb))
         End Do
      End Do
*
      Call mma_deallocate(G2Q)
*
      End Subroutine rddj
