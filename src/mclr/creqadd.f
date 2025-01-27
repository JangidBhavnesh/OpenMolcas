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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine creqadd(q,G2,idSym,MO,Scr,n2)
*
*      Constructs the Q matrix
*
      use MCLR_Data, only: nDens2, nNA, ipMat, nA
      use input_mclr, only: nSym,nAsh,nIsh,nOrb
      Implicit None
      Integer idSym, n2
      Real*8 Q(nDens2),G2(*), MO(n2), Scr(n2)

      Integer iS, jS, kS, lS, ipS, ijS, iAsh, jAsh, kAsh, lAsh, kAA,
     &         lAA, iAA, ikl, ipQ, ipM, iij, ipG
      Real*8 P_ijkl
*                                                                      *
************************************************************************
*                                                                      *
       Integer i,j,itri
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*      Q = (pj|kl)P
*       pi         ijkl
*                                                                      *
************************************************************************
*                                                                      *
       Do iS=1,nSym
          ipS=iEOr(is-1,idsym-1)+1
          if (norb(ips).ne.0) Then
             Do jS=1,nsym
                ijS=iEOR(is-1,js-1)+1
                Do kS=1,nSym
                   ls=iEOr(ijs-1,iEor(ks-1,idsym-1))+1
*                                                                      *
************************************************************************
*                                                                      *
                   Do kAsh=1,nAsh(kS)
                      kAA=kAsh+nIsh(kS)
                      Do lAsh=1,nAsh(lS)
                         lAA=lAsh+nIsh(lS)
                         ikl=nna*(lAsh+nA(lS)-1)+kAsh+nA(kS)
*
*                        Pick up (pj|kl)
*
                         Call Coul(ipS,jS,kS,lS,kAA,lAA,MO,Scr)
*
                         Do iAsh=1,nAsh(iS)
                            iAA=iAsh+nIsh(iS)
                            ipQ=ipMat(ipS,iS)+nOrb(ipS)*(iAA-1)
                            ipM=1+nIsh(jS)*nOrb(ipS)
                            Do jAsh=1,nAsh(jS)
                               iij=nna*(jAsh+nA(jS)-1)+iAsh+nA(iS)
                               ipG=itri(iij,ikl)
                               P_ijkl=G2(ipG)
*
                               Call DaXpY_(nOrb(ipS),P_ijkl,MO(ipM),1,
     &                                     Q(ipQ),1)
                               ipM=ipM+nOrb(ipS)
*
                            End Do
                         End Do
*
                      End Do
                   End Do
*                                                                      *
************************************************************************
*                                                                      *
                End Do  ! kS
             End Do     ! jS
          End If
       End Do           ! iS
*                                                                      *
************************************************************************
*                                                                      *
       End SubRoutine creqadd
