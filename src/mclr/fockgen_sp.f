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
       SubRoutine FockGen_sp(d_0,rDens1,rdens2,Fock,fockout,idsym)
************************************************************************
*                                                                      *
*   Constructs active fockmatrix and Q matrix                          *
*                                                                      *
*   Input: rkappa : Rotation matrix                                    *
*          idsym  : symmetry of perturbation                           *
*                                                                      *
*                                                                      *
*   Output:MO     :MO integrals                                        *
*          Fock   :Fock matrix (one index transformed integrals)       *
*          MOtilde:MO (one index transformed integrals)                *
*                                                                      *
************************************************************************
      use Arrays, only: FIMO
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: nNA, ipCM, ipMat, nA, nDens2
      use input_mclr, only: nSym,nAsh,nIsh,nBas
      Implicit None
      Real*8 d_0
      Integer idSym
      Real*8 Fock(*),fockout(*),rdens2(*),rDens1(*)
!     Real*8 Fock(nDens2),fockout(*),rdens2(*),rDens1(nna*nna)
      Real*8, Allocatable:: MO(:), Scr(:)
      Integer n1, iS, n2, ipS, kS, jS, iB, jA, jAA, kA, kAA, ipM, ipF,
     &        iA, ip1, ip2
      Real*8 rd
*                                                                      *
************************************************************************
*                                                                      *
*
*     Coulomb term: F  =2(pk|ji)d
*                    kp          ij
*
      call dcopy_(nDens2,[0.0d0],0,Fock,1)
*
      n1=0
      Do iS = 1, nSym
         n1=Max(n1,nBas(iS))
      End Do
      n2=n1**2
      Call mma_allocate(MO,n2,Label='MO')
      Call mma_allocate(Scr,n2,Label='Scr')
*
      Do ips=1,nSym
         Do ks=1,nSym
            Do is=1,nSym
               jS=iEor(iEor(ips-1,ks-1),is-1)+1
*
*              Exchange term   F = -(pk|ji)d     (i>j)
*                               pl          kj
*
               If (iEOr(ips-1,iS-1)+1.eq.iDsym .and.
     &             nBas(ipS).gt.0) Then
                  Do iB=1,nIsh(iS)
                     Do jA=1,nAsh(jS)
                        jAA=jA+nIsh(js)
*
                        Call Coul(ipS,kS,iS,jS,iB,jAA,MO,Scr)
*
                        Do kA=1,nAsh(ks)
                           kAA=kA+nIsh(kS)
*
                           ipM=1+(kAA-1)*nBas(ipS)
                           ipF=ipMat(ipS,iS)+nBas(ipS)*(iB-1)
                           rd=rDens1(jA+nA(jS)+(kA+nA(ks)-1)*nna)
                           Call DaXpY_(nBas(ipS),-rd,
     &                                MO(ipM),1,Fock(ipF),1)
                        End Do
                     End Do
                  End Do
               End If
*
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Do iS=1,nSym
         If (nBas(iS).gt.0) Then
            jS=iEOr(is-1,iDSym-1)+1
            Do iA=1,nAsh(is)
               Do jA=1,nAsh(js)
                  rd=rDens1(iA+nA(iS)+(jA+nA(js)-1)*nna)
                  ip1=nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
                  ip2=nBas(iS)*(nIsh(js)+jA-1) +ipmat(is,js)
                 Call DaxPy_(nBAs(iS),Rd,FIMO(ip1),1,Fock(ip2),1)
               End Do
            End Do
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call CreQADD_sp(Fock,rdens2,idsym,MO,Scr,n2)

      Call mma_deallocate(Scr)
      Call mma_deallocate(MO)

*
      Do iS=1,nSym
         js=iEOR(is-1,idsym-1)+1
         If (nbas(is)*nBas(js).ne.0)
     &   Call DGESUB(Fock(ipMat(is,js)),nBas(is),'N',
     &               Fock(ipMat(js,is)),nBas(js),'T',
     %               FockOut(ipMat(is,js)),nBas(is),
     &               nBas(is),nBas(js))
      End Do
      Call DScal_(ndens2,2.0d0,FockOut,1)
      If (idsym.eq.1) Call Add2(Fockout,d_0)
*                                                                      *
************************************************************************
*                                                                      *
      End SubRoutine FockGen_sp
