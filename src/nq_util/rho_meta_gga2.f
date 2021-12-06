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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine Rho_meta_GGA2(nD,mGrid,
     &                         list_s,nlist_s,TabAO,ipTabAO,mAO,nTabAO,
     &                         Fact,mdc,list_bas,Index,nIndex)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
************************************************************************
      use iSD_data
      use k2_arrays, only: DeDe, ipDijS
      use nq_Grid, only: Rho, Sigma, Lapl, Tau
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
      Integer list_s(2,nlist_s), ipTabAO(nlist_s), list_bas(2,nlist_s),
     &        Index(nIndex)
      Real*8 Fact(mdc**2), TabAO(nTabAO)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (Debug) Then
         Write (6,*) 'mAO=',mAO
         Write (6,*) 'mGrid=',mGrid
         Write (6,*) 'nTabAO=',nTabAO
         Write (6,*) 'nlist_s=',nlist_s
         Do iList_s = 1, nList_s
            Write (6,*) 'iList_s=',iList_s
            iSkal = list_s(1,ilist_s)
            iCmp  = iSD( 2,iSkal)
            iBas  = iSD( 3,iSkal)
            mTabAO=iBas*iCmp
            Call RecPrt('Rho_meta_GGA: TabAO',' ',
     &                  TabAO(ipTabAO(iList_s)),mAO,mGrid*mTabAO)
         End Do
      End If
#endif
*
      Rho(:,1:mGrid)=Zero
      Sigma(:,1:mGrid)=Zero
      Tau(:,1:mGrid)=Zero
      Lapl(:,1:mGrid)=Zero
*                                                                      *
************************************************************************
*                                                                      *
      Do ilist_s=1,nlist_s
         iSkal = list_s(1,ilist_s)
         iCmp  = iSD( 2,iSkal)
         iBas  = iSD( 3,iSkal)
         iBas_Eff=list_bas(1,ilist_s)
         kDCRE=list_s(2,ilist_s)
         index_i=list_bas(2,ilist_s)
         mdci  = iSD(10,iSkal)
         iShell= iSD(11,iSkal)
         nFunc_i=iBas*iCmp
*
         Do jlist_s=1,ilist_s
            Fij=Two
            If (jlist_s.eq.ilist_s) Fij=One
            jSkal = list_s(1,jlist_s)
            kDCRR=list_s(2,jlist_s)
            jCmp  = iSD( 2,jSkal)
            jBas  = iSD( 3,jSkal)
            jBas_Eff=list_bas(1,jlist_s)
            index_j =list_bas(2,jlist_s)
            mdcj  = iSD(10,jSkal)
            jShell= iSD(11,jSkal)
            nFunc_j=jBas*jCmp
*
            mDij=nFunc_i*nFunc_j
*
*---------- Get the Density
*
            ijS=iTri(iShell,jShell)
            ip_Tmp=ipDijs
            Call Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ip_Tmp,nD)
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*)
               Write (6,*) 'iS,jS=',iSkal,jSkal
               Write (6,*) 'mDCRij,mDij=',mDCRij,mDij
               Write (6,*) 'ipDij,ipDSij,ipDDij=',ipDij,ipDSij,ipDDij
            End If
#endif
*
            ij = (mdcj-1)*mdc + mdci
*
            iER=iEOr(kDCRE,kDCRR)
            lDCRER=NrOpr(iER)
*
            ip_D_a=ipDij+lDCRER*mDij
            ip_D_b=ip_D_a
            If (nD.ne.1) ip_D_b=ipDSij+lDCRER*mDij
*
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*) 'Rho_meta_GGA2'
               nBB = iBas*jBas
               nCC = iCmp*jCmp
               Write (6,*) 'iShell,jshell=', iShell,jshell
               Write (6,*) 'kDCRE,kDCRR=', kDCRE,kDCRR
               Write (6,*) 'iER,lDCRER=',iER,lDCRER
               Call RecPrt('DAij',' ',DeDe(ip_D_a),nBB,nCC)
               If (nD.ne.1)
     &            Call RecPrt('DBij',' ',DeDe(ip_D_b),nBB,nCC)
            End If
#endif
*
            If (nD.eq.1) Then
               If (iShell.ge.jShell) Then
               Call Do_Rho9a(mGrid,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Fij,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho9a(mGrid,
     &                       DeDe(ip_D_a),                  mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Fij,
     &                       Index(index_i),Index(index_j))
               End If
            Else
               If (iShell.ge.jShell) Then
               Call Do_Rho9_(mGrid,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       Fact(ij)*Fij,
     &                       Index(index_i),Index(index_j))
               Else
               Call Do_Rho9_(mGrid,
     &                       DeDe(ip_D_a),DeDe(ip_D_b),     mAO,
     &                       TabAO(ipTabAO(jList_s)),jBas,jBas_Eff,jCmp,
     &                       TabAO(ipTabAO(iList_s)),iBas,iBas_Eff,iCmp,
     &                       Fact(ij)*Fij,
     &                       Index(index_i),Index(index_j))
               End If
            End If
*
         End Do                      ! jlist_s
      End Do                         ! ilist_s
*
      Return
      End
      Subroutine Do_Rho9a(mGrid,
     &                    DAij,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                        TabAO2,jBas,jBas_Eff,jCmp,
     &                    Fact,Index_i,Index_j)
      use nq_Grid, only: Rho, Sigma, Lapl, Tau
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 DAij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer Index_i(iBas_Eff*iCmp), Index_j(jBas_Eff*jCmp)
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB = Index_j(jCB_Eff)
*

         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB = Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
*
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_12=TabAO1(1,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_13=TabAO1(1,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_14=TabAO1(1,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
*
               Prod_22=TabAO1(2,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_33=TabAO1(3,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_44=TabAO1(4,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
               Prod_1X= TabAO1( 1,iGrid,iCB_Eff) *
     &                 (TabAO2( 5,iGrid,jCB_Eff)
     &                 +TabAO2( 8,iGrid,jCB_Eff)
     &                 +TabAO2(10,iGrid,jCB_Eff))
               Prod_X1=(TabAO1( 5,iGrid,iCB_Eff)
     &                 +TabAO1( 8,iGrid,iCB_Eff)
     &                 +TabAO1(10,iGrid,iCB_Eff))
     &                 *TabAO2( 1,iGrid,jCB_Eff)
*
               Rho(1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho(2,iGrid)=Rho(2,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho(3,iGrid)=Rho(3,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho(4,iGrid)=Rho(4,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho(5,iGrid)=Rho(5,iGrid)
     &                     + (Prod_22+Prod_33+Prod_44)*DAij_
               Rho(6,iGrid)=Rho(6,iGrid) + (Prod_1X
     &                     + Two*(Prod_22+Prod_33+Prod_44)
     &                                   +  Prod_X1) *DAij_

              RhoX=(Prod_21+Prod_12)*DAij_
              RhoY=(Prod_31+Prod_13)*DAij_
              RhoZ=(Prod_41+Prod_14)*DAij_
              Sigma(1,iGrid)=Sigma(1,iGrid)
     &                      +RhoX**2
     &                      +RhoY**2
     &                      +RhoZ**2
              Tau(1,iGrid)=Tau(1,iGrid)
     &                     +(Prod_22+Prod_33+Prod_44)*DAij_
              Lapl(1,iGrid)=Lapl(1,iGrid)
     &                    + (Prod_1X
     &                    + Two*(Prod_22+Prod_33+Prod_44)
     &                    + Prod_X1) *DAij_
            End Do    ! iGrid
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
      Subroutine Do_Rho9_(mGrid,
     &                    DAij,DBij,
     &                    mAO,TabAO1,iBas,iBas_Eff,iCmp,
     &                        TabAO2,jBas,jBas_Eff,jCmp,
     &                    Fact,Index_i,Index_j)
      use nq_Grid, only: Rho, Sigma, Lapl, Tau
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 DAij(iBas*iCmp,jBas*jCmp), DBij(iBas*iCmp,jBas*jCmp),
     &       TabAO1(mAO,mGrid,iBas_Eff*iCmp),
     &       TabAO2(mAO,mGrid,jBas_Eff*jCmp)
      Integer Index_i(iBas_Eff*iCmp), Index_j(jBas_Eff*jCmp)
*                                                                      *
************************************************************************
*                                                                      *
      Do jCB_Eff = 1, jBas_Eff*jCmp
         jCB = Index_j(jCB_Eff)
*
         Do iCB_Eff = 1, iBas_Eff*iCmp
            iCB = Index_i(iCB_Eff)
*
            DAij_=DAij(iCB,jCB)*Fact
            DBij_=DBij(iCB,jCB)*Fact
*
            Do iGrid = 1, mGrid
               Prod_11=TabAO1(1,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_21=TabAO1(2,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_12=TabAO1(1,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_31=TabAO1(3,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_13=TabAO1(1,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_41=TabAO1(4,iGrid,iCB_Eff)*TabAO2(1,iGrid,jCB_Eff)
               Prod_14=TabAO1(1,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
*
               Prod_22=TabAO1(2,iGrid,iCB_Eff)*TabAO2(2,iGrid,jCB_Eff)
               Prod_33=TabAO1(3,iGrid,iCB_Eff)*TabAO2(3,iGrid,jCB_Eff)
               Prod_44=TabAO1(4,iGrid,iCB_Eff)*TabAO2(4,iGrid,jCB_Eff)
               Prod_1X= TabAO1( 1,iGrid,iCB_Eff) *
     &                 (TabAO2( 5,iGrid,jCB_Eff)
     &                 +TabAO2( 8,iGrid,jCB_Eff)
     &                 +TabAO2(10,iGrid,jCB_Eff))
               Prod_X1=(TabAO1( 5,iGrid,iCB_Eff)
     &                 +TabAO1( 8,iGrid,iCB_Eff)
     &                 +TabAO1(10,iGrid,iCB_Eff))
     &                 *TabAO2( 1,iGrid,jCB_Eff)
*
               Rho( 1,iGrid)=Rho(1,iGrid) +     Prod_11      *DAij_
               Rho( 2,iGrid)=Rho(2,iGrid) +     Prod_11      *DBij_
               Rho( 3,iGrid)=Rho(3,iGrid) + (Prod_21+Prod_12)*DAij_
               Rho( 4,iGrid)=Rho(4,iGrid) + (Prod_31+Prod_13)*DAij_
               Rho( 5,iGrid)=Rho(5,iGrid) + (Prod_41+Prod_14)*DAij_
               Rho( 6,iGrid)=Rho(6,iGrid) + (Prod_21+Prod_12)*DBij_
               Rho( 7,iGrid)=Rho(7,iGrid) + (Prod_31+Prod_13)*DBij_
               Rho( 8,iGrid)=Rho(8,iGrid) + (Prod_41+Prod_14)*DBij_
               Rho( 9,iGrid)=Rho( 9,iGrid)
     &                      +(Prod_22+Prod_33+Prod_44)*DAij_
               Rho(10,iGrid)=Rho(10,iGrid)
     &                      +(Prod_22+Prod_33+Prod_44)*DBij_
               Rho(11,iGrid)=Rho(11,iGrid) + (Prod_1X
     &                     + Two*(Prod_22+Prod_33+Prod_44)
     &                                   +  Prod_X1) *DAij_
               Rho(12,iGrid)=Rho(12,iGrid) + (Prod_1X
     &                     + Two*(Prod_22+Prod_33+Prod_44)
     &                                   +  Prod_X1) *DBij_

              RhoX_A=(Prod_21+Prod_12)*DAij_
              RhoY_A=(Prod_31+Prod_13)*DAij_
              RhoZ_A=(Prod_41+Prod_14)*DAij_
              RhoX_A=(Prod_21+Prod_12)*DBij_
              RhoY_A=(Prod_31+Prod_13)*DBij_
              RhoZ_A=(Prod_41+Prod_14)*DBij_
              Sigma(1,iGrid)=Sigma(1,iGrid)
     &                      +RhoX_A**2
     &                      +RhoY_A**2
     &                      +RhoZ_A**2
              Sigma(2,iGrid)=Sigma(2,iGrid)
     &                      +RhoX_A*RhoX_A
     &                      +RhoY_A*RhoY_A
     &                      +RhoZ_A*RhoZ_A
              Sigma(3,iGrid)=Sigma(3,iGrid)
     &                      +RhoX_B**2
     &                      +RhoY_B**2
     &                      +RhoZ_B**2
              Tau(1,iGrid)=Tau(1,iGrid)
     &                     +(Prod_22+Prod_33+Prod_44)*DAij_
              Tau(2,iGrid)=Tau(2,iGrid)
     &                     +(Prod_22+Prod_33+Prod_44)*DBij_
              Lapl(1,iGrid)=Lapl(1,iGrid)
     &                    + (Prod_1X
     &                    + Two*(Prod_22+Prod_33+Prod_44)
     &                    + Prod_X1) *DAij_
              Lapl(2,iGrid)=Lapl(2,iGrid)
     &                    + (Prod_1X
     &                    + Two*(Prod_22+Prod_33+Prod_44)
     &                    + Prod_X1) *DBij_
            End Do    ! iGrid
*
         End Do          ! iCB
      End Do             ! jCB
*
      Return
      End
