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
      Subroutine GugaCtl_dmrg()
*
      use gugx, only: NLEV, A0 => IA0, B0 => IB0, C0 => IC0,
     &                IFCAS, LV1RAS, LV3RAS, LM1RAS, LM3RAS
      Implicit Real*8 (A-H,O-Z)
*
#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "detdim.fh"
#include "spinfo_mclr.fh"
      Integer OrbSym(2*mxBas)
      Integer, Parameter:: iPrint=0
*
      Interface
      SUBROUTINE MKGUGA(NSM,NLEV,NSYM,STSYM,NCSF,Skip_MKSGNUM)
      IMPLICIT None

      Integer NLEV, NSYM, STSYM
      Integer NSM(NLEV)
      Integer NCSF(NSYM)
      Logical, Optional:: Skip_MKSGNUM
      End SUBROUTINE MKGUGA
      End Interface

*
      ntRas1=0
      ntRas2=0
      ntRas3=0
      Do iSym=1,nSym
         ntRas1=ntRas1+nRs1(iSym)
         ntRas2=ntRas2+nRs2(iSym)
         ntRas3=ntRas3+nRs3(iSym)
      End Do
*
      B0=iSpin-1
      A0=(nActEl-B0)/2
      C0=ntASh-A0-B0
      If ( (2*A0+B0).ne.nActEl ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' 2*A0+B0.ne.nActEl '
         Write (6,*)
      End If
      If ( A0.lt.0 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' A0.lt.0'
         Write (6,*)
      End If
      If ( B0.lt.0 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' B0.lt.0'
         Write (6,*)
      End If
      If ( C0.lt.0 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' C0.lt.0'
         Write (6,*)
      End If
*
      iOrb=0
      Do iSym=1,nSym
         Do iBas=1,nRs1(iSym)
            iOrb=iOrb+1
            OrbSym(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs2(iSym)
            iOrb=iOrb+1
            OrbSym(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs3(iSym)
            iOrb=iOrb+1
            OrbSym(iOrb)=iSym
         End Do
      End Do
*
      NLEV=ntASh
      LV1RAS=ntRas1
      LV3RAS=LV1RAS+ntRas2
      LM1RAS=2*LV1RAS-nHole1
      LM3RAS=nActEl-nElec3
!
      IFCAS=1
      Call mkGUGA(OrbSym,NLEV,NSYM,State_Sym,NCSF)
      NCONF=NCSF(State_Sym)

*
      If ( nConf.ne.NCSF(state_sym).and.(nConf.ne.1) ) then
         Write (6,*) "Set nConf=NCSF(state_sym)"
         Write (6,*)
         nConf=NCSF(state_sym)
      End If
*
      Call mkGUGA_Free()

      End Subroutine GugaCtl_dmrg
