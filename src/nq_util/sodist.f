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
      Subroutine SODist(SO_tmp,mAO,nCoor,mBas,mBas_Eff,
     &                  nCmp,nDeg,SOValue,nMOs,iAO)

      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 SO_tmp(mAO*nCoor,mBas,nCmp*nDeg),
     &       SOValue(mAO*nCoor,nMOs)
      Integer   iOff_MO(0:7), iOff_CMO(0:7)
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Character*80 Label
#endif
*
#ifdef _DEBUGPRINT_
      Write (6,*) 'SODist: MO-Coefficients'
      iOff=1
      Do iIrrep = 0, nIrrep-1
         If (nBas(iIrrep).gt.0) Then
            Write (6,*) ' Symmetry Block',iIrrep
            Call RecPrt(' ',' ',CMOs(iOff),nBas(iIrrep),nBas(iIrrep))
         End If
         iOff=iOff+nBas(iIrrep)**2
      End Do
#endif
*
*---- Compute some offsets
*
      itmp1=1
      Do iIrrep = 0, nIrrep-1
         iOff_MO(iIrrep)=itmp1
         itmp1=itmp1+nBas(iIrrep)
      End Do
*
      Do i1 = 1, nCmp
         iDeg=0
         Do iIrrep = 0, nIrrep-1
            iSO=iAOtSO(iAO+i1,iIrrep)
            If (iSO<0) Cycle

            iDeg=iDeg+1
            iOff=(i1-1)*nDeg+iDeg
!
!---------- Distribute contribution to all SO's in this irrep
!
            iMO =iOff_MO(iIrrep)
            Call DaXpY_(mAO*nCoor*mBas,One,SO_tmp(:,:,iOff:),1,
     &                  SOValue(:,iMO+iSO-1:),1)
          End Do
      End Do
*
#ifdef _DEBUGPRINT_
      Write (Label,'(A)')'SODist: SOValue(mAO*nCoor,nMOs)'
      Call RecPrt(Label,' ',SOValue(1,1),mAO*nCoor,nMOs)
#endif
*
      Return
      End
