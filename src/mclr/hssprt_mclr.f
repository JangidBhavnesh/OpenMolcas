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
      SubRoutine HssPrt_MCLR(ideg,Hess,ldisp)
      use stdalloc, only: mma_allocate, mma_deallocate
      use input_mclr, only: nIrrep,nSym,ChIrr
      Implicit None
      integer ideg(*)
      Real*8     Hess(*)
      Integer  ldisp(nsym)
*     local variables
      Integer  kDisp(8)
      Character Title*39
      Real*8, Allocatable:: Temp(:)
      Integer i,iaa,ii,iIrrep,j,jj
*
      integer idisp,jdisp,Ind
      Ind(idisp,jdisp)=idisp*(idisp-1)/2+jdisp
*
      iDisp=0
      Do iIrrep=1,nIrrep
           kDisp(iIrrep)=iDisp
           iDisp=iDisp+lDisp(iIrrep)
           Write (6,*) lDisp(iIrrep)
      End Do
*
      Call mma_allocate(Temp,iDisp**2,Label='Temp')
      iaa=0
      Do iIrrep=1,nIrrep
        If (ldisp(iirrep).ne.0) Then
        Write(title,'(A,A)') 'Hessian in Irrep ',
     &                chirr(iIrrep)

        Do i=1,lDisp(iirrep)
         Do j=1,i
          ii=ind(i,j)
          jj=iaa+ind(i,j)
          Temp(ii)=Hess(jj)*
     &        sqrt(DBLE(ideg(i+kdisp(iirrep))*
     &                    ideg(j+kdisp(iirrep)) ))
         End Do
        End Do
        Call TriPrt(title,' ',Temp,ldisp(iirrep))
        iaa=iaa+ind(ldisp(iirrep),ldisp(iirrep))
        End If
       End Do
       Call mma_deallocate(Temp)

      End SubRoutine HssPrt_MCLR
