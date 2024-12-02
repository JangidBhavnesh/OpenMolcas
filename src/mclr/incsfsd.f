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
      Subroutine InCSFSD(iState,State_sym,GUGA)
      use Str_Info, only: CNSM
      use stdalloc, only: mma_allocate
      use MCLR_Data, only: iALLO,i1,iAnders,lConf,llDet
      use MCLR_Data, only: LuCSF2SD
      Implicit None
      Integer iState, State_sym
      Logical GUGA
      Integer idum(1)
      Integer iSym,iAdr,i,iad
*

*     Place pointer
*
      iSym=iEor(State_Sym-1,iState-1)+1
*
      If (isym.eq.1.and.i1.eq.1) Return
      If (isym.eq.iAnders) Return
*
      iAdr=2
      If (iSym.eq.1) iAdr=1
      iad=0
      Do i=1,iState-1
         Call iDafile(LUCSF2SD,0,idum,lldet,iad)
         Call iDafile(LUCSF2SD,0,idum,lconf,iad)
      End Do
*
      If (iSym.ne.1) Then
         If (iAnders.eq.-9)  Then
             Call mma_allocate(CNSM(2)%icts,lldet,Label='ICTS')
             Call mma_allocate(CNSM(2)%iconf,lConf,Label='ICONF')
             iAllo=1
         End If
         iAnders=isym
      End If
      If (iSym.eq.1) Then
          If (i1.eq.-9) Then
           Call mma_allocate(CNSM(1)%icts,lldet,Label='ICTS')
           Call mma_allocate(CNSM(1)%iconf,lConf,Label='ICONF')
           i1=1
          End If
      End If

!      open(unit=1422,file="det.index") ! yma
!      do i=1,lldet
!        write(1422,*)CNSM(iAdr)%icts(i)
!      end do
!      close(1422)

! calculated from zoo.f, the GUGA number for determinent
      Call iDafile(LUCSF2SD,2,CNSM(iAdr)%icts,lldet,iad)
      Call iDafile(LUCSF2SD,2,CNSM(iAdr)%iconf,lconf,iad)

c Avoid unused argument warnings
      If (.False.) Call Unused_logical(GUGA)
      End Subroutine InCSFSD
