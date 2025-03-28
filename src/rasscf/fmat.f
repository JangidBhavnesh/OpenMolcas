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
* Copyright (C) 1989, Bjorn O. Roos                                    *
*               1989, Per Ake Malmqvist                                *
*               1991,1993,1996, Markus P. Fuelscher                    *
************************************************************************
      Subroutine Fmat(CMO,PUVX,D,D1A,FI,FA)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Update the Fock matrix for the active orbitals and transform     *
*     it to MO basis as well as the matrix FI (Fock matrix) for        *
*     frozen and inactive orbitals).                                   *
*                                                                      *
*     calling arguments:                                               *
*     CMO     : array of real*8                                        *
*               MO-coefficients                                        *
*     PUVX    : array of real*8                                        *
*               two-electron integrals (pu!vx)                         *
*     D       : array of real*8                                        *
*               averaged one-body density matrix                       *
*     D1A     : array of real*8                                        *
*               active one body density matrix in AO-basis             *
*     FI      : array of real*8                                        *
*               inactive Fock matrix. In input is in AO basis.         *
*               In output in MO basis.                                 *
*     FA      : array of real*8                                        *
*               active Fock matrix. In input in AO Basis.              *
*               In output in MO basis. It is also modified by scaling  *
*               exchange part in ExFac .ne. 1.0d0                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     B.O. Roos and P.Aa. Malmqvist                                    *
*     University of Lund, Sweden, 1989                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     - updated for MOLCAS version 2                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1991               *
*     - updated for MOLCAS version 3                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1993               *
*     - updated for integral direct and reaction field calculations    *
*       M.P. Fuelscher, University of Lund, Sweden, 1996               *
*                                                                      *
************************************************************************

      Use RunFile_procedures, Only: Get_dExcdRa
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_global, only: KSDFT, DFTFOCK, ECAS, EMY, ExFac, NAC,
     &                         NewFock, nFint, VIA, VIA_DFT, l_casdft
      use printlevel, only: DEBUG
      use output_ras, only: LF,IPRLOC
      use general_data, only: NSYM,NTOT1,NASH,NBAS,NFRO,NISH,NORB

      Implicit None
      Real*8 CMO(*) , PUVX(*) , D(*) , D1A(*) , FI(*) , FA(*)

      Character(LEN=16), Parameter :: ROUTINE='FMAT    '

      Real*8, Allocatable :: TmpFck(:), Tmp1(:), Tmp2(:), TmpD1A(:)
      Integer iPrLev, iOff, iSym, iBas, i, iFro, ij, iOff1, iOff2,
     &        iOff3, iOrb, ipTmpFck, ipTmpFckA, ipTmpFckI, j, jOrb,
     &        nTmpFck
      Real*8, External:: DDot_

C Local print level (if any)
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE

        write(6,*) repeat('*',65)
        write(6,*) 'Entering FMAT routine called by SXCTL!'
        write(6,*) repeat('*',65)
        write(6,*) 'printing input matrices :'
        write(6,*) repeat('*',65)
        Write(LF,*)
        Write(LF,*) ' CMOs in FMAT'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
         iOff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*)(CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do

         Write(LF,*)
         Write(LF,*) ' PUVX in FMAT'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         call wrtmat(PUVX,1,nFint, 1, nFint)

         Write(LF,*)
         Write(LF,*) ' ---------------------'
        CALL TRIPRT('Averaged one-body density matrix D, in MO in FMAT',
     &              ' ',D,NAC)

         Write(LF,*)
         Write(LF,*) ' D1A in AO basis in FMAT'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         iOff = 1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          call wrtmat(D1A(iOff),iBas,iBas, iBas, iBas)
          iOff = iOff + iBas*iBas
         End DO

         Write(LF,*)
         Write(LF,*) ' FI in AO-basis in FMAT'
         Write(LF,*) ' --------------'
         Write(LF,*)
         iOff = 1
         Do iSym = 1,nSym
           iOrb = nOrb(iSym)
           Call TriPrt(' ',' ',FI(iOff),iOrb)
           iOff = iOff + (iOrb*iOrb+iOrb)/2
         End Do

         Write(LF,*)
         Write(LF,*) ' FA in AO-basis in FMAT'
         Write(LF,*) ' --------------'
         Write(LF,*)
         iOff = 1
         Do iSym = 1,nSym
           iOrb = nOrb(iSym)
           Call TriPrt(' ',' ',FA(iOff),iOrb)
           iOff = iOff + (iOrb*iOrb+iOrb)/2
         End Do
       End If

*     create FA in AO basis
      Call mma_allocate(Tmp1,nTot1,Label='Tmp1')
      Call Fold(nSym,nBas,D1A,Tmp1)
      If(KSDFT.ne.'SCF') NewFock=0
c      If (NewFock.eq.0) Then
c         nBMX=0
c         Do iSym=1,nSym
c            nBMX=Max(nBMX,nBas(iSym))
c         End Do
c         Call FZero(FA,nTot1)
c         Call FTwo_Drv(nSym,nBas,nAsh,nSkipX,
c     &                    Tmp1,D1A,FA,nTot1,
c     &                    ExFac,nBMX,CMO)
c      End If

*     Inactive-active contribution to ECAS
      VIA=dDot_(nTot1,FI,1,Tmp1,1)
      ECAS=EMY+VIA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*) ' Total core energy fmat:       ',EMY
        Write(LF,*) ' inactive-active interaction:  ',VIA
        Write(LF,*) ' CAS energy (core+interaction):',ECAS
      End If
      Call mma_deallocate(Tmp1)

*     print FI and FA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI in AO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FI(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
        Write(LF,*)
        Write(LF,*) ' FA in AO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If

*     transform FI from AO to MO basis
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        If (iBas==0) Cycle
        iOrb = nOrb(iSym)
        If (iOrb==0) Cycle
        iFro = nFro(iSym)
        Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
        Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
        Call Square(FI(iOff1),Tmp1,1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Tmp1,iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Tmp2,iBas)
       Call DGEMM_Tri('T','N',iOrb,iOrb,iBas,
     &                 1.0D0,Tmp2,iBas,
     &                       CMO(iOff2+(iFro*iBas)),iBas,
     &                 0.0D0,FI(iOff3),iOrb)
        Call mma_deallocate(Tmp2)
        Call mma_deallocate(Tmp1)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

*     transform FA from AO to MO basis
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        If (iBas==0) Cycle
        iOrb = nOrb(iSym)
        If (iOrb==0) Cycle
        iFro = nFro(iSym)
        Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
        Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
        Call Square(FA(iOff1),Tmp1,1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Tmp1,iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Tmp2,iBas)
        Call DGEMM_Tri('T','N',iOrb,iOrb,iBas,
     &                 1.0D0,Tmp2,iBas,
     &                       CMO(iOff2+(iFro*iBas)),iBas,
     &                 0.0D0,FA(iOff3),iOrb)
        Call mma_deallocate(Tmp2)
        Call mma_deallocate(Tmp1)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

c**************************************************************************
c              Add DFT part to Fock matrix:                               *
c**************************************************************************
      If(KSDFT(1:3).ne.'SCF'.and.KSDFT(1:3).ne.'PAM'.and.
     &      .not. l_casdft ) Then
        ipTmpFckI=-99999
        ipTmpFckA=-99999
        Call Get_dExcdRa(TmpFck,nTmpFck)
        ipTmpFck = 1
        If(nTmpFck.eq.NTOT1) Then
           ipTmpFckI=ipTmpFck
        Else If(nTmpFck.eq.2*NTOT1) Then
           ipTmpFckI=ipTmpFck
           ipTmpFckA=ipTmpFck+nTot1
        Else
           Write(LF,*) ' Somethings wrong in dim. DFT',nTmpFck
           Call Abend()
        End If
        Call mma_allocate(TmpD1A,nTot1,Label='TmpD1A')
        Call Fold(nSym,nBas,D1A,TmpD1A)
        VIA_DFT=dDot_(nTot1,TmpFck(ipTmpFckI),1,TmpD1A,1)
        Call mma_deallocate(TmpD1A)
*
*          Transform alpha density from AO to MO
*
        iOff1 = 1
        iOff2 = 1
        iOff3 = 1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          If (iBas==0) Cycle
          iOrb = nOrb(iSym)
        If (iOrb==0) Cycle
          iFro = nFro(iSym)
          Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
          Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
          Call Square(TmpFck(ipTmpFckI+iOff1-1),
     &                Tmp1,1,iBas,iBas)
          Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Tmp1,iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Tmp2,iBas)
          Call DGEMM_Tri('T','N',iOrb,iOrb,iBas,
     &                   1.0D0,Tmp2,iBas,
     &                         CMO(iOff2+(iFro*iBas)),iBas,
     &                   0.0D0,TmpFck(ipTmpFckI+iOff3-1),iOrb)
          Call mma_deallocate(Tmp2)
          Call mma_deallocate(Tmp1)
          iOff1 = iOff1 + (iBas*iBas+iBas)/2
          iOff2 = iOff2 + iBas*iBas
          iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
        End Do
*
*          Transform Active DFT Fock from AO to MO
*
        If(ipTmpFckA.ne.-99999) Then
        iOff1 = 1
        iOff2 = 1
        iOff3 = 1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          If (iBas==0) Cycle
          iOrb = nOrb(iSym)
          If (iOrb==0) Cycle
          iFro = nFro(iSym)
          Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
          Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
          Call Square(TmpFck(ipTmpFckA+iOff1-1),
     &                Tmp1,1,iBas,iBas)
          Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Tmp1,iBas,
     &               CMO(iOff2+(iFro*iBas)),iBas,
     &               0.0d0,Tmp2,iBas)
          Call DGEMM_Tri('T','N',iOrb,iOrb,iBas,
     &                   1.0D0,Tmp2,iBas,
     &                         CMO(iOff2+(iFro*iBas)),iBas,
     &                   0.0D0,TmpFck(ipTmpFckA+iOff3-1),iOrb)
          Call mma_deallocate(Tmp2)
          Call mma_deallocate(Tmp1)
          iOff1 = iOff1 + (iBas*iBas+iBas)/2
          iOff2 = iOff2 + iBas*iBas
          iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
        End Do
        End If
*
c        If(DFTFOCK(1:4).ne.'ROKS') Then
c          Write(LF,*) ' Just add a,b to FA,FI',DFTFOCK(1:4)
c        Else
c          Write(LF,*) ' ROKS formula',DFTFOCK(1:4)
c        End If
*
        If(DFTFOCK(1:4).ne.'ROKS') Then
          call daxpy_(NTOT1,1.0D0,TmpFck(ipTmpFckI),1,FI,1)
          If(ipTmpFckA.ne.-99999)
     &    call daxpy_(NTOT1,1.0D0,TmpFck(ipTmpFckA),1,FA,1)
        Else If (DFTFOCK(1:4).eq.'ROKS') Then
           iOff1 = 0
           Do iSym = 1,nSym
              Do iOrb=1,nOrb(iSym)
                 Do jOrb=1,iOrb
                    ij=iOff1+iOrb*(iOrb-1)/2+jOrb
                    If(iOrb.le.nIsh(iSym)) Then
                      FI(ij)=FI(ij)+0.5d0*
     &                  (TmpFck(ipTmpFckI+ij-1)+TmpFck(ipTmpFckA+ij-1))
                    End If
                    If (iOrb.gt.nIsh(iSym).and.
     &                  iOrb.le.nIsh(iSym)+nAsh(iSym)) Then
                       If (jOrb.le.nIsh(iSym)) Then
                          FI(ij)=FI(ij)+TmpFck(ipTmpFckA+ij-1)
                       Else
                          FI(ij)=FI(ij)+0.5d0*(TmpFck(ipTmpFckI+ij-1)+
     &                                         TmpFck(ipTmpFckA+ij-1))
                       End If
                    End If
                    If (iOrb.gt.nIsh(iSym)+nAsh(iSym)) Then
                       If(jOrb.gt.nIsh(iSym).and.
     &                    jOrb.le.nIsh(iSym)+nAsh(iSym)) Then
                          FI(ij)=FI(ij)+TmpFck(ipTmpFckI+ij-1)
                       Else
                          FI(ij)=FI(ij)+0.5d0*(TmpFck(ipTmpFckI+ij-1)+
     &                                         TmpFck(ipTmpFckA+ij-1))
                       End If

                    End If
                 End Do
              End Do
              iOff1 = iOff1 + (nOrb(iSym)*nOrb(iSym)+nOrb(iSym))/2
           End Do
        Else
           Write(LF,*) " Not implemented yet"
        End If
        Call mma_deallocate(TmpFck)
      End If
***************************************************************************
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FA in MO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If
*     update Fock matrix by rescaling exchange term...
      If (NewFock.eq.1) Call Upd_FA(PUVX,FA,D,ExFac)

*     print FI and FA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI in MO-basis in fmat'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FI(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FA in MO-basis in fmat after upd_FA'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If

      End Subroutine Fmat
