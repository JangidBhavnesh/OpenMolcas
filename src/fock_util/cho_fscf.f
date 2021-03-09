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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************

      SUBROUTINE CHO_FSCF(rc,nDen,ipFLT,nForb,nIorb,Porb,ipDLT,ExFac)

**********************************************************************
*  Author : F. Aquilante
*
C *************** INACTIVE AO-BASIS FOCK MATRIX **********************
C
C   F(ab) = sum_J  Lab,J * V(J)  -  sum_Jk  Lka,J * Lkb,J
C
**********************************************************************
C
C      V(J) = sum_gd  Lgd,J * Dtot(gd)
C
C      a,b,g,d:  AO-index
C      k:        MO-index   belonging to (Frozen+Inactive)
C
**********************************************************************
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      use Data_structures, only: CMO_Type, Laq_Type
      use Data_structures, only: Allocate_Laq, Deallocate_Laq
      Implicit Real*8 (a-h,o-z)

      Integer   rc,nDen
      Integer   iSkip(8)
      Integer   ISTLT(8)
      Real*8    tread(2),tcoul(2),texch(2)
      Real*8    FactCI,FactXI,ExFac
      Integer   ipDLT(nDen),ipFLT(nDen)
      Type (CMO_Type)   Porb(nDen)
      Integer   nForb(8,nDen),nIorb(8,nDen)
#ifdef _DEBUGPRINT_
      Logical   Debug
#endif
      Logical   DoRead
      Character*50 CFmt
      Character(LEN=8), Parameter:: SECNAM = 'CHO_FSCF'
#include "chotime.fh"

#include "real.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Real*8, Parameter:: xone = -one

      Logical add
      Character*6 mode

      Real*8, Allocatable:: Lrs(:,:)
      Real*8, Allocatable:: VJ(:)
      Integer:: nAux(8)

      Type (Laq_type) Laq(2)

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
************************************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in CASSCF-debug
#endif

      FactCI = one
      FactXI = xone*ExFac

      DoRead  = .false.
      IREDC= -1  ! unknwn reduced set

      If (nDen.ne.1 .and. nDen.ne.2) then
         write(6,*)SECNAM//'Invalid parameter nDen= ',nDen
         call abend()
      EndIf


        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        ! 1 --> CPU   2 --> Wall
        tread(:) = zero  !time read/transform vectors
        tcoul(:) = zero  !time for computing Coulomb
        texch(:) = zero  !time for computing Exchange

C ==================================================================

c --- Various offsets
c --------------------
        ISTLT(1)=0
      DO ISYM=2,NSYM
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB ! Inactive D and F matrices
      END DO

      iLoc = 3 ! use scratch location in reduced index arrays

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

         If (NumCho(jSym).lt.1) GOTO 1000


C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
         mTvec = 0  ! mem for storing the half-transformed vec

         do l=1,nSym
            k=Muld2h(l,JSYM)
            Mmax = 0
            do jDen=1,nDen
               Mmax = Max(Mmax,nForb(k,jDen)+nIorb(k,jDen))
            end do
            mTvec = mTvec + nBas(l)*Mmax
         end do

         mTvec=Max(mTvec,1)

C ------------------------------------------------------------------
C ------------------------------------------------------------------

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

         Do JRED=JRED1,JRED2

            CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

            If (nVrs.eq.0) GOTO 999  ! no vectors in that (jred,jsym)

            if (nVrs.lt.0) then
               Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
               call abend()
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM//'cho_X_setred non-zero return code.',
     &                        '   rc= ',irc
              call abend()
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            If(JSYM.eq.1)Then

               Call GetMem('rsDtot','Allo','Real',ipDab,nRS)
               Call GetMem('rsFC','Allo','Real',ipFab,nRS)
               Call Fzero(Work(ipDab),nRS)
               Call Fzero(Work(ipFab),nRS)

            EndIf

            Call mma_maxDBLE(LWORK)

            nVec  = Min(LWORK/(nRS+mTvec),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',nRS+mTvec
               WRITE(6,*) 'jsym= ',jsym
               rc = 33
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

            If(JSYM.eq.1)Then
C --- Transform the density to reduced storage
               mode = 'toreds'
               add  = .false.
               nMat=1
               Call move_sto(irc,iLoc,nMat,ipDLT,ipDab,mode,add)
            EndIf

C --- BATCH over the vectors ----------------------------

            nBatch = (nVrs-1)/nVec + 1

            DO iBatch=1,nBatch

               If (iBatch.eq.nBatch) Then
                  JNUM = nVrs - nVec*(nBatch-1)
               else
                  JNUM = nVec
               endif


               JVEC = nVec*(iBatch-1) + iVrs
               IVEC2 = JVEC - 1 + JNUM

               CALL CWTIME(TCR1,TWR1)

               CALL CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,
     &                        NUMV,IREDC,MUSED)

               If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
                  rc=77
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

               If(JSYM.eq.1)Then
C ************ (alpha+beta) COULOMB CONTRIBUTION  ****************
C
C --- Contraction with the density matrix
C ---------------------------------------
C --- V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * DI(rs)
C==========================================================
C
                  CALL CWTIME(TCC1,TWC1)

                  Call mma_allocate(VJ,JNUM,Label='VJ')

                  CALL DGEMV_('T',nRS,JNUM,
     &                 ONE,Lrs,nRS,
     &                 Work(ipDab),1,ZERO,VJ,1)

C --- FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
C===============================================================

                  Fact = dble(min(jVec-iVrs,1))

                  CALL DGEMV_('N',nRS,JNUM,
     &                 FactCI,Lrs,nRS,
     &                 VJ,1,Fact,Work(ipFab),1)


                  CALL CWTIME(TCC2,TWC2)
                  tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                  tcoul(2) = tcoul(2) + (TWC2 - TWC1)

                  Call mma_deallocate(VJ)

               EndIf  ! Coulomb contribution


               iSwap = 2  ! LpJ,b are returned
C *************** EXCHANGE CONTRIBUTIONS  ***********************

               Do jDen=1,nDen

                  nAux(:) =nForb(:,jDen)+nIorb(:,jDen)
                  Call Allocate_Laq(Laq(jDen),nAux,nBas,nVec,JSYM,nSym,
     &                              iSwap)

                  CALL CWTIME(TCR3,TWR3)

                  kMOs = jDen  ! 1--> alpha  2-->beta MOs
                  nMOs = jDen

C --- Set up the skipping flags
C -------------------------------------------------------------
                  Do i=1,nSym

                     k = Muld2h(i,JSYM)
                     iSkip(k) = Min(1,
     &                    nBas(i)*(nForb(k,jDen)+nIorb(k,jDen)))

                  End Do
C -------------------------------------------------------------


C *********************** HALF-TRANSFORMATION  ****************

                  CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                            JSYM,iSwap,IREDC,nMOs,kMOs,POrb,
     &                            Laq,DoRead)

                  CALL CWTIME(TCR4,TWR4)
                  tread(1) = tread(1) + (TCR4 - TCR3)
                  tread(2) = tread(2) + (TWR4 - TWR3)

                  if (irc.ne.0) then
                     rc = irc
                     RETURN
                  endif


                  CALL CWTIME(TCX1,TWX1)

                  Do iSyma=1,nSym

                     iSymk = MulD2h(JSYM,iSyma)

C ---------------------------------------------------------------------
c *** Compute only the LT part of the InActive exchange matrix ********
C
C     FI(ab) = FI(ab) + FactXI * sum_Jk  LkJ,a * LkJ,b
C ---------------------------------------------------------------------
                     NK = nForb(iSymk,jDen) + nIorb(iSymk,jDen)

                     If (iSkip(iSymk).ne.0) Then

                        ISFI = ipFLT(jDen) + ISTLT(iSyma)

                        CALL DGEMM_TRI('T','N',nBas(iSyma),nBas(iSyma),
     &                         NK*JNUM,FactXI,Laq(jDen)%pA(iSymk)%A,
     &                         NK*JNUM,Laq(jDen)%pA(iSymk)%A,NK*JNUM,
     &                         One,Work(ISFI),nBas(iSyma))


                     EndIf

C --------------------------------------------------------------------
                  End Do  !loop over MOs symmetries

                  CALL CWTIME(TCX2,TWX2)
                  texch(1) = texch(1) + (TCX2 - TCX1)
                  texch(2) = texch(2) + (TWX2 - TWX1)


                  Call Deallocate_Laq(Laq(jDen))
               End Do   ! loop over densities

C --------------------------------------------------------------------
C --------------------------------------------------------------------

            END DO  ! end batch loop


            If(JSYM.eq.1)Then
c --- backtransform fock matrix to full storage
               mode = 'tofull'
               add  = .true.
               nMat = nDen
               Call move_sto(irc,iLoc,nMat,ipFLT,ipFab,mode,add)
            EndIf

C --- free memory
            Call mma_deallocate(Lrs)

            If(JSYM.eq.1)Then
              Call GetMem('rsFC','Free','Real',ipFab,nRS)
              Call GetMem('rsDtot','Free','Real',ipDab,nRS)
            EndIf


999         Continue

         END DO   ! loop over red sets

1000     CONTINUE

      END DO   !loop over JSYM

      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1


*
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky SCF timing from '//SECNAM
      Write(6,CFmt)'------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ/TRANSFORM VECTORS           '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB                          '
     &                           //'         ',tcoul(1),tcoul(2)
         Write(6,'(2x,A26,2f10.2)')'EXCHANGE                         '
     &                           //'         ',texch(1),texch(2)
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif


c Print the Fock-matrix
#ifdef _DEBUGPRINT_
      if(Debug) then !to avoid double printing in SCF-debug

      WRITE(6,'(6X,A)')'TEST PRINT FROM '//SECNAM
      WRITE(6,'(6X,A)')
      WRITE(6,'(6X,A)')'***** FOCK MATRIX AO-BASIS ***** '
      Do jDen=1,nDen
        if(nDen.eq.2) Then
          if(jden.eq.1) WRITE(6,'(6X,A)')'******** ALPHA SPIN ******** '
          if(jden.eq.2) WRITE(6,'(6X,A)')'******** BETA SPIN ********* '
        endif
        DO ISYM=1,NSYM
           ISFI=ipFLT(jDen)+ISTLT(ISYM)
           IF( NBAS(ISYM).GT.0 ) THEN
             WRITE(6,'(6X,A)')
             WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
             call TRIPRT('','',Work(ISFI),NBAS(ISYM))
           ENDIF
        END DO
      END DO

      endif

#endif

      rc  = 0


      Return
      END

**************************************************************
**************************************************************



      SUBROUTINE move_sto(irc,iLoc,nDen,ipXLT,ipXab,mode,add)
      use ChoArr, only: iRS2F
      use ChoSwp, only: IndRed
      Implicit Real*8 (a-h,o-z)
      Integer  ISLT(8),cho_isao,nDen
      External cho_isao
      Integer ipXLT(nDen),ipXab
      Logical add
      Character*6 mode

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

************************************************************************
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
************************************************************************


c Offsets to symmetry block in the LT matrix
      ISLT(1)=0
      DO ISYM=2,NSYM
         ISLT(ISYM) = ISLT(ISYM-1)
     &              + NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
      END DO

**************************************************

      jSym = 1 ! only total symmetric density

      xf=0.0d0
      if (add) xf=1.0d0 !accumulate contributions

      If (mode.eq.'toreds') then

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            Do jDen=1,nDen

               kfrom = ipXLT(jDen) + isLT(iSyma) + iab - 1

               Work(ipXab+jRab-1) = xf*Work(ipXab+jRab-1)
     &                            +    Work(kfrom)

            End Do

         End Do  ! jRab loop

      ElseIf (mode.eq.'tofull') then

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            Do jDen=1,nDen

               kto = ipXLT(jDen) + isLT(iSyma) + iab - 1

               Work(kto) = xf*Work(kto)
     &                   +    Work(ipXab+jRab-1)

            End Do


         End Do  ! jRab loop

      Else

         write(6,*)'Wrong input parameter. mode = ',mode
         irc = 66
         Call abend()

      EndIf

      irc = 0

      Return
      End
