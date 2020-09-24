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
* Copyright (C) 1987, Bjorn O. Roos                                    *
*               1992, Per Ake Malmqvist                                *
*               1998, Jun-ya Hasegawa                                  *
************************************************************************
*--------------------------------------------*
* 1987  B. O. ROOS                           *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE TRACTL(iPart)
C  SECOND ORDER TWO-ELECTRON TRANFORMATION PROGRAM. CONTROL SECTION
C
C  THIS SUBROUTINE SETS UP THE MEMORY ALLOCATIONS FOR TRA2 AND LOOPS
C  OVER THE SYMMETRY BLOCKS. TRA2 IS CALLED ONCE FOR EACH SYMMETRY
C  BLOCK OF INTEGRALS. SYMMETRY BLOCKED AO INTEGRALS MUST HAVE BEEN
C  GENERATED BY INTSORT ON UNIT LUINTA=40.
C
C WRITTEN IN GARCHING IN SEPTEMBER 1987
C AUTHOR: BYOERN ROOS
C         DEPARTMENT OF THEORETICAL CHEMISTRY
C         CHEMICAL CENTRE
C         P.O.B. 124
C         S-221 00 LUND SWEDEN     TEL: 46-10 82 51
C     ********** IBM-3090 RELEASE 87 09 14 **********
C 92-12-04 P-AA M: Changed for use with CASPT2 MOLCAS-3 version.
C Also transforms 1-el integrals. -> Actually that part is commented out
C
c 98-09-02 J.Hasegawa Modified for non-squared integrals.
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL IFTEST
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "intgrl.fh"
#include "SysDef.fh"
#include "trafo.fh"
      DIMENSION nBasXX(8),Keep(8)
      Logical iSquar
      Logical DoCholesky

      IFTEST=.FALSE.
#ifdef _DEBUG_
      IfTest=.True.
#endif

C Copy data to common ERI.
      NSYMZ=NSYM
      DO I=1,NSYM
        NORBZ(I)=NORB(I)
        NOSHZ(I)=NOSH(I)
        LUINTMZ=LUINTM
      END DO

C Open temporary files for half-transformed integrals.
C They are closed at end of TRACTL.

C The MO coefficients were allocated and read in STINI. They are
C available at WORK(LCMO).
C
C     RETRIEVE BASE DATA FROM UNIT LUINTA
C
**JHsta
      Call GetOrd(IRC,iSquar,nSymXX,nBasXX,Keep)
      IF ( OUTFMT.EQ.'LONG    ' ) THEN
        If(iSquar)      write(6,*)'TRACTL OrdInt status: squared'
        If(.not.iSquar) write(6,*)'TRACTL OrdInt status: non-squared'
      ENDIF
**JHend
      IF(IRC.NE.0) THEN
        WRITE(6,*)' TRACTL, called to transform the two-electron'
        WRITE(6,*)' integrals, got non-zero return code from'
        WRITE(6,*)' subroutine GETORD. The return code is IRC=',IRC
        WRITE(6,*)' Do you have a valid ORDINT file? If you do,'
        WRITE(6,*)' please inform the MOLCAS group -- this may be'
        WRITE(6,*)' a bug. Anyway, the calculations must stop, sorry.'
        CALL QUIT(_RC_IO_ERROR_READ_)
      END IF
** PAM2007: For unknown reasons, one extra word is needed.
*      lBuf = 1+MAX(255*255,NBMX**2)
* but note that tractl, being a utility, can be called from other
* programs and should not take the value NBMX from caspt2.fh...
* hence correction by AJS below.
c
c     Correction by AJS, Jan. 12, 2009. Defines the value
c     of NBMX
c     ---------------------------------------------------
c
      NBMX=1
      DO I=1,NSYM
c        WRITE(6,'(a,i2,a,i5)') 'ISYM=',I,'   nBas(ISYM)=',nBasXX(I)
        NBMX=MAX(NBMX,nBasXX(i))
      ENDDO
c     ---------------------------------------------------
*      lBuf = MAX(255*255,NBMX**2)
      lBuf = 1+NBMX**2
c      write(6,'(2(A,I10))') 'NBMX=',NBMX,'   lBuf=',lBuf
*
* COMPARE CONTENT OF 1EL and 2EL INTEGRAL FILE
      IERR=0
      IF ( NSYMXX.NE.NSYM ) THEN
        IERR=1
      ELSE
        DO ISYM=1,NSYM
          IF (NBAS(ISYM).NE.NBASXX(ISYM)) IERR=1
        END DO
      END IF
      IF(IERR.NE.0) THEN
        WRITE(6,*)'     *** ERROR IN SUBROUTINE TRACTL ***'
        WRITE(6,*)'          INCOMPATIBLE BASIS DATA'
        WRITE(6,*)
        WRITE(6,*)' JOBIPH NR OF SYMM:', NSYM
        WRITE(6,*)' JOBIPH NR OF BASIS FUNCTIONS/SYMM:'
        WRITE(6,'(1x,8I5)')(NBAS(I),I=1,NSYM)
        WRITE(6,*)
        WRITE(6,*)' ORDINT NR OF SYMM:', NSYMXX
        WRITE(6,*)' ORDINT NR OF BASIS FUNCTIONS/SYMM:'
        WRITE(6,'(1x,8I5)')(NBASXX(I),I=1,NSYMXX)
        CALL ERRTRA
        CALL SYSHALT('TRACTL')
      END IF
C
C     SET ADDRESS FIELD FOR OUTPUT INTEGRAL FILE
C
      LIADUT=3*36*36
      DO I=1,36*36
       IAD2M(1,I)=0
       IAD2M(2,I)=0
       IAD2M(3,I)=0
      END DO
      IAD13=0
      CALL iDAFILE(LUINTM,1,IAD2M,LIADUT,IAD13)
C
C     LOOP OVER QUADRUPLES OF SYMMETRIES (NSP,NSQ,NSR,NSS)
C     NOTE THAT THE INTEGRALS ON LUINTA HAVE TO BE SORTED IN THE SAME
C     ORDER AS THE LOOP STRUCTURE BELOW (USE PROGRAM INTSORT)
C
c Allocate largest possible array as work space:
      If ( IFTEST ) then
        Write(6,*)
     &  ' Symmetry  Basis functions   total orbitals    active orbitals'
        Write(6,*)
     &  ' -------------------------------------------------------------'
      End If
      CALL GETMEM('LW1','MAX','REAL',LW1,MEMX)
      MEMX=MAX(MEMX-1*MEMX/6,0)

C --- Reserve space for Integral generation from Cholesky vectors
      Call DecideOnCholesky(DoCholesky)
      If (DoCholesky) then
         MEMX = MAX(MEMX-MEMX/10,0)
         write(6,*)'Memx= ',MEMX
      End If
C -----------------------------------------------------------------
      CALL GETMEM('LW1','ALLO','REAL',LW1,MEMX)
      NCHAIN=0
      LMOP1=1
      ITP=0
      DO NSP=1,NSYM
       IF(NSP.NE.1) ITP=ITP+NASH(NSP-1)
       NBP=NBAS(NSP)
       IF(NSP.NE.1) LMOP1=LMOP1+NBAS(NSP-1)**2
       LMOP=LMOP1+NBP*NFRO(NSP)
       LMOP2=LMOP
       NOP=NORB(NSP)
       NOCP=NOSH(NSP)
       KEEPP=KEEP(NSP)
       ISP=NSP
       LMOQ1=1
       ITQ=0
       DO NSQ=1,NSP
        IF(NSQ.NE.1) ITQ=ITQ+NASH(NSQ-1)
        NBQ=NBAS(NSQ)
        IF(NSQ.NE.1) LMOQ1=LMOQ1+NBAS(NSQ-1)**2
        LMOQ=LMOQ1+NBQ*NFRO(NSQ)
        LMOQ2=LMOQ
        KEEPQ=KEEP(NSQ)
        NOQ=NORB(NSQ)
        NOCQ=NOSH(NSQ)
        ISQ=NSQ
        NSPQ=MUL(NSP,NSQ)
        LMOR1=1
        ITR=0
**JHsta
        NSymR=NSP
        If(iSquar)NSymR=NSYM
        DO NSR=1,NSymR
**JHend
         IF(NSR.NE.1) ITR=ITR+NASH(NSR-1)
         NBR=NBAS(NSR)
         IF(NSR.NE.1) LMOR1=LMOR1+NBAS(NSR-1)**2
         LMOR=LMOR1+NBR*NFRO(NSR)
         LMOR2=LMOR
         KEEPR=KEEP(NSR)
         NOR=NORB(NSR)
         NOCR=NOSH(NSR)
         NSPQR=MUL(NSPQ,NSR)
         ISR=NSR
         LMOS1=1
         ITS=0
**JHsta
         NSymS=NSR
         If(NSP.eq.NSR)NSymS=NSQ
         If(iSquar)NSymS=NSR
         DO NSS=1,NSymS
**JHend
          IF(NSS.NE.1) ITS=ITS+NASH(NSS-1)
          NBS=NBAS(NSS)
          IF(NSS.NE.1) LMOS1=LMOS1+NBAS(NSS-1)**2
          LMOS=LMOS1+NBS*NFRO(NSS)
          LMOS2=LMOS
          IF(NSPQR.NE.NSS) GO TO 101
          NOS=NORB(NSS)
          NOCS=NOSH(NSS)
          KEEPS=KEEP(NSS)
          ISS=NSS
C
          KEEPT=KEEPP+KEEPQ+KEEPR+KEEPS
          NOCCT=NOCP*NOCQ*NOCR*NOCS
          IF(NOCCT.NE.0.AND.KEEPT.NE.0) GO TO 901
          IF(KEEPT.EQ.0) NCHAIN=NCHAIN+1
          IF(NOP*NOQ*NOR*NOS.EQ.0) GO TO 101
C
C         CALLING SEQUENCE FOR SECOND ORDER TRANSFORMATION TRA2
C         FIRST ALLOCATE AND CHECK MEMORY
C
          NBPQ=NBP*NBQ
          IF(ISP.EQ.ISQ) NBPQ=(NBP**2+NBP)/2
          NBRS=NBR*NBS
          IF(ISR.EQ.ISS) NBRS=(NBR**2+NBR)/2
          NOTU=NOCR*NOCS
          IF(ISR.EQ.ISS) NOTU=(NOCR**2+NOCR)/2

          If ( IFTEST ) then
             Write(6,'(1X,4I2,2X,4I4,2X,4I4,2X,4I4)')
     &            NSP,NSQ,NSR,NSS,
     &            NBP,NBQ,NBR,NBS,
     &            NOP,NOQ,NOR,NOS,
     &            NOCP,NOCQ,NOCR,NOCS
          End If

CPAM      MEMX=MXMEM
CPAM      LW1=1
cJHsta
          Mxx1=max(lBuf,NOP*NBQ,NBP*NOCQ,NOCP*NBQ)
          Mxx2=max(NBR*NBS,NBP*NBQ,NOP*NOR,NOP*NOS,NOQ*NOR,NOQ*NOS)
          Mxx3=max(NOCR*NBS,NBR*NOCS)
          LW2=LW1+Mxx1
          LW3=LW2+Mxx2
          LW4=LW3+Mxx3
          LRUPQ=NBP*NBQ*NBR*NOCS
          LURPQ=NBP*NBQ*NOCR*NBS
          MEMLFT=MEMX-LW4+LW1
* I.E. MEMLFT=MEMX-MXX1-MXX2-MXX3, possibly negative...
          LPQTU=NBP*NBQ*NOCR*NOCS
          LATRU=NOP*NOCQ*NBR*NOCS
          LTARU=NOCP*NOQ*NBR*NOCS
          LATUS=NOP*NOCQ*NOCR*NBS
          LTAUS=NOCP*NOQ*NOCR*NBS
          LTUPQ=max(LPQTU,LATRU,LTARU)
          MEMT=LRUPQ+LURPQ+LTUPQ
          L2=max(LATUS,LTAUS)

          LRUPQM=NBR*NOCS
          If(LRUPQM.ne.0)LRUPQM=max(NBR*NOCS,NBPQ)
          LURPQM=NBS*NOCR
          If(LURPQM.ne.0)LURPQM=max(NBS*NOCR,NBPQ)
          LTUPQM=MAX(NOTU,NOCQ*NOCS*NOP,NOCP*NOCS*NOQ)
          If(LTUPQM.ne.0)  LTUPQM=max(LTUPQM,NBPQ,NBR*NOP,NBR*NOQ)
          L2M=MAX(NOCQ*NOCR*NOP,NOCP*NOCR*NOQ)
          If(L2M.ne.0)L2M=max(NOP*NBS,NOQ*NBS)

          IF(MEMT.GT.MEMLFT.OR.L2.GT.MEMLFT-LURPQ) THEN
C           LRUPQ=INT((1.0D0*MEMLFT*LRUPQ+MEMT-1)/MEMT)
C           LURPQ=INT((1.0D0*MEMLFT*LURPQ+MEMT-1)/MEMT)
            NX=MEMLFT/(LRUPQM+LURPQM+LTUPQM)
            iiPart=1
            if(iPart.gt.0) then
            iiPart=iPart
            endif
            LRUPQ=NX*LRUPQM*iiPart
            LURPQ=NX*LURPQM*iiPart
          ENDIF
          LTUPQ=MAX(0,MEMLFT-LRUPQ-LURPQ)
c          print *,'LRUPQ=',LRUPQ
c          print *,'LURPQ=',LURPQ
c          print *,'LTUPQ=',LTUPQ

          IF(LRUPQ.LT.LRUPQM) GO TO 902
          IF(LURPQ.LT.LURPQM) GO TO 902
          IF(LTUPQ.LT.LTUPQM) GO TO 902
          IF(LRUPQ+LTUPQ.LT.L2M) GO TO 902
          LW5=LW4+LURPQ
          LW6=LW5+LRUPQ

          If(.not.iSquar)then
* Keep addresses LW2, LW3, LW4, and save LTUPQ (in common /TRAFO/),
* for use in the calls to tr2Sq or to tr2NsA.
           LTUPQX=LTUPQ

* Recompute memory requirements, now for the tr2NsB call.
           Mxx1=max(lBuf,NBP*NOCQ,NBR*NOS,NOR*NBS)
           Mxx2=max(NBP*NBQ,NBR*NBS)
C LPQRS, MEMT integers, changed to Re*8. Defined and used in the
C following section only. Named XLPQRS, XMEMT, +small local changes.
           XLPQRS=DBLE(NBRS)
           XLPQRS=XLPQRS*DBLE(NBP*NBQ)
           LTURS=NOCP*NOCQ*NBR*NBS
* Mxx1 words needed...
           LW2B=LW1+Mxx1
* Another Mxx2 words needed...
           LW3B=LW2B+Mxx2
* The next line can also be written ''MEMLFT = MEMX-MXX1-MXX2''
* This could possibly be small or negative.
           MEMLFT=MAX(0,MEMX-LW3+LW1)
           XMEMT=XLPQRS+DBLE(LTURS)
* XMEMT is NBRS*NBP*NBQ + NOCP*NOCQ*NBR*NBS
           IF(XMEMT.GT.DBLE(MEMLFT)) THEN
            LPQRS=INT((DBLE(MEMLFT)*XLPQRS)/XMEMT+0.5d0)
           ELSE
            LPQRS=INT(XLPQRS)
           ENDIF
C XLPQRS, XMEMT not used after this.
           LRSmx=LPQRS/NBPQ
           If(LRSmx.gt.NBRS)LRSmx=NBRS
           Nread=NBRS/LRSmx
           Nrest=mod(NBRS,LRSmx)
           If(Nrest.ne.0)Nread=Nread+1
           LRS=NBRS/Nread
           Nrest=mod(NBRS,Nread)
           If(Nrest.ne.0)LRS=LRS+1
           LPQRS=LRS*NBPQ
           MaxRS=LRS
           LTURS=MEMLFT-LPQRS
           IF(LPQRS.LT.NBPQ) GO TO 903
           LTURSM=NOCP*NOCQ
           If(LTURSM.ne.0) LTURSM=max(LTURSM,NBRS)
           IF(LTURS.LT.LTURSM) GO TO 903
           LW4B=LW3B+LPQRS
           LTUPQ=LTURS
          Endif
          If (iSquar)then
*            TR2Sq(CMO,X1,X2,X3,URPQ,RUPQ,TUPQ,lBuf)
             Call tr2Sq(Work(LCMO),
     &                  Work(LW1),
     &                  Work(LW2),
     &                  Work(LW3),
     &                  Work(LW4),
     &                  Work(LW5),
     &                  Work(LW6),lBuf)
          Else
*            tr2NsA(CMO,X1,X2,X3,pqUs,pqrU,pqTU,lBuf)
c          LW2=LW1+Mxx1
c          LW3=LW2+Mxx2
c          LW4=LW3+Mxx3
c          LW5=LW4+LURPQ
c          LW6=LW5+LRUPQ

           If (IFTEST) then
            write(6,*) 'Calling tr2Nsa'
            write(6,*) 'MEMX=',MEMX
            write(6,*) 'lLW1=',LW2-LW1
            write(6,*) 'lLW2=',LW3-LW2
            write(6,*) 'lLW3=',LW4-LW3
            write(6,*) 'lLW4=',LW5-LW4
            write(6,*) 'lLW5=',LW6-LW5
            write(6,*) 'lLW6=',MEMX-(LW6-LW1)
            write(6,*)
           End If

             LTUPQ=LTUPQX
             Call tr2NsA1(Work(LCMO),
     &                   Work(LW1),LW2-LW1,
     &                   Work(LW2),LW3-LW2,
     &                   Work(LW3),LW4-LW3,
     &                   Work(LW4),LW5-LW4,
     &                   Work(LW5),LW6-LW5,
     &                   Work(LW6),MEMX-(LW6-LW1),lBuf)
             Call tr2NsA2(Work(LCMO),
     &                   Work(LW1),LW2-LW1,
     &                   Work(LW2),LW3-LW2,
     &                   Work(LW5),LW6-LW5,
     &                   Work(LW6),MEMX-(LW6-LW1))
             Call tr2NsA3(Work(LCMO),
     &                   Work(LW1),LW2-LW1,
     &                   Work(LW2),LW3-LW2,
     &                   Work(LW4),LW5-LW4,
     &                   Work(LW5),MEMX-(LW5-LW1))
*            tr2NsB(CMO,X1,X2,pqrs,TUrs,lBuf,MAXRS)
             LTUPQ=LTURS
             Call tr2NsB(Work(LCMO),
     &                   Work(LW1),
     &                   Work(LW2B),
     &                   Work(LW3B),
     &                   Work(LW4B),lBuf,MaxRS)
          End If
  101     CONTINUE
         END DO
        END DO
       END DO
      END DO
      If ( IFTEST ) then
       Write(6,*)
     & ' --------------------------------------------------------------'
      End If
      CALL GETMEM('LW1','FREE','REAL',LW1,MEMX)
C
C     FINALLY WRITE OUT THE DAFILE ADDRESS LIST ON UNIT 13
C
      IAD13=0
      CALL iDAFILE(LUINTM,1,IAD2M,LIADUT,IAD13)
C
CPAM01 Also transform 1-electron integrals, and put CMOs on LUONEM.
CPAM01      CALL TRAONE(WORK(LCMO),KEEP)

      RETURN
C
C     HERE IF INTERPHASE FROM SORT IN ERROR
C
  901 CONTINUE
      WRITE(6,'(/5X,A,8I6)')
     & 'ERROR IN KEEP PARAMETER FROM INTSORT FILE:  ',(KEEP(I),I=1,NSYM)
      WRITE(6,'(/5X,A,8I6)')
     & 'NOT CONSISTENT WITH OCCUPIED ORBITAL SPACE: ',(NOSH(I),I=1,NSYM)
      WRITE(6,'(/5X,A)') 'PROGRAM STOP IN SUBROUTINE TRACTL'
      GOTO 999
C
C     HERE IF NOT ENOUGH CORE SPACE
C
  902 CONTINUE
      WRITE(6,'(/1X,A)')'NOT ENOUGH CORE SPACE FOR SORTING IN TRA2'
      WRITE(6,'(/1X,A,I12)')'TOTAL SORTING SPACE IS',MEMLFT
      WRITE(6,'(/1X,A,I12,A,I12)')'STEP1: AVAILABLE IS',LRUPQ,
     &                                    '  NEEDED IS',LRUPQM
      WRITE(6,'(/1X,A,I12,A,I12)')'STEP2:    ''''         ',
     &                               LTUPQ,'  NEEDED IS',LTUPQM
      WRITE(6,'(/1X,A,I12,A,I12)')'STEP3:    ''''         ',
     &                            LRUPQ+LTUPQ,'  NEEDED IS',L2M
      GOTO 999

  903 CONTINUE
      WRITE(6,'(/1X,A)')'NOT ENOUGH CORE SPACE FOR SORTING IN TRATWO2'
      WRITE(6,'(/1X,A,I12)')'TOTAL SORTING SPACE IS',MEMLFT
      WRITE(6,'(/1X,A,I12,A,I12)')'STEP1: AVAILABLE IS',LPQRS,
     &                                    '  NEEDED IS',NBPQ
      WRITE(6,'(/1X,A,I12,A,I12)')'STEP1:     ''''        ',LTURS,
     &                                    '   ''''        ',LTURSM
      GOTO 999

 999  CONTINUE
      CALL Abend
      END
