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
      SUBROUTINE SYGTOSD(ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,CISYG,CISD,
     &                   detocc,detcoeff,SPTRA)
      use stdalloc, only: mma_allocate, mma_deallocate
#     include "intent.fh"
      IMPLICIT NONE
      INTEGER ICNFTAB(*),ISPNTAB(*),ISSTAB(*),IFSBTAB(*)
      REAL*8 CISYG(*),CISD(*)
      character(len=*), intent(_OUT_) :: detocc(*)
      real(8) :: detcoeff(*), SPTRA(*)

      INTEGER NASPRT
      INTEGER MINOP,MAXOP,NACTEL,NOPEN,NCLSD
      INTEGER NCNF,NSPD,NCPL,ISPART,ISST
C     INTEGER KFSB,IBLK,ISPD,I,IPOS,IORB,ISYM
      INTEGER KFSB,IBLK,ISPD,I,IPOS,IORB
      INTEGER ISPN,NO,IOSTA,IOEND,IMORS,ISBSTR
      INTEGER ICNF,IEL,IEL1,IEL2,IFORM,IFSB,IOCC
      INTEGER IPART,IREST
      INTEGER ISORB,ISPEND,ISPSTA,ISUM,ISYGEND,ISYGSTA,NWRD
      INTEGER IWORD,IWRD,JSST,KCNF,KCNFINF,KGSLIM,KGSORB,KHSHMAP
      INTEGER KMRSSBS,KSPNINF,KSSTARR,KSSTTB,KSPN,LSPTRA
      INTEGER LSYM,MORSBITS,MXBLK,NAPART,NBLK,NHEAD
      INTEGER NHSHMAP,NOCC,NOP,NORB,NSP,NSSTP,NSYM
C     INTEGER IERR,ICPL,KSBSMRS,JMORS,NFSB
      INTEGER IERR,                   NFSB
      INTEGER OCC2MRS
      EXTERNAL OCC2MRS
CC add an occupation array in the usual 0,u,d,2 format
      character(len=1), allocatable :: occ(:)
      integer :: idet

      Real*8, Allocatable:: BLK(:)
      Integer, Allocatable:: OrbArr(:), OccArr(:), STArr(:), DIM(:),
     &                       SSArr(:), SBSET(:)
CC CC

C Unbutton the configuration table:
      NACTEL=ICNFTAB(3)
      NORB  =ICNFTAB(4)
      MINOP =ICNFTAB(5)
      MAXOP =ICNFTAB(6)
      NSYM  =ICNFTAB(7)
      LSYM  =ICNFTAB(8)
      NAPART=ICNFTAB(9)
      IFORM =ICNFTAB(10)
      NHEAD=10
      KGSORB =NHEAD+1
      KGSLIM =KGSORB+(NSYM+1)*(NAPART+1)
      KCNFINF =KGSLIM+2*NAPART
CTEST      write(*,*)' Table at KGSORB:'
CTEST      do ipart=0,napart
CTEST       write(*,'(1x,10i5)')(icnftab(kgsorb+isym+(nsym+1)*ipart)
CTEST     &                                             ,isym=0,nsym)
CTEST      end do
CTEST      write(*,*)' Table at KGSLIM:'
CTEST      write(*,'(1x,10i5)')(icnftab(kgslim+2*(ipart-1)),ipart=1,napart)
CTEST      write(*,'(1x,10i5)')(icnftab(kgslim+1+2*(ipart-1)),ipart=1,napart)
C Unbutton the spin coupling table:
      KSPNINF=9
C Unbutton the Substring Table:
      NASPRT=ISSTAB(5)
      MORSBITS=ISSTAB(6)
      NSSTP =ISSTAB(7)
      KSSTTB=15
CTEST      KSBSMRS=ISSTAB(11)
      KMRSSBS=ISSTAB(12)
C Unbutton the Fock Sector Block table:
      NHEAD=7
      KSSTARR=NHEAD+1
      NFSB   =IFSBTAB(2)
      NHSHMAP=IFSBTAB(6)
      KHSHMAP=IFSBTAB(7)

C MXBLK=Largest individual SYG block of determinants:
      MXBLK=0
      idet=0
      DO NOPEN=MINOP,MAXOP
        NCNF=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP)))
        NCPL=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+2)
        MXBLK=MAX(NCNF*NCPL,MXBLK)
      END DO
C A variety of small temporary arrays used:
      CALL mma_allocate(BLK,MXBLK,Label='BLK')
      CALL mma_allocate(ORBARR,NACTEL,Label='OrbArr')
      CALL mma_allocate(OCCARR,2*NORB,Label='OccArr')
      CALL mma_allocate(STARR,NASPRT,Label='STArr')
      CALL mma_allocate(DIM,NASPRT,Label='Dim')
      CALL mma_allocate(SSARR,NASPRT,Label='SSArr')
      CALL mma_allocate(SBSET,NSSTP,Label='SBSet')
      call mma_allocate(occ,norb,label='occ')
C We will need later the accumulated number of substrings of
C earlier substring types:
      ISUM=0
      DO ISST=1,NSSTP
        SBSET(ISST)=ISUM
        ISUM=ISUM+ISSTAB(KSSTTB+5*(ISST-1))
      END DO
C Loop over nr of open shells.
      ISYGEND=0
      DO NOPEN=MINOP,MAXOP
        NCLSD=(NACTEL-NOPEN)/2
        IF(NCLSD.LT.0) cycle
        IF(2*NCLSD+NOPEN.NE.NACTEL) cycle
        NOCC=NCLSD+NOPEN
        IF(NOCC.GT.NORB) cycle
        NCNF=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP)))
        IF(NCNF.EQ.0) cycle
        KCNF=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+1)
        NWRD=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+2)
        NCPL=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+1)
        NSPD=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+2)
        KSPN=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+4)
        IF(NSPD.EQ.0) cycle
        IF(NCPL.EQ.0) cycle
C ISYGSTA=1st element of each block
        NBLK=NCPL*NCNF
        ISYGSTA=ISYGEND+1
        ISYGEND=ISYGEND+NBLK
C Location of spin coupling coefficients:
        LSPTRA=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+5)
C Matrix multiplication into temporary array:
        CALL  DGEMM_('N','N',NSPD,NCNF,NCPL,1.0D0,
     &               SPTRA(LSPTRA),NSPD,
     &               CISYG(ISYGSTA),NCPL,0.0D0,
     &               BLK,NSPD)

C There is no phase factor in the reorder of orbitals from SYG
C to SD. But there is a fairly lengthy procedure for finding the
C correct position of each CI coefficient.
        IBLK=0
C Loop over configurations
        IWORD = 0 ! dummy initialize
        DO ICNF=1,NCNF
          IF(IFORM.EQ.1) THEN
            DO IEL=1,NOCC
              ORBARR(+IEL)=ICNFTAB(KCNF-1+IEL+NWRD*(ICNF-1))
            END DO
          ELSE IF(IFORM.EQ.2) THEN
            IEL2=0
            IEL1=NCLSD
            DO IORB=1,NORB
              IOCC=ICNFTAB(KCNF-1+IORB+NWRD*(ICNF-1))
              IF(IOCC.EQ.1) THEN
                IEL1=IEL1+1
                ORBARR(+IEL1)=IORB
              ELSE
                IEL2=IEL2+1
                ORBARR(+IEL2)=IORB
              END IF
            END DO
          ELSE IF(IFORM.EQ.3) THEN
            DO IEL=1,NOCC
              IWRD=(3+IEL)/4
              IREST=(3+IEL)-4*IWRD
              IF(IREST.EQ.0) THEN
                IWORD=ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
              END IF
              IORB=MOD(IWORD,256)
              IWORD=IWORD/256
              ORBARR(+IEL)=IORB
            END DO
          ELSE IF(IFORM.EQ.4) THEN
            IEL2=0
            IEL1=NCLSD
            DO IORB=1,NORB
              IWRD=(IORB+14)/15
              IREST=IORB+14-15*IWRD
              IF(IREST.EQ.0) THEN
                IWORD=ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
              END IF
              IOCC=MOD(IWORD,4)
              IWORD=IWORD/4
              IF(IOCC.EQ.1) THEN
                IEL1=IEL1+1
                ORBARR(+IEL1)=IORB
              ELSE
                IEL2=IEL2+1
                ORBARR(+IEL2)=IORB
              END IF
            END DO
          END IF

CTEST      write(*,'(1x,a,10i5)')'Configuration:',
CTEST     &                          (ORBARR(+iel),iel=1,nocc)
CTEST      write(*,*)' Loop over spin determinants.'
C Loop over spin determinants
          DO ISPD=1,NSPD
CTEST      write(*,'(1x,a,20i3)')'Spin determinant:',
CTEST     &                    (ISPNTAB(KSPN-1+I+NOPEN*(ISPD-1)),I=1,nopen)
            IBLK=IBLK+1
C count the determinants
            idet=idet+1
C Construct occupation number array:
            CALL ICOPY(2*NORB,[0],0,OCCARR,1)
            DO IEL=1,NCLSD
              IORB=ORBARR(+IEL)
              OCCARR(+2*IORB-1)=1
              OCCARR(+2*IORB  )=1
            END DO
            DO I=1,NOPEN
C Spin of each electron is coded as 1 for alpha, 0 for beta.
              ISPN=ISPNTAB(KSPN-1+I+NOPEN*(ISPD-1))
              IEL=NCLSD+I
              ISORB=2*ORBARR(+IEL)-ISPN
              OCCARR(+ISORB)=1
            END DO
C Identify substrings:
C Loop over active partitions. Subdivide as needed into subpartitions.
CTEST      write(*,*)' Identify substrings.'
CTEST      write(*,'(1x,a,10i5)')'Occupation array:',
CTEST     &                          (OCCARR(+isorb),isorb=1,2*norb)
CC construct occupation array in 0,u,d,2 format
            do IORB=1,2*norb-1,2
              if ((OCCARR(IORB) == 1)
     &            .and. (OCCARR(1+IORB) == 1)) then
                occ((IORB+1)/2)='2'
              else if ((OCCARR(IORB) == 1)
     &                 .and. (OCCARR(1+IORB) == 0)) then
                occ((IORB+1)/2)='u'
              else if ((OCCARR(IORB) == 0)
     &                 .and. (OCCARR(1+IORB) == 1)) then
                occ((IORB+1)/2)='d'
              else if((OCCARR(IORB) == 0)
     &                .and. (OCCARR(1+IORB) == 0)) then
                occ((IORB+1)/2)='0'
              end if
            end do
CC
            IOEND=0
            ISPEND=0
CTEST      write(*,'(1x,a,10i5)')'NAPART:',NAPART
            DO IPART=1,NAPART
CTEST      write(*,'(1x,a,10i5)')' In loop, IPART:',IPART
             NOP=2*ICNFTAB(KGSORB+(NSYM+1)*IPART)
CTEST      write(*,'(1x,a,10i5)')'            NOP:',NOP
CTEST             IF(NOP.EQ.0) write(*,*)' (Skip it.)'
             IF(NOP.EQ.0) cycle
             NSP=(NOP+MORSBITS-1)/MORSBITS
CTEST      write(*,'(1x,a,10i5)')'            NSP:',NSP
             ISPSTA=ISPEND+1
             ISPEND=ISPEND+NSP
             DO ISPART=ISPSTA,ISPEND
CTEST      write(*,'(1x,a,10i5)')'Loop lims ISPSTA,ISPEND:',ISPSTA,ISPEND
              NO=MIN(NOP,MORSBITS)
              NOP=NOP-NO
              IOSTA=IOEND+1
              IOEND=IOEND+NO
CTEST      write(*,'(1x,a,10i5)')'IOSTA,IOEND:',IOSTA,IOEND
CTEST      write(*,'(1x,a,10i5)')'Occ array:',
CTEST     &                     (OCCARR(+ISORB),ISORB=IOSTA,IOEND)
              IMORS=OCC2MRS(NO,OCCARR(+IOSTA))
CTEST      write(*,'(1x,a,10i5)')'IMORS=',IMORS
C Position in Morsel-to-Substring table:
              IPOS=KMRSSBS+2*(IMORS+(2**MORSBITS)*(ISPART-1))
C Substring ID number
              ISBSTR=ISSTAB(IPOS)
C Test:
CTEST              JMORS=ISSTAB(KSBSMRS+2*(ISBSTR-1))
CTEST              IF(IMORS.NE.JMORS)THEN
CTEST                WRITE(*,*)' Mistranslated morsel!!'
CTEST      write(*,'(1x,a,4i12)')'IMORS->ISBSTR:',IMORS,ISBSTR
CTEST      write(*,'(1x,a,4i12)')'but ISBSTR->IMORS:',ISBSTR,JMORS
CTEST      write(*,'(1x,a,4i12)')'KMRSSBS:',KMRSSBS
CTEST      write(*,'(1x,a,4i12)')'KSBSMRS:',KSBSMRS
CTEST                CALL ABEND()
CTEST              END IF
CTEST      write(*,'(1x,a,10i5)')'ISBSTR:',ISBSTR
C Substring type ISST, nr of such substrings is NDIM
              ISST=ISSTAB(IPOS+1)
CTEST      write(*,'(1x,a,10i5)')'ISST  :',ISST
              STARR(+ISPART)=ISST
              DIM(ISPART)=ISSTAB(KSSTTB+5*(ISST-1))
              SSARR(ISPART)=ISBSTR-SBSET(ISST)
             END DO
            END DO
CTEST      write(*,*)' Finally, substring types and substrings:'
CTEST      write(*,'(1x,a,10i5)')'Substr types:',
CTEST     &                   (STARR(+ISPART),ISPART=1,NASPRT)
CTEST      write(*,'(1x,a,10i5)')'Substrings  :',
CTEST     &                    (SSARR(ISPART),ISPART=1,NASPRT)
CTEST      write(*,'(1x,a,10i5)')'Dimensions  :',
CTEST     &                    (DIM(ISPART),ISPART=1,NASPRT)
C Position within FS block:
            IPOS=(SSARR(NASPRT)-1)
            DO ISPART=NASPRT-1,1,-1
              IPOS=DIM(ISPART)*IPOS+(SSARR(ISPART)-1)
            END DO
            IPOS=IPOS+1
C Identify Fock Sector Block:
CTEST      write(*,*)' Arguments in HSHGET call:'
CTEST      write(*,'(1x,a,10i5)')'Key:',
CTEST     &                   (STARR(+ISPART),ISPART=1,NASPRT)
CTEST      write(*,'(1x,a,10i5)')'Size of key:',NASPRT
CTEST      write(*,'(1x,a,10i5)')'Size of items stored:',NASPRT+2
CTEST      write(*,'(1x,a,10i5)')'Items stored at KSSTARR=',KSSTARR
CTEST      write(*,'(1x,a,10i5)')'      Map size  NHSHMAP=',NHSHMAP
CTEST      write(*,'(1x,a,10i5)')'  Map stored at KHSHMAP=',KHSHMAP
            CALL HSHGET(STARR,NASPRT,NASPRT+2,IFSBTAB(KSSTARR),
     &                            NHSHMAP,IFSBTAB(KHSHMAP),IFSB)
CTEST      write(*,'(1x,a,10i5)')' Map returns index IFSB=',IFSB
CTEST      write(*,'(1x,a,10i5)')' Item stored there is  =',
CTEST     &                 (IFSBTAB(KSSTARR-1+ISPART+(NASPRT+2)*(IFSB-1)),
CTEST     &                                              ISPART=1,NASPRT+2)
C Position of this FS block in SD wave function:
            KFSB=IFSBTAB(KSSTARR+(NASPRT+2)*IFSB-1)
C Temporary check, may be removed later. See that we have picked up
C the correct FS block.
            IERR=0
            DO ISPART=1,NASPRT
              JSST=IFSBTAB(KSSTARR-1+ISPART+(NASPRT+2)*(IFSB-1))
              ISST=STARR(+ISPART)
              IF(ISST.NE.JSST) IERR=1
            END DO
            IF(IERR.NE.0) THEN
              WRITE(6,*) ' SYGTOSD Error:'//
     &                     ' Hash map returned the wrong FS block!'
              WRITE(6,'(1x,a,8I8)')'NOPEN,ICNF,ISPN:',NOPEN,ICNF,ISPN
              WRITE(6,'(1x,a,20I3)')'Configuration:',
     &                               (ORBARR(+IEL),IEL=1,NACTEL)
              WRITE(6,'(1x,a,20I3)')'Determinant:',
     &                            (OCCARR(+ISORB),ISORB=1,2*NORB)
              WRITE(6,'(1x,a,10I5)')'Substring type combination:',
     &                         (STARR(+ISPART),ISPART=1,NASPRT)
              WRITE(6,'(1x,a,10I5)')'Substring combination:',
     &                         (SSARR(ISPART),ISPART=1,NASPRT)
              WRITE(6,'(1x,a,8I8)')'Hash table says IFSB=',IFSB
              IF(IFSB.GT.0 .AND. IFSB.LE.NFSB) THEN
              WRITE(6,'(1x,a,8I8)')'but that FS block would contain',
     & (IFSBTAB(KSSTARR-1+ISPART+(NASPRT+2)*(IFSB-1)),ISPART=1,NASPRT)
              ELSE
              WRITE(6,*)'but there is no such FS block!'
              WRITE(6,*)' The FS block table follows:'
              CALL PRFSBTAB(IFSBTAB)
              END IF
              CALL ABEND()
            END IF
C Finally:
            CISD(KFSB-1+IPOS)=BLK(IBLK)
            detcoeff(idet)=CISD(KFSB-1+IPOS)
            write(detocc(idet),*) occ
C End of spin-determinant loop
          END DO
C End of loop over configurations
        END DO
C End of loop over nr of open shells
      END DO
      CALL mma_deallocate(BLK)
      CALL mma_deallocate(ORBARR)
      CALL mma_deallocate(OCCARR)
      CALL mma_deallocate(STARR)
      CALL mma_deallocate(DIM)
      CALL mma_deallocate(SSARR)
      CALL mma_deallocate(SBSET)
      call mma_deallocate(occ)
      RETURN
      END
