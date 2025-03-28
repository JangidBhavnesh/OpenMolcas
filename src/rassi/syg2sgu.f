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
      SUBROUTINE SYG2SGU(IMODE,SGS,CIS,LSYM,ICNFTAB,ISPNTAB,CIOLD,CINEW)
      use rassi_aux, only: ipglob
      use gugx, only: SGStruct, CIStruct, mxlev
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer IMODE
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Integer LSYM
      INTEGER ICNFTAB(*),ISPNTAB(*)
      Real*8 CIOLD(*),CINEW(*)

      Integer, PARAMETER:: NBUFFER=600,MXCPI=15
      Integer KWALK(NBUFFER)
      Real*8 PHASE(NBUFFER)
      Integer ICNUM(NBUFFER)
      Integer ICASE(400)
      Integer :: IFUP2CS(0:1)=[2,1]
      Integer, Allocatable:: MWS2W(:), OrbArr(:)

C SGS       : Data that define a Split Graph
C qCIS : Data that define a CI array structure
C IMODE=0 transforms a Symmetric Group CI array to SGUGA
C IMODE=1 transforms a Split GUGA CI array to Symm Group
C ...Configuration and Spin Coupling tables, fill this in later.
C CIOLD and CINEW are obvious.

CTEST      write(*,*)' SYG2SGU, LSYM=',LSYM
C Dereference CIS and SGS       for some data:
      NCONF =CIS%NCSF(LSYM)
      NWALK =CIS%nWalk
      CALL mma_allocate(MWS2W,NWALK,Label='MWS2W')
      NSYM  =ICNFTAB(7)
      CALL MSTOW(SGS,CIS,MWS2W,NSYM)
CTEST      write(*,*)' NCONF=',NCONF
C MWS2W is a table which gives the upper or lower walk
C index as function of the MAW sum.

C Inspect the top row of the DRT to find NACTEL and spin:
      NLEV  =SGS%DRT(1,1)
      IF (NLEV.GT.MXLEV) THEN
        WRITE(6,*) ' SYG2SGU: error: number of levels exceeds MXLEV'
        WRITE(6,'(1X,2(A,I4))') ' NLEV = ',NLEV,' MXLEV = ',MXLEV
        CALL AbEnd()
      END IF
      NACTEL=SGS%DRT(1,2)

C Now a good bound on MINOP, the minimum number of open
C shells, would be MLTPLC-1. This is the best bound, and it
C does not depend on any assumed Ms.

CTEST      write(*,*)' Test prints in SYG2SGU.'
C A buffer of packed walks is used:
      NWLKLST=0
      IWLKPOS=1
C Nr of integers used to store each total walk:
      MIPWLK=1+(NLEV-1)/MXCPI
C Max nr of Split-Graph CSF''s, stored as walks in the
C buffer KWALK
      MXWLK=NBUFFER/MIPWLK
C Pick up various data from conf and spin tables
C Unbutton Configuration table:

      MINOP =ICNFTAB(5)
      MAXOP =ICNFTAB(6)
      NSYM  =ICNFTAB(7)
      LSYM  =ICNFTAB(8)
      NAPART=ICNFTAB(9)
      IFORM =ICNFTAB(10)
*PAM07 Statement to fool intel 10.1 compiler do do the right thing:
      IF(MINOP.GT.MAXOP) WRITE(6,*) MINOP,MAXOP

CTEST      write(*,*)' SYG2SGU Test print, Config table:'
CTEST      write(*,'(1x,a,8I8)')'MINOP,MAXOP:',MINOP,MAXOP
CTEST      write(*,'(1x,a,8I8)')'  NSYM:',NSYM
CTEST      write(*,'(1x,a,8I8)')'  LSYM:',LSYM
CTEST      write(*,'(1x,a,8I8)')'NAPART:',NAPART
CTEST      write(*,'(1x,a,8I8)')' IFORM:',IFORM
      NHEAD=10
      KGSORB=NHEAD+1
      KGSLIM =KGSORB+(NSYM+1)*(NAPART+1)
      KCNFINF=KGSLIM+2*NAPART
CTEST      write(*,'(1x,a,8I8)')'KCNFINF:',KCNFINF
CTEST      write(*,*)'   NOPEN    ISYM    NCNF    KCNF    IWRD'
CTEST      do nopen=minop,maxop
CTEST       do isym=1,nsym
CTEST        NCNF=ICNFTAB(KCNFINF+3*(ISYM-1+NSYM*(NOPEN-MINOP)))
CTEST        KCNF=ICNFTAB(KCNFINF+3*(ISYM-1+NSYM*(NOPEN-MINOP))+1)
CTEST        IWRD=ICNFTAB(KCNFINF+3*(ISYM-1+NSYM*(NOPEN-MINOP))+2)
CTEST      write(*,'(1x,8I8)') NOPEN,ISYM,NCNF,KCNF,IWRD
CTEST       end do
CTEST      end do
C Unbutton Spin table:
      KSPNINF=9
C Test prints:
CTEST      write(*,*)' SYG2SGU Test print, spin table:'
CTEST      write(*,'(1x,a,8I8)')' MLTPL:',ISPNTAB(3)
CTEST      write(*,'(1x,a,8I8)')' MS2  :',ISPNTAB(4)
CTEST      write(*,'(1x,a,8I8)')' MINOP:',ISPNTAB(5)
CTEST      write(*,'(1x,a,8I8)')' MAXOP:',ISPNTAB(6)
CTEST      write(*,'(1x,a,8I8)')' IOPEN:',
CTEST     &     (ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+0),NOPEN=MINOP,MAXOP)
CTEST      write(*,'(1x,a,8I8)')' NCL  :',
CTEST     &     (ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+1),NOPEN=MINOP,MAXOP)
CTEST      write(*,'(1x,a,8I8)')' NSPD :',
CTEST     &     (ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+2),NOPEN=MINOP,MAXOP)
CTEST      write(*,'(1x,a,8I8)')'KSPCPL:',
CTEST     &     (ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+3),NOPEN=MINOP,MAXOP)
CTEST      write(*,'(1x,a,8I8)')'KSPDET:',
CTEST     &     (ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+4),NOPEN=MINOP,MAXOP)
CTEST      write(*,'(1x,a,8I8)')'LSPTRA:',
CTEST     &     (ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+5),NOPEN=MINOP,MAXOP)
C End of test prints:
C
      CALL mma_allocate(ORBARR,NACTEL,Label='OrbArr')

C Loop over nr of open shells
C NCSYMG=Nr of Symmetric-Group CSF''s treated so far.
CTEST      write(*,*)' SYG2SGU Loop over NOPEN.'
      NCSYMG=0
      DO NOPEN=MINOP,MAXOP
        NCLSD=(NACTEL-NOPEN)/2
        IF(NCLSD.LT.0) cycle
        IF(2*NCLSD+NOPEN.NE.NACTEL) cycle
        NOCC=NCLSD+NOPEN
        IF(NOCC.GT.NLEV) cycle
CTEST      write(*,'(1x,a,8I8)')'NOPEN,NCLSD,NOCC:',NOPEN,NCLSD,NOCC
        NCNF=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP)))
        KCNF=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+1)
        NWRD=ICNFTAB(KCNFINF+3*(LSYM-1+NSYM*(NOPEN-MINOP))+2)
        NCPL=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+1)
        KCPL=ISPNTAB(KSPNINF+6*(NOPEN-MINOP)+3)
CTEST      write(*,'(1x,a,8I8)')'NCNF,KCNF:',NCNF,KCNF
CTEST      write(*,'(1x,a,8I8)')'NCPL,KCPL:',NCPL,KCPL
CTEST      if(NCNF.EQ.0) write(*,*)'    (Skip)'
CTEST      if(NCNF.GT.0) write(*,*)'         Loop over ICNF.'


* Here follows four if-clauses on the cases of IFORM=1..4.
* The contents in all four cases is almost the same, but it is
* *DELIBERATE* that it has not been rewritten as one single piece of code
* with a few IF statement inside. It seems that some compilers do not
* optimize the code correctly then.
* Also, a number of IF-clauses have been replaced by arithmetic computation
* for the same reason. Pardon the clumsy result...
        IWORD = 0 ! dummy initialize
        IF(IFORM.EQ.1) THEN
************************************************************************
* The IFORM=1 case:
C Long loop over configurations
        DO ICNF=1,NCNF
          DO IEL=1,NOCC
            ORBARR(IEL)=ICNFTAB(KCNF-1+IEL+NWRD*(ICNF-1))
          END DO
CTEST          write(*,'(1x,a,8I8)')'Configuration:',
CTEST     &       (ORBARR(IEL),IEL=1,NOCC)
CTEST          write(*,'(1x,a,8I8)')'Nr of spin coupl NCPL:',NCPL
C Loop over spin couplings
          DO ICPL=1,NCPL
            DO I=1,NLEV
              ICASE(I)=0
            END DO
            DO I=1,NCLSD
              IORB=ORBARR(I)
              ICASE(IORB)=3
            END DO
            DO I=1,NOPEN
              IORB=ORBARR(NCLSD+I)
              IFUP=ISPNTAB(KCPL-1+I+NOPEN*(ICPL-1))
*              ICASE(IORB)=1
*              IF(IFUP.NE.1) ICASE(IORB)=2
              ICASE(IORB)=IFUP2CS(IFUP)

            END DO
C A phase factor will be induced by the reordering.
            PHS=1.0D0
            NODD=0
            DO LEV=1,NLEV
              I=ICASE(LEV)
*              IF(I.EQ.1) NODD=1-NODD
*              IF(I.EQ.2) NODD=1-NODD
              NODD=((2*NODD-1)*I*(I-3))/2+NODD
*              IF(NODD.EQ.1) THEN
*                IF(I.EQ.2) PHS=-PHS
*                IF(I.EQ.3) PHS=-PHS
*              END IF
              PHS=DBLE((3+NODD*I*(I-1)*(2*I-7))/3)*PHS
            END DO
CTEST          write(*,'(1x,a,8I8)')
CTEST     &            '      Walk:',(ICASE(IORB),IORB=1,NLEV)
C Pack the walk and add it to the list.
            CALL PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
            IWLKPOS=IWLKPOS+MIPWLK
            NWLKLST=NWLKLST+1
            PHASE(NWLKLST)=PHS
            IF(NWLKLST.EQ.MXWLK) THEN
C Translate the packed walks to a list of CSF ID-numbers
              CALL W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
              IWLKPOS=1
C Loop over this list.
              ICSYMG=NCSYMG
              IF ( IMODE.EQ.0 ) THEN
                DO I=1,NWLKLST
                  PHS=PHASE(I)
                  ICSPLT=ICNUM(I)
                  ICSYMG=NCSYMG+I
                  CINEW(ICSPLT)=CIOLD(ICSYMG)*PHS
CTEST      write(*,'(1x,a,8I8)')'ICSPLT<-ICSYMG:',ICSPLT,ICSYMG
                END DO
              ELSE
                DO I=1,NWLKLST
                  PHS=PHASE(I)
                  ICSPLT=ICNUM(I)
                  ICSYMG=NCSYMG+I
                  CINEW(ICSYMG)=CIOLD(ICSPLT)*PHS
CTEST      write(*,'(1x,a,8I8)')'ICSYMG<-ICSPLT:',ICSYMG,ICSPLT
                END DO
              END IF
              NCSYMG=ICSYMG
              NWLKLST=0
            END IF

C End of loop over spin couplings
          END DO
C End of loop over configurations
        END DO
        ELSE IF(IFORM.EQ.2) THEN
************************************************************************
* The IFORM=2 case:
C Long loop over configurations
        DO ICNF=1,NCNF
          IEL2=0
          IEL1=NCLSD
          DO IORB=1,NLEV
            IOCC=ICNFTAB(KCNF-1+IORB+NWRD*(ICNF-1))
            IF(IOCC.EQ.1) THEN
              IEL1=IEL1+1
              ORBARR(IEL1)=IORB
            ELSE
              IEL2=IEL2+1
              ORBARR(IEL2)=IORB
            END IF
          END DO
CTEST          write(*,'(1x,a,8I8)')'Configuration:',
CTEST     &       (ORBARR(IEL),IEL=1,NOCC)
CTEST          write(*,'(1x,a,8I8)')'Nr of spin coupl NCPL:',NCPL
C Loop over spin couplings
          DO ICPL=1,NCPL
            DO I=1,NLEV
              ICASE(I)=0
            END DO
            DO I=1,NCLSD
              IORB=ORBARR(I)
              ICASE(IORB)=3
            END DO
            DO I=1,NOPEN
              IORB=ORBARR(NCLSD+I)
              IFUP=ISPNTAB(KCPL-1+I+NOPEN*(ICPL-1))
*              ICASE(IORB)=1
*              IF(IFUP.NE.1) ICASE(IORB)=2
              ICASE(IORB)=IFUP2CS(IFUP)
            END DO
C A phase factor will be induced by the reordering.
            PHS=1.0D0
            NODD=0
            DO LEV=1,NLEV
              I=ICASE(LEV)
*              IF(I.EQ.1) NODD=1-NODD
*              IF(I.EQ.2) NODD=1-NODD
              NODD=((2*NODD-1)*I*(I-3))/2+NODD
*              IF(NODD.EQ.1) THEN
*                IF(I.EQ.2) PHS=-PHS
*                IF(I.EQ.3) PHS=-PHS
*              END IF
              PHS=DBLE((3+NODD*I*(I-1)*(2*I-7))/3)*PHS
            END DO
CTEST          write(*,'(1x,a,8I8)')
CTEST     &            '      Walk:',(ICASE(IORB),IORB=1,NLEV)
C Pack the walk and add it to the list.
            CALL PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
            IWLKPOS=IWLKPOS+MIPWLK
            NWLKLST=NWLKLST+1
            PHASE(NWLKLST)=PHS
            IF(NWLKLST.EQ.MXWLK) THEN
C Translate the packed walks to a list of CSF ID-numbers
              CALL W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
              IWLKPOS=1
C Loop over this list.
              ICSYMG=NCSYMG
              IF ( IMODE.EQ.0 ) THEN
                DO I=1,NWLKLST
                  PHS=PHASE(I)
                  ICSPLT=ICNUM(I)
                  ICSYMG=NCSYMG+I
                  CINEW(ICSPLT)=CIOLD(ICSYMG)*PHS
CTEST      write(*,'(1x,a,8I8)')'ICSPLT<-ICSYMG:',ICSPLT,ICSYMG
                END DO
              ELSE
                DO I=1,NWLKLST
                  PHS=PHASE(I)
                  ICSPLT=ICNUM(I)
                  ICSYMG=NCSYMG+I
                  CINEW(ICSYMG)=CIOLD(ICSPLT)*PHS
CTEST      write(*,'(1x,a,8I8)')'ICSYMG<-ICSPLT:',ICSYMG,ICSPLT
                END DO
              END IF
              NCSYMG=ICSYMG
              NWLKLST=0
            END IF

C End of loop over spin couplings
          END DO
C End of loop over configurations
        END DO
        ELSE IF(IFORM.EQ.3) THEN
************************************************************************
* The IFORM=3 case:
C Long loop over configurations
        DO ICNF=1,NCNF
          DO IEL=1,NOCC
            IWRD=(3+IEL)/4
            IREST=(3+IEL)-4*IWRD
            IF(IREST.EQ.0) THEN
              IWORD=ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
            END IF
            IORB=MOD(IWORD,256)
            IWORD=IWORD/256
            ORBARR(IEL)=IORB
          END DO
CTEST          write(*,'(1x,a,8I8)')'Configuration:',
CTEST     &       (ORBARR(IEL),IEL=1,NOCC)
CTEST          write(*,'(1x,a,8I8)')'Nr of spin coupl NCPL:',NCPL
C Loop over spin couplings
          DO ICPL=1,NCPL
            DO I=1,NLEV
              ICASE(I)=0
            END DO
            DO I=1,NCLSD
              IORB=ORBARR(I)
              ICASE(IORB)=3
            END DO
            DO I=1,NOPEN
              IORB=ORBARR(NCLSD+I)
              IFUP=ISPNTAB(KCPL-1+I+NOPEN*(ICPL-1))
*              ICASE(IORB)=1
*              IF(IFUP.NE.1) ICASE(IORB)=2
              ICASE(IORB)=IFUP2CS(IFUP)
            END DO
C A phase factor will be induced by the reordering.
            PHS=1.0D0
            NODD=0
            DO LEV=1,NLEV
              I=ICASE(LEV)
*              IF(I.EQ.1) NODD=1-NODD
*              IF(I.EQ.2) NODD=1-NODD
              NODD=((2*NODD-1)*I*(I-3))/2+NODD
*              IF(NODD.EQ.1) THEN
*                IF(I.EQ.2) PHS=-PHS
*                IF(I.EQ.3) PHS=-PHS
*              END IF
              PHS=DBLE((3+NODD*I*(I-1)*(2*I-7))/3)*PHS
            END DO
CTEST          write(*,'(1x,a,8I8)')
CTEST     &            '      Walk:',(ICASE(IORB),IORB=1,NLEV)
C Pack the walk and add it to the list.
            CALL PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
            IWLKPOS=IWLKPOS+MIPWLK
            NWLKLST=NWLKLST+1
            PHASE(NWLKLST)=PHS
            IF(NWLKLST.EQ.MXWLK) THEN
C Translate the packed walks to a list of CSF ID-numbers
              CALL W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
              IWLKPOS=1
C Loop over this list.
              ICSYMG=NCSYMG
              IF ( IMODE.EQ.0 ) THEN
                DO I=1,NWLKLST
                  PHS=PHASE(I)
                  ICSPLT=ICNUM(I)
                  ICSYMG=NCSYMG+I
                  CINEW(ICSPLT)=CIOLD(ICSYMG)*PHS
CTEST      write(*,'(1x,a,8I8)')'ICSPLT<-ICSYMG:',ICSPLT,ICSYMG
                END DO
              ELSE
                DO I=1,NWLKLST
                  PHS=PHASE(I)
                  ICSPLT=ICNUM(I)
                  ICSYMG=NCSYMG+I
                  CINEW(ICSYMG)=CIOLD(ICSPLT)*PHS
CTEST      write(*,'(1x,a,8I8)')'ICSYMG<-ICSPLT:',ICSYMG,ICSPLT
                END DO
              END IF
              NCSYMG=ICSYMG
              NWLKLST=0
            END IF

C End of loop over spin couplings
          END DO
C End of loop over configurations
        END DO
        ELSE IF(IFORM.EQ.4) THEN
************************************************************************
* The IFORM=4 case:
C Long loop over configurations
        DO ICNF=1,NCNF
          IEL2=0
          IEL1=NCLSD
          DO IORB=1,NLEV
            IWRD=(IORB+14)/15
            IREST=IORB+14-15*IWRD
            IF(IREST.EQ.0) THEN
              IWORD=ICNFTAB(KCNF-1+IWRD+NWRD*(ICNF-1))
            END IF
            IOCC=MOD(IWORD,4)
            IWORD=IWORD/4
            IF(IOCC.EQ.1) THEN
              IEL1=IEL1+1
              ORBARR(IEL1)=IORB
            ELSE
              IEL2=IEL2+1
              ORBARR(IEL2)=IORB
            END IF
          END DO
CTEST          write(*,'(1x,a,8I8)')'Configuration:',
CTEST     &       (ORBARR(IEL),IEL=1,NOCC)
CTEST          write(*,'(1x,a,8I8)')'Nr of spin coupl NCPL:',NCPL
C Loop over spin couplings
          DO ICPL=1,NCPL
            DO I=1,NLEV
              ICASE(I)=0
            END DO
            DO I=1,NCLSD
              IORB=ORBARR(I)
              ICASE(IORB)=3
            END DO
            DO I=1,NOPEN
              IORB=ORBARR(NCLSD+I)
              IFUP=ISPNTAB(KCPL-1+I+NOPEN*(ICPL-1))
*              ICASE(IORB)=1
*              IF(IFUP.NE.1) ICASE(IORB)=2
              ICASE(IORB)=IFUP2CS(IFUP)
            END DO
C A phase factor will be induced by the reordering.
            PHS=1.0D0
            NODD=0
            DO LEV=1,NLEV
              I=ICASE(LEV)
*              IF(I.EQ.1) NODD=1-NODD
*              IF(I.EQ.2) NODD=1-NODD
              NODD=((2*NODD-1)*I*(I-3))/2+NODD
*              IF(NODD.EQ.1) THEN
*                IF(I.EQ.2) PHS=-PHS
*                IF(I.EQ.3) PHS=-PHS
*              END IF
              PHS=DBLE((3+NODD*I*(I-1)*(2*I-7))/3)*PHS
            END DO
CTEST          write(*,'(1x,a,8I8)')
CTEST     &            '      Walk:',(ICASE(IORB),IORB=1,NLEV)
C Pack the walk and add it to the list.
            CALL PKWLK(NLEV,MIPWLK,1,KWALK(IWLKPOS),ICASE)
            IWLKPOS=IWLKPOS+MIPWLK
            NWLKLST=NWLKLST+1
            PHASE(NWLKLST)=PHS
            IF(NWLKLST.EQ.MXWLK) THEN
C Translate the packed walks to a list of CSF ID-numbers
              CALL W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
              IWLKPOS=1
C Loop over this list.
              ICSYMG=NCSYMG
              IF ( IMODE.EQ.0 ) THEN
                DO I=1,NWLKLST
                  PHS=PHASE(I)
                  ICSPLT=ICNUM(I)
                  ICSYMG=NCSYMG+I
                  CINEW(ICSPLT)=CIOLD(ICSYMG)*PHS
CTEST      write(*,'(1x,a,8I8)')'ICSPLT<-ICSYMG:',ICSPLT,ICSYMG
                END DO
              ELSE
                DO I=1,NWLKLST
                  PHS=PHASE(I)
                  ICSPLT=ICNUM(I)
                  ICSYMG=NCSYMG+I
                  CINEW(ICSYMG)=CIOLD(ICSPLT)*PHS
CTEST      write(*,'(1x,a,8I8)')'ICSYMG<-ICSPLT:',ICSYMG,ICSPLT
                END DO
              END IF
              NCSYMG=ICSYMG
              NWLKLST=0
            END IF

C End of loop over spin couplings
          END DO
C End of loop over configurations
        END DO
        END IF
************************************************************************

C End of loop over NOPEN
      END DO
C
C As above, processing what remains in the KWALK buffer.
      CALL W2SGORD(SGS,CIS,MWS2W,NWLKLST,KWALK,ICNUM)
      IWLKPOS=1
      IF ( IMODE.EQ.0 ) THEN
        DO I=1,NWLKLST
          PHS=PHASE(I)
          ICSPLT=ICNUM(I)
          ICSYMG=NCSYMG+I
          CINEW(ICSPLT)=CIOLD(ICSYMG)*PHS
C      write(*,'(1x,a,8I8)')'ICSPLT<-ICSYMG:',ICSPLT,ICSYMG
        END DO
      ELSE
        DO I=1,NWLKLST
          PHS=PHASE(I)
          ICSPLT=ICNUM(I)
CTEST      IF(ICSPLT.LT.0 .or. ICSPLT.GT.NCONF) STOP 'Ohh'
          ICSYMG=NCSYMG+I
          CINEW(ICSYMG)=CIOLD(ICSPLT)*PHS
C      write(*,'(1x,a,8I8)')'ICSYMG<-ICSPLT:',ICSYMG,ICSPLT
        END DO
      END IF
      NCSYMG=NCSYMG+NWLKLST
      NWLKLST=0

      IF( IPGLOB.GE.5 ) THEN
        WRITE(6,*)
        WRITE(6,*)' CI vector reordered in SYG2SGU'
        IF(IMODE.EQ.0) THEN
          WRITE(6,*)' OLD CI vector in SYMMETRIC GROUP order'
          WRITE(6,'(10F12.8)') (CIOLD(I),I=1,NCONF)
          WRITE(6,*)' NEW CI vector in SPLIT GRAPH UGA order'
          WRITE(6,'(10F12.8)') (CINEW(I),I=1,NCONF)
        ELSE
          WRITE(6,*)' OLD CI vector in SPLIT GRAPH UGA order'
          WRITE(6,'(10F12.8)') (CIOLD(I),I=1,NCONF)
          WRITE(6,*)' NEW CI vector in SYMMETRIC GROUP order'
          WRITE(6,'(10F12.8)') (CINEW(I),I=1,NCONF)
        END IF
      ENDIF
C
      CALL mma_deallocate(OrbArr)
      CALL mma_deallocate(MWS2W)

      END SUBROUTINE SYG2SGU

      SUBROUTINE PKWLK(N,IPWLK,NWALK,IWALK,ICASE)
      Integer IWALK(*),ICASE(N,NWALK)
C PURPOSE: PACK THE GUGA STEP NUMBERS INTO THE ARRAY IWALK.
C EACH OF THE NWALK WALKS HAS N STEP NUMBERS, 2 BITS EACH,
C AT MOST 15 TO AN INTEGER ELEMENT OF IWALK, EACH NEW WALK
C IS ALIGNED ON INTEGERS.
C NOTE: Can be used for upper, lower, or total walks, so the
C number of integers used for each walk (IPWLK) is given as
C call parameter.
      IPOS=0
      DO I=1,NWALK
        LEND=0
        DO J=1,IPWLK
          LSTA=LEND+1
          LEND=MIN(LSTA+14,N)
          IPOS=IPOS+1
          IWORD=0
          DO L=LEND,LSTA,-1
            IWORD=4*IWORD+ICASE(L,I)
          END DO
          IWALK(IPOS)=IWORD
        END DO
      END DO
      END SUBROUTINE PKWLK

      SUBROUTINE UPKWLK(N,IPWLK,NWALK,IWALK,ICASE)
      Integer IWALK(*),ICASE(N,NWALK)
* See companion subroutine PKWLK.
      IPOS=0
      DO I=1,NWALK
        LEND=0
        DO J=1,IPWLK
          LSTA=LEND+1
          LEND=MIN(LSTA+14,N)
          IPOS=IPOS+1
          IWORD=IWALK(IPOS)
          DO L=LSTA,LEND
            NEXT=IWORD/4
            ICASE(L,I)=IWORD-4*NEXT
            IWORD=NEXT
          END DO
        END DO
      END DO
      END SUBROUTINE UPKWLK


      SUBROUTINE W2SGORD(SGS,CIS,MWS2W,
     &                 NLIST,KWALK,ICNUM)
      use gugx, only: SGStruct, CIStruct
      use stdalloc, only: mma_allocate, mma_deallocate
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Integer MWS2W(*), NLIST
      Integer KWALK(*),ICNUM(NLIST)

      Integer, PARAMETER :: MXCPI=15
      Integer, Allocatable:: ICS(:)

C Purpose: Given a list of bit-packed total walks,
C translate this into a list of elements of a CI array.
C MXCPI= Max nr of case numbers packed in one integer.
C Dereference SGS:

      nLev   =SGS%nLev
      nVert  =SGS%nVert
      MidLev =SGS%MidLev
      MVSta  =SGS%MVSta
C Dereference CIS:

      nMidV =CIS%nMidV
      NIPWLK=CIS%nIpWlk
C Nr of integers used to store each total walk:
      MIPWLK=1+(NLEV-1)/MXCPI
C Allocate scratch space for case numbers:
      CALL mma_allocate(ICS,NLEV,Label='ICS')
      CALL W2SGORD1(NLEV,NVERT,NMIDV,NIPWLK,SGS%ISM,MIDLEV,
     &            MVSTA,CIS%IOCSF,CIS%NOW,CIS%IOW,
     &            SGS%DOWN,SGS%MAW,ICS,
     &            MWS2W,MIPWLK,NLIST,KWALK,ICNUM)
      CALL mma_deallocate(ICS)
      END SUBROUTINE W2SGORD

      SUBROUTINE W2SGORD1(NLEV,NVERT,NMIDV,NIPWLK,ISM,MIDLEV,
     &                  MVSTA,IOCSF,NOW,IOW,IDOWN,MAW,ICS,
     &                  MWS2W,MIPWLK,NLIST,KWALK,ICNUM)
      use Symmetry_Info, only: nSym=>nIrrep, MUL
      Integer IOCSF(NSYM,NMIDV,NSYM)
      Integer NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      Integer ISM(NLEV),IDOWN(NVERT,0:3),MAW(NVERT,0:3)
      Integer KWALK(MIPWLK,NLIST),ICNUM(NLIST)
      Integer MWS2W(*)
      Integer ICS(NLEV)
C Purpose: For a wave function in Split GUGA storage structure,
C given KWALK(J,I) with J=1..MIPWLK that contains the
C complete Guga walk, as a packed array of case numbers, construct
C ICNUM(I), which is the index of this configuration in the
C Split-Guga scheme.
C Use MAWS to WLK table, MWS2W
CTEST      write(*,*)' In W2SGORD1. NLIST=',NLIST
CTEST      write(*,*)'    NIPWLK:',NIPWLK
CTEST      write(*,*)'    MIPWLK:',MIPWLK
CTEST      write(*,*)' MAW table:'
CTEST      do I=1,NVERT
CTEST       write(*,'(1x,I4,5x,4i4)')I,(MAW(I,J),J=0,3)
CTEST      end do

      DO ICONF=1,NLIST
C Unpack total walk to ICS()
CTEST        write(*,*)' ICONF=',ICONF
CTEST        write(*,*)' KWALK()=',(KWALK(I,ICONF),I=1,MIPWLK)
        CALL UPKWLK(NLEV,MIPWLK,1,KWALK(1,ICONF),ICS)

CTEST        write(*,*)' ICS:'
CTEST        write(*,'(1x,30I3)')(ICS(I),I=1,NLEV)
C Follow upper walk down to MIDLEV:
        MAWSU=0
        IUV=1
        ISYUP=1
        IDV=-1000000000
        DO LEV=NLEV,MIDLEV+1,-1
          IC=ICS(LEV)
          IF(((IC+1)/2).EQ.1)ISYUP=MUL(ISM(LEV),ISYUP)
          IDV=IDOWN(IUV,IC)
          MAWSU=MAWSU+MAW(IUV,IC)
CTEST          write(*,'(1x,4i6)')iuv,ic,idv,maw(iuv,ic)
          IUV=IDV
        END DO
C We have found the midvertex number:
        MV=IDV+1-MVSTA
C Follow lower walk down to MIDLEV:
        ISYDWN=1
        MAWSD=0
        DO LEV=MIDLEV,1,-1
          IC=ICS(LEV)
          IF(((IC+1)/2).EQ.1)ISYDWN=MUL(ISM(LEV),ISYDWN)
          IDV=IDOWN(IUV,IC)
          MAWSD=MAWSD+MAW(IUV,IC)
CTEST          write(*,'(1x,4i6)')iuv,ic,idv,maw(iuv,ic)
          IUV=IDV
        END DO
CTEST        write(*,*)' MAWSU:',MAWSU
CTEST        write(*,*)' MAWSD:',MAWSD
        IUW=MWS2W(MAWSU)
        IDW=MWS2W(MAWSD)
CTEST        write(*,*)' IUW:',IUW
CTEST        write(*,*)' IDW:',IDW
C Subtract the offsets:
CTEST      write(6,*)' Offset, IOW(1,ISYUP,MV):',IOW(1,ISYUP,MV)
CTEST      write(6,*)' Offset, IOW(2,ISYDWN,MV):',IOW(2,ISYDWN,MV)
CTEST      write(6,*)' NIPWLK:',NIPWLK
CTEST      write(6,*)' MIPWLK:',MIPWLK
        IUW=IUW-IOW(1,ISYUP,MV)/NIPWLK
        IDW=IDW-IOW(2,ISYDWN,MV)/NIPWLK
CTEST      write(6,*)' With offsets subtracted:'
CTEST      write(6,*)' IUW, IDW:',IUW,IDW
C Split-Guga storage scheme: an element in a set of matrices.
C Offset to the matrix we want is IOCSF(ISYUP,MV,ISYCI).
C Leading dimension=nr of upwalks in this block.
        ISYCI=MUL(ISYUP,ISYDWN)
        IOFF=IOCSF(ISYUP,MV,ISYCI)
        LDIM=NOW(1,ISYUP,MV)
        ICNUM(ICONF)=IOFF+IUW+LDIM*(IDW-1)
      END DO
      END SUBROUTINE W2SGORD1

      SUBROUTINE MSTOW(SGS,CIS,MWS2W,nSym)
      use gugx, only: SGStruct, CIStruct
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Integer MWS2W(*)
      Integer nSym

      Integer, allocatable:: ICS(:)

      NLEV  =SGS%nLev
      NVERT =SGS%nVert
      MIDLEV=SGS%MidLev

      NMIDV =CIS%nMidV
      NIPWLK=CIS%nIpWlk
      NWALK =CIS%nWalk
      CALL mma_allocate(ICS,NLEV,Label='ICS')
      CALL MSTOW1(NSYM,NLEV,NVERT,NMIDV,NIPWLK,NWALK,
     &            MIDLEV,ICS,CIS%NOW,CIS%IOW,
     &            CIS%ICase,SGS%UP,SGS%DOWN,
     &            SGS%MAW,MWS2W)
      CALL mma_deallocate(ICS)
      END SUBROUTINE MSTOW

      SUBROUTINE MSTOW1(NSYM,NLEV,NVERT,NMIDV,NIPWLK,NWALK,
     &                  MIDLEV,ICS,NOW,IOW,IWALK,
     &                  IUP,IDOWN,MAW,MWS2W)
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer ICS(NLEV)
      Integer NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      Integer IDOWN(NVERT,0:3),IUP(NVERT,0:3)
      Integer MAW(NVERT,0:3)
      Integer IWALK(NIPWLK*NWALK)
      Integer MWS2W(NWALK)

C Purpose: From the list of packed up- and downwalks, construct
C the table MWS2W, such that MAW sums can be translated to the
C corresponding walks of the Split-GUGA.

      DO MV=1,NMIDV
        DO ISYUP=1,NSYM
          NUP=NOW(1,ISYUP,MV)
          IF(NUP.EQ.0) cycle
          IUOFF=IOW(1,ISYUP,MV)/NIPWLK
          DO IUW=1,NUP
            IUWTOT=IUOFF+IUW
C Unpack upper walk to ICS()
            CALL UPKWLK(NLEV-MIDLEV,NIPWLK,1,IWALK(1+NIPWLK*(IUWTOT-1)),
     &                  ICS(MIDLEV+1))
            MS=0
            IUV=1
            DO LEV=NLEV,MIDLEV+1,-1
              IC=ICS(LEV)
              MS=MS+MAW(IUV,IC)
              IUV=IDOWN(IUV,IC)
            END DO
            MWS2W(MS)=IUWTOT
          END DO
        END DO
      END DO

      DO MV=1,NMIDV
        DO ISYDWN=1,NSYM
          NDWN=NOW(2,ISYDWN,MV)
          IF(NDWN.EQ.0) cycle
          IDOFF=IOW(2,ISYDWN,MV)/NIPWLK
          DO IDW=1,NDWN
            IDWTOT=IDOFF+IDW
C Unpack lower walk to ICS()
            CALL UPKWLK(MIDLEV,NIPWLK,1,IWALK(1+NIPWLK*(IDWTOT-1)),ICS)
            MS=0
            IDV=NVERT
            DO LEV=1,MIDLEV
              IC=ICS(LEV)
              IUV=IUP(IDV,IC)
              MS=MS+MAW(IUV,IC)
              IDV=IUV
            END DO
            MWS2W(MS)=IDWTOT
          END DO
        END DO
      END DO

      END SUBROUTINE MSTOW1
