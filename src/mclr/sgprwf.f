!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE SGPRWF_MCLR(SGS,CIS,LSYM,PRWTHR,NCONF,CI)
!
!     PURPOSE: PRINT THE WAVEFUNCTION (SPIN COUPLING AND OCCUPATIONS)
!
!     NOTE:    THIS ROUTINE USES THE SPLIT GRAPH GUGA CONVENTION, I.E.,
!              CI BLOCKS ARE MATRICES CI(I,J), WHERE THE  FIRST INDEX
!              REFERS TO THE UPPER PART OF THE WALK.
!
      use struct, only: SGStruct, CIStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      Type(SGStruct) SGS
      Type(CIStruct) CIS
      Integer  LSYM
      Real*8 PRWTHR
      Integer  NCONF
      Real*8 CI(NCONF)


      Integer ICS(50)
      Character(LEN=120) Line
!
      Associate(nSym=>SGS%nSym, nLev=>SGS%nLev, MidLev=>SGS%MidLev,
     &           nMidV=>CIS%nMidV, nIpWlk=>CIS%nIpWlk,
     &           NSM=>SGS%ISm, NOCSF=>CIS%NOCSF, IOCSF=>CIS%IOCSF,
     &           NOW=>CIS%NOW, IOW=>CIS%IOW, ICASE=>CIS%ICASE)

      Line(1:16)='      conf/sym  '
      iOff=16
      iSym=nSm(1)
      Do Lev=1,nLev
         If ( nSm(Lev).ne.iSym ) iOff=iOff+1
         Write (Line(iOff+Lev:),'(I1)') nSm(Lev)
         If ( nSm(Lev).ne.iSym ) iSym=nSm(Lev)
      End Do
      iOff=iOff+nLev+3
      Line(iOff:iOff+15)='   Coeff  Weight'
      Write (6,'(A)') Line(1:iOff+15)
      Line=' '
!
!     THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
!     WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
!

      DO MV=1,NMIDV
        DO ISYUP=1,NSYM
          NCI=NOCSF(ISYUP,MV,LSYM)
          IF(NCI.EQ.0) Cycle
          NUP=NOW(1,ISYUP,MV)
          ISYDWN=1+IEOR(ISYUP-1,LSYM-1)
          NDWN=NOW(2,ISYDWN,MV)
          ICONF=IOCSF(ISYUP,MV,LSYM)
          IUW0=1-NIPWLK+IOW(1,ISYUP,MV)
          IDW0=1-NIPWLK+IOW(2,ISYDWN,MV)
          IDWNSV=0
          DO IDWN=1,NDWN
            DO IUP=1,NUP
              ICONF=ICONF+1
              COEF=CI(ICONF)
! -- SKIP OR PRINT IT OUT?
              IF(ABS(COEF).LT.PRWTHR) Cycle
              IF(IDWNSV.NE.IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=ICASE(ICDPOS)
! -- UNPACK LOWER WALK.
                NNN=0
                DO LEV=1,MIDLEV
                  NNN=NNN+1
                  IF(NNN.EQ.16) THEN
                    NNN=1
                    ICDPOS=ICDPOS+1
                    ICDWN=ICASE(ICDPOS)
                  END IF
                  IC1=ICDWN/4
                  ICS(LEV)=ICDWN-4*IC1
                  ICDWN=IC1
                END DO
                IDWNSV=IDWN
              END IF
              ICUPOS=IUW0+NIPWLK*IUP
              ICUP=ICASE(ICUPOS)
! -- UNPACK UPPER WALK:
              NNN=0
              DO LEV=MIDLEV+1,NLEV
                NNN=NNN+1
                IF(NNN.EQ.16) THEN
                  NNN=1
                  ICUPOS=ICUPOS+1
                  ICUP=ICASE(ICUPOS)
                END IF
                IC1=ICUP/4
                ICS(LEV)=ICUP-4*IC1
                ICUP=IC1
              END DO
! -- PRINT IT!
              Write (Line(1:),'(I8)') iConf
              iOff=10
              iSym=nSm(1)
              Do Lev=1,nLev
                 If ( nSm(Lev).ne.iSym ) iOff=iOff+1

                 Select case (ICS(Lev))
                    Case (3)
                       Write (Line(iOff+Lev:),'(A1)') '2'
                    Case (2)
                       Write (Line(iOff+Lev:),'(A1)') 'd'
                    Case (1)
                       Write (Line(iOff+Lev:),'(A1)') 'u'
                    Case (0)
                       Write (Line(iOff+Lev:),'(A1)') '0'
                    Case Default
                       Call Abend()
                 End Select
                 If ( nSm(Lev).ne.iSym ) iSym=nSm(Lev)
              End Do
              iOff=iOff+nLev+3
              Write (Line(iOff:),'(2F8.5)') COEF,COEF**2
              Write (6,'(6X,A)') Line(1:iOff+15)
              Line=' '
            END DO
          END DO
        END DO
      END DO
!
!     EXIT
!
      End Associate
      END
