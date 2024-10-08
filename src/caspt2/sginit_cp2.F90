!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994,2006, Per Ake Malmqvist                           *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
! 2006  PER-AAKE MALMQUIST                   *
!--------------------------------------------*
      SUBROUTINE SGINIT_CP2(nSym,iSpin,nActEl,nHole1,nEle3,nRas1T,nRas2T,nRas3T,SGS,CIS,EXS)
      use stdalloc, only: mma_deallocate
      use gugx, only: SGStruct, CIStruct, EXStruct
      use MkGUGA_mod, only: MKGUGA
      IMPLICIT None
      Integer nSym,iSpin,nActEl,nHole1,nEle3,nRas1T,nRas2T,nRas3T
      Type(SGStruct) SGS
      Type(CIStruct) CIS
      Type(EXStruct) EXS

      SGS%nSym=nSym
      SGS%iSpin=iSpin
      SGS%nActEl=nActEl

      Associate ( LM1RAS=>SGS%LM1RAS, LM3RAS=>SGS%LM3RAS,               &
     &            LV1RAS=>SGS%LV1RAS, LV3RAS=>SGS%LV3RAS,               &
     &            IFRAS=>SGS%IFRAS)

      LV1RAS=NRAS1T
      LV3RAS=nRas1T+NRAS2T
      LM1RAS=2*nRas1T-NHOLE1
      LM3RAS=NACTEL-NELE3
      IF ((NRAS1T+NRAS3T)/=0) Then
         IFRAS=1
      Else
         IFRAS=0
      End If

      CALL MKGUGA(SGS,CIS)

!     FORM VARIOUS OFFSET TABLES:
!     NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!           TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.
!
      CALL MKCOT(SGS,CIS)
!
!     CONSTRUCT THE CASE LIST
!
      Call MKCLIST(SGS,CIS)

! DECIDE MIDLEV AND CALCULATE MODIFIED ARC WEIGHT TABLE.

      CALL MKMAW(SGS)

! THE DAW, UP AND RAW TABLES WILL NOT BE NEEDED ANY MORE:

! CALCULATE SEGMENT VALUES. ALSO, MVL AND MVR TABLES.

      CALL MKSEG(SGS,CIS,EXS)

! NIPWLK: NR OF INTEGERS USED TO PACK EACH UP- OR DOWNWALK.

      CALL NRCOUP(SGS,CIS,EXS)

      CALL MKCOUP(SGS,CIS,EXS)


      Call mma_deallocate(CIS%ISGM)
      Call mma_deallocate(CIS%VSGM)
      Call mma_deallocate(CIS%IVR)

      Call mma_deallocate(SGS%MAW)

      CALL mma_deallocate(SGS%DRT)
      Call mma_deallocate(SGS%DOWN)
      CALL mma_deallocate(SGS%DAW)
      CALL mma_deallocate(SGS%UP)
      CALL mma_deallocate(SGS%RAW)
      Call mma_deallocate(SGS%LTV)

      End Associate

      END SUBROUTINE SGINIT_CP2

