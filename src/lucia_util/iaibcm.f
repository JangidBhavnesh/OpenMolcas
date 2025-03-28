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
* Copyright (C) 1995,1999, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE IAIBCM(ICISPC,IAIB)
*
* obtain allowed combinbation of alpha- and beta- supergroups
* for CI space ICISPC
*
* Master for IAIBCM_GAS
*
*      Jeppe Olsen, august 1995
*                   I_RE_MS2 added, May 99
*
      use lucia_data, only: ICMBSPC,IGSOCCX,LCMBSPC,NGAS
      use lucia_data, only: IPRDIA
      use lucia_data, only: I_RE_MS2_SPACE,I_RE_MS2_VALUE
      use lucia_data, only: IBSPGPFTP,ISPGPFTP,NELFGP
      use lucia_data, only: NOCTYP
      use lucia_data, only: MXPNGAS
      IMPLICIT NONE
      INTEGER ICISPC
*. Output
      INTEGER IAIB(*)

      INTEGER IATP,IBTP,NOCTPA,NOCTPB,IOCTPA,IOCTPB
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
C?    write(6,*) ' IAIB ::::::'
C?    write(6,*) ' LCMBSPC, ICISPC, ICMBSPC '
C?    WRITE(6,*) ICISPC,  LCMBSPC(ICISPC)
C?    WRITE(6,*) (ICMBSPC(II,ICISPC),II=1, LCMBSPC(ICISPC))

      CALL IAIBCM_GAS(LCMBSPC(ICISPC),
     &                ICMBSPC(1,ICISPC),
     &                  IGSOCCX,
     &                   NOCTPA,
     &                   NOCTPB,
*
     &                ISPGPFTP(1,IOCTPA),
     &                ISPGPFTP(1,IOCTPB),NELFGP,MXPNGAS, NGAS, IAIB,
     &                   IPRDIA,I_RE_MS2_SPACE,I_RE_MS2_VALUE)
*
      END SUBROUTINE IAIBCM
