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
      SUBROUTINE TODSCP(A,NDIM,MBLOCK,IFIL)
*
C TRANSFER ARRAY REAL*8  A(LENGTH NDIM) TO DISCFIL IFIL IN
C RECORDS WITH LENGTH NBLOCK.
*
* Packed version : Store only nonzero elements
*. Small elements should be xeroed outside
      use Constants, only:Zero
      use lucia_data, only: IDISK
      IMPLICIT NONE
      INTEGER NDIM,MBLOCK,IFIL
      REAL*8 A(*)
C-jwk-cleanup      INTEGER START,STOP
      REAL*8, EXTERNAL :: INPROD
      INTEGER ISCR(2), IDUMMY(1)
*
      INTEGER, PARAMETER :: LPBLK=50000
      INTEGER IPAK(LPBLK)
      REAL*8 XPAK(LPBLK)
      INTEGER IPACK,IMZERO,MMBLOCK,IELMNT,LBATCH
      REAL*8 XNORM
*
*
C?    write(6,*) ' entering TODSCP, file = ', IFIL
C?    CALL XFLUSH(6)
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = INPROD(A,A,NDIM)
        IF(XNORM.EQ.Zero) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. Packing
        ISCR(2) = 1
C       CALL ITODS(ISCR,2,MMBLOCK,IFIL)
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
*. Loop over packed records of dimension LPBLK
      IELMNT = 0
 1000 CONTINUE
*. The next LPBLK elements
      LBATCH = 0
*. Obtain next batch of elemnts
  999 CONTINUE
       IF(NDIM.GE.1) THEN
       IELMNT = IELMNT+1
       IF(A(IELMNT).NE.ZERO) THEN
         LBATCH=LBATCH+1
         IPAK(LBATCH) = IELMNT
         XPAK(LBATCH) = A(IELMNT)
       END IF
       END IF
       IF(LBATCH.EQ.LPBLK.OR.IELMNT.EQ.NDIM) goto 998
       GOTO 999
*. Send to DISC
 998   CONTINUE
       IDUMMY(1)=LBATCH
       CALL IDAFILE(IFIL,1,IDUMMY,1,IDISK(IFIL))
       IF(LBATCH.GT.0) THEN
         CALL IDAFILE(IFIL,1,IPAK,LBATCH,IDISK(IFIL))
         CALL DDAFILE(IFIL,1,XPAK,LBATCH,IDISK(IFIL))
       END IF
       IF(IELMNT.EQ.NDIM) THEN
         IDUMMY(1)=-1
         CALL IDAFILE(IFIL,1,IDUMMY,1,IDISK(IFIL))
       ELSE
         IDUMMY(1)=0
         CALL IDAFILE(IFIL,1,IDUMMY,1,IDISK(IFIL))
         GOTO 1000
       END IF
*. End of loop over records of truncated elements
c      END IF
 1001 CONTINUE
*
      END SUBROUTINE TODSCP
