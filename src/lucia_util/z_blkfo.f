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
* Copyright (C) 1998, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE Z_BLKFO(ISPC,ISM,    IATP,    IBTP,
     &                    KPCI1BT,  KPCIBT, KPCBLTP,  NBATCH,
     &                     NBLOCK)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Local_Arrays, only: CLBT, CLEBT
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM,IATP,IBTP.
*
* Output is given in the form of pointers to vectors in WORK
* where the info is stored :
*
* CLBT : Length of each Batch ( in blocks)
* CLEBT : Length of each Batch ( in elements)
* KPCI1BT : Length of each block
* KPCIBT  : Info on each block
* KPCBLTP : BLock type for each symmetry
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "cicisp.fh"
#include "stinf.fh"
#include "cstate.fh"
#include "csm.fh"
#include "strbas.fh"
#include "crun.fh"
      Integer, Allocatable:: LCIOIO(:)
      Integer, Allocatable:: SVST(:)
*
* Some dummy initializations
      NTEST = 00
#ifdef _DEBUGPRINT_
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' Output from Z_BLKFO '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
#endif
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Pointers to output arrays
      CALL mma_allocate(CLBT ,MXNTTS,Label='CLBT')
      CALL mma_allocate(CLEBT,MXNTTS,Label='CLEBT')
      CALL GETMEM('CI1BT ','ALLO','INTE',KPCI1BT,MXNTTS)
      CALL GETMEM('CIBT  ','ALLO','INTE',KPCIBT ,8*MXNTTS)
      CALL GETMEM('CBLTP ','ALLO','INTE',KPCBLTP,NSMST)
*.    ^ These should be preserved after exit so put mark for flushing here
*. Info needed for generation of block info
      Call mma_allocate(LCIOIO,NOCTPA*NOCTPB,Label='LCIOIO')
      CALL IAIBCM(ISPC,LCIOIO)
      Call mma_allocate(SVST,1,Label='SVST')
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,IWORK(KPCBLTP),SVST)
      Call mma_deallocate(SVST)
*. Allowed length of each batch
c      IF(ISIMSYM.EQ.0) THEN
        LBLOCK = MXSOOB
c      ELSE
c        LBLOCK = MXSOOB_AS
c      END IF
*
      LBLOCK = MAX(LBLOCK,LCSBLK)
* JESPER : Should reduce I/O
      IF (ENVIRO(1:6).EQ.'RASSCF') THEN
         LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
         IF(PSSIGN.NE.0.0D0) LBLOCK = INT(2.0D0*XISPSM(IREFSM,1))
      ENDIF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' LBLOCK = ', LBLOCK
      END IF
*
*. Batches  of C vector
      CALL PART_CIV2(      IDC,
     &               IWORK(KPCBLTP),
     &               IWORK(KNSTSO(IATP)),
     &               IWORK(KNSTSO(IBTP)),NOCTPA,NOCTPB, NSMST,LBLOCK,
     &               LCIOIO,
*
     &               ISMOST(1,ISM),
     &               NBATCH,
     &               CLBT,
     &               CLEBT,IWORK(KPCI1BT),
     &               IWORK(KPCIBT),0,ISIMSYM)
*. Number of BLOCKS
      NBLOCK = IFRMR(IWORK(KPCI1BT),1,NBATCH)
     &       + IFRMR(CLBT,1,NBATCH) - 1
      IF(NTEST.GE.1) THEN
         WRITE(6,*) ' Number of batches', NBATCH
         WRITE(6,*) ' Number of blocks ', NBLOCK
      END IF
*. Length of each block
      CALL EXTRROW(IWORK(KPCIBT),8,8,NBLOCK,IWORK(KPCI1BT))
*
      Call mma_deallocate(LCIOIO)
      RETURN
      END
