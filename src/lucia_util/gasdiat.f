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
      SUBROUTINE GASDIAT(    DIAG,   LUDIA,   ECORE,  ICISTR,     I12,
     &                      IBLTP,  NBLOCK,  IBLKFO)
      use stdalloc, only: mma_allocate, mma_deallocate
      use strbas, only: NSTSO
      use lucia_data, only: IPRDIA
      use lucia_data, only: PSSIGN
      use lucia_data, only: MXNSTR,I_AM_OUT,N_ELIMINATED_BATCHES
      use lucia_data, only: IDISK
      use lucia_data, only: NTOOB,IREOST,IREOTS,NACOB
      use lucia_data, only: NOCTYP
      use lucia_data, only: NELEC
#ifdef _DEBUGPRINT_
      use lucia_data, only: IBSPGPFTP
#endif
      use csm_data, only: NSMST

*
* CI diagonal in SD basis for state with symmetry ISM in internal
* space ISPC
*
* GAS version, Winter of 95
*
* Driven by table of TTS blocks, May97
*
      IMPLICIT NONE
* =====
*.Input
* =====
*
      INTEGER LUDIA,ICISTR,I12,NBLOCK
      REAL*8 ECORE
      INTEGER IBLTP(*)
      INTEGER IBLKFO(8,NBLOCK)

*
* ======
*.Output
* ======
      REAL*8 DIAG(*)

      Integer, Allocatable:: LASTR(:), LBSTR(:)
      Real*8, Allocatable:: LSCR2(:)
      Real*8, Allocatable:: LJ(:), LK(:), LXB(:), LH1D(:), LRJKA(:)
      INTEGER, EXTERNAL:: IMNMX
      INTEGER NTEST,IATP,IBTP,NAEL,NBEL,NOCTPA,MAXA
#ifdef _DEBUGPRINT_
      Integer NOCTPB,IOCTPA,IOCTPB
#endif
*
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRDIA)
*
** Specifications of internal space
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
*
#ifdef _DEBUGPRINT_
      NOCTPB = NOCTYP(IBTP)
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' GASDIA speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
      END IF
#endif
*
**. Local memory
*
      CALL mma_allocate(LJ   ,NTOOB**2,Label='LJ')
      CALL mma_allocate(LK   ,NTOOB**2,Label='LK')
      Call mma_allocate(LSCR2,2*NTOOB**2,Label='LSCR2')
      CALL mma_allocate(LXB  ,NACOB,Label='LXB')
      CALL mma_allocate(LH1D ,NACOB,Label='LH1D')
*. Space for blocks of strings
      Call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
      Call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
      MAXA = IMNMX(NSTSO(IATP)%I,NSMST*NOCTPA,2)
      CALL mma_allocate(LRJKA,MAXA,Label='LRJKA')
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      CALL GT1DIA(LH1D)
      CALL GTJK(LJ,LK,NTOOB,LSCR2,IREOTS,IREOST)
      IF( LUDIA .GT. 0 ) IDISK(LUDIA)=0
      CALL GASDIAS(NAEL,LASTR,NBEL,LBSTR,
     &             NACOB,DIAG,NSMST,
     &             LH1D,LXB,LJ,LK,
     &             NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &             LUDIA,ECORE,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &             LRJKA,I12,IBLTP,NBLOCK,IBLKFO,
     &             I_AM_OUT,N_ELIMINATED_BATCHES)
*.Flush local memory
      Call mma_deallocate(LJ)
      Call mma_deallocate(LK)
      Call mma_deallocate(LSCR2)
      Call mma_deallocate(LXB)
      Call mma_deallocate(LH1D)
      Call mma_deallocate(LASTR)
      Call mma_deallocate(LBSTR)
      Call mma_deallocate(LRJKA)
*
      END SUBROUTINE GASDIAT
