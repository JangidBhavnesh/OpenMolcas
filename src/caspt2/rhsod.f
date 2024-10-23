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
* Copyright (C) Steven Vancoillie                                      *
************************************************************************
*SVC: compute RHS elements "on demand". If we have access to all the
* Cholesky vectors, we can just instruct a process to compute its own
* block of RHS elements, computing the integrals directly. This is much
* more computationally intensive, but should scale much better since we
* go from a badly scaling scatter algorithm to no communication at all.
* This also eliminates the need for the GA library in creating the RHS.
* FIXME: optimizations needed, remove double computation of integrals

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD(IVEC)
      use caspt2_output, only:iPrGlb
      use PrintLevel, only: verbose
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use EQSOLV
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"


      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,'(1X,A)') ' Using RHS on-demand algorithm'
      END IF

#ifdef _MOLCAS_MPP_
      IF (.NOT.Is_Real_Par()) THEN
        WRITE(6,'(1X,A)') 'RHSOD: error: fake parallel not supported'
        CALL AbEnd()
      END IF
#endif

      CALL RHSOD_A(IVEC)
      CALL RHSOD_B(IVEC)
      CALL RHSOD_C(IVEC)
      CALL RHSOD_D(IVEC)
      CALL RHSOD_E(IVEC)
      CALL RHSOD_F(IVEC)
      CALL RHSOD_G(IVEC)
      CALL RHSOD_H(IVEC)

#ifdef _DEBUGPRINT_
* compute and print RHS fingerprints
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') 'Case','Sym','Fingerprint'
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') '====','===','==========='
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NAS*NIS.NE.0) THEN
            CALL RHS_ALLO (NAS,NIS,lg_W)
            CALL RHS_READ (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
            DNRM2 = RHS_DDOT(NAS,NIS,lg_W,lg_W)
            WRITE(6,'(1X,I4,1X,I3,1X,F18.11)') ICASE,ISYM,DNRM2
          END IF
        END DO
      END DO
#endif


      END


************************************************************************
* SUBROUTINES FOR THE SEPARATE CASES
************************************************************************

*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      SUBROUTINE RHSOD_A(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      use caspt2_output, only:iPrGlb
      use PrintLevel, only: debug
      use caspt2_data, only: FIMO
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER IVEC

      INTEGER IOBRA(8,8), IOKET(8,8)
      REAL*8, ALLOCATABLE:: BRA(:), KET(:)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case A'
      END IF

************************************************************************
* Case A:
C   RHS(tvx,j)=(tj,vx)+FIMO(t,j)*kron(v,x)/NACTEL
************************************************************************

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(1,NBRA,IOBRA)
      CALL CHOVEC_SIZE(2,NKET,IOKET)

      CALL mma_allocate(BRA,NBRA,LABEL='BRA')
      CALL mma_allocate(KET,NKET,LABEL='KET')

      CALL CHOVEC_READ(1,BRA,NBRA)
      CALL CHOVEC_READ(2,KET,NKET)

      ICASE=1
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      NFIMOES=0
      DO ISYM=1,NSYM

        NAS=NTUV(ISYM) !NASUP(ISYM,ICASE)
        NIS=NISH(ISYM) !NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 1

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IJ=IISTA,IIEND
          ISYJ=ISYM
          DO ITVX=IASTA,IAEND ! these are always all elements
            ITVXTOT=ITVX+NTUVES(ISYM)
            ITABS=MTUV(1,ITVXTOT)
            IVABS=MTUV(2,ITVXTOT)
            IXABS=MTUV(3,ITVXTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
            IX  =MTREL(1,IXABS)
            ISYX=MTREL(2,IXABS)
! compute integrals (tiuv)
            NV=NVTOT_CHOSYM(MUL(ISYT,ISYJ)) ! JSYM=ISYT*ISYI=ISYU*ISYV
            ITJ=IT-1+NASH(ISYT)*(IJ-1)
            IVX=IV-1+NASH(ISYV)*(IX-1)
            IOFFTJ=1+IOBRA(ISYT,ISYJ)+NV*ITJ
            IOFFVX=1+IOKET(ISYV,ISYX)+NV*IVX
            TJVX=DDOT_(NV,BRA(IOFFTJ),1,KET(IOFFVX),1)
! A(tvx,j) = (tjvx) + FIMO(t,j)*delta(v,x)/NACTEL
            IF (ISYT.EQ.ISYJ.AND.IVABS.EQ.IXABS) THEN
              ITTOT=IT+NISH(ISYT)
              FTJ=FIMO(NFIMOES+(ITTOT*(ITTOT-1))/2+IJ)
              ATVXJ=TJVX+FTJ/DBLE(MAX(1,NACTEL))
            ELSE
              ATVXJ=TJVX
            END IF
! write element A(tvx,j)
            IDX=ITVX+NAS*(IJ-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=ATVXJ
#else
            GA_Arrays(lg_w)%Array(IDX)=ATVXJ
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 1      CONTINUE

        NFIMOES=NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2

      END DO
************************************************************************

      CALL mma_deallocate(BRA)
      CALL mma_deallocate(KET)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_C(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      use caspt2_output, only:iPrGlb
      use PrintLevel, only: debug
      use caspt2_data, only: FIMO
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER IVEC

      INTEGER IOBRA(8,8), IOKET(8,8)
      REAL*8, ALLOCATABLE:: BRA(:), KET(:)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case C'
      END IF

************************************************************************
* Case C:
C   RHS(tvx,a)=(at,vx)+(FIMO(a,t)-Sum_u(au,ut))*delta(v,x)/NACTEL
************************************************************************

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(3,NBRA,IOBRA)
      CALL CHOVEC_SIZE(2,NKET,IOKET)

      CALL mma_allocate(BRA,NBRA,LABEL='BRA')
      CALL mma_allocate(KET,NKET,LABEL='KET')

      CALL CHOVEC_READ(3,BRA,NBRA)
      CALL CHOVEC_READ(2,KET,NKET)

      ICASE=4
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      NFIMOES=0
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE) !NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE) !NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 4

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IA=IISTA,IIEND
          ISYA=ISYM
          DO ITVX=IASTA,IAEND ! these are always all elements
            ITVXTOT=ITVX+NTUVES(ISYM)
            ITABS=MTUV(1,ITVXTOT)
            IVABS=MTUV(2,ITVXTOT)
            IXABS=MTUV(3,ITVXTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
            IX  =MTREL(1,IXABS)
            ISYX=MTREL(2,IXABS)
! compute integrals (at,vx)
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYT)) ! JSYM=ISYT*ISYI=ISYU*ISYV
            IAT=IA-1+NSSH(ISYA)*(IT-1)
            IVX=IV-1+NASH(ISYV)*(IX-1)
            IOFFAT=1+IOBRA(ISYA,ISYT)+NV*IAT
            IOFFVX=1+IOKET(ISYV,ISYX)+NV*IVX
            ATVX=DDOT_(NV,BRA(IOFFAT),1,KET(IOFFVX),1)

! W(tvx,a) = (at,vx) + (FIMO(a,t)-Sum_u(au,ut))*delta(v,x)/NACTEL
! write element W(tvx,j), only the (at,vx) part
            IDX=ITVX+NAS*(IA-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=ATVX
#else
            GA_arrays(lg_w)%Array(IDX)=ATVX
#endif
          END DO
! now, add in the part with corrections to the integrals
          IATOT=IA+NISH(ISYM)+NASH(ISYM)
          DO IT=1,NASH(ISYM)
            ITTOT=IT+NISH(ISYM)
            FAT=FIMO(NFIMOES+(IATOT*(IATOT-1))/2+ITTOT)
            SUMU=0.0D0
            ITABS=NAES(ISYM)+IT
            DO IUABS=1,NASHT
              IUUT=KTUV(IUABS,IUABS,ITABS)-NTUVES(ISYM)
              IDX=IUUT+NAS*(IA-IISTA)
#ifdef _MOLCAS_MPP_
              SUMU=SUMU+DBL_MB(MW+IDX-1)
#else
              SUMU=SUMU+GA_Arrays(lg_W)%Array(IDX)
#endif
            END DO
            ADDONE=(FAT-SUMU)/DBLE(MAX(1,NACTEL))
            DO IVABS=1,NASHT
              ITVV=KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
              IDX=ITVV+NAS*(IA-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=DBL_MB(MW+IDX-1)+ADDONE
#else
              GA_Arrays(lg_w)%Array(IDX)=GA_Arrays(lg_w)%Array(IDX)
     &                                  +ADDONE
#endif
            END DO
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 4      CONTINUE

        NFIMOES=NFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2

      END DO
************************************************************************

      CALL mma_deallocate(BRA)
      CALL mma_deallocate(KET)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_B(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      use caspt2_output, only:iPrGlb
      use PrintLevel, only: debug
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER IVEC

      INTEGER IOSYM(8,8)
      REAL*8, ALLOCATABLE:: CHOBUF(:)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case B'
      END IF

************************************************************************
* Case B (2,3):
C   Let  W(tv,j,l)=(jt,lv):
C   BP(tv,jl)=((tj,vl)+(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
C   BM(tv,jl)=((tj,vl)-(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
************************************************************************

      SQRTH=SQRT(0.5D0)

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(1,NCHOBUF,IOSYM)

      CALL mma_allocate(CHOBUF,NCHOBUF,LABEL='CHOBUF')

      CALL CHOVEC_READ(1,CHOBUF,NCHOBUF)

      iCASE=2
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 2

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IJGEL=IISTA,IIEND
          IJGELTOT=IJGEL+NIGEJES(ISYM)
          IJABS=MIGEJ(1,IJGELTOT)
          ILABS=MIGEJ(2,IJGELTOT)
          IJ  =MIREL(1,IJABS)
          ISYJ=MIREL(2,IJABS)
          IL  =MIREL(1,ILABS)
          ISYL=MIREL(2,ILABS)
          DO ITGEU=IASTA,IAEND ! these are always all elements
            ITGEUTOT=ITGEU+NTGEUES(ISYM)
            ITABS=MTGEU(1,ITGEUTOT)
            IVABS=MTGEU(2,ITGEUTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integrals (ajcl) and (alcj)
            NV=NVTOT_CHOSYM(MUL(ISYT,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
            ITJ=IT-1+NASH(ISYT)*(IJ-1)
            IVL=IV-1+NASH(ISYV)*(IL-1)
            IOFFTJ=1+IOSYM(ISYT,ISYJ)+NV*ITJ
            IOFFVL=1+IOSYM(ISYV,ISYL)+NV*IVL
            TJVL=DDOT_(NV,CHOBUF(IOFFTJ),1,CHOBUF(IOFFVL),1)

            NV=NVTOT_CHOSYM(MUL(ISYT,ISYL))
            ITL=IT-1+NASH(ISYT)*(IL-1)
            IVJ=IV-1+NASH(ISYV)*(IJ-1)
            IOFFTL=1+IOSYM(ISYT,ISYL)+NV*ITL
            IOFFVJ=1+IOSYM(ISYV,ISYJ)+NV*IVJ
            TLVJ=DDOT_(NV,CHOBUF(IOFFTL),1,CHOBUF(IOFFVJ),1)

! BP(tv,jl)=((tj,vl)+(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
            SCL=0.5D0
            IF (ITABS.EQ.IVABS) SCL=SCL*0.5D0
            IF (ILABS.EQ.IJABS) SCL=SCL*SQRTH
            BPTVJL=SCL*(TJVL+TLVJ)
! write element HP(ac,jl)
            IDX=ITGEU+NAS*(IJGEL-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=BPTVJL
#else
            GA_Arrays(lg_w)%Array(IDX)=BPTVJL
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 2      CONTINUE
      END DO
************************************************************************



      iCASE=3
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 3

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IJGTL=IISTA,IIEND
          IJGTLTOT=IJGTL+NIGTJES(ISYM)
          IJABS=MIGTJ(1,IJGTLTOT)
          ILABS=MIGTJ(2,IJGTLTOT)
          IJ  =MIREL(1,IJABS)
          ISYJ=MIREL(2,IJABS)
          IL  =MIREL(1,ILABS)
          ISYL=MIREL(2,ILABS)
          DO ITGTU=IASTA,IAEND ! these are always all elements
            ITGTUTOT=ITGTU+NTGTUES(ISYM)
            ITABS=MTGTU(1,ITGTUTOT)
            IVABS=MTGTU(2,ITGTUTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integrals (tj,vl) and (tlvj)
            NV=NVTOT_CHOSYM(MUL(ISYT,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
            ITJ=IT-1+NASH(ISYT)*(IJ-1)
            IVL=IV-1+NASH(ISYV)*(IL-1)
            IOFFTJ=1+IOSYM(ISYT,ISYJ)+NV*ITJ
            IOFFVL=1+IOSYM(ISYV,ISYL)+NV*IVL
            TJVL=DDOT_(NV,CHOBUF(IOFFTJ),1,CHOBUF(IOFFVL),1)

            NV=NVTOT_CHOSYM(MUL(ISYT,ISYL))
            ITL=IT-1+NASH(ISYT)*(IL-1)
            IVJ=IV-1+NASH(ISYV)*(IJ-1)
            IOFFTL=1+IOSYM(ISYT,ISYL)+NV*ITL
            IOFFVJ=1+IOSYM(ISYV,ISYJ)+NV*IVJ
            TLVJ=DDOT_(NV,CHOBUF(IOFFTL),1,CHOBUF(IOFFVJ),1)

! BM(tv,jl)=((tj,vl)-(tl,vj))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(j,l))
            SCL=0.5D0
            !IF (ITABS.EQ.IVABS) SCL=SCL*0.5D0
            !IF (ILABS.EQ.IJABS) SCL=SCL*SQRTH
            BMTVJL=SCL*(TJVL-TLVJ)
! write element BM(tv,jl)
            IDX=ITGTU+NAS*(IJGTL-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=BMTVJL
#else
            GA_Arrays(lg_w)%Array(IDX)=BMTVJL
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 3      CONTINUE
      END DO
************************************************************************

      CALL mma_deallocate(CHOBUF)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_F(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      use caspt2_output, only:iPrGlb
      use PrintLevel, only: debug
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER IVEC

      INTEGER IOSYM(8,8)
      REAL*8, ALLOCATABLE:: CHOBUF(:)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case F'
      END IF

************************************************************************
* Case F (8,9):
C FP(tv,ac)=((at,cv)+(av,ct))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(a,c))
C FM(tv,ac)= -((at,cv)-(av,ct))/(2*SQRT(1+Kron(a,c))
************************************************************************

      SQRTH=SQRT(0.5D0)

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(3,NCHOBUF,IOSYM)

      CALL mma_allocate(CHOBUF,NCHOBUF,LABEL='CHOBUF')

      CALL CHOVEC_READ(3,CHOBUF,NCHOBUF)

      iCASE=8
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 8

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IAGEB=IISTA,IIEND
          IAGEBTOT=IAGEB+NAGEBES(ISYM)
          IAABS=MAGEB(1,IAGEBTOT)
          ICABS=MAGEB(2,IAGEBTOT)
          IA  =MAREL(1,IAABS)
          ISYA=MAREL(2,IAABS)
          IC  =MAREL(1,ICABS)
          ISYC=MAREL(2,ICABS)
          DO ITGEU=IASTA,IAEND ! these are always all elements
            ITGEUTOT=ITGEU+NTGEUES(ISYM)
            ITABS=MTGEU(1,ITGEUTOT)
            IVABS=MTGEU(2,ITGEUTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integrals (ta,vc) and (tc,va)
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYT)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAT=IA-1+NSSH(ISYA)*(IT-1)
            ICV=IC-1+NSSH(ISYC)*(IV-1)
            IOFFAT=1+IOSYM(ISYA,ISYT)+NV*IAT
            IOFFCV=1+IOSYM(ISYC,ISYV)+NV*ICV
            ATCV=DDOT_(NV,CHOBUF(IOFFAT),1,CHOBUF(IOFFCV),1)

            NV=NVTOT_CHOSYM(MUL(ISYA,ISYV)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAV=IA-1+NSSH(ISYA)*(IV-1)
            ICT=IC-1+NSSH(ISYC)*(IT-1)
            IOFFAV=1+IOSYM(ISYA,ISYV)+NV*IAV
            IOFFCT=1+IOSYM(ISYC,ISYT)+NV*ICT
            AVCT=DDOT_(NV,CHOBUF(IOFFAV),1,CHOBUF(IOFFCT),1)

! FP(tv,ac)=((at,cv)+(av,ct))*(1-Kron(t,v)/2)/(2*SQRT(1+Kron(a,c))
            SCL=0.5D0
            IF (ITABS.EQ.IVABS) SCL=SCL*0.5D0
            IF (IAABS.EQ.ICABS) SCL=SCL*SQRTH
            FPTVAC=SCL*(ATCV+AVCT)
! write element FP(tv,ac)
            IDX=ITGEU+NAS*(IAGEB-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=FPTVAC
#else
            GA_Arrays(lg_w)%Array(IDX)=FPTVAC
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 8      CONTINUE
      END DO
************************************************************************



      iCASE=9
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 9

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IAGTB=IISTA,IIEND
          IAGTBTOT=IAGTB+NAGTBES(ISYM)
          IAABS=MAGTB(1,IAGTBTOT)
          ICABS=MAGTB(2,IAGTBTOT)
          IA  =MAREL(1,IAABS)
          ISYA=MAREL(2,IAABS)
          IC  =MAREL(1,ICABS)
          ISYC=MAREL(2,ICABS)
          DO ITGTU=IASTA,IAEND ! these are always all elements
            ITGTUTOT=ITGTU+NTGTUES(ISYM)
            ITABS=MTGTU(1,ITGTUTOT)
            IVABS=MTGTU(2,ITGTUTOT)
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integrals (at,cv) and (av,ct)
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYT)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAT=IA-1+NSSH(ISYA)*(IT-1)
            ICV=IC-1+NSSH(ISYC)*(IV-1)
            IOFFAT=1+IOSYM(ISYA,ISYT)+NV*IAT
            IOFFCV=1+IOSYM(ISYC,ISYV)+NV*ICV
            ATCV=DDOT_(NV,CHOBUF(IOFFAT),1,CHOBUF(IOFFCV),1)

            NV=NVTOT_CHOSYM(MUL(ISYA,ISYV)) ! JSYM=ISYA*ISYA=ISYC*ISYC
            IAV=IA-1+NSSH(ISYA)*(IV-1)
            ICT=IC-1+NSSH(ISYC)*(IT-1)
            IOFFAV=1+IOSYM(ISYA,ISYV)+NV*IAV
            IOFFCT=1+IOSYM(ISYC,ISYT)+NV*ICT
            AVCT=DDOT_(NV,CHOBUF(IOFFAV),1,CHOBUF(IOFFCT),1)

! FM(tv,ac)= -((at,cv)-(av,ct))/(2*SQRT(1+Kron(a,c))
            SCL=0.5D0
            !IF (ITABS.EQ.IVABS) SCL=SCL*0.5D0
            !IF (IAABS.EQ.ICABS) SCL=SCL*SQRTH
            FMTVAC=SCL*(AVCT-ATCV)
! write element FM(tv,ac)
            IDX=ITGTU+NAS*(IAGTB-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=FMTVAC
#else
            GA_Arrays(lg_w)%Array(IDX)=FMTVAC
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 9      CONTINUE
      END DO
************************************************************************

      CALL mma_deallocate(CHOBUF)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_H(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      use caspt2_output, only:iPrGlb
      use PrintLevel, only: debug
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER IVEC

      INTEGER IOSYM(8,8)
      REAL*8, ALLOCATABLE:: CHOBUF(:)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case H'
      END IF

************************************************************************
* Case H:
C   WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
C   WM(jl,ac)=((ajcl)-(alcj))*SQRT(3.0D0)
************************************************************************

      SQRT3=SQRT(3.0D0)
      SQRTH=SQRT(0.5D0)

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(4,NCHOBUF,IOSYM)

      CALL mma_allocate(CHOBUF,NCHOBUF,LABEL='CHOBUF')

      CALL CHOVEC_READ(4,CHOBUF,NCHOBUF)

      iCASE=12
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NAGEB(ISYM)
        NIS=NIGEJ(ISYM)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 12

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IJGEL=IISTA,IIEND
          IJGELTOT=IJGEL+NIGEJES(ISYM)
          IJABS=MIGEJ(1,IJGELTOT)
          ILABS=MIGEJ(2,IJGELTOT)
          IJ  =MIREL(1,IJABS)
          ISYJ=MIREL(2,IJABS)
          IL  =MIREL(1,ILABS)
          ISYL=MIREL(2,ILABS)
          DO IAGEB=IASTA,IAEND ! these are always all elements
            IAGEBTOT=IAGEB+NAGEBES(ISYM)
            IAABS=MAGEB(1,IAGEBTOT)
            ICABS=MAGEB(2,IAGEBTOT)
            IA  =MAREL(1,IAABS)
            ISYA=MAREL(2,IAABS)
            IC  =MAREL(1,ICABS)
            ISYC=MAREL(2,ICABS)
! compute integrals (ajcl) and (alcj)
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
            IAJ=IA-1+NSSH(ISYA)*(IJ-1)
            ICL=IC-1+NSSH(ISYC)*(IL-1)
            IOFFAJ=1+IOSYM(ISYA,ISYJ)+NV*IAJ
            IOFFCL=1+IOSYM(ISYC,ISYL)+NV*ICL
            AJCL=DDOT_(NV,CHOBUF(IOFFAJ),1,CHOBUF(IOFFCL),1)

            NV=NVTOT_CHOSYM(MUL(ISYA,ISYL))
            IAL=IA-1+NSSH(ISYA)*(IL-1)
            ICJ=IC-1+NSSH(ISYC)*(IJ-1)
            IOFFAL=1+IOSYM(ISYA,ISYL)+NV*IAL
            IOFFCJ=1+IOSYM(ISYC,ISYJ)+NV*ICJ
            ALCJ=DDOT_(NV,CHOBUF(IOFFAL),1,CHOBUF(IOFFCJ),1)

! HP(ac,jl)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
            SCL=1.0D0
            IF (IAABS.EQ.ICABS) SCL=SCL*SQRTH
            IF (ILABS.EQ.IJABS) SCL=SCL*SQRTH
            HPACJL=SCL*(AJCL+ALCJ)
! write element HP(ac,jl)
            IDX=IAGEB+NAS*(IJGEL-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=HPACJL
#else
            GA_Arrays(lg_w)%Array(IDX)=HPACJL
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 12     CONTINUE
      END DO
************************************************************************



      iCASE=13
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NAGTB(ISYM)
        NIS=NIGTJ(ISYM)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 13

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IJGTL=IISTA,IIEND
          IJGTLTOT=IJGTL+NIGTJES(ISYM)
          IJABS=MIGTJ(1,IJGTLTOT)
          ILABS=MIGTJ(2,IJGTLTOT)
          IJ  =MIREL(1,IJABS)
          ISYJ=MIREL(2,IJABS)
          IL  =MIREL(1,ILABS)
          ISYL=MIREL(2,ILABS)
          DO IAGTB=IASTA,IAEND ! these are always all elements
            IAGTBTOT=IAGTB+NAGTBES(ISYM)
            IAABS=MAGTB(1,IAGTBTOT)
            ICABS=MAGTB(2,IAGTBTOT)
            IA  =MAREL(1,IAABS)
            ISYA=MAREL(2,IAABS)
            IC  =MAREL(1,ICABS)
            ISYC=MAREL(2,ICABS)
! compute integrals (ajcl) and (alcj)
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYJ)) ! JSYM=ISYA*ISYJ=ISYC*ISYL
            IAJ=IA-1+NSSH(ISYA)*(IJ-1)
            ICL=IC-1+NSSH(ISYC)*(IL-1)
            IOFFAJ=1+IOSYM(ISYA,ISYJ)+NV*IAJ
            IOFFCL=1+IOSYM(ISYC,ISYL)+NV*ICL
            AJCL=DDOT_(NV,CHOBUF(IOFFAJ),1,CHOBUF(IOFFCL),1)

            NV=NVTOT_CHOSYM(MUL(ISYA,ISYL))
            IAL=IA-1+NSSH(ISYA)*(IL-1)
            ICJ=IC-1+NSSH(ISYC)*(IJ-1)
            IOFFAL=1+IOSYM(ISYA,ISYL)+NV*IAL
            IOFFCJ=1+IOSYM(ISYC,ISYJ)+NV*ICJ
            ALCJ=DDOT_(NV,CHOBUF(IOFFAL),1,CHOBUF(IOFFCJ),1)

! HP(ac,jl)=((ajcl)-(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
            SCL=SQRT3
            HMACJL=SCL*(AJCL-ALCJ)
! write element HP(ac,jl)
            IDX=IAGTB+NAS*(IJGTL-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=HMACJL
#else
            GA_Arrays(lg_W)%Array(IDX)=HMACJL
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 13     CONTINUE
      END DO
************************************************************************

      CALL mma_deallocate(CHOBUF)

      RETURN
      END


*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_D(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      use caspt2_output, only:iPrGlb
      use caspt2_data, only: FIMO
      use PrintLevel, only: debug
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER IVEC

      INTEGER IOBRA1(8,8), IOKET1(8,8), IOBRA2(8,8), IOKET2(8,8)
      REAL*8, ALLOCATABLE:: BRABUF1(:), KETBUF1(:),
     &                      BRABUF2(:), KETBUF2(:)
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
      DIMENSION NFIMOES(8)

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case D'
      END IF

************************************************************************
* Case D (5,6):
C D1(tv,aj)=(aj,tv) + FIMO(a,j)*Kron(t,v)/NACTEL
C D2(tv,aj)=(tj,av)
************************************************************************

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(4,NBRABUF1,IOBRA1)
      CALL CHOVEC_SIZE(2,NKETBUF1,IOKET1)

      CALL mma_allocate(BRABUF1,NBRABUF1,LABEL='BRABUF1')
      CALL mma_allocate(KETBUF1,NKETBUF1,LABEL='KETBUF1')

      CALL CHOVEC_READ(4,BRABUF1,NBRABUF1)
      CALL CHOVEC_READ(2,KETBUF1,NKETBUF1)

      CALL CHOVEC_SIZE(3,NBRABUF2,IOBRA2)
      CALL CHOVEC_SIZE(1,NKETBUF2,IOKET2)

      CALL mma_allocate(BRABUF2,NBRABUF2,LABEL='BRABUF2')
      CALL mma_allocate(KETBUF2,NKETBUF2,LABEL='KETBUF2')

      CALL CHOVEC_READ(3,BRABUF2,NBRABUF2)
      CALL CHOVEC_READ(1,KETBUF2,NKETBUF2)

      iCASE=5
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      ! set up FIMO access
      ACTINV=1.0D0/DBLE(MAX(1,NACTEL))
      IFIMOES=0
      DO ISYM=1,NSYM
        NFIMOES(ISYM)=IFIMOES
        IFIMOES=IFIMOES+(NORB(ISYM)*(NORB(ISYM)+1))/2
      END DO

      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 8

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

! cases D1, D2 share the RHS along the tu superindex
        NAS1=NAS/2
        IASTA1=IASTA
        IAEND1=IAEND/2
        IASTA2=IAEND1+1
        IAEND2=IAEND

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
        DO IAJ=IISTA,IIEND
          IAJTOT=IAJ+NIAES(ISYM)
          IJABS=MIA(1,IAJTOT)
          IAABS=MIA(2,IAJTOT)
          IA  =MAREL(1,IAABS)
          ISYA=MAREL(2,IAABS)
          IJ  =MIREL(1,IJABS)
          ISYJ=MIREL(2,IJABS)
          DO ITV=IASTA1,IAEND1 ! these are always all elements
            ITABS=MTU(1,ITV+NTUES(ISYM))
            IVABS=MTU(2,ITV+NTUES(ISYM))
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integral (aj,tv)
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYJ))
            IOAJ=IA-1+NSSH(ISYA)*(IJ-1)
            IOTV=IT-1+NASH(ISYT)*(IV-1)
            IOFFAJ=1+IOBRA1(ISYA,ISYJ)+NV*IOAJ
            IOFFTV=1+IOKET1(ISYT,ISYV)+NV*IOTV
            AJTV=DDOT_(NV,BRABUF1(IOFFAJ),1,KETBUF1(IOFFTV),1)

! D1(tv,aj)=(aj,tv) + FIMO(a,j)*Kron(t,v)/NACTEL
! integrals only
            IDX=ITV+NAS*(IAJ-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=AJTV
#else
            GA_Arrays(lg_w)%Array(IDX)=AJTV
#endif
          END DO
! now, dress with FIMO(a,j), only if T==V, so ISYT==ISYV, so if ISYM==1
          IF (ISYM.EQ.1) THEN
            IATOT=IA+NISH(ISYA)+NASH(ISYA)
            FAJ=FIMO(NFIMOES(ISYA)+(IATOT*(IATOT-1))/2+IJ)
            ONEADD=FAJ*ACTINV
            DO IUABS=1,NASHT
              IUU=KTU(IUABS,IUABS)
              IDX=IUU+NAS*(IAJ-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=DBL_MB(MW+IDX-1)+ONEADD
#else
              GA_Arrays(lg_w)%Array(IDX)=GA_Arrays(lg_w)%Array(IDX)
     &                                  +ONEADD
#endif
            END DO
          END IF
          DO ITV=IASTA2,IAEND2 ! these are always all elements
            ITABS=MTU(1,ITV-NAS1+NTUES(ISYM))
            IVABS=MTU(2,ITV-NAS1+NTUES(ISYM))
            IT  =MTREL(1,ITABS)
            ISYT=MTREL(2,ITABS)
            IV  =MTREL(1,IVABS)
            ISYV=MTREL(2,IVABS)
! compute integral (av,tj)
            NV=NVTOT_CHOSYM(MUL(ISYA,ISYV))
            IOAV=IA-1+NSSH(ISYA)*(IV-1)
            IOTJ=IT-1+NASH(ISYT)*(IJ-1)
            IOFFAV=1+IOBRA2(ISYA,ISYV)+NV*IOAV
            IOFFTJ=1+IOKET2(ISYT,ISYJ)+NV*IOTJ
            AVTJ=DDOT_(NV,BRABUF2(IOFFAV),1,KETBUF2(IOFFTJ),1)

! D2(tv,aj)=(av,tj) + FIMO(a,j)*Kron(t,v)/NACTEL
            IDX=ITV+NAS*(IAJ-IISTA)
#ifdef _MOLCAS_MPP_
            DBL_MB(MW+IDX-1)=AVTJ
#else
            GA_Arrays(lg_w)%Array(IDX)=AVTJ
#endif
          END DO
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 8      CONTINUE

      END DO
************************************************************************

      CALL mma_deallocate(BRABUF1)
      CALL mma_deallocate(KETBUF1)

      CALL mma_deallocate(BRABUF2)
      CALL mma_deallocate(KETBUF2)

************************************************************************

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_E(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      use caspt2_output, only:iPrGlb
      use PrintLevel, only: debug
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER IVEC

      INTEGER IOBRA(8,8), IOKET(8,8)
      REAL*8, ALLOCATABLE:: BRABUF(:), KETBUF(:)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case E'
      END IF

************************************************************************
* Case E (6,7):
C EP(v,ajl)=((aj,vl)+(al,vj))/SQRT(2+2*Kron(j,l))
C EM(v,ajl)=((aj,vl)-(al,vj))*SQRT(3/2)
************************************************************************

* -SVC- Case E is slightly special, in that the inactive superindices are
* so large, that it is suboptimal to have a direct translation table for
* them. Instead, the code loops over symmetry blocks of A-JL and figures
* out if the indices on the processor fall within a block or not. Within
* a A-JL symmetry block, NA(ISYA) and NIGEJ(ISYJL) are known, so they can
* be determined by integer division. This could be optimized by combining
* it with loop peeling (on the todo list?).

      SQRTH=SQRT(0.5D0)
      SQRTA=SQRT(1.5D0)

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(4,NBRABUF,IOBRA)
      CALL CHOVEC_SIZE(1,NKETBUF,IOKET)

      CALL mma_allocate(BRABUF,NBRABUF,LABEL='BRABUF')
      CALL mma_allocate(KETBUF,NKETBUF,LABEL='KETBUF')

      CALL CHOVEC_READ(4,BRABUF,NBRABUF)
      CALL CHOVEC_READ(1,KETBUF,NKETBUF)

      iCASE=6
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 6

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYA=1,NSYM
          ISYJL=MUL(ISYA,ISYM)
          ISYV=ISYM

          NA=NSSH(ISYA)
          NJL=NIGEJ(ISYJL)
          ! what is start/end in this block?
          IAJGELSTA=MAX(IISTA-IOFF,1)
          IAJGELEND=MIN(IIEND-IOFF,NA*NJL)

          DO IAJGEL=IAJGELSTA,IAJGELEND
            IJGEL=(IAJGEL-1)/NA+1
            IA=IAJGEL-NA*(IJGEL-1)
            IJGELTOT=IJGEL+NIGEJES(ISYJL)
            IJABS=MIGEJ(1,IJGELTOT)
            ILABS=MIGEJ(2,IJGELTOT)
            IJ  =MIREL(1,IJABS)
            ISYJ=MIREL(2,IJABS)
            IL  =MIREL(1,ILABS)
            ISYL=MIREL(2,ILABS)
            DO IV=IASTA,IAEND ! these are always all elements
! compute integrals (ajvl) and (alvj)
              NV=NVTOT_CHOSYM(MUL(ISYA,ISYJ))
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IVL=IV-1+NASH(ISYV)*(IL-1)
              IOFFAJ=1+IOBRA(ISYA,ISYJ)+NV*IAJ
              IOFFVL=1+IOKET(ISYV,ISYL)+NV*IVL
              AJVL=DDOT_(NV,BRABUF(IOFFAJ),1,KETBUF(IOFFVL),1)

              NV=NVTOT_CHOSYM(MUL(ISYA,ISYL))
              IAL=IA-1+NSSH(ISYA)*(IL-1)
              IVJ=IV-1+NASH(ISYV)*(IJ-1)
              IOFFAL=1+IOBRA(ISYA,ISYL)+NV*IAL
              IOFFVJ=1+IOKET(ISYV,ISYJ)+NV*IVJ
              ALVJ=DDOT_(NV,BRABUF(IOFFAL),1,KETBUF(IOFFVJ),1)

! EP(v,ajl)=((aj,vl)+(al,vj))/SQRT(2+2*Kron(j,l))
              IF (ILABS.EQ.IJABS) THEN
                SCL=0.5D0
              ELSE
                SCL=SQRTH
              END IF
              EP=SCL*(AJVL+ALVJ)
! write element EP
              IDX=IV+NAS*(IAJGEL+IOFF-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=EP
#else
              GA_Arrays(lg_W)%Array(IDX)=EP
#endif
            END DO
          END DO

          IOFF=IOFF+NA*NJL
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 6      CONTINUE
      END DO
************************************************************************



      iCASE=7
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 7

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYA=1,NSYM
          ISYJL=MUL(ISYA,ISYM)
          ISYV=ISYM

          NA=NSSH(ISYA)
          NJL=NIGTJ(ISYJL)
          ! what is start/end in this block?
          IAJGTLSTA=MAX(IISTA-IOFF,1)
          IAJGTLEND=MIN(IIEND-IOFF,NA*NJL)

          DO IAJGTL=IAJGTLSTA,IAJGTLEND
            IJGTL=(IAJGTL-1)/NA+1
            IA=IAJGTL-NA*(IJGTL-1)
            IJGTLTOT=IJGTL+NIGTJES(ISYJL)
            IJABS=MIGTJ(1,IJGTLTOT)
            ILABS=MIGTJ(2,IJGTLTOT)
            IJ  =MIREL(1,IJABS)
            ISYJ=MIREL(2,IJABS)
            IL  =MIREL(1,ILABS)
            ISYL=MIREL(2,ILABS)
            DO IV=IASTA,IAEND ! these are always all elements
! compute integrals (ajvl) and (alvj)
              NV=NVTOT_CHOSYM(MUL(ISYA,ISYJ))
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IVL=IV-1+NASH(ISYV)*(IL-1)
              IOFFAJ=1+IOBRA(ISYA,ISYJ)+NV*IAJ
              IOFFVL=1+IOKET(ISYV,ISYL)+NV*IVL
              AJVL=DDOT_(NV,BRABUF(IOFFAJ),1,KETBUF(IOFFVL),1)

              NV=NVTOT_CHOSYM(MUL(ISYA,ISYL))
              IAL=IA-1+NSSH(ISYA)*(IL-1)
              IVJ=IV-1+NASH(ISYV)*(IJ-1)
              IOFFAL=1+IOBRA(ISYA,ISYL)+NV*IAL
              IOFFVJ=1+IOKET(ISYV,ISYJ)+NV*IVJ
              ALVJ=DDOT_(NV,BRABUF(IOFFAL),1,KETBUF(IOFFVJ),1)

! EM(v,ajl)=((aj,vl)-(al,vj))*SQRT(3/2)
              EM=SQRTA*(AJVL-ALVJ)
! write element EM
              IDX=IV+NAS*(IAJGTL+IOFF-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=EM
#else
              GA_Arrays(lg_W)%Array(IDX)=EM
#endif
            END DO
          END DO

          IOFF=IOFF+NA*NJL
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 7      CONTINUE
      END DO
************************************************************************

      CALL mma_deallocate(BRABUF)
      CALL mma_deallocate(KETBUF)

      RETURN
      END

*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD_G(IVEC)
      USE SUPERINDEX
      USE CHOVEC_IO
      use caspt2_output, only:iPrGlb
      use PrintLevel, only: debug
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
#ifndef _MOLCAS_MPP_
      use fake_GA, only: GA_Arrays
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER IVEC

      INTEGER IOBRA(8,8), IOKET(8,8)
      REAL*8, ALLOCATABLE:: BRABUF(:), KETBUF(:)
*      Logical Incore
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'RHS on demand: case G'
      END IF

************************************************************************
* Case G (10,11):
C GP(v,jac)=((av,cj)+(cv,aj))/SQRT(2+2*Kron(a,b))
C GM(v,jac)=((av,cj)-(cv,aj))*SQRT(3/2)
************************************************************************

* -SVC- Case G is slightly special, in that the inactive superindices are
* so large, that it is suboptimal to have a direct translation table for
* them. Instead, the code loops over symmetry blocks of J-AC and figures
* out if the indices on the processor fall within a block or not. Within
* a J-AC symmetry block, NJ(ISYJ) and NAGEB(ISYAC) are known, so they can
* be determined by integer division. This could be optimized by combining
* it with loop peeling (on the todo list?).

      SQRTH=SQRT(0.5D0)
      SQRTA=SQRT(1.5D0)

************************************************************************
CSVC: read in all the cholesky vectors (need all symmetries)
************************************************************************
      CALL CHOVEC_SIZE(3,NBRABUF,IOBRA)
      CALL CHOVEC_SIZE(4,NKETBUF,IOKET)

      CALL mma_allocate(BRABUF,NBRABUF,LABEL='BRABUF')
      CALL mma_allocate(KETBUF,NKETBUF,LABEL='KETBUF')

      CALL CHOVEC_READ(3,BRABUF,NBRABUF)
      CALL CHOVEC_READ(4,KETBUF,NKETBUF)

      iCASE=10
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 10

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYJ=1,NSYM
          ISYAC=MUL(ISYJ,ISYM)
          ISYV=ISYM

          NJ=NISH(ISYJ)
          NAC=NAGEB(ISYAC)
          ! what is start/end in this block?
          IJAGECSTA=MAX(IISTA-IOFF,1)
          IJAGECEND=MIN(IIEND-IOFF,NJ*NAC)

          DO IJAGEC=IJAGECSTA,IJAGECEND
            IAGEC=(IJAGEC-1)/NJ+1
            IJ=IJAGEC-NJ*(IAGEC-1)
            IAGECTOT=IAGEC+NAGEBES(ISYAC)
            IAABS=MAGEB(1,IAGECTOT)
            ICABS=MAGEB(2,IAGECTOT)
            IA  =MAREL(1,IAABS)
            ISYA=MAREL(2,IAABS)
            IC  =MAREL(1,ICABS)
            ISYC=MAREL(2,ICABS)
            DO IV=IASTA,IAEND ! these are always all elements
! compute integrals (ajvl) and (alvj)
              NV=NVTOT_CHOSYM(MUL(ISYA,ISYV))
              IAV=IA-1+NSSH(ISYA)*(IV-1)
              ICJ=IC-1+NSSH(ISYC)*(IJ-1)
              IOFFAV=1+IOBRA(ISYA,ISYV)+NV*IAV
              IOFFCJ=1+IOKET(ISYC,ISYJ)+NV*ICJ
              AVCJ=DDOT_(NV,BRABUF(IOFFAV),1,KETBUF(IOFFCJ),1)

              NV=NVTOT_CHOSYM(MUL(ISYC,ISYV))
              ICV=IC-1+NSSH(ISYC)*(IV-1)
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IOFFCV=1+IOBRA(ISYC,ISYV)+NV*ICV
              IOFFAJ=1+IOKET(ISYA,ISYJ)+NV*IAJ
              CVAJ=DDOT_(NV,BRABUF(IOFFCV),1,KETBUF(IOFFAJ),1)

C GP(v,jac)=((av,cj)+(cv,aj))/SQRT(2+2*Kron(a,b))
              IF (IAABS.EQ.ICABS) THEN
                SCL=0.5D0
              ELSE
                SCL=SQRTH
              END IF
              GP=SCL*(AVCJ+CVAJ)
! write element EP
              IDX=IV+NAS*(IJAGEC+IOFF-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=GP
#else
              GA_Arrays(lg_w)%Array(IDX)=GP
#endif
            END DO
          END DO

          IOFF=IOFF+NJ*NAC
        END DO
************************************************************************

        CALL RHS_RELEASE_UPDATE (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 10     CONTINUE
      END DO
************************************************************************



      iCASE=11
************************************************************************
* outer loop over symmetry blocks in the RHS
************************************************************************
      DO ISYM=1,NSYM

        NAS=NASUP(ISYM,ICASE)
        NIS=NISUP(ISYM,ICASE)
        NW=NAS*NIS

        IF(NW.EQ.0) GOTO 11

        CALL RHS_ALLO (NAS,NIS,lg_W)
        CALL RHS_ACCESS (NAS,NIS,lg_W,IASTA,IAEND,IISTA,IIEND,MW)
        NW=NAS*(IIEND-IISTA+1)

************************************************************************
* inner loop over RHS elements in symmetry ISYM
************************************************************************
! find start and end block
        IOFF=0
        DO ISYJ=1,NSYM
          ISYAC=MUL(ISYJ,ISYM)
          ISYV=ISYM

          NJ=NISH(ISYJ)
          NAC=NAGTB(ISYAC)
          ! what is start/end in this block?
          IJAGTCSTA=MAX(IISTA-IOFF,1)
          IJAGTCEND=MIN(IIEND-IOFF,NJ*NAC)

          DO IJAGTC=IJAGTCSTA,IJAGTCEND
            IAGTC=(IJAGTC-1)/NJ+1
            IJ=IJAGTC-NJ*(IAGTC-1)
            IAGTCTOT=IAGTC+NAGTBES(ISYAC)
            IAABS=MAGTB(1,IAGTCTOT)
            ICABS=MAGTB(2,IAGTCTOT)
            IA  =MAREL(1,IAABS)
            ISYA=MAREL(2,IAABS)
            IC  =MAREL(1,ICABS)
            ISYC=MAREL(2,ICABS)
            DO IV=IASTA,IAEND ! these are always all elements
! compute integrals (ajvl) and (alvj)
              NV=NVTOT_CHOSYM(MUL(ISYA,ISYV))
              IAV=IA-1+NSSH(ISYA)*(IV-1)
              ICJ=IC-1+NSSH(ISYC)*(IJ-1)
              IOFFAV=1+IOBRA(ISYA,ISYV)+NV*IAV
              IOFFCJ=1+IOKET(ISYC,ISYJ)+NV*ICJ
              AVCJ=DDOT_(NV,BRABUF(IOFFAV),1,KETBUF(IOFFCJ),1)

              NV=NVTOT_CHOSYM(MUL(ISYC,ISYV))
              ICV=IC-1+NSSH(ISYC)*(IV-1)
              IAJ=IA-1+NSSH(ISYA)*(IJ-1)
              IOFFCV=1+IOBRA(ISYC,ISYV)+NV*ICV
              IOFFAJ=1+IOKET(ISYA,ISYJ)+NV*IAJ
              CVAJ=DDOT_(NV,BRABUF(IOFFCV),1,KETBUF(IOFFAJ),1)

C GM(v,jac)=((av,cj)-(cv,aj))*SQRT(3/2)
              GM=SQRTA*(AVCJ-CVAJ)
! write element GM
              IDX=IV+NAS*(IJAGTC+IOFF-IISTA)
#ifdef _MOLCAS_MPP_
              DBL_MB(MW+IDX-1)=GM
#else
              GA_Arrays(lg_W)%Array(IDX)=GM
#endif
            END DO
          END DO

          IOFF=IOFF+NJ*NAC
        END DO
************************************************************************

        CALL RHS_Release_Update (lg_W,IASTA,IAEND,IISTA,IIEND)
        CALL RHS_SAVE (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
        CALL RHS_FREE (lg_W)
 11     CONTINUE
      END DO
************************************************************************

      CALL mma_deallocate(BRABUF)
      CALL mma_deallocate(KETBUF)

      RETURN
      END
