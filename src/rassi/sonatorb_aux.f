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
      SUBROUTINE SONATORB_PLOT (DENS, FILEBASE, CHARTYPE, ASS, BSS)
      use OneDat, only: sNoNuc, sNoOri
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
      Real*8 DENS(6,NBTRI)
      CHARACTER(LEN=*) FILEBASE
      CHARACTER(LEN=8) CHARTYPE
      INTEGER ASS,BSS

      CHARACTER(LEN=25) FNAME
      CHARACTER(LEN=16) KNUM
      CHARACTER(LEN=16) FNUM,XNUM
      CHARACTER(LEN=8) LABEL
      CHARACTER CDIR
      Real*8 Dummy(1)
      Integer iDummy(7,8)
      Real*8, allocatable:: SZZ(:), VEC(:), VEC2(:), DMAT(:), SCR(:)
      Real*8, allocatable:: VNAT(:), EIG(:), OCC(:)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLOTTING SECTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

      NBMX2=NBMX**2

c SZZ  - AO Overlap integral
c VEC  - AO Overlap eigenvectors
c EIG  - AO Overlap eigenvalues
c VEC2 - Eigenvectors of density matrix
c SCR  - Temporary for matrix multiplication
C NOTE: SCR COULD PROBABLY BE SOMETHING LIKE NBMX*(NBMX+1)/2
C       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
C       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
      CALL mma_allocate(SZZ,NBTRI,Label='SZZ')
      CALL mma_allocate(VEC,NBSQ,Label='VEC')
      CALL mma_allocate(VEC2,NBMX2,Label='VEC2')
      CALL mma_allocate(SCR,NBMX2,Label='SCR')
      CALL mma_allocate(EIG,NBST,Label='EIG')
      SZZ(:)=0.0D0
      VEC(:)=0.0D0
      VEC2(:)=0.0D0
      SCR(:)=0.0D0
      EIG(:)=0.0D0

      CALL mma_allocate(VNAT,NBSQ,Label='VNAT')
      VNAT(:)=0.0D0
      CALL mma_allocate(OCC,NBST,Label='OCC')
      OCC(:)=0.0D0

C READ ORBITAL OVERLAP MATRIX.
      IRC=-1

c IOPT=6, origin and nuclear contrib not read
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      LABEL='MLTPL  0'
      CALL RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF


C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=1
      LV=1
      LE=1
      VEC(:)=0.0D0
      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO I=1,NB**2,(NB+1)
          VEC(LV-1+I)=1.0D00
        END DO
        CALL JACOB(SZZ(LS),VEC,NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO I=1,NB
          EIG(LE1)=SZZ(LS1)
          X=1.0D00/SQRT(MAX(SZZ(LS1),1.0D-14))
          CALL DSCAL_(NB,X,VEC(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
        END DO
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
      END DO

      CALL mma_deallocate(SZZ)

      CALL mma_allocate(DMAT,NBMX2,Label='DMAT')

      IF(ITYPE.LE.2) THEN
        ISTART=3
        IEND=3
      ELSE
        ISTART=1
        IEND=3
      END IF

      DO IDIR=ISTART,IEND

        CDIR='?'
        IF(IDIR.EQ.1) CDIR='X'
        IF(IDIR.EQ.2) CDIR='Y'
        IF(IDIR.EQ.3) CDIR='Z'

        INV=1
        II2=0
        IOCC=0
        LV=1
        LE=1
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
          IF(NB.EQ.0) GOTO 1750

C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVALUES OF THE OVERLAP MATRIX:

C expand the triangular matrix for this symmetry to a square matrix
          DMAT(:)=0.0D0
          CALL DCOPY_(NBMX2,[0.0D00],0,SCR,1)
          DO J=1,NB
          DO I=1,J
            II2=II2+1
            IJ=NB*(J-1)+I
            JI=NB*(I-1)+J
            IF(I.NE.J) THEN
              DMAT(IJ)=DENS(IDIR,II2)/2.0d0
              DMAT(JI)=DENS(IDIR,II2)/2.0d0
            ELSE
              DMAT(IJ)=DENS(IDIR,II2)
              DMAT(JI)=DENS(IDIR,II2)
            END IF
          END DO
          END DO

          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 DMAT,NB,VEC(LV),NB,
     &                 0.0D0,SCR,NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 VEC(LV),NB,SCR,NB,
     &                 0.0D0,DMAT,NB)

          ID1=1
          ID2=1
          DO I=1,NB
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID1),NB)
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO


C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          SCR(:)=0.0D0
          ISCR=1
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
c simple averaging
              SCR(ISCR)=(DMAT(IJ)+DMAT(JI))/2.0d0

c add a factor of two to convert spin -> sigma
              IF(ITYPE.GE.3) SCR(ISCR)=SCR(ISCR)*2.0d0
              ISCR=ISCR+1
            END DO
          END DO

C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          CALL DCOPY_(NBMX2,[0.0D0],0,VEC2,1)
          CALL DCOPY_(NB,[1.0D0],0,VEC2,NB+1)

          CALL JACOB(SCR,VEC2,NB,NB)
          CALL JACORD(SCR,VEC2,NB,NB)

C JACORD ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=0
          DO I=1,NB
            II=II+I
            OCC(IOCC+NB+1-I)=SCR(II)
          END DO
          IOCC=IOCC+NB

C REEXPRESS THE EIGENVALUES IN AO BASIS FUNCTIONS. REVERSE ORDER.
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 VEC(LV),NB,VEC2,NB,
     &                 0.0D0,SCR,NB)
          I1=1
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,SCR(I1),1,VNAT(I2),1)
            I1=I1+NB
          END DO
          INV=INV+NB**2
          LV=LV+NB**2
          LE=LE+NB
1750      CONTINUE
        END DO

C WRITE OUT THIS SET OF NATURAL SPIN ORBITALS
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2)') ASS,".",BSS
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".",CDIR
       END IF
       WRITE(FNUM,'(I8)') BSS
       FNUM=ADJUSTL(FNUM)
       IF (ASS.NE.BSS) THEN
         WRITE(XNUM,'(I8,A)') ASS,'_'//TRIM(FNUM)
         FNUM=ADJUSTL(XNUM)
       END IF
       IF (ITYPE.GT.2) FNUM=CDIR//TRIM(FNUM)

       FNAME=FILEBASE//'.'//TRIM(FNUM)
       IF(ITYPE.EQ.1)
     &        WRITE(6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &             VNAT, OCC, Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

c       Test a few values
C        CALL ADD_INFO("SONATORB_PLOT", VNAT, 1, 4)

c    ONLYFOR NATURAL ORBITALS
      if(ITYPE.EQ.1)
     &       CALL ADD_INFO("SONATORB_NO_OCC", OCC, SUM(NBASF), 4)

      END DO

      CALL mma_deallocate(DMAT)
      CALL mma_deallocate(VEC)
      CALL mma_deallocate(VEC2)
      CALL mma_deallocate(SCR)
      CALL mma_deallocate(EIG)
      CALL mma_deallocate(VNAT)
      CALL mma_deallocate(OCC)

      END SUBROUTINE SONATORB_PLOT

      SUBROUTINE SONATORB_CPLOT (DENS, FILEBASE, CHARTYPE, ASS, BSS)
      use OneDat, only: sNoNuc, sNoOri, sOpSiz
      use rassi_aux, only: ipglob
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
      Real*8 DENS(6,NBTRI)
      CHARACTER(LEN=*) FILEBASE
      CHARACTER(LEN=8) CHARTYPE
      INTEGER ASS,BSS

      CHARACTER(LEN=25) FNAME
      CHARACTER(LEN=16) KNUM
      CHARACTER(LEN=16) FNUM,XNUM
      CHARACTER(LEN=8) LABEL
      CHARACTER CDIR
      Real*8 Dummy(1)
      Integer IDUM(1),iDummy(7,8)
      Real*8, Allocatable:: SZZ(:), VEC(:), VEC2(:), VEC2I(:), SCR(:)
      Real*8, Allocatable:: SCRI(:), EIG(:)
      Real*8, Allocatable:: VNAT(:), VNATI(:), OCC(:)
      Real*8, Allocatable:: DMAT(:), DMATI(:)
      Real*8, Allocatable:: SANG(:)
      Real*8, Allocatable:: SANGF(:), SANGTR(:), SANGTI(:)
      Real*8, Allocatable:: SANGTR2(:), SANGTI2(:)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLOTTING SECTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

      NBMX2=NBMX**2

c SZZ  - AO Overlap integral
c VEC  - AO Overlap eigenvectors
c EIG  - AO Overlap eigenvalues
c VEC2 - Eigenvectors of density matrix
c SCR  - Temporary for matrix multiplication
C NOTE: SCR COULD PROBABLY BE SOMETHING LIKE NBMX*(NBMX+1)/2
C       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
C       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
      CALL mma_allocate(SZZ,NBTRI,Label='SZZ')
      SZZ(:)=0.0D0
      CALL mma_allocate(VEC,NBSQ,Label='VEC')
      VEC(:)=0.0D0
      CALL mma_allocate(VEC2,NBMX2,Label='VEC2')
      VEC2(:)=0.0D0
      CALL mma_allocate(VEC2I,NBMX2,Label='VEC2I')
      VEC2I(:)=0.0D0
      CALL mma_allocate(SCR,NBMX2,Label='SCR')
      SCR(:)=0.0D0
      CALL mma_allocate(SCRI,NBMX2,Label='SCRI')
      SCRI(:)=0.0D0
      CALL mma_allocate(EIG,NBST,Label='EIG')
      EIG(:)=0.0D0

      CALL mma_allocate(VNAT,NBSQ,Label='VNAT')
      VNAT(:)=0.0D0
      CALL mma_allocate(VNATI,NBSQ,Label='VNATI')
      VNATI(:)=0.0D0
      CALL mma_allocate(OCC,NBST,Label='OCC')
      OCC(:)=0.0D0

C READ ORBITAL OVERLAP MATRIX.
      IRC=-1

c IOPT=6, origin and nuclear contrib not read
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      LABEL='MLTPL  0'
      CALL RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF


C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=1
      LV=1
      LE=1
      VEC(:)=0.0D0
      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO I=1,NB**2,(NB+1)
          VEC(LV-1+I)=1.0D00
        END DO
        CALL JACOB(SZZ(LS),VEC(LV),NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO I=1,NB
          EIG(LE1)=SZZ(LS1)
          X=1.0D00/SQRT(MAX(SZZ(LS1),1.0D-14))
          CALL DSCAL_(NB,X,VEC(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
        END DO
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
      END DO

      CALL mma_deallocate(SZZ)

      CALL mma_allocate(DMAT,NBMX2,Label='DMAT')
      CALL mma_allocate(DMATI,NBMX2,Label='DMATI')

      IF(ITYPE.LE.2) THEN
        ISTART=3
        IEND=3
      ELSE
        ISTART=1
        IEND=3
      END IF

      DO IDIR=ISTART,IEND

          CDIR='?'
          IF(IDIR.EQ.1) CDIR='X'
          IF(IDIR.EQ.2) CDIR='Y'
          IF(IDIR.EQ.3) CDIR='Z'


cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
C read in ao matrix for angmom or mltpl
      CALL mma_allocate(SANG,NBTRI,Label='SANG')
      SANG(:)=0.0D0

      IRC=-1
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      JOPT=ibset(0,sOpSiz)

      IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
        ICMP=1
        LABEL='MLTPL  0'
        CALL iRDONE(IRC,JOPT,LABEL,ICMP,IDUM,ISYLAB)
        CALL  RDONE(IRC,IOPT,LABEL,ICMP,SANG,ISYLAB)

        IF ( IRC.NE.0 ) THEN
          WRITE(6,*)
          WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
          WRITE(6,*)'      MLTPL0 INTEGRALS ARE NOT AVAILABLE'
          WRITE(6,*)'      IRC:',IRC
          WRITE(6,*)
          CALL ABEND()
        END IF

      ELSE IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
        ICMP=3
        LABEL='ANGMOM'
        CALL iRDONE(IRC,JOPT,LABEL,ICMP,IDUM,ISYLAB)
        CALL  RDONE(IRC,IOPT,LABEL,ICMP,SANG,ISYLAB)

        IF ( IRC.NE.0 ) THEN
          WRITE(6,*)
          WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
          WRITE(6,*)'      ANGMOM INTEGRALS ARE NOT AVAILABLE'
          WRITE(6,*)'      IRC:',IRC
          WRITE(6,*)
          CALL ABEND()
        END IF

      END IF

cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
        INV=1
        II2=0
        IOCC=0
        LV=1
        LE=1
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
          IF(NB.EQ.0) cycle

C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVALUES OF THE OVERLAP MATRIX:

C expand the triangular matrix for this symmetry to a square matrix
          DMAT(:)=0.0D0
          DMATI(:)=0.0D0
          SCR(:)=0.0D0
          SCRI(:)=0.0D0

          DO J=1,NB
          DO I=1,J
            II2=II2+1
            IJ=NB*(J-1)+I
            JI=NB*(I-1)+J
            IF(I.NE.J) THEN
              DMAT(IJ)=DENS(IDIR,II2)/2.0d0
              DMAT(JI)=DENS(IDIR,II2)/2.0d0
              DMATI(IJ)=-DENS(IDIR+3,II2)/2.0d0
              DMATI(JI)= DENS(IDIR+3,II2)/2.0d0
            ELSE
              DMAT(IJ)=DENS(IDIR,II2)
              DMATI(JI)=DENS(IDIR+3,II2)
            END IF
          END DO
          END DO

          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 DMAT,NB,VEC(LV),NB,
     &                 0.0D0,SCR,NB)
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 DMATI,NB,VEC(LV),NB,
     &                 0.0D0,SCRI,NB)



          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 VEC(LV),NB,SCR,NB,
     &                 0.0D0,DMAT,NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 VEC(LV),NB,SCRI,NB,
     &                 0.0D0,DMATI,NB)

          ID1=1
          ID2=1
          DO I=1,NB
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID1),NB)
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID2),1)
            CALL DSCAL_(NB,EIG(LE-1+I),DMATI(ID1),NB)
            CALL DSCAL_(NB,EIG(LE-1+I),DMATI(ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO


C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          SCR(:)=0.0D0
          SCRI(:)=0.0D0

          ISCR=1
          ISCRI=1
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
c simple averaging
              SCR(ISCR)=(DMAT(JI)+DMAT(IJ))/2.0d0
              SCRI(ISCRI)=(DMATI(JI)-DMATI(IJ))/2.0d0
c add a factor of two to convert spin -> sigma
              IF(ITYPE.GE.3) SCR(ISCR)=SCR(ISCR)*2.0d0
              IF(ITYPE.GE.3) SCRI(ISCRI)=SCRI(ISCRI)*2.0d0
              ISCR=ISCR+1
              ISCRI=ISCRI+1
            END DO
          END DO

C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          VEC2(:)=0.0D0
          VEC2I(:)=0.0D0

          CALL CPLOT_DIAG(SCR,SCRI, NB,VEC2,VEC2I)

C LAPACK ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=0
          DO I=1,NB
            II=II+I
            OCC(IOCC+NB+1-I)=SCR(II)
          END DO
          IOCC=IOCC+NB

C REEXPRESS THE EIGENVECTORS IN AO BASIS FUNCTIONS. REVERSE ORDER.
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 VEC(LV),NB,VEC2,NB,
     &                 0.0D0,SCR,NB)
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 VEC(LV),NB,VEC2I,NB,
     &                 0.0D0,SCRI,NB)

          I1=1
          I1I=1
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,SCR(I1),1,VNAT(I2),1)
            CALL DCOPY_(NB,SCRI(I1I),1,VNATI(I2),1)
            I1=I1+NB
            I1I=I1I+NB
          END DO
          INV=INV+NB**2
          LV=LV+NB**2
          LE=LE+NB
        END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC TESTING
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(IPGLOB.GE.4) THEN

      CALL mma_allocate(SANGF,NBMX**2,Label='SANGF')
      SANGF(:)=0.0D0
      CALL mma_allocate(SANGTR,NBMX**2,Label='SANGTR')
      CALL mma_allocate(SANGTI,NBMX**2,Label='SANGTI')
      SANGTR(:)=0.0D0
      SANGTI(:)=0.0D0
      CALL mma_allocate(SANGTR2,NBMX**2,Label='SANGTR2')
      CALL mma_allocate(SANGTI2,NBMX**2,Label='SANGTI2')
      SANGTR2(:)=0.0D0
      SANGTI2(:)=0.0D0

      INV=0
      INV2=0
      II=0
      SUM = 0.0d0
      SUMI = 0.0d0

      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        IF(NB.EQ.0) GOTO 1860

c       Expand integrals for this symmetry to full storage
        SANGF(:)=0.0D0

        DO J=1,NB
        DO I=1,J
          IJ=NB*(J-1)+I-1
          JI=NB*(I-1)+J-1

          SANGF(1+JI) = SANG(1+II)

          IF(I.NE.J) THEN
            IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
              SANGF(1+IJ) = -SANG(1+II)
            ELSE
              SANGF(1+IJ) =  SANG(1+II)
            END IF
          END IF

          II=II+1

        END DO
        END DO

        IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,SANGF,NB,
     &             VNAT(1+INV),NB,0.0d0,SANGTR,NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,SANGF,NB,
     &              VNATI(1+INV),NB,0.0d0,SANGTI,NB)

          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,VNAT(1+INV),NB,
     &             SANGTR,NB,0.0d0,SANGTR2,NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,VNATI(1+INV),NB,
     &             SANGTI,NB,1.0d0,SANGTR2,NB)

          CALL DGEMM_('T','N',NB,NB,NB,-1.0d0,VNATI(1+INV),NB,
     &             SANGTR,NB,0.0d0,SANGTI2,NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,VNAT(1+INV),NB,
     &             SANGTI,NB,1.0d0,SANGTI,NB)

        ELSE IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN

          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,SANGF,NB,
     &             VNAT(1+INV),NB,0.0d0,SANGTI,NB)
          CALL DGEMM_('T','N',NB,NB,NB,-1.0d0,SANGF,NB,
     &             VNATI(1+INV),NB,0.0d0,SANGTR,NB)

          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,VNAT(1+INV),NB,
     &             SANGTR,NB,0.0d0,SANGTR2,NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,VNATI(1+INV),NB,
     &             SANGTI,NB,1.0d0,SANGTR2,NB)

          CALL DGEMM_('T','N',NB,NB,NB,-1.0d0,VNATI(1+INV),NB,
     &             SANGTR,NB,0.0d0,SANGTI2,NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,VNAT(1+INV),NB,
     &             SANGTI,NB,1.0d0,SANGTI2,NB)

        END IF

c Sum over the trace
        DO I = 1,NB
          IJ = I+(I-1)*NB-1
          SUM  = SUM  + OCC(I+INV2) * SANGTR2(1+IJ)
          SUMI = SUMI + OCC(I+INV2) * SANGTI2(1+IJ)
        END DO

1860    CONTINUE

        INV=INV+NB**2
        INV2=INV2+NB

      END DO

        WRITE(6,*) "Ben P TEST for JA:"
        WRITE(6,*) "REAL: ",SUM
        WRITE(6,*) "IMAG: ",SUMI

        CALL mma_deallocate(SANGF)
        CALL mma_deallocate(SANGTR)
        CALL mma_deallocate(SANGTI)
        CALL mma_deallocate(SANGTR2)
        CALL mma_deallocate(SANGTI2)
      END IF ! IPGLOB >= 4

      CALL mma_deallocate(SANG)

C WRITE OUT THIS SET OF NATURAL SPIN ORBITALS
C REAL PART
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".","R"
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A,A,A)')ASS,".",BSS,".",CDIR,".","R"
       END IF
       WRITE(FNUM,'(I8)') BSS
       FNUM=ADJUSTL(FNUM)
       IF (ASS.NE.BSS) THEN
         WRITE(XNUM,'(I8,A)') ASS,'_'//TRIM(FNUM)
         FNUM=ADJUSTL(XNUM)
       END IF
       IF (ITYPE.GT.2) FNUM=CDIR//TRIM(FNUM)

       FNAME=FILEBASE//'.'//TRIM(FNUM)//'.R'
       IF(ITYPE.EQ.1)
     &        WRITE(6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &             VNAT, OCC, Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

C IMAGINARY PART
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".","I"
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A,A,A)')ASS,".",BSS,".",CDIR,".","I"
       END IF

       FNAME=FILEBASE//'.'//TRIM(FNUM)//'.I'
       IF(ITYPE.EQ.1)
     &        WRITE(6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &             VNATI, OCC, Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

c       Test a few values
C        CALL ADD_INFO("SONATORB_CPLOTR", VNAT, 1, 4)
C        CALL ADD_INFO("SONATORB_CPLOTI", VNATI, 1, 4)
C        CALL ADD_INFO("SONATORB_CPLOTO", OCC, 1, 4)

      END DO

      CALL mma_deallocate(DMAT)
      CALL mma_deallocate(DMATI)
      CALL mma_deallocate(VEC)
      CALL mma_deallocate(VEC2)
      CALL mma_deallocate(VEC2I)
      CALL mma_deallocate(SCR)
      CALL mma_deallocate(SCRI)
      CALL mma_deallocate(VNAT)
      CALL mma_deallocate(VNATI)
      CALL mma_deallocate(OCC)

      END SUBROUTINE SONATORB_CPLOT




      SUBROUTINE CPLOT_DIAG(MATR, MATI, DIM, EIGVECR, EIGVECI)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER DIM
      REAL*8 MATR(DIM*(DIM+1)/2),MATI(DIM*(DIM+1)/2)
      REAL*8 EIGVECR(DIM,DIM),EIGVECI(DIM,DIM)

      REAL*8 CEIGVAL(DIM)
      COMPLEX*16 MATFULL((DIM*(DIM+1)/2))
      COMPLEX*16 CEIGVEC(DIM,DIM)
      COMPLEX*16 ZWORK(2*DIM-1)
      REAL*8 RWORK(3*DIM-2)
      INTEGER INFO

      DO J=1,(DIM*(DIM+1)/2)
          MATFULL(J) = CMPLX(MATR(J),MATI(J),kind=8)
c          MATFULL(J) = CMPLX(MATR(J),0.0d0,kind=8)
      END DO


      call zhpev_('V','U',DIM,MATFULL,CEIGVAL,
     &           CEIGVEC,DIM,ZWORK,RWORK,INFO)


      IF(INFO.NE.0) THEN
          WRITE(6,*) "Error in diagonalization"
          WRITE(6,*) "INFO: ",INFO
          CALL ABEND()
      END IF

      DO I=1,DIM
      DO J=1,DIM
          EIGVECR(I,J) = REAL(CEIGVEC(I,J))
          EIGVECI(I,J) = AIMAG(CEIGVEC(I,J))
      END DO
      END DO

      CALL DCOPY_(DIM*(DIM+1)/2,[0.0D00],0,MATR,1)
      CALL DCOPY_(DIM*(DIM+1)/2,[0.0D00],0,MATI,1)

      DO J=1,DIM
         MATR((J*(J-1)/2)+J) = CEIGVAL(J)
      END DO

      END SUBROUTINE CPLOT_DIAG
