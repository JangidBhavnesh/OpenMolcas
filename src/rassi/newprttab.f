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
      Subroutine NEWPRTTAB(NSYM,NFRO,NISH,NRAS1,NRAS2,NRAS3,
     &                     NSSH,NDEL)
      use rassi_global_arrays, only: PART
      use stdalloc, only: mma_allocate
      IMPLICIT NONE
      INTEGER NSYM,NFRO,NISH,NRAS1,NRAS2,NRAS3,NSSH,NDEL
      INTEGER NSIZE,ITYPE,NPART,ISYM,IPART,ISUM
      DIMENSION NFRO(NSYM),NISH(NSYM),NRAS1(NSYM),NRAS2(NSYM)
      DIMENSION NRAS3(NSYM),NSSH(NSYM),NDEL(NSYM)

C Executable statements
      NPART=7
      NSIZE=4+(NSYM+1)*(NPART+1)
      ITYPE=93
      CALL mma_allocate(PART,NSIZE,Label='PART')
      PART(1)=NSIZE
      PART(2)=ITYPE
      PART(3)=NPART
      PART(4)=NSYM
      DO ISYM=1,NSYM
        PART(5+ISYM+(NSYM+1)*1)=NRAS1(ISYM)
        PART(5+ISYM+(NSYM+1)*2)=NRAS2(ISYM)
        PART(5+ISYM+(NSYM+1)*3)=NRAS3(ISYM)
        PART(5+ISYM+(NSYM+1)*4)=NISH(ISYM)
        PART(5+ISYM+(NSYM+1)*5)=NSSH(ISYM)
        PART(5+ISYM+(NSYM+1)*6)=NFRO(ISYM)
        PART(5+ISYM+(NSYM+1)*7)=NDEL(ISYM)
        ISUM=0
        IPART=1
10       CONTINUE
          ISUM=ISUM+PART(5+ISYM+(NSYM+1)*IPART)
         IPART=IPART+1
         IF(IPART.le.7) GO TO 10
c
c VV: original code was unrolled uncorrectly by GCC 3.0.4
c        DO IPART=1,7
c          ISUM=ISUM+PART(5+ISYM+(NSYM+1)*IPART)
c        END DO
c
        PART(5+ISYM+(NSYM+1)*0)=ISUM
      END DO
      DO IPART=0,7
        ISUM=0
        DO ISYM=1,NSYM
          ISUM=ISUM+PART(5+ISYM+(NSYM+1)*IPART)
        END DO
        PART(5+(NSYM+1)*IPART)=ISUM
      END DO

      END Subroutine NEWPRTTAB
