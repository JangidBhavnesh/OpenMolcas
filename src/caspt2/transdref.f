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
      SUBROUTINE TRANSDREF(TORB,NTORB,DREF,NDREF)
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER NDREF
      REAL*8 TORB(NTORB),DREF(NDREF)

      REAL*8, ALLOCATABLE:: DSQ(:), TSQ(:), TMP(:)
* Purpose: given an orbital transformation array
* transform the DREF array (blocked triangular, active)

      NT=0
      NAMX=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        NO=NI+NR1+NR2+NR3+NS
        NA=NR1+NR2+NR3
        NAMX=MAX(NAMX,NA)
        NT=NT+NI**2+NR1**2+NR2**2+NR3**2+NS**2
      END DO
      CALL mma_allocate(DSQ,NAMX**2,LABEL='DSQ')
      CALL mma_allocate(TSQ,NAMX**2,LABEL='TSQ')
      CALL mma_allocate(TMP,NAMX**2,LABEL='TMP')
      IDOFF=0
      ITOFF=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NA=NR1+NR2+NR3
        NS=NSSH(ISYM)
        NO=NI+NA+NS
        IF (NO.eq.0) GOTO 99
* Copy the matrices to square storage: first fill with zeroes.
        TSQ(1:NA**2)=0.0D0
        IOFF=0
        ITOFF=ITOFF+NI**2
        DO I=1,NR1
         II=IOFF+I
         DO J=1,NR1
          JJ=IOFF+J
          TSQ(II+NA*(JJ-1))=TORB(ITOFF+I+NR1*(J-1))
         END DO
        END DO
*---
        ITOFF=ITOFF+NR1**2
        IOFF=IOFF+NR1
        DO I=1,NR2
         II=IOFF+I
         DO J=1,NR2
          JJ=IOFF+J
          TSQ(II+NA*(JJ-1))=TORB(ITOFF+I+NR2*(J-1))
         END DO
        END DO
*---
        ITOFF=ITOFF+NR2**2
        IOFF=IOFF+NR2
        DO I=1,NR3
         II=IOFF+I
         DO J=1,NR3
          JJ=IOFF+J
          TSQ(II+NA*(JJ-1))=TORB(ITOFF+I+NR3*(J-1))
         END DO
        END DO
*--- Finally, the secondary orbitals (non-deleted, virtual).
        ITOFF=ITOFF+NR3**2
        ITOFF=ITOFF+NS**2
* Now transfer the DREF matrix block to square storage:
        IJ=0
        DO I=1,NA
         DO J=1,I
          IJ=IJ+1
          DSQ(J+NA*(I-1))=DREF(IDOFF+IJ)
          DSQ(I+NA*(J-1))=DREF(IDOFF+IJ)
         END DO
        END DO
* Transform, first do DSQ*TSQ -> TMP...
       CALL DGEMM_('N','N',NA,NA,NA,1.0D0,DSQ,NA,TSQ,NA,
     &              0.0D0,TMP,NA)
* ... and then do TSQ(transpose)*TMP -> DSQ...
       CALL DGEMM_('T','N',NA,NA,NA,1.0D0,TSQ,NA,TMP,NA,
     &              0.0D0,DSQ,NA)
* Transfer DSQ values back to D, in triangular storage.
       IJ=0
       DO I=1,NA
        DO J=1,I
         IJ=IJ+1
         DREF(IDOFF+IJ)=DSQ(I+NA*(J-1))
        END DO
       END DO
       IDOFF=IDOFF+(NA*(NA+1))/2
* and repeat, using next symmetry block.
  99   CONTINUE
      END DO
      CALL mma_deallocate(DSQ)
      CALL mma_deallocate(TSQ)
      CALL mma_deallocate(TMP)


      RETURN
      END
