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
Module rassi_data
INTEGER     NFROT,NISHT,NASHT,NOSHT,NSSHT,NDELT,NBST,             &
            NFRO(8),NISH(8),NASH(8),NOSH(8),NSSH(8),NDEL(8),      &
            NBASF(8),                                             &
            NFES(8),NIES(8),NAES(8),NOES(8),NSES(8),NDES(8),      &
            NBES(8)
CHARACTER(LEN=8)    WFTYPE
! NISHT  - NR OF  INACTIVE ORBITALS, TOTAL
! NASHT  - NR OF    ACTIVE ORBITALS, TOTAL
! NOSHT  - NR OF  OCCUPIED ORBITALS, TOTAL
! NSSHT  - NR OF SECONDARY ORBITALS, TOTAL
! NISH(8),...,NSSH(8), AS ABOVE, BUT BY SYMMETRY TYPE.
! NBST, NBAS(8) - SIMILAR, NR OF BASIS FUNCTIONS.
REAL*8 ENUC
INTEGER      NBMX,NBTRI,NBSQ,NBSQPR(8),                           &
             MXCI,NCMO,NTDMAB,NTDMZZ,NTDMS,NTDMA,                 &
             NSXY,NTRA,NCXA,LNILPT,LINILPT
! NBMX   - MAX NR OF BASIS FUNCTIONS OF ANY SPECIFIC SYMMETRY.
! NBTRI  - TOTAL SIZE OF TRIANGULAR SYMMETRY BLOCKS OF BASIS FNCS.
! NBSQ   - D:O, SQUARE SYMMETRY BLOCKS.
! NBSQPR - ACCUMULATED NR OF SQUARE SYMMETRY BLOCKS OF PREVIOUS
!          SYMMETRY TYPES. NBSQPR(1)=0.
! MXCI   - LARGEST NEEDED CI ARRAY FOR A STATE OF ANY SYMMETRY.
! NCMO   - SIZE OF CMO COEFFICIENT ARRAYS, = SUM(NOSH(I)*NBAS(I)).
! NTDMAB - SIZE OF TRANS.D. MATRIX IN BIORTHONORMAL MO BASIS.
! NTDMZZ - SIZE OF TRANS.D. MATRIX IN AO BASIS.
! NSXY   - SIZE OF MO OVERLAP ARRAY.
! NTRA   - SIZE OF TRANSFORMATION COEFFICIENT ARRAY.
! NCXA   - SIZE OF TRANSFORMATION MATRIX.
! LNILPT - WORK(LNILPT) IS A VALID DUMMY FIELD
!! LINILPT- IWORK(LINILPT) IS A VALID DUMMY FIELD
Real*8 ChFracMem
! ChFracMem - fraction of memory for the Cholesky vectors buffer
End Module rassi_data
