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
Module Local_Arrays
use stdalloc, only: mma_allocate, mma_deallocate

Private
Integer, Allocatable, Public:: CLBT(:), CLEBT(:), CI1BT(:), CIBT(:), CBLTP(:)

Public :: Allocate_Local_Arrays, Deallocate_Local_Arrays
Contains

Subroutine Allocate_Local_Arrays(MXNTTS,NSMST)
Integer :: MXNTTS,NSMST
CALL mma_allocate(CLBT ,MXNTTS,Label='CLBT')
CALL mma_allocate(CLEBT,MXNTTS,Label='CLEBT')
CALL mma_allocate(CI1BT,MXNTTS,Label='CI1BT')
CALL mma_allocate(CIBT ,8*MXNTTS,Label='CIBT')
CALL mma_allocate(CBLTP,NSMST,Label='CBLTP')
End Subroutine Allocate_Local_Arrays

Subroutine Deallocate_Local_Arrays()
If(Allocated(CLBT))      CALL mma_deallocate(CLBT)
If(Allocated(CLEBT))     CALL mma_deallocate(CLEBT)
If(Allocated(CI1BT))      CALL mma_deallocate(CI1BT)
If(Allocated(CIBT))      CALL mma_deallocate(CIBT)
If(Allocated(CBLTP))      CALL mma_deallocate(CBLTP)
End Subroutine Deallocate_Local_Arrays
End Module Local_Arrays
