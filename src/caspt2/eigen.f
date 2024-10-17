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
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************
* This subroutine diagonalizes a real symmetric matrix A using
* the Jacobi algorithm
**
      subroutine eigen(A,U,N)
      use stdalloc, only: mma_allocate, mma_deallocate
      implicit None
#include "WrkSpc.fh"

      integer N
      real(8) A(N,N)
      real(8) U(N,N)

      integer:: NSCR, IJ, I, J
      real(8), allocatable:: SCR(:)



      NSCR=(N*(N+1))/2
      call mma_allocate(SCR,NSCR,LABEL='SCR')

      IJ=0
      do I=1,N
        do J=1,I
          IJ=IJ+1
          SCR(IJ)=A(I,J)
        end do
      end do

* Initialize U as the identity matrix
      U=0.0d0
      call dcopy_(N,[1.0D0],0,U,N+1)

* Call Jacobi algorithm
      call JACOB(SCR,U,N,N)

      call mma_deallocate(SCR)

      end subroutine eigen
