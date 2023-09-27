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
      Module k2_setup
      Integer, Parameter:: nDArray=11,nDScalar=9

      Real*8, Allocatable :: Data_k2(:)
      Logical :: k2_processed=.False.
      Integer, Allocatable :: IndK2(:,:)
      Integer nk2, nIndk2
      End Module k2_setup
