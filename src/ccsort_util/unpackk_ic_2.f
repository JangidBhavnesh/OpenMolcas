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
       subroutine unpackk_ic_2 (i,vint,ndimvi,ndimvj,Vic)
!
!     this routine vint(j,k,l) = <i,j|k,l>
!     for given i from incore (reduced) expanded block Vic
!     ie. symp=symr,symq=syms
!
!     i      - value of pivot index (I)
!     vint   - array of integrals (O)
!     ndimvi - (norb(symi)) (I)
!     ndimvj - (norb(symj)) (I)
!     Vic    - incore expanded block of integrals (I)
!
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimvi,ndimvj
       real*8 vint(1:ndimvj,1:ndimvi,1:ndimvj)
       real*8 Vic(1:(ndimvi*(ndimvi+1)/2),1:(ndimvj*(ndimvj+1)/2))
!
!     help variables
!
      integer j,k,l,ik,jl
!
!
        do k=1,ndimvi
        if (i.ge.k) then
          ik=i*(i-1)/2+k
        else
          ik=k*(k-1)/2+i
        end if
          jl=0
          do j=1,ndimvj
          do l=1,j
            jl=jl+1
            vint(j,k,l)=Vic(ik,jl)
            vint(l,k,j)=Vic(ik,jl)
          end do
          end do
        end do
!
       return
       end
