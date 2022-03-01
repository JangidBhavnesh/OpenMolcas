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
       subroutine esb_ic_3 (symp,Vic,dimp,pqind)
!
!     this routine realize expansion of symmetry block
!     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL),
!     (for case symp=symr=symq=syms)
!     provided such integrals exists
!     It found corresponding (IJ|KL) and expand it to
!     matrix vic (prqs)
!

#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
!
           integer symp,dimp
!LD           integer symp,dimp,dimq
               real*8 Vic(1:(dimp*(dimp+1)/2)*((dimp*(dimp+1)/2)+1)/2)
        real*8 val1
!

       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
       integer pqind(1:mbas,1:mbas)
!
!     help variables
       integer i,j,maxx,ik,jl,ikjl
       integer ind(1:4)
#include "tratoc.fh"
       integer INDMAX
       parameter (INDMAX=nTraBuf)
       REAL*8    TWO(INDMAX)
!
!I        calc pqind
!
!LD        if (dimp.ge.dimp) then
          maxx=dimp
!LD        else
!LD          maxx=dimq
!LD        end if
!
        do i=1,maxx
        do j=1,maxx
          if (i.ge.j) then
            pqind(i,j)=i*(i-1)/2+j
          else
            pqind(i,j)=j*(j-1)/2+i
          end if
        end do
        end do
!
!II    get  adress
       idis13=idis(symp,symp,symp)
!
!III.1define order of indices
!
       ni=np(symp,symp,symp)
       nj=nq(symp,symp,symp)
       nk=nr(symp,symp,symp)
       nl=ns(symp,symp,symp)
!
!
!     define NSI,NSJ,NSK,NSL
       ind(ni)=symp
       ind(nj)=symp
       ind(nk)=symp
       ind(nl)=symp
       NSI=ind(1)
       NSJ=ind(2)
       NSK=ind(3)
       NSL=ind(4)
!
       indtemp=indmax+1
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
         if (fullprint.ge.3) write (6,*) ' * K ind ',KOLD
!
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
         if (fullprint.ge.3) write (6,*) ' ** L ind ',LOLD
!
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
         if (fullprint.ge.3) write (6,*) ' *** I ind ',IOLD
!
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
         if (fullprint.ge.3) write (6,*) ' **** J ind ',JOLD
!
!
!*    read block of integrals if neccesarry
!
       if (indtemp.eq.(indmax+1)) then
       indtemp=1
!     read block
       CALL dDAFILE(LUINTM,2,TWO,INDMAX,IDIS13)
       end if
!
!*    write integrals to appropriate possitions
!
       val1=TWO(indtemp)
!
!:1   combination (ij|kl) -> (ij|kl)
!     since yes1 is always 1, if structure is skipped
       ind(1)=iold
       ind(2)=jold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
!
!:2    def iklj
       ik=pqind(i1,k1)
       jl=pqind(j1,l1)
       if (ik.ge.jl) then
         ikjl=ik*(ik-1)/2+jl
       else
         ikjl=jl*(jl-1)/2+ik
       end if
!
       Vic(ikjl)=val1
!
!
        indtemp=indtemp+1
!
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
!
!
       return
       end
