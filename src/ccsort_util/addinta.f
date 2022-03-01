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
       subroutine addinta (wrk,wrksize,                                 &
     & syma,ammap)
!
!     this routine do for all a in syma
!     1- reconstruct #2 <_a,m,p,q> from TEMPDA2 file
!     2- prepair corresponding <_am p q> (like <amef>aaaa) to #3
!     and write it to opened INTA1-4
!     N.B.  this routine use followuing foreign routines:
!     wrtmap
!     wri
!
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer syma
       integer ammap(1:mbas,1:8,1:8)
!
!     help variables
!
       integer lenefaaaa,lenefbaab,lenefbbbb,lenefabab
       integer lenejaaaa,lenejbaab,lenejbaba,lenejbbbb,lenejabab,       &
     & lenejabba
       integer posst,rc,a
!
!*    mapd2 and mapi2 of #2 <_a,m|p,q> are prepaired
!
!*    make required mapd3 and mapi3 and write them to INTA1-4
!     define lengths of this mediates
!
!*1   to INTA1 <m,_a||ef>aaaa, <m,_a||ef>baab
       call ccsort_grc0(3,2,1,3,3,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenefaaaa)
       call dawrtmap (luna1,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,2,3,4,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenefbaab)
       call dawrtmap (luna1,mapd3,mapi3,rc)
!
!*2   to INTA2 <m,_a||ef>bbbb, <m,_a||ef>abab
       call ccsort_grc0(3,2,2,4,4,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenefbbbb)
       call dawrtmap (luna2,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,1,3,4,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenefabab)
       call dawrtmap (luna2,mapd3,mapi3,rc)
!
!*3   to INTA3 <m,_a||ej>aaaa, <m,_a||ej>baab, <m,_a||ej>baba
       call ccsort_grc0(3,0,1,3,1,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejaaaa)
       call dawrtmap (luna3,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,2,3,2,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejbaab)
       call dawrtmap (luna3,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,2,4,1,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejbaba)
       call dawrtmap (luna3,mapd3,mapi3,rc)
!
!*4   to INTA4 <m,_a||ej>bbbb, <m,_a||ej>abba, <m,_a||ej>abab
       call ccsort_grc0(3,0,2,4,2,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejbbbb)
       call dawrtmap (luna4,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,1,4,1,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejabba)
       call dawrtmap (luna4,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,1,3,2,0,syma,poss30,posst,mapd3,mapi3)
       call deflength (mapd3,lenejabab)
       call dawrtmap (luna4,mapd3,mapi3,rc)
!
!
!*    cycle over a
!
       do 1000 a=1,nvb(syma)
!
!*    reconstruct #2 <_a,m,p,q> for given _a
       call mkampq (wrk,wrksize,                                        &
     & a,ammap)
!
!*    get contributions to INTA2 <m,_a||ef>bbbb, <m,_a||ef>abab
!     and wtite it there
!
       if (lenefbbbb.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,2,2,4,4,1,1)
       call dawri (luna2,lenefbbbb,wrk(mapd3(1,1)))
       end if
!
       if (lenefabab.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,0,1,3,4,1,0)
       call dawri (luna2,lenefabab,wrk(mapd3(1,1)))
       end if
!
!*    get contributions to INTA4 <m,_a||ej>bbbb, <m,_a||ej>abba, <m,_a||ej>abab
!     and wtite it there
!
       if (lenejbbbb.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,0,2,4,2,1,1)
       call dawri (luna4,lenejbbbb,wrk(mapd3(1,1)))
       end if
!
       if (lenejabba.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,0,1,4,1,0,1)
       call dawri (luna4,lenejabba,wrk(mapd3(1,1)))
       end if
!
       if (lenejabab.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,0,1,3,2,1,0)
       call dawri (luna4,lenejabab,wrk(mapd3(1,1)))
       end if
!
       if (a.gt.(nvb(syma)-nva(syma))) then
!     contributions to INTA1 and INTA3 only for a-alfa
!
!*    get contributions to INTA1 <m,_a||ef>aaaa, <m,_a||ef>baab if any
!     and wtite it there
!
       if (lenefaaaa.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,2,1,3,3,1,1)
       call dawri (luna1,lenefaaaa,wrk(mapd3(1,1)))
       end if
!
       if (lenefbaab.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,0,2,3,4,0,1)
       call dawri (luna1,lenefbaab,wrk(mapd3(1,1)))
       end if
!
!*    get contributions to INTA3 <m,_a||ej>aaaa, <m,_a||ej>baab, <m,_a||ej>baba
!     and wtite it there
!
       if (lenejaaaa.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,0,1,3,1,1,1)
       call dawri (luna3,lenejaaaa,wrk(mapd3(1,1)))
       end if
!
       if (lenejbaab.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,0,2,3,2,0,1)
       call dawri (luna3,lenejbaab,wrk(mapd3(1,1)))
       end if
!
       if (lenejbaba.gt.0) then
       call expmpq (wrk,wrksize,                                        &
     & syma,0,2,4,1,1,0)
       call dawri (luna3,lenejbaba,wrk(mapd3(1,1)))
       end if
!
       end if
!
 1000   continue
!
       return
       end
