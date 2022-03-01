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
       SUBROUTINE mkadress (NOIPSB)

!FUE
!FUE   It is strictly forbidden to use any construct to fetch disk
!FUE   addresses which do not make use of subroutines provided in
!FUE   the MOLCAS utilities as otherwise transparency and portability
!FUE   is lost.
!FUE

!
!     this routine prepair information arrays
!     for integral blocks per symmetry in primitive form
!     it defines:
!     noipsb(nijkl) - number of integrals per symmetry block
!     idispsb(nijkl)- initial addreses for exach symmetry block
!     typ(pa,qa,ra) - typ of symmetry block
!     types of (ij|kl):
!     1 - si=sk, si=sj, sk=sl
!     2 - si=sk, si=sj, sk>sl
!     3 - si=sk, si>sj, sk=sl
!     4 - si=sk, si>sj, sk>sl
!     5 - si>sk, si=sj, sk=sl
!     6 - si>sk, si=sj, sk>sl
!     7 - si>sk, si>sj, sk=sl
!     8 - si>sk, si>sj, sk>sl
!
!     idis(pa,qa,ra)- initial addreses for given symmetry block
!     np(pa,qa,ra)  - possition of p index in original block (ij|kl) (on tape)
!     nq(pa,qa,ra)  - possition of q index in original block (ij|kl) (on tape)
!     nr(pa,qa,ra)  - possition of r index in original block (ij|kl) (on tape)
!     ns(pa,qa,ra)  - possition of s index in original block (ij|kl) (on tape)
!
!     N.B. typ,idis,np,nq,nr,ns are trensfered thhrough common block /edpand2/
!
!
       IMPLICIT REAL*8(A-H,O-Z)

#include "tratoc.fh"
       integer INDMAX,MXFUNC
       PARAMETER (INDMAX=nTraBuf,MXFUNC=200)


#include "SysDef.fh"
#include "ccsort.fh"
#include "reorg.fh"

       integer NOIPSB(106)
!      integer idispsb(106)
!
!     help variables
!
       integer sense
       integer p,q,r,s,pa,qa,ra
       integer IND,INDT,ISPQRS,NINT,NSLM,idistemp,idishelp
       integer jlow,ilow,iup,jup,kup,lup,iold,jold,kold,lold
       integer norbp,nsi,nsj,nsk,nsl,nsij,nsijk
       real*8 dum(1)

!FUE   - pick the start addresses of each symmetry allowed integral
!FUE     block from the tranformed two electron integral file
        idistemp=0
        Call iDaFile(LUINTM,0,iTraToc,nTraToc,idistemp)
!FUE   - the option 0 in the call to dafile does not make any I/O
!FUE     but just returns the updated disk address
!FUE   - at this point idistemp points to the first record

!
       IND=0
       INDT=0
!FUE   idistemp=1
!FUE   idisadd=150
!
       do 100 NSI=1,NSYM
       do 101 NSJ=1,NSYM
       do 102 NSK=1,NSYM
       typ(NSI,NSJ,NSK)=0
 102    continue
 101    continue
 100    continue
!

       if (fullprint.gt.0) then
       Write(6,'(6X,A)') 'Transformed integral blocks:'
       Write(6,'(6X,A)') '----------------------------'
       Write(6,*)
       Write(6,'(6X,A)')                                                &
     & 'block  symmetry      no. of        no. of '
       Write(6,'(6X,A)')                                                &
     & ' no.    spec.        orbitals     integrals'
       Write(6,'(6X,A)')                                                &
     & '-------------------------------------------'
       end if
!
       ISPQRS=0
       DO 300 NSI=1,NSYM
       DO 301 NSJ=1,NSI
       NSIJ=MUL(NSI,NSJ)
       DO 302 NSK=1,NSI
       NSIJK=MUL(NSK,NSIJ)
       NSLM=NSK
       IF(NSK.EQ.NSI) NSLM=NSJ
       DO 303 NSL=1,NSLM
       IF(NSIJK.NE.NSL) GO TO 303
       NORBP=NORB(NSI)*NORB(NSJ)*NORB(NSK)*NORB(NSL)
       IF(NORBP.EQ.0)GO TO 303
       ISPQRS=ISPQRS+1
!
!     def
!
!     redefine indices from Parr to Dirac notation <p,q|r,s> <-> (i,j|k,l)
!
       p=nsi
       q=nsk
       r=nsj
       s=nsl
!
!     def. sense
!
!     type (ij|kl)
!     1 - si=sk, si=sj, sk=sl
!     2 - si=sk, si=sj, sk>sl
!     3 - si=sk, si>sj, sk=sl
!     4 - si=sk, si>sj, sk>sl
!     5 - si>sk, si=sj, sk=sl
!     6 - si>sk, si=sj, sk>sl
!     7 - si>sk, si>sj, sk=sl
!     8 - si>sk, si>sj, sk>sl
!
       if (nsi.eq.nsk) then
       if (nsi.eq.nsj) then
       if (nsk.eq.nsl) then
       sense=1
       else
       sense=2
       end if
       else
       if (nsk.eq.nsl) then
       sense=3
       else
       sense=4
       end if
       end if
       else
       if (nsi.eq.nsj) then
       if (nsk.eq.nsl) then
       sense=5
       else
       sense=6
       end if
       else
       if (nsk.eq.nsl) then
       sense=7
       else
       sense=8
       end if
       end if
       end if
!
       if (nsijk.ne.nsl) then
       sense=0
       else if (NORB(NSI)*NORB(NSJ)*NORB(NSK)*NORB(NSL).eq.0) then
       sense=0
       end if
!
!
!1:   perm <pq|rs> -> <pq|rs>
!
       pa=p
       qa=q
       ra=r
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=1
       nq(pa,qa,ra)=3
       nr(pa,qa,ra)=2
       ns(pa,qa,ra)=4
!
!2:   perm <pq|rs> -> <rq|ps> 1-3
!
       pa=r
       qa=q
       ra=p
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=2
       nq(pa,qa,ra)=3
       nr(pa,qa,ra)=1
       ns(pa,qa,ra)=4
!
!3:   perm <pq|rs> -> <ps|rq> 2-4
!
       pa=p
       qa=s
       ra=r
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=1
       nq(pa,qa,ra)=4
       nr(pa,qa,ra)=2
       ns(pa,qa,ra)=3
!
!4:   perm <pq|rs> -> <pq|rs> 1-3,2-4
!
       pa=r
       qa=s
       ra=p
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=2
       nq(pa,qa,ra)=4
       nr(pa,qa,ra)=1
       ns(pa,qa,ra)=3
!
!5:   perm <pq|rs> -> <qp|sr> 1-2,3-4
!
       pa=q
       qa=p
       ra=s
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=3
       nq(pa,qa,ra)=1
       nr(pa,qa,ra)=4
       ns(pa,qa,ra)=2
!
!6:   perm <pq|rs> -> <qp|sr> 1-2,3-4 -> <sp|qr> 1-3
!
       pa=s
       qa=p
       ra=q
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=4
       nq(pa,qa,ra)=1
       nr(pa,qa,ra)=3
       ns(pa,qa,ra)=2
!
!7:   perm <pq|rs> -> <qp|sr> 1-2,3-4 -> <qr|sp> 2-4
!
       pa=q
       qa=r
       ra=s
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=3
       nq(pa,qa,ra)=2
       nr(pa,qa,ra)=4
       ns(pa,qa,ra)=1
!
!8:   perm <pq|rs> -> <qp|sr> 1-2,3-4 -> <sr|qp> 1-3,2-4
!
       pa=s
       qa=r
       ra=q
       typ(pa,qa,ra)=sense
       idis(pa,qa,ra)=idistemp
       np(pa,qa,ra)=4
       nq(pa,qa,ra)=2
       nr(pa,qa,ra)=3
       ns(pa,qa,ra)=1
!
!
!      idispsb(ispqrs)=idistemp
       idishelp=0
!
!     ******************************************************************
!
!     LOOP OVER THE USED ORBITALS OF EACH SYMMETRY BLOCK
!     THIS LOOPING IS COPIED FROM THE MANUAL OF THE 4-INDEX
!     TRANSFORMATION PROGRAM
!
!     NINT COUNTS INTEGRAL LABELS IN THE GIVEN SYMMETRY BLOCK
!
       NINT=0
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
!
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
!
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
!
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
!
       IND=IND+1
       INDT=INDT+1
       NINT=NINT+1
       idishelp=idishelp+1
       if (idishelp.gt.INDMAX) then
!FUE     - all integrals in the reord were processed, hence
!FUE       update the disk address by a dummy I/O
         dum(1)=0.0d0
         Call dDaFile(LUINTM,0,dum,INDMAX,idistemp)
!FUE     idistemp=idistemp+idisadd
         idishelp=1
       end if
!
       IF (IND.LT.INDMAX) GO TO 404
       IND=0
!
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
       if (idishelp.gt.0) then
!FUE     - all integrals in the reord were processed, hence
!FUE       update the disk address by a dummy I/O
         dum(1)=0.0d0
         Call dDaFile(LUINTM,0,dum,INDMAX,idistemp)
!FUE   idistemp=idistemp+idisadd
       end if
!
!     WRITING THE LAST RECORD OF LABELS IN THE GIVEN SYMMETRY BLOCK
!     RECORDS ON LUPACK ARE FORMATTED TO 28KB LENGTH
!
       IF(IND.NE.0) THEN
       IND=0
       ENDIF
!
       NOIPSB(ISPQRS)=NINT
!
!     ******************************************************************
!
       if (fullprint.gt.0) then
       Write(6,'(6X,I5,2X,4I2,2X,4I4,2X,I8)')                           &
     &       ISPQRS,                                                    &
     &       NSI,NSJ,NSK,NSL,                                           &
     &       IUP,JUP,KUP,LUP,                                           &
     &       NINT
       end if
!
!     ******************************************************************
!
!
 303    CONTINUE
 302    CONTINUE
 301    CONTINUE
 300    CONTINUE
       if (fullprint.gt.0) then
       Write(6,'(6X,A)')                                                &
     & '-------------------------------------------'
       end if
!
       RETURN
       END
!
