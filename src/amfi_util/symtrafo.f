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
      Subroutine SymTrafo(LUPROP,ip,lOper,nComp,nBas,nIrrep,Label,
     &                    MolWgh,SOInt,LenTot)
cbs
cbs   Purpose: combine SO-integrals from amfi to symmetry-adapted
cbs   integrals on one file AOPROPER_MF_SYM
cbs
cbs
      implicit real*8 (a-h,o-z)
#include "para.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "real.fh"
      parameter(maxorbs=MxOrb)
      parameter(maxcent=MxAtom)
      Real*8 SOInt(LenTot)
      character*8 xa,ya,za,xa2
      character*3 END
CBS   character*20 filename
      logical  EX
CBS   namelist /SYMTRA/ none
      dimension xa(4),ya(4),za(4)
      dimension xa2(4)
      dimension ncent(maxorbs), Lval(maxorbs),mval(maxorbs),
     *nadpt(maxorbs),nphase(8,maxorbs),idummy(8),
     *Lhighcent(maxcent),Lcent(MxCart),Mcent(MxCart),
     *ncontcent(0:Lmax)   ,
     *numballcart(maxcent),ifirstLM(0:Lmax,-Lmax:Lmax,maxcent)
      Integer ip(nComp), nBas(0:nIrrep-1), lOper(nComp), ipC(3,MxAtom)
      Character Label*8
c###############################################################################
      IPNT(I,J)=(max(i,j)*max(i,j)-max(i,j))/2 +min(i,j)
*
      END='   '  ! added due to cray warnings. B.S. 04/10/04
      xa(1)='********'
      ya(1)='********'
      za(1)='********'
      xa(2)='        '
      ya(2)='        '
      Za(2)='        '
      xa(3)='ANTISYMM'
      ya(3)='ANTISYMM'
      Za(3)='ANTISYMM'
      xa(4)='X1SPNORB'
      ya(4)='Y1SPNORB'
      ZA(4)='Z1SPNORB'
c
c     read information from SYMINFO
      isymunit=isfreeunit(58)
      call f_inquire('SYMINFO',EX)
      if (.not.EX) Call SysAbendMsg('systrafo',
     &                  'SYMINFO not present','Sorry')
      call molcas_open(isymunit,'SYMINFO')
      rewind(isymunit)
*define _DEBUG_
#ifdef _DEBUG_
      write(6,*) 'Symmetry adapation of the SO-integrals'
#endif
      read(isymunit,*)
      read(isymunit,*)
      read(isymunit,*)
      numboffunct=0
      do while(END.ne.'END')
         numboffunct=numboffunct+1
         read(isymunit,'(A3)') END
      enddo
#ifdef _DEBUG_
      write(6,*) 'there are totally ',numboffunct,' functions'
#endif
      if (numboffunct.gt.maxorbs)
     &Call SysAbendMsg('symtrafo','increase maxorbs in symtrafo',' ')
      rewind isymunit
      read(isymunit,*)
      read(isymunit,*)
      numbofcent=0
      do irun=1,numboffunct
         read(isymunit,*) index,ncent(irun),lval(irun),
     &                    mval(irun),nadpt(irun),
     &                    (nphase(I,irun),I=1,nadpt(irun))
         numbofcent=max(numbofcent,ncent(irun))
         if (index.ne.irun)
     &      Call SysAbendMsg('symtrafo',
     &                       'weird numbering  on SYMINFO',' ' )
      End Do
      Close(iSymUnit)
#ifdef _DEBUG_
      write(6,*) 'number of unique centres' , numbofcent
#endif
c
c     clean up arrays for new integrals
      numboffunct3=(numboffunct*numboffunct+numboffunct)/2
*
      Call GetMem('AMFII','Allo','Real',ipInt,numboffunct3*3)
      ipX=ipInt-1
      ipY=ipX+numboffunct3
      ipZ=ipY+numboffunct3
      call dcopy_(numboffunct3*3,[Zero],0,Work(ipInt),1)
*
      nSOs=0
      Do iIrrep = 0, nIrrep-1
         nSOs=nSOs+nBas(iIrrep)
      End Do
      Call GetMem('iSO','Allo','Inte',ipiSO,nSOs*2)
      iSO_a=0
      Do iIrrep=0,nIrrep-1
         iSO_r=0
         Do iBas = 1, nBas(iIrrep)
            iSO_a=iSO_a+1
            iSO_r=iSO_r+1
            iWork(ipiSO+(iSO_a-1)*2  )=iIrrep
            iWork(ipiSO+(iSO_a-1)*2+1)=iSO_r
         End Do
      End Do
*
*     loop over unique centres to read integrals and information
*
      iunit=LUPROP
      nSCR=numboffunct3*3
      Call GetMem('SCR','Allo','Real',ipSCR,nSCR)
      Call FZero(Work(ipSCR),nSCR)
      ipSCRX=ipSCR
      ipSCRY=ipSCRX+numboffunct3
      ipSCRZ=ipSCRY+numboffunct3
      length3_tot=0
*
*     In a MPI run not all atomic block will be available in
*     all processes. Make up so we know later if a particular
*     atom is present.
*
      do jcent=1,numbofcent
         ipC(1,jCent)=ip_Dummy
      End Do
      do jcent=1,numbofcent
*
#ifdef _DEBUG_
      write(6,*) 'read integrals and info for centre ',jcent
#endif
*
*        Note that when running in parallel this list is incomplete.
*        Hence, we process the centers which each process host.
*
              read(iunit,END=199)  iCent
              read(iunit)  xa2
     &        ,numbofsym,(idummy(I),
     &        i=1,numbofsym),
     &        numballcart(icent),(Lcent(i),
     &        I=1,numballcart(icent)),
     &        (mcent(i),I=1,numballcart(icent)),
     &        Lhighcent(icent),(ncontcent(I),I=0,Lhighcent(icent))
#ifdef _DEBUG_
              write(6,*) numballcart(icent) ,
     &        'functions on centre ',icent
#endif
              length3=ipnt(numballcart(icent),numballcart(icent))
              ipC(1,iCent)=ipSCRX
              ipC(2,iCent)=ipSCRY
              ipC(3,iCent)=ipSCRZ
              read(iunit) (Work(i),i=ipSCRX,ipSCRX+length3-1)
              read(iunit)  Ya
              read(iunit) (Work(i),i=ipSCRY,ipSCRY+length3-1)
              read(iunit)  Za
              read(iunit) (Work(i),i=ipSCRZ,ipSCRZ+length3-1)
              ipSCRX=ipSCRX+length3
              ipSCRY=ipSCRY+length3
              ipSCRZ=ipSCRZ+length3
              length3_tot=length3_tot+length3
culf
c      check if any L-value is missing
      LLhigh = Lhighcent(icent)
      do i=1,Lhighcent(icent)
         if(ncontcent(I).eq.0) LLhigh=LLhigh-1
      enddo
      Lhighcent(icent)=LLhigh
cbs   determize where the first function of a special type is..
      not_defined=ipnt(numboffunct,numboffunct)+1
      do Lrun=0,Lhighcent(icent)
         do Mrun=-Lrun,Lrun
            ifirstLM(Lrun,Mrun,icent)=not_defined
         enddo
      enddo
      do iorb=1,numballcart(icent)
         Lrun=Lcent(iorb)
         Mrun=Mcent(iorb)
#ifdef _DEBUG_
         write(6,*) 'iorb,Lrun,mrun',iorb,Lrun,mrun
#endif
         ifirstLM(Lrun,Mrun,icent)=min(iorb,ifirstLM(Lrun,Mrun,icent))
      enddo
*
cbs   determined..
cbs   check if all of them were found
      do Lrun=0,Lhighcent(icent)
         do Mrun=-Lrun,Lrun
            if (ifirstLM(Lrun,Mrun,icent).eq.not_defined) then
               write(6,*) 'problems for centre,L,M ',icent,Lrun,Mrun
               Call SysAbendMsg('symtrafo',
     &                          'problems with L- and M-values',' ')
            endif
          enddo
      enddo
      enddo    !end of loop over centres
 199  Continue
#ifdef _DEBUG_
      ipSCRX=ipSCR
      ipSCRY=ipSCRX+numboffunct3
      ipSCRZ=ipSCRY+numboffunct3
      Write (6,*) 'length3_tot=',length3_tot
      Call RecPrt('SCRX',' ',Work(ipSCRX),1,length3_tot)
      Call RecPrt('SCRY',' ',Work(ipSCRY),1,length3_tot)
      Call RecPrt('SCRZ',' ',Work(ipSCRZ),1,length3_tot)
      Do iCent = 1, numbofcent
         Write (6,*) ipC(1,iCent),ipC(2,iCent),ipC(3,iCent)
      End Do
#endif
*     If this process does not have any blocks of integrals proceed
*     directly to the distribution step.
      If (Length3_tot.eq.0) Go To 299
cbs
cbs   Finally the transformation!!!!
cbs
cbs
      icentprev=0
      jcentprev=0
      ilcentprev=-1
      jlcentprev=-1
      imcentprev=20
      jmcentprev=20
      isame=1
      jsame=1
      lauf=0
      do irun=1,numboffunct
*        Skip if center corresponding to this basis function is
*        not available at this process.
         if (ncent(irun).eq.icentprev.AND.ilcentprev.eq.lval(irun)
     &       .AND.imcentprev.eq.mval(irun)) then
            isame=isame+1
         else
            isame=1
            icentprev=ncent(irun)
            ilcentprev=lval(irun)
            imcentprev=mval(irun)
         endif
      do jrun=1,irun
      lauf=lauf+1
      if (ncent(jrun).eq.jcentprev.and.jlcentprev.eq.lval(jrun)
     &    .AND.jmcentprev.eq.mval(jrun)) then
         jsame=jsame+1
      else
         jsame=1
         jcentprev=ncent(jrun)
         jlcentprev=lval(jrun)
         jmcentprev=mval(jrun)
      endif
cbs   check for same centers
      if (ncent(irun).eq.ncent(jrun)) then
      if (lval(irun).eq.lval(jrun).and.lval(irun).gt.0) then
      if (iabs(iabs(mval(irun))-iabs(mval(jrun))).le.1) then
*
*
*
cbs   the only cases  where non-zero integrals occur
      if (nadpt(irun).eq.1) then
         coeff=1d0
      else
         icoeff=0
         do icc=1,nadpt(irun)
            icoeff=icoeff+nphase(icc,irun)*nphase(icc,jrun)
         enddo
         coeff=DBLE(icoeff)
         If (MolWgh.eq.2) Then
            coeff=coeff/DBLE(nadpt(irun))
         Else
            coeff=coeff/DBLE(nadpt(irun)*nadpt(irun))
         End If
      endif
cbs   determine indices of atomic integrals
      indexi=ifirstLM(lval(irun),mval(irun),ncent(irun))+isame-1
      indexj=ifirstLM(lval(irun),mval(jrun),ncent(irun))+jsame-1
      laufalt=ipnt(indexi,indexj)
*
      If (ipC(1,nCent(iRun)).ne.ip_Dummy) Then
      ipSCRX=ipC(1,ncent(irun))-1+laufalt
      ipSCRY=ipC(2,ncent(irun))-1+laufalt
      ipSCRZ=ipC(3,ncent(irun))-1+laufalt
cDebugDebug
c     Write (*,*) 'laufalt=',laufalt
c     Write (*,*) 'ip''s:',ipSCRX,ipSCRY,ipSCRZ
c     Write (*,*) Work(ipSCRX),Work(ipSCRY),Work(ipSCRZ)
cDebugDebug
      if (indexi.gt.indexj) then
         Work(ipX+lauf)=-coeff*Work(ipSCRX)
         Work(ipY+lauf)=-coeff*Work(ipSCRY)
         Work(ipZ+lauf)=-coeff*Work(ipSCRZ)
      else
         Work(ipX+lauf)=coeff*Work(ipSCRX)
         Work(ipY+lauf)=coeff*Work(ipSCRY)
         Work(ipZ+lauf)=coeff*Work(ipSCRZ)
      endif
      End If
*
*
      endif
      endif
      endif
      enddo ! jrun
      enddo ! irun
 299  Continue
*     This test is not valid for parallel execusion, since here
*     we have only incomplete lists.
*     if (lauf.ne.numboffunct3)
*    & Call SysAbendMsg('symtrafo', 'error in numbering ',' ' )
      Do iComp = 1, nComp
         Do iSO=1, numboffunct
            j1=iWork(ipiSO+(iSO-1)*2  )
            iSO_r=iWork(ipiSO+(iSO-1)*2+1)
            Do jSO = 1, iSO
               j2=iWork(ipiSO+(jSO-1)*2  )
               jSO_r=iWork(ipiSO+(jSO-1)*2+1)
               j12=iEor(j1,j2)
               If (iAnd(lOper(iComp),2**j12).eq.0) Go To 99
*
               iOff = iPntSO(Max(j1,j2),Min(j1,j2),lOper(iComp),nBas)
               iOff = iOff + ip(iComp)-1
*
               iOff2=ipInt + (iComp-1)*numboffunct3 -1 + iPnt(iSO,jSO)
*
               If (j1.eq.j2) Then
                  ijSO=iPnt(iSO_r,jSO_r)
                  SOInt(iOff+ijSO)=-Work(iOff2)
               Else
                  nRow=nBas(j1)
                  SOInt(iOff+(jSO_r-1)*nRow+iSO_r)=-Work(iOff2)
               End If
*
 99            Continue
            End Do
         End Do
*
*------- Write out integrals to ONEINT for this specific component of the
*        operator.
*
*
         iOpt=0
         iRC=-1
         iSmLbl=lOper(iComp)
         Call GADSum(SOInt(ip(iComp)),n2Tri(iSmLbl))
         Call WrOne(iRC,iOpt,Label,iComp,SOInt(ip(iComp)),iSmLbl)
         If (iRC.ne.0) then
            Call SysAbendMsg('symtrafo',
     &            '     Error in subroutine ONEEL ',
     &      '     Abend in subroutine WrOne')
         End If
*
      End Do
*
#ifdef _DEBUG_
      Call PrMtrx(Label,lOper,nComp,ip,SOInt)
#endif
*
      Call GetMem('iSO','Free','Inte',ipiSO,nSOs*2)
      Call GetMem('SCR','Free','Real',ipSCR,nSCR)
      Call GetMem('AMFII','Free','Real',ipInt,numboffunct3*3)
CBS   write(6,*) 'Symmetry transformation successfully done'
      Return
      End
