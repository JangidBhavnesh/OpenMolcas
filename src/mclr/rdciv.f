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
      Subroutine RdCIV()
************************************************************************
*                                                                      *
*     Read the contents of the JOBIPH file.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use negpre, only: ERAS, LuCIV
      use gugx, only: SGS, CIS, EXS
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
#include "Input.fh"
#include "Files_mclr.fh"
      Real*8, Allocatable:: OCIvec(:), Tmp(:)
      Integer iDisk, iDisk1, i, Iter
      Real*8 Temp

      Call DaName(LuCIV,'ROOTS')
      iDisk=0
*----------------------------------------------------------------------*
*     Load the CI vector for all roots                            *
*----------------------------------------------------------------------*

      Call mma_Allocate(OCIvec,nConf,Label='OCIvec')

      iDisk=iToc(4)
      idisk1=0
      Do i=1,lroots
       Call dDaFile(LuJob,2,OCIvec,nConf,iDisk)
       Call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,
     &              nRs1,nRs2,nRs3,
     &              SGS,CIS,EXS,OCIvec,1,State_Sym,State_Sym)
       NCSF(1:nSym)=CIS%NCSF(1:nSym)
       NCONF=CIS%NCSF(State_Sym)
       Call mkGuga_Free(SGS,CIS,EXS)

       Call dDafile(LuCIV,1,OCIvec,nconf,iDisk1)
      End Do

      Call mma_deAllocate(OCIvec)
*----------------------------------------------------------------------*
*     Load state energy                                                *
*----------------------------------------------------------------------*

      Call mma_allocate(Tmp,mxRoot*mxIter,Label='Tmp')
      iDisk=iToc(6)
      Call dDaFile(LuJob,2,Tmp,mxRoot*mxIter,iDisk)
      Do i=0,lroots-1
      ERAS(i+1)=0.0d0
      Do  iter=1,mxIter
         Temp=Tmp(iter*mxRoot+i)
         If ( Temp.ne.0.0D0 ) ERAS(i+1)=Temp
      End Do
      End Do
      Call mma_deallocate(Tmp)
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      End Subroutine RdCIV
