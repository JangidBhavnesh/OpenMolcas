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
* Copyright (C) 2020, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine Find_Distance(Ref,Point,Dir,Fact,Dist,nAtom,
     &                         BadConstraint)
      use Slapaf_Info, only: RefGeo
      Implicit None
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "info_slapaf.fh"
#include "sbs.fh"
      Integer, Intent(In) :: nAtom
      Real*8, Intent(In) :: Ref(3,nAtom),Dir(3,nAtom),Fact,Dist
      Real*8, Intent(Out) :: Point(3,nAtom)
      Logical, Intent(Out) :: BadConstraint
      Real*8, Allocatable :: OldRef(:,:),Dummy(:)
      Real*8 :: R,CurFact,PrevR,Correct
      Real*8, Parameter :: Thr = 1.0d-6
      Integer :: nCoor,i
      Logical :: Invar
*                                                                      *
************************************************************************
*                                                                      *
      Invar=(iAnd(iSBS,2**7).eq.0).and.(iAnd(iSBS,2**8).eq.0)
      nCoor=3*nAtom
      Call mma_allocate(OldRef,3,nAtom,Label='OldRef')
      Call mma_allocate(Dummy,nCoor,Label='Dummy')
      OldRef(:,:) = RefGeo(:,:)
      RefGeo(:,:)    = Ref(:,:)

      R=Zero
      CurFact=Zero
      Correct=Fact
      i=0
      Do While (Abs(One-R/Dist).gt.Thr)

*       Add the scaled direction vector
        CurFact=CurFact+Correct
        Point(:,:) = Ref(:,:) + CurFact * Dir(:,:)

*       Align and measure distance
        PrevR=R
        Call Align(Point(:,:),Ref(:,:),nAtom)
        If (MEP_Type.eq.'SPHERE') Then
          Call SphInt(Point,nAtom,ip_Dummy,R,Dummy,
     &                .False.,'dummy   ',Work(ip_Dummy),.False.)
        Else If (MEP_Type.eq.'TRANSVERSE') Then
          Call Transverse(Point,nAtom,R,Dummy,
     &                .False.,'dummy   ',Work(ip_Dummy),.False.)
        End If

*       Stop if too many iterations or if the constraint is moving
*       in the wrong direction
        i=i+1
        If (i.gt.5.or.Correct*(R-PrevR).lt.Zero) Exit
        Correct=(One-R/Dist)*Fact
      End Do
      BadConstraint=(Abs(One-R/Dist).gt.Thr)

      RefGeo(:,:) = OldRef(:,:)
      Call mma_deallocate(OldRef)
      Call mma_deallocate(Dummy)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End Subroutine Find_Distance
