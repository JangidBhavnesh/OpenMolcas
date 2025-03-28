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
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************

#include "compiler_features.h"

#ifndef _HAVE_EXTRA_

* Broadcast a file from the master to the slaves

      Subroutine PFGet_ASCII(FName)
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: mpp_rootid, King
#endif
      Implicit None
      Character (Len=*), Intent(In) :: FName
#ifdef _MOLCAS_MPP_
#include "SysDef.fh"
#include "mafdecls.fh"
      Integer, Parameter :: LBuf=4096
      Character (Len=LBuf) :: Buf
      Integer :: LU, Err, FLen, Pos, Num
      Logical :: Found, Failed
      Integer, External :: IsFreeUnit
      Interface
        Subroutine GA_Brdcst(type,buf,lenbuf,root)
          Integer type,lenbuf,root
          Type(*) buf
        End Subroutine GA_Brdcst
      End Interface

      ! Note that each process opens only one file, so there is a single
      ! unit number LU
      LU=10
      LU=IsFreeUnit(LU)
      ! Check file existence and read size on the master
      If (King()) Then
        Call f_Inquire(FName, Found)
        If (Found) Then
          Call Molcas_Open_Ext2(LU, FName, "stream", "unformatted", Err,
     &                          .False., 0, "old", Failed)
          If (Failed .or. (Err .ne. 0)) Then
            Write(6,*) "Failed to open file ", Trim(FName)
            Call AbEnd()
          End If
          Inquire(LU, Size=FLen)
        Else
          FLen=0
        End If
      End If
      ! Broadcast the file size
      Err=0
      Call GA_Brdcst(MT_INT, FLen, 1*ItoB, mpp_rootid)
      If (FLen .le. 0) Return
      ! Open file for writing in the slaves
      If (.not.King()) Then
        Call Molcas_Open_Ext2(LU, FName, "stream", "unformatted", Err,
     &                        .False., 0, "replace", Failed)
        If (Failed .or. (Err .ne. 0)) Then
          Write(6,*) "Failed to open file ", Trim(FName)
          Call AbEnd()
        End If
      End If
      ! Pass the file content in chunks
      Pos=0
      Do While (Pos .lt. FLen)
        ! Length of this chunk
        Num = Min(LBuf, FLen-Pos)
        ! The master reads the file
        If (King()) Then
          Read(LU, IOStat=Err) Buf(1:Num)
          If (Err .ne. 0) Then
            Write(6,*) "Error reading the file ", Trim(FName)
            Call AbEnd()
          End If
        End If
        Call GA_Brdcst(MT_BYTE, Buf(1:Num), Num, mpp_rootid)
        ! The slaves write the file
        If (.not. King()) Then
          Write(LU, IOStat=Err) Buf(1:Num)
          If (Err .ne. 0) Then
            Write(6,*) "Error writing the file ", Trim(FName)
            Call AbEnd()
          End If
        End If
        Pos = Pos + LBuf
      End Do
      Close(LU)
      Return
#else
      ! Avoid unused argument warnings
      If (.False.) Call Unused_Character(FName)
#endif

      End Subroutine PFGet_ASCII

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#     include "macros.fh"
      subroutine empty_PFGet_ASCII()
      end subroutine empty_PFGet_ASCII

#endif
