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
Module cstate_mclr
#include "cstate_mclr.fh"
save
!
! Stuff from cstate_mclr.fh
!        iRefSm        :        Reference symmetry
!        iRefML        :        ML value
!        iRefPA        :        Parity
!        iRefL         :        L value
!        MS2           :        2*MS
!        MULTS         :        Spin multiplicity
!        nRoot         :        Number of roots
!        IDC           :        Ms combinations
!        PSSIGN        :        Ms combination PS factor
!        PLSIGN        :        Ms combination PL factor
!        IntSel        :        Internal space
!        iAlign        :        Not in use
!        Ethers        :        Threshold (energy cont)
!        Cthres        :        Threshold (coeff)
!        NGenSym       :        Number of reference symmetries
!        iGenSym       :        Reference symmetries
!        InvSym        :        Ger/UnGer inv sym
!        iKram         :        Kramer symmetry
!. CSTATE
!Integer IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS,NROOT,IDC,INTSEL,IALIGN,   &
!        NGENSYM,IGENSYM(100),INVSYM,IKRAM
!Real*8 PSSIGN,PLSIGN, ETHRES,CTHRES
End Module cstate_mclr
