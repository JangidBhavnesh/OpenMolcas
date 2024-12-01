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
Module cicisp_mclr
#include "detdim.fh"
#include "cicisp_mclr.fh"
save
!
! icisps **
! smost ***
!     niCISp                :        Number of CI-spaces (1) *
!     iAStFI                :       Alpha string type for internal CI space *
!     iBStFI                :       Beta string type for internal CI space (iBZTP)*
!     iActi                :        Active/Inactive space  (1)*
!     MnR1IC                :         Minimum number of electrons in RAS1*
!     MxR1IC                :        Maximum    -              "             - RAS1*
!     MnR3IC                :       Minimum    -             "       - RAS3*
!     MxR3IC                :        Maximum    -              "             - RAS3*
!     iZCI                :        Internal zero order space *
!     iRCI                :        Number of zero order space *
!     nElCI                :           Number of electrons per CI space *
!     nAElCI                :           Number of alpha electrons per CI space *
!     nBElCI                :           Number of beta electrons per CI space *
!     xispsm                :       Number of det. for each (symmetry,CI space) **
!     ismost                :       Symmetry operator ASym=ismost(BSym,iTOTSM) ***(istead of ieor)
!     MXSB                 :       Largest symmetry block **
!     MXSOOB                :       Largest block        **
!
!Integer NICISP,IASTFI(MXPICI),IBSTFI(MXPICI),IACTI(MXPICI),MNR1IC(MXPICI),MXR1IC(MXPICI),    &
!        MNR3IC(MXPICI),MXR3IC(MXPICI),IZCI,IRCI(3,7,7),NELCI(MXPICI),NAELCI(MXPICI),NBELCI(MXPICI), &
!        ISMOST(MXPCSM,MXPCSM),MXSB,MXSOOB,ldet,lcsf
!Real*8 XISPSM(MXPCSM,MXPICI)
Private MXPIRR,MXPOBS,MXPR4T,MXINKA,MXPORB,MXPXOT,MXPXST,MXPSHL, &
        MXPL,MXPXT,MXPICI,MXPSTT,MXPCSM,MXPCTP,MXCNSM,MXPWRD, &
        MXNMS,MTYP,MXPNGAS,MXPNSMST,MXPPTSPC

End Module cicisp_mclr
