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
      SUBROUTINE INTCTL2(IF_TRNSF)
      use caspt2_output, only: iPrGlb
      use caspt2_gradient, only: do_grad, nStpGrd
      use caspt2_data, only: CMO, FIMO
      use PrintLevel, only: debug
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "intgrl.fh"
#include "caspt2_grad.fh"

      LOGICAL IF_TRNSF


* Compute using Cholesky vectors.
* Frozen, inactive and active Fock matrix in AO basis:
      Call GetMem('FFAO','ALLO','REAL',LFFAO,NBTRI)
      Call GetMem('FIAO','ALLO','REAL',LFIAO,NBTRI)
      Call GetMem('FAAO','ALLO','REAL',LFAAO,NBTRI)
* tracho2 makes many allocations but should deallocate everything
* before its return.
      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' INTCTL2 calling TRACHO2...'
        CALL XFLUSH(6)
      END IF
      Call TraCho2(CMO,Work(LDREF),
     &             Work(LFFAO),Work(LFIAO),Work(LFAAO),IF_TRNSF)
      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' INTCTL2 back from TRACHO2.'
        CALL XFLUSH(6)
      END IF
* All extra allocations inside tracho2 should now be gone.

* For gradient calculation, it is good to have FIAO and FAAO
      IF (do_grad.or.nStpGrd.eq.2) THEN
        !! FFAO has one-electron Hamiltonian
        CALL DCOPY_(NBTRI,WORK(LFFAO),1,WORK(ipFIMO),1)
        CALL DAXPY_(NBTRI,1.0D+00,WORK(LFIAO),1,WORK(ipFIMO),1)
        CALL DCOPY_(NBTRI,WORK(ipFIMO),1,WORK(ipFIFA),1)
        CALL DAXPY_(NBTRI,1.0D+00,WORK(LFAAO),1,WORK(ipFIFA),1)
      END IF
* Transform them to MO basis:
      CALL DCOPY_(notri,[0.0D0],0,WORK(LHONE),1)
      FIMO(:)=0.0D0
      CALL DCOPY_(notri,[0.0D0],0,WORK(LFAMO),1)
c Compute FIMO, FAMO, ...  to workspace:
      Call FMat_Cho(CMO,Work(LFFAO),Work(LFIAO),Work(LFAAO),
     &              Work(LHONE),FIMO,Work(LFAMO))
      Call GetMem('FFAO','FREE','REAL',LFFAO,NBTRI)
      Call GetMem('FIAO','FREE','REAL',LFIAO,NBTRI)
      Call GetMem('FAAO','FREE','REAL',LFAAO,NBTRI)

      RETURN
      END
