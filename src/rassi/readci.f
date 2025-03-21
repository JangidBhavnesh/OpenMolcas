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
      SUBROUTINE READCI(ISTATE,SGS,CIS,NCI,CI)
      use rassi_aux, only: ipglob
      use rassi_global_arrays, only: JBNUM, LROOT
      use gugx, only: SGStruct, CIStruct
#ifdef _HDF5_
      USE mh5, ONLY: mh5_is_hdf5, mh5_open_file_r, mh5_exists_attr,
     &               mh5_fetch_attr, mh5_fetch_dset, mh5_close_file
      use Cntrl, only: NROOTS
#endif
      use Cntrl, only: NSTATE, PRCI, CITHR, IRREP, JBNAME, MLTPLT
      use cntrl, only: iTOC15, LuIph

      IMPLICIT NONE
#ifdef _HDF5_
! pick up MxRoot
#include "rasdim.fh"
      integer :: refwfn_id
      integer :: root2state(mxroot), IDXCI
#endif

      INTEGER ISTATE
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      INTEGER NCI
      REAL*8 CI(NCI)

      INTEGER I, IAD, IDISK, JOB, LROOT1, LSYM


      IF(ISTATE.LT.1 .OR. ISTATE.GT.NSTATE) THEN
        WRITE(6,*)'RASSI/READCI: Invalid ISTATE parameter.'
        WRITE(6,*)' ISTATE, NSTATE:',ISTATE,NSTATE
        CALL ABEND()
      END IF
      JOB=JBNUM(ISTATE)
      LROOT1=LROOT(ISTATE)

#ifdef _HDF5_
************************************************************************
*
* For HDF5 formatted job files
*
************************************************************************
      If (mh5_is_hdf5(jbname(job))) Then
        refwfn_id = mh5_open_file_r(jbname(job))
        if (mh5_exists_attr(refwfn_id, 'ROOT2STATE')) then
          call mh5_fetch_attr(refwfn_id,'ROOT2STATE',root2state)
          IDXCI = root2state(lroot1)
        else
          IDXCI = lroot1
        end if
        IF (IDXCI.LE.0.OR.IDXCI.GT.NROOTS(JOB)) THEN
          call WarningMessage(2,'Invalid CI array index, abort!')
          call AbEnd
        END IF
        call mh5_fetch_dset(refwfn_id,
     &         'CI_VECTORS',CI,[NCI,1],[0,IDXCI-1])
        call mh5_close_file(refwfn_id)
      Else
#endif
************************************************************************
*
* For JOBIPH/JOBMIX formatted job files
*
************************************************************************
        CALL DANAME(LUIPH,JBNAME(JOB))
        IAD=0
        CALL IDAFILE(LUIPH,2,ITOC15,30,IAD)
        IDISK=ITOC15(4)
        DO I=1,LROOT1-1
          CALL DDAFILE(LUIPH,0,CI,NCI,IDISK)
        END DO
        CALL DDAFILE(LUIPH,2,CI,NCI,IDISK)
        CALL DACLOS(LUIPH)
************************************************************************
#ifdef _HDF5_
      End If
#endif

      IF(IPGLOB.gt.0 .and. PRCI) THEN
        WRITE(6,*)' READCI called for state ',ISTATE
        WRITE(6,*)' This is on JobIph nr.',JOB
        WRITE(6,*)' JobIph file name:',JBNAME(JOB)
        WRITE(6,*)' It is root nr.',LROOT(ISTATE)
        WRITE(6,*)' Its length NCI=',NCI
        WRITE(6,*)' Its symmetry  =',IRREP(JOB)
        WRITE(6,*)' Spin multiplic=',MLTPLT(JOB)
        LSYM=IRREP(JOB)
        CALL PRWF(SGS,CIS,LSYM,CI,CITHR)
      END IF

      END SUBROUTINE READCI
