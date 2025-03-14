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
      INTEGER FUNCTION ISYMCN_MCLR(ICL,IOP,NCL,NOPEN)
*
* Master routine for symmetry of configuration
* with NCL doubly occupied orbitals and NOPEN singly occupied shells
*
      use MCLR_Data, only: ISMFTO
      IMPLICIT None
      INTEGER ICL(*),IOP(*)
      INTEGER NCL,NOPEN

      INTEGER IEL,IVV,JVV,KVV
*
      ISYMCN_MCLR = 1
      DO IEL = 1, NOPEN
        IVV=ISYMCN_MCLR-1
        JVV=ISMFTO(IOP(IEL))-1
        KVV = IEOR(IVV,JVV)
        ISYMCN_MCLR = KVV+1
      END DO
*
c Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_integer_array(ICL)
         CALL Unused_integer(NCL)
      END IF
      END FUNCTION ISYMCN_MCLR
