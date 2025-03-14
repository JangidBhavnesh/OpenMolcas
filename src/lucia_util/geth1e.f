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
      FUNCTION GETH1E(IORB,ITP,ISM,JORB,JTP,JSM)
      use GLBBAS, only: INT1, PINT1, PGINT1, PGINT1A
      use lucia_data, only: IH1FORM
      use symmetry_info, only: MULTD2H => Mul
      use lucia_data, only: IBSO,IOBPTS,IREOTS,NACOBS
      use lucia_data, only: MXPNGAS
*
* One-electron integral for active
* orbitals (IORB,ITP,ISM),(JORB,JTP,JSM)
*
* The orbital symmetries are used to obtain the
* total symmetry of the operator
      IMPLICIT None
      Integer IORB,ITP,ISM,JORB,JTP,JSM
      INTEGER IJSM
      REAL*8 GETH1E
      REAL*8, EXTERNAL:: GTH1ES
*
      IJSM = MULTD2H(ISM,JSM)
      GETH1E = 0.0D0
      IF(IH1FORM.EQ.1) THEN
*. Normal integrals, lower triangular packed
        IF(IJSM.EQ.1) THEN
C?        WRITE(6,*) ' GETH1E, old route '
          GETH1E =
     &    GTH1ES(IREOTS,PINT1,INT1,IBSO,MXPNGAS,
     &              IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
        ELSE
          GETH1E =
     &    GTH1ES(IREOTS,PGINT1(IJSM)%I,INT1,IBSO,MXPNGAS,
     &              IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
        END IF
      ELSE
*. Integrals are in full blocked form
        GETH1E =
     &  GTH1ES(IREOTS,PGINT1A(IJSM)%I,INT1,IBSO,MXPNGAS,
     &         IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,0)
      END IF
*
      END FUNCTION GETH1E
