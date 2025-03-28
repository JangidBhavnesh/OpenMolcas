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
      SUBROUTINE GETINT(    XINT,     ITP,     ISM,     JTP,     JSM,
     &                       KTP,     KSM,     LTP,     LSM,  IXCHNG,
     &                      IKSM,    JLSM,   ICOUL)
      use GLBBAS, only: PINT2, KINH1
      use lucia_data, only: NSMOB
      use lucia_data, only: NOBPTS,NTOOBS
*
* Outer routine for accessing integral block
*
      IMPLICIT NONE
*
      INTEGER  ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,ICOUL
      REAL*8 XINT(*)
*
      INTEGER NTEST,NI,NK,NIK,NJ,NL,NJL
      NTEST = 00
*
      IF(NTEST.GE.1) THEN
c       WRITE(6,*) ' I_USE_SIMTRH in GETINT =', I_USE_SIMTRH
       WRITE(6,*) ' GETINT : ICOUL = ', ICOUL
       WRITE(6,*)       'ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM : '
       WRITE(6,'(8I4)')  ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM
      END IF
*. Read integrals in in RASSCF format
        CALL GETINCN_RASSCF(   XINT,    ITP,    ISM,    JTP,    JSM,
     &                          KTP,    KSM,    LTP,    LSM, IXCHNG,
     &                       IKSM,JLSM,PINT2,NSMOB,KINH1,
     &                        ICOUL)

      IF(NTEST.NE.0) THEN
        IF(ITP.EQ.0) THEN
          NI = NTOOBS(ISM)
        ELSE
          NI = NOBPTS(ITP,ISM)
        END IF
*
        IF(KTP.EQ.0) THEN
          NK = NTOOBS(KSM)
        ELSE
          NK = NOBPTS(KTP,KSM)
        END IF
*
        IF(IKSM.EQ.0) THEN
          NIK = NI * NK
        ELSE
          NIK = NI*(NI+1)/2
        END IF
*
        IF(JTP.EQ.0) THEN
          NJ = NTOOBS(JSM)
        ELSE
          NJ = NOBPTS(JTP,JSM)
        END IF
*
        IF(LTP.EQ.0) THEN
          NL = NTOOBS(LSM)
        ELSE
          NL = NOBPTS(LTP,LSM)
        END IF
*
        IF(JLSM.EQ.0) THEN
          NJL = NJ * NL
        ELSE
          NJL = NJ*(NJ+1)/2
        END IF
        WRITE(6,*) ' 2 electron integral block for TS blocks '
        WRITE(6,*) ' Ixchng :', IXCHNG
        WRITE(6,*) ' After GETINC '
        WRITE(6,'(1X,4(A,I2,A,I2,A))')
     &  '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &  '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
        CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
      END IF
*
      END SUBROUTINE GETINT
