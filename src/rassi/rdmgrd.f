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
      SUBROUTINE RDMGRD(JOB,IDISP,LABEL,STYPE,ISYMP,NARRAY,ARRAY)
      use rassi_aux, only: ipglob
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NJOB, MINAME
      use cntrl, only: LuMck
      use Symmetry_Info, only: nSym=>nIrrep, MUL
      use rassi_data, only: NBASF

      IMPLICIT None
C Purpose: Read in the derivatives of 1-electron integrals
C of some operator, with respect to some displacement IDISP.
C ISYMP is the symmetry irrep label of the derivatives.
      Integer JOB,IDISP,ISYMP,NARRAY
      Real*8 ARRAY(NARRAY)
      CHARACTER(LEN=8) LABEL,STYPE

      Integer ITOFF(8),IAOFF(8)
      Real*8, Allocatable:: TEMP(:)
      INTEGER IRC, IOPT, ISUM, IS, JS, NBI, NBJ, NBIJ, NTEMP, ISCODE,
     &        LT, IA1, J, I, IA2
      REAL*8 F

      IF(JOB.LT.1 .OR. JOB.GT.NJOB) THEN
        WRITE(6,*)' RASSI/RDMGRD: Invalid JOB parameter.'
        WRITE(6,*)' JOB:',JOB
        CALL ABEND()
      END IF

      IF(IPGLOB.GT.3) THEN
        WRITE(6,*)' RDMGRD called for JOB=',JOB
        WRITE(6,*)' perturbed by displacement nr.',IDISP
        WRITE(6,*)' MckInt file name:',MINAME(JOB)
        WRITE(6,*)' Operator name LABEL=',LABEL
        WRITE(6,*)' Symmetry type STYPE=',STYPE
        WRITE(6,*)' Irrep label   ISYMP=',ISYMP
        WRITE(6,*)' Length NARRAY=',NARRAY
      END IF

C Open MCKINT file:
      IRC=-1
      IOPT=0
      CALL OPNMCK(IRC,IOPT,MINAME(JOB),LUMCK)
      IF(IRC.NE.0) THEN
        WRITE(6,*)'RASSI/RDMGRD: Failed to open '//MINAME(JOB)
        WRITE(6,*)'Unit nr LUMCK=',LUMCK
        WRITE(6,*)'Option code IOPT=',IOPT
        WRITE(6,*)'Return code IRC =',IRC
        CALL ABEND()
      END IF

C Addressing integral blocks in the buffer:
      ISUM=0
      DO IS=1,NSYM
       JS=MUL(IS,ISYMP)
       IF(IS.GE.JS) THEN
        ITOFF(IS)=ISUM
        ITOFF(JS)=ISUM
        NBI=NBASF(IS)
        NBJ=NBASF(JS)
        NBIJ=NBI*NBJ
        IF(IS.EQ.JS) NBIJ=(NBIJ+NBI)/2
        ISUM=ISUM+NBIJ
       END IF
      END DO
      NTEMP=ISUM
C Read MCKINT file:
      IOPT=0
      ISCODE=2**(ISYMP-1)
C Get temporary buffer to read data by RDMCK calls
      CALL mma_allocate(TEMP,NTEMP,Label='TEMP')
C Read 1-electron integral derivatives:
      IRC=NTEMP
      CALL dRDMCK(IRC,IOPT,LABEL,IDISP,TEMP,ISCODE)
      IF(IRC.NE.0) THEN
        WRITE(6,*)'RDMGRD: RDMGRD failed to read '//MINAME(JOB)
        WRITE(6,*)'  Displacement IDISP=',IDISP
        WRITE(6,*)'    Option code IOPT=',IOPT
        WRITE(6,*)'    Data label LABEL=',LABEL
        WRITE(6,*)'Symmetry code ISCODE=',ISCODE
        WRITE(6,*)'    Return code IRC =',IRC
        CALL ABEND()
      END IF

C Addressing integral blocks in ARRAY:
      ISUM=0
      DO IS=1,NSYM
       JS=MUL(IS,ISYMP)
       IAOFF(IS)=ISUM
       NBI=NBASF(IS)
       NBJ=NBASF(JS)
       NBIJ=NBI*NBJ
       ISUM=ISUM+NBIJ
      END DO
      IF(ISUM.GT.NARRAY) THEN
        WRITE(6,*)'RASSI/RDMGRD: Output ARRAY has insufficient length.'
        WRITE(6,*)' Input parameter NARRAY=',NARRAY
        WRITE(6,*)' Needed size       ISUM=',ISUM
        CALL ABEND()
      END IF
C Move buffer integrals into ARRAY in proper format:
      DO IS=1,NSYM
        NBI=NBASF(IS)
        IF(NBI.LE.0) GOTO 11
        IF(ISYMP.EQ.1) THEN
         LT=1+ITOFF(IS)
         IA1=1+IAOFF(IS)
         CALL SQUARE(TEMP(LT),ARRAY(IA1),1,NBI,NBI)
         IF(STYPE(1:4).EQ.'ANTI') THEN
          DO J=1,NBI-1
           DO I=J+1,NBI
            ARRAY(IA1-1+I+NBI*(J-1))=-ARRAY(IA1-1+J+NBI*(I-1))
           END DO
          END DO
         END IF
        ELSE
         JS=MUL(IS,ISYMP)
         IF(IS.LT.JS) GOTO 11
         NBJ=NBASF(JS)
         IF(NBJ.LE.0) GOTO 11
         LT=1+ITOFF(IS)
         IA1=1+IAOFF(IS)
         IA2=1+IAOFF(JS)
         CALL DCOPY_(NBI*NBJ,TEMP(LT),1,ARRAY(IA1),1)
         F=1.0D0
         IF(STYPE(1:4).EQ.'ANTI') F=-F
         DO I=1,NBI
          DO J=1,NBJ
           ARRAY(IA2-1+J+NBJ*(I-1))=F*ARRAY(IA1-1+I+NBI*(J-1))
          END DO
         END DO
        END IF
  11    CONTINUE
      END DO
C Get rid of temporary buffer
      CALL mma_deallocate(TEMP)

C Close MCKINT file:
      IRC=-1
      IOPT=0
      CALL CLSMCK(IRC,IOPT)
      IF(IRC.NE.0) THEN
        WRITE(6,*)'RASSI/RDMGRD: Failed to close '//MINAME(JOB)
        WRITE(6,*)'Unit nr LUMCK=',LUMCK
        WRITE(6,*)'Option code IOPT=',IOPT
        WRITE(6,*)'Return code IRC =',IRC
        CALL ABEND()
      END IF

      END SUBROUTINE RDMGRD
