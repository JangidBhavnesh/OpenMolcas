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
      SUBROUTINE MECTL(PROP,OVLP,HAM,ESHFT)
      use rassi_aux, only: ipglob
      use rassi_global_arrays, only: HDIAG
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NSTATE, NPROP, PRXVR, IFHAM, ToFile,
     &                 IfHDia, IfShft, PrMER, IfDCPL, iComp, IPUSED,
     &                 PNAME, PNUC, PORIG
      use cntrl, only: FnEig, LuEig

      IMPLICIT None
      REAL*8 PROP(NSTATE,NSTATE,NPROP),OVLP(NSTATE,NSTATE),
     &       HAM(NSTATE,NSTATE),ESHFT(NSTATE)

      Real*8, Allocatable:: DerCpl(:), NucChg(:)
      Integer nCol, iProp, i, ISTA, IFON, j, iState, iDisk, jState,
     &        NST, nAtom, IEND, IfHD
      REAL*8 X, PLIMIT, PMAX
*

C Print results:
      NCOL=4
      IF(PRXVR.and.IPGLOB.ge.0) THEN
      WRITE(6,*)
      Call CollapseOutput(1,'Expectation values for input states')
      WRITE(6,'(3X,A)')     '-----------------------------------'
      WRITE(6,*)
      WRITE(6,*)' EXPECTATION VALUES OF 1-ELECTRON OPERATORS'
      WRITE(6,*)' FOR THE RASSCF INPUT WAVE FUNCTIONS:'
      WRITE(6,*)
      DO IPROP=1,NPROP
       IF(IPUSED(IPROP).NE.0) THEN

* Skip printing if all the diagonal values are very small
*  (presumed zero for reasons of selection rules)
        PLIMIT=1.0D-10
        PMAX=0.0D0
        DO I=1,NSTATE
         PMAX=MAX(PMAX,ABS(PROP(I,I,IPROP)+PNUC(IPROP)*OVLP(I,I)))
        END DO
        IF(PMAX.GE.PLIMIT) THEN

        DO ISTA=1,NSTATE,NCOL
          IEND=MIN(NSTATE,ISTA+NCOL-1)
          WRITE(6,*)
          WRITE(6,'(1X,A,A8,A,I4)')
     &  'PROPERTY: ',PNAME(IPROP),'   COMPONENT:',ICOMP(IPROP)
          WRITE(6,'(1X,A,3ES17.8)')
     &'ORIGIN    :',(PORIG(I,IPROP),I=1,3)
          WRITE(6,'(1X,A,I8,3I17)')
     &'STATE     :',(I,I=ISTA,IEND)
          WRITE(6,*)
          WRITE(6,'(1X,A,4(1x,G16.9))')
     &'ELECTRONIC:',(PROP(I,I,IPROP),I=ISTA,IEND)
          WRITE(6,'(1X,A,4(1x,G16.9))')
     &'NUCLEAR   :',(PNUC(IPROP)*OVLP(I,I),I=ISTA,IEND)
          WRITE(6,'(1X,A,4(1x,G16.9))')
     &'TOTAL     :',(PROP(I,I,IPROP)+PNUC(IPROP)*OVLP(I,I),I=ISTA,IEND)
          WRITE(6,*)
        END DO

        END IF
       END IF
      END DO
      Call CollapseOutput(0,'Expectation values for input states')
      END IF

      IFON=1
      X=0.0D0
      DO I=2,NSTATE
       DO J=1,I-1
        X=MAX(X,ABS(OVLP(I,J)))
       END DO
      END DO
      IF (X.GE.1.0D-6) IFON=0
      IFHD=1
      X=0.0D0
      DO I=2,NSTATE
       DO J=1,I-1
        X=MAX(X,ABS(HAM(I,J)))
       END DO
      END DO
      IF (X.GE.1.0D-6) IFHD=0
      IF(IFON.eq.0) IFHD=0

      IF(IPGLOB.GE.2) THEN
       IF(IFHAM) THEN
         WRITE(6,*)
         WRITE(6,*)' HAMILTONIAN MATRIX FOR THE ORIGINAL STATES:'
         WRITE(6,*)
         IF(IFHD.eq.1) THEN
          WRITE(6,*)' Diagonal, with energies'
          WRITE(6,'(5(1X,F15.8))')(HAM(J,J),J=1,NSTATE)
         !do J=1,NSTATE
         !WRITE(6,*) HAM(J,J)
         !enddo
         ELSE
          DO ISTA=1,NSTATE,5
            IEND=MIN(ISTA+4,NSTATE)
            WRITE(6,'(10X,5(8X,A3,I4,A3))')
     &         (' | ', I, ' > ',I=ISTA,IEND)
            DO J=1,NSTATE
              WRITE(6,'(A3,I4,A3,5(2X,F16.8))')
     &           ' < ',J,' | ', (HAM(I,J),I=ISTA,IEND)
            END DO
          ENDDO
         END IF
       END IF
      END IF

      IF(IPGLOB.GE.2) THEN
        WRITE(6,*)
        WRITE(6,*)'     OVERLAP MATRIX FOR THE ORIGINAL STATES:'
        WRITE(6,*)
        IF(IFON.eq.1) THEN
         WRITE(6,*)' Diagonal, with elements'
         WRITE(6,'(5(1X,F15.8))')(OVLP(J,J),J=1,NSTATE)
        ELSE
         DO ISTATE=1,NSTATE
          WRITE(6,'(5(1X,F15.8))')(OVLP(ISTATE,J),J=1,ISTATE)
         END DO
        END IF
      END IF

* Addition by A.Ohrn. If ToFile keyword has been specified, we put
* numbers on the auxiliary rassi-to-qmstat file.
      If(ToFile) then
        If(.not.IfHam) then
          Write(6,*)
          Write(6,*)'You ask me to print hamiltonian, but there is '
     &//'none to print!'
          Call Abend()
        Endif
        Call DaName(LuEig,FnEig)
        iDisk=0
        Do iState=1,nState
          Do jState=1,iState
            Call dDaFile(LuEig,1,Ham(iState,jState),1,iDisk)
          Enddo
        Enddo
        Do iState=1,nState
          Do jState=1,iState
            Call dDaFile(LuEig,1,OvLp(iState,jState),1,iDisk)
          Enddo
        Enddo
*-- File is closed in eigctl.
        Call DaClos(LuEig)
      Endif
* End of addition by A.Ohrn.

      IF(IFHAM .AND. (IFHDIA.OR.IFSHFT)) THEN
        DO ISTATE=1,NSTATE
          IF(.NOT.IFSHFT) ESHFT(ISTATE)=0.0D0
          IF(IFHDIA) ESHFT(ISTATE)=ESHFT(ISTATE)+
     &              (HDIAG(ISTATE)-HAM(ISTATE,ISTATE))
        END DO
        DO ISTATE=1,NSTATE
         DO JSTATE=1,NSTATE
          HAM(ISTATE,JSTATE)=HAM(ISTATE,JSTATE)+
     &      0.5D0*(ESHFT(ISTATE)+ESHFT(JSTATE))*OVLP(ISTATE,JSTATE)
         END DO
        END DO
        IFHD=1
        DO I=2,NSTATE
         DO J=1,J-1
          IF(ABS(HAM(I,J)).GE.1.0D-10) IFHD=0
         END DO
        END DO
        IF(IFON.eq.0) IFHD=0
        IF(IPGLOB.GE.2) THEN
         WRITE(6,*)
         WRITE(6,*)' USER-MODIFIED HAMILTONIAN FOR THE ORIGINAL STATES:'
         WRITE(6,*)'(With user shifts, and/or replaced diagonal'
         WRITE(6,*)' elements, including overlap corrections.)'
         IF(IFHD.eq.1) THEN
          WRITE(6,*)' Diagonal, with energies'
          WRITE(6,'(5(1X,F15.8))')(HAM(J,J),J=1,NSTATE)
         !do J=1,NSTATE
         !WRITE(6,*) HAM(J,J)
         !enddo
         ELSE
          DO ISTATE=1,NSTATE
            WRITE(6,'(5(1X,F15.8))')(HAM(ISTATE,J),J=1,ISTATE)
          END DO
         END IF
        END IF
      END IF
CPAM00 End of updated HDIA/SHIFT section.

      IF(IPGLOB.GT.0 .and. PRMER) THEN
      WRITE(6,*)
      Call CollapseOutput(1,'Matrix elements for input states')
      WRITE(6,'(3X,A)')     '--------------------------------'
      WRITE(6,*)
      WRITE(6,*)' MATRIX ELEMENTS OF 1-ELECTRON OPERATORS'
      WRITE(6,*)' FOR THE RASSCF INPUT WAVE FUNCTIONS:'
      WRITE(6,*)
      DO IPROP=1,NPROP
       IF(IPUSED(IPROP).NE.0) THEN
        DO ISTA=1,NSTATE,NCOL
          IEND=MIN(NSTATE,ISTA+NCOL-1)
          WRITE(6,*)
          WRITE(6,'(1X,A,A8,A,I4)')
     &  'PROPERTY: ',PNAME(IPROP),'   COMPONENT:',ICOMP(IPROP)
          WRITE(6,'(1X,A,3ES17.8)')
     &'ORIGIN    :',(PORIG(I,IPROP),I=1,3)
          WRITE(6,'(1X,A,I8,3I17)')
     &'STATE     :',(I,I=ISTA,IEND)
          WRITE(6,*)
          DO J=1,NSTATE
            WRITE(6,'(1X,I5,6X,4(1x,G16.9))')
     & J,(PROP(J,I,IPROP)+PNUC(IPROP)*OVLP(J,I),I=ISTA,IEND)
          END DO
          WRITE(6,*)
        END DO
       END IF
      END DO
      Call CollapseOutput(0,'Matrix elements for input states')
      END IF
cnf
      If (IfDCpl) Then
         Call Get_iScalar('Unique atoms',natom)
         Call mma_allocate(NucChg,natom,Label='NucChg')
         Call Get_dArray('Effective nuclear Charge',NucChg,nAtom)
         nST = nState*(nState+1)/2
         Call mma_allocate(DerCpl,3*natom*nST,Label='DerCpl')
         Call AppDerCpl(natom,nST,NucChg,Prop,DerCpl,HAM)
         Call mma_deallocate(DerCpl)
         Call mma_deallocate(NucChg)
      End If
cnf
      CALL Put_dArray('SFS_HAM' ,HAM,NSTATE**2)
      CALL Put_dArray('SFS_OVLP',OVLP,NSTATE**2)

      WRITE(6,*)

      END SUBROUTINE MECTL
