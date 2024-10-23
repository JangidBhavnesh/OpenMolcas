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
#include "xrhs.fh"
      SUBROUTINE MODRHS(IVEC,FIMO,NFIMO)
      USE SUPERINDEX
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      INTEGER IVEC, NFIMO
      REAL*8 FIMO(NFIMO)

      REAL*8, ALLOCATABLE:: WA(:), WC(:), WD(:)


***************************************************************
* Case A:
      ICASE=1
      IFOFF=0
      DO ISYM=1,NSYM
       IF(NINDEP(ISYM,1).EQ.0) GOTO 100

       NAS=NTUV(ISYM)
       NIS=NISH(ISYM)
       NWA=NAS*NIS
       IF(NWA.EQ.0) GOTO 100
       CALL mma_allocate(WA,NWA,Label='WA')
       CALL RHS_ALLO (NAS,NIS,lg_A)
C Read W from disk:
       CALL RHS_READ (NAS,NIS,lg_A,ICASE,ISYM,IVEC)
       CALL RHS_GET (NAS,NIS,lg_A,WA)
* Insert one-electron contribution to coupling <A|0>:
* WA(tvv,j)=FIMO(t,j)/NACTEL (+two-electron part)
       ISYJ=ISYM
       ISYT=ISYM
       NAT=NASH(ISYT)
       NIT=NISH(ISYT)
       DO IT=1,NAT
        ITTOT=NIT+IT
        ITABS=NAES(ISYT)+IT
        DO IJ=1,NIT
         VALUE=FIMO(IFOFF+(ITTOT*(ITTOT-1))/2+IJ)/DBLE(MAX(1,NACTEL))
         DO IVABS=1,NASHT
          IW1=KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
          IW2=IJ
          WA(IW1+NAS*(IW2-1))=WA(IW1+NAS*(IW2-1))+VALUE
         END DO
        END DO
       END DO
       CALL RHS_PUT (NAS,NIS,lg_A,WA)
C Put W on disk:
       CALL RHS_SAVE (NAS,NIS,lg_A,ICASE,ISYM,IVEC)
       CALL RHS_FREE (lg_A)
       CALL mma_deallocate(WA)

 100   CONTINUE
* End of loop over ISYM.
       NO=NORB(ISYM)
       IFOFF=IFOFF+(NO*(NO+1))/2
      END DO

***************************************************************
* Case C:
      ICASE=4
      IFOFF=0
      DO ISYM=1,NSYM
       IF(NINDEP(ISYM,4).EQ.0) GOTO 200
       NAS=NTUV(ISYM)
       NIS=NSSH(ISYM)
       NWC=NAS*NIS
       IF(NWC.EQ.0) GOTO 200
       CALL mma_allocate(WC,NWC,LABEL='WC')
       CALL RHS_ALLO (NAS,NIS,lg_C)
C Read W from disk:
       CALL RHS_READ (NAS,NIS,lg_C,ICASE,ISYM,IVEC)
       CALL RHS_GET (NAS,NIS,lg_C,WC)
* Insert one-electron contribution to coupling <C|0>:
* WC(xuu,a)=(FIMO(a,x)-sum((ay,yx), y=1,NASHT) )/NACTEL (+ two-el part)
       NIX=NISH(ISYM)
       NAX=NASH(ISYM)
       NSX=NSSH(ISYM)
       DO IX=1,NAX
        IXTOT=NIX+IX
        IXABS=NAES(ISYM)+IX
        DO IA=1,NSX
         IATOT=NIX+NAX+IA
         SUM=FIMO(IFOFF+(IATOT*(IATOT-1))/2+IXTOT)
         DO IYABS=1,NASHT
          IYYW=KTUV(IYABS,IYABS,IXABS)-NTUVES(ISYM)
          IYYWA=IYYW+NAS*(IA-1)
          SUM=SUM-WC(IYYWA)
         END DO
         ONEADD=SUM/DBLE(MAX(1,NACTEL))
         DO IUABS=1,NASHT
          IW1=KTUV(IXABS,IUABS,IUABS)-NTUVES(ISYM)
          IW2=IA
          WC(IW1+NAS*(IW2-1))= WC(IW1+NAS*(IW2-1))+ONEADD
         END DO
        END DO
       END DO
       CALL RHS_PUT (NAS,NIS,lg_C,WC)
C Put W on disk:
       CALL RHS_SAVE (NAS,NIS,lg_C,ICASE,ISYM,IVEC)
       CALL RHS_FREE (lg_C)
       CALL mma_deallocate(WC)

 200   CONTINUE
* End of loop over ISYM.
       NO=NORB(ISYM)
       IFOFF=IFOFF+(NO*(NO+1))/2
      END DO

***************************************************************
* Case D1:
      ICASE=5
      ISYM=1
      IF(NINDEP(ISYM,5).EQ.0) RETURN

      NAS=NASUP(ISYM,5)
      NIS=NISUP(ISYM,5)
      NWD=NAS*NIS
      IF(NWD.EQ.0) RETURN
      CALL mma_allocate(WD,NWD,LABEL='WD')
      CALL RHS_ALLO (NAS,NIS,lg_D)
C Read W from disk:
      CALL RHS_READ (NAS,NIS,lg_D,ICASE,ISYM,IVEC)
      CALL RHS_GET (NAS,NIS,lg_D,WD)

* Insert one-electron contribution to coupling <D1|0>:
* Compute WD1(vv,aj)=FIMO(a,j)/NACTEL (+ two-el part)
      IFOFF=0
      IAJ=0
      DO ISYJ=1,NSYM
       NIJ=NISH(ISYJ)
       NAJ=NASH(ISYJ)
       NSJ=NSSH(ISYJ)
       DO IA=1,NSJ
        IATOT=NIJ+NAJ+IA
        DO IJ=1,NIJ
         ONEADD=FIMO(IFOFF+(IATOT*(IATOT-1))/2+IJ)/DBLE(MAX(1,NACTEL))
         IAJ=IAJ+1
         DO ISYU=1,NSYM
          DO IU=1,NASH(ISYU)
           IUABS=NAES(ISYU)+IU
           IUU=KTU(IUABS,IUABS)-NTUES(ISYM)
           IWD=IUU+NAS*(IAJ-1)
           WD(IWD)=WD(IWD)+ONEADD
          END DO
         END DO
        END DO
       END DO
       NO=NORB(ISYJ)
       IFOFF=IFOFF+(NO*(NO+1))/2
      END DO
      CALL RHS_PUT (NAS,NIS,lg_D,WD)

C Put W on disk:
      CALL RHS_SAVE (NAS,NIS,lg_D,ICASE,ISYM,IVEC)
      CALL RHS_FREE (lg_D)
      CALL mma_deallocate(WD)

      END SUBROUTINE MODRHS
