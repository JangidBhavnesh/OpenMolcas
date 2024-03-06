      SUBROUTINE mkism_cp2()

      use fciqmc_interface, only: DoFCIQMC
      use stdalloc, only: mma_allocate
      use gugx, only: SGS, L2ACT, LEVEL

      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"

#include "SysDef.fh"
      Integer nLev

      INTEGER I,IT,ITABS,ILEV,ISYM, iq

      NLEV=NASHT
      SGS%nLev = NLEV
      Call mma_allocate(SGS%ISM,NLEV,Label='ISM')
C ISM(LEV) IS SYMMETRY LABEL OF ACTIVE ORBITAL AT LEVEL LEV.
C PAM060612: With true RAS space, the orbitals must be ordered
C first by RAS type, then by symmetry.
      ITABS=0
      DO ISYM=1,NSYM
        DO IT=1,NASH(ISYM)
          ITABS=ITABS+1
! Quan: Bug in LEVEL(ITABS) and L2ACT
          if (DoCumulant .or. DoFCIQMC) then
             do iq=1,NLEV
               LEVEL(iq)=iq
               L2ACT(iq)=iq
             enddo
          endif
          ILEV=LEVEL(ITABS)
          SGS%ISM(ILEV)=ISYM
        END DO
      END DO

      END SUBROUTINE mkism_cp2