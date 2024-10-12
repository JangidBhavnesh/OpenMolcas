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
      SUBROUTINE PT2INI()
      USE INPUTDATA, ONLY: INPUT, READIN_CASPT2
      USE REFWFN, ONLY: REFWFN_INIT, REFWFN_INFO, REFWFN_DATA,
     &                  REFWFN_CLOSE
      USE PT2WFN
      use caspt2_gradient, only: do_grad
      use caspt2_data, only: FIMO, FAMO, FIFA, HONE, DREF, PREF, DMIX,
     &                       DWGT
      use stdalloc, only: mma_allocate
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "intgrl.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "compiler_features.h"

      INTEGER LuSpool
C     Cholesky
      Integer iSym, iRC

*
* Probe the RunFile for some basic information
*
      Call Get_cArray('Seward Title',Header,144)
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      Call Get_iScalar('Unique atoms',nUniqAt)
      nbast=sum(nbas(1:nsym))
      nbsqt=sum(nbas(1:nsym)**2)
      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nbast)
      Call DecideOnCholesky(IfChol)
* PAM Feb 2008: The following statement was moved here from
* prpctl, in case prpctl call is bypassed by keyword ''NOPROP''.
      Call Put_cArray('Relax Method','CASPT2  ',8)
*
*     Allocate Input struct, read and process the Input
*
      ALLOCATE(Input)
      LuSpool=21
      Call SpoolInp(LuSpool)
      CALL READIN_CASPT2(LuSpool,nSym)
      Call Close_LuSpool(LuSpool)
* Initialize scratch files.
      CALL OPNFLS_CASPT2()
* Initialize the reference wavefunction file and read basic info
      Call refwfn_init(Input%FILE)
      Call refwfn_info()
* the input processing needs some data from the reference wavefunction
* to be set, hence this phase occurs after reading that data.
      Call ProcInp_Caspt2()
*
* Allocate some arrays that will stay globally allocated:
*

* Finally read the MO and CI data from the refwfn file, and close it as
* we have no more need for it and the same filename might be reused for
* the pt2wfn file. We do this after input processing because we need to
* know which roots to pick up. The MOs are stored on LUONEM, at address
* IAD1M(1), and the CI arrays on LUCIEX at IDCIEX.
      Call refwfn_data()
      Call refwfn_close()

* If necessary, make modifications to the inactive/virtual orbitals. The
* new MOs, if any, will overwrite the originals on LUONEM at IAD1M(1).
      If (Input%modify_correlating_MOs) Then
        Call correlating_orbitals()
      End If

* After possible reconfiguration of inactives/virtuals, the computation
* of total sizes of orbital spaces is done here, before any other code
* (that might rely on these to be correctly set!)
      Call wfnsizes()

* Create the PT2 wavefunction file (formerly JOBMIX). The reference file
* should not be active, as it might be the same file (in which case it
* is overwritten).
      call pt2wfn_init()
      call pt2wfn_data()
*
* Global allocations active throughout the program (NOTRI was computed
* during the call to wfnsizes, so this is done here.)
*
* The total fock matrix (sum of inactive and active contrib.)
      NFIFA=NOTRI
      CALL mma_allocate(FIFA,NFIFA,Label='FIFA')
* The one-electron Hamiltonian
      NHONE=NOTRI
      CALL mma_allocate(HONE,NHONE,Label='HONE')
* The fock matrix with contributions from inactive orbitals, only.
      NFIMO=NOTRI
      CALL mma_allocate(FIMO,NFIMO,Label='FIMO')
* The fock matrix with contributions from active orbitals, only.
      NFAMO=NOTRI
      CALL mma_allocate(FAMO,NFAMO,Label='FAMO')
* Density matrices, active indices.
      CALL mma_allocate(DREF,NDREF,Label='DREF')
      CALL mma_allocate(PREF,NPREF,Label='PREF')
* All the density matrices kept in memory
      CALL mma_allocate(DMIX,NDREF,NSTATE,Label='DMIX')
      CALL mma_allocate(DWGT,NSTATE,NSTATE,Label='DWGT')

* Print input data
      CALL PRINP_CASPT2()

C INITIALIZE SPLIT-GRAPH UGA TABLES.
      CALL POLY0()

C Initialize superindex tables. Check memory needs.
      CALL SIZES()

C Initialize sizes, offsets etc used in equation solver.
      CALL EQCTL1()

      If (IfChol) then
* initialize Cholesky information
        Call Cho_X_init(irc,0.0d0)
        if (irc.ne.0) then
          write(6,*) 'CASPT2: Non-zero rc in Cho_X_init'
          CALL QUIT(irc)
        endif
* import_ch transfers some values from Cholesky module to chocaspt2.fh
        call import_cho(numcho_pt2,infvec_n2_pt2,maxvec_pt2)
        Call setup_cho(nSym,nIsh,nAsh,nSsh,NumCho_pt2,'Allo')
* get unit numbers for Cholesky MO vectors
        Do iSym=1,nSym
          Call Cho_Caspt2_OpenF(0,1,iSym,nIsplit(iSym))
          Call Cho_Caspt2_OpenF(0,2,iSym,nAsplit(iSym))
        End Do
* set up size info for cholesky vector transformation in tracho2
        call trachosz()
      End If

* Allocate global orbital arrays:
      CALL GETMEM('LCMOPT2','ALLO','REAL',LCMOPT2,NCMO)
* Allocate global orbital transformation arrays:
      CALL GETMEM('TORB','ALLO','REAL',LTORB,NTORB)
      CALL GETMEM('TAT','ALLO','REAL',LTAT,NTAT)

! initialize quantities for gradient calculation
      If (do_grad) Then
        CALL GrdIni()
      End If

      END SUBROUTINE PT2INI

      SUBROUTINE PT2CLS()
      USE SUPERINDEX
      USE INPUTDATA, ONLY: CLEANUP_INPUT
      USE PT2WFN
      use gugx, only: SGS, CIS, EXS
      use caspt2_data, only: FIMO, FAMO, FIFA, HONE, DREF, PREF, DMIX,
     &                       DWGT
      use stdalloc, only: mma_deallocate
* NOT TESTED
#if 0
      use OFembed, only: FMaux
#endif
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"

      Integer iSym
C     Cholesky return code
      INTEGER irc
C     size of idsct array
      INTEGER NIDSCT

      If (IfChol) then
*---  Finalize Cholesky information
         Call Cho_X_Final(irc)
         if (irc.ne.0) then
            write(6,*) 'CASPT2: Non-zero rc in Cho_X_Final'
            CALL QUIT(irc)
         endif
*     Delete files of MO Cholesky vectors.
         Do iSym=1,nSym
            Call Cho_Caspt2_OpenF(3,1,iSym,nIsplit(iSym)) !close+delete
            Call Cho_Caspt2_OpenF(3,2,iSym,nAsplit(iSym)) !close+delete
         End Do
*     Deallocate memory used as index arrays in the setup section
         Call setup_cho(nSym,nIsh,nAsh,nSsh,NumCho_pt2,'Free')
* NOT TESTED
#if 0
         Call mma_deallocate(FMaux,safe='*')
#endif
         ! deallocate chovec_io arrays
         call trachosz_free()
      End If

* Deallocate SGUGA tables:
      CALL mkGUGA_Free(SGS,CIS,EXS)

C     Deallocate MAGEB, etc, superindex tables:
      CALL SUPFREE()
* Deallocate global array for Fock matrix, etc:
      CALL mma_deallocate(FIFA)
      CALL mma_deallocate(HONE)
      CALL mma_deallocate(FIMO)
      CALL mma_deallocate(FAMO)
      CALL mma_deallocate(DREF)
      CALL mma_deallocate(PREF)
      CALL mma_deallocate(DMIX)
      CALL mma_deallocate(DWGT)
* Deallocate global orbital transformation arrays:
      CALL GETMEM('TORB','FREE','REAL',LTORB,NTORB)
      CALL GETMEM('TAT','FREE','REAL',LTAT,NTAT)
* Deallocate global orbital arrays:
      CALL GETMEM('LCMOPT2','FREE','REAL',LCMOPT2,NCMO)
* Deallocate global RHS disk offsets (allocated in eqctl1):
      NIDSCT=MXSCT*8*MXCASE*MXVEC
      CALL GETMEM('IDSCT','FREE','INTE',LIDSCT,NIDSCT)

      call pt2wfn_close()
C     Close all files:
      CALL CLSFLS_CASPT2()

C free input struct
      CALL CleanUp_Input()
      End SUBROUTINE PT2CLS
