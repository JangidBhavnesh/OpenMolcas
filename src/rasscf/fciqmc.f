************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2014, Giovanni Li Manni                                *
*               2019, Oskar Weser                                      *
************************************************************************
      module fciqmc
#ifdef _MOLCAS_MPP_
      use mpi
#endif
      use filesystem, only : chdir_, getcwd_, get_errno_, strerror_
      use fortran_strings, only : str
      use stdalloc, only : mma_allocate, mma_deallocate

      use rasscf_data, only : iter, lRoots, nRoots, S, KSDFT, NAC, EMY,
     &    rotmax, ener, iAdr15, iRoot, Weight, nacpr2, nacpar
      use general_data, only : nSym, iSpin, LuInta, nactel,
     &    jobIph, ntot, ntot1, ntot2,
     &    nFro, nIsh, nRs1, nRs2, nRs3, nDel, nAsh, nBas
      use gugx_data, only : IfCAS
      use gas_data, only : ngssh, iDoGas, nGAS, iGSOCCX

      use fcidump_reorder, only : get_P_GAS, get_P_inp,ReOrFlag,ReOrInp
      use fcidump, only : make_fcidumps, transform
      use fciqmc_make_inp, only : make_inp
      use fciqmc_read_RDM, only : read_neci_RDM
      implicit none
      private
      public :: fciqmc_ctl, DoNECI, DoEmbdNECI, cleanup
      logical ::
     &  DoEmbdNECI = .false.,
     &  DoNECI = .false.
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      integer*4 :: error
#endif
      interface
        integer function isfreeunit(iseed)
          integer, intent(in) :: iseed
        end function
      end interface
      save
      contains

!>  @brief
!>    Start and control FCIQMC.
!>
!>  @author Giovanni Li Manni, Oskar Weser
!>
!>  @details
!>  For meaning of global variables NTOT1, NTOT2, NACPAR
!>  and NACPR2, see src/Include/general.inc and src/Include/rasscf.inc.
!>  This routine will replace CICTL in FCIQMC regime.
!>  Density matrices are generated via double-run procedure in NECI.
!>  They are then dumped on arrays DMAT, DSPN, PSMAT, PAMAT to replace
!>  what normally would be done in DavCtl if NECI is not used.
!>  F_In is still generated by SGFCIN... in input contains
!>  only two-electron terms as computed in TRA_CTL2.
!>  In output it contains also the one-electron contribution
!>
!>  @paramin[in] CMO MO coefficients
!>  @paramin[in] DIAF DIAGONAL of Fock matrix useful for NECI
!>  @paramin[in] D1I_MO Inactive 1-dens matrix
!>  @paramin[in] TUVX Active 2-el integrals
!>  @paramin[inout] F_In Fock matrix from inactive density
!>  @paramin[inout] D1S_MO Average spin 1-dens matrix
!>  @paramin[out] DMAT Average 1 body density matrix
!>  @paramin[out] PSMAT Average symm. 2-dens matrix
!>  @paramin[out] PAMAT Average antisymm. 2-dens matrix
!>  @paramin[in] fake_run  If true the NECI run is not performed, but
!>    the RDMs are read from previous runs.
      subroutine fciqmc_ctl(CMO, DIAF, D1I_AO, D1A_AO, TUVX, F_IN,
     &                      D1S_MO, DMAT, PSMAT, PAMAT,
     &                      fake_run)
      use write_orbital_files, only : get_typeidx
      implicit none
#include "output_ras.fh"
#include "rctfld.fh"
#include "timers.fh"
      real*8, intent(in) ::
     &    CMO(nTot2), DIAF(nTot),
     &    D1I_AO(nTot2), D1A_AO(nTot2), TUVX(nAcpr2)
      real*8, intent(inout) :: F_In(nTot1), D1S_MO(nAcPar)
      real*8, intent(out) :: DMAT(nAcpar),
     &    PSMAT(nAcpr2), PAMAT(nAcpr2)
      logical, intent(in), optional :: fake_run
      logical :: fake_run_
      real*8, save :: NECIen
      integer :: iPRLEV, iOff, iSym, iBas, i, j, jRoot,
     &    permutation(sum(nAsh(:nSym))), file_id, typeidx(7, 8)
      character(*), parameter :: filename = 'ORTORB'
      character(len=80) :: orbfile_title
      real*8 :: orbital_E(nTot), folded_Fock(nAcPar), vDummy(2)

      parameter(ROUTINE = 'FCIQMC_clt')

      call qEnter(routine)

      fake_run_ = merge(fake_run, .false., present(fake_run))

! Local print level (if any)
      iprlev = iprloc(1)
      if(iprlev.ge.debug) then
        write(lf,*)
        write(lf,*) ' ===================='
        write(lf,*) ' Entering FCIQMC_Ctl'
        write(lf,*) ' ===================='
        write(lf,*)
        write(lf,*) ' iteration count =', ITER
        write(lf,*) ' IFCAS value     =', IFCAS
        write(lf,*) ' lroots,nroots   =', lroots,nroots
        write(lf,*)
      end if
! set up flag 'IFCAS' for GAS option, which is set up in gugatcl originally.
! IFCAS = 0: This is a CAS calculation
! IFCAS = 1: This is a RAS calculation
! IFCAS = 2: This is a GAS calculation
      if(iprlev.ge.debug) then
        write(lf,*)
        write(lf,*) ' CMO in FCIQMC_CTL'
        write(lf,*) ' ---------------------'
        write(lf,*)
        ioff=1
        do isym = 1,nsym
          ibas = nbas(isym)
          if(ibas.ne.0) then
            write(6,*) 'Sym =', isym
            do i= 1,ibas
              write(6,*) (cmo(ioff+ibas*(i-1)+j),j=0,ibas-1)
            end do
            ioff = ioff + (ibas*ibas)
          end if
        end do
      end if

! SOME DIRTY SETUPS
! TODO(Giovanni): No dirty setups
      S = 0.5D0 * DBLE(ISPIN-1)

      call check_options(lRoots, lRf, KSDFT, iDoGAS, iGSOCCX, nGAS)

! Produce a working FCIDUMP file
! TODO: permutation has to be applied at more places
      select case (ReOrFlag)
        case (2:)
          permutation = get_P_inp(ReOrInp)
        case (-1)
          permutation = get_P_GAS(nGSSH)
      end select

! This call is not side effect free, sets EMY and modifies F_IN
      call transform(iter, CMO, DIAF, D1I_AO, D1A_AO, D1S_MO,
     &      F_IN, orbital_E, folded_Fock)

      if (ReOrFlag /= 0) then
        call make_fcidumps(orbital_E, folded_Fock, TUVX, EMY,
     &                     permutation)
      else
        call make_fcidumps(orbital_E, folded_Fock, TUVX, EMY)
      end if

      file_id = 50
      file_id = isfreeunit(file_id)
      typeidx = get_typeidx(
     &    nFro, nIsh, nRs1, nRs2, nRs3, nBas, nGSSH, nDel)
      orbfile_title = 'Orbitals that are used for FCIQMC.'
      Call WrVec(filename, file_id, 'CIE', nSym, nBas, nBas,
     &           CMO, vDummy, orbital_E, typeidx, orbfile_title)

! Run NECI
      call Timing(Rado_1, Swatch, Swatch, Swatch)
#ifdef _MOLCAS_MPP_
      if (is_real_par()) call MPI_Barrier(MPI_COMM_WORLD, error)
#endif
      if (.not. fake_run_) call run_neci(DoEmbdNECI, NECIen)
! NECIen so far is only the energy for the GS.
! Next step it will be an array containing energies for all the optimized states.
      do jRoot = 1, lRoots
        ENER(jRoot, ITER) = NECIen
      end do


      call get_neci_RDM(D1S_MO, DMAT, PSMAT, PAMAT)

! print matrices
      if (IPRLEV >= DEBUG) then
        call TRIPRT('Averaged one-body density matrix, DMAT',
     &              ' ',DMAT,NAC)
        call TRIPRT('Averaged one-body spin density matrix, DS',
     &              ' ',D1S_MO,NAC)
        call TRIPRT('Averaged two-body density matrix, P',
     &              ' ',PSMAT,NACPAR)
        call TRIPRT('Averaged antisymmetric two-body density matrix,PA',
     &              ' ',PAMAT,NACPAR)
      end if

      if (nAsh(1) /= nac) call dblock(dmat)


      call Timing(Rado_2, Swatch, Swatch, Swatch)
      Rado_2 = Rado_2 - Rado_1
      Rado_3 = Rado_3 + Rado_2

      call qExit(routine)
      end subroutine fciqmc_ctl


      subroutine run_neci(DoEmbdNECI, NECien)
        implicit none
        logical, intent(in) :: DoEmbdNECI
        real*8, intent(out) :: NECIen
        if (DoEmbdNECI) then
          call make_inp(readpops=iter >= 5 .and. abs(rotmax) < 1d-2)
#ifdef _NECI_
          write(6,*) 'NECI called automatically within Molcas!'
          if (myrank /= 0) call chdir_('..')
          call necimain(NECIen)
          if (myrank /= 0) call chdir_('tmp_'//str(myrank))
#else
          call WarningMessage(2, 'EmbdNECI is given in input, '//
     &'so the embedded NECI should be used. Unfortunately MOLCAS was '//
     &'not compiled with embedded NECI. Please use -DNECI=ON '//
     &'for compiling or use an external NECI.')
#endif
        else
          call make_inp()
          if (myrank == 0) then
            call write_ExNECI_message()
            call wait_and_read(NECIen)
            write(6,*) 'I read the following energy:', NECIen
          end if
#ifdef _MOLCAS_MPP_
          call MPI_Bcast(NECIen, 1, MPI_REAL8, 0,MPI_COMM_WORLD,error)
#endif
        end if
      end subroutine run_neci


      subroutine wait_and_read(NECIen)
        implicit none
        real*8, intent(out) :: NECIen
        logical :: newcycle_found
        integer :: LuNewC
        newcycle_found = .false.
        do while(.not. newcycle_found)
          call sleep(1)
          call f_Inquire('NEWCYCLE', newcycle_found)
        end do
        write(6, *) 'NEWCYCLE file found. Proceding with SuperCI'
        LuNewC = 12
        call molcas_open(LuNewC, 'NEWCYCLE')
          read(LuNewC,*) NECIen
        close(LuNewC, status='delete')
      end subroutine wait_and_read


      subroutine abort_(message)
        implicit none
        character(*), intent(in) :: message
        call WarningMessage(2, message)
        call QTrace()
        call Abend()
      end subroutine


      subroutine cleanup()
        use fciqmc_make_inp, only : make_inp_cleanup => cleanup
        use fciqmc_read_RDM, only : read_RDM_cleanup => cleanup
        use fcidump, only : fcidump_cleanup => cleanup
        implicit none
        call make_inp_cleanup()
        call read_RDM_cleanup()
        call fcidump_cleanup()
      end subroutine cleanup


      subroutine check_options(lroots, lRf, KSDFT,
     &      DoGAS, iGSOCCX, nGAS)
        implicit none
        integer, intent(in) :: lroots, iGSOCCX(:, :),nGAS
        logical, intent(in) :: lRf, DoGAS
        character(*), intent(in) :: KSDFT
        logical :: Do_ESPF
        if (lroots > 1) then
          call abort_('FCIQMC does not support State Average yet!')
        end if
        call DecideOnESPF(Do_ESPF)
        if ( lRf .or. KSDFT /= 'SCF' .or. Do_ESPF) then
          call abort_('FCIQMC does not support Reaction Field yet!')
        end if
        if (DoGAS) then
          if (.not. all(iGSOCCX(:nGAS, 1) == iGSOCCX(:nGAS, 2))) then
            call abort_('Only disconnected GAS spaces are '//
     &        'currently supported in FCIQMC.')
          end if
        end if
      end subroutine check_options


      subroutine write_ExNECI_message()
        implicit none
        character(1024) :: h5fcidmp, fcidmp, fcinp, newcycle, WorkDir
        integer :: L, err

        call getcwd_(WorkDir, err)
        if (err /= 0) write(6, *) strerror_(get_errno_())
        call prgmtranslate_master('H5FCIDMP', h5fcidmp, L)
        call prgmtranslate_master('FCIDMP', fcidmp, L)
        call prgmtranslate_master('FCINP', fcinp, L)
        call prgmtranslate_master('NEWCYCLE', newcycle, L)


        write(6,'(A)')'Run NECI externally.'
        write(6,'(A)')'Get the (example) NECI input:'
        write(6,'(4x, A, 1x, A, 1x,  A)')
     &    'cp', trim(fcinp), '$NECI_RUN_DIR'
        write(6,'(A)')'Get the ASCII formatted FCIDUMP:'
        write(6,'(4x, A, 1x, A, 1x,  A)')
     &    'cp', trim(fcidmp), '$NECI_RUN_DIR/FCIDUMP'
        write(6,'(A)')'Or the HDF5 FCIDUMP:'
        write(6,'(4x, A, 1x, A, 1x,  A)')
     &    'cp', trim(h5fcidmp), '$NECI_RUN_DIR'
        write(6, *)
        write(6,'(A)') "When finished do:"
        write(6,'(4x, A)')
     &    'cp TwoRDM_aaaa.1 TwoRDM_abab.1 TwoRDM_abba.1 '//
     &    'TwoRDM_bbbb.1 TwoRDM_baba.1 TwoRDM_baab.1 '//trim(WorkDir)
        write(6,'(4x, A)')'echo $your_RDM_Energy > '//trim(newcycle)
        call xflush(6)
      end subroutine write_ExNECI_message

!> Generate density matrices for Molcas
!>   Neci density matrices are stored in Files TwoRDM_**** (in spacial orbital basis).
!>   I will be reading them from those formatted files for the time being.
!>   Next it will be nice if NECI prints them out already in Molcas format.
      subroutine get_neci_RDM(D1S_MO, DMAT, PSMAT, PAMAT)
        real*8, intent(out) :: D1S_MO(nAcPar), DMAT(nAcpar),
     &      PSMAT(nAcpr2), PAMAT(nAcpr2)
        real*8, allocatable ::
!> ONE-BODY DENSITY
     &    DTMP(:),
!> SYMMETRIC TWO-BODY DENSITY
     &    Ptmp(:),
!> ANTISYMMETRIC TWO-BODY DENSITY
     &    PAtmp(:),
!> ONE-BODY SPIN DENSITY
     &    DStmp(:)
        real*8 :: Scal
        integer :: jRoot, kRoot, iDisk, jDisk

        call mma_allocate(DTMP, nAcPar, label='Dtmp ')
        call mma_allocate(DStmp, nAcPar, label='DStmp')
        call mma_allocate(Ptmp, nAcPr2, label='Ptmp ')
        call mma_allocate(PAtmp, nAcPr2, label='PAtmp')

        call read_neci_RDM(DTMP, DStmp, Ptmp, PAtmp)

! COMPUTE AVERAGE DENSITY MATRICES
        do jRoot = 1, lRoots
          Scal = 0.0d0
          do kRoot = 1, nRoots
            if (iRoot(kRoot) == jRoot) Scal = Weight(kRoot)
          end do
          DMAT(:) = SCAL * DTMP(:)
          D1S_MO(:) = SCAL * PSMAT(:)
          PSMAT(:) = SCAL * Ptmp(:)
          PAMAT(:) = SCAL * PAtmp(:)
! Put it on the RUNFILE
          call Put_D1MO(DTMP,NACPAR)
          call Put_P2MO(Ptmp,NACPR2)
! Save density matrices on disk
          iDisk = IADR15(4)
          jDisk = IADR15(3)
          call DDafile(JOBIPH, 1, DTMP, NACPAR, jDisk)
          call DDafile(JOBIPH, 1, DStmp, NACPAR, jDisk)
          call DDafile(JOBIPH, 1, Ptmp, NACPR2, jDisk)
          call DDafile(JOBIPH, 1, PAtmp, NACPR2, jDisk)
        end do

        call mma_deallocate(DTMP)
        call mma_deallocate(DStmp)
        call mma_deallocate(Ptmp)
        call mma_deallocate(PAtmp)
      end subroutine get_neci_RDM

      function get_typeidx(
     &  nFro, nIsh, nRs1, nRs2, nRs3, nBas, nGSSH, nDel) result(typeidx)
        implicit none
        integer, intent(in) ::  nFro(:), nIsh(:), nRs1(:), nRs2(:),
     &    nRs3(:), nBas(:), nGSSH(:, :), nDel(:)
        integer :: typeidx(7, 8)

        typeidx(:, :nSym) = 0

        typeidx(1, :nSym) = nFro(:nSym)
        typeidx(2, :nSym) = nIsh(:nSym)
        if (.not. iDoGAS) then
          typeidx(3, :nSym) = nRS1(:nSym)
          typeidx(4, :nSym) = nRS2(:nSym)
          typeidx(5, :nSym) = nRS3(:nSym)
        else
          typeidx(4, :nSym) = sum(nGssh(1:nGAS, :nSym), dim=1)
        end if
        typeidx(7, :nSym) = nDel(:nSym)
        typeidx(6, :nSym) = nBas(:nSym) - sum(typeidx(:, :nSym), dim=1)
      end function
      end module fciqmc
