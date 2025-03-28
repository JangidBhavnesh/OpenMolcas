!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
Module Lucia_Data
Implicit None
#include "Molcas.fh"
Private :: MaxBfn,MaxBfn_Aux,MxAO,mxAtom,mxroot,mxNemoAtom,Mxdbsc,lCache,mxact,mxina,mxbas,mxOrb,mxSym,mxGAS, &
           LENIN,LENIN1,LENIN2,LENIN3,LENIN4,LENIN5,LENIN6,LENIN8

! Stuff from mxpdim.fh
Integer, PARAMETER:: MXPIRR = 20
Integer, PARAMETER:: MXPOBS = 20
Integer, PARAMETER:: MXPR4T = 10
Integer, PARAMETER:: MXPORB = 500
Integer, PARAMETER:: MXPICI = 30
Integer, PARAMETER:: MXPSTT = 2500
Integer, PARAMETER:: MXPCSM = 100
!. Note : MXPNGAS = MXPR4T+6 !!
!. Required in order to handle GAS and RAS within /LUCINP/
Integer, PARAMETER:: MXPNGAS = MXPR4T + 6
Integer, PARAMETER:: MXPNSMST = 8
!. Largest allowed division of space for perturbation operator
Integer, PARAMETER:: MXPPTSPC=20
!
!     PARAMETER(MXPTSOB = 100) ! introduced by Kamal Sharkas in commit 91722faa6
!                              ! probably for some large-scale cases
!                              ! changed back to default - stknecht/dec 2015
Integer, PARAMETER:: MXPTSOB = 35


! Stuff from strinp.fh
INTEGER       NSTTYP,NELEC(MXPSTT)

! Stuff from stinf.fh
Integer      ISTAC(MXPSTT,2),NOCTYP(MXPSTT)

! Stuff from orbinp.fh
Integer       NINOB,NACOB,NDEOB,NOCOB,NTOOB,                      &
              NORB1,NORB2,NORB3,                                  &
              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),               &
              NINOBS(MXPOBS),                                     &
              NACOBS(MXPOBS),NOCOBS(MXPOBS),                      &
              NTOOBS(MXPOBS),NDEOBS(MXPOBS),                      &
              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),       &
              ITPFSO(MXPORB),IBSO(MXPOBS),                        &
              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),     &
              NOBPTS(6+MXPR4T,MXPOBS),IOBPTS(6+MXPR4T,MXPOBS),    &
              ITOOBS(MXPOBS),ITPFTO(MXPORB),ISMFSO(MXPORB),       &
              NOBPT(6+MXPR4T),NAOS_ENV(MXPOBS),NMOS_ENV(MXPOBS),  &
              I_IAD(MXPNGAS),I_IADX(MXPNGAS),MXTSOB,MXTOB,        &
              MXTSOB_P,MXTSOB_H

! Stuff from oper.fh
INTEGER     I12,IPERTOP,IAPR,IPART,I_RES_AB

! Stuff from lucinp.fh
Integer       PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,                    &
              INTSPC,EXTSPC,NRSSH(MXPIRR,3),                      &
              MNRS1R,MXRS3R,NACTEL,                               &
              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPIRR,MXPR4T),       &
              MXR4TP, MXHR0,MXER4,                                &
              NINASH(MXPIRR),                                     &
              NDELSH(MXPIRR),MNRS10,MXRS30,                       &
              MNRS1RE,MXRS3RE,                                    &
              MNRS1ZE,MXRS3ZE

! Stuff from ssave.fh
INTEGER NELIS(4), NSTRKS(4)

! Stuff from loff.fh
INTEGER, PARAMETER :: LOFFI=8**6 !SJS

! Stuff from irat.fh
INTEGER IRAT

! Stuff from io_util.fh
INTEGER IDISK(100)

! Stuff from intform.fh
INTEGER IH1FORM,IH2FORM

! Stuff from gasstr.fh
INTEGER       MNGSOC(MXPNGAS),MXGSOC(MXPNGAS),NGPSTR(MXPNGAS),    &
              IBGPSTR(MXPNGAS),NELFGP(MXPSTT),IGSFGP(MXPSTT),     &
              NSTFGP(MXPSTT),MNELFGP(MXPNGAS),MXELFGP(MXPNGAS),   &
              NELFTP(MXPSTT),NSPGPFTP(MXPSTT),IBSPGPFTP(MXPSTT),  &
              ISPGPFTP(MXPNGAS,MXPSTT),NELFSPGP(MXPNGAS,MXPSTT),  &
              NSTFSMSPGP(MXPNSMST,MXPSTT),                        &
              NGRP,NSTTP,MXNSTR,NTSPGP,MAX_STR_SPGP,              &
              MAX_STR_OC_BLK,MXSMCLS,MXSMCLSE,MXSMCLSE1,          &
              NHLFSPGP(MXPSTT),MNHL, NSTFSMGP(MXPNSMST,MXPSTT),   &
              MINMAX_SM_GP(2,MXPSTT),ISTFSMGP(MXPNSMST,MXPSTT),   &
              I_AM_OUT(MXPSTT),N_ELIMINATED_BATCHES,              &
              IELIMINATED_IN_GAS(MXPNGAS),N_ELIMINATED_GAS,       &
              I_ELIMINATE_GAS, I2ELIMINATED_IN_GAS(MXPNGAS),      &
              N_2ELIMINATED_GAS

! Stuff from cstate.fh
INTEGER IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS,NROOT,IROOT(mxroot),IDC,INTSEL
REAL*8        PSSIGN,PLSIGN

! Stuff from crun.fh
CHARACTER(LEN=6) PROPER(20),ENVIRO,CCFORM
CHARACTER(LEN=8) RESP_OP(2,20), CSEQCI(MXPICI,MXPICI), AVE_OP(20),  &
            ITRACI_CR,ITRACI_CN
REAL*8                  THRES_E,XLAMBDA,E_THRE,C_THRE,              &
            E_CONV,C_CONV,RESP_W(20)
INTEGER     MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV,         &
            ICISTR,NOCSF,IDIAG,NOINT,IDMPIN,MXINKA,ICJKAIB,         &
            INIREF,IRESTRF,MOCAA,MOCAB,IPERT,NPERT,                 &
            IAPRREF,IAPRZER,IIDUM,                                  &
            NSEQCI(MXPICI),ISEQCI(MXPICI,MXPICI),IEXTKOP,           &
            IE0AVEX,IC1DSC,IH0SPC,NPTSPC,                           &
            IOCPTSPC(2,MXPNGAS,MXPPTSPC),IH0INSPC(MXPPTSPC),        &
            IRFROOT,NH0EXSPC,IH0EXSPC(MXPPTSPC),INIDEG,             &
            LCSBLK,NOMOFL,IFINMO,NPSSPC,                            &
            NPSSH(MXPIRR,MXPNGAS),ICLSSEL,IDENSI,IPTEKT,            &
            NPTEKT,IH0ROOT,IRST2,ISKIPEI,                           &
            IXYZSYM(3),NPROP,                                       &
            ITRAPRP,NEXCSTATE,IEXCSYM,IRESPONS,NRESP,               &
            MAXORD_OP(2,20),MXITLE,N_AVE_OP, IROOTHOMING,           &
            IUSE_PH, IADVICE,ITRACI,                                &
            IUSE_PA,IPTFOCK,NSXE,                                   &
            ITRA_FI,                                                &
            ITRA_IN,MULSPC,LPAT,IFMULSPC, IRELAX,                   &
            NCNV_RT,                                                &
            I_RE_MS2_SPACE,I_RE_MS2_VALUE,                          &
            ISIMSYM,                                                &
            INOCALC,ISAVE_EXP,IEXPAND

! Stuff from cprnt.fh
Integer      IPRSTR,IPRCIX,IPRORB,IPRDIA,IPRXT,IPROCC,IPRDEN,     &
             IPRRSP,IPRNCIV,IPRPRO,IPRCC

! Stuff from clunit.fh
INTEGER       LU1INT,LU2INT,LUPRP,                                &
              LUDIA,LUC,LUHC,                                     &
              LUSC1,LUSC2,LUSC3,                                  &
              LUSC34,LUSC35,LUSC36,LUSC37,LUSC38,LUSC39,          &
              LUSC40,                                             &
              LUCIVO,LUMOIN,LUMOUT,LUEXC,                         &
              LUSC51,LUSC52,LUSC53

! Stuff from cintfo.fh
INTEGER       I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2,NINT2_NO_CCSYM

! Stuff from cicisp.fh
Integer       IDUMMY,NICISP,                                      &
              NELCI(MXPICI),                                      &
              ISMOST(MXPCSM,MXPCSM),MXSB,MXSOOB,                  &
              NBLKIC(MXPCSM,MXPICI),LCOLIC(MXPCSM,MXPICI),        &
              MXNTTS,MXSOOB_AS
REAL*8        XISPSM(MXPCSM,MXPICI)

! Stuff from cgas.fh
INTEGER     IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),                    &
            NGSOB(MXPOBS,MXPNGAS),                                &
            NGSOBT(MXPNGAS),IGSOCC(MXPNGAS,2),IGSINA,IGSDEL,      &
            IGSOCCX(MXPNGAS,2,MXPICI),NCISPC,                     &
            NCMBSPC, LCMBSPC(MXPICI),ICMBSPC(MXPSTT,MXPICI),      &
            NMXOCCLS,IPHGAS(MXPNGAS),                             &
            IPHGAS1(MXPNGAS)

! Stuff from cecore.fh
REAL*8        ECORE,ECORE_ORIG,ECORE_HEX

! Stuff from spinfo_lucia.fh
!./SPINFO/ new of Nov. 2001
! NCONF_TOT added april 29
Integer MINOP,MAXOP,NCONF_PER_SYM(MXPCSM),                   &
        NCONF_PER_OPEN(MXPORB+1,MXPCSM),                     &
        NPDTCNF(MXPORB+1),NPCSCNF(MXPORB+1),                 &
        NPCMCNF(MXPORB+1),                                   &
        IB_CONF_REO(MXPORB+1),IB_CONF_OCC(MXPORB+1),         &
        IB_SD_FOR_OPEN(MXPORB+1),                            &
        NCSF_PER_SYM(MXPCSM), NSD_PER_SYM(MXPCSM),           &
        NCONF_ALL_SYM, NCONF_ALL_SYM_FOR_OCCLS(MXPCSM),      &
        IBCONF_ALL_SYM_FOR_OCCLS(MXPCSM),nconf_tot,nCSF_HEXS

End Module Lucia_data
