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
!
c***********************************************************************
c Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
c***********************************************************************
      SUBROUTINE POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,
     1  NCN,CNN,IPOTL,PPAR,QPAR,NSR,NLR,DSCM,RREF,PARM,MMLR,VLIM1,
     2  CMM,NCMM,IVSR,IDSTT,RHOAB)
      IMPLICIT NONE
      INTEGER NBOB
      PARAMETER (NBOB=20)
      INTEGER  I,J,IBOB,IAN1,IAN2,IMN1,IMN2,MN1R,MN2R,IORD,IORDD,IPOTL,
     1  PAD,QAD,PNA,NU1,NU2,NT1,NT2,NCMAX,PPAR,QPAR,NCN,NSR,NLR,NVARB,
     2  NPP,LNPT,GNS,GEL, NCMM,IVSR,LVSR,IDSTT,KDER,MM1, MMLR(3)
      CHARACTER*2 NAME1,NAME2
      REAL*8  A0,A1,A2,A3,ALFA,Asw,Rsw,BETA,BINF,B1,B2,BT,CSAV,U1INF,
     1 U2INF,T1INF,T2INF,YPAD,YQAD,YQADSM,YPNA,YPNASM,ABUND,CNN,
     2 DSCM,DX,DX1,FCT,FC1,FC2,FG1,FG2,MASS1,MASS2,RMASS1,RMASS2,REQ,
     3 RREF,Rinn,Rout,SC1,SC2,SG1,SG2,VLIM,VMIN,XDF,X1,XS,XL,XP1,ZZ,ZP,
     4 ZQ,ZME,Aad1,Aad2,Ana1,Ana2,Rad1,Rad2,Rna1,Rna2,fad1e,fad2e,FSW,
     5 ULR,ULRe,RHOAB,REQP,DM(3),DMP(3),DMPP(3),CMM(3),T0,C6adj,C9adj,
     6 RM3,RET,RETsig,RETpi,RETp,RETm,BFCT,PPOW,PVSR,
     7 U1(0:NBOB),U2(0:NBOB),T1(0:NBOB),T2(0:NBOB),PARM(50),
     8 XO(NPP),VV(NPP),RM2(NPP), bTT(-1:2),cDS(-2:0),bDS(-2:0),VLIM1
      SAVE IBOB,IORD,IORDD,PNA,NVARB
      SAVE REQ,U1,U2,T1,T2,CSAV,BINF,ALFA,Rsw,ZME,
     2 Aad1,Aad2,Ana1,Ana2,Rad1,Rad2,Rna1,Rna2,Rinn,Rout,ULR,ULRe
c** Damping function parameters for printout .....
      DATA bTT/2.44d0,2.78d0,3.126d0,3.471d0/
      DATA bDS/3.3d0,3.69d0,3.95d0/
      DATA cDS/0.423d0,0.40d0,0.39d0/
      SAVE bTT, bDS, cDS
c** Electron mass, as per 2006 physical constants
      DATA ZME/5.4857990943d-4/
      LNPT = 1
      PPAR = 5
      QPAR = 3
      NLR = 3
      IORD = NLR
      VLIM = VLIM1
      REQ = 4.1700105834769996
      DSCM = 3.337678701485D+02     
      RREF =8
      PARM(1) =-5.156803528943D-01 
      PARM(2) =-9.585070416286D-02 
      PARM(3) = 1.170797201140D-01 
      PARM(4) =-2.282814434665D-02
      MMLR(1) = 6
      MMLR(2) = 8
      MMLR(3) = 10
      CMM(1)  = 6719000
      CMM(2)  = 112635000.00000000     
      CMM(3)  = 2786940000.0000000
      RHOAB   = 0.54
      IVSR    = -2
      IDSTT   = 1
      NCMM    = 3
      WRITE(6,*) ''
c=======================================================================
c** Generate an MLR potential as per Dattani & Le Roy J.Mol.Spec. 2011
c=======================================================================
      WRITE(6,*) 'IPOTL=',IPOTL,'Beginning to process MLR potential!'
      WRITE(6,*) 'MMLR=',MMLR
      WRITE(6,*) 'CMM=',CMM
      WRITE(6,*) ''
      IF(IPOTL.EQ.4) THEN
          IF(LNPT.GT.0) THEN
            NCN= MMLR(1)
            CNN= CMM(1)
            ULRe= 0.d0
            IF((NCMM.GE.3).AND.(MMLR(2).LE.0)) THEN
                RET= 9.36423830d-4*REQ
            ELSE
            KDER= 0
            CALL dampF(REQ,RHOAB,NCMM,MMLR,IVSR,IDSTT,KDER,DM,DMP,DMPP)
                  DO  J= 1,NCMM
                          ULRe= ULRe + DM(J)*CMM(J)/REQ**MMLR(J)
                  ENDDO
            WRITE(6,*) 'Made it out of dampF'      
            ENDIF
            BINF= DLOG(2.d0*DSCM/ULRe)
            WRITE(6,602) NCN,PPAR,QPAR,DSCM,REQ
            WRITE(6,607) PPAR,PPAR,QPAR,NSR,NLR,IORD+1,
     1                                       (PARM(J),J= 1,IORD+1)
            WRITE(6,613) RREF
              IF(RHOAB.GT.0.d0) THEN
                  PVSR= 0.5d0*IVSR
                  IF(IDSTT.GT.0) THEN
                      PVSR= 0.5d0*IVSR
                      WRITE(6,664) RHOAB,PVSR,bDS(IVSR),cDS(IVSR),PVSR
                  ELSE
                      LVSR= IVSR/2
                      WRITE(6,666) RHOAB,LVSR,bTT(LVSR)
                  ENDIF
              ENDIF
            WRITE(6,617) BINF,MMLR(1),CMM(1),MMLR(1) 
              IF(NCMM.GT.1) THEN
                  DO  I= 2,NCMM !Removed IF stmnt that prints C10 nicely
                    WRITE(6,619) MMLR(I),CMM(I),MMLR(I)
                  ENDDO
              ENDIF
          ENDIF
c  Loop over distance array XO(I)
          WRITE(6,*) XO
          WRITE(6,*) 'PPAR=',PPAR
          WRITE(6,*) 'REQ=',REQ
          WRITE(6,*) 'XO=',XO(1)
          DO  I= 1,NPP
              ZZ= (XO(i)**PPAR- REQ**PPAR)/(XO(i)**PPAR+ REQ**PPAR)
              ZP= (XO(i)**PPAR-RREF**PPAR)/(XO(i)**PPAR+RREF**PPAR)
              ZQ= (XO(i)**QPAR-RREF**QPAR)/(XO(i)**QPAR+RREF**QPAR)
              BETA= 0.d0
              IF(ZZ.GT.0) IORDD= NLR
              IF(ZZ.LE.0) IORDD= NSR
              DO  J= IORDD,0,-1
                  BETA= BETA*ZQ+ PARM(J+1)
              ENDDO
c  Calculate MLR exponent coefficient
              BETA= BINF*ZP + (1.d0- ZP)*BETA
              ULR= 0.d0
c** Calculate local value of uLR(r)
              IF((NCMM.GE.3).AND.(MMLR(2).LE.0)) THEN
                  RM3= 1.d0/XO(I)**3
                ELSE
c                 IVSR gets corrupted so make sure it's -2.                
c                 IVSR=-2
                  WRITE(6,*) 'IVSR=',IVSR
                  IF(RHOAB.GT.0.d0) CALL dampF(XO(I),RHOAB,NCMM,MMLR,
     1                                    IVSR,IDSTT,KDER,DM,DMP,DMPP)
                  DO  J= 1,NCMM
                      IF(RHOAB.LE.0.d0) THEN
                          ULR= ULR + CMM(J)/XO(I)**MMLR(J)
                      ELSE
                          ULR= ULR + DM(J)*CMM(J)/XO(I)**MMLR(J)
                      ENDIF
                  ENDDO
              ENDIF
              BETA= (ULR/ULRe)*DEXP(-BETA*ZZ)
              VV(I)= DSCM*(1.d0 - BETA)**2 - DSCM + VLIM
              WRITE(6,*) 'VLIM=',VLIM
              WRITE(6,*) 'DSCM=',DSCM
              WRITE(6,*) 'ZZ=',ZZ
              WRITE(6,*) 'ULRe=',ULRe
              WRITE(6,*) 'ULR=',ULR
              WRITE(6,*) 'BETA=',BETA
              WRITE(6,*) 'VV=',VV
          ENDDO
      ENDIF
      WRITE(6,*) 'Finishing MLR generation...'
c
      RETURN
  602 FORMAT(/' MLR(n=',i1,'; p=',I1,', q=',I1,') Potential with:   De='
     1 ,F10.3,'[cm-1]    Re=',F12.8,'[A]')
  607 FORMAT('   with exponent coefficient   beta(r)= beta{INF}*y',I1,
     1  ' + [1-y',i1,']*Sum{beta_i*y',i1,'^i}'/6x,'exponent coefft. powe
     2r series orders',I4,' for  R < Re  and',I4,' for  R > Re'/6x,
     3  'and',i3,' coefficients:',1PD16.8,2D16.8:/(10x,4D16.8:))
  613 FORMAT(6x,'with radial variables  y_p & y_q  defined w.r.t.',
     1 '  RREF=',F10.7)
  615 FORMAT(6x,'radial variables  y_p & y_q  defined w.r.t.',
     1 '  RREF= Re=' F10.7)
  614 FORMAT(/' Potential is an SPF expansion in  (r-Re)/(',F5.2,
     1  '* r)  with   Re=',f12.9/5x,'De=',g18.10,'   b0=',
     2  1PD16.9,'   and',i3,'  b_i  coefficients:'/(5D16.8))
  616 FORMAT(/' Potential is an O-T expansion in  (r-Re)/[',f5.2,
     1  '*(r+Re)]  with   Re=',f12.9/5x,'De=',G18.10,
     2  '   c0=',1PD16.9,'   and',i3,'  c_i coefficients:'/(5D16.8))
  617 FORMAT('   while betaINF=',f12.8,'  & uLR defined by  C',i1,' =',
     1  1PD13.6,'[cm-1 Ang','^',0P,I1,']')
  618 FORMAT(/' Potential is a general GPEF expansion in  (r**',i1,
     1  ' - Re**',i1,')/(',SP,F5.2,'*r**',SS,i1,SP,F6.2,'*Re**',SS,i1,
     2  ')'/5x,'with   Re=',f12.9,'   De=',g18.10,'   g0=',1PD16.9/
     3  5x,'and',i3,'  g_i coefficients:  ',3D16.8/(5D16.8:))
  619 FORMAT(50x,'C',I1,' =',1PD13.6,'[cm-1 Ang','^',0P,I1,']')
  621 FORMAT(50x,'C',I2,'=',1PD13.6,'[cm-1 Ang','^',0P,I2,']')
  620 FORMAT(/' Potential is a power series in  r  of  order',i3,
     1 ' with   V(r=0)=',f11.4/3x,'& coefficients (from linear term):',
     2 1P2d16.8:/(5x,4D16.8:))
  626 FORMAT('   De=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]   and'/
     1  '     Damping function  D(r)= exp[ -',0Pf6.4,'*(',f7.4,
     1  '/X -1.0)**',f5.2,']')
  636 FORMAT(3x,'where   fsw(r) = 1/[1 - exp{',f7.4,'*(r -',f7.4,')}]')
  642 FORMAT(' where for  r < Rinn=',F7.4,'   V=',SP,F12.4,1x,1PD13.6,
     1  '/R**12' )
  644 FORMAT('  and  for  r > Rout=',F7.3,'   V= VLIM ',
     1 (SP,1PD14.6,'/r**',SS,I2):/(39x,SP,1PD14.6,'/r**',SS,I2))
  652 FORMAT(6x,'where the radial variable   y_',I1,'= (r**',I1,' - RREF
     4**',i1,')/(r**',I1,' + RREF**',i1, ')')
  654 FORMAT(10x,'is defined w.r.t.   RREF=',F11.8)
  656 FORMAT(10x,'is defined w.r.t.   RREF= Re= ',F11.8)
  664 FORMAT(4x,'uLR inverse-power terms incorporate DS-type damping wit
     1h   RHOAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',SS,f4.1,'}'/8x,'Dm(r)= [1 - exp(-',f5.2,
     3 '(RHOAB*r)/m -',f6.3,'(RHOAB*r)^2/sqrt{m})]^{m',SP,F4.1,'}')
  666 FORMAT(4x,'uLR inverse-power terms incorporate TT-type damping wit
     1h   RHOAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',I2,'}'/8x,'Dm(r)= [1 - exp(-bTT*r)*SUM{(bTT*r
     3)^k/k!}]   where   bTT=',f6.3,'*RHOAB')
      END

c***********************************************************************
      SUBROUTINE dampF(r,RHOAB,NCMM,MMLR,IVSR,IDSTT,KDER,DM,DMP,DMPP)
      INTEGER NCMM,NCMMax,MMLR(NCMM),IVSR,IDSTT,KDER,IDFF,FIRST,
     1  Lsr,m,MM,MMAX
      REAL*8 r,RHOAB,bTT(-2:2),cDS(-4:0),bDS(-4:0),aTT,br,XP,YP,
     1  TK, DM(NCMM),DMP(NCMM),DMPP(NCMM),SM(-3:25),
     2  bpm(20,-2:0), cpm(20,-2:0),ZK
       DATA bTT/2.10d0,2.44d0,2.78d0,3.126d0,3.471d0/
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.40d0,0.39d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
      WRITE(6,*) 'Made it inside of dampF! IVSR=',IVSR
      WRITE(6,*) 'IDSTT=',IDSTT
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  IDFF= -2,0
                      bpm(m,IDFF)= bDS(IDFF)/DFLOAT(m)
                      cpm(m,IDFF)= cDS(IDFF)/DSQRT(DFLOAT(m))
                  ENDDO
              ENDDO
              FIRST= 0
          ENDIF
          br= RHOAB*r
          WRITE(6,*) 'NCMM=',NCMM
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,IVSR) + cpm(MM,IVSR)*br)*br)
              YP= 1.d0 - XP
              ZK= MM-1.d0
              DM(m)= YP**(MM-1)
c... Actually ...  DM(m)= YP**(MM + IVSR/2)  :  set it up this way to
c   avoid taking exponential of a logarithm for fractional powers (slow)
              IF(IVSR.EQ.-4) THEN
                  ZK= ZK- 1.d0
                  DM(m)= DM(m)/YP
                  ENDIF
              IF(IVSR.EQ.-3) THEN
                  ZK= ZK- 0.5d0
                  DM(m)= DM(m)/DSQRT(YP)
                  ENDIF
              IF(IVSR.EQ.-1) THEN
                  ZK= ZK+ 0.5d0
                  DM(m)= DM(m)*DSQRT(YP)
                  ENDIF
              IF(IVSR.EQ.0) THEN
                  ZK= MM
                  DM(m)= DM(m)*YP
                  ENDIF
          ENDDO
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  IDSTT=',i3,'   IVSR=',i3,'  no dampi
     1ng function is defined')
  602 FORMAT( /,' ***ERROR ***  RHOAB=', F7.4,'  yields an invalid Dampi
     1ng Function definition')
      END
