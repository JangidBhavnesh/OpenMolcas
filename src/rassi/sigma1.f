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
      SUBROUTINE SIGMA_1(ISGSTRUCT,SGS,ICISTRUCT,IXSTRUCT,
     &                 NMIDV,MXEO,NVTAB,NICOUP,ISM,
     &                 IP,IQ,CPQ,ISYCI,CI,SGM,NOCSF,IOCSF,NOW,
     &                 IOW,NOCP,IOCP,ICOUP,VTAB,MVL,MVR)
      use Struct, only: nSGSize, nCISize, nXSize, SGStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ISM(*)
      DIMENSION NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION NOCP(MXEO,NSYM,NMIDV),IOCP(MXEO,NSYM,NMIDV)
      DIMENSION VTAB(NVTAB)
      DIMENSION ICOUP(3,NICOUP)
      DIMENSION MVL(NMIDV,2),MVR(NMIDV,2)
      DIMENSION CI(*),SGM(*)
#include "symmul.fh"
      Dimension iSGStruct(nSGSize)
      Type (SGStruct) SGS
      Dimension iCIStruct(nCISize)
      Dimension iXStruct (nXSize)
#include "WrkSpc.fh"
      SAVE ICALL
      DATA ICALL /0/




      NLEV  =ISGSTRUCT(2)
      MIDLEV=ISGSTRUCT(8)
      NLEV  =SGS%nLev
      MIDLEV=SGS%MidLev
      NIPWLK=ICISTRUCT(2)
      LICASE=ICISTRUCT(9)
      NT1MX =IXSTRUCT(10)
      NT2MX =IXSTRUCT(11)
      NT3MX =IXSTRUCT(12)
      NT4MX =IXSTRUCT(13)

      IF(ABS(CPQ).LT.1.0D-12) RETURN
      ICALL=ICALL+1
C SYMMETRY OF ORBITALS:
      ISYP=ISM(IP)
      ISYQ=ISM(IQ)
      ISYPQ=MUL(ISYP,ISYQ)
C SYMMETRY OF SIGMA ARRAY:
      ISYSGM=MUL(ISYPQ,ISYCI)
      IF(IP.GT.IQ) THEN
C THEN THIS IS AN EXCITING OPERATOR.
        IF(IP.LE.MIDLEV) GOTO 1600
        IF(IQ.GT.MIDLEV) GOTO 1700
        GOTO 1800
      ELSE IF(IQ.GT.IP) THEN
C THEN THIS IS A DEEXCITING OPERATOR.
        IF(IQ.LE.MIDLEV) GOTO 1300
        IF(IP.GT.MIDLEV) GOTO 1400
        GOTO 1500
      END IF
      IF(IP.GT.MIDLEV) GOTO 1200
C SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
C IP=IQ < MIDLEV.
      DO 100 MVSGM=1,NMIDV
        DO 101 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 101
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOLW=IOW(2,ISYDSG,MVSGM)
          IPSHFT=2*(IP-1)
          LLW=LICASE+IOLW-NIPWLK+IPSHFT/30
          IPSHFT=MOD(IPSHFT,30)
          IPPOW=2**IPSHFT
          DO 102 J=1,NDWNSG
            JC=IWORK(LLW+J*NIPWLK)
CPAM96            ICS=IAND(ISHFT(JC,-IPSHFT),3)
            ICS=MOD(JC/IPPOW,4)
            IF(ICS.EQ.0) GOTO 102
            X=CPQ*DBLE((1+ICS)/2)
            JSTA=ISGSTA+NUPSG*(J-1)
            CALL DAXPY_(NUPSG,X,CI(JSTA),1,SGM(JSTA),1)
102       CONTINUE
101     CONTINUE
100   CONTINUE
      RETURN
1200  CONTINUE
C SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
C IP=IQ>MIDLEV
      DO MVSGM=1,NMIDV
        DO ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) CYCLE
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOUW=IOW(1,ISYUSG,MVSGM)
          IPSHFT=2*(IP-1-MIDLEV)
          LUW=LICASE+IOUW-NIPWLK+IPSHFT/30
          IPSHFT=MOD(IPSHFT,30)
          IPPOW=2**IPSHFT
          DO I=1,NUPSG
            IC=IWORK(LUW+I*NIPWLK)
CPAM96            ICS=IAND(ISHFT(IC,-IPSHFT),3)
            ICS=MOD(IC/IPPOW,4)
            IF(ICS.EQ.0) cycle
            X=CPQ*DBLE((1+ICS)/2)
            ISTA=ISGSTA-1+I
            CALL DAXPY_(NDWNSG,X,CI(ISTA),NUPSG,SGM(ISTA),NUPSG)
          END DO
        END DO
      END DO
      RETURN
1300  CONTINUE
C DEEXCITING OPERATOR, IP<IQ<=MIDLEV.
      DO 300 MVSGM=1,NMIDV
        DO 301 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 301
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYDC=MUL(ISYPQ,ISYDSG)
          NDWNC=NOW(2,ISYDC,MVSGM)
          IF(NDWNC.EQ.0) GOTO 301
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUSG,MVSGM,ISYCI)
          INDEO=2*NLEV+(IQ*(IQ-1))/2+IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
C CASE IS: LOWER HALF, DEEXCITE:
            CALL DEX1 (CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA),
     *             NCP,ICOUP(1,LICP),VTAB)
          END IF
301     CONTINUE
300   CONTINUE
      RETURN
1400  CONTINUE
C DEEXCITING OPERATOR, MIDLEV<IP<IQ
      DO 400 MVSGM=1,NMIDV
        DO 401 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 401
          ISYUC=MUL(ISYPQ,ISYUSG)
          NUPC=NOW(1,ISYUC,MVSGM)
          IF (NUPC.EQ.0) GOTO 401
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOC=IOCSF(ISYUC,MVSGM,ISYCI)
          INDEO=2*NLEV+(IQ*(IQ-1))/2+IP
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
C CASE IS: UPPER HALF, DEEXCITE:
            CALL DEX2 (CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA),
     *             NCP,ICOUP(1,LICP),VTAB)
          END IF
401     CONTINUE
400   CONTINUE
      RETURN
1500  CONTINUE
C DEEXCITING CASE, IP<=MIDLEV<IQ.
C ALLOCATE TEMPORARY WORK AREA:
        NTMPMX=MAX(NT1MX,NT2MX)
        CALL GETMEM('TMP   ','ALLO','REAL',LTMP,NTMPMX)
      DO 500 MVSGM=1,NMIDV
        MV4=MVR(MVSGM,1)
        MV5=MVR(MVSGM,2)
        DO 501 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 501
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYUC=MUL(ISYQ,ISYUSG)
          ISYDC=MUL(ISYP,ISYDSG)
          IF(MV4.EQ.0) GOTO 499
          NUPC=NOW(1,ISYUC,MV4)
          IF(NUPC.EQ.0) GOTO 499
          NDWNC=NOW(2,ISYDC,MV4)
          IF(NDWNC.EQ.0) GOTO 499
          INDEO=IQ
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.EQ.0) GOTO 499
          NTMP=NUPSG*NDWNC
          CALL FZERO(WORK(LTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUC,MV4,ISYCI)
C CASE IS: UPPER HALF, DEEXCITE:
          CALL DEX2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LTMP),
     *             NCP,ICOUP(1,LICP),VTAB)
          INDEO=IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
C CASE IS: LOWER HALF, DEEXCITE:
            X=1.0D00
            CALL DEX1 (X,NUPSG,WORK(LTMP),SGM(ISGSTA),
     *               NCP,ICOUP(1,LICP),VTAB)
          END IF
499       CONTINUE
          IF(MV5.EQ.0) GOTO 501
          NUPC=NOW(1,ISYUC,MV5)
          IF(NUPC.EQ.0) GOTO 501
          NDWNC=NOW(2,ISYDC,MV5)
          IF(NDWNC.EQ.0) GOTO 501
          INDEO=NLEV+IQ
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.EQ.0) GOTO 501
          NTMP=NUPSG*NDWNC
          CALL FZERO(WORK(LTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUC,MV5,ISYCI)
C CASE IS: UPPER HALF, DEEXCITE:
          CALL DEX2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LTMP),
     *             NCP,ICOUP(1,LICP),VTAB)
          INDEO=NLEV+IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
C CASE IS: LOWER HALF, DEEXCITE:
            X=1.0D00
            CALL DEX1 (X,NUPSG,WORK(LTMP),SGM(ISGSTA),
     *               NCP,ICOUP(1,LICP),VTAB)
          END IF
501     CONTINUE
500   CONTINUE
      CALL GETMEM('TMP   ','FREE','REAL',LTMP,NTMPMX)
      RETURN
1600  CONTINUE
C EXCITING CASE, IQ<IP<=MIDLEV.
      DO 600 MVSGM=1,NMIDV
        DO 601 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 601
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYDC=MUL(ISYPQ,ISYDSG)
          NDWNC=NOW(2,ISYDC,MVSGM)
          IF(NDWNC.EQ.0) GOTO 601
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUSG,MVSGM,ISYCI)
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          NCP=NOCP(INDEO,ISYDC,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDC,MVSGM)
C CASE IS: LOWER HALF, EXCITE:
            CALL EXC1 (CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA),
     *             NCP,ICOUP(1,LICP),VTAB)
          END IF
601     CONTINUE
600   CONTINUE
      RETURN
1700  CONTINUE
C EXCITING CASE, MIDLEV<IQ<IP
      DO 700 MVSGM=1,NMIDV
        DO 701 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 701
          ISYUC=MUL(ISYPQ,ISYUSG)
          NUPC=NOW(1,ISYUC,MVSGM)
          IF (NUPC.EQ.0) GOTO 701
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOC=IOCSF(ISYUC,MVSGM,ISYCI)
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          NCP=NOCP(INDEO,ISYUC,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYUC,MVSGM)
C CASE IS: UPPER HALF, EXCITE:
            CALL EXC2 (CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA),
     *             NCP,ICOUP(1,LICP),VTAB)
          END IF
701     CONTINUE
700   CONTINUE
      RETURN
1800  CONTINUE
C EXCITING CASE, IQ<=MIDLEV<IP
C ALLOCATE TEMPORARY WORK AREAS:
      NTMPMX=MAX(NT3MX,NT4MX)
      CALL GETMEM('TMP   ','ALLO','REAL',LTMP,NTMPMX)
      DO 800 MVSGM=1,NMIDV
        MV1=MVL(MVSGM,2)
        MV2=MVL(MVSGM,1)
        IF((MV1.EQ.0).AND.(MV2.EQ.0)) GOTO 800
        DO 801 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 801
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYUC=MUL(ISYP,ISYUSG)
          ISYDC=MUL(ISYQ,ISYDSG)
          IF(MV2.EQ.0) GOTO 799
          NUPC=NOW(1,ISYUC,MV2)
          IF(NUPC.EQ.0) GOTO 799
          NDWNC=NOW(2,ISYDC,MV2)
          IF(NDWNC.EQ.0) GOTO 799
          INDEO=IP
          NCP=NOCP(INDEO,ISYUC,MV2)
          IF(NCP.EQ.0) GOTO 799
          NTMP=NUPSG*NDWNC
          CALL FZERO(WORK(LTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUC,MV2)
          IOC=IOCSF(ISYUC,MV2,ISYCI)
C CASE IS: UPPER HALF, EXCITE:
          CALL EXC2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LTMP),
     *             NCP,ICOUP(1,LICP),VTAB)
          INDEO=IQ
          NCP=NOCP(INDEO,ISYDC,MV2)
          IF(NCP.EQ.0) GOTO 799
          LICP=1+IOCP(INDEO,ISYDC,MV2)
C CASE IS: LOWER HALF, EXCITE:
          X=1.0D00
          CALL EXC1 (X,NUPSG,WORK(LTMP),SGM(ISGSTA),
     *               NCP,ICOUP(1,LICP),VTAB)
799       CONTINUE
          IF(MV1.EQ.0) GOTO 801
          NUPC=NOW(1,ISYUC,MV1)
          IF(NUPC.EQ.0) GOTO 801
          NDWNC=NOW(2,ISYDC,MV1)
          IF(NDWNC.EQ.0) GOTO 801
          INDEO=NLEV+IP
          NCP=NOCP(INDEO,ISYUC,MV1)
          IF(NCP.EQ.0) GOTO 801
          NTMP=NUPSG*NDWNC
          CALL FZERO(WORK(LTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUC,MV1)
          IOC=IOCSF(ISYUC,MV1,ISYCI)
C CASE IS: UPPER HALF, EXCITE:
          CALL EXC2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LTMP),
     *             NCP,ICOUP(1,LICP),VTAB)
          INDEO=NLEV+IQ
          NCP=NOCP(INDEO,ISYDC,MV1)
          IF(NCP.EQ.0) GOTO 801
          LICP=1+IOCP(INDEO,ISYDC,MV1)
C CASE IS: LOWER HALF, EXCITE:
          X=1.0D00
          CALL EXC1 (X,NUPSG,WORK(LTMP),SGM(ISGSTA),
     *               NCP,ICOUP(1,LICP),VTAB)
801   CONTINUE
800   CONTINUE
      CALL GETMEM('TMP   ','FREE','REAL',LTMP,NTMPMX)

      END SUBROUTINE SIGMA_1
