c=======================================================================
      PROGRAM MAIN
c=======================================================================
c-------------------  Last updated 9 February 2004 ---------------------
c=======================================================================
ccc   INCLUDE 'arrsizes.h'
c-----------------------------------------------------------------------
c  Utility routine to summarize dimensioning of arrays
c-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,
     1        mxsets,mxfreq
      REAL*8 CCM,PI
c-----------------------------------------------------------------------
c  mxdata - maximum number of input data points
c  mxisp  - maximum number of points for initial state potential array
c           (also used for number of points in transition moment array)
c  mxnj    - maxiumum value of j quantum number allowed
c  mxfsp  - maximum number of points for final state potential array
c  mxnp   - maximum number of parameters total
c  mxntp  - maximum number of turning points to be read in
c  mxprm  - maximum number of parameters for final state pot'l or TMF
c  mxv    - largest value for the v quantum number
c  mxfs   - maximum number of final states allowed
c  mxisot - maximum number of isotopomers allowed
c  mxsets - maximum number of data sets allowed
c  mxfreq - maximum number of data points allowed in a given set
c-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
c-----------------------------------------------------------------------
      INTEGER AN1,AN2,C,CHARGE,ERR,FITIT,FREQYN,MN1(mxisot),MN2(mxisot),
     1  GEGS1,GEGS2,GNS1,GNS2, I,IFR,IFS,ILRF,IP,IR2F,IROUND,ISET,
     2  ISO,IVJ,IWRSCH,IWROVR,J,JFRPW,JREF,K,LNPT,LPDER,LPFS,LPRINT,
     3  LPTMF,LPRWF,M,MAXV,MCALC, N,NCNF,NCNI,NFS,NIN,NISTP,NPRS,NPRF,
     4  NPRM,NDATA,NBEG,NEND,NPTOT,NPPFREE,NSETS,NUSEF, PRINTYN,RPD,V,
     5  NGPRND, BOLTZ(mxsets),DTYPE(mxsets),IFRPW(mxsets),ISOT(mxsets),
     6  J1ST(mxsets),VMAX(mxsets),JMAX(mxsets),NFREQ(mxsets),NJ(mxsets),
     7  NVJ(mxsets),PUNITS(mxsets),PQR(mxsets),SCALE(mxsets),
     7  GFS(mxfs),
     8  CD(mxsets,mxfs),CN(mxsets,mxfs),STYPE(mxfs),FSVAR(mxprm,mxfs),
     9  FSTYPE(mxfs),NFSPRM(mxfs),NTPFS(mxfs),OMEGA(0:mxfs),OTMF(mxfs),
     a  TMFTYP(mxfs),TMFVAR(0:mxprm-1,mxfs),V1ST(mxsets),VFIX(mxfreq),
     b  JFIX(mxfreq),XCOORD(mxfs), INNER(0:mxv,mxisot),INNR(0:mxv)
c-----------------------------------------------------------------------
      REAL*8 AB1,AB2,ADD,AVALUE(mxfs),BFCT(mxisot),BVALUE(mxfs),
     1   CM(mxnp,mxnp),CNNF,DFREQ,dIdT(0:mxprm-1,mxfreq,mxfs,mxsets),
     2   DFACT,DSE,EPS,FACTOR(mxsets),FCT(mxfsp),FREQ(mxfreq,mxsets),
     3   FREQ1,FSPRM(mxprm,mxfs),MASS1(mxisot),MASS2(mxisot),
     4   GV(0:mxv),RCNST(7),MCI(0:mxv,0:7,mxisot),MU(mxisot),
     5   OBS(mxfreq,mxsets),OVRCRT,OUTPUT(mxfreq,mxsets),POPCRT,
     6   PS(mxnp),PU(mxnp),PV(mxnp),RAD(mxfsp),REXFS,REXTMF,RFACTF,RH,
     7   RMAX,RMAXOUT,RMIN,RMS(mxsets),RTP(2,mxfs),RTPF(mxfsp),
     8   SF(mxsets),
     8   SQRTMU(mxisot),TMFPRM(0:mxprm-1,mxfs),TEMP(mxsets),TSTPS,TSTPU,
     9   UNC(mxfreq,mxsets),VFACTF,VF(mxfsp,mxfs),VI(mxisp,mxisot),
     a   VIS(mxisp),VLIMI,VLIMF(mxfs),VSHFTF(mxfs),VTP(2,mxfs),
     b   RM2(mxisp),VICD(mxisp,mxisot),VTPF(mxfsp),WAVL(mxfreq,mxsets),
     c   XM2(mxfsp),YD(mxdata),YO(mxdata),YU(mxdata),zfs(mxfsp,mxfs),
     d   WF1(mxisp),ztmf(mxisp,mxfs),TWOPIC, PMAX,BvWN,GAMA, UCUTOFF
c-----------------------------------------------------------------------
      CHARACTER*2 N1,N2,NAME1(mxisot),NAME2(mxisot)
      CHARACTER*75 TITLE
      CHARACTER*70 INFO(mxsets)
c-----------------------------------------------------------------------
      COMMON /MF/ BFCT,EPS,FACTOR,FREQ,MCI,OBS,UNC,OVRCRT,POPCRT,TEMP,
     1   VI,VICD,VLIMI,WAVL,BOLTZ,CD,CN,DTYPE,FITIT,INNER,ISOT,IWRSCH,
     2   IWROVR,J1ST,JFIX,JFRPW,JMAX,MCALC,NJ,OMEGA,PQR,printyn,
     3   V1ST,VFIX,VMAX
      COMMON /MD/ DFACT,FSVAR,LPDER,NPPFREE
      COMMON /MGt/ REXTMF
      COMMON /MGf/ AVALUE,BVALUE,REXFS,RTP,VTP,zfs,FSTYPE,XCOORD
      COMMON /MFD/ dIdT,OUTPUT,SF,IFRPW,NFREQ,NFS,NSETS,NVJ,RPD,SCALE,
     1             TMFVAR,INFO
      COMMON /MFGf/ VF,VLIMF
      COMMON /MFGt/ XM2,ztmf,LPTMF,NIN,TMFTYP
      COMMON /MDGf/ FSPRM,NFSPRM
      COMMON /MFGtGf/ RAD
      COMMON /MFDGt/ TMFPRM,OTMF,GFS
c?????
      real*8 xtmf
c?????
c-----------------------------------------------------------------------
c  LNPT needed in prepot.f (must be greater that zero here)
c  NPRS needed in prepot.f subroutine genint
c  MCALC  - specifies how radial overlap calculations are performed
c            <= 0 to use delta function approximation (UNAVAILABLE)
c            >  0 for exact quantal calculation
c  EPS is converged precision (in cm-1) of calculated eigenvalues
c  POPCRT is fraction of initial-state population to be included when
c         the thermal sums terminate in direct sums over the rotational
c         or vibrational populations.
c-----------------------------------------------------------------------
      DATA LNPT/1/,NPRS/1/,NDATA/0/,NGPRND/1/
      EPS= 1.d-05
      POPCRT= 0.9999d0
      MCALC= 1
      TWOPIC= 2.d0*PI*CCM
c=======================================================================
c TITLE is a title or output header of up to 75 characters, read on a
c   single line enclosed between single quotes: e.g.  'title of problem'
c=======================================================================
      READ(5,*) TITLE
c=======================================================================
c  AN1   - atomic number of atom 1
c  AN2   - atomic number of atom 2
c  CHARGE - +/-(integer) charge on the molecule
c  NISTP  - number of isotopomers to be considered
c  NFS    - number of final state potentials considered
c  NSETS  - the number of data sets to be predicted or fitted to
c  FITIT  - if > 1 do a fit to experimental data
c           if = 0 do forward calculation only
c=======================================================================
      READ(5,*) AN1, AN2, CHARGE, NISTP, NFS, NSETS, FITIT
c=======================================================================
      IF(NSETS.gt.mxsets) THEN
         WRITE(6,617) mxsets
         STOP
      ENDIF
      IF(NISTP.gt.mxisot) THEN
         WRITE(6,619) mxisot
         STOP
      ENDIF
      IF(NFS.gt.mxfs) THEN
         WRITE(6,621) mxfs
         STOP
      ENDIF
      WRITE(6,600) TITLE,NISTP
      DO  iso= 1,NISTP
c** Loop over isotopomers: read mass numbers & calcuate reduced masses
c  Read data for lightest isotomer first and repeat for subsequent ones
c=======================================================================
c  MN1(iso)   - mass number of atom 1 of isotopomer 'iso'
c  MN2(iso)   - mass number of atom 2 of isotopomer 'iso'
c=======================================================================
          READ(5,*) MN1(iso), MN2(iso)
c=======================================================================
          IF((AN1.GT.0).AND.(AN1.LE.109)) THEN
              CALL MASSES(AN1,MN1(iso),NAME1(iso),GEGS1,GNS1,
     1                                                 MASS1(iso),AB1)
            ELSE
c ... If particle-1 is not a normal stable atomic species, read in a
c  particle name (2-characters surrounded by 's) and mass in [u]
c=======================================================================
              READ(5,*) NAME1(iso), MASS1(iso)
c=======================================================================
            ENDIF
          IF((AN2.GT.0).AND.(AN2.LE.109)) THEN
              CALL MASSES(AN2,MN2(iso),NAME2(iso),GEGS2,GNS2,
     1                                                 MASS2(iso),AB2)
            ELSE
c ... If particle-2 is not a normal stable atomic species, read in a
c  particle name (2-characters surrounded by 's) and mass in [u]
c=======================================================================
              READ(5,*) NAME2(iso), MASS2(iso)
c=======================================================================
            ENDIF
c-----------------------------------------------------------------------
          MU(iso)= MASS1(iso)*MASS2(iso)/(MASS1(iso)+MASS2(iso))
          SQRTMU(iso)= DSQRT(MU(iso))
          WRITE(6,602) NAME1(iso),MN1(iso),NAME2(iso),MN2(iso),
     1                                   MASS1(iso),MASS2(iso),MU(iso)
          ENDDO
c=======================================================================
c  RH     - radial mesh size (in angstroms) for numerical integration
c  RMIN   - lower limit (in angstroms) for numerical integration
c  RMAX   - upper limit (in angstroms) for numerical integration is
c           internally set to the smaller of the read-in RMAX, or the
c           largest distance allowed by array dimensions
c  OVRCRT - convergence criterion for exact quantal calculation of
c           continuum wavefunctions
c=======================================================================
      READ(5,*) RH, RMIN, RMAX, OVRCRT
c=======================================================================
      IF(FITIT.GT.0) THEN
          WRITE(6,606) NSETS
        ELSE
          FITIT= 0
          WRITE(6,610) NSETS
        ENDIF
      IF(MCALC.GT.0) THEN
          WRITE(6,626) OVRCRT
        ELSE
          WRITE(6,628)
          STOP
        ENDIF
      NIN= (RMAX - RMIN)/RH + 1
c-----------------------------------------------------------------------
c   NIN  is the number of initial-state potential mesh points
c-----------------------------------------------------------------------
      IF(NIN.GE.mxisp) NIN = mxisp
      RMAXOUT= RMIN + NIN*RH
      WRITE(6,652) RMIN,RMAXOUT,RH
c======================================================================
c  IWRSCH controls print level in matrix element subroutines in SCHRQ
c     < 0  allows only printing of errors/warnings [normal setting]
c     = 0  prevents all printing there
c     > 0  prints results; > 1 more details; ; > 2 even more details; ...
c IWROVR  controls print level inside OVRLAP and OVRPD
c     < 0  allows only printing of errors/warnings [normal setting]
c     = 0  no print
c     = 1  prints overlap integrals and converged amplitudes
c     > 1  prints additional details
c=======================================================================
      READ(5,*) IWRSCH, IWROVR
c=======================================================================
      IF(FITIT.GT.0) THEN
c=======================================================================
c IROUND controls rounding in nllssrr subroutine
c    not equal 0 causes sequential rounding and re-fitting to be
c      performed, with each parameter being rounded at the |iround|th
c      significant digit of its local uncertainty
c    > 0 rounding selects in turn remaining parameter w/ largest
c      relative uncertainty
c    < 0 rounds parameters sequentially from last to first
c    = 0 simply stops after full convergence (no rounding)
c LPDER : if greater than zero, partial derivative array is written
c            to channel 10
c UCUTOFF : ignore input data with uncertaities > UCUTOFF
c DFACT : (real) scaling factor in calculation of derivatives by
c        differences in dyidpj subroutine; typically set to 1.0 (this
c        is still under development)
c LPRINT : specifies the level of printing inside nllssrr.f
c     = 0  no print except for failed convergence.
c     < 0  only converged, unrounded parameters, PU & PS's
c     >= 1  print converged parameters, PU & PS's
c     >= 2  also print parameter change each rounding step
c     >= 3  also indicate nature of convergence
c     >= 4  also print convergence tests on each cycle
c     >= 5  also parameters changes & uncertainties, each cycle
c=======================================================================
          READ(5,*) IROUND, LPDER, UCUTOFF, DFACT, LPRINT
c=======================================================================
          IF(IROUND.NE.0) WRITE(6,607) IABS(IROUND)
          IF(IROUND.GT.0) WRITE(6,608)
          IF(IROUND.LT.0) WRITE(6,609)
          ENDIF
      DO 10 iset= 1,NSETS
c=======================================================================
c INFO  is a name for the data set, of up to 70-characters, read in on a
c    single line, enclosed in single quotes: e.g.  'name for data set'
c ISOT  specifies the isotopomer to be considered for this data set
c BOLTZ > 0 for calculation with Boltzman population of initial states
c        = 0 does calculation for a specific (v,J) level
c DTYPE = 1  if property is (a sum) of final-state intensities
c       = 2  if property is ratio of (sums of) final-state intensities
c IFRPW  identifies type of data in this set
c            =  0 for predissociation
c            =  1 for (decadic) molar absorption coefficients
c            =  3 for spontaneous emission
c            = -1 for constant freqency factor in absorption
c            = -3 for constant freqency factor in emission
c PQR   = 0  if using the Q-branch approximation (J'= J")
c        = 1  if the P, Q, and R braches are calculated, weighted by the
c             Honl-London factor and summed over.
c=======================================================================
          READ(5,*) INFO(iset)
          READ(5,*) ISOT(iset), BOLTZ(iset), DTYPE(iset), IFRPW(iset),
     1                                                       PQR(iset)
c=======================================================================
          WRITE(6,637) iset,INFO(iset)
          IF(NFS.LE.1) THEN
              CN(iset,1)= 1
            ELSE
              IF(DTYPE(iset).EQ.1) THEN
c=======================================================================
c  If property is sum of intensities for multiple final states, read
c       integer final-state weight coefficients here
c  E(tot)= CN1*E(fs1) + CN2*E(fs2) + ... + CNnfs*E(nfs); CNi = 1 or 0
c=======================================================================
                  READ(5,*) (CN(iset,ifs), ifs= 1,NFS)
c=======================================================================
                  WRITE(6,629) NFS,(CN(iset,ifs),ifs= 1,NFS)
                ELSEIF(DTYPE(iset).EQ.2) THEN
c=======================================================================
c  If property is ratio of intensities, read (integer) numerator (CN)
c         and denominator (CD) final-state sum weight coefficients here.
c     E(ratio) = {CN(1)*E(fs1) + CN(2)*E(fs2) + ... + CN(nfs)*E(nfs)}/
c                {CD(1)*E(fs1) + CD(2)*E(fs2) + ... + CD(nfs)*E(nfs)}
c  CN(i) = 1 or 0
c  CD(i) = 1 or 0
c=======================================================================
                  READ(5,*) (CN(iset,ifs), ifs= 1,NFS)
                  READ(5,*) (CD(iset,ifs), ifs= 1,NFS)
c=======================================================================
                  WRITE(6,630) iset,(CN(iset,ifs),ifs= 1,NFS)
                  WRITE(6,631) (CD(iset,ifs),ifs= 1,NFS)
                ELSEIF((DTYPE(iset).LE.0).OR.(DTYPE(iset).GE.3)) THEN
                  WRITE(6,633) iset,DTYPE(iset)
                  STOP
                ENDIF
            ENDIF
          TEMP(iset)= 0.d0
          JFRPW= IABS(IFRPW(iset))
          IF((JFRPW.NE.0).AND.(JFRPW.NE.1).AND.(JFRPW.NE.3)) THEN
              WRITE (6,604) IFRPW(iset)
              STOP
              ENDIF
          IF(JFRPW.EQ.1) WRITE(6,632)
          IF(JFRPW.EQ.3) WRITE(6,634)
          BFCT(ISOT(iset))= MU(ISOT(iset))*RH*RH/16.85762908d0
          IF(JFRPW.EQ.0) FACTOR(iset)= 9.17555390D10*SQRTMU(ISOT(iset))
          IF(JFRPW.EQ.1) FACTOR(iset)= 8.4397201D0*SQRTMU(ISOT(iset))
          IF(JFRPW.EQ.3) FACTOR(iset)= 2.4313849D-8*SQRTMU(ISOT(iset))
          IF(PQR(iset).GT.0) WRITE(6,613)
          IF(PQR(iset).LE.0) WRITE(6,615)
c
c** Data/control parameter input for fit/prediction of predissociation
c-----------------------------------------------------------------------
          IF(IFRPW(iset).EQ.0) THEN
              BOLTZ(iset)= 0
              NJ(iset)= 0
c=======================================================================
c  NVJ  is number of (v,J) states for which predissociation calculated
c    In forward calculation, if NFV=0, perform for specified (v,J) range
c=======================================================================
              READ(5,*) NVJ(iset)
c=======================================================================
              IF(FITIT.LE.0) THEN
                  IF(NVJ(iset).GT.0) THEN
c** For No-FIT predictions: either read in (v,J)'s for predissociating
c   levels (if  NVJ > 0), ... OR (if NVJ=0) read integers to specify a
c   range of initial-state levels V1ST.le.v.le.VMAX; J1ST.le.J.le.JMAX
c=======================================================================
c If  NVJ > 0  read in list of initial (v,J) levels.
c   VFIX - specific value of v for level (v,J) in forward calculation
c   JFIX - specific value of J for level (v,J)
c=======================================================================
                      READ(5,*) (VFIX(ivj), JFIX(ivj), ivj= 1,NVJ(iset))
c=======================================================================
                      VMAX(iset)= 0
                      DO  ivj= 1, NVJ(iset)
                          VMAX(iset)= MAX0(VMAX(iset),VFIX(ivj))
                          ENDDO
                      WRITE(6,645) NVJ(iset),iset,NAME1(ISOT(iset)),
     1                MN1(ISOT(iset)),NAME2(ISOT(iset)),MN2(ISOT(iset)),
     2                (VFIX(ivj),JFIX(ivj),ivj= 1,NVJ(iset))
                    ELSEIF(NVJ(iset).LE.0) THEN
c=======================================================================
c If read-in  NVJ=0 ,  calculate predissociation rates for all levels
c   from  v= V1ST  up to  VMAX,  for  J= J1ST  up to  JMAX
c=======================================================================
                      READ(5,*) V1ST(iset),VMAX(iset), J1ST(iset),
     1                                                      JMAX(iset)
c=======================================================================
                      WRITE(6,649) V1ST(iset),J1ST(iset),VMAX(iset),
     1                                                      JMAX(iset)
                    ENDIF
c
                ELSEIF(FITIT.GT.0) THEN
                  IF(NVJ(iset).LE.0) THEN
                      WRITE(6,603) iset,NVJ(iset),FITIT
                      STOP
                      ENDIF
c=======================================================================
c  PUNITS > 0 for predissociation data input as widths (FWHM in cm-1)
c         = 0 for predissociation data input as lifetimes (seconds)   
c         < 0 for predissociation data input as rates (inverse seconds)
c  [Predissociation internally treated as rates (in s-1), so convert
c  values if in other units.  For not fitting, PUNITS is dummy variable]
c       FWHM(cm-1) = RATE(s-1)/2*pi*c = 1/TAU(s)*2*pi*c
c       ** 1/(2*pi*c) = 0.5308837458D-11               
c=======================================================================
                  READ(5,*) PUNITS(iset)
c=======================================================================
                  ivj= 1
                  DO  i= 1,nvj(iset)
c=======================================================================
c For a FIT, need to read each predissociation datum  OBS  and its
c   uncertainty  UNC  for each initial-state level  (v,J)= (VFIX,JFIX)
c   Note: PUNITS (above) specifies units of input OBS and UNC values;
c   calculation done as rates (s-1) so OBS and UNC converted internally
c=======================================================================
                      READ(5,*) VFIX(ivj), JFIX(ivj), OBS(ivj,iset),
     1                                                   UNC(ivj,iset)
c=======================================================================
                      IF(UNC(ivj,iset).LT.UCUTOFF) ivj= ivj+1
                      ENDDO
                      nvj(iset)= ivj-1
                  IF(PUNITS(iset).LT.0) THEN
                      WRITE(6,643) NVJ(iset),iset,(VFIX(ivj),JFIX(ivj),
     1                   OBS(ivj,iset),UNC(ivj,iset),ivj= 1,NVJ(iset))
                    ELSEIF(PUNITS(iset).EQ.0) THEN
                      WRITE(6,647) NVJ(iset),iset,(VFIX(ivj),JFIX(ivj),
     1                   OBS(ivj,iset),UNC(ivj,iset),ivj= 1,NVJ(iset))
                      DO  ivj= 1,NVJ(iset)
                          OBS(ivj,iset)= 1.d0/OBS(ivj,iset)
                          UNC(ivj,iset)= 1.d0/UNC(ivj,iset)
                          ENDDO
                    ELSEIF(PUNITS(iset).GT.0) THEN
                      WRITE(6,648) NVJ(iset),iset,(VFIX(ivj),JFIX(ivj),
     1                   OBS(ivj,iset),UNC(ivj,iset),ivj= 1,NVJ(iset))
                      DO ivj= 1,NVJ(iset)
                          OBS(ivj,iset)= TWOPIC*OBS(ivj,iset)
                          UNC(ivj,iset)= TWOPIC*UNC(ivj,iset)
                          ENDDO
                    ENDIF
                  VMAX(iset)= 0
cc                JMAX(iset)= 0
                  DO  ivj= 1,nvj(iset)
c ... now - finally save OBS in correct units for fitting.
                      NDATA= NDATA+1
                      YO(NDATA)= OBS(ivj,iset)
                      YU(NDATA)= UNC(ivj,iset)
                      VMAX(iset)= MAX(VMAX(iset),VFIX(ivj))
cc                    JMAX(iset)= MAX(JMAX(iset),JFIX(ivj))
                      ENDDO
                      NFREQ(iset)= NVJ(iset)
                ENDIF
c??? (not clear why Geoff set these ... but ...
              NFREQ(iset)= 1
              FREQ(1,iset)= 0.d0
              WAVL(1,iset)= 0.d0
              ENDIF
c ... end of input for case of a predissociation dataset
c
          IF(IFRPW(iset).NE.0) THEN
c** Case of absorption or emission data sets for either thermal initial
c    state or transitions (normally emission) from one initial (v,J)
              IF(BOLTZ(iset).GT.0) THEN         
c-----------------------------------------------------------------------
c** Data/control parameter input for a "thermal" property:  BOLTZ > 0
c=======================================================================
c  TEMP  is the Kelvin temperature of the current input data set
c  VMAX  is the upper bound cutoff v value for the thermal sum over v
c  NJ  specifies how sum over a thermal rotational population is done
c      < 0 for direct sum from J= 0 to J= |NJ|
c      = 0 perform all calculations with J= 0
c      > 0 sums over nj average J values in nj equally weighted
c          segments of the rotational population for that v
c=======================================================================
                  READ(5,*) TEMP(iset), VMAX(iset), NJ(iset)
c=======================================================================
                  WRITE(6,625) TEMP(iset),NAME1(ISOT(iset)),
     1    MN1(ISOT(iset)),NAME2(ISOT(iset)),MN2(ISOT(iset)),VMAX(iset)
                  WRITE(6,624) POPCRT
                  IF(NJ(iset).GT.mxnj) THEN
                      WRITE(6,611) mxnj,NJ(iset),NJ(iset)
                      NJ(iset)= -NJ(iset)
                      ENDIF
                  IF(NJ(iset).LE.0) THEN
                      JMAX(iset)= IABS(NJ(iset))
                      IF(JMAX(iset).GT.0) WRITE(6,618) JMAX(iset)
                      IF(JMAX(iset).EQ.0) WRITE(6,620)
                    ELSE
                      JMAX(iset)= NJ(iset)-1
                      WRITE(6,622) NJ(iset)
                    ENDIF
                ELSEIF(BOLTZ(iset).LE.0) THEN
c** Control parameters for non-thermal absorption/emission from the
c  particular initial-state level:  v= V1ST, J= J1ST
c=======================================================================
                  READ(5,*) V1ST(iset), J1ST(iset)
c=======================================================================
                  WRITE(6,623) NAME1(ISOT(iset)),MN1(ISOT(iset)),
     1       NAME2(ISOT(iset)),MN2(ISOT(iset)),V1ST(iset),J1ST(iset)
                  VMAX(iset)= V1ST(iset)
                  JMAX(iset)= J1ST(iset)
                  NVJ(iset)= 0
                ENDIF
c
              IF(FITIT.LE.0) THEN
c** For Non-FIT forward calculation, specify frequencies for calculations
c=======================================================================
c* For forward calculation, specify desired transition energy mesh
c  NFREQ  is the number of transition energies (in cm-1) for this set
c  FREQ(i)= FREQ1 + (i-1)*DFREQ    ; values negative for emission
c=======================================================================
                  READ(5,*) NFREQ(iset), FREQ1, DFREQ
c=======================================================================
                  IF(NFREQ(iset).GT.mxfreq) THEN
c** If number of requested frequencies exceeds array size - STOP!)
                        WRITE(6,888)
  888 FORMAT(/' *** ERROR *** input  NFREQ=',i10,'  > array dimension',
     1  '   mxfreq=',i6)
                        STOP
                        ENDIF
                  DO ifr= 1,NFREQ(iset)
                      FREQ(ifr,iset)= FREQ1 + (ifr-1)*DFREQ
                      WAVL(ifr,iset)= 1.d7/FREQ(ifr,iset)
                      ENDDO
                  WRITE(6,640) NFREQ(iset),(FREQ(ifr,iset),
     1                                             ifr= 1,NFREQ(iset))
                  ENDIF
c
              IF(FITIT.GT.0) THEN
c** For a FIT, read experimental absorption/emiss. intensities or ratios
c=======================================================================
c NFREQ    is the number of data to be read in the current set
c FREQYN > 0  if read-in ordinate values are energies in (cm-1) 
c        = 0  if read-in ordinate values are wavelengths in (nm) 
c=======================================================================
                  READ(5,*) NFREQ(iset), FREQYN
c=======================================================================
                  IF(NFREQ(iset).GT.mxfreq) THEN
                      WRITE(6,616) mxfreq
                      STOP
                      ENDIF
                  ifr= 1
                  DO  i= 1, NFREQ(iset)
c=======================================================================
c  FREQ  input transition energies (cm-1) [or wavelengths (nm)]
c        (for emission, values should be negative)
c  OBS   experimentally observed intensity values for current set
c  UNC   uncertainties in the experimental values
c=======================================================================
                      READ(5,*) FREQ(ifr,iset), OBS(ifr,iset),
     1                                                   UNC(ifr,iset)
c=======================================================================
                      IF(UNC(ifr,iset).LT.UCUTOFF) ifr= ifr+1
                      ENDDO
                      NFREQ(iset)= ifr-1
                  DO ifr= 1,NFREQ(iset)
                      NDATA= NDATA+1
                      YO(NDATA)= OBS(ifr,iset)
                      YU(NDATA)= UNC(ifr,iset)
                      ENDDO
                  IF(FREQYN.GT.0) THEN
c** As required ... generate wavelengths/frequencies for calculation
                      DO ifr= 1,NFREQ(iset)
                          WAVL(ifr,iset)= 1.d7/FREQ(ifr,iset)
                          ENDDO
                    ELSE
                      DO ifr= 1,NFREQ(iset)
                          WAVL(ifr,iset)= FREQ(ifr,iset)
                          FREQ(ifr,iset)= 1.d7/WAVL(ifr,iset)
                          ENDDO
                    ENDIF
                  WRITE(6,638) NFREQ(iset),iset,(FREQ(ifr,iset),
     1  WAVL(ifr,iset),OBS(ifr,iset),UNC(ifr,iset),ifr= 1,NFREQ(iset))
                  ENDIF
              ENDIF
          IF(IFRPW(iset).LT.0) WRITE(6,644) JFRPW,FREQ(1,iset),JFRPW
   10     CONTINUE
c=======================================================================
c ... end of loop to input data/cases for performing calculations
c=======================================================================
      DO iset= 1,NSETS
          SCALE(iset)= 0
          SF(iset)= 1.d0
          ENDDO
      IF((FITIT.GT.0).AND.(NSETS.GT.1)) THEN
c
c** For a fit to multiple data sets, may wish to apply a global scaling
c  to certain data subsets to allow for (say) uncertainties in
c  concentration measurements for different isotopomers.
c=======================================================================
c  SF(iset)  is a multiplicative factor for scaling 2'nd, 3'rd, ... etc.
c     data set values relative to 1'st data set.
c  SCALE is a integer flag that tells whether or not scaling factors for
c     2'nd, 3'rd, ...) data sets are to be held fixed or varied in a fit.
c      = 0  does not vary scaling factor for the set
c      = 1  all scaling factor to vary in the fit.
c=======================================================================
          READ(5,*) (SF(iset), iset= 2, NSETS)
          READ(5,*) (SCALE(iset), iset= 2, NSETS)
c=======================================================================
          ENDIF
c=======================================================================
c              Preparation of Initial-State Potential
c=======================================================================
c  Prepare distance array  R = RAD(i)  defining potential mesh points up
c  to the upper dimensioned limit mxfsp
c  Also prepare array XM2(i) = 1/R^2
c-----------------------------------------------------------------------
      DO  i= 1,mxfsp
          RAD(i)= RMIN + (i-1)*RH
          XM2(i)= 1.d0/RAD(i)**2
          ENDDO   
c
      WRITE(6,656)
c** Allow different isotopomers to have different initial-state potentials
      DO  20 iso= 1,NISTP
          DO  i= 1,NIN
              RM2(i)= XM2(i)
              ENDDO
          CALL PREPOT(LNPT,AN1,AN2,MN1(iso),MN2(iso),NIN,OMEGA(0),
     1                                         RAD,RM2,VLIMI,VIS,NCNI)
c-----------------------------------------------------------------------
c  prepot inputs:  LNPT,MASS1(iso),MASS2(iso),VLIMI,NIN,RAD(i),RM2(i)
c  prepot outputs: VIS(i) - the generated initial-state potential array
c                          in units cm-1
c                  NCNI - read inside subroutine regarding potential tail
c                  RM2(i) - array may change to account for corrections
c=======================================================================
c in PREPOT:  READ(5,*) NTPI, LPPOTI, OMEGA(0), VLIMI
c=======================================================================
c  NTPI    > 0 if reading in potential points and interpolating over
c             and extrapolating beyond them to get desired potential
c          = 0 to generate an analytic potential to be specified below
c  LPPOTI  > 0 to print generated potential and its derivatives-by-
c             differences at every |LPPOTI|th point
c          = 0 to prevent such printing
c  OMEGA  the (integer) total electronic angular momentum projection
c         quantum number (required for proper rotational intensities)
c  VLIMI   - the absolute energy (in cm-1) at large R asymptote
c=======================================================================
c  For pointwise potential (NTPI > 0)
c in PREPOT:  READ(5,*) NUSEI, IR2I, ILRI, NCNI, CNNI
c in PREPOT:  READ(5,*) RFACTI, EFACTI, VSHIFTI
c in PREPOT:  READ(5,*) (XI(I), YI(I),I=1,NTPI)
c=======================================================================
c  NUSEI > 0 to use NUSEI-point piecewise polynomial interpolation for
c            read in potential points (typically 6,8,or 10)
c        = 0 to use cubic spline interpolation scheme
c
c  IR2I  > 0 causes numerical interpolation to be performed over YI*XI^2
c            rather than over the potential YI themselves; this may
c            improve interpolation for a steep repulsive wall.
c
c  ILRI  specifies how to extrapolate beyond outer turning points
c        < 0 fits last 3 to:  VLIMI - A exp[-b(R-R0)^2]
c        = 0 fits last 3 to:  VLIMI - A R^p exp(-bR)
c        = 1 fits last 2 to:  VLIMI - A/(R^B)
c        = 2 or 3  fit last 2 or 3 points to:
c                     VLIMI - sum[C(NCNI+2m)/R^(NCNI+2m)]; m = 0,ILRI-1
c        > 3  fits last few points to:
c                     VLIMI - sum[C(NCNI+m)/R^(NCNI+m)]; m = 0,ILRI-1
c
c  NCNI  is used (if ILR > 1) to specify limiting inverse-power behavior
c        of function tail;  see ILRI.
c
c  For ILRI > 1 cases, if CNNI.NE.0 fix limiting coefficient C(NCNI) at
c  this read-in value, rather than from fit to outermost turning points.
c
c  RFACTI - factor to convert read-in turning point x-values to angstroms
c  EFACTI - factor to convert read-in turning point y-values to cm-1
c  VSHIFTI - the energy (in cm-1) which must be added to the read-in
c         potential points to make then consistent with the VLIMI value.
c
c  (Xi,Yi) - are the NTPI pairs of potential points
c=======================================================================
c  Generate analytical initial-state potential (NTPI.LE.0) in PREPOT
c* Potentials generated in cm-1 with potential asymptote at energy VLIM
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** IPOTL specifies the type of potential function to be generated.
c** MPAR & NPAR are integers for specifying potential types.
c** NVARB is number of (real*8) potential parameters read in.
c** IBOB specifies whether (if > 0) or not (if .le. 0) atomic mass
c      dependent Born-Oppenheimer breakdown corrections will be included
c** For all functions considered, well depth and equilibrium distance
c  are read as  DSCM (cm-1)  and REQ (Angstroms), respectively.
c* [Most read-in parameters are dimensionless (scaled by DSCM & REQ).]
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** If IPOTL=1  generate an L.J.(MPAR,NPAR) potential.
c** If IPOTL=2  generate an MLJ(NPAR) potential [JCP 112, 3949 (2000)]
c      If MPAR.ge.0  exponent parameter is polynomial of order (NVARB-1)
c           in z=(R-Re)/(R+Re), with the NVARB coefficients  PARM(j)
c      If MPAR < 0  exponent polynomial in  z  has order (NVARB-4) with
c           coefficients PARM(i) (i= 1,NVARB-3), & includes a switching
c           function with exponent coefficient  ALPHA= PARM(NVARB)  and
c           RSW= PARM(NVARB-1),  defined to yield limiting inverse-power
c           potential coefficient  Cn= PARM(NVARB-2).
c** If IPOTL=3  generate a Morse or Extended Morse Oscillator potential
c       with exponent factor "beta" defined as a power series of order
c       (NVARB-1) in  z=(R-Re)/(R+Re)  with NVARB coefficients PARM(i).
c       Set  NVARB= 1 for conventional "simple" Morse potential.
c*   Special option #1: set  MPAR= -1  to produce Wei Hua's 4-parameter
c      modified Morse function with  b= PARM(1)  and C= PARM(2).
c*   Special option #2: set  MPAR= -2  to produce Coxon's "Generalized
c      Morse Oscillator" potential with exponent expansion in (R-Re)]
c ...  otherwise, set  MPAR.ge.0
c** If IPOTL=4  use Seto's modification of Surkus' GPEF expansion in
c       z = [R^NPAR - Re^NPAR]/[a*R^NPAR + b*Re^NPAR] where
c       a=PARM(NVARB-1) & b=PARM(NVARB), which incorporates Dunham, SPF,
c       O-T and other forms: V(z) = c_0 z^2 [1 + c_1 z + c_2 z^2 + ...]
c       where  c_0 [cm-1] is read in as DSCM, and the first (NVARB-2)
c       PARM(i)'s are the  c_i  (i > 0).  [MPAR is dummy parameter here]
c  * For Dunham case:  NPAR=1, PARM(NVARB-1)= 0.0, PARM(NVARB)= 1.0
c  * For SPF case:  NPAR=1, PARM(NVARB-1)= 1.0, PARM(NVARB)= 0.0
c  * For Ogilvie-Tipping:  NPAR=1, PARM(NVARB-1)= 0.5 = PARM(NVARB)
c  * NOTE that for Surkus NPAR < 0 case:  z(NPAR,a,b)= z(|NPAR|,-b,-a)
c      Generate & return the  D_e  value implied by these coefficients.
c** If IPOTL=5  generate generalized HFD(NPAR,6,8,10,12,14) potential.
c       PARM(1-3) are the parameters defining the HFD damping function
c       D(x)=exp[-pparm(1)*(PARM(2)/x - 1)**PARM(3)] {for x < PARM(2)}
c       PARM(4) the quadratic coefficient in the exponent, and
c       PARM(5) is the power of  x=R/Req  multiplying the repulsive term
c              AREP*x**PARM(5) *exp[-beta*x - PARM(4)*x**2] ;
c       PARM(6-11)  are the reduced C_NPAR, C_6, C_8, C_10, C_12 and C14
c       parameters (NPAR < 6), while  AREP and  beta  are defined
c       by having the potential minimum at  x=1.  For NVARB < 11, higher
c       C_m coefficients automatically zero;  necessarily  NVARB.ge.7 .
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** IBOB > 0, add atomic-mass-dependent Born-Openheimer breakdown
c  correction functions to rotationless and/or centrifugal potential(s).
c  Both expressed as power series in  z= (R-Re)/(R+Re) starting with the
c  constant term, using the mass shift convention of Le Roy [J.Mol.Spec.
c  194, 189 (1999)].  Adiabatic B-O-B potential correction fx. defined
c  by polynomials of order NC1 with (NC1+1) coefficients {CA1(i)} for
c  atom-1  and order NC2 with (NC2+1) coefficients {CA2(i)} for atom-2,
c  while centrifugal correction fx. defined polynomial of order NG1 with
c  (NG1+1) coefficients {GA1(i)} for atom-1 a nd order NG2 with (NG2+1)
c  coefficients {GA2(i)} for atom-2.
c** Input parameters IANi & IMNi are the atomic & mass number of atom-i
c  (i=1,2), while integers RMN1 & RMN2 read here are the mass numbers of
c  the reference isotopes defining the B-O-B correction functions.
c** NC1 & NC2 are orders of polynomials DELTA(V,atom-i) defining
c  'adiabatic' corrections to the rotationless potential for atoms 1 & 2
c  DELTA(V)= (1-M1ref/M1)*DELTA(V,atom-1) + (1-M2ref/M2)*DELTA(V,atom-2)
c** NG1 & NG2 are orders of polynomials q1(z) & q2(z) defining B-O-B
c   correction to the centrifugal potential:
c      V(centrifugal)= [1 + (M1ref/M1)*q1(z) + (M2ref/M2)*q2(z)]/R**2
c ... to omit a particular correction set associated NCi or NGi .lt.0
c** RX > 0.0  invokes Coxon's (older) expansions in (R-Re) for potential
c     correction and in  [(R-Rx)**j - (Re-Rx)**j] for centrifugal corrn.
c ... OTHERWISE (to use Le Roy B-O-B formalism) set  RX.le.0.d0 !!
c-----------------------------------------------------------------------
c** Read inside subroutine POTGEN
c         IF(LNPT.GT.0) THEN
c             READ(5,*) IPOTL, MPAR, NPAR, NVARB, IBOB, DSCM, REQ
c             IF(NVARB.GT.0) READ(5,*) (PARM(I), I=1,NVARB)
c             IF(IBOB.GT.0) THEN
c                 READ(5,*) RMN1, RMN2, NC1, NC2, NG1, NG2, RX
c                 IF(NC1.GE.0) READ(5,*) (CA1(I), I=0,NC1)
c                 IF(NC2.GE.0) READ(5,*) (CA2(I), I=0,NC2)
c                 IF(NG1.GE.0) READ(5,*) (GA1(I), I=0,NG1)
c                 IF(NG2.GE.0) READ(5,*) (GA2(I), I=0,NG2)
c             ENDIF
c         ENDIF
c-----------------------------------------------------------------------
          DO  i= 1,NIN
c* Need to scale VIS(i) potential to internal units before entering ALF
              VIS(i)= BFCT(iso)*VIS(i)
              VI(i,iso)= VIS(i)
              VICD(i,iso)= RM2(i)
              ENDDO
          MAXV= 0
          DO  iset= 1,NSETS
              IF(VMAX(iset).GT.MAXV) MAXV= VMAX(iset)
              ENDDO
c-----------------------------------------------------------------------
c  Automatic Level Finder ALF:
c    inputs:  NIN,RMIN,RH,VIS,VLIMI,VMAX,ZMU,EPS,NCNI
c    outputs: GV(v) calculated eigenvalue array
c         INNR(v) associates level v with (=1) inner vs. outer (=0) well
c         ERR - value representing degree of subroutine success 
c-----------------------------------------------------------------------
          CALL ALF(NIN,RMIN,RH,VIS,WF1,VLIMI,MAXV,ERR,MU(iso),EPS,
     1                                            NCNI,GV,INNR,IWRSCH)
c** Now call CDJOEL to get CDC's needed for eigenvalue predictions
          JREF= 0
          BvWN= RH**2/BFCT(iso)
          DO  v=0,MAXV
              MCI(v,0,iso)= GV(v)
              INNER(v,iso)= INNR(v)
              CALL SCHRQ(v,JREF,GV(v),GAMA,PMAX,VLIMI,VIS,WF1,BFCT(iso),
     1                 EPS,RMIN,RH,NIN,NBEG,NEND,INNR(v),IWRSCH,LPRWF)
              CALL CDJOEL(GV(v),NBEG,NEND,BvWN,RH,IWRSCH,VIS,WF1,RM2,
     1                                                          RCNST)
              DO  c= 1,7
                  MCI(v,c,iso)= RCNST(c)
                  ENDDO
              ENDDO
c-----------------------------------------------------------------------
          WRITE(6,662) NAME1(iso),MN1(iso),NAME2(iso),MN2(iso),MAXV,
     1                             (v,(MCI(v,c,iso),c= 0,5),v= 0,MAXV)
          WRITE(6,699)
   20     CONTINUE
c=======================================================================
c  REXFS  - reference distance about which the final-state potentials
c           are expanded
c  REXTMF - reference point about which the transition moment functions
c           are expanded
c  LPFS > 0 causes every LPFS'th point of the final-state potential
c           arrays to be written to channel-9;  otherwise no print-out.
c  LPTMF > 0 causes every LPFS'th point of transition moment function
c           arrays to be written to channel-9;  otherwise no print-out.
c        For predissociation with derivative operator coupling, the
c        separate dW/dR and dPSIc/dR components are collected and
c        written WHEN FITIT = 0 (forward calc).
c=======================================================================
      READ(5,*) REXFS, REXTMF, LPFS, LPTMF
c=======================================================================
c Preparation of Final-State Potentials and Transition Moment Functions
c-----------------------------------------------------------------------
      WRITE(6,664) NFS
      NPTOT= 0
      NPPFREE= 0
      DO 60 ifs= 1,NFS
c=======================================================================
c* FSTYPE = 1 for exponential repulsive potential
c    VF = VLIMF + A * exp{-(R-REXFS)*[B + C*z + D*z^2 + E*z^3 + F*z^4};
c    Values of A-F read as parameters FSPRM(i,ifs), i= 1-6 respectively
c
c* FSTYPE = 3 for Extended Morse Oscillator potential
c   VF = VLIMF + A1*[exp{-(R-A2)*(A3 + A4*z + ...)} -1]**2 - A1
c   Values of A1-A6 read as parameters FSPRM(i,ifs), i= 1-6 respectively
c
c* FSTYPE = 2 for pot'l defined by NTPFS turning pts. [RTPF(i),VTPF(i)]
c     with a repulsive exponential inner wall attached to the 2
c     innermost points. Interpolate to get potential array at required
c     mesh with NUSEF-point piecewise polynomials or cubic splines.
c     In this option, also add energy VSHFTF (cm-1) to the read-in
c     potential points to make them consistent with the stated VLIMF
c     value (often VSHFTF=Te). Read-in turning point distances and
c     energies are multiplied by RFACT and VFACT, respectively,
c     converting units to Angstroms and cm-1.
c
c     Repulsive inner wall defined by fitting the 2 innermost read-in
c     turning points to the form:
c
c     VF(R) = A + B exp[-(R-REXFS)*(b0 + b1*z + b2*z^2 + ... + b5*z^5)
c
c     where A and B are determined by the first 2 turning points;
c     parameters b(j)= FSPRM[(j+1),nfsprm]
c
c  OMEGA  - is the total electronic angular momentum progection quantum
c           number.  Required when calculating Honl-London factor needed
c           for P/Q/R intensity calculation for option 'PQR(iset)=1'
c
c  NFSPRM - total number of final state parameters, fixed and free
c
c  VLIMF  - absolute energy (cm-1) at potential asymptote (large R)
c
c  XCOORD is integer specifying expansion coordinate for the potential
c          = p (p=1-9)  for  zfs = (R^p - REXFS^p)/(R^p + REXFS^p)
c          = 10 for zfs = (R-REXFS)/R
c          = 11 for zfs = (R-REXFS)/(REXFS)
c=======================================================================
          READ(5,*) FSTYPE(ifs), OMEGA(ifs), NFSPRM(ifs), VLIMF(ifs),
     1                                                      XCOORD(ifs)
c=======================================================================
          IF(NFSPRM(ifs).GT.6) THEN
              NFSPRM(ifs)= mxprm
              WRITE(6,666) NFSPRM(ifs)
              ENDIF
c=======================================================================
c  FSPRM(i) - array of nfsprm parameters defining the final state PES
c=======================================================================
          READ(5,*) (FSPRM(j,ifs), j=1, NFSPRM(ifs))
c=======================================================================
          IF (FITIT.GT.0) THEN
c=======================================================================
c  FSVAR(i) - specifies the variability of fsprm(i)
c              = 1 if parameter allowed to vary
c              = 0 if parameter fixed
c=======================================================================
              READ(5,*) (FSVAR(j,ifs), j=1, NFSPRM(ifs))
c=======================================================================
              DO j= 1,NFSPRM(ifs)
                  IF((DABS(FSPRM(j,ifs)).LT.1.d-5).AND.
     1                                       (FSVAR(j,ifs).GE.1)) THEN
                      FSPRM(j,ifs)= 1.d-3
                      WRITE(6,698) j,ifs
                      ENDIF
                  ENDDO
              ENDIF
          WRITE(6,667) ifs,OMEGA(ifs)
          IP= XCOORD(ifs)
          IF(XCOORD(ifs).EQ.1) WRITE(6,668) REXFS
          IF((XCOORD(ifs).GE.2).AND.(XCOORD(ifs).LE.9))
     1                                  WRITE(6,669) IP,IP,IP,IP,REXFS
          IF(XCOORD(ifs).EQ.10) WRITE(6,670) REXFS
          IF(XCOORD(ifs).EQ.11) WRITE(6,671) REXFS
          IF((XCOORD(ifs).LE.0).OR.(XCOORD(ifs).GE.12)) THEN
              WRITE(6,672) XCOORD(ifs)
              ENDIF
          IF(FSTYPE(ifs).EQ.1) THEN
              CALL GENFS(ifs)
              WRITE(6,674) VLIMF(ifs),(i,FSPRM(i,ifs),i=1,NFSPRM(ifs))
            ELSEIF(FSTYPE(ifs).GE.3) THEN
              CALL GENFS(ifs)
              WRITE(6,675) VLIMF(ifs),(i,FSPRM(i,ifs),i=1,NFSPRM(ifs))
c-----------------------------------------------------------------------
c  GENFS has ability to calculate and re-calculate new final state
c  potentials for fitting iterations or forward calculation
c-----------------------------------------------------------------------
            ELSEIF(FSTYPE(ifs).EQ.2) THEN
c=======================================================================
c  NTPFS   - number of data points KNOWN for potential:  see description
c        in PREPOT (above) for meanings of NUSEF, IR2F, ILRF, NCNF, CNNF
c  RFACTF  - conversion factor so that radial array is in Angstroms
c  VFACTF  - conversion factor so that potential array is in cm-1
c  VSHFTF  - the energy (cm-1) which must be added to the potential to
c            make it consistent with the chosen VLIMF value
c  RTPF(i) - radial array of turning points
c  VTPF(i) - potential array of turning points
c=======================================================================
              READ(5,*) NTPFS(ifs)
              READ(5,*) NUSEF, IR2F, ILRF, NCNF, CNNF
              READ(5,*) RFACTF, VFACTF, VSHFTF(ifs)
              READ(5,*) (RTPF(i), VTPF(i), i= 1,NTPFS(ifs))
c=======================================================================
              IF(NUSEF.GT.0) WRITE(6,676) VLIMF(ifs),NUSEF,NTPFS(ifs)
              IF(NUSEF.LE.0) WRITE(6,677) VLIMF(ifs),NTPFS(ifs)
              IF(IR2F.GT.0) WRITE(6,678)
              IF((ILRF.GT.1).AND.(DABS(CNNF).GT.0.D0))
     1                                          WRITE(6,679) CNNF,NCNF
              WRITE(6,680) VSHFTF(ifs),RFACTF,VFACTF,
     1                                (RTPF(i),VTPF(i),i=1,NTPFS(ifs))
              WRITE(6,681)
              DO i= 1,NTPFS(ifs)
                  RTPF(i)= RTPF(i)*RFACTF
                  VTPF(i)= VTPF(i)*VFACTF + VSHFTF(ifs)
                  ENDDO
c-----------------------------------------------------------------------
c  RTP(i) and VTP(i) arrays are 1st and 2nd radial turning points
c  needed later in determination of inner repulsive wall.
c-----------------------------------------------------------------------
              DO i=1,2
                  RTP(i,ifs)= RTPF(i)
                  VTP(i,ifs)= VTPF(i)
                  ENDDO
              IF(IR2F.GT.0) THEN
                  DO  i= 1,NTPFS(ifs)
                      VTPF(i)= VTPF(i)*RTPF(i)**2
                      ENDDO
                  ENDIF
              NPRF= mxfsp
              CALL GENINT(LNPT,mxfsp,RAD,FCT,NUSEF,IR2F,NTPFS(ifs),
     1                  RTPF,VTPF,VLIMF(ifs),ILRF,NCNF,CNNF,NPRS,NPRF)
              WRITE(6,681)
c-----------------------------------------------------------------------
c  GENINT subroutine used to interpolate over the entire range and
c  extrapolate for both inner and outer regions of final state 2.  The
c  inner region will later be overwritten as it is fit to experimental
c  data.  Points found using GENINT will remain fixed for the mesh from
c  RTPF(2) and outward but inner points must be variable (to fit data);
c  mesh here will be determined by the added repulsive inner wall.
c
c  GENINT inputs: LNPT, mxfsp, RAD(i), NUSEF, IR2F, {RTPF(i),VTPF(i)},
c             VLIMF(ifs), ILRF, NCNF, CNNF, NPRS, NPRF
cgtk  NPRS       - no purpose here (to be removed ?????)
cgtk  NPRF       - no purpose here (to be removed ?????)
c
c  GENINT outputs:  FCT is 1-dimensional final-state potential array
c-----------------------------------------------------------------------
c  Now need to get Final State array into 2 dimensional form. 
c-----------------------------------------------------------------------
              DO i=1,mxfsp
                  VF(i,ifs)= FCT(i)
                  ENDDO
              CALL GENFS(ifs)
              IF(FSTYPE(ifs).EQ.2)  WRITE(6,601) AVALUE(ifs),
     1              BVALUE(ifs),(j,FSPRM(j,ifs),j=1,NFSPRM(ifs))
              ENDIF
c-----------------------------------------------------------------------
c  Prepare Transition Moment Function arrays here
c=======================================================================
c  GFS  is a positive  integer defining the electronic degeneracy of the
c       transition associated with this state.
c  TMFTYP specifies coordinate ztmf(R) in power series expansion
c         < 0 for predisociation case where the TMF is not a power
c             series but an operator instead
c                         P= -hbar^2/2*mu {dW/dr + 2W*d/dr}
c                         where W(r)= a/{4a^2 + (r - Rc)^2}
c             [TWO CASES: -1 for one Lorentzian, -2 for two Lorentzians]
c            * for the 1 Lorentzian case read in a, Rc as tmf parameters
c            * for the 2 Lorentzian case read in as a1, Rc1, a2, Rc2
c         = 0  for TMF defined by read-in array of points
c         = 1  for  ztmf = (r - Rextmf)/(r + Rextmf)
c         = p = 2-10  for  ztmf = (r^p - Rextmf^p)/(r^p + Rextmf^p)
c         = 11  for  ztmf = (r - Rextmf)/Rextmf      [Dunham]
c         = 12  for  ztmf = 1/r^2
c         = 13  for  ztmf = r         [the distance R itself]
c
c  OTMF   - the order transition moment function power series in {ztmf}
c         - this is a DUMMY value for TMFTYP < 0
c=======================================================================
          READ(5,*) GFS(ifs), TMFTYP(ifs), OTMF(ifs)
c=======================================================================
          IF(TMFTYP(ifs).LT.0) OTMF(ifs)= 2*IABS(TMFTYP(ifs)) - 1
          IF(OTMF(ifs)-1.GT.mxprm) THEN
              WRITE(6,682) OTMF(ifs),mxprm-1
              STOP
              ENDIF
c=======================================================================
c  TMFPRM - coefficients of the transition moment function power series
c         - note that for TMFTYP = 0, these coefficients create a power
c           series in function defined by read-in points (usually want
c           linear function with leading coefficient fixed at zero)
c         - for TMFTYP = -1, enter the "a" and "Rc" as
c           parameters 0 and 1, respectively
c         - for TMFTYP = -2, enter the "a" and "Rc" values for the 1st
c           Lorentzian followed by the "a" and "Rc" values for the 2nd
c=======================================================================
          READ(5,*) (TMFPRM(m,ifs), m=0, OTMF(ifs))
c=======================================================================
          IF(FITIT.GT.0) THEN
c=======================================================================
c  TMFVAR(i) - specifies whether the TMFPRM(i) are fixed or free in fit
c               = 1 if parameter allowed to vary
c               = 0 if parameter fixed
c=======================================================================
              READ(5,*) (TMFVAR(m,ifs), m=0, OTMF(ifs))
c=======================================================================
            ELSE
              DO  m= 0,OTMF(ifs)
                  TMFVAR(m,ifs)= 0
                  ENDDO
            ENDIF
c=======================================================================
c  GENTMF SUMMARY
c  inputs - ifs, TMFTYP
c  outputs  ztmf(i,ifs) - radial array for transition moment function
c                         when TMFTYP .ge. 0
c                       - transition moment array itself when TMFTYP < 0
c=======================================================================
c  IF TMFTYP(ifs)= 0 , GENTMF implements the following reads statements:
c       READ(5,*) NPTMF, TMFLIM
c       READ(5,*) NUSETMF, ILRTMF, NCNTMF, CNNTMF
c       READ(5,*) RFACTMF, MFACTMF
c       READ(5,*) (Xi(i),Yi(i),i=1,NPTMF)
c=======================================================================
c     RFACTMF - factor to convert read-in distances to angstroms
c     MFACTMF - factor to convert read-in y-values to moment fct. units
c=======================================================================
          IF(TMFTYP(ifs).GE.0) THEN
              CALL GENTMF(ifs)
              WRITE(6,686)(TMFPRM(i,ifs),i=0,OTMF(ifs))
              ENDIF
          IF(TMFTYP(ifs).EQ.-1) WRITE(6,695) (TMFPRM(i,ifs),i= 0,1)
c
          IF(FITIT.GT.0) THEN
              DO i=1,NFSPRM(ifs)
                  IF(FSVAR(i,ifs).GT.0) THEN
                      NPPFREE= NPPFREE+ 1
                      NPTOT= NPTOT + 1
                      PV(NPTOT)= FSPRM(i,ifs)
                      ENDIF
                  ENDDO
              DO m= 0,OTMF(ifs)
                  IF(TMFVAR(m,ifs).GT.0) THEN
                      NPTOT= NPTOT + 1
                      PV(NPTOT)= TMFPRM(m,ifs)
                      ENDIF
                  ENDDO
              ENDIF
   60     CONTINUE
      WRITE(6,654) (GFS(i),i= 1,NFS)
      WRITE(6,687)
      DO ifs= 1,NFS
          WRITE(6,688) (i,ifs,FSPRM(i,ifs),FSVAR(i,ifs),i=1,NFSPRM(ifs))
          ENDDO
      WRITE(6,689)
      DO ifs= 1,NFS
          WRITE(6,690) (m,ifs,TMFPRM(m,ifs),TMFVAR(m,ifs),m=0,OTMF(ifs))
          ENDDO
      IF((FITIT.GT.0).AND.(NSETS.GT.1)) THEN
          WRITE(6,691)
          DO iset= 1,NSETS
              IF(SCALE(iset).GE.1) THEN
                  NPTOT= NPTOT + 1
                  PV(NPTOT)= SF(iset)
                  ENDIF
              WRITE(6,692) iset,SF(iset),SCALE(iset)
              ENDDO
          ENDIF
      WRITE(6,699)
      IF(FITIT.LE.0) THEN
          printyn= 1
          CALL FORWARD
          ENDIF
      IF(FITIT.GT.0) THEN
          printyn= 0
          CALL NLLSSRR(NDATA,NPTOT,mxnp,IROUND,NGPRND,LPRINT,YO,YU,YD,
     1                                    PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c-----------------------------------------------------------------------
c  nllssrr input  - ndata,nptot,mxnp,iround,lprint,yo,yu,pv
c  nllssrr output - yd,pv,pu,ps,cm,tstps,tstpu,dse
c
c  write out the correlation matrix below
c-----------------------------------------------------------------------
          IF(NPTOT.GT.1) THEN
              WRITE(6,693) CM(1,1)
              DO i= 2,NPTOT
                  WRITE(6,694) i,(CM(i,k),k= 1,i)
                  ENDDO
              ENDIF
c-----------------------------------------------------------------------
c Before doing final calculation, map final parameters from nllssrr
c (pv values) back onto internal variables of bcont
c-----------------------------------------------------------------------
          nprm= 0
          DO ifs= 1,NFS
              DO n= 1,NFSPRM(ifs)
                  IF(FSVAR(n,ifs).GT.0) THEN
                      nprm= nprm + 1
                      FSPRM(n,ifs)= PV(nprm)
                      ENDIF
                  ENDDO
              DO m= 0,OTMF(ifs)
                  IF(TMFVAR(m,ifs).GT.0) THEN
                      nprm= nprm + 1
                      TMFPRM(m,ifs)= PV(nprm)
                      ENDIF
                  ENDDO
              CALL GENFS(ifs)
              IF(FSTYPE(ifs).EQ.2)  WRITE(6,601) AVALUE(ifs),
     1                   BVALUE(ifs),(j,FSPRM(j,ifs),j= 1,NFSPRM(ifs))
              ENDDO
          DO iset= 1,NSETS
              IF(SCALE(iset).EQ.1) THEN
                  nprm= nprm + 1
                  SF(iset)= PV(nprm)
                  ENDIF
              ENDDO
c-----------------------------------------------------------------------
c Now do final calculation of intensities
c-----------------------------------------------------------------------
          printyn= 1
          CALL FORWARD
          ENDIF
      IF(FITIT.GT.0) THEN
          i= 0
          DO iset= 1,NSETS
              RMS(iset)= 0.d0
              IF(BOLTZ(iset).GT.0) THEN
c                 WRITE(8,802) iset,INFO(iset)
                  DO ifr= 1,NFREQ(iset)
                      i= i + 1
                      ADD= (OUTPUT(ifr,iset)-YO(i))/YU(i)
                      ADD= ADD*ADD
                      RMS(iset)= RMS(iset) + ADD
c                     WRITE(8,803) FREQ(ifr,iset),YO(i),YU(i),
c    1                        OUTPUT(ifr,iset),OUTPUT(ifr,iset)-YO(i),
c    2                        (OUTPUT(ifr,iset)-YO(i))/YU(i)
                      ENDDO
                ELSE
c                 WRITE(8,804) iset,INFO(iset)
                  DO ivj= 1,NVJ(iset)
                      i= i + 1
                      ADD= (OUTPUT(ivj,iset)-YO(i))/YU(i)
                      ADD= ADD*ADD
                      RMS(iset)= RMS(iset) + ADD
c                     WRITE(8,805) VFIX(ivj),JFIX(ivj),YO(i),YU(i),
c    1                        OUTPUT(ivj,iset),OUTPUT(ivj,iset)-YO(i),
c    2                        (OUTPUT(ivj,iset)-YO(i))/YU(i)
                      ENDDO
                ENDIF
              IF(BOLTZ(iset).GT.0) RMS(iset)= RMS(iset)/NFREQ(iset)
              IF(BOLTZ(iset).LE.0) RMS(iset)= RMS(iset)/NVJ(iset)
              RMS(iset)= DSQRT(RMS(iset))
c             WRITE(8,806) NFREQ(iset),iset,RMS(iset)
              WRITE(6,806) NFREQ(iset),iset,RMS(iset)
              ENDDO
          ENDIF
c** If desired ... print potential & TMF arrays compactly to channel-9
      DO  ifs= 1, NFS
          IF (LPFS.GT.0) THEN
              WRITE(9,900) ifs,LPFS
              WRITE(9,902) (RAD(i),VF(i,ifs),i=1,NIN,LPFS)
              ENDIF
          IF (LPTMF.GT.0) THEN
              WRITE(9,904) ifs,LPTMF
              DO i= 1,NIN,LPTMF
                  ztmf(i,5)= 0.d0
                  xtmf= 1.d0
                  DO m= 0,OTMF(ifs)
                      zTMF(i,5)= zTMF(i,5) + TMFPRM(m,ifs)*xtmf
                      xtmf= xtmf*ztmf(i,ifs)
                      ENDDO
                  ENDDO
              WRITE(9,902) (RAD(i),zTMF(i,5), i= 1,NIN,LPTMF)
              ENDIF
          ENDDO
      STOP
c-----------------------------------------------------------------------
  600 FORMAT(/1x,a75/1x,25('===')/' Consider',I3,' isotopomer(s)'/1x,
     1 69('-')/6x,'Isotopomer',7x,'Mass of atom-1   Mass of atom-2    Re
     2duced mass'/2x,'----------------- ',3('   --------------'))
  601 FORMAT(/" Potential is type-2, so PREPOT's inward extrapolation is
     1 replaced by"/10x,'a variable inner exponential wall ',
     2 'of the form:'/'   V(r)= A + B*exp{-(r - REXFS)*[A1 + A2*y + ',
     3 'A3*y^2 + ...]};'/16x,'A =',f14.6/16x,'B =',F14.6/
     4 (15x,'A',i1,' =',F14.8))
  602 FORMAT(2x,A2,'(',i3,') - ',A2,'(',I3,')',3(3x,F14.9))
  603 FORMAT(/' *** ERROR *** Read-in  NVJ(iset=',i2,')=',i3,
     1  '  when  FITIT=',i2)
  604 FORMAT(/' **BCONT ERROR** ',I3,' is an invalid value for IFRPW.')
  606 FORMAT(2x,69('-')//' Perform a Fit to',i3,' Set(s) of Experimental
     1 Data.'/1x,48('='))
  607 FORMAT(/' Apply "Sequential Rounding & Refitting" at digit-',
     1  i1,' of the (local) parameter')
  608 FORMAT(3x,'uncertainty, selecting remaining parameter with largest
     1 relative uncertainty.')
  609 FORMAT(3x,'uncertainty, proceeding sequentially from the LAST para
     1meter to the FIRST.')
  610 FORMAT(2x,69('-')//' Perform Forward Calculations for',i3,' Cases'
     1  /1x,52('='))
  611 FORMAT(/' Value of nj above upper limit, mxnj =',i3,'. Instead ',
     1 'of dividing the population'/1x,'into',i3,' equal segments, do ',
     2 'direct sum from J = 0 to J =',i3,'.')     
  612 FORMAT(' * Assume Boltzman population of initial states.')
  613 FORMAT(' * Sum over P, Q and R branches for final-state rotational
     1 levels')
  615 FORMAT(' * Use Q-branch approximation for final state rotational l
     1evels')
  616 FORMAT(/' ** DIMENSIONING ERROR ** number of input frequencies ',
     1 'cannot exceed ',i4)
  617 FORMAT(/' ** DIMENSIONING ERROR ** number of input data sets ',
     1 'cannot exceed ',i3)
  618 FORMAT(' * Directly sum over thermally weighted rotational levels'
     1 ,' from   J = 0  to',I4)
  619 FORMAT(/' ** DIMENSIONING ERROR ** number of input isotopes ',
     1 'cannot exceed ',i2)
  620 FORMAT(/' Fix J = 0 rather than summing over rotational ',
     1 'distribution.')
  621 FORMAT(/' ** DIMENSIONING ERROR ** number of final states ',
     1 'cannot exceed ',i2)
  622 FORMAT(/' Divide rotational distbn for each vib. level into',I3,
     1 ' equally weighted segments.')
  623 FORMAT(' Transitions from initial-state ',A2,'(',i3,')-',A2,'(',
     1  i3,') level   v=',i3,'    J=',i4)
  624 FORMAT(' * Initial-state thermal sum includes',F8.5,' of the vibra
     1tional population')
  625 FORMAT(' * Thermal  T=',F7.2,'K  initial state sum for ',A2,'(',
     1  i3,') - ',A2,'(',i3,')  has   VMAX=',i3)
  626 FORMAT( ' Calculations use exact (numerical) final-state continuum
     1 wave functions'/' Amplitude convergence assumed when constant to'
     2  ,1PD8.1,' at 3 successive maxima')
  628 FORMAT(2x,69('-')/' Perform calculations using delta function appr
     1oximation for final-state'/'   continuum wave functions (reflectio
     2n method).'//' *** NB. ONLY EXACT QUANTAL CALCULATION AVAILABLE ',
     3 'AT THIS TIME ***')
  629 FORMAT(' Property is sum of intensities over all',i2,
     1  ' final states with weight factors:'/(5x,10I4))
  630 FORMAT(' Case',i2,' property is ratio of final-state intensities w
     1ith'/8x,'Numerator weight coefficients:',10I4:/(38x,10I4))
  631 FORMAT(6x,'& Denomator weight coefficients:',10I4)
  632 FORMAT(' * Units for calculated Molar Absorption Coefficients: ',
     1 '(L/mol*cm)/(cm-1)')
  633 FORMAT(/' *** ERROR *** For case',i2,' read-in  DTYPE=',i2,' INVAL
     1ID  (Not equal 1 or 2)!')
  634 FORMAT(' * Units of calculated Einstein Emission Intensity ',
     1 'Coefficients: (sec-1)/(cm-1)')
  637 FORMAT(/' Set-',i2,': ',A70/1x,35('=='))
  638 FORMAT(/i4,' Set #',i2,'  Experimental Intensity Coefficients:'/
     1 3x,69('-')/2('   FREQ/cm-1   WAVL/nm    OBS   UNC ')/
     2 2(3x,11('---'))/(2(f12.2,f10.4,f8.2,f6.2)))
  640 FORMAT(/' Predict',i4,'  Intensities at Frequencies (cm-1):'/
     1  (8f10.2:))
  643 FORMAT(/i4,' Experimental Predisociation Rates (s-1) in Data Set #
     1',i2/1x,36('--')/2('   v   J  RATE(sec-1)     UNC   ',8x)/
     2                       2('   -   -  -----------  ---------',8x)/
     3 (2(i4,i4,1PD13.4,D11.3,8x)))
  644 FORMAT(/' **NOTE** For modelling purposes, replace  (frequency)^',
     1 I1,'  factor'/' with the constant:  (',F8.2,')^',I1)
  645 FORMAT(/' Consider ',i3,' specific (v,J) Levels for Set #',i2,';',
     1 2x,A2,'(',i3,') - ',A2,'(',I3,')'/
     2 1x,39('--')/(8(2x:'(',i2,',',i2,')')):)
  647 FORMAT(/1x,i3,' Experimental Predisociation Lifetimes for Set ',
     1 i2/1x,36('--')/2('   v   J  LIFETIME(s)     UNC   ',8x)/
     2            2('   -   -  -----------  ---------',8x)/
     3 (2(i4,i4,1PD13.4,D11.3,8x)))
  648 FORMAT(/1x,i3,' Experimental Predisociation Linewidths (FWHM) for'
     1 ,' Set ',i2/1x,36('--')/2('   v   J  WIDTH(cm-1)     UNC   ',8x)/
     2                       2('   -   -  -----------  ---------',8x)/
     3 (2(i4,i4,1PD13.4,D11.3,8x)))
  649 FORMAT(/' Perform calculations for ALL (v,J) levels ',
     1 'between (',I3,',',I3,') and (',I3,',',I3,')')
  652 FORMAT(' Integrate initial-state wavefunctions from  RMIN =',f6.3,
     1 ' to  RMAX =',f7.3/5x,'with mesh  RH =',f10.7,' [Angstroms].')
  654 FORMAT(/' Calculations use electronic state degeneracy factors:',
     1    8i4)
  656 FORMAT(/' For the initial electronic state:'/1x,17('=='))
  662 FORMAT(/' Data for ',A2,'(',i3,')-',A2,'(',I3,
     1 ') uses initial-state levels up to  v = ',I2/1x,12('=='),14x,
     2 'whose molecular constants are:'/2x,'v',7x,'Gv',
     3  9x,'Bv',9x,'-Dv',11x,'Hv',10x,'Lv',11x,'Mv'/1x,39('--')/
     4  (I3,0PF13.5,F11.7,1PD13.5,D13.5,D13.5,D13.5))
  664 FORMAT(/1x,25('==')/' Consider transitions to',i3,' separate final
     1 state(s)'/1x,25('=='))
  666 FORMAT(/' **CAUTION** Value for NFSPRM was too large, ',
     1 'internally reduced to',I3)
  667 FORMAT(/' Final State',i2,'   with   OMEGA=',i2/1x,13('='))
  668 FORMAT(' Radial expansion variable is   y = (r - REXFS)/(r + REXFS
     1)  &  REXFS=',F9.6)
  669 FORMAT(' Radial expansion variable is   y = (r^',i1,' - REXFS^',
     1 i1,')/(r^',i1,' + REXFS^',i1,')'/11x,'where   REXFS=',F9.6)
  670 FORMAT(' Radial expansion variable is   y = (r - REXFS)/r',
     1 '  &   REXFS=',F9.6)
  671 FORMAT(' Radial expansion variable is   y = (r - REXFS)/REXFS',
     1 '  &   REXFS=',F9.6)
  672 FORMAT(/' *** Value for  XCOORD(ifs)=',i3,'  INVALID **** (must be
     1 .ge.1  and  .le.11).')
  674 FORMAT(' Potential is type-1,  an exponential of the form:'/7x,
     1 'V(r) = VLIMF + A1*exp[-(r - REXFS)*(A2 + A3*y + A4*y^2 + ...)]'/
     2 6x,'VLIMF =',f17.8/(9x,'A',i1,' =',F17.8))
  675 FORMAT(' Potential is type-3,  an Extended Morse Oscillator functi
     1on:'/7x,'V(r) = VLIMF + A1*{exp[-(r - A2)*(A3 + A4*y + A5*y^2 + ..
     2.)] - 1}**2 - A1'/6x,'VLIMF =',f17.8/(9x,'A',i1,' =',F17.8))
  676 FORMAT('- Absolute energy at potential asymptote:',T45,'VLIMF =',
     1 F12.5,' cm-1'/'- Perform',I3,'-point piecewise polynomial interpo
     2lation over',I5,' input points')
  677 FORMAT('- Absolute energy at potential asymptote:',T45,'VLIMF =',
     1 F12.5,' cm-1'/'- Perform cubic spline interpolation over the',
     2 I5,' input points')
  678 FORMAT('- Interpolation actually performed over modified input arr
     1ay:   V(i) * R(i)**2')
  679 FORMAT('- Beyond read-in points extrapolate to limiting asymptotic
     1 behaviour:'/20x,'V(R)  =  V(lim) - (',D16.7,')/R**',I2)
  680 FORMAT('- To make input points V(i) consistent with  V(lim),  add'
     1 ,'  V(shift)=',F12.4/'- Scale input points:  (distance)*',
     2 1PD16.9,'  &  (energy)*',D16.9/13x,'to get required internal unit
     3s  [Angstroms & cm-1 for potentials]'/
     4 3('      R(i)         Y(i)  ')/3(3X,11('--'))/(3(0PF13.8,F12.4)))
  681 FORMAT(1x,38('--'))
  682 FORMAT(/' *** Fatal Error in TMF ***  OTMF =',i3,' unallowed - ',
     1 'must be less than ',I2)
  686 FORMAT(' Coefficients of power series expansion for transition mom
     1ent function are:'/(3x,6(f12.8):))
  687 FORMAT(/' Parameter Summary'/1x,17('=')/' i ifs     FSPRM(i,ifs)',
     1 '   FSVAR(i,ifs)'/' - --- -----------------  ------------')
  688 FORMAT(1x,i1,2x,i1,1x,f18.10,9x,i1)
  689 FORMAT(/' m ifs     TMFPRM(m,ifs)  TMFVAR(m,ifs)'/
     1 ' - --- -----------------  -------------')
  690 FORMAT(1x,i1,2x,i1,1x,f18.10,9x,i1)
  691 FORMAT(/' iset     SCALING FACTOR      SCALE '/
     1 ' ----  -----------------  -------------')
  692 FORMAT(1x,i2,13x,f8.6,9x,i1)
  693 FORMAT(/14x,'Correlation Matrix'/'  1',f7.3,4x,9('--'))
  694 FORMAT(i3,20(f7.3))
  695 FORMAT(' Use Lorentzian Function for TMF; W= a/{4a^2 + (R-Rc)^2}'
     1  /6x,'a= ',f12.8,/5x,'Rc= ',f12.8)
  698 FORMAT(/' *WARNING* Variable Final State Param. cannot be zero, so
     1 set  FSPRM(',i1,',',i1,')= 0.001')
  699 FORMAT(1x,39('--'))
  802 FORMAT(/' Summary for data set #',i1/1x,a70/56x,'CALC-OBS'/
     1  '  FREQUENCY',5x,'Y(obs)    UNC     Y(calc)    CALC-OBS',
     2  5x,'/UNC'/1x,31('--'))
  803 FORMAT(2(f11.3),f8.3,3(f11.4))
  804 FORMAT(/' Summary for data set #',i1/1x,a70/71x,'CALC-OBS'/
     1  2x,' v   J   ',5x,'Y(obs)    UNC     Y(calc)    CALC-OBS',
     2  5x,'/UNC'/1x,31('--'))
  805 FORMAT(2(i4),5(1PD13.5),f13.4)
  806 FORMAT(' RMS deviation for the',i4,' data of set-',i2,' = ',f8.3)
  900 FORMAT(/1x,'Final State',I2,'  potential array written at every ',
     1 i3,'th mesh point.')
  902 FORMAT((4(F8.3,F12.4)):)
  904 FORMAT(/' Transition Moment Function for Final State',i2,
     1 ' written every ',i3,'th mesh point.')
c-----------------------------------------------------------------------
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE FORWARD
c=======================================================================
c--------------------  Last updated 9 February 2004 --------------------
c=======================================================================
ccc   INCLUDE 'arrsizes.h'
c-----------------------------------------------------------------------
c  Utility routine to summarize dimensioning of arrays
c-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,
     1        mxsets,mxfreq
      REAL*8 CCM,PI
c-----------------------------------------------------------------------
c  mxdata - maximum number of input data points
c  mxisp  - maximum number of points for initial state potential array
c           (also used for number of points in transition moment array)
c  mxnj    - maxiumum value of j quantum number allowed
c  mxfsp  - maximum number of points for final state potential array
c  mxnp   - maximum number of parameters total
c  mxntp  - maximum number of turning points to be read in
c  mxprm  - maximum number of parameters for final state pot'l or TMF
c  mxv    - largest value for the v quantum number
c  mxfs   - maximum number of final states allowed
c  mxisot - maximum number of isotopomers allowed
c  mxsets - maximum number of data sets allowed
c  mxfreq - maximum number of data points allowed in a given set
c-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
c=======================================================================
      INTEGER c,FITIT,i,ifr,ifs,INNER(0:mxv,mxisot),iset,IV,ivj,IWRSCH,
     1 IWROVR,j,JDP,JFRPW,JNNER,JMIN,jp,jpmax,jpmin,JPWR,KV,LPRWF,LPTMF,
     2 m,MCALC,MISS,NBEG,NEND,NFS,NIN,NSETS,printyn,RPD,v,VMIN,VLAST,
     3  BOLTZ(mxsets),CD(mxsets,mxfs),CN(mxsets,mxfs),DTYPE(mxsets),
     4  GFS(mxfs),IFRPW(mxsets),ISOT(mxsets),J1ST(mxsets),JFIX(mxfreq),
     5  JM(mxnj),JMAX(mxsets),NFREQ(mxsets),NJ(mxsets),NVJ(mxsets),
     6  OMEGA(0:mxfs),OTMF(mxfs),PQR(mxsets),SCALE(mxsets),TMFTYP(mxfs),
     7  TMFVAR(0:mxprm-1,mxfs),V1ST(mxsets),VFIX(mxfreq),VMAX(mxsets)
c=======================================================================
      REAL*8 CHNG,EVIN,ESAV,FACT,FACT1,FRQFCT,FTST,GAMA,HLFACT,JSCALE,
     1  EDIF,EFN,EJTST,EPS,OVR,OVRCRT,QTOT,QTOTI,QTST,RDF,RH,RMIN,SUM,
     2  TCM,VLIMI,ZKLIM,
     3  BFCT(mxisot),DER(0:mxprm-1),dIdT(0:mxprm-1,mxfreq,mxfs,mxsets),
     4  ESEGM(mxnj),FACTOR(mxsets),FREQ(mxfreq,mxsets),
     5  MCI(0:mxv,0:7,mxisot),OBS(mxfreq,mxsets),UNC(mxfreq,mxsets),
     6  OUTPUT(mxfreq,mxsets),
     6  POPCRT,POPF,PSI(mxfsp),QVB(0:mxv),RAD(mxfsp),SF(mxsets),
     7  TMFPRM(0:mxprm-1,mxfs),TEMP(mxsets),TOTFS(mxfreq,mxfs,mxsets),
     8  UMAX,UJ(mxfsp),VF(mxfsp,mxfs),VI(mxisp,mxisot),
     9  VICD(mxisp,mxisot),VLIMF(mxfs),WAVL(mxfreq,mxsets),XM2(mxfsp),
     a  ztmf(mxisp,mxfs),zin(mxisp)
      CHARACTER*70 INFO(mxsets)
c=======================================================================
      COMMON /MF/ BFCT,EPS,FACTOR,FREQ,MCI,OBS,UNC,OVRCRT,POPCRT,TEMP,
     1   VI,VICD,VLIMI,WAVL,BOLTZ,CD,CN,DTYPE,FITIT,INNER,ISOT,IWRSCH,
     2   IWROVR,J1ST,JFIX,JFRPW,JMAX,MCALC,NJ,OMEGA,PQR,printyn,
     3   V1ST,VFIX,VMAX
      COMMON /MFD/ dIdT,OUTPUT,SF,IFRPW,NFREQ,NFS,NSETS,NVJ,RPD,SCALE,
     1             TMFVAR,INFO
      COMMON /MFGf/ VF,VLIMF
      COMMON /MFGt/ XM2,ztmf,LPTMF,NIN,TMFTYP
      COMMON /MFGtGf/ RAD
      COMMON /MFDGt/ TMFPRM,OTMF,GFS
c=======================================================================
      DATA LPRWF/0/
c=======================================================================
c  LPRWF specifies print level for initial state wavefuntion (chan.8)
c    > 0  print wavefunction every LPRWF-th point.
c    < 0  print every |LPRWF|-th point of wave function starting at
c         R(NBEG) with step size |LPRWF|*RH
c    = 0  to avoid printing
c-----------------------------------------------------------------------
      RMIN= RAD(1)
      RH= RAD(2)-RAD(1)
      ivj= 0
c-----------------------------------------------------------------------
c begin by zeroing state-intensity and partial derivative arrays
c-----------------------------------------------------------------------
      DO iset = 1,NSETS
          DO ifr = 1,mxfreq
              DO ifs = 1,NFS
                  TOTFS(ifr,ifs,iset)= 0.d0
                  IF(RPD.GT.0) THEN
                      DO m=0,OTMF(ifs)
                          dIdT(m,ifr,ifs,iset)= 0.d0
                          ENDDO
                      ENDIF
                  ENDDO
              ENDDO
          ENDDO
c-----------------------------------------------------------------------
c  Begin actual calculation ... consider NSETS of data
c-----------------------------------------------------------------------
      IF(FITIT.GT.0) WRITE(12,1200)
 1200 FORMAT(' VARIABLES = FREQ,OBS,unc,calc,int1,int2,int3,int4')
 1201 FORMAT(/'zone t="',a70,'"')
      DO 1000 iset = 1,NSETS
          IF(printyn.EQ.1) THEN
              WRITE(6,697)
              WRITE(6,690) INFO(iset)
              ENDIF
          IF((IFRPW(iset).EQ.0).AND.(printyn.EQ.1)) THEN
c ... for predissociation calculations
              WRITE(6,600) (ifs,ifs= 1,NFS)
              WRITE(6,699)
              IF(FITIT.GT.0) THEN
                  WRITE(11,1105) (ifs,ifs= 1,NFS)
                  WRITE(11,699)
                  WRITE(12,1201) INFO(iset)
                  ENDIF
              ENDIF
          IF(BOLTZ(iset).GT.0) THEN
c ... for thermal absorption/emission calculations
              VMIN = 0
              JMIN = 0
              TCM= 0.6950387D0*TEMP(iset)
            ELSE
c ... for a transition from a specified state  V1ST,J1ST
              TCM= 9.d9
              IF(NVJ(iset).GT.0) THEN
c ... if that transition is predissociation
                  VMIN= 1
                  VMAX(iset)= 1
                  JMIN= 1
                  JMAX(iset)= NVJ(iset)
                ELSE
c ... if that transition is absorption/emission
                  VMIN= V1ST(iset)
                  VMAX(iset)= V1ST(iset)
                  JMIN= J1ST(iset)
                  JMAX(iset)= J1ST(iset)
                ENDIF
            ENDIF   
          QTOT= 0.d0
c-----------------------------------------------------------------------
c  Specify rotational energy at which to cut off (direct) thermal J-sum
c-----------------------------------------------------------------------
          IF((NJ(iset).LT.0).AND.(BOLTZ(iset).GT.0))
     1                                   EJTST= -TCM*DLOG(1.d0-POPCRT)
c-----------------------------------------------------------------------
c  Begin loop over vibrational levels
c-----------------------------------------------------------------------
          DO 250 v = VMIN,VMAX(iset)
              QVB(v)= 0.d0
              IF((NJ(iset).GE.0).AND.(BOLTZ(iset).GT.0))
     1               CALL JAVGE(v,NJ(iset),JM,MCI(v,1,ISOT(iset)),TCM)
c-----------------------------------------------------------------------
c** JAVGE uses equally weighted segments of the rotational population
c         intead of doing complete sum over all J" values.
c  Inputs:  v,NJ,{MCI(v,1,isot)= Bv},TCM
c  Output:  JM(NJ)- the average J for each of the nj equally weighted
c                  segments of the rotational population for vibrational
c                  level whose rotational constant is BV
c-----------------------------------------------------------------------
c  For predissociation,  calculations, the counters below are used to
c  represent the jth predissociating state. j runs from 1 to NVJ(iset)
c-----------------------------------------------------------------------
              DO 200 j = JMIN,JMAX(iset)
c-----------------------------------------------------------------------
c  j is the counter for loop over lower state J" (=JDP)
c  Loop to incorporate population scaling factor for absolute result
c-----------------------------------------------------------------------
                  ivj= ivj + 1
                  IF((BOLTZ(iset).GT.0).OR.(NVJ(iset).LE.0)) THEN
                      IV= v
                      IF(NJ(iset).LE.0) JDP= j
                      IF(NJ(iset).GT.0) JDP= JM(j+1)
                    ELSE
                      IV= VFIX(j)
                      JDP= JFIX(j)
                    ENDIF
                  EVIN= 0.d0
                  JPWR= 1
                  DO c= 0,7
                      EVIN= EVIN + MCI(IV,c,ISOT(iset))*JPWR
                      JPWR= JPWR*JDP*(JDP+1)
                      ENDDO
                  JSCALE= DBLE(RH*RH*(JDP*(JDP+1)-OMEGA(0)**2))
c-----------------------------------------------------------------------
c  now create centrifugally-distorted potential UJ(i) for initial state
c-----------------------------------------------------------------------
                  DO  i= 1,NIN
                      UJ(i)= VI(i,ISOT(iset))+ JSCALE*VICD(i,ISOT(iset))
                      ENDDO
                  MISS= 0
c-----------------------------------------------------------------------
c  UJ(i) above is used here for the initial state potential array
c  WARNING - definition changes later and is used for final state
c-----------------------------------------------------------------------
   50             KV= IV
c-----------------------------------------------------------------------
c  SCHRQ Summary
c  inputs  KV,JDP,EVIN,VLIMI,UJ,BFCT,EPS,RMIN,RH,NIN,INNER,IWRSCH,LPRWF
c  outputs  KV   - value of v may change inside program
c           EVIN - which is now more precisely calculated
c           GAMA -
c           UMAX - maximum value of potential barrier (real)
c           PSI  - wavefunction array
c           NBEG - beginning mesh point for wavefunction array
c           NEND - final mesh point for wavefunction array
c-----------------------------------------------------------------------
c  INNER specifies wave func'n matching (& initiation) condition (schrq)
c    = 0  match inward & outward solutions at outermost wave function
c         maximum; otherwise match at inner edge of classically allowed
c         region
c    < 0  uses zero slope inner boundary condition.
c   
c   For most cases set INNER = 0 but to find "inner-well-dominated"
c    solutions of an asymmetric double minimum potential, set INNER > 0
c   To find symmetric eigenfunctions of a symmetric potential, set
c    INNER < 0  & start integration (set RMIN) at potential mid point
c-----------------------------------------------------------------------
                  JNNER= INNER(KV,ISOT(iset))
                  ESAV= EVIN
                  CALL SCHRQ(KV,JDP,EVIN,GAMA,UMAX,VLIMI,UJ,PSI,
     1                    BFCT(ISOT(iset)),EPS,RMIN,RH,NIN,NBEG,NEND,
     2                    JNNER,IWRSCH,LPRWF)
                  IF(KV.NE.IV) THEN
c** If find the WRONG level, try to fix it up!
                      MISS= MISS+1
                      IF(MISS.EQ.1) THEN
c ... first, assuming a double-well problem, switch the value of INNER
                          IF(INNER(KV,ISOT(iset)).GT.0) JNNER= 0
                          IF(INNER(KV,ISOT(iset)).LE.0) JNNER= 1
                        ELSE
c ... otherwise ... just jiggle up or down with original INNER
                          JNNER= INNER(KV,ISOT(iset))
                          IF(KV.GT.IV) THEN
                              CHNG= (MCI(KV,0,ISOT(iset))-
     1                        MCI(IV,0,ISOT(iset)))*(KV-IV)/ABS(KV-IV+1)
                              WRITE(6,602) IV,JDP,ESAV,KV,EVIN,ESAV-CHNG
                              EVIN= ESAV - CHNG
                            ELSE
                              CHNG= (MCI(IV,0,ISOT(iset))-
     1                        MCI(KV,0,ISOT(iset)))*(IV-KV)/ABS(IV-KV+1)
                              WRITE(6,602) IV,JDP,ESAV,KV,EVIN,ESAV+CHNG
                              EVIN= ESAV + CHNG
                            ENDIF
                        ENDIF
                      IF(MISS.GT.4) STOP
                      GOTO 50
                      ENDIF
                  IF((NJ(iset).GE.0).AND.(BOLTZ(iset).GT.0))
     1                                                ESEGM(j+1)= EVIN
c-----------------------------------------------------------------------
c ESEGM(J") is the energy of the Jth equally weighted rotational segment
c-----------------------------------------------------------------------
                  IF(BOLTZ(iset).GT.0) THEN
                      IF(NJ(iset).LT.0) THEN
                          POPF=(2.d0*JDP + 1)*
     1                           DEXP(-(EVIN-MCI(0,0,ISOT(iset)))/TCM)
                        ELSE
                          POPF= DEXP(-(MCI(v,0,ISOT(iset))-
     1                                  MCI(0,0,ISOT(iset)))/TCM)*TCM/
     2                         (MAX(NJ(iset),1)*MCI(v,1,ISOT(iset)))
                        ENDIF
                    ELSEIF(BOLTZ(iset).LE.0) THEN
                      POPF= 1.d0
                    ENDIF
                  QVB(v)= QVB(v) + POPF
                  QTOT= QTOT + POPF
                  FACT1= FACTOR(iset)*POPF
c-----------------------------------------------------------------------
c  Now loop over contributions from excited electronic final states
c-----------------------------------------------------------------------
                  DO 150 ifs = 1,NFS
c-----------------------------------------------------------------------
c  Change TMF array dimensions from 2 to 1 (for ovrlap subroutine)
c-----------------------------------------------------------------------
                      DO i= 1,NIN
                          zin(i)= ztmf(i,ifs)
                          ENDDO
c-----------------------------------------------------------------------
c  For model predissociation calculations, calculate bound-state
c  expectation value of transition moment function.
c-----------------------------------------------------------------------
cgtk - please check if this entire IF loop is even necessary
cgtk               IF(IFRPW(iset).EQ.0) THEN
cgtk                   SUM= 0.d0
cgtk                   DO i= NBEG,NEND
cgtk                       RDF= 1.D0
cgtk                       DO m= 0,OTMF(ifs)
cgtk                           SUM= SUM + (TMFPRM(m,ifs)*RDF)*PSI(i)**2
cgtk                           RDF= RDF*ztmf(i,ifs)
cgtk                           ENDDO
cgtk                       ENDDO
cgtk                   SUM= SUM*RH
cgtk                   IF(printyn.EQ.1) WRITE(6,604) KV,JDP,EVIN,SUM
cgtk                   ENDIF
c-----------------------------------------------------------------------
c  If PQR = 0, program collapses sum over final-state rotational
c    quantum numbers to the single (Q-branch) term
c  If PQR = 1, program allows for P,Q, and R branches with appropriate
c    Honl-London factors
c-----------------------------------------------------------------------
                      IF(PQR(iset).LE.0) THEN
                          jpmin= JDP
                          jpmax= JDP
                        ELSE
                          jpmin= JDP-1
                          jpmax= JDP+1
                        ENDIF
                      DO 125 jp= jpmin,jpmax
c-----------------------------------------------------------------------
c  jp = J', the upper state value of J
c-----------------------------------------------------------------------
                          CALL HONL(HLFACT,JDP,jp,OMEGA,PQR(iset),ifs)
c-----------------------------------------------------------------------
c  HONL inputs:  JDP,jp,OMEGA(i=0,ifs),PQR(iset),ifs
c       output:  HLFACT - the Honl-London factor, divided by (2J+1), for
c                       a specified branch for a given delta-omega value
c-----------------------------------------------------------------------
                          IF(HLFACT.LE.0.d0) GOTO 125
c** Include Electronic degeneracy factor ...
                          HLFACT= HLFACT*GFS(ifs)
c-----------------------------------------------------------------------
c  Add centrifugal term to final-state potential before overlap calcn.
c  Note that here UJ(i) is used for the final state potential array
c  since the initial potential is no longer needed - wavefunction array
c  already calculated and stored.
c  Note that final-state potential is in (scaled) internal units
c-----------------------------------------------------------------------
                          JSCALE= DBLE(RH*RH*jp*(jp+1))
                          DO i= 1,mxfsp
                              UJ(i)= VF(i,ifs)*BFCT(ISOT(iset)) +
     1                                                   JSCALE*XM2(i)
                              ENDDO
c-----------------------------------------------------------------------
                          DO 100 ifr= 1,NFREQ(iset)
                              IF(JFRPW.EQ.3) THEN
                                  EFN= EVIN - FREQ(ifr,iset)
                                ELSE
                                  EFN= EVIN + FREQ(ifr,iset)
                                ENDIF
                              IF((EFN.LE.VLIMF(ifs)).AND.
     1                                        (IFRPW(iset).EQ.0)) THEN
                                  IF(printyn.EQ.1)
     1                                    WRITE(6,606) IV,JDP,EFN,0.d0
                                  OUTPUT(ivj,iset)= 0.d0
                                  TOTFS(ivj,ifs,iset)= 0.d0
                                  GOTO 200
                                  ENDIF
                              IF(EFN.LE.VLIMF(ifs)) GOTO 100
c-----------------------------------------------------------------------
c  EFN must be above the lower asymptote for upper potential
c-----------------------------------------------------------------------
                              FRQFCT= 1.d0
                              IF(IFRPW(iset).GT.0)
     1                             FRQFCT= FREQ(ifr,iset)**IFRPW(iset)
                              IF(IFRPW(iset).LT.0)
     1                                     FRQFCT= FREQ(1,iset)**JFRPW
                              ZKLIM= DSQRT(EFN - VLIMF(ifs))
                              FACT= FACT1*HLFACT*FRQFCT/ZKLIM
c-----------------------------------------------------------------------
c  Perform overlap integral calculation in OVRLAP or OVRPD
c  **OVRPD is the special case for PreDissociation with TMF as operator
cgtk OVRDLT not available at this time
c-----------------------------------------------------------------------
cgtk                          IF (MCALC.EQ.0) CALL OVRDLT()
c-----------------------------------------------------------------------
                              IF (MCALC.GT.0) THEN
                                  IF(TMFTYP(ifs).GE.0) THEN
      CALL OVRLAP(BFCT(ISOT(iset)),DER,EFN,OVR,OVRCRT,PSI,RH,RMIN,
     1   TMFPRM,UJ,VLIMF(ifs),zin,ifs,IWROVR,jp,NEND,OTMF(ifs),TMFTYP)
                                    ELSE
      CALL OVRPD(BFCT(ISOT(iset)),DER,EFN,OVR,OVRCRT,RAD,PSI,TMFPRM,UJ,
     1     VLIMF(ifs),FITIT,ifs,IWROVR,jp,LPTMF,NEND,OTMF(ifs),TMFVAR)
                                    ENDIF
                                  ENDIF
                              OVR= FACT*OVR
                              IF(RPD.GT.0) THEN
                                  DO m= 0,OTMF(ifs)
                                      IF(IFRPW(iset).EQ.0) THEN
                                          dIdT(m,ivj,ifs,iset)=
     1                                            FACT*DER(m)*SF(iset)
c  note that for predissociation, there is no accumulation of TOTFS
c  since dealing with only one frequency... easiest to add scaling
c  factor here for predissociation - scale abs/emmission cases later
                                        ELSE
                                          dIdT(m,ifr,ifs,iset)=
     1                               dIdT(m,ifr,ifs,iset)+ FACT*DER(m)
                                        ENDIF
                                      ENDDO
                                  ENDIF
                              IF(IFRPW(iset).EQ.0) THEN
c???
c??? Geoff Had SF scaling here????  TOTFS(ivj,ifs,iset)= OVR*SF(iset)
c???
                                  TOTFS(ivj,ifs,iset)= OVR
                                ELSE
                                  TOTFS(ifr,ifs,iset)=
     1                                         TOTFS(ifr,ifs,iset)+OVR
                                ENDIF
  100                         CONTINUE
c ... end of loop over frequencies
c-----------------------------------------------------------------------
  125                     CONTINUE
c ... end of loop over J'
c-----------------------------------------------------------------------
  150                 CONTINUE
c ... end of loop over final electronic states
c-----------------------------------------------------------------------
                  IF(IFRPW(iset).EQ.0) THEN
                      CALL ADD(iset,ivj,NFS,DTYPE(iset),CN,CD,TOTFS,
     1                                           RPD,OTMF,dIdT,OUTPUT)
c  subroutine "ADD" simply adds the contribution of each final state
                      IF(printyn.EQ.1) THEN
                          WRITE(6,606) KV,JDP,EFN,OUTPUT(ivj,iset),
     1                                 (TOTFS(ivj,ifs,iset),ifs=1,NFS)
                          IF(FITIT.GT.0) THEN
                              WRITE(11,1106) KV,JDP,EFN,OBS(ivj,iset),
     1                 OUTPUT(ivj,iset),(TOTFS(ivj,ifs,iset),ifs=1,NFS)
                              WRITE(12,1106) KV,JDP,EFN,OBS(ivj,iset),
     1                 OUTPUT(ivj,iset),(TOTFS(ivj,ifs,iset),ifs=1,NFS)
                              ENDIF
                          ENDIF
                      ENDIF
                  IF((NJ(iset).LT.0).AND.(BOLTZ(iset).GT.0)) THEN
c-----------------------------------------------------------------------
c  In appropriate case, truncate direct sum over initial-state
c  rotational population when fraction POPCRT of the rotational
c  population for that vibrational level is accounted for.
c-----------------------------------------------------------------------
                      IF((EVIN-MCI(v,0,ISOT(iset))).GT.EJTST) THEN
                          IF(printyn.EQ.1)
     1                             WRITE(6,608) KV,TEMP(iset),JDP,EVIN
                          GOTO 210
                          ENDIF
                    ELSEIF((BOLTZ(iset).LE.0).AND.(IFRPW(iset).NE.0))
     1                                                            THEN
                      WRITE(6,609) V1ST(iset),J1ST(iset)
                      ENDIF
  200             CONTINUE
c ... end of loop over initial-state J"
c-----------------------------------------------------------------------
  210         IF((IFRPW(iset).NE.0).AND.(printyn.EQ.1).AND.
     1  (NJ(iset).GT.0)) WRITE(6,610) KV,(JM(i),ESEGM(i),i=1,NJ(iset))
c-----------------------------------------------------------------------
c** As appropriate, truncate thermal vibrational sum when fraction
c  POPCRT of possible levels have been accounted for.  (Missing levels'
c  effect estimated using a harmonic approximation).
c-----------------------------------------------------------------------
              IF((v.GT.0).AND.(BOLTZ(iset).GT.0)) THEN
                  IF (v.LT.VMAX(iset)) THEN
                      EDIF= MCI(v+1,0,ISOT(iset))-MCI(v,0,ISOT(iset))
                    ELSE
                      EDIF= MCI(v,0,ISOT(iset))-MCI(v-1,0,ISOT(iset))
                    ENDIF
                  FTST= DEXP(-(EDIF+MCI(v,0,ISOT(iset))-
     1              MCI(0,0,ISOT(iset)))/TCM)*TCM/MCI(v,1,ISOT(iset))/
     2                                          (1.d0-DEXP(-EDIF/TCM))
                  QTST= QTOT/(QTOT+FTST)
                  VLAST= v
                  IF(QTST.GT.POPCRT) GOTO 252
              ENDIF
  250     CONTINUE
c-----------------------------------------------------------------------
c ... 250 ends loop over initial-state v"
c  at this point multiply by scaling factor and boltzman factor for
c  absorption and emmission cases
c-----------------------------------------------------------------------
cc  252     IF(BOLTZ(iset).GT.0) THEN   
  252     CONTINUE
          IF(BOLTZ(iset).GT.0) THEN
              DO i= 0,VLAST
                  QVB(i)= QVB(i)/(QTOT+FTST)
                  ENDDO
              IF(printyn.EQ.1) WRITE(6,620) TEMP(iset),VLAST+1,QTST,
     1                                           (i,QVB(i),i= 0,VLAST)
cc              ENDIF
              QTOTI= 1.d0/QTOT
              DO ifs = 1,NFS
                  DO ifr= 1,NFREQ(iset)
                      TOTFS(ifr,ifs,iset)= QTOTI*TOTFS(ifr,ifs,iset)
     1                                                       *SF(iset)
c???
c??? Omit SF here ... property of OBS not of CALC !!
c???                  TOTFS(ifr,ifs,iset)= QTOTI*TOTFS(ifr,ifs,iset)
c???
                      IF(RPD.gt.0) THEN
                          DO m= 0,OTMF(ifs)
                              dIdT(m,ifr,ifs,iset)= QTOTI*
     1                                   dIdT(m,ifr,ifs,iset)*SF(iset)
                              ENDDO
                          ENDIF
                      ENDDO
                  ENDDO
              ENDIF
          CALL ADD(iset,NFREQ(iset),NFS,DTYPE(iset),CN,CD,TOTFS,RPD,
     1                                               OTMF,dIdT,OUTPUT)
          IF((IFRPW(iset).NE.0).AND.(printyn.EQ.1)) THEN
              IF(SF(iset).NE.1.d0) WRITE(6,617) SF(iset)
              WRITE(6,699)
              IF(DTYPE(iset).EQ.1) THEN
                  WRITE(6,614) iset,TEMP(iset),(ifs,ifs=1,NFS)
                  WRITE(6,699)
                  DO ifr=1,NFREQ(iset)
                      WRITE(6,616) FREQ(ifr,iset),WAVL(ifr,iset),
     1                OUTPUT(ifr,iset),(TOTFS(ifr,ifs,iset),ifs=1,NFS)
                      ENDDO
                  IF(FITIT.GT.0) THEN
                      WRITE(11,1100) iset,TEMP(iset),
     1                                      INFO(iset),(ifs,ifs=1,NFS)
                      WRITE(11,699)
                      WRITE(12,1201) INFO(iset)
                      ENDIF
                  ENDIF
              IF(DTYPE(iset).EQ.2) THEN
                  WRITE(6,615) iset,TEMP(iset),
     1                                   (CN(iset,ifs),ifs,ifs= 1,NFS)
                  WRITE(6,623) (CD(iset,ifs),ifs,ifs= 1,NFS)
                  WRITE(6,622) (ifs,ifs=1,NFS)
                  WRITE(6,699)
                  DO ifr=1,NFREQ(iset)
                      WRITE(6,618) FREQ(ifr,iset),WAVL(ifr,iset),
     1                OUTPUT(ifr,iset),(TOTFS(ifr,ifs,iset),ifs=1,NFS)
                      ENDDO
                  IF(FITIT.GT.0) THEN
                      WRITE(11,1101) iset,TEMP(iset),INFO(iset),
     1                                                 (ifs,ifs=1,NFS)
                      WRITE(11,699)
                      WRITE(12,1201) INFO(iset)
                      ENDIF
                  ENDIF
              WRITE(6,698)
              IF(FITIT.GT.0) WRITE(11,699)
              DO ifr=1,NFREQ(iset)
                  IF(FITIT.GT.0) THEN
                      WRITE(11,1102) FREQ(ifr,iset),
     1                   OBS(ifr,iset),UNC(ifr,iset),OUTPUT(ifr,iset),
     2                 (OUTPUT(ifr,iset)-OBS(ifr,iset))/UNC(ifr,iset),
     3                                 (TOTFS(ifr,ifs,iset),ifs=1,NFS)
                      WRITE(12,1102) FREQ(ifr,iset),
     1                   OBS(ifr,iset),UNC(ifr,iset),OUTPUT(ifr,iset),
     2                                 (TOTFS(ifr,ifs,iset),ifs=1,NFS)
                      ENDIF
                  ENDDO
              IF(FITIT.GT.0) WRITE(11,698)
              ENDIF
 1000     CONTINUE
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT(' Predissociation Rates [s-1]'//3x,'v',3x,'J"',3x,
     1 'E(v,J)',4x,'TOTAL RATE',5x,'State-',I1,4(8x:'State-',I1))
  602 FORMAT(' *** WARNING ... while searching for  v=',I3,', J=',i3,
     1 '  with  EVIN=',f9.2/5x,'actually find  v=',I3,'  at  E=',F9.2,
     2 '  Reset  EVIN=',f10.2,'  & try again.')
  604 FORMAT(' For   E(v=',I2,' J=',I3,') =',F11.4,',  <coupling ',
     1 'function> =',G16.8)
  606 FORMAT(1x,i3,1x,i3,f11.3,6(1PD14.6))
  608 FORMAT(1x,'Truncate  v=',I2,'  thermal  T=',f7.2,'K  rotational su
     1m at  E(v,J=',i3,')=',f8.2,'[wn]')
  609 FORMAT(' Calculation is for initial-state level   v=',i3,'   J=',
     1  I3)
  610 FORMAT(' Initial-state level v =',I2,' has rotational sublevels'/
     1(4(3x:'E(',I3,')=',F10.3)))
  614 FORMAT(' Calculated Transition Intensity Coefficients for set #',
     1  I2,' [T = ',F6.1,']'/' FREQ/cm-1  WAVL/nm      TOTAL',
     2  5(4x:'State-',I1)/)
  615 FORMAT(' Calculated Transition Intensity Ratios for set #',I2,
     1  ' [Temp = ',F6.1,']'/11x,I2,'*I(fs',SS,I1,')',
     2  5(SP,I3,'*I(fs',SS,I1,')':))
  623 FORMAT('   Ratio = ',60('-')/11x,I2,'*I(fs',SS,I1,')',
     1  5(SP,I3,'*I(fs',SS,I1,')':))
  622 FORMAT(' FREQ/cm-1  WAVL(nm)     Ratio',5(4x,'State-',I1:))
  616 FORMAT(f9.2,f10.4,f12.2,6(f11.5))
  617 FORMAT( ' Calculated absorption coefficients multiplied by scaling
     1 factor:',f9.5)
  618 FORMAT(f9.2,f10.4,f12.6,6(f11.5))
  620 FORMAT(' At T=',f7.2,'K population fraction factors for',I3,
     1 ' vibrational levels containing'/5x,'population fraction',f11.8,
     2 ' are:   Fvib(v=',I2,')=',f12.9/(43x,'Fvib(v=',I2,')=',f12.9))
  690 FORMAT(1x,a70)
  697 FORMAT(/40('=='))
  698 FORMAT(40('=='))
  699 FORMAT(40('--'))
 1100 FORMAT(/' Calculated Transition Intensity Coefficients for  ',
     1 'Data Set #',I2,' [Temp = ',F6.1,']'/1x,a70/' FREQ(cm-1) ',
     2 '   OBS     UNC     CALC  [c-o]/u',5(2x:'State-',I1)/)
 1101 FORMAT(/' Calculated Transition Intensity Ratios for  ',
     1 'Data Set #',I2,' [Temp = ',F6.1,']'/1x,a70/' FREQ(cm-1) ',
     2 '   OBS     UNC     CALC  [c-o]/u',5(2x:'State-',I1)/)
 1102 FORMAT(f10.2,f9.3,f8.3,f9.3,f8.3,6(f9.3:))
 1105 FORMAT(/1x,79('-')/1x,'Calculated Predissociation Rates'//3x,
     1 'v',3x,'J"',3x,'E(v,J)',4x,' OBS RATE ',4x,'CALC RATE',6x,
     2 'State-',I1,4(8x:'State-',I1))
 1106 FORMAT(1x,i3,1x,i3,f11.3,7(1PD14.6))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c***********************************************************************
      SUBROUTINE ADD(iset,n,nfs,type,cn,cd,totfs,rpd,otmf,dIdT,output)
c-----------------------------------------------------------------------
c Add individual intensity contributions to get total sum or ratio value
c-----------------------------------------------------------------------
ccc   INCLUDE 'arrsizes.h'
c-----------------------------------------------------------------------
c  Utility routine to summarize dimensioning of arrays
c-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,
     1        mxsets,mxfreq
      REAL*8 CCM,PI
c-----------------------------------------------------------------------
c  mxdata - maximum number of input data points
c  mxisp  - maximum number of points for initial state potential array
c           (also used for number of points in transition moment array)
c  mxnj    - maxiumum value of j quantum number allowed
c  mxfsp  - maximum number of points for final state potential array
c  mxnp   - maximum number of parameters total
c  mxntp  - maximum number of turning points to be read in
c  mxprm  - maximum number of parameters for final state pot'l or TMF
c  mxv    - largest value for the v quantum number
c  mxfs   - maximum number of final states allowed
c  mxisot - maximum number of isotopomers allowed
c  mxsets - maximum number of data sets allowed
c  mxfreq - maximum number of data points allowed in a given set
c-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
      INTEGER i,ifs,m,n,nfs,cd(mxsets,mxfs),cn(mxsets,mxfs),rpd,
     1        iset,otmf(mxfs),type
      REAL*8 dIdT(0:mxprm-1,mxfreq,mxfs,mxsets),
     1   output(mxfreq,mxsets),sumnum,sumden,
     2   totfs(mxfreq,mxfs,mxsets)
c
      DO i= 1,n
          sumnum= 0.d0
          sumden= 0.d0
          DO ifs= 1,nfs
              IF(type.EQ.1) THEN
                  sumnum= sumnum + cn(iset,ifs)*totfs(i,ifs,iset)
                ELSEIF(type.EQ.2) THEN
                  sumnum= sumnum + cn(iset,ifs)*totfs(i,ifs,iset)
                  sumden= sumden + cd(iset,ifs)*totfs(i,ifs,iset)
                ENDIF
              ENDDO
          IF((type.EQ.2).AND.(rpd.GT.0)) THEN
              DO ifs= 1,nfs
                  DO m= 0,otmf(ifs) 
                      dIdT(m,i,ifs,iset)= dIdT(m,i,ifs,iset)*
     1     (cn(iset,ifs)*sumden - cd(iset,ifs)*sumnum)/(sumden*sumden)
                      ENDDO
                  ENDDO
              ENDIF
          IF(type.EQ.1) THEN
              output(i,iset)= sumnum
            ELSEIF(type.EQ.2) THEN
              output(i,iset)= sumnum/sumden
            ENDIF
          ENDDO
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c  last updated 5 Oct'01 rjl
c=======================================================================
      SUBROUTINE GENFS(ifs)
c=======================================================================
c** Subroutine to prepare final-state potential: currently 3 options:
c* FSTYPE= 1 :  Repulsive exponential with variable exponent
c          V = VLIMF + A1*exp{-(r- Rexfs)*[A2 + A3*y + A4*y^2 ...]}
c* FSTYPE= 2 :  Numerical potential@ longer range with the 2 innermost
c          points defining B0 & B1 in expansion
c          V = B0 + B1*exp{-(r- Rexfs)*[A2 + A3*y +A4*y^2 ...]}
c* FSTYPE= 3 :  EMO potential
c          V = VLIMF + A1*[1 + exp{-(r- A2))*[A3 + A4*y + ...]}]^2 - A1
c-----------------------------------------------------------------------
c  inputs  FSPRM(j,ifs) - final state parameters (one set per final
c                         state)
c          FSTYPE(ifs)  - specifies the type of final state potential
c          ifs          - final state counter
c          mxfsp        - number of final state potential points
c          NFSPRM(ifs)  - number of final state parameters
c          RAD(i)       - radial distance array
c          REXFS        - point about which potential is expanded
c          RTP(i,ifs)   - (for FSTYPE=2) radial distance for turning
c                                        points i= 1 or 2
c          VTP(i,ifs)   - (for FSTYPE=2) potential value for turning
c                                        points i= 1 or 2
c          XCOORD(ifs)  - expansion COORDinate y for final-state potenl.
c                         XCOORD = p (p=1-9)
c                            for y= zfs= (r^p - REXFS^p)/(r^p + REXFS^p)
c                         XCOORD = 10 for y= zfs= (r-REXFS)/r
c                         XCOORD = 11 for y= zfs= (r-REXFS)/(REXFS)
c
c  outputs  VF(i,ifs)   - final state potential array
c-----------------------------------------------------------------------
ccc   INCLUDE 'arrsizes.h'
c-----------------------------------------------------------------------
c  Utility routine to summarize dimensioning of arrays
c-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,
     1        mxsets,mxfreq
      REAL*8 CCM,PI
c-----------------------------------------------------------------------
c  mxdata - maximum number of input data points
c  mxisp  - maximum number of points for initial state potential array
c           (also used for number of points in transition moment array)
c  mxnj    - maxiumum value of j quantum number allowed
c  mxfsp  - maximum number of points for final state potential array
c  mxnp   - maximum number of parameters total
c  mxntp  - maximum number of turning points to be read in
c  mxprm  - maximum number of parameters for final state pot'l or TMF
c  mxv    - largest value for the v quantum number
c  mxfs   - maximum number of final states allowed
c  mxisot - maximum number of isotopomers allowed
c  mxsets - maximum number of data sets allowed
c  mxfreq - maximum number of data points allowed in a given set
c-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
c-----------------------------------------------------------------------
      INTEGER FSTYPE(mxfs),i,ifs,IP,j,NFSPRM(mxfs),NUMIMP,XCOORD(mxfs)
c-----------------------------------------------------------------------
      REAL*8 AVALUE(mxfs),BETA,EF,ES,EXPON1,EXPON2,BVALUE(mxfs),
     1    FSPRM(mxprm,mxfs),numer,R1,R2,RAD(mxfsp),REXFS,REQFS,RH,RMIN,
     2    RTP(2,mxfs),VF(mxfsp,mxfs),VLIMF(mxfs),VTP(2,mxfs),
     3    zfs(mxfsp,mxfs),z1,z2,ZZ,ZZ1,ZZ2
c-----------------------------------------------------------------------
      COMMON /MGf/ AVALUE,BVALUE,REXFS,RTP,VTP,zfs,FSTYPE,XCOORD
      COMMON /MFGf/ VF,VLIMF
      COMMON /MDGf/ FSPRM,NFSPRM
      COMMON /MFGtGf/ RAD
c-----------------------------------------------------------------------
      RMIN= RAD(1)
      RH= RAD(2)-RAD(1)
c
      IP= XCOORD(ifs)
      IF(FSTYPE(ifs).EQ.1) THEN
c-----------------------------------------------------------------------
c  Final state type 1 has a purely repulsive analytic potential form:
c      VF = VLIMF + A1 * exp{-(R-REXFS)*[A2 + A3*y + ... + A6*y^4};
c  Parameters A1-A6 are read in as FSPRM(i), i= 1-6 respectively.
c  Note that VLIMF is added to this potential in main program.
c-----------------------------------------------------------------------
          DO  i=1, mxfsp
              numer= RAD(i)-REXFS
              IF(IP.EQ.0) zfs(i,ifs)= numer/REXFS       
              IF(IP.EQ.1) zfs(i,ifs)= numer/(RAD(i)+REXFS)
              IF(IP.GE.2) zfs(i,ifs)= (RAD(i)**IP - REXFS**IP)/
     1                                        (RAD(i)**IP + REXFS**IP)
              BETA= 0.d0
              ZZ= 1.d0
              DO j= 2,NFSPRM(ifs)
                  BETA= BETA + FSPRM(j,ifs)*ZZ
                  ZZ= ZZ*zfs(i,ifs)
                  ENDDO
              VF(i,ifs)= VLIMF(ifs) + FSPRM(1,ifs)*DEXP(-(numer)*BETA)
              ENDDO
          ENDIF
      IF(FSTYPE(ifs).EQ.2) THEN
c-----------------------------------------------------------------------
c  Final State 2 is defined by NTPFS turning points [RTPF(i),VTPF(i)]
c  with a repulsive exponential inner wall attached to the 2 innermost
c  points.  Outer portion was previously  obtained by interpolating (in
c  GENINT) to get potential array with piecewise polynomials or cubic
c  splines.  VSHFTF (cm-1) was added to the read-in potential points to
c  make them consistent with the stated VLIMF value.
c  Units for distance and energy are Angstroms and cm-1 repectively.
c-----------------------------------------------------------------------
c  Repulsive inner wall is defined by fitting the 2 innermost read-in
c  turning points to the form:
c    VF(R) = A + B exp[-(R-REXFS)*(b0 + b1*y + b2*y^2 + ... + b5*y^5)
c            if XCOORD = 1   y = zfs(R) = (R-REXFS)/(R+REXFS)
c            if XCOORD = p(=1-9)   y = zfs(R) = (R^p - REXFS^p)/
c                                                (R^p + REXFS^p)
c            if XCOORD = 10  y = zfs(R) = (R-REXFS)/R
c            if XCOORD = 11  y = zfs(R) = (R-REXFS)/REXFS
c  where A and B are determined by the 2 turning points
c  and parameters b(j)= FSPRM[(j+1),ifs]
c-----------------------------------------------------------------------
          R1= RTP(1,ifs) - REXFS
          R2= RTP(2,ifs) - REXFS
          IF (XCOORD(ifs).EQ.10) THEN
              z1= R1/RTP(1,ifs)
              z2= R2/RTP(2,ifs)
            ELSEIF (XCOORD(ifs).EQ.11) THEN
              z1= R1/REXFS
              z2= R2/REXFS
            ELSEIF (XCOORD(ifs).EQ.1) THEN
              z1= R1/(RTP(1,ifs) + REXFS)
              z2= R2/(RTP(2,ifs) + REXFS)
            ELSEIF ((XCOORD(ifs).GE.2).AND.(XCOORD(ifs).LE.9)) THEN
              z1=(RTP(1,ifs)**IP -REXFS**IP)/(RTP(1,ifs)**IP +REXFS**IP)
              z2=(RTP(2,ifs)**IP -REXFS**IP)/(RTP(2,ifs)**IP +REXFS**IP)
            ENDIF
c-----------------------------------------------------------------------
c  Now solve two equations & two unknowns to get A and B parameters for
c  inner wall of final state 2 (force exponential to go through the
c  first two turning points).  A= AVALUE and B= BVALUE below.
c-----------------------------------------------------------------------
          ZZ1= 1.d0
          ZZ2= 1.d0
          EXPON1= 0.d0
          EXPON2= 0.d0
          DO j=1,NFSPRM(ifs)
              EXPON1= EXPON1 + FSPRM(j,ifs)*ZZ1
              EXPON2= EXPON2 + FSPRM(j,ifs)*ZZ2
              ZZ1= ZZ1*z1
              ZZ2= ZZ2*z2
              ENDDO
          ES= DEXP(-R1*EXPON1)
          EF= DEXP(-R2*EXPON2)
          BVALUE(ifs)= (VTP(1,ifs)-VTP(2,ifs))/(ES-EF)
          AVALUE(ifs)= VTP(1,ifs)-BVALUE(ifs)*ES
c-----------------------------------------------------------------------
c  NUMIMP is NUMber of Inner Mesh Points (for repulsive inner wall)
c-----------------------------------------------------------------------
          NUMIMP= IDNINT((RTP(2,ifs) - RMIN)/RH + 1.d0)
          IF(NUMIMP.LE.0) NUMIMP= 1
c-----------------------------------------------------------------------
c  Now add on the exponential inner wall in region up to RTP(2,ifs)
c-----------------------------------------------------------------------
          DO  i=1,NUMIMP
              numer= RAD(i)-REXFS
              IF(XCOORD(ifs).EQ.1) zfs(i,ifs)= numer/(RAD(i)+REXFS)
              IF((XCOORD(ifs).GE.2).AND.(XCOORD(ifs).GE.9))
     1     zfs(i,ifs)= (RAD(i)**IP -REXFS**IP)/(RAD(i)**IP +REXFS**IP)
              IF(XCOORD(ifs).EQ.10) zfs(i,ifs)= numer/RAD(i)
              IF(XCOORD(ifs).EQ.11) zfs(i,ifs)= numer/REXFS
              BETA= 0.d0
              ZZ= 1.d0
              DO j= 1,NFSPRM(ifs)
                  BETA= BETA + FSPRM(j,ifs)*ZZ
                  ZZ= ZZ*zfs(i,ifs)
                  ENDDO
              VF(i,ifs)= AVALUE(ifs) + BVALUE(ifs)*DEXP(-(numer)*BETA)
              ENDDO
          ENDIF
c
      IF(FSTYPE(ifs).EQ.3) THEN
c-----------------------------------------------------------------------
c** Generate an Extended Morse Oscillator final-state potential
c   V = VLIMF + A1*[1 + exp{-(R- A2)*[A3 + A4*y + A5*y^2 +...]}]^2 - A1
c* Parameters A1-A6 are read in as FSPRM(i), i= 1-6 respectively, while
c  VLIMF is added to this potential in main program.
c-----------------------------------------------------------------------
          DO  i=1, mxfsp
              REQFS= FSPRM(2,ifs)
              numer= (RAD(i) - REQFS)
              zfs(i,ifs)= (RAD(i)**IP - REXFS**IP)/
     1                                        (RAD(i)**IP + REXFS**IP)
              BETA= 0.d0
              ZZ= 1.d0
              DO j= 3,NFSPRM(ifs)
                  BETA= BETA + FSPRM(j,ifs)*ZZ
                  ZZ= ZZ*zfs(i,ifs)
                  ENDDO
              VF(i,ifs)= VLIMF(ifs) + FSPRM(1,ifs)*
     1                         ((DEXP(-numer*BETA)- 1.d0)**2 - 1.d0)
              ENDDO
          ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE GENTMF(ifs)
c=======================================================================
c   Last modified   6 October 2001
c=======================================================================
ccc   INCLUDE 'arrsizes.h'
c-----------------------------------------------------------------------
c  Utility routine to summarize dimensioning of arrays
c-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,
     1        mxsets,mxfreq
      REAL*8 CCM,PI
c-----------------------------------------------------------------------
c  mxdata - maximum number of input data points
c  mxisp  - maximum number of points for initial state potential array
c           (also used for number of points in transition moment array)
c  mxnj    - maxiumum value of j quantum number allowed
c  mxfsp  - maximum number of points for final state potential array
c  mxnp   - maximum number of parameters total
c  mxntp  - maximum number of turning points to be read in
c  mxprm  - maximum number of parameters for final state pot'l or TMF
c  mxv    - largest value for the v quantum number
c  mxfs   - maximum number of final states allowed
c  mxisot - maximum number of isotopomers allowed
c  mxsets - maximum number of data sets allowed
c  mxfreq - maximum number of data points allowed in a given set
c-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
      INTEGER I,ifs,IP,IR2TMF,J,LNPT,LPTMF,ILRTMF,m,NCNTMF,NIN,NPRS,
     1    NPRF,NPTMF,NUSETMF,NROW,  GFS(mxfs),OTMF(mxfs),TMFTYP(mxfs)
      REAL*8 CNNTMF,FCT(mxisp),MFACTMF,RAD(mxfsp),REXTMF,RFACTMF,
     1       TMF(mxfsp),TMFLIM,TMFPRM(0:mxprm-1,mxfs),Xi(mxntp),xtmf,
     2       XM2(mxfsp),Yi(mxntp),ztmf(mxisp,mxfs)
c-----------------------------------------------------------------------
      COMMON /MGt/ REXTMF
      COMMON /MFGt/ XM2,ztmf,LPTMF,NIN,TMFTYP
      COMMON /MFGtGf/ RAD
      COMMON /MFDGt/ TMFPRM,OTMF,GFS
c-----------------------------------------------------------------------
      DATA LNPT/1/,NPRS/1/,IR2TMF/0/
c-----------------------------------------------------------------------
      IF(TMFTYP(ifs).EQ.0) THEN
c-----------------------------------------------------------------------
c** Read   NPTMP  points of transition moment function with asymptotic
c                 value of  TMFLIM
c* Interpolate with NUSETMF-point piecewise polynomials (or splines for
c    NUSETMF.le.0), which are extrapolated to the asymptote as specified by
c    parameters ILRTMF, CNC & CNN (see read #20).
c RFACT - factor converts read-in distances to angstroms
c MFACT - factor converts read-in moment values to debye (for absorption
c         or emission or cm-1 for predissociation).
c=======================================================================
          READ(5,*) NPTMF, TMFLIM
          READ(5,*) NUSETMF, ILRTMF, NCNTMF, CNNTMF
          READ(5,*) RFACTMF, MFACTMF
          READ(5,*) (XI(i), YI(i), i=1, NPTMF)
c=======================================================================
          WRITE(6,610) NPTMF, TMFLIM
          IF(NUSETMF.GT.0) WRITE(6,616) NUSETMF, NPTMF
          IF(NUSETMF.LE.0) WRITE(6,618) NPTMF
          IF((ILRTMF.GT.1).AND.(DABS(CNNTMF).GT.0.D0))
     1                                     WRITE(6,610) CNNTMF, NCNTMF
          WRITE(6,622) RFACTMF, MFACTMF
          NROW= (NPTMF+2)/3
          DO  J= 1,NROW
              WRITE(6,624) (XI(I), YI(I), I=J, NPTMF, NROW)
              ENDDO
          DO  I= 1,NPTMF
              XI(I)= XI(I)*RFACTMF
              YI(I)= YI(I)*MFACTMF
              ENDDO
  616 FORMAT(' Perform',I3,'-point piecewise polynomial interpolation ov
     1er',I5,' input points' )
  618 FORMAT(' Perform cubic spline interpolation over the',I5,' input p
     1oints' )
  622 FORMAT(' Scale input points:  (distance)*',1PD16.9,'   &  (moment)
     1*',D16.9/4x,'to get required units  [Angstroms & debye (or cm-1 fo
     2r predissociation)]'/
     3  3('      X(i)         Y(i)  ')/3(3X,11('--')))
  624 FORMAT((3(F12.6,F13.6)))
          NPRF= NIN
          CALL GENINT(LNPT,NIN,RAD,FCT,NUSETMF,IR2TMF,NPTMF,XI,YI,
     1                          TMFLIM,ILRTMF,NCNTMF,CNNTMF,NPRS,NPRF)
          DO i=1,NIN
              ztmf(i,ifs)= FCT(i)
              ENDDO
        ELSEIF((TMFTYP(ifs).GE.1).AND.(TMFTYP(ifs).LE.9)) THEN
          IP= TMFTYP(ifs)
          WRITE(6,614) IP,IP,IP,IP,REXTMF
          DO i=1,NIN
              ztmf(i,ifs)= (RAD(i)**IP - REXTMF**IP)/
     1                                       (RAD(i)**IP + REXTMF**IP)
              ENDDO
        ELSEIF(TMFTYP(ifs).EQ.10) THEN
          WRITE(6,623) REXTMF
          DO i=1,NIN
              ztmf(i,ifs)= (RAD(i) - REXTMF)/RAD(i)
              ENDDO
        ELSEIF(TMFTYP(ifs).EQ.11) THEN
          WRITE(6,613) REXTMF
          DO i=1,NIN
              ztmf(i,ifs)= (RAD(i) - REXTMF)/REXTMF
              ENDDO
        ELSEIF(TMFTYP(ifs).EQ.12) THEN
          WRITE(6,612)
          DO i= 1,NIN
              ztmf(i,ifs)= XM2(i)
              ENDDO
        ELSEIF(TMFTYP(ifs).EQ.13) THEN
          WRITE(6,611)
          DO i= 1,NIN
              ztmf(i,ifs)= RAD(i)
              ENDDO
        ELSEIF(TMFTYP(ifs).EQ.14) THEN
          WRITE(6,621)
          DO i= 1,NIN
              ztmf(i,ifs)= 2.d0
              ENDDO
        ELSE
          WRITE(6,601) TMFTYP(ifs)
          STOP
        ENDIF
      RETURN
c-----------------------------------------------------------------------
  601 FORMAT(/' *Fatal Error in GENTMF* TMFTYP =',i3,' unallowed ',
     1 '- must be .LE. 14')
  610 FORMAT(' Transition moment function defined by interpolating over'
     1  ,I4,' read-in points'/5x,'and approaching the asymptotic value',
     2  f12.6)
  621 FORMAT(' Transition moment function is a power series in the CONST
     1ANT 2.d0 {for testing}')
  611 FORMAT(' Transition moment function is a power series in the dista
     1nce  R [Angstroms]')
  612 FORMAT(' Transition moment function is a power series in the varia
     1ble   z = 1/R**2')
  613 FORMAT(' Transition moment function is a power series in the Dunha
     1m variable'/4x,'z = (R - REXTMF)/REXTMF   with   REXTMF =',F8.5,
     2 ' [Angstroms]')
  623 FORMAT(' Transition moment function is a power series in the SPF v
     1ariable'/4x,'z = (R - REXTMF)/R   with   REXTMF=',F8.5,
     2 ' [Angstroms]')
  614 FORMAT(' Transition moment function is a power series in the Surku
     1s variable'/4x,'z = (R^',i1,' - REXTMF^',i1,')/(R^',i1,
     2 ' + REXTMF^',i1,')   with   REXTMF=',F8.5,' [Angstroms]')
  820 FORMAT(' Coefficients of power series expansion for transition mom
     1ent function are:'/(3x,6(f12.8):))
c-----------------------------------------------------------------------
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DYIDPJ(i,NDATA,NPTOT,YC,PV,PD,PS,RMSR)
c=======================================================================
c   Last modified  29 April 2001
c=======================================================================
c  i     - datum whose derivative is required
c  nptot - total number of parameters to be varied
c  yc    - current Ycalc value
c  pv    - parameter array, pv(n)
c  pd    - partial derivative array, pd(i)
c  ps    - parameter sensitivities, ps(i)
c  rmsr  - root mean square residuals
c-----------------------------------------------------------------------
ccc   INCLUDE 'arrsizes.h'
c-----------------------------------------------------------------------
c  Utility routine to summarize dimensioning of arrays
c-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,
     1        mxsets,mxfreq
      REAL*8 CCM,PI
c-----------------------------------------------------------------------
c  mxdata - maximum number of input data points
c  mxisp  - maximum number of points for initial state potential array
c           (also used for number of points in transition moment array)
c  mxnj    - maxiumum value of j quantum number allowed
c  mxfsp  - maximum number of points for final state potential array
c  mxnp   - maximum number of parameters total
c  mxntp  - maximum number of turning points to be read in
c  mxprm  - maximum number of parameters for final state pot'l or TMF
c  mxv    - largest value for the v quantum number
c  mxfs   - maximum number of final states allowed
c  mxisot - maximum number of isotopomers allowed
c  mxsets - maximum number of data sets allowed
c  mxfreq - maximum number of data points allowed in a given set
c-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
      INTEGER cycle,i,ifr,ifs,iset,m,n,idata,LPDER,NDATA,NFS,nprm,
     1  NPTOT,NPPFREE,NSETS,RPD, GFS(mxfs),
     2  FSVAR(mxprm,mxfs),IFRPW(mxsets),NFREQ(mxsets),NFSPRM(mxfs),
     3  NVJ(mxsets),OTMF(mxfs),SCALE(mxsets),TMFVAR(0:mxprm-1,mxfs),
     4  UPPER(mxsets)
      REAL*8 DELTAP(mxnp),dIdT(0:mxprm-1,mxfreq,mxfs,mxsets),
     1       DTIVE(mxdata,mxnp),FSPRM(mxprm,mxfs),
     2       OUTPUT(mxfreq,mxsets),PD(mxnp),DFACT,PS(mxnp),
     3       PV(mxnp),RMSR,SF(mxsets),TMFPRM(0:mxprm-1,mxfs),YC,
     4       YCALC(mxdata)
      CHARACTER*70 INFO(mxsets)
c-----------------------------------------------------------------------
      COMMON /MD/ DFACT,FSVAR,LPDER,NPPFREE
      COMMON /MFD/ dIdT,OUTPUT,SF,IFRPW,NFREQ,NFS,NSETS,NVJ,RPD,SCALE,
     1             TMFVAR,INFO
      COMMON /MDGf/ FSPRM,NFSPRM
      COMMON /MFDGt/ TMFPRM,OTMF,GFS
c-----------------------------------------------------------------------
      SAVE DTIVE,YCALC,cycle
      DATA cycle/0/
c-----------------------------------------------------------------------
c  When called for first datum, create complete partial derivative array
c  to will be used in the subsequent calls for rest of the data
c-----------------------------------------------------------------------
      IF(i.EQ.1) THEN
          cycle= cycle + 1
          IF(NPPFREE.GT.0) WRITE(6,600)
          DO idata= 1,NDATA
              DO nprm= 1,NPTOT
                  DTIVE(idata,nprm)= 0.d0
                  ENDDO
              ENDDO
          nprm= 0
c-----------------------------------------------------------------------
c  First copy current PV values onto internal potential and transition
c  moment function parameter arrays
c-----------------------------------------------------------------------
          DO ifs= 1,NFS
              DO n= 1,NFSPRM(ifs)
                  IF(FSVAR(n,ifs).GT.0) THEN
                      nprm= nprm + 1
                      FSPRM(n,ifs)= PV(nprm)
                      IF(RMSR.GT.0.d0) THEN
                          DELTAP(nprm)= DFACT*PS(nprm)*NPTOT/RMSR
c    1                      *DSQRT(DBLE(NDATA-NPTOT)/DBLE(NDATA))/RMSR
                        ELSE
                          DELTAP(nprm)= 0.0010d0*FSPRM(n,ifs)
c-----------------------------------------------------------------------
c  DELTAP - step size used to define derivatives-by-differences
c-----------------------------------------------------------------------
                        ENDIF
                      WRITE(6,601) nprm, DELTAP(nprm)
                      IF(LPDER.GT.0) WRITE(10,601) nprm, DELTAP(nprm)
                      ENDIF
                  ENDDO
              CALL GENFS(ifs)
              DO m= 0,OTMF(ifs)
                  IF(TMFVAR(m,ifs).GT.0) THEN
                      nprm= nprm + 1
                      TMFPRM(m,ifs)= PV(nprm)
                      ENDIF
                  ENDDO
              ENDDO
          DO iset= 1,NSETS
              IF(SCALE(iset).EQ.1) THEN
                  nprm= nprm + 1
                  SF(iset)= PV(nprm)
                  ENDIF
              ENDDO
          RPD= 1
c-----------------------------------------------------------------------
c  RPD is a flag that tells FORWARD whether or not to return transition
c  moment function partial derivatives (RPD= 1 to return  derivatives)
c-----------------------------------------------------------------------
          CALL FORWARD
          RPD= 0
          IDATA= 0
          DO iset= 1,NSETS
              IF(IFRPW(iset).EQ.0) UPPER(iset)= NVJ(iset)
              IF(IFRPW(iset).NE.0) UPPER(iset)= NFREQ(iset)
              DO ifr= 1,UPPER(iset)
                  idata= idata + 1
                  YCALC(idata)= OUTPUT(ifr,iset)
                  ENDDO
              ENDDO
c-----------------------------------------------------------------------
c  begin partial derivative calculation for final state potential
c  parameters here - do SYMMETRIC partial derivatives-by-differences
c-----------------------------------------------------------------------
          nprm= 0
          DO ifs= 1,NFS
              DO n= 1,NFSPRM(ifs)
                  IF(FSVAR(n,ifs).GT.0) THEN
                      nprm= nprm + 1
                      FSPRM(n,ifs)= PV(nprm) + DELTAP(nprm)
                      CALL GENFS(ifs)
                      CALL FORWARD
                      idata= 0
                      DO iset= 1,NSETS
                          DO ifr= 1, UPPER(iset)
                              idata= idata + 1
                              DTIVE(idata,nprm)= (OUTPUT(ifr,iset)-
     1                                      YCALC(idata))/DELTAP(nprm)
                              ENDDO
                          ENDDO
                      FSPRM(n,ifs)= PV(nprm) - 2*DELTAP(nprm)
                      CALL GENFS(ifs)
                      CALL FORWARD
                      idata= 0
                      DO iset= 1,NSETS
                          DO ifr= 1, UPPER(iset)
                              idata= idata + 1
                              DTIVE(idata,nprm)=0.5d0*(DTIVE(idata,nprm)
     1                 + (YCALC(idata)-OUTPUT(ifr,iset))/DELTAP(nprm))
                              ENDDO
                          ENDDO
                      FSPRM(n,ifs)= PV(nprm)
                      ENDIF
                  ENDDO
              CALL GENFS(ifs)
              DO m= 0,OTMF(ifs)         
                  IF(TMFVAR(m,ifs).GT.0) THEN
                      nprm= nprm + 1
                      idata= 0
                      DO iset= 1,NSETS
                          DO ifr= 1,UPPER(iset)
                              idata= idata + 1
                              DTIVE(idata,nprm)= dIdT(m,ifr,ifs,iset)
                              ENDDO
                          ENDDO
                      ENDIF
                  ENDDO
              ENDDO
          idata= 0
          DO iset= 1,NSETS
              IF(SCALE(iset).EQ.0) THEN
                  idata= idata + UPPER(iset)
                ELSE
                  nprm= nprm + 1
                  DO ifr= 1,UPPER(iset)
                      idata= idata + 1
                      DTIVE(idata,nprm)= OUTPUT(ifr,iset)/SF(iset)
                      ENDDO
                  ENDIF
              ENDDO
          IF(LPDER.GT.0) THEN
              WRITE(10,100) cycle
              DO idata= 1,NDATA
                  WRITE(10,101) YCALC(idata),(DTIVE(idata,nprm),
     1                                                  nprm= 1,NPTOT)
                  ENDDO
              ENDIF
          ENDIF     
c-----------------------------------------------------------------------
c  For each datum i, collect yalculated values and partial derivatives
c  below
c-----------------------------------------------------------------------
      YC= YCALC(i)
      DO nprm= 1,NPTOT
          PD(nprm)= DTIVE(i,nprm)
          ENDDO
      RETURN
  600 FORMAT(/' Change in Potential Parameter Values (Derivatives-',
     1 'by-Differences)')
  601 FORMAT(1x,'DELTAP(',i2,')= ',f12.6)
  100 FORMAT(/5x,'YCALC',10x,'Partial Derivatives for fitting cycle',i3)
  101 FORMAT(1x,30(1PD14.6))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE OVRLAP(BFCT,DER,EFN,OVR,OVRCRT,PSI,RH,RMIN,TMFPRM,VJ,
     1                  VLIM,z,ifs,IWR,JP,NEND,OTMF,TMFTYP)
c=======================================================================
c        Routine by R.J. Le Roy;  Last Modified 19 August 2003
c=======================================================================
c Calculate overlap integral Franck-Condon Moment FCM(i) array 
c between the given bound state wave function PSI(i) (which is zero
c for  i > NEND) and the J' = JP continuum final state wave function
c (asymptotocally normalized to unit amplitude) at energy EFN on the
c effective potential VJ(I) with asymptote VLIM, with input array
c z(i).  z(i) is the input (radial) array whose moments are being taken.
c NOTE that z(i) = zin(i) = ztmf
c
c On entry, energy units for EFN and VLIM are (cm-1), while VJ(I)
c incorporates the factor BFCT (i.e., VJ/BFCT has units cm-1).
c
c Convergence of asymptotic wave function normalization defined by
c requirement that W.K.B. fits to 3 successive maxima must agree
c relatively to within OVRCRT. 
c-----------------------------------------------------------------------
ccc   INCLUDE 'arrsizes.h'
c-----------------------------------------------------------------------
c  Utility routine to summarize dimensioning of arrays
c-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,
     1        mxsets,mxfreq
      REAL*8 CCM,PI
c-----------------------------------------------------------------------
c  mxdata - maximum number of input data points
c  mxisp  - maximum number of points for initial state potential array
c           (also used for number of points in transition moment array)
c  mxnj    - maxiumum value of j quantum number allowed
c  mxfsp  - maximum number of points for final state potential array
c  mxnp   - maximum number of parameters total
c  mxntp  - maximum number of turning points to be read in
c  mxprm  - maximum number of parameters for final state pot'l or TMF
c  mxv    - largest value for the v quantum number
c  mxfs   - maximum number of final states allowed
c  mxisot - maximum number of isotopomers allowed
c  mxsets - maximum number of data sets allowed
c  mxfreq - maximum number of data points allowed in a given set
c-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
      INTEGER i,ifs,MESH1,MESH2,MESH3,step,first,last,IWR,JP,TURNPT,m,
     1        OTMF,NAMP,NEND,TMFTYP(mxfs)
      REAL*8 AMP1,AMP2,AMP3,AMP4,BFCT,DER(0:mxprm-1),DI,FCFACT,
     1       EFN,ELIM,ER,FCM(0:MXPRM-1),EDIFF1,EDIFF2,EDIFFi,HALF,HARG,
     2       NFACT,
     2       RH,RMIN,OVR,OVRCRT,PSI(NEND),S0,S1,S2,SG1,SG2,SGi,Si,
     3       SNARG,SQKINF,
     3       THIRD,TMFPRM(0:mxprm-1,mxfs),VLIM,VJ(mxfsp),VV,XIITH,XX,
     4       Y1,Y2,Y3,z(mxisp),
     4       ZTST,ZZ0,ZZ1
c-----------------------------------------------------------------------
      HALF=  1.D0/2.D0
      XIITH= 1.D0/12.D0
      THIRD= 1.D0/3.D0
      ER= EFN*BFCT
      ELIM= VLIM*BFCT
      SQKINF= DSQRT(ER-ELIM)
      AMP1= 1.D0
      AMP2= 2.D0
      AMP3= 0.D0
      AMP4= 0.D0
      DO m= 0,OTMF
         DER(m)= 0.D0
         FCM(m)= 0.d0
      ENDDO
c-----------------------------------------------------------------------
c** Locate first turning point and use Airy function to estimate
c  appropriate integration starting point such that  PSI(1) .LE. 1.D-10
c-----------------------------------------------------------------------
      MESH1= 1
      EDIFF1= VJ(MESH1)-ER
      step= DINT(0.05d0/RH)
      IF(step.LT.1) step= 1
      first= step+1
      DO i= first,mxfsp,step
         TURNPT= i
         EDIFF2= VJ(i)-ER
         IF(EDIFF2.LE. 0.D0) GOTO 4
         MESH1= i
         EDIFF1= EDIFF2
      ENDDO     
      IF(IWR.NE.0) THEN
         WRITE(6,607) JP,EFN
         OVR= 0.d0
         RETURN
      ENDIF
    4 MESH2= TURNPT
      TURNPT= MESH1+(MESH2-MESH1)*EDIFF1/(EDIFF1-EDIFF2)
      IF(IABS(TURNPT-MESH2).LE.1) GOTO 6
      IF((TURNPT.LE.0).OR.(TURNPT.GT.mxfsp)) THEN
         IF(IWR.NE.0) WRITE(6,601) JP,EFN
         OVR= 0.d0
         RETURN
      ENDIF
      MESH1= MESH2
      EDIFF1= EDIFF2
      EDIFF2= VJ(TURNPT)-ER
      GOTO 4
    6 DI= 10.D0/(VJ(TURNPT-1)-VJ(TURNPT))**THIRD
      step= DINT(DI)
      MESH1= MAX0(1,TURNPT-step)
      IF(MESH1.GE.NEND) THEN
         OVR= 0.D0
         RETURN
      ENDIF
    8 EDIFF1= VJ(MESH1)-ER
      IF(EDIFF1.LT.10.D0) GOTO 10
c-----------------------------------------------------------------------
c** Adjust starting point outward to ensure integration scheme stability
c-----------------------------------------------------------------------
      MESH1= MESH1+1
      IF((MESH1-mxfsp).LT.0) GOTO 8
      IF((MESH1-mxfsp).GE.0) THEN
         OVR= 0.D0
         RETURN
      ENDIF
   10 MESH2= MESH1+1
      MESH3= MESH2+1
c-----------------------------------------------------------------------
c** WKB starting condition for wave function
c-----------------------------------------------------------------------
      S0= 1.D0
      EDIFF1= VJ(MESH1)-ER
      EDIFFi= VJ(MESH2)-ER
      IF((EDIFF1.GT. 0.D0).AND.(EDIFFi.GT. 0.D0)) THEN
         SG1= DSQRT(EDIFF1)
         SGi= DSQRT(EDIFFi)
         Si= S0*DSQRT(SG1/SGi)*DEXP((SG1+SGi)/2.D0)
         IF(Si.LE.S0) S0= 0.D0
      ELSE
         VV= VJ(MESH2)/BFCT
         IF(IWR.NE.0) WRITE(6,608) JP,EFN,VV,MESH2
         S0= 0.D0
         Si= 1.D0
      ENDIF   
c-----------------------------------------------------------------------
c  notationally speaking, all Yi's refer to values used in the Numerov
c  Algorithm for wavefunciont propagation.
c-----------------------------------------------------------------------
      Y1= S0*(1.D0-XIITH*EDIFF1)
      Y2= Si*(1.D0-XIITH*EDIFFi)
c-----------------------------------------------------------------------
c  Use trapezoid rule for numerical integration.  Initialize FCM(m)
c  values using first section of area.   
c-----------------------------------------------------------------------
      ZZ0= 1.D0
      ZZ1= 1.D0
      DO m= 0,OTMF
         FCM(m)= HALF*S0*PSI(MESH1)*ZZ0 + Si*PSI(MESH2)*ZZ1
         ZZ0= ZZ0*z(MESH1)
         ZZ1= ZZ1*z(MESH2)
      ENDDO
      S2= S0
c-----------------------------------------------------------------------
c** Integrate outward to first turning point.  NOTE that Airy-estimated
c  initialization minimizes need for renormalizations.
c-----------------------------------------------------------------------
      DO 16 i= MESH3,TURNPT
         Y3= Y2+Y2-Y1+EDIFFi*Si
         Y1= Y2
         Y2= Y3
         EDIFFi= VJ(I)-ER
         S1= S2
         S2= Si
         Si= Y3/(1.D0-XIITH*EDIFFi)
c-----------------------------------------------------------------------
c** If bound wavefx. non-negligible, accumulate overlap moments
c  NOTE that FCFACT is the Franck-Condon factor
c-----------------------------------------------------------------------
         IF(I.LE.NEND) THEN
            FCFACT= Si*PSI(i)
            DO m= 0,OTMF
               FCM(m)= FCM(m) + FCFACT
               FCFACT= FCFACT*z(i)
            ENDDO
         ENDIF
c-----------------------------------------------------------------------
c** If wavefuntion too large in forbidden region, renormalize it ...
c-----------------------------------------------------------------------
         IF((Si.GE.1.D32).OR.(i.EQ.TURNPT)) THEN
            NFACT= 1.D0/Si
            Si= 1.D0
            IF(S0.GT.1.D-30) S0= S0*NFACT
            DO m= 0,OTMF
               FCM(m)= FCM(m)*NFACT
            ENDDO
            Y1= Y1*NFACT
            Y2= Y2*NFACT
         ENDIF
   16 CONTINUE
      IF((IWR.NE.0).AND.(S0/SI.GT.1.D-8))
     1                           WRITE(6,602)JP,EFN,MESH1,S0/Si,TURNPT
      MESH2= TURNPT+1
c-----------------------------------------------------------------------
c** If turning point NOT past end of range for bound state wavefx., then
c   integrate from turning point to end of bound-state wave function
c-----------------------------------------------------------------------
      IF(TURNPT.LT.NEND) THEN
         DO i= MESH2,NEND
            Y3= Y2 + Y2 - Y1 + EDIFFi*Si
            Y1= Y2
            Y2= Y3
            EDIFFi= VJ(i)-ER
            S1= S2
            S2= Si
            Si= Y3/(1.D0 - XIITH*EDIFFi)
            FCFACT= Si*PSI(i)
            DO m= 0,OTMF
               FCM(m)= FCM(m) + FCFACT
               FCFACT= FCFACT*z(i)
            ENDDO
         ENDDO   
         MESH2= NEND+1
      ENDIF
c-----------------------------------------------------------------------
c** Continue wave function propagation until amplitude converges
c-----------------------------------------------------------------------
      NAMP= 0
      DO 30 i= MESH2,mxfsp
         Y3= Y2 + Y2 - Y1 + EDIFFi*Si
         Y1= Y2
         Y2= Y3
         EDIFF2= EDIFFi
         EDIFFi= VJ(i)-ER
         S1= S2
         S2= Si
         Si= Y3/(1.D0-XIITH*EDIFFi)
         IF((Si.GE.S2).OR.(S1.GT.S2)) GOTO 30
c-----------------------------------------------------------------------
c** At successive maxima, fit solution to W.K.B. form to determine
c  apparent asymptotic amplitude.
c-----------------------------------------------------------------------
         SG2= DSQRT(-EDIFF2)
         SGi= DSQRT(-EDIFFi)
         HARG= (SG2 + SGi)/2.D0
         SNARG= 1.D0/DSQRT(1.D0 + ((DSQRT(SG2/SGi)*S2/Si- DCOS(HARG))
     1          /DSIN(HARG))**2)
         NAMP= NAMP+1
         AMP4= AMP3
         AMP3= AMP2
         AMP2= AMP1
         AMP1= Si*DSQRT(SGi/SQKINF)/SNARG
         XX= RMIN + (i-1)*RH
         IF(IWR.GE.3) WRITE(6,604) JP,EFN,XX,AMP1
         last= i
c-----------------------------------------------------------------------
c** Test successive amplitudes for convergence
c-----------------------------------------------------------------------
         ZTST= OVRCRT*AMP1
         IF((DABS(AMP1-AMP2).LT.ZTST).AND.(DABS(AMP2-AMP3).LT.ZTST))
     1      GOTO 35
   30 CONTINUE
      IF(IWR.NE.0) WRITE(6,603) JP,EFN,AMP1,AMP2,AMP3,AMP4
   35 OVR= 0.d0
      DO m= 0,OTMF
          FCM(m)= FCM(m)*RH/AMP1
          OVR= OVR + TMFPRM(m,ifs)*FCM(m)
          ENDDO
      IF(IWR.GE.1) WRITE(6,605) JP,EFN,last,XX,(m,FCM(m),m=0,OTMF)
      IF(IWR.GE.2) WRITE(6,606) S0,AMP1,AMP2,AMP3,AMP4
      DO m= 0,OTMF
          DER(m)= 2.d0*OVR*FCM(m)
          ENDDO
      OVR= OVR*OVR
      RETURN
c-----------------------------------------------------------------------
  601 FORMAT(' *** OVRLAP BOMBed ***  For   J = ',I3,'   EFN = ',
     1  F10.2,'   never got to first turning point')
  602 FORMAT(' ** WARNING **  For   J = ', I3 ,'   EFN = ',F10.2,
     1 '  starting wavefunction is PSI(',I4,') = ',D10.3,' and  I(turn.p
     2t.) = ',I4)
  603 FORMAT(' ** WARNING ** For J= ',I3,' EFN= ',F9.2,' amplitude not',
     1 ' converged by end of range'/3x,'Last four values are ',
     2 4(1PD14.6))
  604 FORMAT(' At  J=',I3,'   E=',F9.2,'   R=',F6.3,'  apparent asymptot
     1ic amplitude',1PD14.6)
  605 FORMAT(' At  J=',I3,'   E= ',F12.2,'   R(end)= R(',I5,')=',
     1  F7.4,'    FCM(',I1,')=',F12.8:/
     2  (4x,3(5x,'FCM(',I1,')=',F12.8:)))
  606 FORMAT(5X,'S0= ',1PD10.3,'   & last 4 amplitudes are',2D14.6/
     1    45x,2D14.6)
  607 FORMAT(' *** ERROR ***   At   J = ',I3,'   EFN = ',F10.2,
     1  '  have  V .GT. E  everywhere.' )
  608 FORMAT(' *** Caution ***  For   J = ',I3,'   (EFN= ',F10.2,
     1  ') .GE. (V = ',F10.2,')   at  I = ',I4,' ,so initialize with a n
     2ode.')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE OVRPD(BFCT,DER,EFN,OVR,OVRCRT,RAD,PSI,TMFPRM,
     1                 VJ,VLIM,FITIT,ifs,IWR,JP,LPTMF,NEND,OTMF,TMFVAR)
c-----------------------------------------------------------------------
c ** USE THIS ROUTINE FOR TMFTYP < 0
c           * TMFTYP = -1 for predisociation case where the TMF is not
c                         a power series but an operator instead
c                         P= -hbar^2/2*mu {dW/dR + 2W*d/dR}
c                         where W(R)= a/{4a^2 + (R-Rc)^2}
c-----------------------------------------------------------------------
c Calculate overlap integral, OVR,
c between the given bound state wave function PSI(i) (which is zero
c for  i > NEND) and the J' = JP continuum final state wave function
c (asymptotocally normalized to unit amplitude) at energy EFN on the
c effective potential VJ(I) with asymptote VLIM, with input array
c
c On entry, energy units for EFN and VLIM are (cm-1), while VJ(I)
c incorporates the factor BFCT (i.e., VJ/BFCT has units cm-1).
c
c Convergence of asymptotic wave function normalization defined by
c requirement that W.K.B. fits to 3 successive maxima must agree
c relatively to within OVRCRT. 
c-----------------------------------------------------------------------
c  Original routine by R.J. Le Roy; Modified by G.T. Kraemer; 
c                 Current version 29 April 2001
c ** updating to take 5 points of continuum wavefunction
c             for obtaining derivative, dPSIc/dR
c-----------------------------------------------------------------------
ccc   INCLUDE 'arrsizes.h'
c-----------------------------------------------------------------------
c  Utility routine to summarize dimensioning of arrays
c-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,
     1        mxsets,mxfreq
      REAL*8 CCM,PI
c-----------------------------------------------------------------------
c  mxdata - maximum number of input data points
c  mxisp  - maximum number of points for initial state potential array
c           (also used for number of points in transition moment array)
c  mxnj    - maxiumum value of j quantum number allowed
c  mxfsp  - maximum number of points for final state potential array
c  mxnp   - maximum number of parameters total
c  mxntp  - maximum number of turning points to be read in
c  mxprm  - maximum number of parameters for final state pot'l or TMF
c  mxv    - largest value for the v quantum number
c  mxfs   - maximum number of final states allowed
c  mxisot - maximum number of isotopomers allowed
c  mxsets - maximum number of data sets allowed
c  mxfreq - maximum number of data points allowed in a given set
c-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
      INTEGER FITIT,i,ifs,j,m,MESH1,MESH2,MESH5,NL,step,first,last,IWR,
     1        JP,LPTMF,TURNPT,OTMF,NAMP,NEND,NLORZ,
     2        TMFVAR(0:mxprm-1,mxfs)
      REAL*8 ADD1,ADD2,
     1       AI1(3,mxisp),AI2(3,mxisp),ACCUM1(mxisp),ACCUM2(mxisp),
     2       AMP1,AMP2,AMP3,AMP4,AVAL(3),BFCT,DA(3),DRc(3),
     1       DEN,DER(0:mxprm-1),
     1       dIda1,dIda2,dIda3,dIda4,dIdRc1,dIdRc2,dIdRc3,DI,dSdR,EFN,
     2       ELIM,ER,EDIFF(5),FACT1,FACT2,HALF,HARG,NFACT,OVR,
     3       OVR1,OVR1a,OVR1b,OVR2,OVRCRT,OVRLAP(3),RAD(mxfsp),RC(3),RH,
     4       RMIN,RR(mxfsp),PSI(NEND),PSIc(5),PSIc0,SIXTH,SG1,SG2,SNARG,
     5       SQKINF,THIRD,TMFPRM(0:mxprm-1,mxfs),
     6       TOT,VLIM,VJ(mxfsp),VV,XIIth,XX,Y1,Y2,Y3,W(mxisp),ZTST

      REAL*8 gtk,gtk1,gtk2,gtkfact,AIgtk(3,mxisp),ACCUMgtk(mxisp),
     1       OVRgtk,gtkLAP(3)
c-----------------------------------------------------------------------
c  RAD(i)        - radial distance array
c  NEND -  last mesh point for bound state wavefunction
      NLORZ= (OTMF+1)/2
      i= 0
      DO NL= 1,NLORZ
         AVAL(NL)= TMFPRM(i,ifs)
         RC(NL)=   TMFPRM(i+1,ifs)
         i= 2*NL
      ENDDO
      last= 0
      HALF=  1.D0/2.D0
      THIRD= 1.D0/3.D0
      SIXTH= 1.D0/6.D0
      XIIth= 1.D0/12.D0
      RMIN= RAD(1)
      RH= RAD(2)-RAD(1)
      ER= EFN*BFCT
      ELIM= VLIM*BFCT
      SQKINF= DSQRT(ER-ELIM)
      AMP1= 1.D0
      AMP2= 2.D0
      AMP3= 0.D0
      AMP4= 0.D0
      DO m= 0,OTMF
         DER(m)= 0.D0
      ENDDO
c-----------------------------------------------------------------------
c** Locate first turning point and use Airy function to estimate
c  appropriate integration starting point such that  PSI(1) .LE. 1.D-10
c-----------------------------------------------------------------------
      MESH1= 1
      EDIFF(1)= VJ(MESH1)-ER
      step= DINT(0.2/RH)
      IF(step.LT.1) step= 1
      first= step+1
      DO i= first,mxfsp,step
         TURNPT= i
         EDIFF(2)= VJ(i)-ER
         IF(EDIFF(2).LE. 0.D0) GOTO 4
         MESH1= i
         EDIFF(1)= EDIFF(2)
      ENDDO     
      IF(IWR.NE.0) THEN
         WRITE(6,607) JP,EFN
         OVR= 0.D0
         RETURN
      ENDIF
    4 MESH2= TURNPT
      TURNPT= MESH1+(MESH2-MESH1)*EDIFF(1)/(EDIFF(1)-EDIFF(2))
      IF(IABS(TURNPT-MESH2).LE.1) GOTO 6
      IF((TURNPT.LE.0).OR.(TURNPT.GT.mxfsp)) THEN
         IF(IWR.NE.0) WRITE(6,601) JP,EFN
         STOP
      ENDIF
      MESH1= MESH2
      EDIFF(1)= EDIFF(2)
      EDIFF(2)= VJ(TURNPT)-ER
      GOTO 4
    6 DI= 10.D0/(VJ(TURNPT-1)-VJ(TURNPT))**THIRD
      step= DINT(DI)
      MESH1= MAX0(1,TURNPT-step)
      IF(MESH1.GE.NEND) THEN
         OVR= 0.D0
         RETURN
      ENDIF
    8 EDIFF(1)= VJ(MESH1)-ER
      IF(EDIFF(1).LT.10.D0) GOTO 10
c-----------------------------------------------------------------------
c** Adjust starting point outward to ensure integration scheme stability
c-----------------------------------------------------------------------
      MESH1= MESH1+1
      IF((MESH1-mxfsp).LT.0) GOTO 8
      IF((MESH1-mxfsp).GE.0) THEN
         OVR= 0.D0
         RETURN
      ENDIF
   10 MESH2= MESH1+1
c-----------------------------------------------------------------------
c** WKB starting condition for wave function
c-----------------------------------------------------------------------
      DO 100 NL= 1,NLORZ
         PSIc(1)= 1.D0
         DO i= 1,4     
            EDIFF(i)= VJ(MESH1-1+i)-ER
         ENDDO
         IF((EDIFF(1).GT. 0.D0).AND.(EDIFF(2).GT. 0.D0)) THEN
            SG1= DSQRT(EDIFF(1))
            SG2= DSQRT(EDIFF(2))
            PSIc(2)= PSIc(1)*DSQRT(SG1/SG2)*DEXP((SG1+SG2)/2.D0)
            IF(PSIc(2).LE.PSIc(1)) PSIc(1)= 0.D0
         ELSE
            VV= VJ(MESH2)/BFCT
            IF(IWR.NE.0) WRITE(6,608) JP,EFN,VV,MESH2
            PSIc(1)= 0.D0
            PSIc(2)= 1.D0
         ENDIF
         PSIc0= PSIc(2)
c-----------------------------------------------------------------------
c  notationally speaking, all Yi's refer to values used in the Numerov
c  Algorithm for wavefunction propagation.
c-----------------------------------------------------------------------
         Y1= PSIc(1)*(1.D0-XIIth*EDIFF(1))
         Y2= PSIc(2)*(1.D0-XIIth*EDIFF(2))
         DO i= 3,4
            Y3= Y2 + Y2 - Y1 + EDIFF(i-1)*PSIc(i-1)
            PSIc(i)= Y3/(1.D0-XIIth*EDIFF(i))
            Y1= Y2
            Y2= Y3
         ENDDO
c-----------------------------------------------------------------------
c  Use trapezoid rule for numerical integration.  Initialize OVR
c  values using first section of area.  Note that intensity is
c  proportional to 2*PSIb*dW/dR*PSIc + 4*PSIb*W*dPSIc/dR   
c     W= a/{4a^2 + (R-Rc)^2}
c     dW/da= W/a - 8W^2
c     dW/dR= -2(R-Rc)W^2/a
c  use Lagrangian interpolation scheme for determining derivative at
c  mid-point of 5-equally-spaced wavefunction points...
c     dPSIc/dR= {PSIc(i-2) + 8[PSIc(i+1)-PSIc(i-1)] - PSIc(i+2)}/12*RH
c  assume dPSIc/dR = 0 for first and last mesh points (MESH1 & NEND)
c  assume dPSIc/dR = [PSIc(mesh3)-PSIc(mesh1)]/2RH @ MESH2
c  assume dPSIc/dR = [PSIc(NEND)-PSIc(NEND-2)]/2RH @ NEND-1
c  PSI(i) is the stored bound state wavefunction array
c  PSIc(i) is the stored continuum state wavefunction
c
c  see end of program for multiplicative factors for integration
c-----------------------------------------------------------------------
c  Initialize overlap integrals
c  NOTE that the AIi arrays are cumulative integrands at a given mesh
c-----------------------------------------------------------------------
         DO i= 1,NEND
            RR(i)= RAD(i)-RC(NL)
            DEN= 4*AVAL(NL)*AVAL(NL) + RR(i)*RR(i)
            W(i)= AVAL(NL)/DEN
         ENDDO
         OVR1a= HALF*PSI(MESH1)*RR(MESH1)*W(MESH1)*W(MESH1)*PSIc(1)
         OVR1b=      PSI(MESH2)*RR(MESH2)*W(MESH2)*W(MESH2)*PSIc(2)
         gtk1= HALF*PSI(MESH1)*PSIc(1)
         gtk2= PSI(MESH2)*PSIc(2)
         gtk= gtk1 + gtk2
         OVR1= OVR1a + OVR1b
         OVR2= PSI(MESH2)*W(MESH2)*6.d0*(PSIc(3)-PSIc(1))
c3pt     OVR2= PSI(MESH2)*W(MESH2)*(PSIc(3)-PSIc(1))
         AI1(NL,MESH1)= OVR1a
         AI1(NL,MESH2)= OVR1b
         AI2(NL,MESH1)= 0.d0
         AI2(NL,MESH2)= OVR2
         AIgtk(NL,MESH1)= gtk1
         AIgtk(NL,MESH2)= gtk2
c-----------------------------------------------------------------------
c  for portion of operator involving the derivative of final-state wave
c  function, assume dPSI/dR @ MESH1 = 0
c  use dPSI/dR @ MESH2 = (PSI3-PSI1)/2RH... below as 6*(PSI3-PSI1)
c  since entire sum is (later) divided by 12
c-----------------------------------------------------------------------
         IF(TMFVAR(2*NL-2,ifs).GT.0) THEN
            dIda2= OVR1a*W(MESH1) + OVR1b*W(MESH2)
            dIda4= OVR2*W(MESH2)
         ENDIF
         IF(TMFVAR(2*NL-1,ifs).GT.0) THEN
            dIdRc1= OVR1a/RR(MESH1) + OVR1b/RR(MESH2)
            dIdRc2= OVR1a*RR(MESH1)*W(MESH1) + OVR1b*RR(MESH2)*W(MESH2)
            dIdRc3= OVR2*RR(MESH2)*W(MESH2)
         ENDIF
c-----------------------------------------------------------------------
c  Note that there is no contribution at mesh1 for dI/da(4) or dI/dRc(4)
c  since the derivative of the final state wave function is zero here
c-----------------------------------------------------------------------
c** Integrate outward to first turning point.  NOTE that Airy-estimated
c  initialization minimizes need for renormalizations.
c-----------------------------------------------------------------------
         MESH5= MESH2 + 3
         Y3= Y2 + Y2 - Y1 + EDIFF(4)*PSIc(4)
         DO 16 i= MESH5,TURNPT
            EDIFF(5)= VJ(i)-ER
            PSIc(5)= Y3/(1.D0-XIIth*EDIFF(5))
            Y3= Y2 + Y2 - Y1 + EDIFF(5)*PSIc(5)
            Y1= Y2
            Y2= Y3
c----------------------------------------------------------------------
c NOW, If bound wavefx. non-negligible, accumulate overlap integrals
c **Do ALL integration calculations 2 steps behind the propagating front
c   of the wavefunction (since MUST do this for dPSIc/dR)
c-----------------------------------------------------------------------
            IF(i.LE.NEND) THEN
               ADD1= PSI(i-2)*RR(i-2)*W(i-2)*W(i-2)*PSIc(3)
               dSdR= PSIc(1) + 8.d0*(PSIc(4) - PSIc(2)) - PSIc(5)
c3pt           dSdR= PSIc(4) - PSIc(2)
               ADD2= PSI(i-2)*W(i-2)*dSdR
               OVR1= OVR1 + ADD1
               OVR2= OVR2 + ADD2
               gtk= gtk + PSI(i-2)*PSIc(3)
               AI1(NL,i-2)= OVR1
               AI2(NL,i-2)= OVR2
               AIgtk(NL,i-2)= gtk
               IF(TMFVAR(2*NL-2,ifs).GT.0) THEN
                  dIda2= dIda2 + ADD1*W(i-2)
                  dIda4= dIda4 + ADD2*W(i-2)
               ENDIF
               IF(TMFVAR(2*NL-1,ifs).GT.0) THEN
                  dIdRc1= dIdRc1 + ADD1/RR(i-2)
                  dIdRc2= dIdRc2 + ADD1*RR(i-2)*W(i-2)
                  dIdRc3= dIdRc3 + ADD2*RR(i-2)*W(i-2)
               ENDIF
            ENDIF
c-----------------------------------------------------------------------
c  If wavefuntion too large in forbidden region, renormalize it ...
c  Also renormalize overlap integrals etc. here
c-----------------------------------------------------------------------
            IF((PSIc(5).GE.1.D32).OR.(i.EQ.TURNPT)) THEN
               NFACT= 1.D0/PSIc(5)
               PSIc(5)= 1.D0
               DO j= 1,4
                  PSIc(j)= PSIc(j)*NFACT
               ENDDO
               IF(PSIc0.GT.1.D-30) PSIc0= PSIc0*NFACT
               OVR1= OVR1*NFACT
               OVR2= OVR2*NFACT
               gtk= gtk*NFACT
               DO j= MESH1,i-2
                  AI1(NL,j)= AI1(NL,j)*NFACT
                  AI2(NL,j)= AI2(NL,j)*NFACT
                  AIgtk(NL,j)= AIgtk(NL,j)*NFACT
               ENDDO
               Y1= Y1*NFACT
               Y2= Y2*NFACT
               Y3= Y3*NFACT
               IF(TMFVAR(2*NL-2,ifs).GT.0) THEN
                  dIda2= dIda2*NFACT
                  dIda4= dIda4*NFACT         
               ENDIF
               IF(TMFVAR(2*NL-1,ifs).GT.0) THEN
                  dIdRc1= dIdRc1*NFACT       
                  dIdRc2= dIdRc2*NFACT             
                  dIdRc3= dIdRc3*NFACT               
               ENDIF
            ENDIF
            DO j= 1,4
               PSIc(j)= PSIc(j+1)
            ENDDO
   16    CONTINUE
         IF((IWR.NE.0).AND.(PSIc0/PSIc(5).GT.1.D-8))
     1      WRITE(6,602)JP,EFN,MESH1,PSIc0/PSIc(5),TURNPT
c-----------------------------------------------------------------------
c If turning point NOT past end of range for bound state wavefx., then
c integrate from turning point to end of bound-state wave function
c  NOTE that summations will cease at mesh NEND-2 because of the 2-step
c  lag in the calculation.  This should have essentially no effect on
c  integration total
c-----------------------------------------------------------------------
         IF(TURNPT.LT.NEND) THEN
            DO i= TURNPT+1,NEND
               EDIFF(5)= VJ(i)-ER
               PSIc(5)= Y3/(1.D0-XIIth*EDIFF(5))
               Y3= Y2 + Y2 - Y1 + EDIFF(5)*PSIc(5)
               Y1= Y2
               Y2= Y3
c-----------------------------------------------------------------------
c  NOTE that the derivative term is lagging two steps behind the counter
c-----------------------------------------------------------------------
               ADD1= PSI(i-2)*RR(i-2)*W(i-2)*W(i-2)*PSIc(3)
               dSdR= PSIc(1) + 8.d0*(PSIc(4) - PSIc(2)) - PSIc(5)
c3pt           dSdR= PSIc(4) - PSIc(2)
               ADD2= PSI(i-2)*W(i-2)*dSdR
               OVR1= OVR1 + ADD1
               OVR2= OVR2 + ADD2
               gtk= gtk + PSI(i-2)*PSIc(3)
               AI1(NL,i-2)= OVR1
               AI2(NL,i-2)= OVR2
               AIgtk(NL,i-2)= gtk
               IF(TMFVAR(2*NL-2,ifs).GT.0) THEN
                  dIda2= dIda2 + ADD1*W(i-2)
                  dIda4= dIda4 + ADD2*W(i-2)
               ENDIF
               IF(TMFVAR(2*NL-1,ifs).GT.0) THEN
                  dIdRc1= dIdRc1 + ADD1/RR(i-2)
                  dIdRc2= dIdRc2 + ADD1*RR(i-2)*W(i-2)
                  dIdRc3= dIdRc3 + ADD2*RR(i-2)*W(i-2)
               ENDIF
               DO j= 1,4
                  PSIc(j)= PSIc(j+1)
               ENDDO
            ENDDO   
         ENDIF
c-----------------------------------------------------------------------
c Continue wave function propagation until amplitude converges, no
c longer integrating since past the end of bound state range
c-----------------------------------------------------------------------
         NAMP= 0
         DO i= NEND+1,mxfsp
            EDIFF(4)= EDIFF(5)
            EDIFF(5)= VJ(i)-ER
            PSIc(5)= Y3/(1.D0-XIIth*EDIFF(5))
            Y3= Y2 + Y2 - Y1 + EDIFF(5)*PSIc(5)
            Y1= Y2
            Y2= Y3
            IF((PSIc(5).LT.PSIc(4)).AND.(PSIc(3).LT.PSIc(4))) THEN
c-----------------------------------------------------------------------
c At successive maxima, fit solution to W.K.B. form to determine
c apparent asymptotic amplitude.
c-----------------------------------------------------------------------
               SG1= DSQRT(-EDIFF(4))
               SG2= DSQRT(-EDIFF(5))
               HARG= HALF*(SG1+SG2)
               SNARG= 1.D0/DSQRT(1.D0+((DSQRT(SG1/SG2)*PSIc(4)/PSIc(5) -
     1                DCOS(HARG))/DSIN(HARG))**2)
               NAMP= NAMP+1
               AMP4= AMP3
               AMP3= AMP2
               AMP2= AMP1
               AMP1= PSIc(5)*DSQRT(SG2/SQKINF)/SNARG
               XX= RMIN + (i-1)*RH
               IF(IWR.GT.1) WRITE(6,604) JP,EFN,XX,AMP1
               last= i
c-----------------------------------------------------------------------
c Test successive amplitudes for convergence
c-----------------------------------------------------------------------
               ZTST= OVRCRT*AMP1
               IF((DABS(AMP1-AMP2).LT.ZTST).AND.
     1           (DABS(AMP2-AMP3).LT.ZTST)) GOTO 35
            ENDIF
            DO j= 1,4
               PSIc(j)= PSIc(j+1)
            ENDDO
         ENDDO
         IF(IWR.NE.0) WRITE(6,603) JP,EFN,AMP1,AMP2,AMP3,AMP4
   35    FACT1= 2.D0*RH*RH*RH/(AVAL(NL)*BFCT)
         FACT2= -SIXTH*RH*RH/BFCT
c3pt     FACT2= -*RH*RH/BFCT
         OVR1= OVR1*FACT1
         OVR2= OVR2*FACT2
         gtkfact= 1.d0
         gtk= gtk*gtkfact
         DO i= 1,NEND
            AI1(NL,i)= AI1(NL,i)*FACT1/AMP1
            AI2(NL,i)= AI2(NL,i)*FACT2/AMP1
            AIgtk(NL,i)= AIgtk(NL,i)*gtkfact/AMP1
         ENDDO
         OVRLAP(NL)= (OVR1 + OVR2)/AMP1
         gtkLAP(NL)= gtk/AMP1
         IF(TMFVAR(2*NL-2,ifs).GT.0) THEN
            dIda1= OVR1/AVAL(NL)
            dIda2= -16.d0*FACT1*dIda2
            dIda3= OVR2/AVAL(NL)
            dIda4= -8.d0*FACT2*dIda4
            DA(NL)= dIda1 + dIda2 + dIda3 + dIda4
         ENDIF
         IF(TMFVAR(2*NL-1,ifs).GT.0) THEN
            dIdRc1= -FACT1*dIdRc1
            dIdRc2= 4.d0*FACT1*dIdRc2/AVAL(NL)
            dIdRc3= 2.d0*FACT2*dIdRc3/AVAL(NL)
            DRc(NL)= dIdRc1 + dIdRc2 + dIdRc3
         ENDIF
c-----------------------------------------------------------------------
c  note that total intensity is the square of the overlap integral
c  and for total derivative need 2*(overlap integral)*(dI/dp).
c-----------------------------------------------------------------------
  100 CONTINUE
      OVR= 0.d0
      OVRgtk= 0.d0
      DO i= 1,NEND
         ACCUM1(i)= 0.d0
         ACCUM2(i)= 0.d0
         ACCUMgtk(i)= 0.d0
      ENDDO
      DO NL= 1,NLORZ
         OVR= OVR + OVRLAP(NL)
         OVRgtk= OVRgtk + gtkLAP(NL)
         DO i= 1,NEND-2
            ACCUM1(i)= ACCUM1(i) + AI1(NL,i)
            ACCUM2(i)= ACCUM2(i) + AI2(NL,i)
            ACCUMgtk(i)= ACCUMgtk(i) + AIgtk(NL,i)
         ENDDO
      ENDDO
      IF(IWR.GE.1) WRITE(6,605) JP,EFN,last,XX,OVR
      IF(IWR.GE.2) WRITE(6,606) PSIc0,AMP1,AMP2,AMP3,AMP4
      NL= 1
      DO m= 0,OTMF-1,2
         DER(m)= 2.d0*OVR*DA(NL)
         DER(m+1)= 2.d0*OVR*DRc(NL)
         NL= NL+1
      ENDDO
      OVR= OVR*OVR
      OVRgtk= OVRgtk*OVRgtk
      IF(LPTMF.GT.0) THEN
         WRITE(10,1000)
         DO i= 1,NEND-2,LPTMF
            TOT= ACCUM1(i) + ACCUM2(i)   
            WRITE(10,1001) RAD(i),ACCUM1(i),ACCUM2(i),TOT,ACCUMgtk(i)
         ENDDO
      ENDIF
      RETURN
c-----------------------------------------------------------------------
  601 FORMAT(' *** OVRLAP BOMBed *** '/ 'For  J = ',I3,'   EFN = ',
     1  F10.2,'   never got to first turning point')
  602 FORMAT(' ** WARNING **  For   J = ', I3 ,'   EFN = ',F10.2,
     1 '  starting wavefunction is PSI(',I4,') = ',D10.3,' and  I(turn.p
     2t.) = ',I4)
  603 FORMAT(' ** WARNING ** For J= ',I3,' EFN= ',F9.2,' amplitude not',
     1 ' converged by end of range'/3x,'Last four values are ',
     2 4(1PD14.6))
  604 FORMAT(' At J =',I3,'  EFN = ',F10.2,'  R = ',F7.4,
     1 ', apparent asymptotic amplitude is',G15.8)
  605 FORMAT(/' At JP=',I3,' EFN= ',F9.2,'  Range to  R(',I4,')=',
     1 F7.4/' OVRPD gives overlap integral= ',F12.8)
  606 FORMAT(5x,'where  PSIc(initial)= ',D10.3,'  & last 4 amplitudes',
     1 ' are',1P4D15.7)
  607 FORMAT(' *** ERROR ***   At   J = ',I3,'   EFN = ',F10.2,
     1  '  have  V .GT. E  everywhere.' )
  608 FORMAT(' *** Caution ***  For   J = ',I3,'   (EFN= ',F10.2,
     1  ') .GE. (V = ',F10.2,')   at  I = ',I4,' ,so initialize with a n
     2ode.')
cgtk0 FORMAT(' Accumulated Integral Lorentzian type TMF'//
cgtk 1 '   RAD   dW/dR term    dPSIc/dR term    TOTAL')
 1000 FORMAT(' Accumulated Integral Lorentzian type TMF'//
     1 '   RAD   dW/dR term    dPSIc/dR term    TOTAL       PSIb*PSIc')
cgtk0 FORMAT((F7.3,3(1PD14.6)))
 1001 FORMAT((F7.3,4(1PD14.6)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE HONL(HLFACT,JPP,JP,OMEGA,PQR,ifs)
c=======================================================================
c  Routine to calculate Honl-London factors (HLFACT)
c  Last modified 9 April 2007 to take account of Hansson & Watson(2005)
c=======================================================================
ccc   INCLUDE 'arrsizes.h'
c-----------------------------------------------------------------------
c  Utility routine to summarize dimensioning of arrays
c-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,
     1        mxsets,mxfreq
      REAL*8 CCM,PI
c-----------------------------------------------------------------------
c  mxdata - maximum number of input data points
c  mxisp  - maximum number of points for initial state potential array
c           (also used for number of points in transition moment array)
c  mxnj    - maxiumum value of j quantum number allowed
c  mxfsp  - maximum number of points for final state potential array
c  mxnp   - maximum number of parameters total
c  mxntp  - maximum number of turning points to be read in
c  mxprm  - maximum number of parameters for final state pot'l or TMF
c  mxv    - largest value for the v quantum number
c  mxfs   - maximum number of final states allowed
c  mxisot - maximum number of isotopomers allowed
c  mxsets - maximum number of data sets allowed
c  mxfreq - maximum number of data points allowed in a given set
c-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
c=======================================================================
      INTEGER DOMEGA,ifs,JPP,JP,OMEGA(0:mxfs),PQR
      REAL*8  HLFACT,J,L
c
      HLFACT= 1.d0
      IF(PQR.LE.0) GOTO 100
      HLFACT= 0.d0
      IF(JP.LT.0) GOTO 100
      IF((JP.EQ.JPP).AND.(JP.EQ.0)) GOTO 100
      J= DBLE(JPP)
      L= DBLE(OMEGA(0))
      DOMEGA= OMEGA(ifs) - OMEGA(0)
      IF(DOMEGA.EQ.0) THEN
         IF(JP.EQ.JPP+1) HLFACT= (J+1.d0+L)*(J+1.d0-L)/(J+1.d0)
         IF(JP.EQ.JPP)   HLFACT= (2.d0*J+1.d0)*L*L/(J*(J+1.d0))
         IF(JP.EQ.JPP-1) HLFACT= (J+L)*(J-L)/J
      ELSEIF(DOMEGA.EQ.1) THEN
         IF(JP.EQ.JPP+1) HLFACT= (J+2.d0+L)*(J+1.d0+L)/(2.d0*(J+1.d0))
         IF(JP.EQ.JPP)   HLFACT= (J+1.d0+L)*(J-L)*(2.d0*J+1.d0)/
     1                           (2.d0*J*(J+1.d0))
         IF(JP.EQ.JPP-1) HLFACT= (J-1.d0-L)*(J-L)/(2.d0*J)
         IF(MIN(OMEGA(0),OMEGA(ifs)).EQ.0) HLFACT= HLFACT + HLFACT
      ELSEIF(DOMEGA.EQ.-1) THEN
         IF(JP.EQ.JPP+1) HLFACT= (J+2.d0-L)*(J+1.d0-L)/(2.d0*(J+1.d0))
         IF(JP.EQ.JPP)   HLFACT= (J+1.d0-L)*(J+L)*(2.d0*J+1.d0)/
     1                           (2.d0*J*(J+1.d0))
         IF(JP.EQ.JPP-1) HLFACT= (J-1.d0+L)*(J+L)/(2.d0*J)
         IF(MIN(OMEGA(0),OMEGA(ifs)).EQ.0) HLFACT= HLFACT + HLFACT
      ENDIF
c-----------------------------------------------------------------------
c  output HLFACT divided by (2J+1) of its actual value in order to keep
c  proper population cutoff intact (see POPF in forward.f)
c-----------------------------------------------------------------------
      HLFACT= HLFACT/(2.d0*J+1.d0)
  100 RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12


c***********************************************************************
      SUBROUTINE JAVGE(KV,NJ,JM,BV,TCM)
c=======================================================================
c  This subroutine calculates JM(M), the average J for each of the NJ
c  equally weighted  segments of the rotational population for
c  vibrational level whose rotational constant is  Bv, all at
c  temperature  TCM (in cm-1). 
c
c    JM(NJ,M)=-0.5+FJ1(M)*(8*DSQRT(K*T/Bv)+DSQRT(Bv/K*T))-FJ2(M)* ...
c 
c  Taken from Eq.(21) of  Le Roy et al. (J.Chem.Phys. 65, 1485 (1976))
c-----------------  Last updated  23 February 2001 ---------------------
c=======================================================================
      INTEGER   I,IFIRST,JM(20),KV,M,NJ,NJ1
      REAL*8    A1,A2,BV,DERF,E1,E2,F1,F2,FJ1(20),FJ2(20),TCM,TWR,X,ZN   
      DATA IFIRST/-1/
      SAVE FJ2,FJ1,IFIRST
      DERF(X) = 1.D0 - (1.D0 + X*(.278393D0 + X*(.230389D0 + X*(.972D-3
     1          + X*.078108D0))))**(-4)
c=======================================================================
      IF (NJ.NE.IFIRST) THEN
          IF (NJ.LE.0) THEN
              FJ1(1) = 0.d0
              FJ2(1) = 0.d0
              IFIRST= NJ
cc            WRITE(6,600)
              GOTO 100
              ENDIF
          ZN = NJ
          F1 = ZN*0.2215567314D0
          A1 = 0.D0
          E1 = 0.D0
          IF (NJ.GT.1) THEN
              NJ1 = NJ - 1
              DO  M=1,NJ1
                  A2 = A1
                  A1 = DSQRT(DLOG(ZN/DBLE(NJ-M)))
                  E2 = E1
                  E1 = DERF(A1)
                  FJ2(M) = -(NJ-M)*A1+(NJ-M+1)*A2
                  FJ1(M) = F1*(E1-E2)
                  ENDDO
              ENDIF
          A2 = A1
          E2 = E1
          FJ2(NJ) = A2
          FJ1(NJ) = F1*(1.D0-E2)
cc        WRITE(6,601) NJ,(I,FJ1(I),I,FJ2(I),I=1,NJ)
          IFIRST= NJ
          ENDIF
  100 IF (TCM.LE. 0.D0) NJ = 0
      IF (NJ.GT.0) THEN
          F2 = DSQRT(TCM/BV)
          F1 = 4.D0*F2+1.D0/F2
          DO M=1,NJ
              JM(M)=F1*FJ1(M)+F2*FJ2(M)
              ENDDO       
          TWR = TCM/.6950387d0
cc        WRITE(6,603) KV,BV,TWR,NJ,(JM(I),I=1,NJ)
        ELSE
          JM(1) = 0
        ENDIF
      RETURN
c-----------------------------------------------------------------------
cc600 FORMAT(/' Fix J=0 rather than sum over a rotational distribution')
cc601 FORMAT(/' Divide rotational distbn for each vib. level into',I3,
cc   1 ' equally weighted segments.'/
cc   2 (2(4x:'FJ1(',I2,')=',F8.5,'   FJ2(',I2,')=',F8.5)))
cc603 FORMAT(/' For  v =',I2,' with  Bv =',F9.6,' at  T =',
cc   1  F7.1,'  the',I3,' equally weighted J-s are'/(16I5))
c-----------------------------------------------------------------------
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MASSES(IAN,IMN,NAME,GELGS,GNS,MASS,ABUND)
c***********************************************************************
c** For isotope with (input) atomic number IAN and mass number IMN,
c  return (output):  (i) as the right-adjusted 2-character variable NAME
c  the alphabetic symbol for that element,  (ii) the ground state
c  electronic degeneracy GELGS, (iii) the nuclear spin degeneracy GNS,
c  (iv) the atomic mass MASS [amu], and  (v) the natural isotopic
c  abundance ABUND [in percent].   GELGS values based on atomic states
c  in Moore's "Atomic Energy Level" tables, the isotope masses are taken
c  from the 2003 mass table [Audi, Wapstra & Thibault, Nucl.Phys. A729,
c  337-676 (2003)] and other quantities from Tables 6.2 and 6.3 of
c  "Quantities, Units and Symbols in Physical Chemistry", by Mills et
c  al. (Blackwell, 2'nd Edition, Oxford, 1993).
c** If the input value of IMN does not equal one of the tabulated values
c  for atomic species IAN, return the abundance-averaged standard atomic
c  weight of that atom and set GNS=-1 and ABUND=-1.
c                          COPYRIGHT 2005
c** By R.J. Le Roy (with assistance from G.T. Kraemer & J.Y. Seto).
c                      Last modified  1 June 2005
c***********************************************************************
      REAL*8 zm(123,0:10),mass,ab(123,10),abund
      INTEGER i,ian,imn,gel(123),nmn(123),mn(123,10),ns2(123,10),
     1        gelgs, gns
      CHARACTER*2 NAME,AT(123)
c
      DATA  at(1),gel(1),nmn(1),(mn(1,i),i=1,3)/' H',2,3,1,2,3/
      DATA  (zm(1,i),i=0,3)/1.00794d0, 1.00782503207d0, 2.0141017778d0,
     1                      3.0160492777d0/
      DATA  (ns2(1,i),i=1,3)/1,2,1/
      DATA  (ab(1,i),i=1,3)/99.985d0,0.015d0,0.d0/
c
      DATA  at(2),gel(2),nmn(2),(mn(2,i),i=1,2)/'He',1,2,3,4/
      DATA  (zm(2,i),i=0,2)/4.002602d0, 3.0160293191d0, 4.00260325415d0/
      DATA  (ns2(2,i),i=1,2)/1,0/
      DATA  (ab(2,i),i=1,2)/0.000137d0,99.999863d0/
c
      DATA  at(3),gel(3),nmn(3),(mn(3,i),i=1,2)/'Li',2,2,6,7/
      DATA  (zm(3,i),i=0,2)/6.941d0, 6.015122795d0, 7.01600455d0/
      DATA  (ns2(3,i),i=1,2)/2,3/
      DATA  (ab(3,i),i=1,2)/7.5d0,92.5d0/
c
      DATA  at(4),gel(4),nmn(4),(mn(4,i),i=1,1)/'Be',1,1,9/
      DATA  (zm(4,i),i=0,1)/9.012182d0, 9.0121822d0/
      DATA  (ns2(4,i),i=1,1)/3/
      DATA  (ab(4,i),i=1,1)/100.d0/
c
      DATA at(5),gel(5),nmn(5),(mn(5,i),i=1,2)/' B',2,2,10,11/
      DATA (zm(5,i),i=0,2)/10.811d0, 10.0129370d0, 11.0093054d0/
      DATA  (ns2(5,i),i=1,2)/6,3/
      DATA  (ab(5,i),i=1,2)/19.9d0,80.1d0/
c
      DATA at(6),gel(6),nmn(6),(mn(6,i),i=1,3)/' C',1,3,12,13,14/
      DATA (zm(6,i),i=0,3)/12.011d0, 12.d0, 13.0033548378d0,
     1                      14.003241989d0/
      DATA  (ns2(6,i),i=1,3)/0,1,0/
      DATA  (ab(6,i),i=1,3)/98.90d0,1.10d0, 0.d0/
c
      DATA at(7),gel(7),nmn(7),(mn(7,i),i=1,2)/' N',4,2,14,15/
      DATA (zm(7,i),i=0,2)/14.00674d0, 14.0030740048d0, 15.0001088982d0/
      DATA (ns2(7,i),i=1,2)/2,1/
      DATA (ab(7,i),i=1,2)/99.634d0,0.366d0/
c
      DATA at(8),gel(8),nmn(8),(mn(8,i),i=1,3)/' O',5,3,16,17,18/
      DATA (zm(8,i),i=0,3)/15.9994d0, 15.99491461956d0, 16.99913170d0,
     1                      17.9991610d0/
      DATA (ns2(8,i),i=1,3)/0,5,0/
      DATA (ab(8,i),i=1,3)/99.762d0, 0.038d0, 0.200d0/
c
      DATA at(9),gel(9),nmn(9),(mn(9,i),i=1,1)/' F',4,1,19/
      DATA (zm(9,i),i=0,1)/18.9984032d0, 18.99840322d0/
      DATA (ns2(9,i),i=1,1)/1/
      DATA (ab(9,i),i=1,1)/100.d0/
c
      DATA at(10),gel(10),nmn(10),(mn(10,i),i=1,3)/'Ne',1,3,20,21,22/
      DATA (zm(10,i),i=0,3)/20.1797d0, 19.9924401754d0, 20.99384668d0,
     1                       21.991385114d0/
      DATA (ns2(10,i),i=1,3)/0,3,0/
      DATA (ab(10,i),i=1,3)/90.48d0, 0.27d0, 9.25d0/
c
      DATA at(11),gel(11),nmn(11),(mn(11,i),i=1,1)/'Na',2,1,23/
      DATA (zm(11,i),i=0,1)/22.989768d0, 22.9897692809d0/
      DATA (ns2(11,i),i=1,1)/3/
      DATA (ab(11,i),i=1,1)/100.d0/
c
      DATA at(12),gel(12),nmn(12),(mn(12,i),i=1,3)/'Mg',1,3,24,25,26/
      DATA (zm(12,i),i=0,3)/24.3050d0, 23.985041700d0, 24.98583692d0,
     1                       25.982592929d0/
      DATA (ns2(12,i),i=1,3)/0,5,0/
      DATA (ab(12,i),i=1,3)/78.99d0, 10.00d0, 11.01d0/
c
      DATA at(13),gel(13),nmn(13),(mn(13,i),i=1,1)/'Al',2,1,27/
      DATA (zm(13,i),i=0,1)/26.981539d0, 26.98153863d0/
      DATA (ns2(13,i),i=1,1)/5/
      DATA (ab(13,i),i=1,1)/100.d0/
c
      DATA at(14),gel(14),nmn(14),(mn(14,i),i=1,3)/'Si',1,3,28,29,30/
      DATA (zm(14,i),i=0,3)/28.0855d0, 27.9769265325d0, 28.976494700d0,
     1                       29.97377017d0/
      DATA (ns2(14,i),i=1,3)/0,1,0/
      DATA (ab(14,i),i=1,3)/92.23d0, 4.67d0, 3.10d0/

      DATA at(15),gel(15),nmn(15),(mn(15,i),i=1,1)/' P',4,1,31/
      DATA (zm(15,i),i=0,1)/30.973762d0, 30.97376163d0/
      DATA (ns2(15,i),i=1,1)/1/
      DATA (ab(15,i),i=1,1)/100.d0/
c
      DATA at(16),gel(16),nmn(16),(mn(16,i),i=1,4)/' S',5,4,32,33,34,36/
      DATA (zm(16,i),i=0,4)/32.066d0, 31.97207100d0, 32.97145876d0,
     1                       33.96786690d0, 35.96708076d0/
      DATA (ns2(16,i),i=1,4)/0,3,0,0/
      DATA (ab(16,i),i=1,4)/95.02d0, 0.75d0, 4.21d0, 0.02d0/
c
      DATA at(17),gel(17),nmn(17),(mn(17,i),i=1,2)/'Cl',4,2,35,37/
      DATA (zm(17,i),i=0,2)/35.4527d0, 34.96885268d0, 36.96590259d0/
      DATA (ns2(17,i),i=1,2)/3,3/
      DATA (ab(17,i),i=1,2)/75.77d0, 24.23d0/
c
      DATA at(18),gel(18),nmn(18),(mn(18,i),i=1,3)/'Ar',1,3,36,38,40/
      DATA (zm(18,i),i=0,3)/39.948d0, 35.967545106d0, 37.9627324d0,
     1                       39.9623831225d0/
      DATA (ns2(18,i),i=1,3)/0,0,0/
      DATA (ab(18,i),i=1,3)/0.337d0, 0.063d0, 99.600d0/
c
      DATA at(19),gel(19),nmn(19),(mn(19,i),i=1,3)/' K',2,3,39,40,41/
      DATA (zm(19,i),i=0,3)/39.0983d0, 38.96370668d0, 39.96399848d0,
     1                       40.96182576d0/
      DATA (ns2(19,i),i=1,3)/3,8,3/
      DATA (ab(19,i),i=1,3)/93.2581d0, 0.0117d0, 6.7302d0/

      DATA at(20),gel(20),nmn(20),(mn(20,i),i=1,6)/'Ca',1,6,40,42,43,44,
     1                                              46,48/
      DATA (zm(20,i),i=0,6)/40.078d0, 39.96259098d0, 41.95861801d0,
     1         42.9587666d0, 43.9554818d0, 45.9536926d0, 47.952534d0/
      DATA (ns2(20,i),i=1,6)/0,0,7,0,0,0/
      DATA (ab(20,i),i=1,6)/96.941d0, 0.647d0, 0.135d0, 2.086d0,
     1                      0.004d0, 0.187d0/
c
      DATA at(21),gel(21),nmn(21),(mn(21,i),i=1,1)/'Sc',4,1,45/
      DATA (zm(21,i),i=0,1)/44.955910d0, 44.9559119d0/
      DATA (ns2(21,i),i=1,1)/7/
      DATA (ab(21,i),i=1,1)/100.d0/
c
      DATA at(22),gel(22),nmn(22),(mn(22,i),i=1,5)/'Ti',5,5,46,47,48,49,
     1                                              50/
      DATA (zm(22,i),i=0,5)/47.88d0, 45.9526316d0, 46.9517631d0,
     1         47.9479463d0, 48.9478700d0, 49.9447912d0/
      DATA (ns2(22,i),i=1,5)/0,5,0,7,0/
      DATA (ab(22,i),i=1,5)/8.0d0, 7.3d0, 73.8d0, 5.5d0, 5.4d0/
c
      DATA at(23),gel(23),nmn(23),(mn(23,i),i=1,2)/' V',4,2,50,51/
      DATA (zm(23,i),i=0,2)/50.9415d0, 49.9471585d0, 50.9439595d0/
      DATA (ns2(23,i),i=1,2)/12,7/
      DATA (ab(23,i),i=1,2)/0.250d0, 99.750d0/
c
      DATA at(24),gel(24),nmn(24),(mn(24,i),i=1,4)/'Cr',7,4,50,52,53,54/
      DATA (zm(24,i),i=0,4)/51.9961d0, 49.9460442d0, 51.9405075d0,
     1                       52.9406494d0, 53.9388804d0/
      DATA (ns2(24,i),i=1,4)/0,0,3,0/
      DATA (ab(24,i),i=1,4)/4.345d0, 83.789d0, 9.501d0, 2.365d0/
c
      DATA at(25),gel(25),nmn(25),(mn(25,i),i=1,1)/'Mn',6,1,55/
      DATA (zm(25,i),i=0,1)/54.93805d0, 54.9380451d0/
      DATA (ns2(25,i),i=1,1)/5/
      DATA (ab(25,i),i=1,1)/100.d0/
c
      DATA at(26),gel(26),nmn(26),(mn(26,i),i=1,4)/'Fe',9,4,54,56,57,58/
      DATA (zm(26,i),i=0,4)/55.847d0, 53.9396105d0, 55.9349375d0,
     1                       56.9353940d0, 57.9332756d0/
      DATA (ns2(26,i),i=1,4)/0,0,1,0/
      DATA (ab(26,i),i=1,4)/5.8d0, 91.72d0, 2.2d0, 0.28d0/
c
      DATA at(27),gel(27),nmn(27),(mn(27,i),i=1,1)/'Co',10,1,59/
      DATA (zm(27,i),i=0,1)/58.93320d0, 58.9331950d0/
      DATA (ns2(27,i),i=1,1)/7/
      DATA (ab(27,i),i=1,1)/100.d0/
c
      DATA at(28),gel(28),nmn(28),(mn(28,i),i=1,5)/'Ni',9,5,58,60,61,62,
     1                                              64/
      DATA (zm(28,i),i=0,5)/58.69d0, 57.9353429d0, 59.9307864d0,
     1         60.9310560d0, 61.9283451d0, 63.9279660d0/
      DATA (ns2(28,i),i=1,5)/0,0,3,0,0/
      DATA (ab(28,i),i=1,5)/68.077d0,26.223d0,1.140d0,3.634d0,0.926d0/
c
      DATA at(29),gel(29),nmn(29),(mn(29,i),i=1,2)/'Cu',2,2,63,65/
      DATA (zm(29,i),i=0,2)/63.546d0, 62.9295975d0,64.9277895d0/
      DATA (ns2(29,i),i=1,2)/3,3/
      DATA (ab(29,i),i=1,2)/69.17d0, 30.83d0/
c
      DATA at(30),gel(30),nmn(30),(mn(30,i),i=1,5)/'Zn',1,5,64,66,67,68,
     1                                              70/
      DATA (zm(30,i),i=0,5)/65.40d0, 63.9291422d0, 65.9260334d0,
     1         66.9271273d0, 67.9248442d0, 69.9253193d0/
      DATA (ns2(30,i),i=1,5)/0,0,5,0,0/
      DATA (ab(30,i),i=1,5)/48.6d0, 27.9d0, 4.1d0, 18.8d0, 0.6d0/
c
      DATA at(31),gel(31),nmn(31),(mn(31,i),i=1,2)/'Ga',2,2,69,71/
      DATA (zm(31,i),i=0,2)/69.723d0, 68.9255736d0, 70.9247013d0/
      DATA (ns2(31,i),i=1,2)/3,3/
      DATA (ab(31,i),i=1,2)/60.108d0, 39.892d0/
c
      DATA at(32),gel(32),nmn(32),(mn(32,i),i=1,5)/'Ge',1,5,70,72,73,74,
     1                                              76/
      DATA (zm(32,i),i=0,5)/72.61d0, 69.9242474d0, 71.9220758d0,
     1         72.9234589d0, 73.9211778d0, 75.9214026d0/
      DATA (ns2(32,i),i=1,5)/0,0,9,0,0/
      DATA (ab(32,i),i=1,5)/21.23d0, 27.66d0, 7.73d0, 35.94d0, 7.44d0/
c
      DATA at(33),gel(33),nmn(33),(mn(33,i),i=1,1)/'As',4,1,75/
      DATA (zm(33,i),i=0,1)/74.92159d0, 74.9215965d0/
      DATA (ns2(33,i),i=1,1)/3/
      DATA (ab(33,i),i=1,1)/100.d0/
c
      DATA at(34),gel(34),nmn(34),(mn(34,i),i=1,6)/'Se',5,6,74,76,77,78,
     1                                              80,82/
      DATA (zm(34,i),i=0,6)/78.96d0, 73.9224764d0, 75.9192136d0,
     1         76.9199140d0, 77.9173091d0, 79.9165213d0, 81.9166994d0/
      DATA (ns2(34,i),i=1,6)/0,0,1,0,0,0/
      DATA (ab(34,i),i=1,6)/0.89d0, 9.36d0, 7.63d0, 23.78d0, 49.61d0,
     1                      8.73d0/
c
      DATA at(35),gel(35),nmn(35),(mn(35,i),i=1,2)/'Br',4,2,79,81/
      DATA (zm(35,i),i=0,2)/79.904d0, 78.9183371d0, 80.9162906d0/
      DATA (ns2(35,i),i=1,2)/3,3/
      DATA (ab(35,i),i=1,2)/50.69d0, 49.31d0/
c
      DATA at(36),gel(36),nmn(36),(mn(36,i),i=1,6)/'Kr',1,6,78,80,82,83,
     1                                              84,86/
      DATA (zm(36,i),i=0,6)/83.80d0, 77.9203648d0, 79.9163790d0,
     1          81.9134836d0, 82.914136d0, 83.911507d0, 85.91061073d0/
      DATA (ns2(36,i),i=1,6)/0,0,0,9,0,0/
      DATA (ab(36,i),i=1,6)/0.35d0, 2.25d0, 11.6d0, 11.5d0, 57.0d0,
     1                      17.3d0/
c
      DATA at(37),gel(37),nmn(37),(mn(37,i),i=1,2)/'Rb',2,2,85,87/
      DATA (zm(37,i),i=0,2)/85.4678d0, 84.911789738d0, 86.909180527d0/
      DATA (ns2(37,i),i=1,2)/5,3/
      DATA (ab(37,i),i=1,2)/72.165d0, 27.835d0/
c
      DATA at(38),gel(38),nmn(38),(mn(38,i),i=1,4)/'Sr',1,4,84,86,87,88/
      DATA (zm(38,i),i=0,4)/87.62d0, 83.913425d0, 85.9092602d0,
     1                      86.9088771d0, 87.9056121d0/
      DATA (ns2(38,i),i=1,4)/0,0,9,0/
      DATA (ab(38,i),i=1,4)/0.56d0, 9.86d0, 7.00d0, 82.58d0/
c
      DATA at(39),gel(39),nmn(39),(mn(39,i),i=1,1)/' Y',4,1,89/
      DATA (zm(39,i),i=0,1)/88.90585d0, 88.9058483d0/
      DATA (ns2(39,i),i=1,1)/1/
      DATA (ab(39,i),i=1,1)/100.d0/
c
      DATA at(40),gel(40),nmn(40),(mn(40,i),i=1,5)/'Zr',5,5,90,91,92,94,
     1                                              96/
      DATA (zm(40,i),i=0,5)/91.224d0, 89.9047044d0, 90.9056458d0,
     1                      91.9050408d0, 93.9063152d0, 95.9082734d0/
      DATA (ns2(40,i),i=1,5)/0,5,0,0,0/
      DATA (ab(40,i),i=1,5)/51.45d0, 11.22d0, 17.15d0, 17.38d0, 2.80d0/
c
      DATA at(41),gel(41),nmn(41),(mn(41,i),i=1,1)/'Nb',2,1,93/
      DATA (zm(41,i),i=0,1)/92.90638d0, 92.9063781d0/
      DATA (ns2(41,i),i=1,1)/9/
      DATA (ab(41,i),i=1,1)/100.d0/
c
      DATA at(42),gel(42),nmn(42),(mn(42,i),i=1,7)/'Mo',7,7,92,94,95,96,
     1                                              97,98,100/
      DATA (zm(42,i),i=0,7)/95.94d0, 91.906811d0, 93.9050883d0,
     1        94.9058421d0, 95.9046795d0, 96.9060215d0, 97.9054082d0,
     2        99.907477d0/
      DATA (ns2(42,i),i=1,7)/0,0,5,0,5,0,0/
      DATA (ab(42,i),i=1,7)/14.84d0, 9.25d0, 15.92d0, 16.68d0, 9.55d0,
     1                      24.13d0, 9.63d0/
c
      DATA at(43),gel(43),nmn(43),(mn(43,i),i=1,1)/'Tc',6,1,98/
      DATA (zm(43,i),i=0,1)/97.907215d0, 97.907216d0/
      DATA (ns2(43,i),i=1,1)/12/
      DATA (ab(43,i),i=1,1)/100.d0/
c
      DATA at(44),gel(44),nmn(44),(mn(44,i),i=1,7)/'Ru',11,7,96,98,99,
     1                                              100,101,102,104/
      DATA (zm(44,i),i=0,7)/101.07d0, 95.907598d0, 97.905287d0,
     1     98.9059393d0, 99.9042195d0, 100.9055821d0, 101.9043493d0,
     2     103.905433d0/
      DATA (ns2(44,i),i=1,7)/0,0,5,0,5,0,0/
      DATA (ab(44,i),i=1,7)/5.52d0, 1.88d0, 12.7d0, 12.6d0, 17.0d0,
     1                      31.6d0, 18.7d0/
c
      DATA at(45),gel(45),nmn(45),(mn(45,i),i=1,1)/'Rh',10,1,103/
      DATA (zm(45,i),i=0,1)/102.90550d0, 102.905504d0/
      DATA (ns2(45,i),i=1,1)/1/
      DATA (ab(45,i),i=1,1)/100.d0/
c
      DATA at(46),gel(46),nmn(46),(mn(46,i),i=1,6)/'Pd',1,6,102,104,105,
     1                                              106,108,110/
      DATA (zm(46,i),i=0,6)/106.42d0, 101.905609d0, 103.904036d0,
     1       104.905085d0, 105.903486d0, 107.903892d0, 109.905153d0/
      DATA (ns2(46,i),i=1,6)/0,0,5,0,0,0/
      DATA (ab(46,i),i=1,6)/1.02d0, 11.14d0, 22.33d0, 27.33d0, 26.46d0,
     1                      11.72d0/
c
      DATA at(47),gel(47),nmn(47),(mn(47,i),i=1,2)/'Ag',2,2,107,109/
      DATA (zm(47,i),i=0,2)/107.8682d0, 106.905097d0, 108.904752d0/
      DATA (ns2(47,i),i=1,2)/1,1/
      DATA (ab(47,i),i=1,2)/51.839d0, 48.161d0/
c
      DATA at(48),gel(48),nmn(48),(mn(48,i),i=1,8)/'Cd',1,8,106,108,110,
     1                                             111,112,113,114,116/
      DATA (zm(48,i),i=0,8)/112.411d0, 105.906459d0, 107.904184d0,
     1       109.9030021d0, 110.9041781d0, 111.9027578d0, 112.9044017d0,
     2       113.9033585d0, 115.904756d0/
      DATA (ns2(48,i),i=1,8)/0,0,0,1,0,1,0,0/
      DATA (ab(48,i),i=1,8)/1.25d0, 0.89d0, 12.49d0, 12.80d0, 24.13d0,
     1                      12.22d0, 28.73d0, 7.49d0/
c
      DATA at(49),gel(49),nmn(49),(mn(49,i),i=1,2)/'In',2,2,113,115/
      DATA (zm(49,i),i=0,2)/114.818d0, 112.904058d0, 114.903878d0/
      DATA  (ns2(49,i),i=1,2)/9,9/
      DATA (ab(49,i),i=1,2)/4.3d0, 95.7d0/
c
      DATA at(50),gel(50),nmn(50),(mn(50,i),i=1,10)/'Sn',1,10,112,114,
     1                                 115,116,117,118,119,120,122,124/
      DATA (zm(50,i),i=0,10)/118.710d0, 111.904818d0, 113.902779d0,
     1     114.903342d0, 115.901741d0, 116.902952d0, 117.901603d0,
     2     118.903308d0, 119.9021947d0, 121.9034390d0, 123.9052739d0/
      DATA (ns2(50,i),i=1,10)/0,0,1,0,1,0,1,0,0,0/
      DATA (ab(50,i),i=1,10)/0.97d0, 0.65d0, 0.34d0, 14.53d0, 7.68d0,
     1                       24.23d0, 8.59d0, 32.59d0, 4.63d0, 5.79d0/
c
      DATA at(51),gel(51),nmn(51),(mn(51,i),i=1,2)/'Sb',4,2,121,123/
      DATA (zm(51,i),i=0,2)/121.757d0, 120.9038157d0, 122.9042140d0/
      DATA (ns2(51,i),i=1,2)/5,7/
      DATA (ab(51,i),i=1,2)/57.36d0, 42.64d0/
c
      DATA at(52),gel(52),nmn(52),(mn(52,i),i=1,8)/'Te',5,8,120,122,123,
     1                                             124,125,126,128,130/
      DATA (zm(52,i),i=0,8)/127.60d0, 119.904020d0, 121.9030439d0,
     1    122.9042700d0, 123.9028179d0, 124.9044307d0, 125.9033117d0,
     2    127.9044631d0, 129.9062244d0/
      DATA (ns2(52,i),i=1,8)/0,0,1,0,1,0,0,0/
      DATA (ab(52,i),i=1,8)/0.096d0, 2.603d0, 0.908d0, 4.816d0,
     1                      7.139d0, 18.95d0, 31.69d0, 33.80d0/
c
      DATA at(53),gel(53),nmn(53),(mn(53,i),i=1,2)/' I',4,2,127,129/
      DATA (zm(53,i),i=0,2)/126.90447d0, 126.904473d0, 128.904988d0/
      DATA (ns2(53,i),i=1,2)/5,7/
      DATA (ab(53,i),i=1,2)/100.d0,0.d0/
c
      DATA at(54),gel(54),nmn(54),(mn(54,i),i=1,9)/'Xe',1,9,124,126,128,
     1                                          129,130,131,132,134,136/
      DATA (zm(54,i),i=0,9)/131.29d0, 123.9058930d0, 125.904274d0,
     1    127.9035313d0, 128.9047794d0, 129.9035080d0, 130.9050824d0,
     2    131.9041535d0, 133.9053945d0, 135.907219d0/
      DATA (ns2(54,i),i=1,9)/0,0,0,1,0,3,0,0,0/
      DATA (ab(54,i),i=1,9)/0.10d0, 0.09d0, 1.91d0, 26.4d0, 4.1d0,
     1                      21.2d0, 26.9d0, 10.4d0, 8.9d0/
c
      DATA at(55),gel(55),nmn(55),(mn(55,i),i=1,1)/'Cs',2,1,133/
      DATA (zm(55,i),i=0,1)/132.90543d0, 132.905451933d0/
      DATA (ns2(55,i),i=1,1)/7/
      DATA (ab(55,i),i=1,1)/100.d0/
c
      DATA at(56),gel(56),nmn(56),(mn(56,i),i=1,7)/'Ba',1,7,130,132,134,
     1                                             135,136,137,138/
      DATA (zm(56,i),i=0,7)/137.327d0, 129.9063208d0, 131.9050613d0,
     1    133.9045084d0, 134.9056886d0, 135.9045759d0, 136.9058274d0,
     2    137.9052472d0/
      DATA (ns2(56,i),i=1,7)/0,0,0,3,0,3,0/
      DATA (ab(56,i),i=1,7)/0.106d0, 0.101d0, 2.417d0, 6.592d0,
     1                      7.854d0, 11.23d0, 71.70d0/
c
      DATA at(57),gel(57),nmn(57),(mn(57,i),i=1,2)/'La',4,2,138,139/
      DATA (zm(57,i),i=0,2)/138.9055d0, 137.907112d0, 138.9063533d0/
      DATA (ns2(57,i),i=1,2)/10,7/
      DATA (ab(57,i),i=1,2)/0.0902d0, 99.9098d0/
c
      DATA at(58),gel(58),nmn(58),(mn(58,i),i=1,4)/'Ce',9,4,136,138,140,
     1                                             142/
      DATA (zm(58,i),i=0,4)/140.115d0, 135.907172d0, 137.905991d0,
     1    139.9054387d0, 141.909244d0/
      DATA (ns2(58,i),i=1,4)/0,0,0,0/
      DATA (ab(58,i),i=1,4)/0.19d0, 0.25d0, 88.48d0, 11.08d0/
c
      DATA at(59),gel(59),nmn(59),(mn(59,i),i=1,1)/'Pr',10,1,141/
      DATA (zm(59,i),i=0,1)/140.90765d0, 140.9076528d0/
      DATA (ns2(59,i),i=1,1)/5/
      DATA (ab(59,i),i=1,1)/100.d0/
c
      DATA at(60),gel(60),nmn(60),(mn(60,i),i=1,7)/'Nd',9,7,142,143,144,
     1                                             145,146,148,150/
      DATA (zm(60,i),i=0,7)/144.24d0, 141.9077233d0, 142.9098143d0,
     1    143.9100873d0, 144.9125736d0, 145.9131169d0, 147.916893d0,
     2    149.920891d0/
      DATA (ns2(60,i),i=1,7)/0,7,0,7,0,0,0/
      DATA (ab(60,i),i=1,7)/27.13d0, 12.18d0, 23.80d0, 8.30d0, 17.19d0,
     1                       5.76d0, 5.64d0/
c
      DATA at(61),gel(61),nmn(61),(mn(61,i),i=1,1)/'Pm',6,1,145/
      DATA (zm(61,i),i=0,1)/144.912743d0, 144.912749d0/
      DATA (ns2(61,i),i=1,1)/5/
      DATA (ab(61,i),i=1,1)/100.d0/
c
      DATA at(62),gel(62),nmn(62),(mn(62,i),i=1,7)/'Sm',1,7,144,147,148,
     1                                             149,150,152,154/
      DATA (zm(62,i),i=0,7)/150.36d0, 143.911999d0, 146.9148979d0,
     1    147.9148227d0, 148.9171847d0, 149.9172755d0, 151.9197324d0,
     2    153.9222093d0/
      DATA (ns2(62,i),i=1,7)/0,7,0,7,0,0,0/
      DATA (ab(62,i),i=1,7)/3.1d0, 15.0d0, 11.3d0, 13.8d0, 7.4d0,
     1                      26.7d0, 22.7d0/
c
      DATA at(63),gel(63),nmn(63),(mn(63,i),i=1,2)/'Eu',8,2,151,153/
      DATA (zm(63,i),i=0,2)/151.965d0, 150.9198502d0, 152.9212303d0/
      DATA (ns2(63,i),i=1,2)/5,5/
      DATA (ab(63,i),i=1,2)/47.8d0, 52.2d0/
c
      DATA at(64),gel(64),nmn(64),(mn(64,i),i=1,7)/'Gd',5,7,152,154,155,
     1                                              156,157,158,160/
      DATA (zm(64,i),i=0,7)/157.25d0, 151.9197910d0, 153.92086560,
     1    154.9226220d0, 155.9221227d0, 156.9239601d0, 157.9241039d0,
     2    159.9270541d0/
      DATA (ns2(64,i),i=1,7)/0,0,3,0,3,0,0/
      DATA (ab(64,i),i=1,7)/0.20d0, 2.18d0, 14.80d0, 20.47d0, 15.65d0,
     1                      24.84d0, 21.86d0/
c
      DATA at(65),gel(65),nmn(65),(mn(65,i),i=1,1)/'Tb',16,1,159/
      DATA (zm(65,i),i=0,1)/158.92534d0, 158.9253468d0/
      DATA (ns2(65,i),i=1,1)/3/
      DATA (ab(65,i),i=1,1)/100.d0/
c
      DATA at(66),gel(66),nmn(66),(mn(66,i),i=1,7)/'Dy',17,7,156,158,
     1                                           160,161,162,163,164/
      DATA (zm(66,i),i=0,7)/162.50d0, 155.924283d0, 157.924409d0,
     1    159.9251975d0, 160.9269334d0, 161.9267984d0, 162.9287312d0,
     2    163.9291748d0/
      DATA (ns2(66,i),i=1,7)/0,0,0,5,0,5,0/
      DATA (ab(66,i),i=1,7)/0.06d0, 0.10d0, 2.34d0, 18.9d0, 25.5d0,
     1                      24.9d0, 28.2d0/
c
      DATA at(67),gel(67),nmn(67),(mn(67,i),i=1,1)/'Ho',16,1,165/
      DATA (zm(67,i),i=0,1)/164.93032d0, 164.9303221d0/
      DATA (ns2(67,i),i=1,1)/7/
      DATA (ab(67,i),i=1,1)/100.d0/

      DATA at(68),gel(68),nmn(68),(mn(68,i),i=1,6)/'Er',13,6,162,164,
     1                                            166,167,168,170/
      DATA (zm(68,i),i=0,6)/167.26d0, 161.928778d0, 163.929200d0,
     1    165.9302931d0, 166.9320482d0, 167.9323702d0, 169.9354643d0/
      DATA (ns2(68,i),i=1,6)/0,0,0,7,0,0/
      DATA (ab(68,i),i=1,6)/0.14d0, 1.61d0, 33.6d0, 22.95d0, 26.8d0,
     1                      14.9d0/
c
      DATA at(69),gel(69),nmn(69),(mn(69,i),i=1,1)/'Tm',8,1,169/ 
      DATA (zm(69,i),i=0,1)/168.93421d0, 168.9342133d0/
      DATA (ns2(69,i),i=1,1)/1/
      DATA (ab(69,i),i=1,1)/100.d0/
c
      DATA at(70),gel(70),nmn(70),(mn(70,i),i=1,7)/'Yb',1,7,168,170,171,
     1                                            172,173,174,176/
      DATA (zm(70,i),i=0,7)/173.04d0, 167.933897d0, 169.9347618d0,
     1    170.936323580, 171.9363815d0, 172.9382108d0, 173.9388621d0,
     2    175.9425717d0/
      DATA (ns2(70,i),i=1,7)/0,0,1,0,5,0,0/
      DATA (ab(70,i),i=1,7)/0.13d0, 3.05d0, 14.3d0, 21.9d0, 16.12d0,
     1                      31.8d0, 12.7d0/
c
      DATA at(71),gel(71),nmn(71),(mn(71,i),i=1,2)/'Lu',4,2,175,176/
      DATA (zm(71,i),i=0,2)/174.967d0, 174.9407718d0, 175.9426863d0/
      DATA (ns2(71,i),i=1,2)/7,14/
      DATA (ab(71,i),i=1,2)/97.41d0, 2.59d0/
c
      DATA at(72),gel(72),nmn(72),(mn(72,i),i=1,6)/'Hf',5,6,174,176,177,
     1                                             178,179,180/
      DATA (zm(72,i),i=0,6)/178.49d0, 173.940046d0, 175.9414086d0,
     1    176.9432207d0, 177.9436988d0, 178.9458161d0, 179.9465500d0/
      DATA (ns2(72,i),i=1,6)/0,0,7,0,9,0/
      DATA (ab(72,i),i=1,6)/0.162d0, 5.206d0, 18.606d0, 27.297d0,
     1                      13.629d0, 35.100d0/
c
      DATA at(73),gel(73),nmn(73),(mn(73,i),i=1,2)/'Ta',4,2,180,181/
      DATA (zm(73,i),i=0,2)/180.9479d0, 179.9474648d0, 180.9479958d0/
      DATA (ns2(73,i),i=1,2)/16,7/
      DATA (ab(73,i),i=1,2)/0.012d0, 99.988d0/
c
      DATA at(74),gel(74),nmn(74),(mn(74,i),i=1,5)/' W',1,5,180,182,183,
     1                                             184,186/
      DATA (zm(74,i),i=0,5)/183.84d0, 179.946704d0, 181.9482042d0,
     1    182.9502230d0, 183.9509312d0, 185.9543641d0/
      DATA (ns2(74,i),i=1,5)/0,0,1,0,0/
      DATA (ab(74,i),i=1,5)/0.13d0, 26.3d0, 14.3d0, 30.67d0, 28.6d0/
c
      DATA at(75),gel(75),nmn(75),(mn(75,i),i=1,2)/'Re',6,2,185,187/
      DATA (zm(75,i),i=0,2)/186.207d0, 184.9529550d0, 186.9557531d0/
      DATA (ns2(75,i),i=1,2)/5,5/
      DATA (ab(75,i),i=1,2)/37.40d0, 62.60d0/
c
      DATA at(76),gel(76),nmn(76),(mn(76,i),i=1,7)/'Os',9,7,184,186,187,
     1                                             188,189,190,192/
      DATA (zm(76,i),i=0,7)/190.23d0, 183.9524891d0, 185.9538382d0,
     1    186.9557505d0, 187.9558382d0, 188.9581475d0, 189.9584470d0,
     2    191.9614807d0/
      DATA (ns2(76,i),i=1,7)/0,0,1,0,3,0,0/
      DATA (ab(76,i),i=1,7)/0.02d0, 1.58d0, 1.6d0, 13.3d0, 16.1d0,
     1                      26.4d0, 41.0d0/
c
      DATA at(77),gel(77),nmn(77),(mn(77,i),i=1,2)/'Ir',10,2,191,193/
      DATA (zm(77,i),i=0,2)/192.22d0, 190.9605940d0, 192.9629264d0/
      DATA (ns2(77,i),i=1,2)/3,3/
      DATA (ab(77,i),i=1,2)/37.3d0, 62.7d0/
c
c
      DATA at(78),gel(78),nmn(78),(mn(78,i),i=1,6)/'Pt',7,6,190,192,194,
     1                                            195,196,198/
      DATA (zm(78,i),i=0,6)/195.08d0, 189.959932d0, 191.9610380d0,
     1    193.9626803d0, 194.9647911d0, 195.9649515d0, 197.967893d0/
      DATA (ns2(78,i),i=1,6)/0,0,0,1,0,0/
      DATA (ab(78,i),i=1,6)/0.01d0,0.79d0,32.9d0,33.8d0,25.3d0,7.2d0/
c
      DATA at(79),gel(79),nmn(79),(mn(79,i),i=1,1)/'Au',2,1,197/
      DATA (zm(79,i),i=0,1)/196.96654d0, 196.9665687d0/
      DATA (ns2(79,i),i=1,1)/3/
      DATA (ab(79,i),i=1,1)/100.d0/
c
      DATA at(80),gel(80),nmn(80),(mn(80,i),i=1,7)/'Hg',1,7,196,198,199,
     1                                            200,201,202,204/
      DATA (zm(80,i),i=0,7)/200.59d0, 195.965833d0, 197.9667690d0,
     1    198.9682799d0, 199.9683260d0, 200.9703023d0, 201.9706430d0,
     2    203.9734939d0/
      DATA (ns2(80,i),i=1,7)/0,0,1,0,3,0,0/
      DATA (ab(80,i),i=1,7)/0.15d0, 9.97d0, 16.87d0, 23.10d0, 13.18d0,
     1                      29.86d0, 6.87d0/
c
      DATA at(81),gel(81),nmn(81),(mn(81,i),i=1,2)/'Tl',2,2,203,205/
      DATA (zm(81,i),i=0,2)/204.3833d0, 202.9723442d0, 204.9744275d0/
      DATA (ns2(81,i),i=1,2)/1,1/
      DATA (ab(81,i),i=1,2)/29.524d0, 70.476d0/
c
      DATA at(82),gel(82),nmn(82),(mn(82,i),i=1,4)/'Pb',1,4,204,206,207,
     1                                             208/
      DATA (zm(82,i),i=0,4)/207.2d0, 203.9730436d0, 205.9744653d0,
     1    206.9758969d0, 207.9766521d0/
      DATA (ns2(82,i),i=1,4)/0,0,1,0/
      DATA (ab(82,i),i=1,4)/1.4d0, 24.1d0, 22.1d0, 52.4d0/
c
      DATA at(83),gel(83),nmn(83),(mn(83,i),i=1,1)/'Bi',4,1,209/
      DATA (zm(83,i),i=0,1)/208.98037d0, 208.9803987d0/
      DATA (ns2(83,i),i=1,1)/9/
      DATA (ab(83,i),i=1,1)/100.d0/
c
      DATA at(84),gel(84),nmn(84),(mn(84,i),i=1,1)/'Po',5,1,209/
      DATA (zm(84,i),i=0,1)/208.982404d0, 208.9824304d0/
      DATA (ns2(84,i),i=1,1)/1/
      DATA (ab(84,i),i=1,1)/100.d0/
c
      DATA at(85),gel(85),nmn(85),(mn(85,i),i=1,1)/'At',-1,1,210/
      DATA (zm(85,i),i=0,1)/209.987126d0, 209.987148d0/
      DATA (ns2(85,i),i=1,1)/10/
      DATA (ab(85,i),i=1,1)/100.d0/
c
      DATA at(86),gel(86),nmn(86),(mn(86,i),i=1,1)/'Rn',1,1,222/
      DATA (zm(86,i),i=0,1)/222.017571d0, 222.0175777d0/
      DATA (ns2(86,i),i=1,1)/0/
      DATA (ab(86,i),i=1,1)/100.d0/
c
      DATA at(87),gel(87),nmn(87),(mn(87,i),i=1,1)/'Fr',-1,1,223/
      DATA (zm(87,i),i=0,1)/223.019733d0, 223.0197359d0/
      DATA (ns2(87,i),i=1,1)/3/
      DATA (ab(87,i),i=1,1)/100.d0/
c
      DATA at(88),gel(88),nmn(88),(mn(88,i),i=1,1)/'Ra',1,1,226/
      DATA (zm(88,i),i=0,1)/226.025403d0, 226.0254098d0/
      DATA (ns2(88,i),i=1,1)/0/
      DATA (ab(88,i),i=1,1)/100.d0/
c
      DATA at(89),gel(89),nmn(89),(mn(89,i),i=1,1)/'Ac',4,1,227/
      DATA (zm(89,i),i=0,1)/227.027750d0, 227.0277521d0/
      DATA (ns2(89,i),i=1,1)/3/
      DATA (ab(89,i),i=1,1)/100.d0/
c
      DATA at(90),gel(90),nmn(90),(mn(90,i),i=1,1)/'Th',-1,1,232/
      DATA (zm(90,i),i=0,1)/232.038d0, 232.0380553d0/
      DATA (ns2(90,i),i=1,1)/0/
      DATA (ab(90,i),i=1,1)/100.d0/
c
      DATA at(91),gel(91),nmn(91),(mn(91,i),i=1,1)/'Pa',-1,1,231/
      DATA (zm(91,i),i=0,1)/231.03588d0, 231.0358840d0/
      DATA (ns2(91,i),i=1,1)/3/
      DATA (ab(91,i),i=1,1)/100.d0/
c
      DATA at(92),gel(92),nmn(92),(mn(92,i),i=1,4)/' U',-1,4,233,234,
     1                                             235,238/
      DATA (zm(92,i),i=0,4)/238.0289d0, 233.0396352d0, 234.0409521d0,
     1    235.0439299d0, 238.0507882d0/
      DATA (ns2(92,i),i=1,4)/5,0,7,0/
      DATA (ab(92,i),i=1,4)/0.d0, 0.0055d0, 0.7200d0, 99.2745d0/
c
      DATA at(93),gel(93),nmn(93),(mn(93,i),i=1,1)/'Np',-1,1,237/
      DATA (zm(93,i),i=0,1)/237.0481678d0, 237.0481734d0/
      DATA (ns2(93,i),i=1,1)/5/
      DATA (ab(93,i),i=1,1)/100.d0/
c
      DATA at(94),gel(94),nmn(94),(mn(94,i),i=1,1)/'Pu',-1,1,244/
      DATA (zm(94,i),i=0,1)/244.064199d0, 244.064204d0/
      DATA (ns2(94,i),i=1,1)/0/
      DATA (ab(94,i),i=1,1)/100.d0/
c
      DATA at(95),gel(95),nmn(95),(mn(95,i),i=1,1)/'Am',-1,1,243/
      DATA (zm(95,i),i=0,1)/243.061375d0, 243.0613811d0/
      DATA (ns2(95,i),i=1,1)/5/
      DATA (ab(95,i),i=1,1)/100.d0/
c
      DATA at(96),gel(96),nmn(96),(mn(96,i),i=1,1)/'Cm',-1,1,247/
      DATA (zm(96,i),i=0,1)/247.070347d0, 247.070354d0/
      DATA (ns2(96,i),i=1,1)/9/
      DATA (ab(96,i),i=1,1)/100.d0/
c
      DATA at(97),gel(97),nmn(97),(mn(97,i),i=1,1)/'Bk',-1,1,247/
      DATA (zm(97,i),i=0,1)/247.070300d0, 247.070307d0/
      DATA (ns2(97,i),i=1,1)/3/
      DATA (ab(97,i),i=1,1)/100.d0/
c
      DATA at(98),gel(98),nmn(98),(mn(98,i),i=1,1)/'Cf',-1,1,251/
      DATA (zm(98,i),i=0,1)/251.079580d0, 251.079587d0/
      DATA (ns2(98,i),i=1,1)/1/
      DATA (ab(98,i),i=1,1)/100.d0/
c
      DATA at(99),gel(99),nmn(99),(mn(99,i),i=1,1)/'Es',-1,1,252/
      DATA (zm(99,i),i=0,1)/252.082944d0, 252.082980d0/
      DATA (ns2(99,i),i=1,1)/10/
      DATA (ab(99,i),i=1,1)/100.d0/
c
      DATA at(100),gel(100),nmn(100),(mn(100,i),i=1,1)/'Fm',-1,1,257/
      DATA (zm(100,i),i=0,1)/257.095099d0, 257.095105d0/
      DATA (ns2(100,i),i=1,1)/9/
      DATA (ab(100,i),i=1,1)/100.d0/
c
      DATA at(101),gel(101),nmn(101),(mn(101,i),i=1,1)/'Md',-1,1,258/
      DATA (zm(101,i),i=0,1)/258.09857d0, 258.098431d0/
      DATA (ns2(101,i),i=1,1)/16/
      DATA (ab(101,i),i=1,1)/100.d0/
c
      DATA at(102),gel(102),nmn(102),(mn(102,i),i=1,1)/'No',-1,1,259/
      DATA (zm(102,i),i=0,1)/259.100931d0, 259.101030d0/
      DATA (ns2(102,i),i=1,1)/9/
      DATA (ab(102,i),i=1,1)/100.d0/
c
      DATA at(103),gel(103),nmn(103),(mn(103,i),i=1,1)/'Lr',-1,1,260/
      DATA (zm(103,i),i=0,1)/260.105320d0, 260.105500d0/
      DATA (ns2(103,i),i=1,1)/-1/
      DATA (ab(103,i),i=1,1)/100.d0/
c
      DATA at(104),gel(104),nmn(104),(mn(104,i),i=1,1)/'Rf',-1,1,261/
      DATA (zm(104,i),i=0,1)/261.10869d0, 261.108770d0/
      DATA (ns2(104,i),i=1,1)/-1/
      DATA (ab(104,i),i=1,1)/100.d0/
c
      DATA at(105),gel(105),nmn(105),(mn(105,i),i=1,1)/'Db',-1,1,262/
      DATA (zm(105,i),i=0,1)/262.11376d0, 262.114080d0/
      DATA (ns2(105,i),i=1,1)/-1/
      DATA (ab(105,i),i=1,1)/100.d0/
c
      DATA at(106),gel(106),nmn(106),(mn(106,i),i=1,1)/'Sg',-1,1,263/
      DATA (zm(106,i),i=0,1)/263.11822d0, 263.118320d0/
      DATA (ns2(106,i),i=1,1)/-1/
      DATA (ab(106,i),i=1,1)/100.d0/
c
      DATA at(107),gel(107),nmn(107),(mn(107,i),i=1,1)/'Bh',-1,1,262/
      DATA (zm(107,i),i=0,1)/262.12293d0, 262.122890d0/
      DATA (ns2(107,i),i=1,1)/-1/
      DATA (ab(107,i),i=1,1)/100.d0/
c
      DATA at(108),gel(108),nmn(108),(mn(108,i),i=1,1)/'Hs',-1,1,265/
      DATA (zm(108,i),i=0,1)/265.13016d0, 265.130090d0/
      DATA (ns2(108,i),i=1,1)/-1/
      DATA (ab(108,i),i=1,1)/100.d0/
c
      DATA at(109),gel(109),nmn(109),(mn(109,i),i=1,1)/'Mt',-1,1,266/
      DATA (zm(109,i),i=0,1)/266.13764d0, 266.137300d0/
      DATA (ns2(109,i),i=1,1)/-1/
      DATA (ab(109,i),i=1,1)/100.d0/
c
      IF((IAN.LE.0).OR.(IAN.GT.109)) THEN
          MASS= 0.d0
          NAME= 'XX'
          IMN= 0
          WRITE(6,601) IAN
          RETURN
        ELSE
          NAME= AT(IAN)
        ENDIF
      IF((IAN.EQ.1).AND.(IMN.NE.1)) THEN
c** Special case: insert common name for deuterium or tritium
          IF(IMN.EQ.2) NAME=' D'
          IF(IMN.EQ.3) NAME=' T'
          ENDIF
      GELGS= GEL(IAN)
      MASS= -1.d0
      GNS= -1
        ABUND = -1.d0
      DO  I= 1,NMN(IAN)
          if(i.gt.10)  write(6,606) ian,imn,nmn(ian)
  606  format(3i9)
          IF(IMN.EQ.MN(IAN,I)) THEN
              MASS= ZM(IAN,I)
              GNS= NS2(IAN,I)+1
              ABUND = AB(IAN,I)
              ENDIF
          ENDDO
      IF(MASS.LT.0.d0) THEN
          MASS= ZM(IAN,0)
          IF(IMN.NE.0) WRITE(6,602) AT(IAN),IMN
          IMN= 0
          ENDIF
      RETURN
  601 FORMAT(' *** MASSES Data base does not include Atomic Number=',i4)
  602 FORMAT(' *** MASSES Data base does not include ',A2,'(',i3,
     1 '), so use average atomic mass.')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PREPOT(LNPT,IAN1,IAN2,IMN1,IMN2,NPP,OMEGA,RR,RM2,VLIM,
     1                                                        VV,NCN)
c** Driver subroutine of package to read parameters and/or generate
c  values of a potential V(I) at the NPP input distances RR(I).
c=================== Version of  9 February 2004 =======================
c**** Subroutine Input:
c----------------------
c  LNPT  is an integer specifying the operational mode:
c      *  LNPT > 0  : for a new case for which all potential-defining
c                     parameters are read in & a description printed
c      *  LNPT.le.0 : if potential points are to be generated in exactly
c                     the same manner as on preceding call, but at
c                     different distances RR(I) (no reads or writes)
c  IAN1 & IAN2 are the atomic numbers and IMN1 & IMN2 the mass numbers
c        of atoms #1 & 2, used (if needed) to specify isotope masses for
c        calculating adiabatic and/or non-adiabatic B-O-B correction fx.
c  NPP (integer) is the number of input distances  RR(i) (in Angstroms)
c        at which potential values  VV(i) (in cm-1) are to be generated
c  RR  (real array) is set of NPP distances where potential calculated
c  RM2 (real array) on input is the (centrifugal) array of  1/RR(i)**2
c----------------------
c**** Subroutine Output:
c----------------------
c  OMEGA the (integer) elelectronic angular momentum projection q.no.
c  VLIM (cm-1) is the absolute energy at the potential asymptote
c  VV (real array) is the set of function values generated (in cm-1)
c  RM2 values returned may (if appropriate) be modified to include B-O-B
c      corrections to the (centrifugal) potential  1/RR(i)**2
c  NCN is an integer power defining the asymptotically-dominant
c      inverse-power long-range potential tail:  CNN/R**NCN
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Calls GENINT (which calls PLYINTRP, SPLINT & SPLINE) ,  or POTGEN ++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set maximum array dimension for the input function values to be
c  interpolated over & extrapolated beyong
      INTEGER NTPMX
      PARAMETER (NTPMX= 1600)
      INTEGER I,J,IAN1,IAN2,IMN1,IMN2,INPTS,ILR,IR2,JWR,LNPT,LPPOT,LWR,
     1  NCN,NLIN,NPP,NROW,NTP,NUSE,NPRS,NPRF, OMEGA
      REAL*8  RFACT,EFACT,RH,RMIN,VLIM,VSHIFT,VV(NPP),RR(NPP),RM2(NPP),
     1  XI(NTPMX),YI(NTPMX),RWR(20),VWR(20),VWRB(3),D1V(3),D1VB(3),
     2  D2V(3),CNN
c
c** Save variables needed for 'subsequent' LNPT.le.0 calls
      SAVE ILR,IR2,LPPOT,NTP,NUSE
      SAVE CNN,VSHIFT,XI,YI
c
      DATA VWRB/3*0.D0/,D1VB/3*0.D0/
      LPPOT= 0
      NPRS= 1
      NPRF= NPP
c
      IF(LNPT.GT.0) THEN
c** If NTP > 0 :  define potential by interpolation over & extrapolation
c        beyond the NTP read-in turning points using subroutine GENINT.
c   If NTP.le.0 : generate a (fully analytic) potential in POTGEN.
c** If LPPOT > 0 : at every |LPPOT|-th point, print potential and
c      derivatives-by-differences. ***  If  LPPOT < 0  write potential
c      at every |LPPOT|-th point to channel-8 in a compact format **
c  OMEGA  the (integer) total elextronic angular momentum projection
c            quantum number (required for proper rotational intensities)
c** VLIM (cm-1) is the energy associated with the potential asymptote.
c-----------------------------------------------------------------------
          READ(5,*) NTP, LPPOT, OMEGA, VLIM
c-----------------------------------------------------------------------
          WRITE(6,600) OMEGA,VLIM
          IF(NTP.GT.0) THEN
c** For a pointwise potential (NTP > 0), now read points & parameters
c  controlling how the interpolation/extrapolation is to be done.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** NTP (read above) is number of turning points (XI,YI) to be read in.
c** If NUSE > 0  interpolate with NUSE-point piecewise polynomials
c    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE.LE.0)
c    interpolate with cubic spline instead of local polynomials.
c** If IR2 > 0 , interpolate over  YI*XI**2 ; otherwise on  YI  itself
c   This may help if interpolation has trouble on steep repulsive wall.
c** ILR specifies how to extrapolate beyond largest input distance XI(i)
c  If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c  If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c  If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
c  inverse-power terms beginning with  1/R**NCN}. *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  [(cm-1)(Angstroms)**'NCN'].
c* If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c* If ILR > 3 : successive higher power terms differ by factor  1/R
c-----------------------------------------------------------------------
              READ(5,*) NUSE, IR2, ILR, NCN, CNN
c-----------------------------------------------------------------------
              IF(NTP.GT.NTPMX) THEN
                  WRITE(6,602) NTP,NTPMX
                  STOP
                  ENDIF
              IF(NUSE.GT.0) WRITE(6,604) NUSE,NTP
              IF(NUSE.LE.0) WRITE(6,606) NTP
              IF(IR2.GT.0) WRITE(6,608)
              IF((ILR.GT.1).AND.(DABS(CNN).GT.0.D0))WRITE(6,610)CNN,NCN
c** Read in turning points to be interpolated over
c** RFACT & EFACT are factors required to convert units of input turning
c       points (XI,YI) to Angstroms & cm-1, respectively (may = 1.d0)
c** Turning points (XI,YI) must be ordered with increasing XI(I)
c** Energy VSHIFT (cm-1) is added to the input potential points to
c   make their absolute energy consistent with VLIM (often VSHIFT=Te).
c-----------------------------------------------------------------------
              READ(5,*) RFACT, EFACT, VSHIFT
              READ(5,*) (XI(I), YI(I), I= 1,NTP)
c-----------------------------------------------------------------------
              WRITE(6,612) VSHIFT, RFACT, EFACT
              NROW= (NTP+2)/3
              DO  J= 1,NROW
                  IF(EFACT.LE.10.D0) THEN
                      WRITE(6,614) (XI(I),YI(I),I= J,NTP,NROW)
                    ELSE
                      WRITE(6,616) (XI(I),YI(I),I= J,NTP,NROW)
                    ENDIF
                  ENDDO
                  WRITE(6,624)
              DO  I= 1,NTP
                  YI(I)= YI(I)*EFACT+ VSHIFT
                  XI(I)= XI(I)*RFACT
                  ENDDO
              IF(IR2.GT.0) THEN
                  DO  I= 1,NTP
                      YI(I)= YI(I)*XI(I)**2
                      ENDDO
                  ENDIF
              ENDIF
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(NTP.GT.0) THEN
          CALL GENINT(LNPT,NPP,RR,VV,NUSE,IR2,NTP,XI,YI,VLIM,ILR,
     1                                              NCN,CNN,NPRS,NPRF)
        ELSE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Generate a fully analytic potential in subroutine POTGEN ***********
c* Potentials generated in cm-1 with potential asymptote at energy VLIM
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** IPOTL specifies the type of potential function to be generated.
c** MPAR & NPAR are integers for specifying potential types.
c** NVARB is number of (real*8) potential parameters read in.
c** IBOB specifies whether (if > 0) or not (if .le. 0) atomic mass
c      dependent Born-Oppenheimer breakdown corrections will be included
c** For all functions considered, well depth and equilibrium distance
c  are read as  DSCM (cm-1)  and REQ (Angstroms), respectively.
c* [Most read-in parameters are dimensionless (scaled by DSCM & REQ).]
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** If IPOTL=1  generate an L.J.(MPAR,NPAR) potential.
c** If IPOTL=2  generate an MLJ(NPAR) potential [JCP 112, 3949 (2000)]
c      If MPAR > 0  exponent parameter is polynomial of order (NVARB-1)
c           in y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR),
c           with the NVARB coefficients  PARM(j)
c      If MPAR.le.0  exponent polynomial in  y_{1} of order (NVARB-4) with
c           coefficients PARM(i) (i= 1,NVARB-3), & includes a switching
c           function with exponent coefficient  ALPHA= PARM(NVARB)  and 
c           RSW= PARM(NVARB-1),  defined to yield limiting inverse-power
c           potential coefficient  Cn= PARM(NVARB-2).
c** If IPOTL=3  generate a Morse or Extended Morse Oscillator potential
c      with exponent factor "beta" defined as a power series of order
c      (NVARB-1) in  y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR)
c      with NVARB coefficients PARM(i).    [!! MPAR .ge.1 !!]
c    * For conventional "simple" Morse potential,  NVARB=1 & MPAR dummy
c*  Special option #1: set  MPAR= -1  to produce Wei Hua's 4-parameter
c      modified Morse function with  b= PARM(1)  and C= PARM(2).
c*  Special option #2: set  MPAR= -2  to produce Coxon's "Generalized
c      Morse Oscillator" potential with exponent expansion in (R-Re)]
c ...  otherwise, set  MPAR.ge.0
c** If IPOTL=4  use Seto's modification of Surkus' GPEF expansion in
c       z = [R^NPAR - Re^NPAR]/[a*R^NPAR + b*Re^NPAR] where
c       a=PARM(NVARB-1) & b=PARM(NVARB), which incorporates Dunham, SPF,
c       O-T and other forms: V(z) = c_0 z^2 [1 + c_1 z + c_2 z^2 + ...]
c       where  c_0 [cm-1] is read in as DSCM, and the first (NVARB-2)
c       PARM(i)'s are the  c_i  (i > 0).  [MPAR is dummy parameter here]
c  * For Dunham case:  NPAR=1, PARM(NVARB-1)= 0.0, PARM(NVARB)= 1.0
c  * For SPF case:  NPAR=1, PARM(NVARB-1)= 1.0, PARM(NVARB)= 0.0
c  * For Ogilvie-Tipping:  NPAR=1, PARM(NVARB-1)= 0.5 = PARM(NVARB)
c  * NOTE that for Surkus NPAR < 0 case:  z(NPAR,a,b)= z(|NPAR|,-b,-a)
c      Generate & return the  D_e  value implied by these coefficients.
c  * NPAR= 0 generates potential as power series of order-NVARB in R
c    with constant term = (read-in VLIM) & NVARB read-in coefficients
c** If IPOTL=5  generate generalized HFD(NPAR,6,8,10,12,14) potential.
c       PARM(1-3) are the parameters defining the HFD damping function
c       D(x)=exp[-pparm(1)*(PARM(2)/x - 1)**PARM(3)] {for x < PARM(2)}
c       PARM(4) the quadratic coefficient in the exponent, and
c       PARM(5) is the power of  x=R/Req  multiplying the repulsive term
c              AREP*x**PARM(5) *exp[-beta*x - PARM(4)*x**2] ;
c       PARM(6-11)  are the reduced C_NPAR, C_6, C_8, C_10, C_12 and C14
c       parameters (NPAR < 6), while  AREP and  beta  are defined
c       by having the potential minimum at  x=1.  For NVARB < 11, higher
c       C_m coefficients automatically zero;  necessarily  NVARB.ge.7 .
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** IBOB > 0, add atomic-mass-dependent Born-Openheimer breakdown
c  correction functions to rotationless and/or centrifugal potential(s).
c  Both expressed as power series in  z= (R-Re)/(R+Re) starting with the
c  constant term, using the mass shift convention of Le Roy [J.Mol.Spec.
c  194, 189 (1999)].  Adiabatic B-O-B potential correction fx. defined
c  by polynomials of order NC1 with (NC1+1) coefficients {CA1(i)} for
c  atom-1  and order NC2 with (NC2+1) coefficients {CA2(i)} for atom-2,
c  while centrifugal correction fx. defined polynomial of order NG1 with
c  (NG1+1) coefficients {GA1(i)} for atom-1 and order NG2 with (NG2+1)
c  coefficients {GA2(i)} for atom-2.
c** Input parameters IANi & IMNi are the atomic & mass number of atom-i
c  (i=1,2), while integers RMN1 & RMN2 read here are the mass numbers of
c  the reference isotopes defining the B-O-B correction functions.
c** NC1 & NC2 are orders of polynomials DELTA(V,atom-i) defining
c  'adiabatic' corrections to the rotationless potential for atoms 1 & 2
c  DELTA(V)= (1-M1ref/M1)*DELTA(V,atom-1) + (1-M2ref/M2)*DELTA(V,atom-2)
c** NG1 & NG2 are orders of polynomials q1(z) & q2(z) defining B-O-B
c   correction to the centrifugal potential:
c      V(centrifugal)= [1 + (M1ref/M1)*q1(z) + (M2ref/M2)*q2(z)]/R**2
c ... to omit a particular correction set associated NCi or NGi .lt.0
c** RX > 0.0  invokes Coxon's (older) expansions in (R-Re) for potential
c     correction and in  [(R-Rx)**j - (Re-Rx)**j] for centrifugal corrn.
c ... OTHERWISE (to use Le Roy B-O-B formalism) set  RX.le.0.d0 !!
c-----------------------------------------------------------------------
c** Read inside subroutine POTGEN
c         IF(LNPT.GT.0) THEN
c             READ(5,*) IPOTL, MPAR, NPAR, NVARB, IBOB, DSCM, REQ
c             IF(NVARB.GT.0) READ(5,*) (PARM(I), I=1,NVARB)
c             IF(IBOB.GT.0) THEN
c                 READ(5,*) RMN1, RMN2, NC1, NC2, NG1, NG2, RX
c                 IF(NC1.GE.0) READ(5,*) (CA1(I), I=0,NC1)
c                 IF(NC2.GE.0) READ(5,*) (CA2(I), I=0,NC2)
c                 IF(NG1.GE.0) READ(5,*) (GA1(I), I=0,NG1)
c                 IF(NG2.GE.0) READ(5,*) (GA2(I), I=0,NG2)
c             ENDIF
c         ENDIF
c-----------------------------------------------------------------------
          NCN= 99
          CALL POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,RR,RM2,VV,
     1                                                        NCN,CNN)
        ENDIF
      IF(LPPOT.NE.0) THEN
c** If desired, on the first pass (i.e. if LNPT > 0) print the potential
          RH= RR(2)-RR(1)
          INPTS= IABS(LPPOT)
          IF(LPPOT.LT.0) THEN
c** Option to write resulting function compactly to channel-8.
              RMIN= RR(1)
              NLIN= NPP/INPTS+ 1
              WRITE(8,800) NLIN,VLIM
              WRITE(8,802) (RR(I),VV(I),I= 1,NPP,INPTS)
            ELSE
c** Option to print potential & its 1-st three derivatives, the latter
c  calculated by differences, assuming equally spaced RR(I) values.
              WRITE(6,620)
              NPRS= MAX(1,(NPRS- 15*INPTS))
              NPRF= MIN(NPP,(NPRF+ 15*INPTS))
              NLIN= (NPRF-NPRS+1)/(2*INPTS)+1
              RH= INPTS*RH
              DO  I= 1,NLIN
                  LWR= NPRS+INPTS*(I-1)
                  DO  J= 1,2
                      JWR= LWR+(J-1)*NLIN*INPTS
                      IF(JWR.LE.NPP) THEN
                          RWR(J)= RR(JWR)
                          VWR(J)= VV(JWR)
                          D1V(J)= (VWR(J)-VWRB(J))/RH
                          VWRB(J)= VWR(J)
                          D2V(J)= (D1V(J)-D1VB(J))/RH
                          D1VB(J)= D1V(J)
                        ELSE
                          RWR(J)= 0.d0
                          VWR(J)= 0.d0
                        ENDIF
                      IF(I.LE.2) THEN
                          D2V(J)= 0.d0
                          IF(I.EQ.1) D1V(J)= 0.d0
                          ENDIF
                      ENDDO
                  WRITE(6,622) (RWR(J),VWR(J),D1V(J),D2V(J),J= 1,2)
                  ENDDO
            ENDIF
          ENDIF
      IF(LNPT.GT.0) WRITE(6,624)
      RETURN
  600 FORMAT(' State has  OMEGA=',i2, '   and energy asymptote:   Y(lim)
     1=',F12.4,'(cm-1)')
  602 FORMAT(/' **** ERROR in dimensioning of arrays required'
     1 ,' by GENINT;   No. input points ',I5,' > NTPMX =',I4)
  604 FORMAT('- Perform',I3,'-point piecewise polynomial interpolation o
     1ver',I5,' input points' )
  606 FORMAT('- Perform cubic spline interpolation over the',I5,
     1  ' input points' )
  608 FORMAT('- Interpolation actually performed over modified input arr
     1ay:   Y(I) * R(I)**2')
  610 FORMAT('- Beyond read-in points extrapolate to limiting asymptotic
     1 behaviour:'/20x,'Y(R)  =  Y(lim) - (',D16.7,')/R**',I2)
  612 FORMAT('- To make input points Y(i) consistent with  Y(lim),  add'
     1 ,'  Y(shift)=',F12.4/'- Scale input points:  (distance)*',
     2 1PD16.9,'  &  (energy)*',D16.9/13x,'to get required internal unit
     3s  [Angstroms & cm-1 for potentials]'/
     4  3('      R(i)         Y(i)  ')/3(3X,11('--')))
  614 FORMAT((3(F13.8,F12.4)))
  616 FORMAT((3(F12.6,F13.8)))
  620 FORMAT(/'  Function and first 2 derivatives by differences'/
     1  2('     R       Y(R)     d1Y/dR1    d2Y/dR2')/2(2X,19('--')))
  622 FORMAT(2(0PF8.3,F11.3,1PD11.3,D10.2))
  624 FORMAT(1x,38('--'))
  800 FORMAT(I7,' function values with asymptotic value:',F14.6)
  802 FORMAT((1X,3(F12.8,F14.6)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE GENINT(LNPT,NPP,XX,YY,NUSE,IR2,NTP,XI,YI,VLIM,ILR,
     1                                              NCN,CNN,NPRS,NPRF)
c** GENINT produces a smooth function YY(i) at the NPP input distances
c  XX(i) by performing numerical interpolation over the range of the
c  NTP input function values YI(j) at the distances XI(j), and using
c  analytic functions to extrapolate beyond their range to with an
c  exponential at short range and a form specified by ILR, NCN & CNN
c** ILR specifies how to extrapolate beyond largest given turning pts
c   If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c   If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c   If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c* If(ILR.ge.2) fit last turning points to:  VLIM - sum(of ILR
c  inverse-power terms beginning with  1/R**NCN). *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  ((cm-1)(Angstroms)**'NCN').
c  If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c  If ILR > 3 : this factor is  1/R .
c=== Calls subroutines PLYINTRP, SPLINT & SPLINE ==================
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,J,IFXCN,IDER,IR2,ILR,ISR,LNPT,MBEG,MFIN,MINNER,
     1  NN,NPP,NUSE,NUST,NORD,NCN,NCN2,NCN4,NPRS,NPRF,NTP,
     2  IMX1,NMX,JR2,JMAX,MI(10),MF(10)
      REAL*8  ASR,BSR,CSR,ALR,BLR,CLR,DCSR,ADCSR,PDCSR,VRAT,
     1  DX1,DX2,DX3,EX1,EX2,EX3,CNN,VLIM,X1,X2,X3,Y1,Y2,Y3,
     1  XX(NPP),YY(NPP),XI(NTP),YI(NTP),XJ(20),YJ(20),DUMM(20)
c
      SAVE ASR,BSR,CSR,ISR,ALR,BLR,CLR,IMX1,NMX,JR2,JMAX
c
      NUST= NUSE/2
      IF(NUSE.LE.0) NUST= 2
      IDER= 0
      NPRS= 1
      NPRF= NPP
c** Determine if/where need to begin extrapolation beyond input data
c  XX(MI(J))  is the 1-st mesh point past turning point  XI(J) .
c  XX(MF(J))  is the last mesh point before turning point  XI(NTP+1-J)
      DO 6 J = 1,NUST
          MI(J)= 1
          MF(J)= 0
          DO  I= 1,NPP
              IF(XX(I).LE.XI(J)) MI(J)= I+ 1
              IF(XX(I).GE.XI(NTP+1-J)) GO TO 6
              MF(J)= I
              ENDDO
    6     CONTINUE
      IF(NUST.LT.2) THEN
          MFIN= MI(1)-1
        ELSE
          MFIN= MI(2)-1
        ENDIF
      IF(LNPT.GT.0) THEN
c-----------------------------------------------------------------------
c** For a new case determine analytic functions for extrapolating beyond
c  the range of the input points (if necessary) on this or later calls.
c** Try to fit three innermost turning points to  V(R)=A+B*DEXP(-C*R).
c** If unsatisfactory, extrapolate inward with inverse power function
          IF(IR2.LE.0) THEN
              DO  I= 1,4
                  YJ(I)= YI(I)
                  ENDDO
            ELSE
              DO  I= 1,4
                  YJ(I)= YI(I)/XI(I)**2
                  ENDDO
            ENDIF
          X1= XI(1)
          X2= XI(2)
          X3= XI(3)
          Y1= YJ(1)
          Y2= YJ(2)
          Y3= YJ(3)
          IF((Y1-Y2)*(Y2-Y3).LE.0.d0) THEN
c** If 3 innermost points not monotonic, use A+B/X inward extrapoln.
              ISR= 0
              WRITE(6,600)
            ELSE
c** Use cubic through innermost points to get initial trial exponent
c  from ratio of derivatives,  Y''/Y'
              IDER= 2
              ISR= 4
              CALL PLYINTRP(XI,YJ,ISR,X2,XJ,ISR,IDER)
              CSR= XJ(3)/XJ(2)
              DCSR= DABS(CSR*X2)
              IF(DCSR.GT.1.5D+2) THEN
c** If exponential causes overflows, use inverse power inward extrapoln.
                  ISR= 0
                  WRITE(6,602) CSR
                  GO TO 20
                  ENDIF
c** Prepare parameters for inward exponential extrapolation
              VRAT= (Y3- Y2)/(Y1- Y2)
              DX1= X1- X2
              DX3= X3- X2
              EX2= 1.D0
              ADCSR= 1.d99
c** Now iterate (with actual point) to get exact exponent coefficient
              DO  J= 1,15
                  PDCSR= ADCSR
                  EX1= DEXP( CSR*DX1)
                  EX3= DEXP( CSR*DX3)
                  DCSR= (VRAT- (EX3- EX2)/(EX1- EX2)) /
     1   ((X3*EX3- X2 - (X1*EX1- X2)*(EX3-EX2)/(EX1- EX2))/(EX1- EX2))
                  ADCSR= ABS(DCSR)
                  IF((ADCSR.GT.PDCSR).AND.(ADCSR.LT.1.d-8)) GO TO 12
                  IF(ADCSR.LT.1.d-12) GO TO 12
                  CSR= CSR+ DCSR
                  ENDDO
              WRITE(6,604) DCSR
   12         BSR= (Y1-Y2)/(EX1-EX2)
              ASR= Y2-BSR*EX2
              BSR= BSR*DEXP(-CSR*X2)
              WRITE(6,606) X2,ASR,BSR,CSR
            ENDIF
   20     IF(ISR.LE.0) THEN
              IF((X1*X2).LE.0.d0) THEN
c** If 1'st two mesh points of opposite sign, extrapolate linearly
                  ISR= -1
                  ASR= Y2
                  BSR= (Y2- Y1)/(X2- X1)
                  CSR= X2
                  WRITE(6,608) X2,ASR,BSR,CSR
                ELSE
c** For inward extrapolation as inverse power through 1'st two points ..
                  BSR= (Y1-Y2)* X1*X2/(X2- X1)
                  ASR= Y1-BSR/X1
                  CSR= X2
                  WRITE(6,610) X2,ASR,BSR
                ENDIF
              ENDIF
          ENDIF
  600 FORMAT('  ** CAUTION ** Exponential inward extrapolation fails'/
     1 16x,'since first 3 points not monotonic, ... so ...')
  602 FORMAT(' *** CAUTION ** inward extrapolation exponent coefficient
     1   C=',D12.4/10x,'could cause overflows, ... so ...')
  604 FORMAT(' *** CAUTION ** after 15 tries inward extrap. exponent coe
     1fft change is',1PD9.1)
  606 FORMAT(' Extrapolate to   X .le.',F7.4,'  with'/'   Y=',F13.3,
     1  SP,1PD15.6,' * exp(',SS,D13.6,'*X)')
  608 FORMAT(' Extrapolate to   X .le.',F8.4,'   with'/'   Y=',F13.3,
     1  SP,1PD16.7,' * [X - (',SS,F8.4,')]')
  610 FORMAT(' Extrapolate to  X .le.',F8.4,'   with   Y=',F12.3,
     1  SP,1PD15.6,')/X**1')
c
      IF(MFIN.GT.0) THEN
c** If needed, calculate function in inner extrapolation region
          IF(ISR.GT.0) THEN
c ... either as an exponential
              DO  I= 1,MFIN
                  EX1= CSR*XX(I)
                  IF(DABS(EX1).GT.1.D+2) EX1= 1.D+2*DSIGN(1.d0,EX1)
                  YY(I)= ASR+BSR*DEXP(EX1)
                  ENDDO
            ELSEIF(ISR.EQ.0) THEN
c ... or if that fails, as an inverse power
              DO  I= 1,MFIN
                  YY(I)= ASR+BSR/XX(I)
                  ENDDO
            ELSEIF(ISR.LT.0) THEN
c ... or if X changes sign, extrapolate inward linearly
              DO  I= 1,MFIN
                  YY(I)= ASR+ BSR*(XX(I)- CSR)
                  ENDDO
            ENDIF
          ENDIF
c** End of inward extrapolation procedure
c-----------------------------------------------------------------------
      MINNER= MFIN
      IF(NUST.GT.2) THEN
c** If(NUSE.gt.5) minimize spurious behaviour by interpolating with
c  order less than NUSE on intervals near inner end of range
          DO  J= 3,NUST
              NORD= 2*(J-1)
              MBEG= MI(J-1)
              MFIN= MI(J)-1
              IF(MFIN.GE.MBEG) THEN
                  DO  I=  MBEG,MFIN
                      CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
                      YY(I)= DUMM(1)
                      ENDDO
                  ENDIF
              ENDDO
          ENDIF
c** Main interpolation step begins here
c=======================================================================
      MBEG= MI(NUST)
      MFIN= MF(NUST)
      IF(MFIN.GE.MBEG) THEN
          IF(NUSE.LE.0) THEN
c** Either ... use cubic spline for main interpolation step
              CALL SPLINT(LNPT,NTP,XI,YI,MBEG,MFIN,XX,YY)
            ELSE
c ... or use piecewise polynomials for main interpolation step
              DO  I= MBEG,MFIN
                  CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NUSE,IDER)
                  YY(I)= DUMM(1)
                  ENDDO
            ENDIF
          ENDIF
      IF(MFIN.LT.NPP) THEN
          IF(NUST.LE.2) THEN
c** If(NUSE.gt.5) minimize spurious behaviour by interpolating with
c  order less than NUSE on intervals near outer end of range
              MBEG= MF(NUST)+1
            ELSE
              NN= NUST-2
              DO  J= 1,NN
                  NORD= 2*(NUST-J)
                  MBEG= MF(NUST-J+1)+1
                  MFIN= MF(NUST-J)
                  IF(MFIN.GE.MBEG) THEN
                      DO  I= MBEG,MFIN
                          CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
                          YY(I)= DUMM(1)
                          ENDDO
                      END IF
                  ENDDO
            ENDIF
          ENDIF
      MBEG= MFIN+1
      IF((MFIN.GT.MINNER).AND.(IR2.GT.0)) THEN
c** In (IR2.gt.0) option, now remove X**2 from the interpolated function
          DO  I= MINNER+1,MFIN
              YY(I)= YY(I)/XX(I)**2
              ENDDO
          ENDIF
c** Print test of smoothness at join with analytic inward extrapolation
c     IF(LNPT.GT.0) THEN
c         MST= MAX0(MINNER-4,1)
c         NPRS= MST
c         MFN= MST+8
c         IF(MFN.GT.NPP) MFN= NPP
c         IF(MFN.GT.MFIN) MFN= MFIN
c         IF(MINNER.GT.0) WRITE(6,611) X2,((XX(I),YY(I),I= J,MFN,3),
c    1        J= MST,MST+2)
c 611 FORMAT('     Verify smoothness of inner join at   X=',F9.5/
c    1  (3X,3(F10.5,G15.7)))
c         ENDIF
c-----------------------------------------------------------------------
c** To extrapolate potential beyond range of given turning points ...
      IF(LNPT.GT.0) THEN
c** On first entry, calculate things needed for extrapolation constants
          Y1= YI(NTP)
          Y2= YI(NTP-1)
          Y3= YI(NTP-2)
          X1= XI(NTP)
          X2= XI(NTP-1)
          X3= XI(NTP-2)
          IF(IR2.GT.0) THEN
              Y1= Y1/X1**2
              Y2= Y2/X2**2
              Y3= Y3/X3**2
              ENDIF
          ENDIF
c** Check inverse-power tail power ...
      IF(NCN.LE.0) NCN= 6
      IF(ILR.LT.0) THEN
          IF(LNPT.GT.0) THEN
C** For  ILR.lt.0  use  Y = VLIM - ALR * exp[-CLR*(X - BLR)**2]
              EX1= DLOG((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
              EX2= DLOG((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
              BLR= (X1+X2 - (X2+X3)*EX1/EX2)/(2.d0- 2.d0*EX1/EX2)
              CLR= -EX1/(X1+X2-2.d0*BLR)
              ALR= (VLIM-Y1)*DEXP(CLR*(X1-BLR)**2)
              WRITE(6,614) X2,VLIM,ALR,CLR,BLR
              IF(CLR.LT.0.d0) THEN
c ... but replace it by an inverse power of exponent constant negative
                  WRITE(6,612)
                  ILR= 1
                  GO TO 50
                  ENDIF
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM- ALR*DEXP(-CLR*(XX(I) - BLR)**2)
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.EQ.0) THEN
c** For ILR.le.0  use  Y = VLIM - ALR * X**p * exp(-CLR*X)
          IF(LNPT.GT.0) THEN
              EX1= DLOG((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
              EX2= DLOG((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
              DX1= DLOG(X1/X2)/(X1-X2)
              DX2= DLOG(X2/X3)/(X2-X3)
              BLR= (EX1-EX2)/(DX1-DX2)
              CLR= BLR*DX1- EX1
              ALR= (VLIM-Y1)* DEXP(CLR*X1)/X1**BLR
              WRITE(6,616) X2,VLIM,ALR,BLR,CLR
              IF(CLR.LT.0.d0) THEN
c ... but replace it by an inverse power of exponent constant negative
                  WRITE(6,612)
                  ILR= 1
                  GO TO 50
                  ENDIF
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM- ALR*XX(I)**BLR *DEXP(-CLR*XX(I))
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
   50 IF(ILR.EQ.1) THEN
c** For  ILR=1 ,  use     Y = VLIM + ALR/X**BLR
          IF(LNPT.GT.0) THEN
              BLR= DLOG((VLIM-Y2)/(VLIM-Y1))/DLOG(X1/X2)
              ALR= (Y1- VLIM)*X1**BLR
              NCN= BLR
              WRITE(6,618) X2,VLIM,ALR,BLR,NCN
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM+ ALR/XX(I)**BLR
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
c** Set constants for long-range extrapolation
      IFXCN= 0
      IF((CNN.GT.0.d0).OR.(CNN.LT.0.d0)) IFXCN= 1
      NCN2= NCN+2
      IF(ILR.EQ.2) THEN
c** For ILR=2 ,  use   Y = VLIM - CNN/X**NCN - BSR/X**(NCN+2)
c*  If CNN held fixed need ILR > 2  to prevent discontinuity
          IF(LNPT.GT.0) THEN
              IF(IFXCN.LE.0) THEN
                  CNN= ((VLIM-Y1)*X1**NCN2 -
     1                 (VLIM-Y2)*X2**NCN2)/(X1**2-X2**2)
                  ENDIF
              ALR= CNN
              BLR= (VLIM-Y1)*X1**NCN2 - CNN*X1**2
              WRITE(6,620) X2,VLIM,CNN,NCN,BLR,NCN2
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM-(ALR+BLR/XX(I)**2)/XX(I)**NCN
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.EQ.3) THEN
c** For ILR=3 , use   Y = VLIM - (CN + CN2/X**2 + CN4/X**4)/X**NCN
          IF(LNPT.GT.0) THEN
              NCN4= NCN+4
              IF(IFXCN.GT.0) THEN
                  ALR= CNN
                  BLR= (((VLIM-Y1)*X1**NCN-ALR)*X1**4-((VLIM-Y2)
     1                     *X2**NCN-ALR)*X2**4)/(X1**2-X2**2)
                  CLR= ((VLIM-Y1)*X1**NCN-ALR-BLR/X1**2)*X1**4
                ELSE
                  EX1= X1**2
                  EX2= X2**2
                  EX3= X3**2
                  DX1= (VLIM-Y1)*X1**NCN4
                  DX2= (VLIM-Y2)*X2**NCN4
                  DX3= (VLIM-Y3)*X3**NCN4
                  BLR= (DX1-DX2)/(EX1-EX2)
                  ALR= (BLR-(DX2-DX3)/(EX2-EX3))/(EX1-EX3)
                  BLR= BLR-ALR*(EX1+EX2)
                  CLR= DX1-(ALR*EX1+BLR)*EX1
                ENDIF
              WRITE(6,622) X2,VLIM,ALR,NCN,BLR,NCN2,CLR,NCN4
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  EX2= 1.d0/XX(I)**2
                  YY(I)= VLIM-(ALR+EX2*(BLR+EX2*CLR))/XX(I)**NCN
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.GE.4) THEN
c** For ILR.ge.4,   Y = VLIM-SUM(BB(K)/X**K) , (K=NCN,NMX=NCN+ILR-1)
          IF(LNPT.GT.0) THEN
              IF(NCN.LE.0) NCN= 1
              IMX1= ILR-1
              NMX= NCN+IMX1
              JR2= 0
              IF(IR2.GT.0) JR2= 2
              IDER= 0
              JMAX= ILR
              IF(IFXCN.GT.0) JMAX= IMX1
              WRITE(6,624) X2,ILR,NCN,VLIM
              IF(IFXCN.GT.0) WRITE(6,626) NCN,CNN
              ENDIF
c** Actually extrapolate with polynomial fitted to the last JMAX
c  values of  (VLIM - YI(I))*XI(I)**NMX  , & then convert back to  YY(I).
          IF(MBEG.LE.NPP) THEN
              J= NTP- JMAX
              DO  I= 1,JMAX
                  J= J+1
                  XJ(I)= XI(J)
                  YJ(I)= (VLIM-YI(J)/XI(J)**JR2)*XI(J)**NMX
                  IF(IFXCN.GT.0) YJ(I)= YJ(I)- CNN*XI(J)**IMX1
                  ENDDO
              DO  I= MBEG,NPP
                  CALL PLYINTRP(XJ,YJ,JMAX,XX(I),DUMM,JMAX,IDER)
                  YY(I)= DUMM(1)
                  IF(IFXCN.GT.0) YY(I)= YY(I)+ CNN*XX(I)**IMX1
                  YY(I)= VLIM-YY(I)/XX(I)**NMX
                  ENDDO
              ENDIF
          ENDIF
c** Finished extrapolation section.
   90 CONTINUE
c** Test smoothness at outer join to analytic extrapolation function
c     IF((LNPT.GT.0).AND.(MBEG.LE.NPP)) THEN
c         MST= MBEG-5
c         IF(MST.LT.1) MST= 1
c         MFN= MST+8
c         IF(MFN.GT.NPP) MFN= NPP
c         WRITE(6,627) X2,((XX(I),YY(I),I= J,MFN,3),J= MST,MST+2)
c         NPRF= MFN
c         ENDIF
c 627 FORMAT('     Verify smoothness of outer join at   X=',F9.5/
c    1  (3X,3(F10.5,G15.7)))
      RETURN
  612 FORMAT('  *** BUT *** since exponent has positive coefficient, swi
     1tch form ...')
  614 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1  F12.4,' - (',1PD13.6,') * exp{-',0PF10.6,' * (R -',F9.6,')**2}')
  616 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1 F12.4,' - (',1PD13.6,') * R**',0PF10.6,'  * exp{-(',F11.6,'*R)}')
  618 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,SP,1PD15.6,'/X**(',SS,D13.6,')] ,  yielding   NCN=',I3)
  620 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',SS,I1,']')
  622 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/
     1  '   Y=',F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',
     2  SS,I1,SP,D14.6,'/X**',SS,I2,']')
  624 FORMAT(' Function for  X .GE.',F7.3,'  generated by',I3,
     1 '-point inverse-power interpolation'/'   with leading term  1/R**
     2',I1,'  relative to dissociation limit   YLIM=',F11.3)
  626 FORMAT('   and (dimensionless) leading coefficient fixed as   C',
     1  I1,'=',G15.8)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PLYINTRP(XI,YI,NPT,RR,C,NCFT,IDER)
c* From the NPT known mesh points (XI,YI) ,given in order of increasing
c  or decreasing XI(I), select the NCFT points (XJ,YJ) surrounding the
c  given point RR, and by fitting an (NCFT-1)-th degree polynomial through
c  them, interpolate to find the function CC(1) and its first IDER
c  derivatives (CC(I+1),I=1,IDER) evaluated at RR.
c* Adapted by  R.J. Le Roy  from algorithm #416,Comm.A.C.M.;  27/02/1988
c=======================================================================
      INTEGER  I,J,K,I1,I2,IFC,IM,IDER,J1,NH,NPT,NCFT
      REAL*8  RR,XX,XI(NPT),YI(NPT),C(NCFT),XJ(20),YJ(20)
c
      IF((NCFT.GT.20).OR.(NCFT.GT.NPT)) GO TO 101
      NH= NCFT/2
c** First locate the known mesh points (XJ,YJ) bracketing RR
      I1= 1
      I2= NCFT
      IF(NCFT.NE.NPT) THEN
          IF(XI(NPT).LE.XI(1)) THEN
              DO  I= 1,NPT
                  IM= I
                  IF(XI(I).LT.RR) GO TO 20
                  ENDDO
            ELSE
              DO  I= 1,NPT
                  IM= I
                  IF(XI(I).GT.RR) GO TO 20
                  ENDDO
            ENDIF
   20     I1= IM-NH
          IF(I1.LE.0) I1= 1
          I2= I1+NCFT-1
          IF(I2.GT.NPT) THEN
              I2= NPT
              I1= I2-NCFT+1
              ENDIF
          ENDIF
      J= 0
      DO  I= I1,I2
          J= J+1
          XJ(J)= XI(I)-RR
          YJ(J)= YI(I)
          ENDDO
c** Now determine polynomial coefficients C(I).
      DO  I= 2,NCFT
          I1= I-1
          K= I1+1
          DO  J= 1,I1
              K= K-1
              YJ(K)= (YJ(K+1)-YJ(K))/(XJ(I)-XJ(K))
              ENDDO
          ENDDO
      C(1)= YJ(1)
      DO  I= 2,NCFT
          XX= XJ(I)
          C(I)= C(I-1)
          IF(I.NE.2) THEN
              I1= I-1
              K= I1+1
              DO  J= 2,I1
                  K= K-1
                  C(K)= -XX*C(K)+C(K-1)
                  ENDDO
              ENDIF
          C(1)= YJ(I)-XX*C(1)
          ENDDO
c** Finally, convert polynomial coefficients to derivatives at RR.
      IFC= 1
      IF(IDER.GE.NCFT) IDER= NCFT-1
      IF(IDER.LE.1) GO TO 99
      DO  I= 2,IDER
          J= I+1
          IFC= IFC*I
          C(J)= C(J)*IFC
          ENDDO
      IF(J.LT.NCFT) THEN
          J1= J+1
          DO  I= J1,NCFT
              C(I)= 0.D+0
              ENDDO
          ENDIF
   99 RETURN
  101 WRITE(6,601) NCFT,NCFT,NPT
      STOP
  601 FORMAT(/' *** Dimensioning ERROR in PLYINTRP :  either   (NCFT=',
     1  I2,' .GT. 20)   or   (NCFT=',I2,' .GT. NPT=',I3,')')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************
      SUBROUTINE SPLINT(LNPT,NTP,R1,V1,MBEG,MEND,XX,YY)
c** Subroutine to generate (if LNPT.ge.0) 4*NTP coefficients CSP(J)
c  of a cubic spline passing through the NTP points (R1(J),V1(J))
c  and to then calculate values of the resulting function YY(I) at the
c  entering abscissae values XX(I) for  I=MBEG to MEND.
c** If LNPT < 0 , generate function values at the given XX(I) using
c  the coefficients CSP(J) obtained and SAVEd on a preceding call.
c** Assumes both R1(J) & XX(I) are monotonic increasing.
c+++++ Calls only subroutine SPLINE +++++++++++++++++++++++++++++++++++
c======================================================================
      INTEGER MAXSP
      PARAMETER (MAXSP=6400)
      INTEGER  I,IER,I1ST,IDER,JK,K,KK,LNPT,N2,N3,NIPT,NTP,MBEG,MEND
      REAL*8 EPS,R2,RI,RRR,TTMP,R1(NTP),V1(NTP),CSP(MAXSP),
     1  YY(MEND),XX(MEND)
      SAVE CSP
c
      IF(4*NTP.GT.MAXSP) THEN
          WRITE(6,602) MAXSP,NTP
          STOP
          ENDIF
      EPS= 1.D-6*(R1(2)-R1(1))
      N2= 2*NTP
      N3= 3*NTP
      IF(LNPT.GT.0) THEN
c** On first pass for a given data set, generate spline function
c  coefficients in subroutine SPLINE
c** Start by using a cubic polynomial at each end of the range to get
c  the first derivative at each end for use in defining the spline.
          IDER= 1
          NIPT= 4
          I1ST= NTP-3
          CALL PLYINTRP(R1(I1ST),V1(I1ST),NIPT,R1(NTP),CSP,NIPT,IDER)
          TTMP= CSP(2)
          CALL PLYINTRP(R1,V1,NIPT,R1(1),CSP,NIPT,IDER)
          CSP(1)= CSP(2)
          CSP(2)= TTMP
c** Now call routine to actually generate spline coefficients
          CALL SPLINE(R1,V1,NTP,3,CSP,MAXSP,IER)
          IF(IER .NE. 0) THEN
              WRITE(6,604)
              STOP
              ENDIF
          ENDIF
      IF(MEND.LT.MBEG) GO TO 99
c** Now, use spline to generate function at desired points XX(I)
      DO  I= MBEG,MEND
          RI= XX(I)
          RRR= RI-EPS
          KK= 1
c** For a monotonic increasing distance array XX(I),  this statement
c  speeds up the search for which set of cubic coefficients to use.
          IF(I.GT.MBEG) THEN
              IF(XX(I).GT.XX(I-1)) KK= JK
              ENDIF
          DO  K= KK,NTP
              JK= K
              IF(R1(K).GE.RRR) GO TO 64
              ENDDO
   64     CONTINUE
          JK= JK-1
          IF(JK.LT.1) JK= 1
          R2= RI-R1(JK)
          YY(I)= CSP(JK)+R2*(CSP(NTP+JK)+R2*(CSP(N2+JK)+R2*CSP(N3+JK)))
          ENDDO
   99 RETURN
  602 FORMAT(' *** ERROR in SPLINT ***  Array dimension  MAXSP=',I4,
     1  ' cannot contain spline coefficients for  NTP=',I4)
  604 FORMAT(' *** ERROR in generating spline coefficients in SPLINE')
      END
c**********************************************************************
      SUBROUTINE SPLINE(X,Y,N,IOPT,C,N4,IER)
c** Subroutine for generating cubic spline coefficients
c  C(J), (J=1,N4=4*N) through the N points X(I), Y(I).
c** C(I+M*N), M=0-3  are the coefficients of order  0-3  of cubic
c  polynomial expanded about X(I) so as to describe the interval:
c             -  X(I) to X(I+1)  , if  X(I)  in increasing order
c             -  X(I-1) to X(I)  , if  X(I)  in decreasing order.
c** IOPT indicates boundary conditions used in creating the  spline .
c*  If (IOPT=0)  second derivatives = zero at both ends of range.
c*  If (IOPT=1)  1st derivative at first point X(1) fixed at C(1),
c                and 2nd derivative at X(N) = zero.
c*  If (IOPT=2)  1st derivative at last point X(N) fixed at C(2),
c                and 2nd derivative at X(1) = zero.
c*  If (IOPT=3)  constrain first derivatives at end points to have
c                (read in) values  C(1)  at  X(1)  &  C(2)  at  X(N)
c** IER is the error flag.  IER=0  on return if routine successful.
c-----------------------------------------------------------------------
      INTEGER I,II,IER,IOH,IOL,IOPT,J,J1,J2,J3,NER,N,N4,JMP
      REAL*8  A,H,R,DY2,DYA,DYB,XB,XC,YA,YB, X(N),Y(N),C(N4)
c
      JMP= 1
      NER= 1000
      IF(N.LE.1) GO TO 250
c** Initialization
      XC= X(1)
      YB= Y(1)
      H= 0.D0
      A= 0.D0
      R= 0.D0
      DYB= 0.D0
      NER= 2000
c
c  IOL=0 - given derivative at firstpoint
c  IOH=0 - given derivative at last point
c
      IOL= IOPT-1
      IOH= IOPT-2
      IF(IOH.EQ.1) THEN
          IOL= 0
          IOH= 0
          ENDIF
      DY2= C(2)
c
c  Form the system of linear equations
c  and eliminate subsequentially
c
      J= 1
      DO  I= 1,N
          J2= N+I
          J3= J2+N
          A= H*(2.D0-A)
          DYA= DYB+H*R
          IF(I.GE.N) THEN
c
c  set derivative dy2 at last point
c
              DYB= DY2
              H= 0.D0
              IF(IOH.EQ.0) GOTO 200
              DYB= DYA
              GOTO 220
              ENDIF
          J= J+JMP
          XB= XC
          XC= X(J)
          H= XC-XB
c
c  II= 0 - increasing abscissae
c  II= 1 - decreasing abscissae
c
          II= 0
          IF(H.LT.0) II= 1
          IF(H.EQ.0) GO TO 250
          YA= YB
          YB= Y(J)
          DYB= (YB-YA)/H
          IF(I.LE.1) THEN
              J1= II
              IF(IOL.NE.0) GO TO 220
              DYA= C(1)
              ENDIF
200       IF(J1.NE.II) GO TO 250
          A= 1.D0/(H+H+A)
220       R= A*(DYB-DYA)
          C(J3)= R
          A= H*A
          C(J2)= A
          C(I)= DYB
          ENDDO
c
c  back substitution of the system of linear equations
c     and computation of the other coefficients
c
      A= 1.D0
      J1= J3+N+II-II*N
      I= N
      DO  IOL= 1,N
          XB= X(J)
          H= XC-XB
          XC= XB
          A= A+H
          YB= R
          R= C(J3)-R*C(J2)
          YA= R+R
          C(J3)= YA+R
          C(J2)= C(I)-H*(YA+YB)
          C(J1)= (YB-R)/A
          C(I)= Y(J)
          A= 0.D0
          J= J-JMP
          I= I-1
          J2= J2-1
          J3= J3-1
          J1= J3+N+II
          ENDDO
      IER= 0
      RETURN
250   IER= NER
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,
     1                                                        NCN,CNN)
c** Generate analytic potential  VV(i)  as specified by the choice
c  of parameter IPOTL (see comments in PREPOT (& in main program))
c** All potentials generated in units cm-1 with absolute asymptote at
c  (input) energy VLIM for distance array  X0(i) Angstroms.
c** Return with NCN equal to power of asymptotically dominant inverse
c  power term in long range part of potential
c** Born-Oppenheimer correction functions in IPOTL=3 option may have up
c  to NBOB+1 terms.
c-----------------------------------------------------------------------
      INTEGER NBOB
      PARAMETER (NBOB=20)
      INTEGER  I,J,M,IBOB,IAN1,IAN2,IMN1,IMN2,RMN1,RMN2,IORD,IPOTL,
     1  NC1,NC2,NG1,NG2,NCMAX,NPAR,MPAR,NVARB,NPP,LNPT,NCN,GNS,GEL
      CHARACTER*2 NAME1,NAME2
      REAL*8  A0,A1,A2,A3,ALFA,BETA,BINF,B1,B2,CSAV,
     1  ABUND,CNN,DSCM,DX,DX1,FCT,
     2  FC1,FC2,MASS1,MASS2,RMASS1,RMASS2,RC6,RC8,RC10,RC12,RC14,
     3  RCNPAR,RD,RDIF,REQ,RPOW,RX,SC1,SC2,SG1,SG2,VLIM,VMIN,
     4  XDF,X1,XM,XN,XM2C,XP1,ZZ,ZP, CA1(0:NBOB),CA2(0:NBOB),
     5  GA1(0:NBOB),GA2(0:NBOB),PARM(20),XO(NPP),VV(NPP),RM2(NPP)
c
      SAVE IBOB,IPOTL,IORD,MPAR,NPAR,NVARB,DSCM,REQ,PARM,CA1,CA2,GA1,
     1     GA2,RX,CSAV
c
      IF(LNPT.GT.0) THEN
c** Parameter definitions listed preceeding CALL in subroutine PREPOT
c-----------------------------------------------------------------------
          READ(5,*) IPOTL, MPAR, NPAR, NVARB, IBOB, DSCM, REQ
          IF(IPOTL.EQ.1) NVARB= 0
          IF((IPOTL.EQ.3).AND.(NPAR.EQ.-1)) NVARB= 2
          IF(NVARB.GT.0)  READ(5,*) (PARM(I), I=1,NVARB)
          IF(IBOB.GT.0) THEN
              READ(5,*) RMN1, RMN2, NC1, NC2, NG1, NG2, RX
c-----------------------------------------------------------------------
              NCMAX= MAX0(NC1,NC2,NG1,NG2)
              IF(NCMAX.LT.0) THEN
                  IBOB= 0
                ELSE
c** If appropriate, read parameters & prepare to add mass-dep. BOB corrn
                  CALL MASSES(IAN1,IMN1,NAME1,GEL,GNS,MASS1,ABUND)
                  CALL MASSES(IAN1,RMN1,NAME1,GEL,GNS,RMASS1,ABUND)
                  CALL MASSES(IAN2,IMN2,NAME2,GEL,GNS,MASS2,ABUND)
                  CALL MASSES(IAN2,RMN2,NAME2,GEL,GNS,RMASS2,ABUND)
                  WRITE(6,628)
c  For simplicity, first zero out all correction function coefficients
                  DO  I=0,NCMAX
                      CA1(I)= 0.d0
                      CA2(I)= 0.d0
                      GA1(I)= 0.d0
                      GA2(I)= 0.d0
                      ENDDO
                  FC1= 0.d0
                  FC2= 0.d0
c=======================================================================
c** Read actual B-O-B polynomial expansion coefficients
c=======================================================================
                  IF(NC1.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (CA1(I), I=0,NC1)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,630) 1,MASS1,NC1,NAME1,RMN1,NAME1,
     1                                    IMN1,NC1+1,(CA1(I),I= 0,NC1)
                          FC1= 1.d0 - RMASS1/MASS1
                        ELSE
                          WRITE(6,632) 1,MASS1,NC1,NAME1,IMN1,NC1+1,
     1                                               (CA1(I),I= 0,NC1)
                          FC1= 1.d0/MASS1
                        ENDIF
                      ENDIF
                  IF(NC2.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (CA2(I), I=0,NC2)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,630) 2,MASS2,NC2,NAME2,RMN2,NAME2,
     1                                    IMN2,NC2+1,(CA2(I),I= 0,NC2)
                          FC2= 1.d0 - RMASS2/MASS2
                        ELSE
                          WRITE(6,632) 2,MASS2,NC2,NAME2,IMN2,NC2+1,
     1                                                (CA2(I),I=0,NC2)
                          FC2= 1.d0/MASS2
                        ENDIF
                      ENDIF
                  IF(NG1.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (GA1(I), I=0,NG1)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,634) 1,MASS1,NG1,NAME1,RMN1,NAME1,
     1                                    IMN1,NG1+1,(GA1(I),I= 0,NG1)
                        ELSE
                          WRITE(6,636) 1,MASS1,NG1,NAME1,IMN1,RX,
     1                                         NG1+1,(GA1(I),I= 0,NG1)
                        ENDIF                     
                      ENDIF
                  IF(NG2.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (GA2(I), I=0,NG2)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,634) 2,MASS2,NG2,NAME2,RMN2,NAME2,
     1                                    IMN2,NG2+1,(GA2(I),I= 0,NG2)
                        ELSE
                          WRITE(6,636) 2,MASS2,NG2,NAME2,IMN2,RX,
     1                                         NG2+1,(GA2(I),I= 0,NG2)
                        ENDIF
                      ENDIF
                  DO  I=0,NCMAX
                      CA1(I)= CA1(I)*FC1
                      CA2(I)= CA2(I)*FC2
                      IF(RX.LE.0.d0) THEN
                          GA1(I)= GA1(I)*(1.d0-FC1)
                          GA2(I)= GA2(I)*(1.d0-FC2)
                        ELSE
                          GA1(I)= GA1(I)*FC1
                          GA2(I)= GA2(I)*FC2
                        ENDIF
                      ENDDO
                ENDIF
              ENDIF
          ENDIF
      IF(IPOTL.EQ.1) THEN
c=======================================================================
c** Generate a  Lennard-Jones(MPAR,NPAR)  potential here.
c=======================================================================
          XM= MPAR
          XN= NPAR
          XDF= DSCM/(XM-XN)
          IF(LNPT.GE.0) WRITE(6,600) MPAR,NPAR,DSCM,REQ
          NCN= NPAR
          CNN= XM*XDF*REQ**NPAR
          DO  I= 1,NPP
              VV(I)= (XN*(REQ/XO(I))**MPAR - XM*(REQ/XO(I))**NPAR)*XDF
     1                  +VLIM
              ENDDO
          ENDIF
      IF(IPOTL.EQ.2) THEN
c=======================================================================
c** Generate an MLJ potential [as per JCP 112, 3949 (2000)] here ...
c=======================================================================
          NCN= NPAR
          IORD= NVARB-1
          IF(MPAR.LE.0) THEN
c  If appropriate, prepare to calculate switching function
              IORD= IORD- 3
              CNN= PARM(NVARB-2)
              BINF= DLOG(2.d0*DSCM*REQ**NPAR/CNN)
              ALFA= PARM(NVARB- 1)
              RX= PARM(NVARB)
            ELSE
c  Generate limiting Cn value for non-switching case from BINF
              BINF= 0.D0
              DO  I= 1,NVARB
                  BINF= BINF+ PARM(I)
                  ENDDO
              CNN= 2.d0*DSCM*REQ**NPAR *DEXP(-BINF)
            ENDIF
          IF(LNPT.GT.0) THEN
              WRITE(6,602) DSCM,REQ
              IF(MPAR.GT.0) WRITE(6,603) IORD,MPAR,MPAR,MPAR,MPAR,
     1                                    IORD+1,(PARM(J),J= 1,IORD+1)
              IF(MPAR.LE.0) WRITE(6,603) IORD,1,1,1,1,
     1                                    IORD+1,(PARM(J),J= 1,IORD+1)
              IF(MPAR.LE.0) WRITE(6,604) NPAR,NPAR,NPAR,PARM(NVARB-2),
     1                                       PARM(NVARB-1),PARM(NVARB)
              ENDIF
          NCN= NPAR
c  Loop over distance array XO(I)
          DO  I= 1,NPP
              ZZ= (XO(I)- REQ)/(XO(I)+ REQ)
              IF(MPAR.GT.1) ZZ= (XO(i)**MPAR - REQ**MPAR)/
     1                          (XO(i)**MPAR  + REQ**MPAR)
              BETA= 0.d0
              DO  J= IORD,0,-1
                  BETA= BETA*ZZ+ PARM(J+1)
                  ENDDO
c  Calculate and apply switching function to MLJ exponent coefficient
              IF(MPAR.LE.0) BETA= BINF+ (BETA- BINF)/
     1                                 (1.d0+ DEXP(ALFA*(XO(I) - RX)))
              VV(I)= DSCM*(1.d0 - (REQ/XO(I))**NPAR *DEXP(-BETA*ZZ))**2
     1                                                    - DSCM+ VLIM
              ENDDO
          ENDIF
      IF(IPOTL.EQ.3) THEN
c=======================================================================
c** Generate a simple Morse, or Extended (EMO) Morse potential, or as
c  special cases, Coxon's GMO or Wei Hua's generalized Morse
c=======================================================================
          BETA= PARM(1)
          NCN= 99
          IF(LNPT.GE.0) THEN
              IF(MPAR.EQ.-1) THEN
c** Option to generate Wei Hua's extended 4-parameter Morse-type potl.
                  CSAV= PARM(2)
                  WRITE(6,605) DSCM,REQ,CSAV,BETA
                ELSE
                  IF(NVARB.LE.1) WRITE(6,606) DSCM,REQ,BETA
                  IF(NVARB.GT.1) THEN
                      IF(MPAR.GT.0) WRITE(6,608) DSCM,REQ,NVARB-1,
     1                  MPAR,MPAR,MPAR,MPAR,NVARB,(PARM(i),i= 1,NVARB)
                      IF(MPAR.EQ.-2) WRITE(6,610) DSCM,REQ,NVARB-1,
     1                                            (PARM(i),i= 1,NVARB)
                      ENDIF
                ENDIF
              ENDIF
c  Loop over distance array XO(I)
          DO  I= 1,NPP
c ... for Wei Hua's extended Morse function ...
              IF(MPAR.EQ.-1) THEN
                  VV(I)= DSCM*((1.d0 - DEXP(-BETA*(XO(I)-REQ)))/(1.d0
     1                - CSAV*DEXP(-BETA*(XO(I)-REQ))))**2 - DSCM+ VLIM
                ELSE
                  IF(NVARB.GT.1) THEN
                      ZZ= (XO(I)- REQ)/(XO(I)+ REQ)
                      IF(MPAR.GT.1) ZZ= (XO(i)**MPAR - REQ**MPAR)/
     1                                  (XO(i)**MPAR + REQ**MPAR)
c ... for Coxon-Hajigeorgiou "GMO" potential
                      IF(MPAR.EQ.-2) ZZ= (XO(I)- REQ)
                      BETA= 0.d0
                      DO  J= NVARB,1,-1             
                          BETA= BETA*ZZ+ PARM(J)
                          ENDDO
                      ENDIF
                  VV(I)=  DSCM*(1.d0 - DEXP(-BETA*(XO(I)-REQ)))**2
     1                                                    - DSCM+ VLIM
                ENDIF
              ENDDO
          ENDIF
      IF(IPOTL.EQ.4) THEN
c=======================================================================
c** Generate Seto-modified form of Surkus' GPEF function which includes
c  Dunham, SPF and OT forms as special cases.
c=======================================================================
          VMIN= VLIM
          VLIM= 1.d9
          A0= DSCM
          IORD= NVARB-2
          X1= 1.d0
          FCT= PARM(NVARB-1)
          IF((NPAR.NE.0).AND.(DABS(FCT).GT.0.d0)) THEN
              FCT= 1.d0/PARM(NVARB-1)
              DO  J=1,IORD
                  X1= X1+ PARM(J)*FCT**J
                  ENDDO
              DSCM= DSCM*X1*FCT**2 + VMIN
              ENDIF
          IF(NPAR.EQ.1) THEN
c  Cases with power =1 (including Dunham, SPF & O-T expansions).
              IF(DABS(PARM(NVARB-1)).LE.0.d0) THEN
c ... print for Dunham expansion ...
                  WRITE(6,612) PARM(NVARB),REQ,VMIN,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= -99
                  CNN= 0.d0
                  ENDIF
              IF(DABS(PARM(NVARB)).LE.0.d0) THEN
c ... print for Simons-Parr-Finlan expansion ...
                  WRITE(6,614) PARM(NVARB-1),REQ,DSCM,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= 1
                  ENDIF
              IF(DABS(PARM(NVARB)-PARM(NVARB-1)).LE.0.d0) THEN
c ... print for Ogilvie-Tipping expansion ...
                  WRITE(6,616) PARM(NVARB),REQ,DSCM,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= 1
                  ENDIF
              ENDIF
          IF((NPAR.NE.0).AND.((NPAR.NE.1).OR.
     1                ((DABS(PARM(NVARB)-PARM(NVARB-1)).GT.0.d0).AND.
     2               (DABS(PARM(NVARB)*PARM(NVARB-1)).GT.0.d0)))) THEN
c ... print for general GPEF expansion variable ...
              IF(NPAR.LT.0) THEN
c ... for negative NPAR, convert to equivalent positive NPAR case
                  NPAR= -NPAR
                  A1= PARM(NVARB)
                  PARM(NVARB)= -PARM(NVARB-1)
                  PARM(NVARB-1)= -A1
                  ENDIF
              WRITE(6,618) NPAR,NPAR,PARM(NVARB-1),NPAR,PARM(NVARB),
     1                 NPAR,REQ,DSCM,A0,NVARB-2,(PARM(I),I= 1,NVARB-2)
              NCN= NPAR
              ENDIF
          IF(NPAR.EQ.0) THEN
c** For case of simple power series in  R  itself
              WRITE(6,620) NVARB,(PARM(I),I= 1,NVARB)
              DO  I= 1, NPP
                  ZP= 1.d0
                  A1= 0.d0
                  DO  J= 1,NVARB
                      A1= A1+ PARM(J)*ZP
                      ZP= ZP*XO(I)
                      ENDDO
                  VV(I)= A1+ VMIN
                  ENDDO
              VLIM= VV(NPP)
              RETURN
              ENDIF
c ... otherwise - generate potential as a GPEF-type expansion
          DO  I= 1, NPP
              ZZ= (XO(I)**NPAR - REQ**NPAR)/(PARM(NVARB-1)*XO(I)**NPAR
     1                                       + PARM(NVARB)*REQ**NPAR)
              A1= 1.d0
              ZP= 1.d0
              DO  J=1, NVARB-2
                  ZP= ZP*ZZ
                  A1= A1+ PARM(J)*ZP
                  ENDDO
              VV(I)= A0*ZZ*ZZ*A1 + VMIN
              ENDDO
          ENDIF
      IF(IPOTL.EQ.5) THEN
c=======================================================================
c** For generalized  H.F.D.(NPAR,6,8,10,12,14)  potential with reduced
c  form   VBAR = ALFA*x**PARM(5) * exp[-BETR*x - PARM(4)*x**2] - D(x)*
c       [PARM(6)/x**NPAR + PARM(7)/x**6 + PARM(8)/x**8 + PARM(9)/x**10
c       + PARM(10)/X**12 + PARM(11)/X**14]   where   x=r/R_e , 
c  VBAR= V/epsilon   and   D(x)= exp[-PARM(1)*(PARM(2)/x - 1)**PARM(3)]
c  for  x < PARM(2)
c=======================================================================
          A1= PARM(1)
          A2= PARM(2)
          A3= PARM(3)
          B2= PARM(4)
          RC8= 0.d0
          RC10= 0.d0
          RC12= 0.d0
          RC14= 0.d0
          RCNPAR= PARM(6)
          NCN= 6
          IF(RCNPAR.GT.0.d0) NCN= NPAR
          RC6= PARM(7)
          IF(NVARB.ge.8)  RC8= PARM(8)
          IF(NVARB.ge.9)  RC10= PARM(9)
          IF(NVARB.ge.10) RC12= PARM(10)
          IF(NVARB.ge.11) RC14= PARM(11)
          DX= 1.d0
          DX1= 0.d0
          IF(A2.GT.1.d0) THEN
              DX= DEXP(-A1*(A2- 1.d0)**A3)
              DX1= A1*A2*A3*DX*(A2- 1.d0)**(A3- 1.d0)
              ENDIF
          ALFA= -1.D0+ (RCNPAR+ RC6+ RC8+ RC10+ RC12+ RC14)*DX
          IF(ALFA.LE.0.d0) THEN
              WRITE(6,622) RCNPAR,RC6,RC8,RC10,RC12,RC14,ALFA
              STOP
              ENDIF
          B1= ((NPAR*RCNPAR+6.D0*RC6+8.D0*RC8+10.D0*RC10+12.d0*RC12+
     1      14.d0*RC14)*DX - (RCNPAR+RC6+RC8+RC10+RC12+RC14)*DX1)/ALFA
     2      + PARM(5) - 2.D0*B2
          ALFA= ALFA*DEXP(B1+B2)
          IF(LNPT.GE.0) WRITE(6,624) NPAR,PARM(5),B1,B2,ALFA*DSCM,
     1                 RCNPAR,RC6,RC8,RC10,RC12,RC14,A1,A2,A3,DSCM,REQ
          DO  I= 1,NPP
              X1= XO(I)/REQ
              XP1= 0.0D0
              IF((B1*X1+ B2*X1**2).LT.170.D0) XP1= DEXP(-X1*(B1+ B2*X1))
              XP1= XP1*X1**PARM(5)
              FC1= 1.D0
              IF(A2.GT.X1) FC1= DEXP(-A1*(A2/X1- 1.d0)**A3)
              XM2C= (REQ/XO(I))**2
              VV(I)= DSCM*(ALFA*XP1- FC1*(((((RC14*XM2C+RC12)*XM2C+RC10)
     1         *XM2C+ RC8)*XM2C+ RC6)*XM2C**3 + RCNPAR/X1**NPAR)) + VLIM
              ENDDO
          ENDIF
      IF(IBOB.GT.0) THEN
c=======================================================================
c** If appropriate, generate Born-Oppenheimer breakdown correction
c  functions to rotationless and/or centrifugal potential(s).
c [Special "Coxon" option: if  RX > 0.0, expand as per older Coxon work]
c=======================================================================
          IF(RX.GE.0.D0) RDIF= REQ-RX
          DO  I=1,NPP
              IF(RX.LE.0.d0) THEN
                  ZZ= (XO(I)-REQ)/(XO(I)+REQ)
                ELSE
                  ZZ= XO(I)- REQ
                  RD= XO(I)- RX
                ENDIF
              SC1= 0.d0
              SC2= 0.d0
              SG1= 0.d0
              SG2= 0.d0
              RPOW= 1.d0
              DO  J= 0,NCMAX
                  SC1= SC1+ RPOW*CA1(J)
                  SC2= SC2+ RPOW*CA2(J)
                  IF(RX.LE.0.d0) THEN
                      SG1= SG1+ RPOW*GA1(J)
                      SG2= SG2+ RPOW*GA2(J)
                    ELSE
                      M= J-1
                      SG1= SG1+ (RD**J -RDIF**J)*GA1(J)
                      SG2= SG2+ (RD**J -RDIF**J)*GA2(J)
                    ENDIF
                  RPOW= RPOW*ZZ
                  ENDDO
              RM2(I)= (1.d0+ SG1+ SG2)/XO(i)**2
              VV(I)= VV(I) + SC1 + SC2
              ENDDO
          ENDIF
      RETURN
  600 FORMAT(/' Lennard-Jones(',I2,',',I2,') potential with   De=',
     1  F10.3,'(cm-1)   Re =',F10.6,'(A)')
  602 FORMAT(/' Use an MLJ potential with   De =',F10.3,
     1  '(cm-1)    Re =',F12.8,'(A)')
  603 FORMAT(3x,'with parameter BETA an order-',i2,' polynomial in   y =
     1 (R^',i1,' - Re^',i1,')/(R^',i1,' + Re^',i1,')'/'   with',i3,
     2  ' coefficients:',1PD16.8,2D16.8:/(8x,4D16.8:))
  604 FORMAT(' & exponent switching function yielding limiting C',i1,
     1 '/R^',i1,' with   C_',i1,'=',1PD13.6/10x,'defined by   ALPHA_s=',
     2  0Pf9.6,'   R_s=',f10.6)
  605 FORMAT(/' Potential is a Hua-Wei 4-parameter Morse type function w
     1ith   De =',F11.4/11x,'Re =',F12.9,'   C=',f7.4,'   &   beta=',
     1  F13.10,' [1/Angstroms]')
  606 FORMAT(/' Potential is a simple Morse function with   De =',F11.4,
     1  '    Re =',F12.9/39x,'and   beta =',F13.10,' [1/Angstroms]')
  608 FORMAT(/' Potential is Extended Morse Oscillator with   De=',
     1  F11.4,'    Re=',F12.9/3x,'Exponent factor "beta" is order-',i2,
     2 ' power series in  y=(R^',i1,' -Re^',i1,')/(R^',i1,' +Re^',i1,')'
     3 /'   with',I3,' coefficients:',1x,1PD18.9,2D18.9:/(7X,4D18.9:))
  610 FORMAT(/' Potential is Generalized Morse Oscillator with   De=',
     1 F10.3,'   Re=',F11.8/4x,'Exponent factor "beta" is',i3,' order po
     2wer series in (R-Re) with coefficients:'/4x,1PD18.9,3D18.9:/
     3 (4X,4D18.9:))
  612 FORMAT(/' Potential is a Dunham expansion in  (R-Re)/(',f5.2,
     1  ' * Re)  with   Re=',f12.9/'  V(Re)=',f12.4,'    a0=',1PD16.9,
     2  '   and',i3,'  a_i coefficients:'/(5D16.8))
  614 FORMAT(/' Potential is an SPF expansion in  (R-Re)/(',F5.2,
     1  '* R)  with   Re=',f12.9/5x,'De=',g18.10,'   b0=',
     2  1PD16.9,'   and',i3,'  b_i  coefficients:'/(5D16.8))
  616 FORMAT(/' Potential is an O-T expansion in  (R-Re)/[',f5.2,
     1  '*(R+Re)]  with   Re=',f12.9/5x,'De=',G18.10,
     2  '   c0=',1PD16.9,'   and',i3,'  c_i coefficients:'/(5D16.8))
  618 FORMAT(/' Potential is a general GPEF expansion in  (R^',i1,
     1  ' - Re^',i1,')/(',SP,F5.2,'*R^',SS,i1,SP,F6.2,'*Re^',SS,i1,')'/
     2  5x,'with   Re=',f12.9,'   De=',g18.10,'   g0=',1PD16.9/
     3  5x,'and',i3,'  g_i coefficients:  ',3D16.8/(5D16.8:))
  620 FORMAT(/' Potential is an',i3,'-term power series in  R  with coef
     1ficients (starting from power=0):'/(5D16.8))
  622 FORMAT(/' *** ERROR in generating HFD potential *** C',i1,
     1 ', C6, C8, C10, C12, C14  =',6G15.7/10X,'yield    ALFA =',G15.7)
  624 FORMAT(/' Potential is Generalized HFD(',i1,',6,8,10,12,14) with',
     1 '  gamma=',f9.6/'    beta1=',f12.8,'    beta2=',f9.6,'    A=',
     2 1PD16.9/"    reduced {Cn's}:",3D14.6/19x,3D14.6/'    Damping func
     3tion  D(R)= exp[ -',0Pf6.4,'*(',f7.4,'/X -1.0)**',f5.2,']' /
     4 '  & DSCM=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]')
  628 FORMAT(' ')
  630 FORMAT(' B-O-B correction to rotationless potential for atom-',
     1  I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in {(R-R
     2e)/(R+Re)}] * [1- MASS(',A2,i3,')/MASS(',A2,I3,')]'/5x,'with',i3,
     3  ' coefficients:',3G18.10:/(8x,4G18.10:))
  632 FORMAT(' B-O-B correction to rotationless potential for atom-',
     1 I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in (R-Re)
     2]/[MASS(',A2,I3,')]   with',i3,' coefficients:'/(5x,4G18.10:))
  634 FORMAT(' B-O-Breakdown correction to centrifugal term for atom-',
     1 I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in {(R-Re
     2/(R+Re)}] * [MASS(',A2,I3,')/MASS(',A2,I3,')]'/5x,'with',i3,' coef
     3ficients:',3G18.10:/(8x,4G18.10:))
  636 FORMAT(' B-O-Breakdown correction to centrifugal term for atom-',
     1 I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in {(R-Rx
     2)**i - (Re-Rx)**i}] / [MASS(',A2,I3,')]'/5x,'with  Rx=',f6.3,
     3 '  &',i3,' coefficients: ',1PD18.9,D18.9:/(5x,4D18.9:))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c** Version 0.9s dated Mar 7, 2000.
c***********************************************************************
      SUBROUTINE ALF(NDP,RMIN,RH,V,SWF,VLIM,KVMAX,AFLAG,ZMU,EPS,NCN,GV,
     1               INNR,IWR)
c***********************************************************************
c** The subroutine ALF (Automatic vibrational Level Finder) will
c   automatically generate the eigenvalues from the first vibrational
c   level (v=0) to a user specified level (v=KVMAX) or the highest
c   allowed vibrational level of a given smooth single (or double)
c   minimum potential (V). These energies are stored and returned to the
c   calling program in the molecular constants array GV(v=0-KVMAX).
c** For any errors that cannot be resolved within the subroutine, ALF
c   returns AFLAG with a value that defines which error had occured.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++   COPYRIGHT 1998 - 1999  by  Jenning Seto and Robert J. Le Roy   +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c     of it without the express written permission of the authors.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Please inform me of any bugs, by phone at: (519)888-4567, ext. 4099 +
c++++++++ by e-mail to: jyseto@uwaterloo.ca , or write me at: ++++++++++
c+++ Dept. of Chemistry, Univ. Waterloo, Waterloo, Ontario  N2L 3G1 ++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Based on the automatic level finding routine found in LEVEL 6.0
c   written by Robert J. Le Roy
c** Uses the Schrodinger solver subroutine SCHRQ.
c
c** On entry:
c    NDP    is the number of datapoints used for the potential.
c    RMIN   is the inner radial distance of the potential (ang).
c    RH     is the meshvalue (ang).
c           NDP, RMIN, and RH define the radial range over which the
c           potential is defined.
c    V(i)   is the scaled input potential (cm-1).
c           The scaling factor BFCT is (2*mu/hbar^2)*RH^2.
c    VLIM   is the potential asymptote (cm-1).
c    KVMAX  is the maximum vibrational level for which we wish to find.
c    AFLAG  is the rotational state of the potential.
c    ZMU    is the reduced mass of the diatom (amu).
c    EPS    is the energy convergence criterion (cm-1).
c    NCN    is the near dissociation limit radial exponential.
c    IWR    specifies the level of printing inside SCHRQ
c           <> 0 : print error & warning descriptions.
c           >= 1 : also print final eigenvalues & node count.
c           >= 2 : also show end-of-range wave function amplitudes.
c           >= 3 : print also intermediate trial eigenvalues, etc.
c
c** On exit:
c    KVMAX   returns the highest allowed vibrational quantum number if
c            less than the inputed KVMAX.
c    AFLAG   returns calculation outcome to calling program.
c            >=  0 : Subroutine has functioned normally.
c             = -1 : KVMAX larger than number of allowed levels.
c             = -2 : Initial trial energy is unusable.
c             = -3 : Calculated trial energy is unusable.
c             = -4 : Cannot find first vibrational level.
c             = -5 : Calculated trial energy too low.
c             = -6 : Calculated trial energy too high.
c             = -7 : An impossible situation occured.
c             = -8 : Potential found to have a second minimum.
c    GV(v)   contains the vibrational energy level spacings and
c            rotational constants in cm-1 for each level.
c    INNR(v) labels each level as belonging to the inner (INNR = 1) or
c            outer (INNR = 0) well.
c
c** Flags: Modify only when debugging.
c    AWO   specifies the level of printing inside ALF
c          <> 0 : print error & warning descriptions.
c          >  0 : also print intermediate ALF messages.
c    MCO   specifies the level of printing of molecular constants.
c          >  0 : print out vibrational energies to channel-21.
c    INNER specifies wave function matching (& initiation) conditions.
c          = 0 : Match inward & outward solutions at outermost wave
c                function maximum
c          <>0 : Match at inner edge of classically allowed region.
c          < 0 : uses zero slope inner boundary condition.
c          For most normal cases set INNER = 0,  but ......
c            To find "inner-well-dominated" solutions of an asymmetric
c            double minimum potential, set  INNER > 0.
c            To find symmetric eigenfunctions of a symmetric potential,
c            set INNER < 0  & start integration (set RMIN) at potential
c            mid point.
c    LPRWF specifies option of printing out generated wavefunction
c          > 0 : print wave function every LPRWF-th  point.
c          < 0 : compactly write to channel-7 every |LPRWF|-th wave
c                function value.
c          A lead "card" identifies the level, gives the position of
c          1-st point and radial mesh, & states No. of  points.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The dimensioning parameters must be consistant with the sizes of the
c   arrays used in the calling program.
c
c    NVIBMX  is the maximum number of vibrational levels considered.
c            Note: NVIBMX should be larger than KVMAX.
c
      INTEGER NVIBMX
      PARAMETER (NVIBMX = 400)
c
c** NF counts levels found in automatic search option
c
c** OWL holds the vibrational levels that are contained in the outer
c   well.
c** IWL holds the vibrational levels that are contained in the inner
c   well (if present).
c
      INTEGER NDP,KVMAX,NCN,KV,AFLAG,NF,NBEG,NEND,INNR(0:KVMAX),IWR,
     1  I,IZPE,IVDIF,IVCOR,IQT,IEG,LTRY,AWO,MCO,INNER,LPRWF,JROT,
     2  NPMIN, NPMAX, NIWL,IWL(0:NVIBMX),NOWL,OWL(0:NVIBMX)
c
      REAL*8 RMIN,RMAX,RH,V(NDP),SWF(NDP),VLIM,EO,ZMU,EPS,
     1  GV(0:KVMAX),BV(0:NVIBMX),BVDOUT,BVDIN,AO,VD,
     2  BZ,BFCT,PW,PWI,GAMA,VMIN,VMAX,RE,PMAX,VDMV,VDL,VDU,
     3  VPMIN(10), RPMIN(10), VPMAX(10), RPMAX(10),
     3  ZQ, ZH, Z1, Z2
c
      DATA AWO/0/,MCO/0/,LPRWF/0/
c
      DATA ZQ/0.25D0/,ZH/0.5D0/,Z1/1.D0/,Z2/2.D0/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Check if the dimensions are adequate.
c
      IF (KVMAX.GT.NVIBMX) THEN
        WRITE(6,610)
        WRITE(6,613) KVMAX, NVIBMX
        STOP
      ENDIF
c
c** Initialize level counters for each well.
c
      DO I = 0,KVMAX
        INNR(I) = 0
        IWL(I) = 0
        OWL(I) = 0
      END DO
c
c** Initialize the remaining variables and flags.
c
      NF = 0
      NIWL = 0
      NOWL = 0
      KV = 0
      INNER = 0
      LTRY = 0
      CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** Store rotational quantum number.
c
      JROT = AFLAG
c
c** Numerical factor  16.85762908  based on 1998 physical constants.
c
      BZ = ZMU/16.85762908d0
      BFCT = BZ*RH*RH
c
c** RMAX is the outer radial distance over which the potential is
c   defined.
c
      RMAX = RMIN + DBLE(NDP-1)*RH
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Locate the potential minima.
c
      NPMIN = 0
      DO I = 2,NDP-1
        IF ((V(I).LT.V(I-1)).AND.(V(I).LT.V(I+1))) THEN
          NPMIN = NPMIN + 1
          RPMIN(NPMIN) = RMIN + DBLE(I-1)*RH
          VPMIN(NPMIN) = V(I) / BFCT
          IF (NPMIN.EQ.10) GOTO 100
        END IF
      END DO
c
c** If a minimum cannot be found, then print a warning and exit.
c
  100 IF (NPMIN.EQ.0) THEN
        WRITE(6,610)
        WRITE(6,614)
        STOP
      END IF
c
c** If more than two minima are found, then print a warning and exit.
c
      IF (NPMIN.GT.2) THEN
        WRITE(6,605)
        WRITE(6,615) NPMIN, 'minima'
*        STOP
      END IF
c
c** Locate the potential maxima (if it exists).
c
      NPMAX = 0
      DO I = 2,NDP-1
        IF ((V(I).GT.V(I-1)).AND.(V(I).GT.V(I+1))) THEN
          NPMAX = NPMAX + 1
          RPMAX(NPMAX) = RMIN + DBLE(I-1)*RH
          VPMAX(NPMAX) = V(I) / BFCT
          IF (NPMAX.EQ.10) GOTO 150
        END IF
      END DO
c
c** If no maxima were found, then set the maximum to be the value at the
c   end of the range.
c
  150 IF (NPMAX.EQ.0) THEN
        NPMAX = 1
        RPMAX(NPMAX) = RMAX
        VPMAX(NPMAX) = V(NDP) / BFCT
      END IF
c
c** If more than three maxima are found, then print a warning and exit.
c
      IF (NPMAX.GT.3) THEN
        WRITE(6,605)
        WRITE(6,615) NPMAX, 'maxima'
*        STOP
      END IF
c
c** If there is no rotationless barrier to assiciation, then set the
c   final VPMAX to be the value at the end of the range.
c
      IF (RPMAX(NPMAX).LT.RPMIN(NPMIN)) THEN
        NPMAX = NPMAX + 1
        RPMAX(NPMAX) = RMAX
        VPMAX(NPMAX) = V(NDP) / BFCT
      END IF
c
c** If a maxima occurs before a minima, then the potential turns over in
c   short range region and should not be used. Print a warning and exit.
c
      IF (RPMAX(1).LT.RPMIN(1)) THEN
        WRITE(6,610)
        WRITE(6,616) RPMAX(1)
*        STOP
      END IF
c
c** Now find the absolute potential minimum.
c
      VMIN = VPMIN(1)
      RE = RPMIN(1)
      DO I = 2,NPMIN
        IF (VMIN.GT.VPMIN(I)) THEN
          VMIN = VPMIN(I)
          RE = RPMIN(I)
        END IF
      END DO
c
c** Now find the absolute potential minimum.
c
      VMAX = VPMAX(1)
      DO I = 2,NPMAX
        IF (VMAX.LT.VPMAX(I)) VMAX = VPMAX(I)
      END DO
c
c** If the absolute potential maximum is lower than the absolute
c   potential minimum, then print out an error statement and quit.
c
      IF (VMAX.LE.VMIN) THEN
        WRITE(6,610)
        WRITE(6,617)
        STOP
      END IF
c
c** Otherwise, print out the results if desired.
c
      IF (AWO.GT.0) THEN
        WRITE(6,650) NPMIN, VMIN
        WRITE(6,651) NPMAX, VMAX
      END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Calculate 2*NCN/(NCN - 2) for use when calculating trial energies.
c
      PW = 20.0d0
      IF ((NCN.GT.0).AND.(NCN.NE.2)) PW = Z2*DBLE(NCN)/(DBLE(NCN)-Z2)
      IF (VMAX.GT.VLIM) PW = Z2
      PWI = Z1/PW
c
c** Use Lennard-Jones estimation of zero point energy to determine the
c   initial trial energy.
c                            _____________________________
c            vD + 0.5 = ao \/ZMU * De * Re^2 / 16.85762908
c
c                 De = A (vD - v)^3 = A (vD + 0.5)^3
c
c              E(v=0) = VMIN + A [(vD + 0.5)^3 - vD^3]
c
c** Choose AO to have a value of 0.25.
c
      AO = ZQ
      VD = AO * DSQRT(BZ*(VMAX-VMIN)) * RE - ZH
      AO = (VMAX-VMIN)*(Z1 - (VD/(VD+ZH))**3)
      EO = VMIN + AO
c
c** If desired, write out energy level information.
c
      IF (MCO.GE.1) THEN
         WRITE(21,2100)
         WRITE(21,2110) RMIN, RMAX, RH, BZ, ZMU
         WRITE(21,2111) EPS
         WRITE(21,2112)
         WRITE(21,2101)
      END IF
c=========== Begin Actual Eigenvalue Calculation Loop Here =============
c** Compute eigenvalues ... etc up to the KVMAXth vibrational level.
c** When attempts to find the next eigenvalue fails, then perhaps the
c   next level is located in a second (inner) well. If so, then the
c   subroutine will set INNER = 1, and attempt to find that level.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine SCHRQ to find eigenvalue EO and eigenfunction SWF(I).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   10 IF (AWO.GT.0) THEN
        WRITE(6,601)
        IF (INNER.EQ.0) WRITE(6,602)
        IF (INNER.EQ.1) WRITE(6,603)
      END IF
      CALL SCHRQ(KV,JROT,EO,GAMA,PMAX,VLIM,V,SWF,BFCT,EPS,RMIN,RH,NDP,
     1  NBEG,NEND,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The SCHRQ error condition is KV < 0.
c   There are three possible situations to consider:
c     EO > VMAX : Trial energy greater than potential maximum
c     NF = 0 : Looking for the first vibrational level (v = 0)
c     NF > 0 : Looking for the other vibrational levels (v > 0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the case when the next trial energy is higher than the potential
c   maximum, try one last ditch attempt to find the highest bound level
c   (quasi or otherwise) in the potential.
c
      IF ((KV.LT.0).AND.(EO.GT.VMAX)) THEN
        IF (LTRY.LT.1) THEN
          LTRY = 1
          KV = 999
          EO = VMAX - 1.0d-2
c
c** If unsucessful, then print out a warning and exit.
c
        ELSE
          IF (AWO.NE.0) THEN
            WRITE(6,605)
            WRITE(6,606) NF, EO, VMAX
          END IF
          AFLAG = -1
          GOTO 200
        END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If trying to find the first vibrational level (v=0), then double the
c   zero point energy estimation (AO).
c
c                         E(v=0) = VMIN + IQT*AO
c
      ELSEIF ((KV.LT.0).AND.(NF.EQ.0)) THEN
        IF (IQT.GT.1) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,611)
            WRITE(6,620) IQT, EO
          END IF
c
c** If this fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -2
            GOTO 200
          END IF
        END IF
        IQT = IQT + 1
        EO = VMIN + DBLE(IQT)*AO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If trying to find other vibrational levels (v > 0) then switch to
c   use of differences for estimating spacing.
c
      ELSEIF ((KV.LT.0).AND.(NF.GT.0)) THEN
        IF (IVDIF.GT.0) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,612)
            WRITE(6,621) NF,IVDIF
          END IF
c
c** If differences fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -3
            GOTO 200
          END IF
        END IF
        IVDIF = 1
        IF (INNER.EQ.0) THEN
          CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
        ELSE
          CALL DTENG(IEG,NF,NIWL,IWL,NVIBMX,VMIN,GV,EO)
        END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If first level found isn't v=0, try up to 3 times to 'harmonically'
c   estimate improved trial ground state energy.
c
c           E(v=0) = E(v=KV) - (E(v=KV) - VMIN)/(1 + KV/2)
c
      ELSEIF ((KV.GT.0).AND.(NF.EQ.0)) THEN
        IF (IZPE.GT.3) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,611)
            WRITE(6,622) IZPE,GV(0),KV,EO
          END IF
c
c** If differences fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -4
            GOTO 200
          END IF
        END IF
        IZPE = IZPE + 1
        EO = EO - (EO-VMIN)/(Z1+ZH/KV)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the next three cases, KV >= 0 and NF > 0.
c** If the calculated vibrational level is less than the next expected
c   level, then the estimated trial energy is too low.
c** Perhaps the difference in energy between vibrational levels v and
c   v-1 is much greater than the energy between levels v-1 and v-2.
c
c                 E(v) - E(v-1) >> E(v-1) - E (v-2)
c
c   In which case (most likely a potential with a shelf), try twice to
c   estimate a higher trial energy.
c
c   E(v) = E(v-1) + (1+IEG/2) * (2*(E(v-1)-E(v-2)) - (E(v-2)-E(v-3)))
c
      ELSEIF (KV.LT.NF) THEN
        IF (IEG.GT.1) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,612)
            WRITE(6,623) NF, KV
          END IF
c
c** If this fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -5
            GOTO 200
          END IF
        END IF
        IEG = IEG + 1
c
c** If a second minimum is present, then the next vibrational level may
c   be in the inner well. If so, use the inner well vibrational levels
c   to estimate the next trial energy.
c
        IF (INNER.EQ.0) THEN
          CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
        ELSE
          CALL DTENG(IEG,NF,NIWL,IWL,NVIBMX,VMIN,GV,EO)
        END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated vibrational level is the next expected level, then
c   continue.
c
      ELSEIF (KV.EQ.NF) THEN
        NF = NF + 1
        GV(KV) = EO
        INNR(KV) = INNER
        LTRY = 0
        CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c-----------------------------------------------------------------------
c** To ease confusion when using a potential with a second minimum, keep
c   track of levels that are in the outer well seperate from the levels
c   in the inner well.
c** First, calculate the rotational constant for this level (Bv).
c
        BV(KV) = ZH*((SWF(NBEG)/(RMIN+DBLE(NBEG-1)*RH))**2
     1           + (SWF(NEND)/(RMIN+DBLE(NEND-1)*RH))**2)
        DO I = NBEG+1,NEND-1
          BV(KV) = BV(KV) + (SWF(I)/(RMIN+DBLE(I-1)*RH))**2
        END DO
        BV(KV) = BV(KV)*RH/BZ
c
c** Double check that the calculated level is in fact located in the
c   correct well. This can be done (for v <> 0) by comparing the Bv
c   value of the new level and with the Bv values in each well. If the
c   difference is greater than 1.5 times the difference in the other
c   well, then the calculated level is probably in the wrong well.
c
        IF (NOWL.GT.0) THEN
          BVDOUT = DABS(BV(KV) - BV(OWL(NOWL-1)))
        ELSE
          BVDOUT = 9999.9d0
        END IF
        IF (NIWL.GT.0) THEN
          BVDIN = DABS(BV(KV) - BV(IWL(NIWL-1)))
        ELSE
          BVDIN = 9999.9d0
        END IF
        IF (INNER.EQ.0) THEN
          IF ((NOWL.GT.0).AND.(BVDOUT.GT.(1.5d0*BVDIN))) THEN
            IF (MCO.GE.1) THEN
              WRITE(21,2113) KV,'Inner',NIWL,GV(KV)-VMIN,BV(KV)
            END IF
            IWL(NIWL) = KV
            NIWL = NIWL + 1
          ELSE
            IF (MCO.GE.1) THEN
              WRITE(21,2113) KV,'Outer',NOWL,GV(KV)-VMIN,BV(KV)
            END IF
            OWL(NOWL) = KV
            NOWL = NOWL + 1
          END IF
        ELSE
          IF ((NIWL.GT.0).AND.(BVDIN.GT.(1.5d0*BVDOUT))) THEN
            IF (MCO.GE.1) THEN
              WRITE(21,2113) KV,'Outer',NOWL,GV(KV)-VMIN,BV(KV)
            END IF
            OWL(NOWL) = KV
            NOWL = NOWL + 1
          ELSE
            IF (MCO.GE.1) THEN
              WRITE(21,2113) KV,'Inner',NIWL,GV(KV)-VMIN,BV(KV)
            END IF
            IWL(NIWL) = KV
            NIWL = NIWL + 1
          END IF
          INNER = 0
        END IF
c-----------------------------------------------------------------------
c** Now estimate trial energy for next higher vibrational energy level
c   by using the Near-Dissociation Theory result that:
c
c                  (binding energy)**((NCN-2)/(2*NCN))
c
c   is (at least locally) linear in vibrational quantum number.
c
        IF (NF.EQ.1) THEN
          VDMV = ZH/(((VMAX-VMIN)/(VMAX-GV(0)))**PWI - Z1)
        ELSE
          VDMV = Z1/(((VMAX-GV(NF-2))/(VMAX-GV(NF-1)))**PWI - Z1)
        END IF
c
c** If unable to calculate the next trial energy, see if all of the
c   desired levels have been calculated. If not then turn on the warning
c   flag and quit, otherwise print out success message and quit.
c
        IF ((VDMV.LT.Z1).AND.(NCN.GT.2)) THEN
          IF (NF.LE.KVMAX) THEN
            AFLAG = -1
            WRITE(6,640) JROT, KV + VDMV
          ELSEIF (AWO.GT.0) THEN
            WRITE(6,630) KVMAX
          END IF
          GOTO 200
        END IF
c
c** Now calculate the next trial energy.
c
        EO = VMAX - (VMAX-GV(NF-1))*(Z1-Z1/VDMV)**PW
c
c** However, if the level is above the dissociation limit (for
c   potentials with barriers) then use differences to calculate the
c   next trial energy.
c
        IF (EO.GT.VMAX) CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated vibrational level is higher then the next expected
c   level, then try thrice to interpolate harmonically for the missed
c   level.
c
c                   E(v) = E(v-1) + (E(KV) - E(v-1)) / 2
c
      ELSEIF (KV.GT.NF) THEN
        IF (IVCOR.GT.2) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,612)
            WRITE(6,624) IVCOR,KV,EO,(NF-1),GV(NF-1)
          END IF
c
c** If interpolation fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   missing level.
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -6
            GOTO 200
          END IF
        END IF
        IVCOR = IVCOR + 1
c
c** Use NDE theory to determine the missing level.
c
        IF (NPMIN.EQ.1) THEN
          VDU = (VPMAX(1)-EO)**PWI
          VDL = (VPMAX(1)-GV(OWL(NOWL-1)))**PWI
          EO = VPMAX(1) - (VDL + (VDU - VDL) / DBLE(KV - NF + 1))**PW
        ELSEIF (((INNER.EQ.1).AND.(NIWL.GT.0)).OR.
     1          ((INNER.EQ.0).AND.(NOWL.EQ.0))) THEN
          VDU = (VPMAX(1)-EO)**PWI
          VDL = (VPMAX(1)-GV(IWL(NIWL-1)))**PWI
          EO = VPMAX(1) - (VDL + (VDU - VDL) / DBLE(KV - NF + 1))**PW
        ELSEIF (((INNER.EQ.0).AND.(NOWL.GT.0)).OR.
     1          ((INNER.EQ.1).AND.(NIWL.EQ.0))) THEN
          VDU = (VPMAX(2)-EO)**PWI
          VDL = (VPMAX(2)-GV(OWL(NOWL-1)))**PWI
          EO = VPMAX(2) - (VDL + (VDU - VDL) / DBLE(KV - NF + 1))**PW
        END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If an unknown case occurs (quite impossible but don't quote me on
c   it) then write out an error message and exit.
c
      ELSE
        IF (AWO.NE.0) THEN
          WRITE(6,610)
          WRITE(6,666) KV,NF
        END IF
        AFLAG = -7
        GOTO 200
      END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set KV to the next vibrational level to be found unless looking for
c   the highest vibrational level.
c
      IF (KV.NE.999) KV = NF
c
c** If still haven't found all of the vibrational levels then
c   look for the next vibrational level.
c
      IF ((KV.LE.KVMAX).OR.(KV.EQ.999)) GOTO 10
c
c** Otherwise, print out a message saying that all is well.
c
      IF ((KV.GT.KVMAX).AND.(AWO.GT.0)) WRITE(6,630) KVMAX
c
c** If the potential has levels in a second minimum, then print out a
c   list of those levels to channel-21 if desired.
c
      IF ((NIWL.GT.0).AND.(NOWL.GT.0)) THEN
        IF (MCO.GE.1) WRITE(21,2114) NIWL, NOWL
        IF (AWO.NE.0) THEN
          WRITE(6,605)
          WRITE(6,607)
        END IF
        AFLAG = -8
      END IF
c
c** If an error has occured, then set KVMAX to the quantum number of the
c   highest vibrational level found and print out that quantum number
c   and the energy of that level.
c
  200 IF (AFLAG.LT.0) THEN
        KVMAX = NF - 1
        IF (AWO.NE.0) WRITE(6,626) KVMAX, GV(KVMAX)
      END IF
      IF (MCO.GE.1) WRITE(21,2100)
      RETURN
c-----------------------------------------------------------------------
  601 FORMAT(/' Solve by matching inward and outward solutions at')
  602 FORMAT(' the outermost wave function maximum, S(max), where  R = R
     1R(M)')
  603 FORMAT(' the innermost turning point  R1 = R(M)')
  605 FORMAT(/'  *** ALF WARNING ***')
  606 FORMAT(4X,'Next estimated trial energy  E(v=',I3,') =',G15.8/4X,
     1'is greater than the potential maximum  VMAX =',G15.8)
  607 FORMAT(4X,'Potential found to have a second minimum.')
  610 FORMAT(/'  *** ALF ERROR ***')
  611 FORMAT(4X,'Attempt to find zero point level fails!')
  612 FORMAT(4X,'Attempt to find next higher vibrational level fails!')
  613 FORMAT(4X,'Number of vib levels requested=',i4,' exceeds internal
     1ALF array dimension  NVIBMX=',i4)
  614 FORMAT(4X,'Unable to find a potential minimum.')
  615 FORMAT(4X,'There are',I3,'  potential ',A6,' in this potential.')
  616 FORMAT(4X,'The potential turns over in the short range region at R
     1 = ',G15.8)
  617 FORMAT(4X,'VMAX =',G15.8,' found to be less than VMIN =',G15.8)
  620 FORMAT(4X,'Use of energy ',I1,'0% up the potential well (E =',
     1G15.8,')'/4X,' fails to produce a viable vibrational eigenstate.')
  621 FORMAT(4X,'Use of differences to estimate the energy for the next'
     1/4X,' vibrational level (v=',I3,') failed after',I3,'  attempt.')
  622 FORMAT(4X,'After',I3,' tries to harmonically estimate the zero-poi
     1nt energy,'/4X,' initial trial energy',G15.8,'   had yielded   E(v
     2=',I3,') =',G15.8)
  623 FORMAT(4X,'Expecting to find level (v=',I3,') but found level (v='
     1,I3,')')
  624 FORMAT(4X,'After',I3,' tries, failed to interpolate trial energy b
     1etween'/4X,'E(v=',I3,') =',G15.8,'   and   E(v=',I3,') =',G15.8)
  626 FORMAT(4X,'The highest calculated level is  E(v=',I3,') =',G15.8)
  630 FORMAT(/' ALF successfully finds all vibrational levels up to KVMA
     1X=',I3)
  640 FORMAT(/' ALF finds all  J=',i3,'  vib. levels below  vD=',F7.3,
     1  '  estimated by N-D theory')
  650 FORMAT(/' There were',I3,'  potential minima found with the absolu
     1te minimum'/4X,'VMIN =',G15.8,'  cm-1.')
  651 FORMAT(/' There were',I3,'  potential maxima found with the absolu
     1te maximum'/4X,'VMAX =',G15.8,'  cm-1.')
  666 FORMAT(4X,'Undefined case for automatic search.'/,4X,'Values of KV
     1 =',I3,'  and NF =',I3)
 2100 FORMAT(/1X,39('=='))
 2101 FORMAT(/1X,39('--'))
 2110 FORMAT(/' Limits and increment of integration (in Angstroms):'
     1 /'    RMIN =',F6.3,'    RMAX =',F7.3,'    RH =',F9.6,
     2 //' Generate    BZ =',G19.12,' ((1/cm-1)(1/Angstroms**2))'
     3 /' from ZMU:',F15.11,' (amu)')
 2111 FORMAT(/' Eigenvalue convergence criterion is   EPS =',G11.4,'(cm-
     11)')
 2112 FORMAT(/' Calculating properties of the potential described above.
     1 '/' Use Airy function at 3-rd turning point as outer boundary'
     2 /' condition for quasibound levels.')
 2113 FORMAT(' v=',I3,4X,'v(',A5,')=',I3,4X,'Gv=',F16.9,4X,'Bv=',F16.12)
 2114 FORMAT(/' Found',I4,' level(s) in the inner well and',I4,' level(s
     1) in the outer well.')
      END
c***********************************************************************
      SUBROUTINE INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c***********************************************************************
c** This subroutine reinitializes the condition flags when considering a
c   new case (found next vibrational level or finding level in inner
c   well - INNER = 1).
c
c** On entry and exit:
c    IQT     Case when KV < 0 and NF = 0
c            determines the value used for the initial trial energy.
c    IVDIF   Case when KV < 0 and NF > 0
c            is the flag denoting the use of differences to calculate
c            trial energies.
c    IZPE    Case when KV > 0 and NF = 0
c            is the number of times the zero point energy (v = 0) has
c            been estimated harmonically.
c    IEG     Case when KV < NF and NF > 0
c            are the number of times that a larger trial energy is used
c            to find the next level.
c    IVCOR   Case when KV > NF and NF > 0
c            are the number of times that a smaller trial energy is used
c            to find the next level.
c
      INTEGER IZPE,IVDIF,IVCOR,IQT,IEG
c
      IQT = 1
      IVDIF = 0
      IZPE = 0
      IEG = 0
      IVCOR = 0
      RETURN
      END
c***********************************************************************
      SUBROUTINE DTENG(IEG,NF,NVEL,VEL,NVIBMX,VMIN,GV,EO)
c***********************************************************************
c** This subroutine calculates the next trial energy using differences.
c
c** On entry:
c    IEG     factor by which a larger trial energy should be calculated:
c              NVEL = 2 : Increase correction by increments of 25%
c              NVEL > 2 : Increase correction by increments of 50%
c    NF      is the highest calculated vibrational level.
c    NVEL    is the number of levels found in the potential well.
c    VEL(v)  keeps track of all levels in the potential well.
c    NVIBMX  is the maximum number of vibrational levels (dimension).
c    VMIN    is the absolute value of the potential minimum (cm-1).
c    GV(v)   contains the vibrational energy level spacings
c            and rotational constants for each level (cm-1).
c
c** On exit:
c    EO      is the calculated trial energy.
c
      INTEGER IEG,NF,NVEL,NVIBMX,VEL(0:NVIBMX)
c
      REAL*8 VMIN,GV(0:NVIBMX),EO,ZQ,ZH,Z1,Z2
c
      DATA ZQ/0.25D0/,ZH/0.5D0/,Z1/1.D0/,Z2/2.D0/
c
c** If determining the first (non-zero point energy) level in the well,
c   then use the last determined level in the other well plus a larger
c   than harmonic correction that becomes smaller with each new
c   itteration.
c
c           E(v=0) = E(v=NF-1) + (E(v=NF-1)-VMIN)/(NF-1+IEG/4)
c
      IF (NVEL.EQ.0) THEN
        EO = GV(NF-1) + (GV(NF-1) - VMIN)/(NF - 1 + ZQ*DBLE(IEG))
c
c** Try to get v = 1 using smaller-than-harmonic spacing.
c
c                 E(v=1) = E(v=0) + 1.3*(E(v=0)-VMIN)
c
      ELSEIF (NVEL.EQ.1) THEN
        EO = GV(VEL(0)) + 1.3d0*(GV(VEL(0))-VMIN)
c
c** Try to get v = 2 using a sequentially increasing correction.
c
c             E(v=2) = E(v=1) + (0.8+IEG/4)*(E(v=1)-E(v=0))
c
      ELSEIF (NVEL.EQ.2) THEN
        EO = GV(VEL(1)) + (0.8d0+DBLE(IEG)*ZQ)*(GV(VEL(1))-GV(VEL(0)))
c
c** Try to get v > 2 using a sequentially increasing correction.
c
c        E(v) = E(v-1) + (1.0+IEG/2)*(2.0*E(v-1)-3.0*E(v-2)+E(v-3))
c
      ELSE
        EO = GV(VEL(NVEL-1)) + (Z1+DBLE(IEG)*ZH)
     1       *(Z2*GV(VEL(NVEL-1))-3.0d0*GV(VEL(NVEL-2))+GV(VEL(NVEL-3)))
      END IF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c****** R.J. Le Roy  subroutine SCHRQ, last updated  16 May 2000 *******
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2000  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** SCHRQ solves radial Schrodinger equation in dimensionless form
c  d2WF/dR2 = - (E-V(R))*WF(R) ,  where WF(I) is the wave function.
c** Integrate by Numerov method over N mesh points with increment
c  H=RH across range beginning at RMIN .
c** Input trial energy EO, eigenvalue convergence criterion EEPS
c  potential asymptote VLIM, and all returned energies (EO, GAMA & VMAX)
c  have units (cm-1).
c** On entry, the input potential V(I) must include the centrifugal
c  term and the factor:  'BFCT'=2*mu*(2*pi*RH/hPLANCK)**2  (1/cm-1) ,
c  which is also internally incorporated into EO, VLIM & EEPS.
c* Note that these reduced quantities (& the internal eigenvalue E)
c  contain a factor of the squared integration increment  RH**2 .
c  This saves arithmetic work in the innermost loop of the algorithm.
c** For energy in (cm-1), BFCT=ZMU(u)*H(Angst)**2/16.85762908 (1/cm-1)
c** INNER specifies wave function matching (& initiation) condition *
c*  For INNER = 0 , match inward & outward solutions at outermost
c  wave function maximum;  otherwise match at inner edge of classically
c  allowed region. ** INNER<0  uses zero slope inner boundary condition.
c** For most normal cases set INNER=0 ,  but ......
c* to find "inner-well-dominated" solutions of an asymmetric double
c    minimum potential, set  INNER > 0 .
c* To find symmetric eigenfunctions of a symmetric potential, set
c    INNER < 0  & start integration (set RMIN) at potential mid point.
c----------------------------------------------------------------------
      SUBROUTINE SCHRQ(KV,JROT,EO,GAMA,VMAX,VLIM,V,WF,BFCT,EEPS,RMIN,
     1                                 RH,N,NBEG,NEND,INNER,IWR,LPRWF)
c----------------------------------------------------------------------
c** Output vibrational quantum number KV, eigenvalue EO, normalized
c  wave function WF(I), and range, NBEG .le. I .le. NEND  over
c  which WF(I) is defined. *** Have set  WF(I)=0  outside this range.
c* (NBEG,NEND), defined by requiring  abs(WF(I)) < RATST=1.D-9  outside.
c** If(LPRWF.gt.0) print wavefunction WF(I) every LPRWF-th point.
c* If(LPRWF.lt.0) "punch" (i.e., WRITE(10,XXX)) every |LPRWF|-th point
c  of the wave function on disk starting at R(NBEG) with step size
c  of  IPSIQ=|LPRWF|*RH.
c** For energies above the potential asymptote VLIM, locate quasibound
c  levels using Airy function boundary condition and return the level
c  width GAMA and barrier height VMAX, as well as EO.
c** ERROR condition on return is  KV < 0 ; usually KV=-1, but return
c  KV=-2 if error appears to arise from too low trial energy.
c** If(IWR.ne.0) print error & warning descriptions
c  If (IWR.gt.0) also print final eigenvalues & node count.
c  If (IWR.ge.2) also show end-of-range wave function amplitudes
c  If (IWR.ge.3) print also intermediate trial eigenvalues, etc.
c** If input KV.ge.998 , tries to find highest bound level, and
c  trial energy should be only slightly less than VLIM.
c** If input KV < -10 , use log-derivative outer boundary condition at
c  mesh point |KV| , based on incoming value of wave function WF(|KV|)
c  and of the wavefunction derivative at that point, SPNEND, which is
c  brought in as WF(|KV|-1).  For a hard wall condition at mesh point
c  |KV|, set WF(|KV|)=0 and WF(|KV|-1)= -1 before entry.
c----------------------------------------------------------------------
c++ "SCHRQ" calls subroutineas "QBOUND" and "WIDTH", and the latter
c++ calls "LEVQAD" .
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,IBEGIN,ICOR,IJ,IJK,INNER,IPSID,IQTST,IT,ITER,ITP1,
     1         ITP1P,ITP3,IWR,J,JJ,J1,J2,JPSIQ,JQTST,JROT,
     2         KKV,KV,KVIN,LPRWF,M,MS,MSAVE,
     3         N,NBEG,NBEGB,NBEG2,NDN,NEND,NENDCH,NLINES,NPR
      REAL*8  BFCT,DE,DEP,DEPRN,DF,DOLD,DSOC,
     2        E,EEPS,EO,EPS,F,FX,GAMA,GI,GN,H,H2,HT,PROD,PPROD,
     3        RATIN,RATOUT,RATST,RH,RINC,RMIN,RMINN,RR,RSTT,RWR(20),
     4        WF(N),SB,SI,SM,SN,SNEND,SPNEND,SRTGI,SRTGN,SWR(20),
     5        V(N),VLIM,VMAX,VMX,VPR,
     6        WKBTST,XEND,XPR,XPW,DXPW,Y1,Y2,Y3,YIN,YM,YOUT,
     7        Z0,Z1,Z2,ZH
      DATA Z0/0.D0/,ZH/0.5D0/,Z1/1.D0/,Z2/2.D0/
      DATA RATST/1.D-9/,XPW/20.72d0/
      DATA NDN/15/
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DXPW= (XPW+ 2.30d0)/NDN
      KVIN= KV
      KV= -1
      RMINN= RMIN-RH
      GAMA= Z0
      VMAX= VLIM
      VMX= VMAX*BFCT
      H= RH
      H2= H*H
      HT= Z1/12.D+0
      E= EO*BFCT
      EPS= EEPS*BFCT
      DSOC= VLIM*BFCT
      DE= Z0
      RATIN= Z0
      RATOUT= Z0
      IF(IWR.GT.2) THEN
          IF(KVIN.GE.998) then
              WRITE(6,610) EO
            ELSE
              WRITE(6,601) KVIN,JROT,EO,INNER
            ENDIF
          WRITE(6,602)
        ENDIF
      NEND= N
      IF(KVIN.LT.-10) THEN
          NEND= -KVIN
          SNEND= WF(NEND)
          SPNEND= WF(NEND-1)
          ENDIF
      JQTST = 0
c** Start iterative loop; try to converge for up to 15 iterations.
      DO 90 IT= 1,15
          ITER= IT
          IF(INNER.NE.0) GO TO 38
   10     IF(KVIN.LT.-10) THEN
c** If desired, (KVIN < -10) outer boundary set at NEND=|KVIN| and
c  initialize wavefunction with log-derivative condition based on value
c  WF(NEND) & derivative SPNEND at that mesh point (brought in in CALL)
              GN= V(NEND)-E
              GI= V(NEND-1)-E
              SB= SNEND
              SI= SB*(Z1+ ZH*GN)- RH*SPNEND
              GO TO 24
              END IF
          IF(E.GE.DSOC) THEN
c** For quasibound levels, initialize wave function in "QBOUND"
              CALL QBOUND(KVIN,JROT,E,EO,VMX,DSOC,V,RMIN,H,GN,GI,
     1                                 SB,SI,N,ITP3,IWR,IQTST,BFCT,IT)
              NEND= ITP3
              VMAX= VMX/BFCT
              IF(IQTST.GT.0) GO TO 24
              IF(IQTST.LT.0) THEN
                  JQTST = JQTST+IQTST
                  IF((JQTST.LE.-2).OR.(VMAX.LT.VLIM)) GO TO 999
c** Try up to once to find level using trial value just below maximum
                  EO = VMAX-0.1D0
                  E = EO*BFCT
                  GO TO 90
                  ENDIF
              GO TO 20
              ENDIF
c** For  E < DSOC  begin inward integration by using JWKB to estimate
c  optimum (minimum) inward starting point which will still give
c  RATOUT < RATST = exp(-XPW) (ca. 1.d-9) [not needed after 1'st 2 ITER]
          IF(ITER.LE.2) THEN
              NEND= N
c ... first do rough inward search for outermost turning point
              DO  M= N,1,-NDN
                  MS= M
                  GI= V(M)- E
                  IF(GI.LE.0.D0) GO TO 12
                  GN= GI
                  ENDDO
              IF(IWR.NE.0) WRITE(6,611)
              GO TO 999
   12         IF(MS.GE.N) GO TO 998
              FX= GN/(GI-GN)
              SM= ZH*(Z1+ FX)*DSQRT(GN)
              MS= MS+ 2*NDN
              IF(MS.GE.N) GO TO 20
c ... now integrate exponent till JWKB wave fx. would be negligible
              DO  M= MS,N,NDN
                  NEND= M
                  SM= SM+ DSQRT(V(M)- E)
                  IF(SM.GT.DXPW) GO TO 18
                  ENDDO
   18         IF(NEND.LT.N) NEND= NEND+ NDN
              ENDIF
c** For truly bound state initialize wave function as 1-st order WKB
c   solution increasing inward
   20     GN= V(NEND)- E
          GI= V(NEND-1)- E
          MS= NEND-1
          IF(GI.LT.0.d0) GO TO 998
          SRTGN= DSQRT(GN)
          SRTGI= DSQRT(GI)
          SB= Z1
          SI= SB*DSQRT(SRTGN/SRTGI)*DEXP((SRTGN+SRTGI)/Z2)
          IF(SB.GT.SI) THEN
c WOOPS - JWKB gives inward DEcreasing solution, so initialize with node
              IF(IWR.NE.0) WRITE(6,618) JROT,EO,SB/SI
              SI= Z1
              SB= Z0
              ENDIF
   24     M= NEND-1
          Y1= (Z1-HT*GN)*SB
          Y2= (Z1-HT*GI)*SI
          WF(NEND)= SB
          WF(NEND-1)= SI
          MS= NEND
          NENDCH= NEND
          IBEGIN= 3
          IF(INNER.NE.0) IBEGIN= ITP1+2
c** Actual inward integration loop starts here
          DO  I= IBEGIN,NEND
              M= M-1
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(M)-E
              SB= SI
              SI= Y3/(Z1-HT*GI)
              WF(M)= SI
              IF(DABS(SI).GE.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I)  in classically
c  forbidden region where  (V(I) .gt. E)
                  SI= Z1/SI
                  DO  J= M,MS
                      WF(J)= WF(J)*SI
                      ENDDO
                  NENDCH= MS
                  MS= M
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SB= SB*SI
                  SI= Z1
                  ENDIF
              Y1= Y2
              Y2= Y3
c** Test for outermost maximum of wave function.
cc            IF((INNER.EQ.0).AND.(SI.LE.SB)) GO TO 32
c** Test for outer well turning point
              IF((INNER.EQ.0).AND.(GI.lt.0.d0)) GO TO 32
              ENDDO
          IF(INNER.EQ.0) THEN
c** Error mode ... find no wave function maximum.
              KV= -2
              IF(IWR.NE.0) WRITE(6,616) KV,JROT,EO
              GO TO 999
              ENDIF
c** Scale outer part of wave function before proceding
   32     SI= Z1/SI
          MSAVE= M
          RR= RMINN+MSAVE*H
          YIN= Y1*SI
          RATOUT= WF(NEND)*SI
          NEND= NENDCH
          DO  J= MSAVE,NEND
              WF(J)= WF(J)*SI
              ENDDO
          IF(INNER.NE.0) GO TO 70
c-------------------------------------------------------------------
c** Set up to prepare for outward integration **********************
   38     NBEG= 1
          IF(INNER.LT.0) THEN
c** Option to initialize with zero slope at beginning of the range
              SB= Z1
              GN= V(1)-E
              Y1= SB*(Z1-HT*GN)
              Y2= Y1+GN*SB/Z2
              GI= V(2)-E
              SI= Y2/(Z1-HT*GI)
            ELSE
c** Initialize outward integration with a node at beginning of range
   40         GN= V(NBEG)-E
              IF(GN.GT.10.D0) THEN
c** If potential has [V(1)-E] so high that H is (locally) much too
c  large, then shift inner starting point outward.
                  NBEG= NBEG+1
                  IF(NBEG.LT.N) GO TO 40
                  IF(IWR.NE.0) WRITE(6,613)
                  GO TO 999
                  ENDIF
              IF((ITER.LE.1).AND.(IWR.NE.0)) THEN
                  IF(NBEG.GT.1) WRITE(6,609) JROT,EO,NBEG
                  IF(GN.LE.Z0) WRITE(6,604) JROT,EO,E,V(NBEG),NBEG
                  ENDIF
c** Initialize outward wave function with a node:  WF(NBEG) = 0.
              SB= Z0
              SI= Z1
              GI= V(NBEG+1)-E
              Y1= SB*(Z1-HT*GN)
              Y2= SI*(Z1-HT*GI)
            ENDIF
c
          WF(NBEG)= SB
          WF(NBEG+1)= SI
          NBEGB= NBEG
          NBEG2= NBEG+2
          IF(INNER.NE.0) MSAVE= N
c** Actual outward integration loops start here
          DO  I= NBEG2,MSAVE
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(I)-E
              SI= Y3/(Z1-HT*GI)
              WF(I)= SI
              IF(DABS(SI).GE.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I)  in classically forbidden
c  region where  V(I) .gt. E
                  SI= Z1/SI
                  NBEG= NBEGB
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
                  NBEGB= I
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SI= Z1
                  ENDIF
              Y1= Y2
              Y2= Y3
              ITP1= I
c** Exit from this loop at onset of classically allowed region
              IF(GI.LE.Z0) GO TO 52
              ENDDO
          MS= MSAVE
          IF((INNER.EQ.0).AND.(GN.LE.Z0)) GO TO 60
          IF(IWR.NE.0) WRITE(6,612) KVIN,JROT,EO,MSAVE
          GO TO 999
   52     ITP1P= ITP1+1
          MS= ITP1
          IF(INNER.NE.0) GO TO 60
          DO  I= ITP1P,MSAVE
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(I)-E
              SI= Y3/(Z1-HT*GI)
              WF(I)= SI
              IF(DABS(SI).GT.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I) , as needed.
                  SI= Z1/SI
                  NBEG= NBEGB
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
                  NBEGB= I
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SI= Z1
                  ENDIF
              Y1= Y2
              Y2= Y3
              ENDDO
          MS= MSAVE
c** Finished outward integration.  Normalize w.r.t. WF(MSAVE)
   60     SI= Z1/SI
          YOUT= Y1*SI
          YM= Y2*SI
          RATIN= WF(NBEG+1)*SI
          DO  I= NBEG,MS
              WF(I)= WF(I)*SI
              ENDDO
          IF(INNER.NE.0) GO TO 10
c----- Finished numerical integration ... now correct trial energy
c** DF*H  is the integral of  (WF(I))**2 dR
   70     DF= Z0
          DO  J= NBEG,NEND
              DF= DF+WF(J)**2
              ENDDO
c** Add edge correction to DF assuming wave function dies off as simple
c  exponential past R(NEND);  matters only if WF(NEND) unusually large.
          IF((E.LE.DSOC).AND.(WF(NEND).NE.0)) THEN
              IF((KVIN.GE.-10).AND.(WF(NEND-1)/WF(NEND).GT.Z1))
     1              DF= DF+ WF(NEND)**2/(Z2*DLOG(WF(NEND-1)/WF(NEND)))
              ENDIF
          F= (-YOUT-YIN+Z2*YM+GI)
          DOLD= DE
          IF(DABS(F).LE.1.D+30) THEN
              DE= F/DF
            ELSE
              F= 9.9D+30
              DF= F
              DE= DABS(0.01D+0 *(DSOC-E))
            ENDIF
          IF(IWR.GT.2) THEN
              DEPRN = DE/BFCT
              XEND= RMINN+NEND*H
c** RATIN & RATOUT  are wave fx. amplitude at inner/outer ends of range
c  relative to its value at outermost extremum.
              WRITE(6,603) IT,EO,F,DF,DEPRN,MSAVE,RR,RATIN,RATOUT,
     1                                                  XEND,NBEG,ITP1
              ENDIF
c** Test trial eigenvalue for convergence
          IF(DABS(DE).LE.DABS(EPS)) GO TO 100
          E= E+DE
c** KV.ge.998  Option ... Search for highest bound level.  Adjust new
c  trial energy downward if it would have been above dissociation.
          IF((KVIN.GE.998).AND.(E.GT.VMX)) E= VMX- 2.d0*(VMX-E+DE)
          EO= E/BFCT
          IF((IT.GT.4).AND.(DABS(DE).GE.DABS(DOLD)).AND.
     1                                         ((DOLD*DE).LE.Z0)) THEN
c** Adjust energy increment if having convergence difficulties.  Not
c  usually needed except for some quasibounds extremely near  VMAX .
              ICOR= ICOR+1
              DEP= DE/BFCT
              IF(IWR.NE.0) WRITE(6,617) IT,DEP
              DE= ZH*DE
              E= E-DE
              EO= E/BFCT
              ENDIF
   90     CONTINUE
c** End of iterative loop which searches for eigenvalue ************
c-------------------------------------------------------------------*
c** Convergence fails, so return in error condition
      E= E-DE
      EO= E/BFCT
      DEPRN= DE/BFCT
      IF(IWR.NE.0) WRITE(6,620) KVIN,JROT,ITER,DEPRN
      GO TO 999
  100 IF(IWR.NE.0) THEN
          IF(IWR.GE.3) WRITE(6,619)
          IF((DABS(RATIN).GT.RATST).AND.(INNER.GE.0))
     1                                      WRITE(6,614) JROT,EO,RATIN
          IF((E.LT.DSOC).AND.(DABS(RATOUT).GT.RATST)) THEN
              WKBTST= ZH*DABS(V(NEND)-V(NEND-1))/DSQRT((V(NEND)-E)**3)
              IF(WKBTST.GT.1.d-3)WRITE(6,615)JROT,EO,RATOUT,RATST,WKBTST
              ENDIF
          ENDIF
      KKV = 0
c** Perform node count on converged solution
      PROD= WF(ITP1)*WF(ITP1-1)
      J1= ITP1+1
      J2= NEND-1
      DO  J= J1, J2
          PPROD= PROD
          PROD= WF(J)*WF(J-1)
          IF((PPROD.LE.Z0).AND.(PROD.GT.Z0)) KKV= KKV+1
          ENDDO
      KV = KKV
c** Normalize & find interval (NBEG,NEND) where WF(I) is non-negligible
      SN= Z1/DSQRT(H*DF)
      DO  I= NBEG,NEND
          WF(I)= WF(I)*SN
          ENDDO
      IF(ITP1.LE.1) GO TO 122
      J= ITP1P
      DO  I= 1,ITP1
          J= J-1
          IF(DABS(WF(J)).LT.RATST) GO TO 119
          ENDDO
  119 NBEG= J
      IF(NBEG.LE.1) GO TO 122
      J= J-1
      DO  I= 1,J
          WF(I)= Z0
          ENDDO
  122 IF(KVIN.GE.-10) THEN
c** For "non-wall" cases, move NEND inward to where wavefunction
c  "non-negligible"
          J= NEND-1
          DO  I= NBEG,NEND
              IF(DABS(WF(J)).GT.RATST) GO TO 126
              J= J-1
              ENDDO
  126     NEND= J+1
          END IF
      IF(NEND.LT.N) THEN
c** Zero out wavefunction array at distances past NEND
          DO  I= NEND+1,N
              WF(I)= Z0
              ENDDO
          ENDIF
      IF(LPRWF.LT.0) THEN
c** If desired, write every |LPRWF|-th point of the wave function
c  to a file on channel-10, starting at the NBEG-th mesh point.
          JPSIQ= -LPRWF
          NPR= 1+(NEND-NBEG)/JPSIQ
          RINC= RH*JPSIQ
          RSTT= RMINN+NBEG*RH
c** Write every JPSIQ-th point of the wave function for level  v=KV
c  J=JROT , beginning at mesh point NBEG & distance RSTT where
c  the NPR values written separated by mesh step RINC=JPSIQ*RH
          WRITE(10,701) KV,JROT,EO,NPR,RSTT,RINC,NBEG,JPSIQ
          WRITE(10,702) (RMINN+I*RH,WF(I),I=NBEG,NEND,JPSIQ)
          GO TO 140
          ENDIF
c** Print solutions every  LPRWF-th  point, 6 to a line, in columns.
      IF(LPRWF.GT.0) THEN
          NLINES= ((1+(NEND-NBEG)/LPRWF)+3)/4
          IPSID= LPRWF*NLINES
          WRITE(6,605) KV,JROT,EO
          DO  J= 1,NLINES
              JJ= NBEG+(J-1)*LPRWF
              IJK= 0
              DO  IJ= JJ,NEND,IPSID
                  IJK= IJK+1
                  RWR(IJK)= RMINN+IJ*H
                  SWR(IJK)= WF(IJ)
                  ENDDO
              WRITE(6,606) (RWR(I),SWR(I),I= 1,IJK)
              ENDDO
          ENDIF
  140 IF(IWR.EQ.1) WRITE(6,607) KV,JROT,EO
cc
cc    IF(IWR.NE.0) WRITE(6,699)  rminn+itp1*rh,eO,rminn+msave*rh,eo
cc699 format('    & turning points:',2(f8.5,f11.4))
cc
      IF(IWR.GE.2) WRITE(6,607) KV,JROT,EO,ITER,RR,RATIN,RATOUT
c** For quasibound levels, calculate width in subroutine "WIDTH"
      IF((E.GT.DSOC).AND.(KVIN.GT.-10)) CALL WIDTH(KV,JROT,E,EO,DSOC,
     1  V,WF,VMX,RMIN,H,BFCT,IWR,ITP1,ITP3,INNER,N,GAMA)
      RETURN
c** ERROR condition if  E.gt.V(R)  at outer end of integration range.
  998 XPR= RMINN+MS*H
      VPR= V(MS)/BFCT
      IF(IWR.NE.0) WRITE(6,608) EO,MS,VPR,XPR,IT
c** Return in error mode
  999 KV= -1
      RETURN
  601 FORMAT(/' Solve for  v=',I3,'   J=',I3,'   ETRIAL=',1PD15.7,
     1   '  INNER=',i2,'   WF(1st) WF(NEND)' )
  602 FORMAT(' ITER    ETRIAL',8X,'F(E)      DF(E)     D(E)',
     1 5X,'M    R(M)  /WF(M)   /WF(M)  R(NEND) NBEG ITP1'/
     2  1X,96('-'))
  603 FORMAT(I4,1PD15.7,3D10.2,I5,0PF7.3,1P2D9.1,0PF8.2,I4,I5)
  604 FORMAT('      NOTE:  for   J =',I3,'   EO =',F12.4,'   E=',D13.6,
     1 ' .ge. V(R)=',D13.6,'   at initial mesh point',I6)
  605 FORMAT(/' Solution of radial Schr. equation for   E(v=',I3,',J=',
     1  I3,') =',F15.7/2x,4('    R(I)   WF(I)   ')/2X,38('--') )
  606 FORMAT(2X,4(F8.3,F11.7))
  607 FORMAT('E(v=',I3,',J=',I3,')=',F11.4,1x,I3,' Iterations',
     1  '   R(M)=',F6.3,'  WF(NBEG)/WF(M)=',1PD8.1/
     2  57x,'WF(NEND)/WF(M)=',D8.1)
  608 FORMAT(' *** SCHRQ Error:  E=',F9.2,' > V(',I5,')=',F9.2,
     1  '  at  Rmax=',F6.2,'  for  IT=',I2)
  609 FORMAT(' *** For  J=',I3,'   E=',1PD15.7,"  integration can't",
     1 ' start till past mesh'/37x,'point',I5,',  so RMIN smaller than n
     2eeded')
  610 FORMAT(/' Attempt to find the highest bound level starting from',
     1 '   ETRIAL =',1PD9.2)
  611 FORMAT(' *** SCHRQ Error:  inward search at   E=',f9.2,
     1  ' finds no classical region')
  612 FORMAT(/' *** ERROR *** for   v =',I3,'   J =',I3,'   E =',
     1  F12.4,'  Innermost turning point not found by   M = MSAVE =',I5)
  613 FORMAT(/' *** ERROR in potential array ... V(I) everywhere',
     1 ' too big to integrate with given  increment')
  614 FORMAT(' *** CAUTION *** For  J=',I3,'  E=',G15.8/16x,
     1 'WF(first)/WF(Max)=',D9.2,'  suggests  RMIN  may be too large')
  615 FORMAT(' ** CAUTION ** For  J=',I3,'  E=',1PD13.6,
     1 '  WF(NEND)/WF(Max)=',D8.1,' >',D8.1/4X,'& initialization ',
     2 'quality test ',1PD8.1,' > 1.D-3   so RMAX may be too small')
  616 FORMAT(' ** WARNING *** For  v=',I2,', J=',I3,' at  E=',
     1 G14.7,'  WF always has negative slope ... Energy too low or poten
     2tial too weak' )
  617 FORMAT(' *** SCHRQ has a convergence problem, so for  IT=',I2,
     1 '  cut  DE=',1PD10.2,'  in HALF' )
  618 FORMAT(' *** For  J=',I3,'  E=',F9.2,'  JWKB start gives  SB/SI=',
     1  1PD10.3,'  so use a node.')
  619 FORMAT(1X,96('-'))
  620 FORMAT(' *** CAUTION for  v=',I3,'  J=',I3,"  SCHRQ doesn't conver
     1ge by  ITER=",I2,'  DE=',1PD9.2)
  701 FORMAT(/2x,'Level  v=',I3,'   J=',I3,'   E=',F12.4,' ,  wave funct
     1ion at',I6,' points.'/7x,'R(1-st)=',F12.8,'   mesh=',F12.8,
     2  '   NBEG=',I4,'   |LPRWF|=',I3)
  702 FORMAT((1X,4(f9.4,f10.6)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c*******************************************************************
      SUBROUTINE QBOUND(KV,JROT,E,EO,VMX,DSOC,V,RMIN,H,GB,GI,SB,SI,N,
     1  ITP3,IWR,IQTST,BFCT,IT)
c*******************************************************************
c** Subroutine to initialize quasibound level wave function as Airy
c  function at third turning point (if possible). For the relevant
c  theory see: Le Roy & Bernstein, J.Chem.Phys. 54, 5114 (1971) and
c  Le Roy & Liu, J.Chem.Phys.69,3622-31 (1978).
c----------------------------------------------------------------------
c** IQTST  is error flag. *** If (IQTST.lt.0) initialization fails
c  so eigenvalue calculation aborts *** (IQTST.gt.0) for successful
c  Airy function initialization. *** (IQTST=0) if Airy function
c  initialization prevented because 3-rd turning point beyond
c  range, so that WKB initialization is used.
c----------------------------------------------------------------------
      INTEGER I,II,IQTST,IT,ITP3,IWR,J,JROT,K,KV,N
      REAL*8  A1,A2,A13,A23,BFCT,
     1        C1A,C2A,DF,DSOC,E,EO,FBA,FIA,FJ,GB,GBA,GI,GIA,H,
     2        RMIN,RMINN,SB,SI,SL,V(N),VMX,VMXPR,XJ1, Z1,Z3,Z6
      DATA Z1/1.D0/,Z3/3.D0/
     1  Z6/6.D0/,C1A/0.355028053887817D0/,C2A/0.258819403792807D0/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IQTST=1
      RMINN=RMIN-H
c** Start by searching for third turning point.
      J=N
      IF(V(N).GT.E) GO TO 22
      DO  I=2,N
          J=J-1
          IF(V(J).GT.E) GO TO 10
          ENDDO
      GO TO 14
c** Check that there is a classically allowed region inside this point
c  and determine height of barrier maximum.
   10 II=J
      VMX=DSOC
      DO  I=2,J
          II=II-1
          IF(V(II).LE.E) GO TO 16
          IF(V(II).GT.VMX) VMX=V(II)
          ENDDO
c** Energy too high ... find no more than one turning point.
   14 XJ1=RMINN+J*H
c ... Search outward for barrier height to facilitate energy correction
      IF(J.EQ.1) J= 2
      K=J-1
      DO  I=J,N
          IF(V(I).GT.V(K)) GO TO 120
          K=I
          ENDDO
      VMX=V(K)
      GO TO 130
  120 K=K+2
      J=K-1
      DO  I=K,N
          IF(V(I).LT.V(J)) GO TO 126
          J=I
          ENDDO
  126 VMX=V(J)
  130 VMXPR=VMX/BFCT
      IF(IWR.NE.0) WRITE(6,608) JROT,EO,VMXPR,XJ1
      ITP3= J
      IQTST=-1
      GO TO 100
   16 ITP3= J+1
c** ITP3 is the first mesh point outside classically forbidden region
      GB=V(ITP3)-E
      GI=V(ITP3-1)-E
      FJ=GI/(GI-GB)
c** Treat quasibound levels as bound using outer boundary condition
c  of Airy function at third turning point ... as discussed by
c  R.J.Le Roy and R.B.Bernstein  in  J.Chem.Phys. 54,5114(1971).
      SL=(GI-GB)**(Z1/Z3)/H
      IF((SL*H).LT.Z1) THEN
          A1=GI/(SL*H)**2
          A2=GB/(SL*H)**2
          A13=A1*A1*A1
          A23=A2*A2*A2
          FIA=Z1+A13*(A13*(A13+72.D0)+2160.D0)/12960.D0
          GIA=A1+A1*A13*(A13*(A13+90.D0)+3780.D0)/45360.D0
          FBA=Z1+A23*(A23*(A23+72.D0)+2160.D0)/12960.D0
          GBA=A2+A2*A23*(A23*(A23+90.D0)+3780.D0)/45360.D0
c** Airy function  Bi(X)  at points straddling 3-rd turning point
          SI=C1A*FIA+C2A*GIA
          SB=C1A*FBA+C2A*GBA
          GO TO 100
          ENDIF
c** If Airy function expansion unreliable, use zero slope at third
c  turning point as quasibound outer boundary condition.
      DF=GI-GB
      SI=Z1+DF*FJ**3/Z6
      SB=Z1-DF*(Z1-FJ)**3/Z6
      IF(IWR.NE.0) WRITE(6,606) KV,JROT,EO,IT
      GO TO 100
c** If 3-rd turning point beyond range start with WKB wave function
c  at end of range.
   22 IF(IWR.NE.0) WRITE(6,607) JROT,EO
      ITP3= N
      IQTST=0
      GB=V(ITP3)-E
      GI=V(ITP3-1)-E
      VMX=V(ITP3)
      II=ITP3
      DO  I=2,ITP3
          II=II-1
          IF(V(II).LT.VMX) GO TO 100
          VMX=V(II)
          ENDDO
      IF(IWR.NE.0) WRITE(6,604)
c** End of quasibound level initialization schemes.
      IQTST=-9
  100 RETURN
  604 FORMAT(" **** QBOUND doesn't work ... no classically allowed regio
     1n accessible at this energy.")
  606 FORMAT(' *** CAUTION ***  v=',I3,'   J=',I3,'   E=',1PD13.6,
     1 '   IT=',I2/5x,'Airy initialization unstable so use  zero slope',
     2 'at  R(3-rd)' )
  607 FORMAT(' *** For  J=',I3,'  E=',F9.2,
     1  '  R(3-rd) > RMAX  & E < V(N)  so try WKB B.C. @ RMAX')
  608 FORMAT(' For J=',I3,'  ETRY=',F11.4,' > VMAX=',F11.4,
     1  '  find onee turn point:  R=',F6.2)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c** Subroutine to calculates quasibound level tunneling lifetime/width
c** For relevant theory see Le Roy & Liu [J.Chem.Phys.69,3622-31(1978)]
c  and Connor & Smith [Mol.Phys. 43, 397 (1981)].
c** Final level width calculation from Eq.(4.5) of Connor & Smith.
c-----------------------------------------------------------------------
      SUBROUTINE WIDTH(KV,JROT,E,EO,DSOC,V,S,VMX,RMIN,H,BFCT,IWR,ITP1,
     1  ITP3,INNER,N,GAMA)
c++ "WIDTH" calls subroutine "LEVQAD" ++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,IMM,INNER,IRM,ITP1,ITP1P,ITP1P1,ITP2,ITP2M,ITP2M2,
     1         ITP2P1,ITP2P2,ITP3,IWR,JROT,KV,KVI,KVO,
     2         M,M2,N,NN,NST
      REAL*8  AA,ANS1,ANS2,ARG,BFCT,COR,
     1        D1,D2,D3,DFI,DSGB,DSGN,DSOC,DWEB,OMEGJC,
     2        E,EO,EMSC,EMV,FJNLC,G1,G2,G3,GA,GAMA,GAMALG,
     3        H,H2,HBW,HBWB,PI,PMX,RMIN,RMINN,RMX,RT,RT1,RT2,
     4        S(N),SM,TAU,TAULG,TI,TUN0,TNUM,U1,U2,V(N),VMAX,VMX,
     7        XJ,XX,Z0,Z1,Z2,Z4,ZH
      CHARACTER*5 LWELL(2)
      DATA Z0/0.D0/,ZH/0.5D0/,Z1/1.D0/,Z2/2.D0/,Z4/4.D0/,
     1  PI/3.141592653589793D0/
      DATA LWELL/'INNER','OUTER'/
      RMINN=RMIN-H
      H2=H*H
c** ITP1 is first mesh point to right of innermost turning point.
   40 ITP1P=ITP1+1
      ITP1P1=ITP1P+1
      IRM=ITP1-1
c** Calculate JWKB tunneling probability from quadrature over barrier
c** First must locate 2-nd turning point.
      DO  I=ITP1P1,ITP3
          ITP2=I
          IF(V(I).GT.E) GO TO 202
          ENDDO
      GAMA=Z0
      GO TO 250
  202 ITP2P1=ITP2+1
      ITP2P2=ITP2+2
c** ITP2M is the last mesh point before the 2-nd turning point.
      ITP2M=ITP2-1
      ITP2M2=ITP2-2
      G1=V(ITP2M)-E
      G2=V(ITP2)-E
      GA=V(ITP2P1)-E
c** Quadrature over barrier starts here.
      CALL LEVQAD(G1,G2,GA,H,RT,ANS1,ANS2)
      SM=ANS2/H
      IF(GA.LT.Z0) GO TO 218
      SM=SM+ZH*DSQRT(GA)
      PMX=VMX
      M2=ITP2P2
  204 DO  I=M2,ITP3
          M=I
          GA=V(I)-E
          IF(V(I).GT.PMX) PMX=V(I)
          IF(GA.LT.Z0) GO TO 210
          SM=SM+DSQRT(GA)
          ENDDO
      IF(V(M).GT.V(M-1)) THEN
          IF(IWR.NE.0) WRITE(6,602) KV,JROT
          GO TO 250
          ENDIF
      RMX=RMINN+M*H
      U1=DSQRT(GA/(V(M)-DSOC))
      U2=DSQRT((E-DSOC)/(V(M)-DSOC))
      SM=SM-ZH*DSQRT(GA) + (DLOG((Z1+U1)/U2)-U1)*RMX*DSQRT(V(M)-DSOC)/H
      XJ=(DSQRT(Z1+Z4*(V(M)-DSOC)*(RMX/H)**2)-Z1)/Z2
      IF(IWR.NE.0) WRITE(6,603) JROT,EO,XJ,RMX
      GO TO 218
  210 IF(M.LT.ITP3) THEN
c** If encounter a double-humped barrier, take care here.
          IF(IWR.NE.0) WRITE(6,609) KV,JROT,EO,M
          KVO=0
          DSGN=DSIGN(Z1,S(M-1))
c** Find the effective quantum number for the outer well
          DO  I=M,ITP3
              DSGB=DSGN
              DSGN=DSIGN(Z1,S(I))
              IF((DSGN*DSGB).LT.Z0) KVO=KVO+1
              ENDDO
          KVI=KV-KVO
          IF(INNER.EQ.0) THEN
c** For levels of outer well, get correct width by changing ITP1
              ITP1=M
              IF(IWR.GT.0) WRITE(6,610) KVO,LWELL(2)
              GO TO 40
              ENDIF
          IF(IWR.GT.0) WRITE(6,610) KVI,LWELL(1)
c** For "inner-well" levels, locate outer barrier
          DO  I=M,ITP3
              M2=I
              GA=V(I)-E
              IF(GA.GE.Z0) GO TO 204
              ENDDO
          GO TO 218
          ENDIF
      G3=V(M-2)-E
      G2=V(M-1)-E
      CALL LEVQAD(GA,G2,G3,H,RT,ANS1,ANS2)
      SM= SM- ZH*DSQRT(G3)-DSQRT(G2) + ANS2/H
  218 EMSC= -SM/PI
      IF(INNER.NE.0) VMX= PMX
      VMAX= VMX/BFCT
c** Tunneling factors calculated here ** TUN0 is simple WKB result
c  as in Child's eqs.(57c) & (59).
      TUN0= Z0
      IF(DABS(EMSC).LT.25.D0) TUN0= ZH*DEXP(Z2*PI*EMSC)
c ... for permeability calculate Connor-Smith's Eq.(3.7) \omega=OMEGJC
      FJNLC= DSQRT(Z1+ Z2*TUN0)
      OMEGJC= FJNLC- Z1
      IF(TUN0.LT.1.D-6) OMEGJC= TUN0
      OMEGJC= OMEGJC/(FJNLC+ 1.d0)
c** Quadrature for JWKB calculation of vibrational spacing in well HBW
      D1=E-V(IRM)
      D2=E-V(ITP1)
      D3=E-V(ITP1P)
      CALL LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
      RT1=RT
      SM=ANS1/H
      IF(D3.LT.Z0) GO TO 228
      SM=SM+ZH/DSQRT(D3)
      DO  I=ITP1P1,ITP2M2
          IMM=I
          EMV=E-V(I)
          IF(EMV.LT.Z0) GO TO 222
          SM=SM+Z1/DSQRT(EMV)
          ENDDO
      D3=E-V(ITP2M2)
      D2=E-V(ITP2M)
      D1=E-V(ITP2)
      GO TO 226
c** If encounter a double-minimum well, take care here.
  222 D1=EMV
      D2=E-V(IMM-1)
      D3=E-V(IMM-2)
      IF(IWR.NE.0) WRITE(6,605) KV,JROT,EO
  226 CALL LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
      RT2=RT
      SM=SM-ZH/DSQRT(D3) + ANS1/H
c** Get HBW in same energy units (1/cm) associated with BFCT
  228 HBW=Z2*PI/(BFCT*SM)
c** HBW fix up suggested by Child uses his eqs.(48)&(62) for HBW
      AA= DLOG(DABS(EMSC))
      TNUM= Z1/EMSC
c** Derivative of complex gamma function argument calculated as
c  per eq.(6.1.27) in Abramowitz and Stegun.
      ARG= Z2/((ZH/EMSC)**2+Z1)
      NST= DABS(EMSC)*1.D2
      NST= MAX0(NST,4)
      DO  I= 1,NST
          NN= I
          XX= (ZH+NN)/EMSC
          TI= TNUM/(XX*(XX**2+Z1))
          ARG= ARG+TI
          IF(DABS(TI).LT.1.D-10) GO TO 233
          ENDDO
  233 COR= ZH*(EMSC/(NN+Z1))**2
      ARG= ARG+COR-COR**2
      DWEB= (EO-VMAX)*BFCT/(H2*EMSC)
      DFI= (AA + 1.96351002602134D0-ARG)*BFCT/(H2*DWEB)
      HBWB= Z1/(Z1/HBW + DFI/(Z2*PI))
c** Width from formula (4.5) of  Connor & Smith, Mol.Phys.43,397(1981)
c [neglect time delay integral past barrier in their Eq.(4.16)].
      IF(EMSC.GT.-25.D0) THEN
          GAMA = (HBWB/(Z2*PI))*Z4*OMEGJC
          TAU= 0.D0
          IF(GAMA.GT.1.D-60) TAU= 5.308837457D-12/GAMA
c** GAM0 = TUN0*HBW/PI  is the simple WKB width GAMMA(0) discussed by
c  Le Roy & Liu in J.C.P.69,3622(1978).
          IF(IWR.GT.0) WRITE(6,601) TAU,GAMA,HBWB,VMAX
          GO TO 250
          ENDIF
      GAMALG= DLOG10(HBWB/(Z2*PI))+Z2*PI*EMSC/2.302585093D0
      TAULG= DLOG10(5.308837457D-12)-GAMALG
      IF(IWR.GT.0) WRITE(6,611) TAULG,GAMALG,HBWB,VMAX
  250 RETURN
  601 FORMAT('    Lifetime=',1PD10.3,'(s)   Width=',D10.3,'   dG/dv=',
     1 0PF7.2,'   V(max)=',F9.2)
  602 FORMAT(' *** WARNING ***  For   v =',I3,'   J =',I3,'   cannot cal
     1culate width since barrier maximum beyond range')
  603 FORMAT(' *** For  J=',I3,'  E=',F9.2,'  R(3-rd) beyond range so tu
     1nneling calculation uses'/8X,'pure centrifugal potential with  J(a
     2pp)=',F7.2,'  for  R > R(max)=',F7.2)
  605 FORMAT(' **** CAUTION *** Width estimate only qualitative, as have
     1 a double-minimum well for   E(v=',I3,', J=',I3,')=',F15.7/15X,
     2 'a more stable result may be obtained by searching for the quasib
     3ound levels using option: INNER > 0 .')
  609 FORMAT(' *** CAUTION - Permeability estimate not exact as have a d
     1ouble-humped barrier:  E(v=',I3,', J=',I3,') =',G15.8,I6)
  610 FORMAT(16X,'(NOTE: this has the node count of a   v=',I3,2X,A5,
     1 '-well level')
  611 FORMAT(12X,'Log10(lifetime/sec)=',F10.5,' ;   Log10(width/cm-1)=',
     1 F10.5,'   Spacing=',G12.5,'   V(max)=',G14.7,'(cm-1)')
      END
c**********************************************************************
      SUBROUTINE LEVQAD(Y1,Y2,Y3,H,RT,ANS1,ANS2)
c** Subroutine "LEVQAD" fits quadratic  Y = A + B*X + C*X**2  through
c  function values  Y1, Y2, Y3  at equally spaced points separated by
c  distance H, where  Y1 < 0  and (Y2,Y3 .ge.0), locates the function
c  zero (at RT, relative to  X1 < X2 = 0) between points X1 & X2, and
c  evaluates the integral from RT to R3 of   1/sqrt(Y)  , called
c  ANS1, and the integral (same range) of  sqrt(Y) , which is ANS2
c** Alternately, if Y1 & Y3 both  < 0  and only the middle point
c  Y2.ge.0 ,   fit the points to:  Y = A - B*(X-X0)**2 , locate the
c  turning points between which  Y(X) > 0  and evaluate these integrals
c  on this interval.  *************************************************
c----------------------------------------------------------------------
      REAL*8  A,ANS1,ANS2,B,C,CQ,H,HPI,R1,R2,RCQ,RR,RT,SL3,SLT,
     1        X0,Y1,Y2,Y3,Z0,Z1,Z2,Z4,ZT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF((Y1.GE.0).OR.(Y2.LT.0)) GO TO 99
      DATA Z0/0.D0/,Z1/1.D0/,Z2/2.D0/,Z4/4.D0/,HPI/1.570796326794896D0/
      IF(Y3.LT.Z0) GO TO 50
c** Here treat case where both 'Y2' & 'Y3' are positive
      IF(DABS((Y2-Y1)/(Y3-Y2) -1.D0).LT.1.d-10) THEN
c ... special case of true (to 1/10^10) linearity ...
          RT= -H*Y2/(Y2-Y1)
          ANS1= 2.d0*(H-RT)/DSQRT(Y3)
          ANS2= ANS1*Y3/3.D0
          RETURN
          ENDIF
      C=(Y3-Z2*Y2+Y1)/(Z2*H*H)
      B=(Y3-Y2)/H-C*H
      A=Y2
      CQ=B**2-Z4*A*C
      RCQ=DSQRT(CQ)
      R1=(-B-RCQ)/(Z2*C)
      R2=R1+RCQ/C
      IF((R2.LE.Z0).AND.(R2.GE.-H)) RT=R2
      IF((R1.LE.Z0).AND.(R1.GE.-H)) RT=R1
      SL3=Z2*C*H+B
      SLT=Z2*C*RT+B
      IF(C.LT.Z0) GO TO 10
      ANS1=DLOG((Z2*DSQRT(C*Y3)+SL3)/SLT)/DSQRT(C)
      GO TO 20
   10 ANS1=-(DASIN(SL3/RCQ)-DSIGN(HPI,SLT))/DSQRT(-C)
   20 ANS2=(SL3*DSQRT(Y3)-CQ*ANS1/Z2)/(Z4*C)
      IF(RT.GE.H) WRITE(6,601) H,R1,R2
  601 FORMAT(' *** CAUTION *** in LEVQAD, turning point not between poin
     1ts 1 & 2.   H =',F9.6,'   R1 =',F9.6,'   R2 =',F9.6)
      RETURN
c** Here treat case when only 'Y2' is non-negative
   50 RR=(Y2-Y1)/(Y2-Y3)
      X0=H*(RR-Z1)/((RR+Z1)*Z2)
      B=(Y2-Y1)/(H*(Z2*X0+H))
      A=Y2+B*X0**2
      ZT=DSQRT(A/B)
      RT=X0-ZT
      ANS1=Z2*HPI/DSQRT(B)
      ANS2=ANS1*A/Z2
      RETURN
   99 WRITE(6,602) Y1,Y2
  602 FORMAT(' *** ERROR in LEVQAD *** No turning point between 1-st two
     1 points as   Y1=',D10.3,'   Y2=',D10.3)
      ANS1=Z0
      ANS2=Z0
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE CDJOEL(EO,NBEG,NEND,BvWN,RH,WARN,V,WF0,RM2,RCNST)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Subroutine solving the linear inhomogeneous differential equations
c  formulated by J.M. Hutson [J.Phys.B14, 851 (1982)] for treating
c  centrifugal distortion as a perturbation, to determine centrifugal
c  distortion constants of a diatomic molecule.  Uses the algorithm of
c  J. Tellinghuisen [J.Mol.Spectrosc. 122, 455 (1987)].  The current
c  version calculates Bv, Dv, Hv, Lv, Mv, Nv and Ov and writes them out,
c  but does not return values to the calling program.
c
c** On entry:   EO    is the eigenvalue (in units [cm-1])
c               NBEG & NEND  the mesh point range over which the input
c wavefunction  WF0  (in units 1/sqrt(Ang))  has non-negligible values
c               BvWn  is the numerical factor (hbar^2/2mu) [cm-1 Ang^2]
c               RH    is the integration stepsize (in units [Ang])
c               WARN  is an integer flag: > 0 print internal warnings,
c               V(i)  is the effective potential (including centrifugal
c                     term if calculation performed at  J > 0) in
c                     'internal' units, including the factor  RH**2/BvWN
c               RM2(i) is the array  1/(distance**2) in units [1/Ang**2]
c** On exit:    RCNST(i)  is the set of 7 rotational constants: Bv, -Dv,
c                       Hv, Lv, Mv, Nv & Ov
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1994  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Authors: R.J. Le Roy & J. Tellinghuisen         Version of 30/09/1999
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Dimension:  potential arrays  and  vib. level arrays.
      INTEGER NDMINT
      PARAMETER (NDMINT= 40001)
      INTEGER I,M,IPASS,M1,M2,NBEG,NEND,WARN
      REAL*8 V(NEND),WF0(NEND),RM2(NEND),P(NDMINT),WF1(NDMINT),
     1                                            WF2(NDMINT),RCNST(7)
      REAL*8 BvWN,DV,DVV,HVV,HV2,LVV,LV2,MVV,MV2,NVV,OVV,EO,E,RH,RHSQ,
     1  ZTW,AR,R2IN,G2,G3,P0,P1,P2,P3,PI,PIF,PRS,PRT,V1,V2,V3,Y1,Y2,Y3,
     2  TSTHv,TSTLv,TSTMv,AMB,AMB1,AMB2,
     3  OV,OV01,OV02,OV03,OV11,OV12,OV13,OV22,OV23,OV33,
     4  PER01,PER02,PER03,PER11,PER12,PER13,PER22,PER23,PER33
c
      IF(NEND.GT.NDMINT) THEN
          WRITE(6,602) NEND,NDMINT
          RETURN
          ENDIF
      ZTW= 1.D0/12.d0
      RHSQ = RH*RH
      DV = RHSQ/12.D0
      E= EO*RHSQ/BvWN
      IPASS = 1
      OV01 = 0.D0
      OV02 = 0.D0
      OV03 = 0.D0
      OV11 = 0.D0
      OV22 = 0.D0
      OV12 = 0.D0
      OV33 = 0.D0
      OV23 = 0.D0
      OV13 = 0.D0
      PER01 = 0.D0
      PER02 = 0.D0
      PER03 = 0.D0
      PER11 = 0.D0
      PER12 = 0.D0
      PER13 = 0.D0
      PER22 = 0.D0
      PER23 = 0.D0
      PER33 = 0.D0
c** First, calculate the expectation value of  1/r**2  and hence Bv
      R2IN= 0.5D0*(RM2(NBEG)*WF0(NBEG)**2 + RM2(NEND)*WF0(NEND)**2)
      DO   I= NBEG+1, NEND-1
         R2IN= R2IN+ RM2(I)*WF0(I)**2
         ENDDO
      R2IN = R2IN*RH
      RCNST(1)= R2IN*BvWN
c
c** On First pass  IPASS=1  and calculate first-order wavefx., Dv & Hv
c  On second pass  IPASS=2  and calculate second-order wavefx., Lv & Mv
c  On third pass   IPASS=3  and calculate third-order wavefx., Nv & Ov
c
   10 P1= 0.D0
      P2= 0.D0
c
c     P1= WF0(NEND)
c     P2= WF0(NEND-1)
c
      P(NEND) = P1
      P(NEND-1) = P2
      V1 = V(NEND) - E
      V2 = V(NEND-1) - E
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*(RM2(NEND) - R2IN)*WF0(NEND)
          G2 = (RM2(NEND-1) - R2IN)*WF0(NEND-1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NEND) - R2IN)*WF1(NEND)
     1                                                - DVV*WF0(NEND))
          G2 = (RM2(NEND-1) - R2IN)*WF1(NEND-1) - DVV*WF0(NEND-1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NEND) - R2IN)*WF2(NEND)
     1                                - DVV*WF1(NEND) - HVV*WF0(NEND))
          G2 = (RM2(NEND-1) - R2IN)*WF2(NEND-1) - DVV*WF1(NEND-1)
     1                                               - HVV*WF0(NEND-1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      M= NEND-1
c** Now - integrate inward from outer end of range
      DO  I = NBEG+2,NEND
          M = M-1
          Y3 = Y2 + Y2 - Y1 + RHSQ*G2 + V2*P2
          IF(IPASS.EQ.1) G3 = (RM2(M) - R2IN)*WF0(M)
          IF(IPASS.EQ.2) G3 = (RM2(M) - R2IN)*WF1(M) - DVV*WF0(M)
          IF(IPASS.EQ.3) G3 = (RM2(M) - R2IN)*WF2(M) - DVV*WF1(M)
     1                                                    - HVV*WF0(M)
          V3 = V(M) - E
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          IF(V3.LT.0.D0)  GO TO 32
          P(M) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          ENDDO
      GO TO 90
c** Escaped loop at outer turning point:  initialize outward integration
   32 PRS = P3
      PRT = P(M+1)
      P1 = 0.D0
      P2 = 0.D0
c
c     P1 = WF0(NBEG)
c     P2 = WF0(NBEG+1)
c
      P(NBEG) = P1
      P(NBEG+1) = P2
      V1 = V(NBEG) - E
      V2 = V(NBEG+1) - E
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*(RM2(NBEG) - R2IN)*WF0(NBEG)
          G2 = (RM2(NBEG+1) - R2IN)*WF0(NBEG+1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NBEG) - R2IN)*WF1(NBEG)
     1                                                - DVV*WF0(NEND))
          G2 = (RM2(NBEG+1) - R2IN)*WF1(NBEG+1) - DVV*WF0(NBEG+1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NBEG) - R2IN)*WF2(NBEG)
     1                                - DVV*WF1(NEND) - HVV*WF0(NEND))
          G2 = (RM2(NBEG+1) - R2IN)*WF2(NBEG+1) - DVV*WF1(NBEG+1)
     2                                               - HVV*WF0(NBEG+1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      AR = 0.D0
      M1 = M+1
c** Now ... integrate outward from inner end of range
      DO  I = NBEG+2,M1
          Y3 = Y2 + Y2 - Y1 + RHSQ*G2 + V2*P2
          P0 = WF0(I)
          IF(IPASS.EQ.1) G3 = (RM2(I) - R2IN)*P0
          IF(IPASS.EQ.2) G3 = (RM2(I)-R2IN)*WF1(I) - DVV*P0
          IF(IPASS.EQ.3) G3 = (RM2(I)-R2IN)*WF2(I) - DVV*WF1(I) - HVV*P0
          V3 = V(I) - E
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          P(I) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          AR = AR + P0*P3
          ENDDO
c** Average for 2 adjacent mesh points to get Joel's "(a-b)"
      AMB2 = (P3-PRT)/P0
      AMB1 = (P(M)-PRS)/WF0(M)
      AMB = (AMB1+AMB2)*0.5D0
      M2 = M+2
c** Find the rest of the overlap with zero-th order solution ...
      DO  I = M2,NEND
          P0 = WF0(I)
          PI = P(I) + AMB*P0
          P(I) = PI
          AR = AR + PI*P0
          ENDDO
      OV = AR*RH
      DO  I = NBEG,NEND
          P0 = WF0(I)
c ... and project out contribution of zero'th-order part of solution
          PI = P(I) - OV*P0
          PIF = PI*RM2(I)
          IF(IPASS.EQ.1) THEN
c** Now - on first pass accumulate integrals for Dv and Hv
              WF1(I) = PI
              OV01 = OV01 + PI*P0
              OV11 = OV11 + PI*PI
              PER01 = PER01 + PIF*P0
              PER11 = PER11 + PI*PIF
            ELSEIF(IPASS.EQ.2) THEN
c ... and on next pass, accumulate integrals for Lv and Mv
              WF2(I) = PI
              P1 = WF1(I)
              OV02 = OV02 + PI*P0
              OV12 = OV12 + PI*P1
              OV22 = OV22 + PI*PI
              PER02 = PER02 + PIF*P0
              PER12 = PER12 + PIF*P1
              PER22 = PER22 + PI*PIF
            ELSEIF(IPASS.EQ.3) THEN
c ... and on next pass, accumulate integrals for Nv and Ov
              P1 = WF1(I)
              P2 = WF2(I)
              OV03 = OV03 + PI*P0
              OV13 = OV13 + PI*P1
              OV23 = OV23 + PI*P2
              OV33 = OV33 + PI*PI
              PER03 = PER03 + PIF*P0
              PER13 = PER13 + PIF*P1
              PER23 = PER23 + PIF*P2
              PER33 = PER33 + PIF*PI
            ENDIF
          ENDDO
      IF(IPASS.EQ.1) THEN
          DVV = RH*PER01
          HVV = RH*(PER11 - R2IN*OV11)
          IPASS = 2
          RCNST(2) = DVV*BvWN
          RCNST(3) = HVV*BvWn
          GO TO 10
        ELSEIF(IPASS.EQ.2) THEN
          HV2 = RH*PER02*BvWN
          LVV = RH*(PER12 - R2IN*OV12 - DVV*OV11)
          MVV = RH*(PER22 - R2IN*OV22 - 2.D0*DVV*OV12 - HVV*OV11)
          IPASS = 3
          RCNST(4) = LVV*BvWN
          RCNST(5) = MVV*BvWN
          GO TO 10
        ELSEIF(IPASS.EQ.3) THEN
          LV2 = RH*PER03*BvWN
          MV2 = RH*(PER13 - R2IN*OV13 - DVV*OV12 - HVV*OV11)*BvWN
          NVV = RH*(PER23 - R2IN*OV23 - DVV*(OV13 + OV22)
     1                                     - 2.D0*HVV*OV12 - LVV*OV11)
          OVV = RH*(PER33 - R2IN*OV33 - 2.D0*DVV*OV23
     1             - HVV*(2.D0*OV13+ OV22) - 2.D0*LVV*OV12 - MVV*OV11)
          RCNST(6) = NVV*BvWN
          RCNST(7) = OVV*BvWN
        ENDIF
      IF(WARN.GT.0) THEN
          IF(DMAX1(DABS(OV01),DABS(OV02),DABS(OV01)).GT.1.D-9)
     1                                     WRITE(6,604) OV01,OV02,OV03
          TSTHV= dabs(RCNST(3)/HV2-1.D0)
          TSTLV= dabs(RCNST(4)/LV2-1.D0)
          TSTMV= dabs(RCNST(5)/MV2-1.D0)
          IF(DMAX1(TSTHV,TSTLV,TSTMV).GT.1.d-5)
     1                                  WRITE(6,603) TSTHV,TSTLV,TSTMV
          ENDIF
      RETURN
   90 WRITE(6,601) EO
      RETURN
  601 FORMAT(' *** ERROR in CDJOEL *** for input energy  E =',f12.4,
     1   '  never reach outer turning point')
  602 FORMAT(/' *** Dimensioning PROBLEM in CDJOEL ***   NEND=',i6,
     1  ' > NDMINT=',i6)
  603 FORMAT(' ** CAUTION ** Comparison tests for Hv, Lv & Mv give:',
     1 3(1Pd9.1))
  604 FORMAT(' ** CAUTION ** CDJOEL orthogonality tests OV01,OV02 & OV03
     1:',3(1Pd9.1))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE NLLSSRR(NDATA,NPTOT,NPMAX,IROUND,NGPRND,LPRINT,YO,YU,
     1                                 YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c**  Program for performing linear or non-linear least-squares fits and
c  (if desired) automatically using sequential rounding and refitting
c  to minimize the numbers of parameter digits which must be quoted [see
c  R.J. Le Roy, J.Mol.Spectrosc. 191, 223-231 (1998)].          16/11/00
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c             COPYRIGHT 1998-2000  by  Robert J. Le Roy                +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Program uses orthogonal decomposition of the "design" (partial
c  derivative) matrix for the core locally linear (steepest descent)
c  step, following a method introduced (to me) by Dr. Michael Dulick.
c** If no parameters are free (NPTOT=0), simply return RMS(residuals) as
c  calculated from the input parameter values {PV(j)}.
c** A user MUST SUPPLY subroutine  DYIDPJ  to generate the predicted
c  value of each datum and the partial derivatives of each datum w.r.t.
c  each parameter (see below) from the current trial parameters.
c
c** On entry:
c    NDATA  is the number of data to be fitted
c    NPTOT  the number of parameters to be varied (.le.NPMAX).
c           If NPTOT.le.0 , assume  YD(i)=YO(i)  and calculate the (RMS
c           dimensionless deviation)=DSE  from them & YU(i)
c    NPMAX is the maximum number of free parameters allowed by current
c          external array sizes.  Should set internal NPINTMX = NPMAX
c          (may be freely changed by the user).
c    IROUND .ne. 0  causes Sequential Rounding & Refitting to be
c             performed, with each parameter being rounded at the
c            |IROUND|'th sig. digit of its local incertainty.
c        > 0  rounding selects in turn remaining parameter with largest
c             relative uncertainy
c        < 0  round parameters sequentially from last to first
c        = 0  simply stops after full convergence (without rounding).
c    NGPRND in the unusual case when one has MANY parameters and wants
c         to round a number off all at once, setting NGPRDN > 1  causes
c         the last NGPRND parameters to be rounded in the initial step.
c         Use only if necessary:  should NORMALLY set  NGPRND.le.1 .
c    LPRINT  specifies the level of printing inside NLLSSRR
c          if: =  0, no print except for failed convergence.
c               < 0  only converged, unrounded parameters, PU & PS's
c              >= 1  print converged parameters, PU & PS's
c              >= 2  also print parameter change each rounding step
c              >= 3  also indicate nature of convergence
c              >= 4  also print convergence tests on each cycle
c              >= 5  also parameters changes & uncertainties, each cycle
c              >= 6  also print correlation matrix on each cycle
c    YO(i)  are the NDATA 'observed' data to be fitted 
c    YU(i)  are the uncertainties in these YO(i) values
c    PV(j)  are initial trial parameter values (for non-linear fits); 
c           should be set at zero for initially undefined parameters.
c
c** On Exit:   
c    YD(i)  is the array of differences  [Ycalc(i) - YO(i)]
c    PV(j)  are the final converged parameter values
c    PU(j)  are 95% confidence limit uncertainties in the PV(j)'s
c    PS(j)  are 'parameter sensitivities' for the PV(j)'s, defined such
c           that the RMS displacement of predicted data  due to rounding
c           off parameter-j by PS(j) is .le. DSE/10*NPTOT
c    CM(j,k)  is the correlation matrix obtained by normalizing variance
c           /covariance matrix:  CM(j,k) = CM(j,k)/SQRT[CM(j,j)*CM(k,k)]
c    TSTPS = max{|delta[PV(j)]/PS(j)|}  is the parameter sensitivity
c          convergence test:  delta[PV(j)] is last change in parameter-j
c    TSTPU = max{|delta[PV(j)]/PU(j)|}  is the parameter uncertainty
c          convergence test:  delta[PV(j)] is last change in parameter-j
c    DSE    is the predicted (dimensionless) standard error of the fit
c
c  NOTE that the squared 95% confidence limit uncertainty in a property
c  F({PV(j)}) defined in terms of the fitted parameters {PV(j)} (where
c  the L.H.S. involves  [row]*[matrix]*[column]  multiplication) is:
c  [D(F)]^2 = [PU(1)*dF/dPV(1), PU(2)*dF/dPV(2), ...]*[CM(j,k)]*
c                              [PU(2)*dF/dPV(1), PU(2)*dF/dPV(2), ...]
c
c** Externally dimension:  YO, YU and YD  .ge. NDATA
c             PV, PU  and  PS  .ge.  NPTOT (say as NPMAX),
c             CM   as a square matrix with column & row length  NPMAX
c***********************************************************************
      INTEGER NPINTMX
      PARAMETER (NPINTMX=2000)
      INTEGER I,J,K,L,IDF,ITER,NITER,IROUND,JROUND,LPRINT,NDATA,
     1        NGPRND,NPTOT,NPARM,NPMAX,QUIT,KFIX,JFIX,IFXP(NPINTMX)
      REAL*8  YO(NDATA), YU(NDATA), YD(NDATA), PV(NPTOT), PU(NPTOT),
     1  PS(NPTOT),PSS(NPINTMX),PC(NPINTMX),PX(NPINTMX),PY(NPINTMX),
     2  CM(NPMAX,NPMAX), F95(10),
     3  RMSR, RMSRB, DSE, TSTPS, TSTPSB, TSTPU, TFACT, S, UU
      DATA F95/12.7062D0,4.3027D0,3.1824D0,2.7764D0,2.5706D0,2.4469D0,
     1  2.3646D0,2.3060D0,2.2622D0,2.2281D0/
c
      IF((NPTOT.GT.NPMAX).OR.(NPTOT.GT.NPINTMX)
     1                    .OR.(NPTOT.GT.NDATA)) THEN
c** If array dimensioning inadequate, print warning & then STOP
          WRITE(6,602) NPTOT,NPINTMX,NPMAX,NDATA
          STOP
          ENDIF
      NPARM= NPTOT
      TSTPS= 0.d0
      RMSR= 0.d0
      NITER= 0
      QUIT= 0
      DO J= 1,NPTOT
          PS(J)= 0.d0
          IFXP(J)= 0
          ENDDO
      JROUND= IABS(IROUND)
c=======================================================================
c** Beginning of loop to perform rounding (if desired).  NOTE that in
c  sequential rounding, NPARM is the current (iteratively shrinking)
c  number of free parameters.
    6 IF(NPARM.GT.0) TSTPS= 9.d99
c** TFACT  is 95% student t-value for (NDATA-NPARM) degrees of freedom.
c [Approximate expression for (NDATA-NPARM).GT.10 accurate to ca. 0.002]
      TFACT= 0.D0
      IF(NDATA.GT.NPARM) THEN
          IDF= NDATA-NPARM
          IF(IDF.GT.10) TFACT= 1.960D0*DEXP(1.265D0/DBLE(IDF))
          IF(IDF.LE.10) TFACT= F95(IDF)
        ELSE
          TFACT= 0.D0
        ENDIF
c======================================================================
c** Begin iterative convergence loop:  try for up to 20 cycles
      DO 50 ITER= 1,20
          NITER= NITER+ 1
          DSE= 0.d0
          TSTPSB= TSTPS
          RMSRB= RMSR
c** Zero out various arrays
          IF(NPARM.GT.0) THEN
              DO  I = 1,NPARM
c** PSS is the array of Saved Parameter Sensitivities from previous run
c   to be carried into dyidpj subroutine - used in predicting increment
c   for derivatives by differences.
                  PSS(I)= PS(I)
                  PS(I) = 0.D0
                  PU(I) = 0.D0
                  PX(I) = 0.D0
                  PY(I) = 0.D0
                  DO  J = 1,NPARM
                      CM(I,J) = 0.D0
                      ENDDO
                  ENDDO
              ENDIF
c
c========Beginning of core linear least-squares step====================
c
c** Begin by forming the Jacobian Matrix from partial derivative matrix
          DO  I = 1,NDATA
c** User-supplied subroutine DYIDPJ uses current (trial) parameter
c  values {PV} to generate predicted datum # I [y(calc;I)=UU] and its
c  partial derivatives w.r.t. each of the parameters, returning the
c  latter in 1-D array PC.  See dummy sample version at end of listing.
c* [NOTE: if desired, could write DYIDPJ such that the y(calc) values
c     and derivatives for all data are prepared at the same time (when
c     I=1), but only returned here one datum at a time (for I > 1). 
c     However, this would be inappropriate for very large data sets.]
              CALL DYIDPJ(I,NDATA,NPTOT,UU,PV,PC,PSS,RMSR)
              IF((NPARM.LT.NPTOT).AND.(IROUND.GT.0)) THEN
c** For sequential rounding, collapse partial derivative array here
                  DO  J= NPTOT,1,-1
                      IF((IFXP(J).GT.0).AND.(J.LT.NPTOT)) THEN
                          DO  K= J,NPTOT-1
                              PC(K)= PC(K+1)
                              ENDDO
                          PC(NPTOT)= 0.d0
                          ENDIF
                      ENDDO
                  ENDIF
              S = 1.D0 / YU(I)
              YD(I)= UU - YO(I)
              UU = - YD(I) * S
              DSE= DSE+ UU*UU
              IF(NPARM.GT.0) THEN
                  DO  J = 1,NPARM
                      PC(J) = PC(J)*S
                      PS(J) = PS(J)+ PC(J)**2
                      ENDDO
                  CALL QROD(NPARM,NPMAX,NPMAX,CM,PC,PU,UU,PX,PY)
                  ENDIF
              ENDDO
          RMSR= DSQRT(DSE/NDATA)
          IF(NPARM.LE.0) GO TO 60
c
c** Compute the inverse of  CM
          CM(1,1) = 1.D0 / CM(1,1)
          DO  I = 2,NPARM
              L = I - 1
              DO  J = 1,L
                  S = 0.D0
                  DO  K = J,L
                      S = S + CM(K,I) * CM(J,K)
                      ENDDO
                  CM(J,I) = -S / CM(I,I)
                  ENDDO
              CM(I,I) = 1.D0 / CM(I,I)
              ENDDO
c
c** Solve for parameter changes  PC(j)
          DO 26 I = 1,NPARM
              J = NPARM - I + 1
              PC(J) = 0.D0
              DO 24 K = J,NPARM
   24             PC(J) = PC(J) + CM(J,K) * PU(K)
   26         CONTINUE
c
c** Get (upper triangular) "dispersion Matrix" [variance-covarience
c  matrix  without the sigma^2 factor].
          DO 30 I = 1,NPARM
              DO 30 J = I,NPARM
                  UU = 0.D0
                  DO 28 K = J,NPARM
   28                 UU = UU + CM(I,K) * CM(J,K)
   30             CM(I,J) = UU
c** Generate core of Parameter Uncertainties  PU(j) and (symmetric)
c   correlation matrix  CM
          DO 36 J = 1,NPARM
              PU(J) = DSQRT(CM(J,J))
              DO 32 K= J,NPARM
   32             CM(J,K)= CM(J,K)/PU(J)
              DO 34 K= 1,J
                  CM(K,J)= CM(K,J)/PU(J)
   34             CM(J,K)= CM(K,J)
   36         CONTINUE
c** Option to print correlation matrix on first cycle ...
          IF((ITER.EQ.1).AND.(LPRINT.GE.6)) THEN
              WRITE(6,693) CM(1,1)
              DO i= 2,NPTOT
                  WRITE(6,694) i,(CM(i,k),k= 1,i)
                  ENDDO
              ENDIF
c
c** Generate standard error  DSE = sigma^2,  and prepare to calculate
c  Parameter Sensitivities PS
          IF(NDATA.GT.NPARM) THEN
              DSE= DSQRT(DSE/(NDATA-NPARM))
            ELSE
              DSE= 0.d0
            ENDIF
c** Use DSE to get final (95% confid. limit) parameter uncertainties PU
c** Calculate 'parameter sensitivities', changes in PV(j) which would
c  change predictions of input data by an RMS average of  DSE*0.1/NPARM
          UU= DSE*0.1d0/DBLE(NPARM)
          S= DSE*TFACT
          DO 40 J = 1,NPARM
              PU(J)= S* PU(J)
   40         PS(J)= UU*DSQRT(NDATA/PS(J))
c========End of core linear least-squares step==========================
c ... early exit if Rounding cycle finished ...
          IF(QUIT.GT.0) GO TO 60
c
c** Next test for convergence
          TSTPS= 0.D0
          TSTPU= 0.D0
          DO  J= 1,NPARM
              TSTPS= MAX(TSTPS,DABS(PC(J)/PS(J)))
              TSTPU= MAX(TSTPU,DABS(PC(J)/PU(J)))
              ENDDO
          IF(LPRINT.GE.4) WRITE(6,604) ITER,RMSR,TSTPS,TSTPU
c** Now ... update parameters (careful about rounding)
          DO  J= 1,NPTOT
              IF(IFXP(J).GT.0) THEN
c** If parameter held fixed (by rounding process), shift values of
c   change, sensitivity & uncertainty to correct label.
                  DO  I= NPTOT,J+1,-1
                      PC(I)= PC(I-1)
                      PS(I)= PS(I-1)
                      PU(I)= PU(I-1)
                      ENDDO
                  PC(J)= 0.d0
                  PS(J)= 0.d0
                  PU(J)= 0.d0
                ELSE
                  PV(J)= PV(J)+ PC(J)
                ENDIF
              ENDDO
          IF(LPRINT.GE.5) WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),
     1                                    J=1,NPTOT)
          IF(NITER.GT.1) THEN
c** Test for convergence:  for every parameter desire:
c  |parameter change| < |parameter sensitivity|   But STOP iterating
c  if  Max{|change/sens.|} increases AND Max{|change/unc.|} < 0.01
              IF(TSTPS.GT.1.d0) THEN
                  IF((RMSR.GT.RMSRB).AND.(ITER.GT.5)) THEN
                      IF((TSTPU.LT.1.d-2).OR.((TSTPU.LT.0.5d0).AND.
     1                                             (ITER.GT.10))) THEN
                          IF(LPRINT.GE.3) WRITE(6,606) ITER,TSTPU,RMSR
                          GO TO 54
                          ENDIF
                      ENDIF
                ELSE
                  IF(LPRINT.GE.3) WRITE(6,608) ITER,TSTPS,RMSR
                  GO TO 54
                ENDIF
              ENDIF
ccc       CALL FLUSH(6)
   50     CONTINUE
      WRITE(6,610) NPARM,NDATA,ITER,RMSR,TSTPS,TSTPU
c** End of iterative convergence loop for (in general) non-linear case.
c======================================================================
c
   54 IF(NPARM.EQ.NPTOT) THEN
          IF(LPRINT.NE.0) THEN
c** If desired, print unrounded parameters and fit properties
              WRITE(6,616) NDATA,NPARM,RMSR
              WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPARM)
              ENDIF
          IF(IROUND.EQ.0) RETURN
          IF(NGPRND.GT.1) THEN
c** For special case when wish to round off the last NGPRND parameters
c  in a single step ... (sometimes feasible for linear parameters)
              IF(NGPRND.GE.NPTOT) NGPRND= NPTOT-1
              CALL GPROUND(JROUND+1,NPTOT,NPARM,NPMAX,NGPRND,LPRINT,
     1                                                     IFXP,PV,PU)
              GO TO 6
              ENDIF
          ENDIF
c** Automated 'Sequential Rounding and Refitting' section:  round
c  current last parameter, fix it, and return (above) to repeat fit.
      IF(IROUND.LT.0) THEN
c ... if IROUND < 0, sequentially round off 'last' remaining parameter
          JFIX= NPARM
        ELSE
c ... if IROUND > 0, sequentially round off remaining parameter with
c                    largest relative uncertainty.
c ... First, select parameter with the largest relative uncertainty
          K= 0
          TSTPS= 0.d0
          DO  J= 1,NPTOT
              IF(IFXP(J).LE.0) THEN
                  K= K+1
                  TSTPSB= DABS(PU(J)/PV(J))
                  IF(TSTPSB.GT.TSTPS) THEN
                      JFIX= J
                      KFIX= K
                      TSTPS= TSTPSB
                      ENDIF
                  ENDIF
              ENDDO
c** Now redistribute correlation matrix elements for use by ROUND
          DO  J= 1,NPTOT
              IF(IFXP(J).GT.0) THEN 
                  DO  I= NPTOT,J+1,-1
                      CM(KFIX,I)= CM(KFIX,I-1)
                      ENDDO
                  ENDIF
              ENDDO
        ENDIF
      UU= PV(JFIX)
      CALL ROUND(JROUND,NPMAX,NPARM,NPTOT,JFIX,PV,PU,PS,CM)
      IFXP(JFIX)= 1
      IF(LPRINT.GE.2)
     1            WRITE(6,614) JFIX,UU,PU(JFIX),PS(JFIX),PV(JFIX),RMSR
      NPARM= NPARM-1
      IF(NPARM.EQ.0) THEN
c** After rounding complete, make one more pass with all parameters free
c  to get full correct corelation matrix, uncertainties & sensitivities
          NPARM= NPTOT
          QUIT= 1
          DO  J= 1,NPTOT
              IFXP(J)= 0
              ENDDO
c ... reinitialize for derivative-by-differences calculation
          RMSR= 0.d0
          ENDIF
      GO TO 6
c
c** If no parameters varied or sequential rounding completed - simply
c   calculate DSE from RMS residuals and return.
   60 DSE= 0.d0
      IF(NDATA.GT.NPTOT) THEN
          DSE= RMSR*DSQRT(DBLE(NDATA)/DBLE(NDATA-NPTOT))
        ELSE
          DSE= 0.d0
        ENDIF
      IF(NPTOT.GT.0) THEN
          IF(LPRINT.GT.0) THEN
c** Print final rounded parameters with original Uncert. & Sensitivities
              WRITE(6,616) NDATA,NPTOT,RMSR
              WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
              ENDIF
          ENDIF
      RETURN
c
  602 FORMAT(/' *** NLLSSRR problem:  [NPTOT=',i4,'] > min{NPINTMX=',
     1  i4,' NPMAX=',i4,', NDATA=',i6,'}')
  604 FORMAT(' After Cycle #',i2,':   RMSR=',1PD10.3,'    test(PS)=',
     1  1PD8.1,'    test(PU)=',D8.1)
  606 FORMAT(' Effective',i3,'-cycle Cgce:  MAX{|change/unc.|}=',1PD8.1,
     1  ' < 0.01   RMSR=',D10.3)
  608 FORMAT(' Full',i3,'- cycle convergence:  Max{|change/sens.|}=',
     1  1PD8.1,' < 1   RMSR=',D10.2)
  610 FORMAT(' !! CAUTION !! fit of',i4,' parameters to',I6,' data not c
     1onverged after',i3,' Cycles'/5x,'RMS(residuals)=',1PD10.3,
     2 '    test(PS) =',D9.2,'    test(PU) =',D9.2/1x,30('**'))
  612 FORMAT((4x,'PV(',i4,') =',1PD22.14,' (+/-',D8.1,')    PS=',d8.1,
     1  '   PC=',d8.1))
  614 FORMAT(' =',39('==')/' Round Off   PV(',i3,')=',1PD21.13,' (+/-',
     1 D9.2,')    PS=',d9.2/11x,'fix it as ',D21.13,'  & refit:  RMS(res
     2iduals)=',D10.3)
  616 FORMAT(/' Fit of',i6,' data to',i5,' parameters yields   RMS(resid
     1uals)=',G11.4)
  693 FORMAT(/14x,'Correlation Matrix'/'  1',f7.3,4x,9('--'))
  694 FORMAT(i3,20(f7.3))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
        SUBROUTINE QROD(N,NR,NC,A,R,F,B,GC,GS)
C** Performs ORTHOGONAL DECOMPOSITION OF THE LINEAR LEAST-SQUARES   
C            EQUATION J * X = F TO A * X = B(TRANSPOSE) * F WHERE   
C            J IS THE JACOBIAN IN WHICH THE FIRST N ROWS AND COLUMNS
C            ARE TRANSFORMED TO THE UPPER TRIANGULAR MATRIX A     
C            (J = B * A), X IS THE INDEPENDENT VARIABLE VECTOR, AND
C            F IS THE DEPENDENT VARIABLE VECTOR. THE TRANSFORMATION
C            IS APPLIED TO ONE ROW OF THE JACOBIAN MATRIX AT A TIME.
C  PARAMETERS :                                                   
C      N   -  (INTEGER) DIMENSION OF A TO BE TRANSFORMED.       
C      NR  -  (INTEGER) ROW DIMENSION OF A DECLARED IN CALLING PROGRAM.
C      NC  -  (INTEGER) Column DIMENSION OF F DECLARED IN CALLING PROGRAM.
C      A   -  (REAL*8 ARRAY OF DIMENSIONS .GE. N*N) UPPER TRIANGULAR
C             TRANSFORMATION MATRIX.                               
C      R   -  (REAL*8 LINEAR ARRAY OF DIMENSION .GE. N) ROW OF   
C             JACOBIAN TO BE ADDED.                             
C      F   -  (REAL*8 LINEAR ARRAY .GE. TO THE ROW DIMENSION OF THE
C             JACOBIAN) TRANSFORMED DEPENDENT VARIABLE MATRIX.   
C      B   -  (REAL*8) VALUE OF F THAT CORRESPONDS TO THE ADDED 
C             JACOBIAN ROW.                                     
C     GC   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS COSINE TRANSFORMATIONS.
C     GS   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS SINE TRANSFORMATIONS.
C--------------------------------------------------------------------
C  AUTHOR : MICHAEL DULICK, Department of Chemistry,
C           UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO N2L 3G1
C--------------------------------------------------------------------
      INTEGER  I,J,K,N,NC,NR
      REAL*8 A(NR,NC), R(N), F(NR), GC(N), GS(N), B, Z(2)
      DO 10 I = 1,N
          Z(1) = R(I)
          J = I - 1
          DO  K = 1,J
              Z(2) = GC(K) * A(K,I) + GS(K) * Z(1)
              Z(1) = GC(K) * Z(1) - GS(K) * A(K,I)
              A(K,I) = Z(2)
              ENDDO
          GC(I) = 1.D0
          GS(I) = 0.D0
          IF(Z(1) .EQ. 0.D0) GOTO 10
          IF(DABS(A(I,I)) .LT. DABS(Z(1))) THEN
              Z(2) = A(I,I) / Z(1)
              GS(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
              GC(I) = Z(2) * GS(I)
            ELSE
              Z(2) = Z(1) / A(I,I)
              GC(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
              GS(I) = Z(2) * GC(I)
            ENDIF
          A(I,I) = GC(I) * A(I,I) + GS(I) * Z(1)
          Z(2) = GC(I) * F(I) + GS(I) * B
          B = GC(I) * B - GS(I) * F(I)
          F(I) = Z(2)
   10     CONTINUE
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ROUND(IROUND,NPMAX,NPARM,NPTOT,IPAR,PV,PU,PS,CM)
c** Subroutine to round off parameter # IPAR with value PV(IPAR) at the
c  |IROUND|'th significant digit of:  [its uncertainty  PU(IPAR)] .
c** On return, the rounded value replaced the initial value  PV(IPAR).
c** Then ... use the correlation matrix CM and the uncertainties PU(I)
c  in the other (NPTOT-1) [or (NPARM-1) free] parameters to calculate
c  the optimum compensating changes PV(I) in their values.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1998  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER    IROUND,NPMAX,NPARM,NPTOT,IPAR,I,IRND,KRND
      REAL*8  PU(NPMAX),PS(NPMAX),PV(NPMAX),CM(NPMAX,NPMAX),CNST,
     1        CRND,XRND,FCT,Z0
      DATA Z0/0.d0/
      CNST= PV(IPAR)
      XRND= DLOG10(PU(IPAR))
c** If appropriate, base last rounding step on sensitivity (not uncert.)
      IF((NPARM.EQ.1).AND.(PS(IPAR).LT.PU(IPAR))) XRND= DLOG10(PS(IPAR))
c** First ... fiddle with log's to perform the rounding
      IRND= INT(XRND)
      IF(XRND.GT.0) IRND=IRND+1
      IRND= IRND- IROUND
      FCT= 10.D0**IRND
      CRND= PV(IPAR)/FCT
      XRND= Z0
c ... if rounding goes past REAL*8 precision, retain unrounded constant
      IF(DABS(CRND).GE.1.D+16) THEN
          WRITE(6,601) IROUND,IPAR
           RETURN
           ENDIF
      IF(DABS(CRND).GE.1.D+8) THEN
c ... to avoid problems from overflow of I*4 integers ...
          KRND= NINT(CRND/1.D+8)
          XRND= KRND*1.D+8
          CRND= CRND-XRND
          XRND= XRND*FCT
          END IF
      IRND= NINT(CRND)
      CNST= IRND*FCT+ XRND
c** Now ... combine rounding change in parameter # IPAR, together with
c  correlation matrix CM and parameter uncertainties PU to predict
c  changes in other parameters to optimally compensate for rounding off
c  of parameter-IPAR.  Method pointed out by Mary Thompson (Dept. of
c  Statistics, UW),
      IF(IPAR.GT.1) THEN
          XRND= (CNST-PV(IPAR))/PU(IPAR)
          DO  I= 1,NPTOT
              IF(I.NE.IPAR) THEN
                  PV(I)= PV(I)+ CM(IPAR,I)*PU(I)*XRND
                  ENDIF
              ENDDO
          ENDIF
      PV(IPAR)= CNST
      RETURN
  601 FORMAT(' =',39('==')/' Caution:',i3,'-digit rounding of parameter-
     1',i2,' would exceed (assumed) REAL*8'/' ********   precision overf
     2low at 1.D+16, so keep unrounded constant')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE GPROUND(IROUND,NPTOT,NPARM,NPMAX,NGPRND,LPRINT,
     1                                                     IFXP,PV,PU)
c** Subroutine to round off the last NGPRND of the NPTOT parameters
c  PV(i) at the |IROUND|'th significant digit of the smallest of their
c  uncertainties  min{U(i)}.  This procedure does NOT attempt to correct
c  the remaining parameters to compensate for these changes (as ROUND
c  does) and so is not appropriate for nonlinear parameters.
c** On return, the rounded values replaces the initial values of  PV(i).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2000  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,IROUND,NGPRND,NPMAX,NPTOT,NPARM,IPAR,IRND,KRND,LPRINT
      INTEGER  IFXP(NPTOT)
      REAL*8  PU(NPMAX),PV(NPMAX),CNST,CRND,XRND,FCT,XX,YY
c
c** Now, loop over & round off the last NGPRDN parameters,
      IF(LPRINT.GE.2) WRITE(6,602)  NGPRND,NPTOT
      NPARM= NPTOT+ 1
      DO  I= 1,NGPRND
          NPARM= NPARM- 1
c** First ... fiddle with log's to perform the rounding
          XRND= DLOG10(PU(NPARM))
          IRND= INT(XRND)
          IF(XRND.GT.0) IRND=IRND+1
          IRND= IRND- IROUND
          FCT= 10.D0**IRND
          CNST= PV(NPARM)
          YY= CNST
          CRND= PV(NPARM)/FCT
          XRND= 0.d0
c ... if rounding goes past REAL*8 precision, retain unrounded constant
          IF(DABS(CRND).GE.1.D+16) THEN
              WRITE(6,600) IROUND,IPAR
               RETURN
               ENDIF
          IF(DABS(CRND).GE.1.D+8) THEN
c ... to avoid problems from overflow of I*4 integers ...
              KRND= NINT(CRND/1.D+8)
              XRND= KRND*1.D+8
              CRND= CRND-XRND
              XRND= XRND*FCT
              END IF
          IRND= NINT(CRND)
          CNST= IRND*FCT+ XRND
          PV(NPARM) = CNST
          IFXP(NPARM)= 1
          IF(LPRINT.GE.2) WRITE(6,604) NPARM,YY,PV(NPARM)
  604 FORMAT(5x,'Round parameter #',i4,' from',G20.12,'  to',G20.12)
          ENDDO
          NPARM= NPARM- 1
      RETURN
  600 FORMAT(' =',39('==')/' Caution:',i3,'-digit rounding of parameter-
     1',i2,' would exceed (assumed) REAL*8'/' ********   precision overf
     2low at 1.D+16, so keep unrounded constant')
  602 FORMAT(' Rounding off the last',i5,' of',i5,' parameters:')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c     SUBROUTINE DYIDPJ(I,NDATA,NPTOT,UU,PV,PC,PS,RMSR)
c** Illustrative dummy version of DYIDPJ for the case of a fit to a
c  power series of order (NPTOT-1) in X(i). ***  For datum number-i,
c  calculate and return  PD(j)=[partial derivatives of datum-i] w.r.t.
c  each of the free polynomial coefficients varied in the fit
c  (for j=1 to NPTOT).  ** 
c* NOTE that  NDATA, PS and RMSR are useful for cases in which
c  derivatives-by-differences are generated (as for BCONT).
c=====================================================================
c** Use COMMON block(s) to bring in values of the independent variable
c  [here XX(i)] and any other parameters or variables needeed to
c  calculate YC and the partial derivatives.
c=====================================================================
c     INTEGER  I,J,NDATA,NPTOT,MXDATA
c     PARAMETER  (MXDATA= 501)
c     REAL*8  RMSR,YC,PV(NPTOT),PD(NPTOT),PS(NPTOT),POWER,XX(MXDATA)
c     COMMON /DATABLK/XX
c=====================================================================
c** NOTE BENE(!!) for non-linear fits, need to be sure that the
c  calculations of YC and PD(j) are based on the current UPDATED PV(j)
c  values.  If other (than PV) parameter labels are used internally
c  in the calculations, UPDATE them whenever (say)  I = 1 .
c=====================================================================
c     POWER= 1.D0
c     YC= PV(1)
c     PD(1)= POWER
c     DO 10 J= 2,NPTOT
c         POWER= POWER*XX(I)
c         YC= YC+ PV(J)*POWER
c         PD(J)= POWER
c  10     CONTINUE
c     RETURN
c     END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
