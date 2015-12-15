      INTEGER REST,XC,NU,ND,NLU,NLD
      DOUBLE PRECISION RESID_0,REF,CRITERION,MIX 
C***********************************************************************
      PARAMETER(REST=0, XC=0,
     &          RESID_0=1.D-06, REF=5.0, CRITERION=1.D-9, MIX=0.50,
     &          NU = 3, ND = 3,
     &          NLU = 3, NLD = 3)
C***********************************************************************
C**
C**  XC=0: only xc scalar potential.
C**        This case applies for zero total angular momentum, i.e.
C**        for closed shell, where the paramegnetic current is zero.
C**
C**  XC=1: include vector xc potential.
C**        Choose this option for all cases of open shells, i.e. 
C**        when there is a finite current even in the absence
C**        of an external magnetic field.
C**
C**  NU: total number of spin-up electrons 
C**  ND: total number of spin-down electrons 
C**
C**  NLU: number of different L-values occurring in spin-up
C**  NLD: number of different L-values occurring in spin-down
C**
C**  example: if the configuration is 1s 1p (up) and 1s 2s (down), then
C**  NU=2, ND=2, NLU=2, NLD=1
C**
C**  NU,ND,NLU,NLD refer to the orbitals considered in the calculation.
C**  In general, we want to determine more orbitals than actually occupied.
C**
C**

