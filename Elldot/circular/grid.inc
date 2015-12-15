      INTEGER NR,N0
      DOUBLE PRECISION DR,E0,M0,A0,OMEGA0,B,PI
C**---------------------------------------------------------------------
C**   these are the parameters for GaAs:
C**
C**---------------------------------------------------------------------
      PARAMETER(E0 = 12.4D0, M0 = 0.067D0,  A0 = 0.529177*E0/M0,
     &          OMEGA0 = 5.00D0, B=0.00D0, PI = 3.141592653589793)
C      PARAMETER(E0 = 13.1D0, M0 = 0.0665D0,  A0 = 0.529177*E0/M0,
C     &          OMEGA0 = 3.97D0, B=0.0D0, PI = 3.141592653589793)
C      PARAMETER(E0 = 15.15D0, M0 = 0.08D0,  A0 = 0.529177*E0/M0,
C     &          OMEGA0 = 45.0D0, B=12.0D0, PI = 3.141592653589793)
C**---------------------------------------------------------------------
C**
C**   E0 = dielectric constant
C**
C**   M0 = effective mass
C**
C**   A0 = effective Bohr radius (in Angstrom)
C**
C**   OMEGA0 = curvature of the parabolic potential (in meV)
C**
C**   B = external magnetic field in Tesla
C**
C**
C**
C**   Here are the parameters for the grid:
C**
C**   NR = number of radial grid points
C**
C**
C**---------------------------------------------------------------------
      PARAMETER(NR=400, DR=0.017D0, N0=300)
C      PARAMETER(NR=400, DR=0.018D0, N0=310)
C      PARAMETER(NR=400, DR=0.02D0, N0=245)
C      PARAMETER(NR=600, DR=0.02D0, N0=450)
C      PARAMETER(NR=850, DR=0.01D0, N0=540)
C      PARAMETER(NR=900, DR=0.01D0, N0=600)
