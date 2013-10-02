c para.h
c contains simulation parameter list for a barium release

c simulation domain dimensions
      PARAMETER (nx = 3, ny = 3, nz = 401)

c magnetic field and mass info for determining time step

      real b0_init,mproton,q,nf_init
      PARAMETER (b0_init = 1700e-9)    !was 0.2
      PARAMETER (mproton = 16*1.67e-27)      
      PARAMETER (q = 1.6e-19)        !electron charge (Coulombs)
      PARAMETER (nf_init = 3000e15)   !was 0.01

      real lambda_i
      PARAMETER (lambda_i = (3e8/
     x            sqrt((nf_init/1e9)*q*q/(8.85e-12*mproton)))/1e3)

c grid parameters

      PARAMETER (dx = lambda_i/2, dy = lambda_i/2)   !units in km
      PARAMETER (delz = lambda_i/2)          !dz at release coordinates
                                       !this is not the final dz
                                       !d_lamba determines dz
c time stepping parameters
      PARAMETER (nt = 10000)        !number of time steps
      !PARAMETER (dtsub_init = 0.005) !subcycle time step 
      PARAMETER (ntsub = 10.0)   !number of subcycle time steps
      !PARAMETER (dt = dtsub_init*ntsub)     !main time step
      PARAMETER (dt = 0.2*mproton/(q*b0_init))     !main time step
      PARAMETER (dtsub_init = dt/ntsub) !subcycle time step 
      PARAMETER (nout = 2)      !number of steps to diagnosic output 

c output directory
      character(50) out_dir
      PARAMETER (out_dir='./tmp/')

c logical variable for restart
      logical restart
      PARAMETER (restart = .false.)
      PARAMETER (mrestart = 6000)      ! use -1 for no save
      PARAMETER (mbegin = 0)      !mstart

c neutral cloud expansion characteristics
      real vtop,vbottom,vth,vsw
c      PARAMETER(vo = 20.0)
c      PARAMETER(vth = 10.0)
      PARAMETER(vsat = 0.0)
      PARAMETER(vsw = 0.0*57.0)
      PARAMETER(vtop = vsw)
      PARAMETER(vbottom = -vsw)


c max number of ion particles to be produced.  This parameter
c is used to dimension the particle arrays.
      integer*4 Ni_max
      PARAMETER (Ni_max = 3000000)
 
c misc constants
      real mu0,epsilon,pi,rtod,mO,mBa,O_to_Ba,km_to_m,kboltz,melec
      real tempf0,m_pu
      PARAMETER (pi = 3.14159)
      PARAMETER (rtod = 180.0/pi)  !radians to degreesc
      PARAMETER (mu0 = pi*4.0e-7)  !magnetic permeability of free space
      PARAMETER (epsilon = 8.85e-12) !dielectric constant
c      PARAMETER (melec = 9.1e-31)  !mass of electron (kg)
      PARAMETER (mO = mproton)    !mass of H (kg)
c      PARAMETER (mO = 2.3e-25)    !mass of Ba (kg)
      PARAMETER (m_pu = 64.0)
      PARAMETER (mBa = m_pu*mO)    !mass of Ba (kg)
      PARAMETER (O_to_Ba = mO/mBa) !convert E and B for particle move
      PARAMETER (km_to_m = 1e3)    !km to meters conversion
      PARAMETER (kboltz = 1.38e-29)   !kg km^2 / s^2 / K
      PARAMETER (tempf0 = 50*11600.)     !K

      real nn_coef,np_top,np_bottom
      real b0_top,b0_bottom,Lo,vth_top,vth_bottom,vth_max
      real m_top, m_bottom,m_heavy,np_bottom_proton

      PARAMETER (m_heavy = 1.0)
      PARAMETER (np_top = nf_init)
      PARAMETER (np_bottom = nf_init/m_heavy)
      PARAMETER (f_proton_top = 0.5) !fraction relative to top
      PARAMETER (b0_top = 1.0*b0_init)
      PARAMETER (b0_bottom = b0_init)
      PARAMETER (vth_top = 30.00)
      PARAMETER (vth_bottom = 30.00)
      PARAMETER (vth_max = 3*30.0)
      PARAMETER (m_top = mproton)
      PARAMETER (m_bottom = mproton)
      PARAMETER (nn_coef = 1e5)
      PARAMETER (Lo = 4.0*dx) !gradient scale length of boundary


c electron ion collision frequency
      real nu_init, eta_init,lww1,lww2
      PARAMETER (nu_init = 0.0*q*b0_init/mproton)
      PARAMETER (eta_init = 0.0)
      PARAMETER (lww2 = 1.0)    !must be less than 1.0
      PARAMETER (lww1 = (1-lww2)/6.0)  !divide by six for nearest neighbor

c density scaling parameter, alpha, and ion particle array dims
       
      real alpha  
c      PARAMETER (alpha = 1.9263418e-20) !mH...determines particle scaling
      PARAMETER (alpha = (mu0/1e3)*q*(q/mproton)) !mH...determines particle scaling







