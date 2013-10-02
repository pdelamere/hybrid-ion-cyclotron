      PROGRAM MAIND
     
c----------------------------------------------------------------------
c maind.f
c Parallel version with no ion fluid, Nov 24, 2004
c----------------------------------------------------------------------

      include 'incurv.h'

c----------------------------------------------------------------------
c Listing of all declared variables
c
c Note that the format for specifying variables at different time
c levels is based upon the number of 1/2 time steps that the varible
c is behind the current time level.  For example, 
c uf2 is 1 time step behind uf, ufp2 is 1 time step ahead of uf,
c and b12 is 1 full time step behind b1 (not 12 1/2 steps behind b 
c since b does not exist...right). b1p2 is an exception...this is a
c temporary holder for b1 at m+1 in the predictor/corrector update
c of the magnetic field.
c----------------------------------------------------------------------
c      integer time, t1, t2    !keep track of run time
c      external time

      real b0(nx,ny,nz,3),            !ambient magnetic field
     x     b1(nx,ny,nz,3),    !1st order magnetic field
     x     b12(nx,ny,nz,3),   !b1 at previous time step
     x     b1p2(nx,ny,nz,3),  !temporary b1 at time level m+1
     x     bt(nx,ny,nz,3),    !total magnetic field..mc covarient
     x     btmf(nx,ny,nz,3),  !main cell contravarient bt field
     x     btc(nx,ny,nz,3),   !btmf at cell center for particle move
c     x     bdp(nx,ny,nz,3),   !dipole magnetic field
c     x     nf(nx,ny,nz),      !ambient fixed fluid density
c     x     nf1(nx,ny,nz),     !nf at n-1/2
c     x     nf3(nx,ny,nz),     !nf at n-3/2
c     x     nfp1(nx,ny,nz),    !nf at n+1/2
c     x     nn(nx,ny,nz),      !neutral cloud density
c     x     nnd(nx,ny,nz),     !neutral cloud density decrement
     x     np(nx,ny,nz),      !particle ion den at time level n, n+1/2
     x     np3(nx,ny,nz,3),
     x     vp(Ni_max,3),      !particle velocity at t level n+1/2
     x     vp1(Ni_max,3),     !particle velocity at t level n
     x     vplus(Ni_max,3),   !v+ used in velocity update
     x     vminus(Ni_max,3),  !v- used in velocity update
     x     up(nx,ny,nz,3),    !particle flow at time level n, n+1/2
     x     xp(Ni_max,3),      !coordinates of ion particles
c     x     uf(nx,ny,nz,3),    !fluid velocity
c     x     uf2(nx,ny,nz,3),   !fluid velcity at time level n-1
c     x     ufp1(nx,ny,nz,3),  !fluid velocity at time level n+1/2
c     x     ufp2(nx,ny,nz,3),  !fluid velocity at time level n+1
c     x     ui(nx,ny,nz,3),    !total ion flow velocity
     x     aj(nx,ny,nz,3),    !curlB/(alpha*n) 
     x     nu(nx,ny,nz),      !collision frequency
c     x     nuin(nx,ny,nz),    !ion-neutral collision frequency
     x     Ep(Ni_max,3),      !Ion particle electric field
c     x     Ef(nx,ny,nz,3),    !fluid electric field
     x     E(nx,ny,nz,3)     !electric field from electron mom eqn
c     x     uplus(nx,ny,nz,3), !u plus used in velocity update
c     x     uminus(nx,ny,nz,3),!u minus used in velocity update
c     x     pf(nx,ny,nz),      !fluid pressure at n
c     x     pf1(nx,ny,nz)      !fluid pressure at n-1/2

      real temp_p(nx,ny,nz)
      real mnp(nx,ny,nz) !mass density

      real Evp,       !total particle kinetic energy
     x     Euf,       !total fluid kinetic energy
     x     EB1,       !total magnetic field energy
     x     EB1x,      !total b1x energy
     x     EB1y,      !total b1y energy
     x     EB1z,      !total b1z energy
     x     EE,        !total electric field energy
     x     EeP        !total electron pressure energy

      real pup(3),      !total particle momentum
     x     puf(3),      !total fluid momentum
     x     peb(3),      !total momentum carried by E and B fields
     x     input_p(3)   !input momentum

      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      real np_t(nx,ny,nz)
      real np_b(nx,ny,nz)
      real up_t(nx,ny,nz,3)
      real up_b(nx,ny,nz,3)

      real chex_rate
      real bill_rate
      real satnp
c      real gradP(nx,ny,nz,3)
c      real etemp(nx,ny,nz)
c      real ugradu(nx,ny,nz,3)
c      real minnf,maxnf
c      real divu(nx,ny,nz)
      real mindt
      integer t1,t2,cnt_rt
      real time
      integer ierr
      
c      integer seed
      integer seedsize
      integer, dimension(:), allocatable :: seeder

      character(2) filenum(16) !max 16 processors

      filenum = (/'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ',
     x     '9 ','10','11','12','13','14','15','16'/)

c----------------------------------------------------------------------

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, procnum, ierr)


      call system_clock(t1,cnt_rt)
c      seed = float(t1)




c----------------------------------------------------------------------
c Initialize all variables
c----------------------------------------------------------------------
      write(*,*) 'initializing variables...'

      Ni_tot = 140000
      Ni_tot_sw = Ni_tot
      Ni_tot_sys = Ni_tot*procnum
      print *,'Ni_tot_sys, Ni_tot..',Ni_tot_sys,Ni_tot,Ni_tot_sw

   !----------------------------------------------------------------------
   !   check input parameters
   !----------------------------------------------------------------------

      write(*,*) 'alpha...',alpha
      write(*,*) 'c/wpi...',lambda_i,dx,dy,delz
      write(*,*) 'dt......',dt,dtsub_init

      write(*,*) ' '
      write(*,*) 'Bottom parameters...'
      write(*,*) ' '
      va =  b0_init/sqrt(mu0*m_bottom*np_bottom/1e9)/1e3
      write(*,*) 'Alfven velocity.......',va
      write(*,*) 'Thermal velocity......',vth_top
      write(*,*) 'Mach number...........',vbottom/(va + vth_bottom)

      write(*,*) 'Thermal gyroradius..',m_bottom*vth_bottom/(q*b0_init),
     x            m_bottom*vth_bottom/(q*b0_init)/dx
      cwpi = 3e8/sqrt((np_bottom/1e9)*q*q/(epsilon*m_bottom))
      write(*,*) 'Ion inertial length...',cwpi/1e3,cwpi/1e3/dx

      write(*,*) 'Particles per cell....',Ni_tot_sys/(nx*nz)

      write(*,*) ' '
      write(*,*) 'Top parameters...'
      write(*,*) ' '

      va =  b0_init/sqrt(mu0*m_top*np_top/1e9)/1e3
      write(*,*) 'Alfven velocity.......',va
      write(*,*) 'Thermal velocity......',vth_top
      write(*,*) 'Mach number...........',vtop/(va + vth_top)

      write(*,*) 'Thermal gyroradius....',m_top*vth_top/(q*b0_init),
     x            m_top*vth_top/(q*b0_init)/dx
      cwpi = 3e8/sqrt((np_top/1e9)*q*q/(epsilon*m_top))
      write(*,*) 'Ion inertial length...',cwpi/1e3,cwpi/1e3/dx

      write(*,*) 'Particles per cell....',Ni_tot_sys/(nx*nz)
      write(*,*) ' '
      
c      stop

   !----------------------------------------------------------------------



c      Ni_tot = 6
      mstart = 0
      ndiag = 0
      prev_Etot = 1.0
      nuei = 0.0

c initialize seed for each processor

      seed = 123456789 +my_rank*100
      call random_initialize(seed)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)            
      write(*,*) 'rand init ...',ranf(),my_rank


c      call random_seed
c      call random_seed(size = seedsize)
c      allocate(seeder(seedsize))
c      do n = 0,procnum-1 
c         if (my_rank .eq. n) then 
c            call random_seed(get=seeder)
c            call random_seed(put=seeder)
c         endif
c      enddo



      if (.not.(restart)) then
         do 66 i=1,nx
            do 66 j=1,ny
               do 66 k=1,nz
c                  pf(i,j,k) = nf_init*0.05*kboltz*tempf0
c                  pf1(i,j,k) = nf_init*0.05*kboltz*tempf0
c                  nf(i,j,k) = nf_init*0.0
c                  nf1(i,j,k) = nf_init*0.05  
c                  nf3(i,j,k) = nf_init*0.05 
c                  nfp1(i,j,k) = nf_init*0.05  
                  input_E = 0.0
                  input_p = 0.0
                  input_chex = 0.0
                  input_bill = 0.0
 66               continue
         endif

c      do 68 i = 1,nx
c         do 68 j = 1,ny
c            do 68 k = 1,nz
c               uf(i,j,k,1) = -vsw
c               uf2(i,j,k,1) = -vsw
c               ufp1(i,j,k,1) = -vsw
c               ufp2(i,j,k,1) = -vsw
c 68            continue

c      Ni_tot = 4000000

c      m_arr(1:Ni_tot) = mproton
c      m_arr(Ni_tot+1:) = m_pu*mproton  !mass N_2+ = 28.0
c      mrat(1:Ni_tot) = 1.0
c      mrat(Ni_tot+1:) = 1.0/m_pu  !mass N_2+ = 28.0

c add 10% heavy ions

c      m_arr(1:Ni_tot*0.1) = mproton*16.0
c      mrat(1:Ni_tot*0.1) = 1.0/16.0

c      beta = (Ni_tot/(nx*dx*ny*dy*nz*delz))/nf_init

c      write(*,*) 'beta...',beta
  
      call grd7()
      call grd6_setup(b0,bt,b12,b1,b1p2,nu)

      call get_beta()

      input_E = 0.0
c      do i = 1,nx
c         do j = 1,ny
c            do k = 1,nz
c               input_E = input_E + 
c     x          0.5*dx*dy*dz_grid(k)*nf_init*0.01*mO*(vsw*km_to_m)**2
c            enddo
c         enddo
c      enddo

      write(*,*) 'SW particle setup maxwl 1...'

c      call sw_part_setup_temp(np,vp,vp1,xp,input_p,up)
c      call sw_part_setup_maxwl(np,vp,vp1,xp,input_p,up,np_t_flg,
c     x                         np_b_flg)

      call sw_part_setup_maxwl_mix(np,vp,vp1,xp,input_p,up,np_t_flg,
     x                         np_b_flg)
c      call part_setup_maxwl_p(np,vp,vp1,xp,input_p,up,np_t_flg,
c     x                         np_b_flg)
c      call part_setup_maxwl_h(np,vp,vp1,xp,input_p,up,np_t_flg,
c     x                         np_b_flg)
c      write(*,*) 'SW particle setup complete...',Ni_tot,
c     x                 mrat(Ni_tot:Ni_tot+1)


      Ni_tot_sys = Ni_tot*procnum
      write(*,*) 'Ni_tot_sys...top and bottom...',Ni_tot_sys,Ni_tot
      write(*,*) 'Particles per cell....',Ni_tot_sys/(nx*nz)
c      stop

      call f_update_tlev(b1,b12,b1p2,bt,b0)

c----------------------------------------------------------------------




c----------------------------------------------------------------------
c check for restart flag
c----------------------------------------------------------------------
      write(*,*) 'restart status....',restart
      if (restart) then 
         write(*,*) 'opening restart.vars......'
         open(210,file='restart.vars',status='unknown',
     x            form='unformatted')
         write(*,*) 'reading restart.vars......'
         read(210) b1,b12,b1p2,bt,btmf,nn,np,nf,vp,vp1,vplus,vminus,
     x            up,xp,uf,uf2,ufp2,aj,Ep,Ef,E,uplus,uminus,Evp,Euf,
     x            EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot,
     x            ijkp,mstart,input_p,input_EeP,prev_Etot,nf1,nf3,nfp1,
     x            input_chex,input_bill,pf,pf1,mrat,m_arr
         write(*,*) 'restarting hybrid.....'

         if (my_rank .gt. 0) then 
          open(211,file='restart.part'//trim(filenum(my_rank)),
     x            status='unknown',form='unformatted')
          read(211) vp,vp1,vplus,vminus,xp,Ep,input_E,Ni_tot,
     x              ijkp,input_p,mrat,m_arr
         endif
      endif
      
      close(211)
c----------------------------------------------------------------------


c----------------------------------------------------------------------
c write para.h file

      if (my_rank .eq. 0) then

         open(109, file=trim(out_dir)//'para.dat',
     x        status='unknown',form='unformatted')
         
         write(109) nx,ny,nz,dx,dy,delz
         write(109) nt,dtsub_init,ntsub,dt,nout
         write(109) out_dir
c         write(109) model_choice
c         write(109) nf_init,b0_init
c         write(109) nu_init,lww2,lww1
c         write(109) Mdot,Mdot_part
         write(109) vtop,vbottom
         write(109) Ni_max
         write(109) mproton,m_pu,m_heavy
         write(109) np_top,np_bottom
         write(109) b0_top,b0_bottom
         write(109) vth_top,vth_bottom
c         write(109) RIo
         write(109) alpha,beta
c         write(109) comm_sz
         close(109)

      endif
 

c----------------------------------------------------------------------


c----------------------------------------------------------------------
c Initialize diagnostic output files
c----------------------------------------------------------------------
      if (my_rank .eq. 0) then 
         open(110,file=trim(out_dir)//'c.np.dat',status='unknown',
     x        form='unformatted')
         
         open(115,file=trim(out_dir)//'c.np_b.dat',status='unknown',
     x        form='unformatted')
         
         open(120,file=trim(out_dir)//'c.mixed.dat',status='unknown',
     x         form='unformatted')
         
         open(130,file=trim(out_dir)//'c.b1.dat',status='unknown',
     x        form='unformatted')
         
         open(140,file=trim(out_dir)//'c.aj.dat',status='unknown',
     x        form='unformatted')
         
         open(150,file=trim(out_dir)//'c.E.dat',status='unknown',
     x        form='unformatted')
         
         open(160,file=trim(out_dir)//'c.energy.dat',status='unknown',
     x        form='unformatted')
         
c     open(170,file='c.chex.dat',status='unknown',
c     x         form='unformatted')
         
c     open(172,file='c.bill.dat',status='unknown',
c     x         form='unformatted')
         
c     open(175,file='c.satnp.dat',status='unknown',
c     x         form='unformatted')
         
         open(180,file=trim(out_dir)//'c.up.dat',status='unknown',
     x        form='unformatted')
         
         open(190,file=trim(out_dir)//'c.momentum.dat',
     x        status='unknown',
     x        form='unformatted')
         
         open(192,file=trim(out_dir)//'c.p_conserve.dat',
     x        status='unknown',               
     x        form='unformatted')                 
         
         open(300,file=trim(out_dir)//'c.temp_p.dat',status='unknown',
     x         form='unformatted')
         
         open(305,file=trim(out_dir)//'c.xp_0.dat',status='unknown',
     x        form='unformatted')
         
         open(310,file=trim(out_dir)//'c.vp_0.dat',status='unknown',
     x        form='unformatted')

         open(315,file=trim(out_dir)//'c.mrat_0.dat',status='unknown',
     x        form='unformatted')

         open(320,file=trim(out_dir)//'c.np_wake.dat',
     x        status='unknown',form='unformatted')
         
         open(330,file=trim(out_dir)//'c.up_t.dat',status='unknown',
     x         form='unformatted')

         open(340,file=trim(out_dir)//'c.up_b.dat',status='unknown',
     x         form='unformatted')
         
c     open(330,file=trim(out_dir)//'c.ufp2.dat',status='unknown',
c     x         form='unformatted')
         
      open(342,file=trim(out_dir)//'c.test_part.dat',status='unknown',
     x         form='unformatted')
         
      open(350,file=trim(out_dir)//'c.mnp.dat',status='unknown',
     x         form='unformatted')
      endif
      
      
      if (my_rank .gt. 0) then
         open(305,file=trim(out_dir)//'c.xp_'//trim(filenum(my_rank))
     x        //'.dat',
     x        status='unknown',form='unformatted')
         
         open(310,file=trim(out_dir)//'c.vp_'//trim(filenum(my_rank))
     x        //'.dat',
     x        status='unknown',form='unformatted')

         open(315,file=trim(out_dir)//'c.mrat_'//trim(filenum(my_rank))
     x        //'.dat',
     x        status='unknown',form='unformatted')
c         write(*,*) 'filename...',
c     x          my_rank,'c.xp_'//trim(filenum(my_rank))//'.dat'

      endif
c----------------------------------------------------------------------




c======================================================================
c  MAIN LOOP!
c======================================================================

      vp(1,:) = 0.0
c      xp(1,1) = qx(nx/2)+dx/2
c      xp(l,2) = qy(ny/2)+dx/2
c      xp(l,3) = qz(nz/2)+dx/2

      do 1 m = mstart+1, nt

c        write(*,*) ' '
         write(*,*) 'time...', m, m*dt,my_rank

         write(*,*) 'm1...',m_arr(1),mrat(1)

         !Calculate neutral density

         !Ionize cloud and calculate ion density
c         if (Ni_tot .lt. Ni_max) then
cc         call Ionize_pluto(np,vp,vp1,xp,m,input_p,up)
c         call Ionize_pluto_mp(np,vp,vp1,xp,m,input_p,up)
         call Ionize_sw_mp(np,vp,vp1,xp,m,input_p,up)
c         endif

         write(*,*) 'Ni_tot...',Ni_tot,dNi,my_rank !/beta
c         write(*,*) ' '

         call get_interp_weights(xp)
         call update_np(np)             !np at n+1/2
         call update_up(vp,np,up)       !up at n+1/2


         !energy diagnostics
         call get_bndry_Eflux(b1,E)
         call Energy_diag(vp,b0,b1,E,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,
     x                    EeP,etemp,nu,up,np)

         call curlB(bt,np,aj)
         call cov_to_contra(bt,btmf)
         call face_to_center(btmf,btc)       !interp bt to cell center
         call extrapol_up(up,vp,vp1,np)
         call get_Ep(Ep,aj,np,up,btc,nu)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call improve_up(vp1,vplus,vminus,up,np)

         call get_Ep(Ep,aj,np,up,btc,nu)
         call get_vplus_vminus(Ep,btc,vp,vplus,vminus)
         call get_vp_final(Ep,vp,vp1,vplus)

         call move_ion_half(xp,vp,vp1,input_p) !1/2 step ion move to n+1/2
c         call check_min_den_boundary(np,xp,vp,up)

         call get_interp_weights(xp)
         call update_np(np)             !np at n+1/2
         call update_up(vp,np,up)       !up at n+1/2

c         write(*,*) 'np min k=1...',minval(np(:,:,1))/np_bottom

c         write(*,*) 'np max...',maxval(np(:,:,:))
         
c**********************************************************************
c SUBCYCLING LOOP!
c**********************************************************************

         dtsub = dtsub_init
         ntf = ntsub
         mindt = dtsub_init/50.
c check time step
c         write(*,*) 'checking time step...',ntf
c         do i = 1,nx
c            do j = 1,ny
c               do k = 1,nz
c                  ak = 2./dx
c                  btot = sqrt(bt(i,j,k,1)**2 + bt(i,j,k,2)**2 + 
c     x                 bt(i,j,k,3)**2)
c                  a1 = ak**2*Btot/(alpha*(nf(i,j,k)+np(i,j,k)))
c                  a2 = (ak*Btot)**2/(alpha*(nf(i,j,k)+np(i,j,k)))
c                  womega = 0.5*(a1 + sqrt(a1**2 + 4*a2))
c                  phi = womega/ak
c                  deltat = dx/phi
c                  if(deltat .le. 2.0*dtsub) then 
c           write(*,*) 'time stepping error...',deltat,nf(i,j,k)
c                     if (mindt .gt. deltat) then
c                        deltat = mindt
c                        write(*,*) 'mindt...',mindt
c                     endif
c                     do while (2.0*dtsub .gt. deltat)
c                        dtsub = dtsub/2.0
c                        ntf = ntf*2.0
c      write(*,*) 'Changing subcycle time step...',dtsub,deltat,ntf
c                     enddo
c                  endif
c               enddo
c            enddo
c         enddo
         

      do 2 n = 1, ntf

c         write(*,*) 'subcycle step...',n,ntf

         !convert main cell covarient bt to main cell contravarient
         call cov_to_contra(bt,btmf) 
         call curlB(bt,np,aj)     

         !update fluid velocity, uf 

c only need predict_uf when calculating ugradu

cc         call trans_nf_Lax(nf,nf1,nfp1,uf) 
c         call trans_nf_LaxWend1(nf,nf1,nfp1,uf)
c         call trans_pf_LaxWend1(pf,pf1,uf)

c         call get_nuin(nuin,nn,uf)
c         call predict_uf(Ef,b0,b1,b12,uf,uf2,ufp2,nu,np,nf1,uplus, 
c     x                   uminus,ugradu,up,gradP,nuin,bdp,pf1)

c         call predict_nf(nf,nf1,nf3,nfp1,uf,divu,b1)  

c         call get_nuin(nuin,nn,uf)
c         call correct_uf(Ef,btmf,uf,uf2,ufp2,nu,np,nf,uplus,uminus, 
c     x                   ugradu,aj,up,ufp1,gradP,nuin,pf)

c         call trans_nf_LaxWend2(nf,nf1,nfp1,ufp1)
c         call trans_pf_LaxWend2(pf,pf1,ufp1)

         !update magnetic field, b1
c         call predict_B(b1,b12,b1p2,bt,btmf,E,aj,up,uf,uf2,np,nf,nu,
c     x                  gradP) 

         call predict_B(b0,b1,b12,b1p2,bt,btmf,E,aj,up,np,nu) 

c         call correct_nf(nf,nf1,ufp1)

c         call correct_B(b0,b1,b1p2,E,aj,up,uf,np,nfp1,nu,gradP,bdp)
         call correct_B(b0,b1,b1p2,E,aj,up,np,nu)

c         call f_update_tlev(uf,uf2,b1,b12,b1p2,bt,b0,bdp)
         call f_update_tlev(b1,b12,b1p2,bt,b0)


c         call Momentum_diag(up,uf,np,nf,E,b1,pup,puf,peb,input_p)
c         call check_momentum(uf,nf,bt,ugradu)
c         write(192) m
c         write(192) n
c         write(192) surf_tot, graduu_tot, ugradu_tot



 2     continue
c**********************************************************************

         call move_ion_half(xp,vp,vp1,input_p)  !final ion move to n+1
c         call check_min_den_boundary(np,xp,vp,up)

c         call check_min_den(np,xp,vp,vp1,up,bt)

c         call res_chex(xp,vp)

c         write(*,*) 'Momentum conservation...'
c         write(*,*) '  Particles.............',pup(1),pup(2),pup(3)
c         write(*,*) '  Fluid.................',puf(1),puf(2),puf(3)
c         write(*,*) '  ExB...................',peb(1),peb(2),peb(3)
c         write(*,*) '  Normalized............',
c     x                        (pup(1)+puf(1)+peb(1))/input_p(1),
c     x                        (pup(2)+puf(2)+peb(2))/input_p(2),
c     x                        (pup(3)+puf(3)+peb(3))/input_p(3)

         call get_np3(np,np3)
         call update_mixed

         if (my_rank .eq.0) then
            write(342) xp(1:100,:)
         endif
            
c         write(175) b1(nx/2,ny/2,nz/2,:)

c----------------------------------------------------------------------
c diagnostic output
c----------------------------------------------------------------------

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         if (my_rank .eq. 0) then
            write(160) m
            write(160) input_E,input_EeP,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,
     x           EeP,input_chex,input_bill
            write(190) m
            write(190) pup, puf, peb, input_p
            write(320) np(ri-20,rj,rk),np(ri-40,rj,rk),
     x                 np(ri-40,rj,rk+50),np(ri+5,rj,rk)
         endif

         ainput_E = input_E
         abndry_Eflux = bndry_Eflux

         ndiag = ndiag + 1
         if (ndiag .eq. nout) then

c            call separate_np(np_t,np_t_flg)
c            call separate_np(np_b,np_b_flg)
c            call separate_up(vp,np,np_t_flg,up_t)
c            call separate_up(vp,np,np_b_flg,up_b)
            call get_temperature(xp,vp,np,temp_p)
            call update_rho(mnp)
            if (my_rank .eq. 0) then

c     if (m .ge. 275) then
               write(110) m
               write(110) np
c     write(110) nn
               write(115) m
               write(115) np_b 
               write(120) m
               write(120) mixed
               write(130) m
c     write(130) b1
               write(130) bt
               write(140) m
               write(140) aj*alpha*np3
               write(150) m
               write(150) E
c     write(150) Ef
               write(180) m
               write(180) up
               write(300) m
               write(300) temp_p/1.6e-19
c               write(305) m
c               write(305) xp
c               write(310) m
c               write(310) vp
c               write(315) m
c               write(315) mrat
c               write(330) m
               write(330) up_t
               write(340) m
               write(340) up_b
c     write(320) m
c     write(320) uf2
c     write(330) m
c     write(330) ufp2
c     write(340) m
c     write(340) etar
               write(350) m
c     write(350) (pf/nf)/1.6e-25  !saves in units of eV
               write(350) mnp
               ndiag = 0
            endif
 
            if (my_rank .gt. 0) then
c               write(305) m
c               write(305) xp
c               write(310) m
c               write(310) vp
c               write(315) m
c               write(315) mrat
               ndiag = 0
            endif
         endif

         input_E = ainput_E
         bndry_Eflux = abndry_Eflux
         

c----------------------------------------------------------------------


c----------------------------------------------------------------------
c Write restart file
c----------------------------------------------------------------------
         if (my_rank .eq. 0) then
            if (m .eq. mrestart) then
               open(220,file='restart.vars.new',status='unknown',
     x              form='unformatted')
           write(220) b1,b12,b1p2,bt,btmf,nn,np,nf,vp,vp1,vplus,vminus,
     x          up,xp,uf,uf2,ufp2,aj,Ep,Ef,E,uplus,uminus,Evp,Euf,
     x          EB1,EB1x,EB1y,EB1z,EE,EeP,input_E,Ni_tot,
     x          ijkp,mrestart,input_p,input_EeP,prev_Etot,nf1,nf3,nfp1,
     x          input_chex,input_bill,pf,pf1,mrat,m_arr  
            endif
         endif
         
         if (my_rank .gt. 0) then
            if (m .eq. mrestart) then
          open(221,file='restart.part'//trim(filenum(my_rank))//'.new',
     x              status='unknown',form='unformatted')
               write(221) vp,vp1,vplus,vminus,xp,Ep,input_E,Ni_tot,
     x              ijkp,input_p,mrat,m_arr
            endif
         endif
c----------------------------------------------------------------------

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)


 1     continue
c======================================================================

c       if(my_rank .eq. 0) then

          close(110)
          close(115)
          close(120)
          close(130)
          close(140)
          close(150)
          close(160)
          close(170)
          close(172)
          close(175)
          close(180)
          close(190)
          close(192)
          close(210)
c          close(211)
          close(220)
          close(221)
          close(300)
          close(305)
          close(310)
          close(315)
          close(320)
          close(330)
          close(340)
          close(342)
          close(350)

c       endif

       call system_clock(t2,cnt_rt)
       time = (real(t2) - real(t1))/real(cnt_rt)
       write(*,*) 
       write(*,*) 'Elapsed time....',time,' sec'
       write(*,*)

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       call MPI_FINALIZE(ierr)

       stop 
       end



















