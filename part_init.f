c----------------------------------------------------------------------
      SUBROUTINE Energy_diag(vp,b0,b1,E,Evp,Euf,EB1,EB1x,EB1y,EB1z,
     x                       EE,EeP,etemp,nu,up,np)
c----------------------------------------------------------------------
      include 'incurv.h'

      real vp(Ni_max,3),
c     x     uf(nx,ny,nz,3),
c     x     nf(nx,ny,nz),
     x     b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     etemp(nx,ny,nz),
     x     nu(nx,ny,nz),
     x     up(nx,ny,nz,3),
     x     np(nx,ny,nz)

      real mO_q
      parameter(mO_q = mO/q)

      real Evp                  !kinetic energy of particles
      real Euf                  !kinetic energy of fluid flow
      real EB1,EB1x,EB1y,EB1z   !Magnetic field energy 
      real EE                   !Electric field energy 
      real EeP                  !Electron pressure energy
      real total_E              !total energy
      real aveEvp               !average particle energy
      real norm_E               !normalized energy
      real vol                  !volume of cell
      real denf                 !fluid density

      real recvbuf
      integer count
      count = 1

      Euf = 0.0
      EB1 = 0.0
      EB1x = 0.0
      EB1y = 0.0
      EB1z = 0.0
      EE = 0.0
      EeP = 0.0                
      do 10 i=1,nx-1
         j = 2
c         do 10 j=1,ny-1
            do 10 k=1,nz-1
               vol = dx*dy*dz_cell(k)*km_to_m**3
               EB1x = EB1x + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,1))**2 
               EB1y = EB1y + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,2))**2 
               EB1z = EB1z + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,3))**2 
c               EeP = EeP + kboltz*etemp(i,j,k)
               do 10 m=1,3
                  denf = np(i,j,k)/(km_to_m**3)
                  Euf = Euf + 0.5*mO*denf*vol*(up(i,j,k,m)*km_to_m)**2
                  EB1 = EB1 + 
     x              (vol/(2.0*mu0))*(mO_q*(b1(i,j,k,m)-b0(i,j,k,m)))**2
                  EE = EE + (epsilon*vol/2.0)*
     x                      (mO_q*E(i,j,k,m)*km_to_m)**2
 10               continue

c      input_EeP = input_EeP + EeP

c      write(*,*) 'Energy diag...',Ni_tot,m_arr(2000000)
 
      Evp = 0.0
      do 15 l=1,Ni_tot
         do 15 m=1,3
            Evp = Evp + 0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
 15   continue

c      write(*,*) 'Energy diag 2...',Ni_tot,m_arr(2000000)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(Evp,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_Evp = recvbuf

      call MPI_ALLREDUCE(input_E,recvbuf,count,
     x     MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      S_input_E = (recvbuf)


      total_E = S_Evp+EE+EB1
      aveEvp = S_Evp/S_input_E

c      write(*,*) 'Input energy (J).............',S_input_E
cc      write(*,*) 'Input EeP energy (J).........',input_EeP
c      write(*,*) 'Total vp energy (J)..........',S_Evp
c      write(*,*) 'Total up energy (J)..........',Euf
c      write(*,*) 'Total B energy (J)...........',EB1/S_input_E
c      write(*,*) 'Total E energy (J)...........',EE/S_input_E
cc      write(*,*) 'Total EeP energy (J).........',EeP
c      write(*,*) 'Total energy (J).............',total_E
cc      write(*,*) 'Total energy w/ eP (J).......',total_E+EeP
c      write(*,*) 'Energy thru boundaries.......',bndry_Eflux/S_input_E
c      write(*,*) 'Normalized particle energy...',aveEvp
      write(*,*) 'Normalized energy............',total_E/S_input_E,
     x   my_rank
      write(*,*) 'Normalized energy (bndry)....',
     x                (total_E)/(S_input_E+bndry_Eflux)
c      write(*,*) 'Normalized energy (no b1z)...',(S_Evp+Euf+EE+EB1x+
c     x                                            EB1y)/S_input_E
cc      write(*,*) 'Normalized energy (w/ eP)....',
cc     x                             (total_E+EeP)/(input_E + input_EeP)
c      write(*,*) ' '

      norm_E = total_E/S_input_E

c      if (prev_Etot .eq. 0.0) then prev_Etot = norm_E
c      do 20 i=1,nx 
c         do 20 j=1,ny
c            do 20 k=1,nz
c               nu(i,j,k) = nu(i,j,k) + 
c     x                 nu(i,j,k)*2.0*((norm_E - prev_Etot)/norm_E)
c 20            continue
      prev_Etot = norm_E

      return
      end
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup(np,vp,vp1,xp,input_p,up)
c----------------------------------------------------------------------
      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)

c      Ni_tot = 120000

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx)-qx(1))
         xp(l,2) = qy((ny/2)-8)+(1.0-pad_ranf())*(qy((ny/2)+8)-
     x                                        qy((ny/2)-8))
         xp(l,3) = qz((nz/2)-8)+(1.0-pad_ranf())*(qz((nz/2)+8)-
     x                                        qz((nz/2)-8))

c         i = nint(nx*pad_ranf())
c         j = nint((ny/2) + 16.*(0.5-pad_ranf()))
c         k = nint((nz/2) + 16.*(0.5-pad_ranf()))
cc         write(*,*) 'l...',l,i,j,k

c         xp(l,1) = qx(i)+dx*(0.5-pad_ranf())
c         xp(l,2) = qy(j)+dy*(0.5-pad_ranf())
c         xp(l,3) = qz(k)+dz_grid(k)*(0.5-pad_ranf())

         vp(l,1) = -vsw
         vp(l,2) = 0.0
         vp(l,3) = 0.0


         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*mproton*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + mproton*vp(l,m) / beta
 45      continue
 10      continue

        
 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_temp(np,vp,vp1,xp,input_p,up)
c----------------------------------------------------------------------
      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand

c      Ni_tot = 120000

c      do n = 0,procnum-1
c      if (my_rank .eq. n) then
      do 10 l = 1,Ni_tot
c         write(*,*) 'procnum, random number...',n,pad_ranf()

         phi = 2.0*pi*pad_ranf()
         flg = 0
         do 30 while (flg .eq. 0)
            theta = pi*pad_ranf()
            f = sin(theta)
            rnd = pad_ranf()
            if (f .ge. rnd) flg = 1
 30      continue


         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))


         flg = 0
         do 40 while (flg .eq. 0)
            v = (100*pad_ranf())
c            f = (vth**2/exp(1.0))*v**2*exp(-(v)**2 / vth**2)
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vp(l,1) = vsw + v*cos(phi)*sin(theta)
               vp(l,2) = v*sin(phi)*sin(theta)
               vp(l,3) = v*cos(theta)
            endif

c         vp(l,1) = -vsw
c         vp(l,2) = 0.0
c         vp(l,3) = 0.0

 40      continue


         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*mproton*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + mproton*vp(l,m) / beta
 45      continue
 10      continue
c      endif
c      enddo

 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl(np,vp,vp1,xp,input_p,up,np_t_flg,
     x                               np_b_flg)
c----------------------------------------------------------------------
      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      integer Ni_tot_1

      integer flg
      real nprat

      v1 = 1.0

      np_t_flg(:) = 0
      np_b_flg(:) = 0


      nprat = np_bottom/np_top
      Ni_tot_heavy = Ni_tot

c add heavies to bottom first

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)
            if (pad_ranf() .le. rnd) flg = 1
c            if (xp(l,3) .ge. qz(nz/2)) then
c               np_t_flg(l) = 1
c               m_arr(l) = m_top
c               mrat(l) = mproton/m_top
c            endif
c            if (xp(l,3) .lt. qz(nz/2)) then
            np_b_flg(l) = 1
            m_arr(l) = m_bottom
            mrat(l) = mproton/m_bottom
c            endif
            
 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 30 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 30      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = vth_bottom

c         vth = 0.5*(vth_top + vth_bottom) + 
c     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 44      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue


c add 10% initial Ni_tot for uniform proton background

         Ni_tot_1 = Ni_tot + 1
         Ni_tot = Ni_tot + nint(Ni_tot*0.1)


         do 49 l = Ni_tot_1,Ni_tot
            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))         
            np_b_flg(l) = 0
            m_arr(l) = mproton
            mrat(l) = 1.0

            ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
            ijkp(l,2) = nint(xp(l,2)/dy)
            
            k=1
            do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
               ijkp(l,3) = k    !grid
               k=k+1
 50         continue
            k=ijkp(l,3)
            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
               ijkp(l,3) = k+1
            endif
            
            vth = vth_top
            
            flg = 0
            do 52 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vx = v
               endif
 52         continue
            flg = 0
            do 54 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vy = v
               endif
 54         continue
            flg = 0
            do 56 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vz = v
               endif
 56         continue
            
            ii = ijkp(l,1)
            kk = ijkp(l,3)
            dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x           tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
            dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x           (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
            
            vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
            vp(l,2) = vy 
            vp(l,3) = vz        !+dvz 
            
            if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
            if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
            do  m=1,3
               vp1(l,m) = vp(l,m)
               input_E = input_E + 
     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
               input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
            enddo
            

 49      enddo


c add protons to top half

         Ni_tot_1 = Ni_tot+1
         Ni_tot = Ni_tot + Ni_tot_heavy*np_top/np_bottom

         do 60 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 70 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
             if (pad_ranf() .le. rnd) flg = 1
c             if (xp(l,3) .ge. qz(nz/2)) then
             np_t_flg(l) = 1
             m_arr(l) = m_top
             mrat(l) = mproton/m_top
c             endif
c             if (xp(l,3) .lt. qz(nz/2)) then
c                np_b_flg(l) = 1
c                m_arr(l) = m_bottom
c                mrat(l) = mproton/m_bottom
c             endif

 70      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 100 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 100     continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = vth_top

c         vth = 0.5*(vth_top + vth_bottom) + 
c     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 90 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 90      continue
         flg = 0
         do 92 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 92      continue
         flg = 0
         do 94 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 94      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 95 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 95      continue
 60   continue
c      endif
c      enddo




 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl_mix(np,vp,vp1,xp,input_p,up,
     x     np_t_flg,np_b_flg)
c----------------------------------------------------------------------
      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      integer Ni_tot_1

      integer flg
      real nprat

      v1 = 1.0

      np_t_flg(:) = 0
      np_b_flg(:) = 0

c      nprat = np_bottom/np_top
c      Ni_tot_heavy = Ni_tot

c add cold populations first

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
c         flg = 0
c         do 20 while (flg .eq. 0)
         xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
cc            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)  !for bottom
cc            rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0 !for top
c            if (pad_ranf() .le. rnd) flg = 1
cc            if (xp(l,3) .ge. qz(nz/2)) then
cc               np_t_flg(l) = 1
cc               m_arr(l) = m_top
cc               mrat(l) = mproton/m_top
cc            endif
cc            if (xp(l,3) .lt. qz(nz/2)) then
c            np_b_flg(l) = 1
         m_arr(l) = mproton
         mrat(l) = 1.0
c            endif
            
c 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 30 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 30      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = vth_bottom

c         vth = 0.5*(vth_top + vth_bottom) + 
c     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 44      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
c         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
c     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
c         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
c     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
c         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,1) = vsw + vx
         vp(l,2) = vy 
         vp(l,3) = vz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue


c add hot and fast population

         Ni_tot_1 = Ni_tot + 1

         Ni_tot = Ni_tot + 2*Ni_tot


         do 49 l = Ni_tot_1,Ni_tot

            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

c            flg = 0
c            do 22 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
c               rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0) !for bottom
c               rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0 !for top
c               if (pad_ranf() .le. rnd) flg = 1
c               np_t_flg(l) = 1
            m_arr(l) = mproton
            mrat(l) = 1.0
c 22         continue

            ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
            ijkp(l,2) = nint(xp(l,2)/dy)
            
            k=1
            do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
               ijkp(l,3) = k    !grid
               k=k+1
 50         continue
            k=ijkp(l,3)
            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
               ijkp(l,3) = k+1
            endif
            
            vth = vth_top
            
            flg = 0
            do 52 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vx = v
               endif
 52         continue
            flg = 0
            do 54 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vy = v
               endif
 54         continue
            flg = 0
            do 56 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vz = v
               endif
 56         continue
            
            ii = ijkp(l,1)
            kk = ijkp(l,3)
c            dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
c     x           tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
c            dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
c     x           (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
            
c            vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
            vp(l,1) = 1.0*vsw + vx !+dvx
            vp(l,2) = vy 
            vp(l,3) = vz        !+dvz 
            
            if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
            if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
            do  m=1,3
               vp1(l,m) = vp(l,m)
               input_E = input_E + 
     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
               input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
            enddo
            

 49      enddo


c add pickup distribution

         Ni_tot_1 = Ni_tot + 1

         Ni_tot = Ni_tot + 0.0*Ni_tot


         do 69 l = Ni_tot_1,Ni_tot

            xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
            xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

c            flg = 0
c            do 22 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz-1)-qz(1))
c               rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0) !for bottom
c               rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0 !for top
c               if (pad_ranf() .le. rnd) flg = 1
c               np_t_flg(l) = 1
            m_arr(l) = mproton
            mrat(l) = 1.0
c 22         continue

            ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
            ijkp(l,2) = nint(xp(l,2)/dy)
            
            k=1
            do 70 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
               ijkp(l,3) = k    !grid
               k=k+1
 70         continue
            k=ijkp(l,3)
            if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
               ijkp(l,3) = k+1
            endif
            
            vth = vth_top
            
            flg = 0
            do 72 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vx = v
               endif
 72         continue
            flg = 0
            do 74 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vy = v
               endif
 74         continue
            flg = 0
            do 76 while (flg .eq. 0)
               v = (2*vth_max*pad_ranf())-vth_max
               f = exp(-(v)**2 / vth**2)
               rnd = pad_ranf()
               if (f .ge. rnd) then 
                  flg = 1
                  vz = v
               endif
 76         continue
            
            ii = ijkp(l,1)
            kk = ijkp(l,3)

            theta = pad_ranf()*2*PI
            vp(l,1) = 1.0*vsw*cos(theta) + vx !+dvx
            vp(l,2) = vy 
            vp(l,3) = 1.0*vsw*sin(theta) + vz        !+dvz 
            
            if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
            if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
            do  m=1,3
               vp1(l,m) = vp(l,m)
               input_E = input_E + 
     x              0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
               input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
            enddo
            

 69      enddo



         
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end
c----------------------------------------------------------------------





c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl_eq(np,vp,vp1,xp,input_p,up,
     x                               np_t_flg,np_b_flg)
c----------------------------------------------------------------------
      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      integer Ni_tot_1
      real n1,n2,mr,nr,p0

      integer flg
      real nprat

      v1 = 1.0

      np_t_flg(:) = 0
      np_b_flg(:) = 0


      nprat = np_bottom/np_top

      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)
            if (pad_ranf() .le. rnd) flg = 1
c            if (xp(l,3) .ge. qz(nz/2)) then
c               np_t_flg(l) = 1
c               m_arr(l) = m_top
c               mrat(l) = mproton/m_top
c            endif
c            if (xp(l,3) .lt. qz(nz/2)) then
            np_b_flg(l) = 1
            m_arr(l) = m_bottom
            mrat(l) = mproton/m_bottom
c            endif
            
 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = vth_bottom

c         vth = 0.5*(vth_top + vth_bottom) + 
c     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 44      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue
      

         Ni_tot_1 = Ni_tot+1
         Ni_tot = Ni_tot + Ni_tot*np_top/np_bottom



         do 60 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 70 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
             if (pad_ranf() .le. rnd) flg = 1
c             if (xp(l,3) .ge. qz(nz/2)) then
             np_t_flg(l) = 1
             m_arr(l) = m_top
             mrat(l) = mproton/m_top
c             endif
c             if (xp(l,3) .lt. qz(nz/2)) then
c                np_b_flg(l) = 1
c                m_arr(l) = m_bottom
c                mrat(l) = mproton/m_bottom
c             endif

 70      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 100 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 100     continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif


c update np and then determine eq vth for protons (vth_top)


c      write(*,*) 'get interp weights...'
c      call get_interp_weights(xp)
c      write(*,*) 'update_np...'
c      call update_np(np)
c      write(*,*) 'update_up...'


c 60   continue

c      do 62 l = Ni_tot_1,Ni_tot


c         ii = ijkp(l,1)
c         jj = ijkp(l,2)
         kk = ijkp(l,3)

         vth = vth_top

         n1 = np_bottom*(1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
         n2 = np_top*(1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
      
         p0 = m_bottom*np_bottom*vth_bottom**2

         if ((n1 .gt. 0) .and. (n2 .gt. 0)) then 
         vth = sqrt((p0 - n1*m_bottom*vth_bottom*vth_bottom)/(n2*m_top))
c            write(*,*) 'vth_top...',l,vth,xp(l,3)
         endif

c         if ((n2 .gt. 0) .and. (n1 .gt. 0)) then 
c            vth = (n1*m_bottom*vth_bottom**2/(n2*m_top))**(0.5)
c            write(*,*) 'vth_top...',l,vth,xp(l,3)
c         endif



c         vth = 0.5*(vth_top + vth_bottom) + 
c     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 90 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 90      continue
         flg = 0
         do 92 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 92      continue
         flg = 0
         do 94 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 94      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 95 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 95      continue
 60   continue
c      endif
c      enddo



 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end
c----------------------------------------------------------------------





c----------------------------------------------------------------------
      SUBROUTINE part_setup_maxwl_p(np,vp,vp1,xp,input_p,up,np_t_flg,
     x                               np_b_flg)
c----------------------------------------------------------------------
      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      integer Ni_tot_1

      integer flg
      real nprat

      v1 = 1.0

      np_t_flg(:) = 0
      np_b_flg(:) = 0


      nprat = np_bottom/np_top


      do 10 l = 1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)
            if (pad_ranf() .le. rnd) flg = 1
c            if (xp(l,3) .ge. qz(nz/2)) then
c               np_t_flg(l) = 1
c               m_arr(l) = m_top
c               mrat(l) = mproton/m_top
c            endif
c            if (xp(l,3) .lt. qz(nz/2)) then
            np_b_flg(l) = 1
            m_arr(l) = m_bottom
            mrat(l) = mproton/m_bottom
c            endif
            
 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 44      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue
      

         Ni_tot_1 = Ni_tot+1
         Ni_tot = Ni_tot + Ni_tot*np_top/np_bottom

         do 60 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 70 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
             if (pad_ranf() .le. rnd) flg = 1
c             if (xp(l,3) .ge. qz(nz/2)) then
             np_t_flg(l) = 1
             m_arr(l) = m_top
             mrat(l) = mproton/m_top
c             endif
c             if (xp(l,3) .lt. qz(nz/2)) then
c                np_b_flg(l) = 1
c                m_arr(l) = m_bottom
c                mrat(l) = mproton/m_bottom
c             endif

 70      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 100 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 100     continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 90 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 90      continue
         flg = 0
         do 92 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 92      continue
         flg = 0
         do 94 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 94      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 95 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 95      continue
 60   continue
c      endif
c      enddo



 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE part_setup_maxwl_h(np,vp,vp1,xp,input_p,up,np_t_flg,
     x                               np_b_flg)
c----------------------------------------------------------------------
      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)
      integer Ni_tot_1

      integer flg
      real nprat

      v1 = 1.0

c      np_t_flg(:) = 0
c      np_b_flg(:) = 0


      nprat = np_bottom/np_top

      write(*,*) 'Ni_tot...',Ni_tot

      dNi = Ni_tot/2
      Ni_tot_1 = Ni_tot
      Ni_tot = Ni_tot + dNi/m_heavy

      write(*,*) 'Ni_tot...',Ni_tot,my_rank,m_heavy


      do 10 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))
         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
            rnd = ((1.0-tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0)
            if (pad_ranf() .le. rnd) flg = 1
c            if (xp(l,3) .ge. qz(nz/2)) then
c               np_t_flg(l) = 1
c               m_arr(l) = m_top
c               mrat(l) = mproton/m_top
c            endif
c            if (xp(l,3) .lt. qz(nz/2)) then
            np_b_flg(l) = 1
            m_arr(l) = m_heavy*mproton
            mrat(l) = 1.0/m_heavy
c            endif
            
 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 40      continue
         flg = 0
         do 42 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 42      continue
         flg = 0
         do 44 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 44      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 45      continue
 10      continue
      
         write(*,*) 'Ni_tot...',Ni_tot

         Ni_tot_1 = Ni_tot+1
c         Ni_tot = Ni_tot + Ni_tot*np_top/np_bottom
         Ni_tot = Ni_tot + dNi/m_heavy

         write(*,*) 'Ni_tot...',Ni_tot,my_rank


         do 60 l = Ni_tot_1,Ni_tot

         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 70 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0
             if (pad_ranf() .le. rnd) flg = 1
c             if (xp(l,3) .ge. qz(nz/2)) then
             np_t_flg(l) = 1
             m_arr(l) = m_heavy*mproton
             mrat(l) = 1.0/m_heavy
c             endif
c             if (xp(l,3) .lt. qz(nz/2)) then
c                np_b_flg(l) = 1
c                m_arr(l) = m_bottom
c                mrat(l) = mproton/m_bottom
c             endif

 70      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 100 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 100     continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 90 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vx = v
            endif
 90      continue
         flg = 0
         do 92 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vy = v
            endif
 92      continue
         flg = 0
         do 94 while (flg .eq. 0)
            v = (2*vth_max*pad_ranf())-vth_max
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .ge. rnd) then 
               flg = 1
               vz = v
            endif
 94      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 95 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*m_arr(l)*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + m_arr(l)*vp(l,m) / beta
 95      continue
 60   continue
c      endif
c      enddo



 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE sw_part_setup_maxwl_1(np,vp,vp1,xp,input_p,up,np_t_flg,
     x                               np_b_flg)
c----------------------------------------------------------------------
      include 'incurv.h'

      real np(nx,ny,nz)
      real vp(Ni_max,3)
      real vp1(Ni_max,3)
      real xp(Ni_max,3)
      real input_p(3)
      real up(nx,ny,nz,3)
      real phi,theta,rnd,f,v
      real rand
      real vx,vy,vz
      real dvx,dvz,v1
      integer np_t_flg(Ni_max)
      integer np_b_flg(Ni_max)

      integer flg
      real nprat

      v1 = 1.0

      np_t_flg(:) = 0
      np_b_flg(:) = 0


      nprat = np_bottom/np_top

c      Ni_tot = 120000

c      do n = 0,procnum-1
c      if (my_rank .eq. n) then
      Ni_tot = Ni_tot + Ni_tot*nint((1/nprat)-1)/2

      do 10 l = 1,Ni_tot
c         write(*,*) 'procnum, random number...',n,pad_ranf()

c         phi = 2.0*pi*pad_ranf()
c         flg = 0
c         do 30 while (flg .eq. 0)
c            theta = pi*pad_ranf()
c            f = sin(theta)
c            rnd = pad_ranf()
c            if (f .ge. rnd) flg = 1
c 30      continue


         xp(l,1) = qx(1)+(1.0-pad_ranf())*(qx(nx-1)-qx(1))
         xp(l,2) = qy(1)+(1.0-pad_ranf())*(qy(ny-1)-qy(1))

         
         flg = 0
         do 20 while (flg .eq. 0)
            xp(l,3) = qz(1)+(1.0-pad_ranf())*(qz(nz)-qz(1))
             rnd = (1-nprat)* 
     x             ((1.0+tanh((xp(l,3)-qz(nz/2))/(Lo)))/2.0) +
     x             nprat
c            write(*,*) 'rnd...',rnd
             if (pad_ranf() .le. rnd) flg = 1
             if (xp(l,3) .ge. qz(nz/2)) np_t_flg(l) = 1
             if (xp(l,3) .lt. qz(nz/2)) np_b_flg(l) = 1

 20      continue

         ijkp(l,1) = nint(xp(l,1)/dx) !particle grid location index
         ijkp(l,2) = nint(xp(l,2)/dy)
         
         k=1
         do 50 while(xp(l,3) .gt. qz(k)) !find k on non-uniform 
            ijkp(l,3) = k       !grid
            k=k+1
 50      continue
         k=ijkp(l,3)
         if (xp(l,3) .gt. (qz(k)+(dz_grid(k)/2))) then
            ijkp(l,3) = k+1
         endif

         vth = 0.5*(vth_top + vth_bottom) + 
     x     0.5*(vth_top - vth_bottom)*tanh((qz(ijkp(l,3))-qz(nz/2))/Lo)

         flg = 0
         do 40 while (flg .eq. 0)

            vy = (400*pad_ranf())-200
            vz = (400*pad_ranf())-200
            vx = (400*pad_ranf())-200
            
            v = sqrt(vx**2 + vy**2 + vz**2)
            f = exp(-(v)**2 / vth**2)
            rnd = pad_ranf()
            if (f .gt. rnd) then 
               flg = 1
            endif
 40      continue

         ii = ijkp(l,1)
         kk = ijkp(l,3)
         dvx = -2.0*v1*cos(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x    tanh((qz(kk)-qz(nz/2))/Lo)*(cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         dvz = v1*sin(PI*(qx(ii)-qx(nx/2))/qx(nx/2))*
     x       (cosh((qz(kk)-qz(nz/2))/Lo))**(-2)
         
         vp(l,1) = vsw*(tanh((qz(k)-qz(nz/2))/(Lo))) + vx !+dvx
         vp(l,2) = vy 
         vp(l,3) = vz !+dvz 

         if (xp(l,3) .gt. qz(nz/2)) mix_ind(l) = 1
         if (xp(l,3) .le. qz(nz/2)) mix_ind(l) = 0
            
         do 45 m=1,3
            vp1(l,m) = vp(l,m)
            input_E = input_E + 
     x           0.5*mproton*(vp(l,m)*km_to_m)**2 /beta
            input_p(m) = input_p(m) + mproton*vp(l,m) / beta
 45      continue
 10      continue
c      endif
c      enddo



 
      write(*,*) 'get interp weights...'
      call get_interp_weights(xp)
      write(*,*) 'update_np...'
      call update_np(np)
      write(*,*) 'update_up...'
      call update_up(vp,np,up)
      write(*,*) 'update_up complete...'

   

      return
      end
c----------------------------------------------------------------------











