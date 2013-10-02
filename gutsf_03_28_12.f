c----------------------------------------------------------------------
      SUBROUTINE f_update_tlev(b1,b12,b1p2,bt,b0)
c loops run 1 to n since values are only being copied
c----------------------------------------------------------------------
      include 'incurv.h'
      
      real b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     b0(nx,ny,nz,3)
c     x     bdp(nx,ny,nz,3)

      do 10 i=1,nx
         do 10 j=1,ny
            do 10 k=1,nz
               bt(i,j,k,1) = b1p2(i,j,k,1)! + b0(i,j,k,1)
               bt(i,j,k,2) = b1p2(i,j,k,2)! + b0(i,j,k,2)
               bt(i,j,k,3) = b1p2(i,j,k,3)! + b0(i,j,k,3) 
               do 10 m=1,3
c                  uf2(i,j,k,m) = uf(i,j,k,m)
                  b12(i,j,k,m) = b1(i,j,k,m)
                  b1(i,j,k,m) = b1p2(i,j,k,m)
 10               continue
      
      return
      end
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      SUBROUTINE crossf(aa,bbmf,cc)
c The cross product is formed at the main cell center.  a is assumed
c be main cell contravarient (cell face) and b is assumed to be
c main cell covarient (cell edge).  The result is main cell
c contravarient (cell face).

c Can only center vectors on main cell for loops 3 to n...and can
c only extrapolate back for loops 2 to n-1.  Must handle other cases
c separately.

c The magnetic field does not require extrapolation to cell centers
c on boundaries since dB/dx = 0 is the boundary condition.  That is
c just copy the interior values to the boundary.
c----------------------------------------------------------------------
CVD$R VECTOR

      include 'incurv.h'

      real aa(nx,ny,nz,3)        !main cell contravarient vector 
      real bbmf(nx,ny,nz,3)      !main cell contravarient vector
      real cc(nx,ny,nz,3)        !cross product result, main cell
                                 !contravarient (cell face)

      real ax,ay,az,bx,by,bz    !dummy vars
      real temp                 !used to vectorize loop
      real zfrc(nz)             !0.5*dz_grid(k)/dz_cell(k)
c      real ct(nx,ny,nz,3)       !temp main cell center cross product
      real aac(3),bbc(3)


c extrapolate(/interpolate) to main cell center and do cross product


      call periodic(aa)
      call periodic(bbmf)


      do 5 k=1,nz
         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
 5       continue



      do 10 k=2,nz-1      
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               im = i-1         !assume daa/dxyz = 0 at boundary
               jm = j-1         !bbmf is given on boundary
               km = k-1

               ax = 0.5*(aa(i,j,k,1) + aa(im,j,k,1))
               bx = 0.5*(bbmf(i,j,k,1) + bbmf(im,j,k,1))

               ay = 0.5*(aa(i,j,k,2) + aa(i,jm,k,2))
               by = 0.5*(bbmf(i,j,k,2) + bbmf(i,jm,k,2))

               az = zfrc(k)*(aa(i,j,k,3) - aa(i,j,km,3)) + aa(i,j,km,3)
               bz = zfrc(k)*(bbmf(i,j,k,3) - bbmf(i,j,km,3))
     x                     + bbmf(i,j,km,3)

               ct(i,j,k,1) = ay*bz - az*by
               ct(i,j,k,2) = az*bx - ax*bz
               ct(i,j,k,3) = ax*by - ay*bx

 10            continue

       call periodic(ct)

c extrapolate back to main cell contravarient positions.
c ...just average across cells since cell edges are centered
c about the grid points.
      
      do 60 k=2,nz-1
         do 60 j=2,ny-1
            do 60 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1

               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

 60            continue

      call periodic(cc)


      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE cov_to_contra(bt,btmf)
c Converts total magnetic field from main cell covarient positions
c to main cell contravarient positions.  This is then used in the
c fluid velocity update routines.  This routine assumes that cell 
c edges and cell centers are "topologically centered".  So the grid
c points do not reside at the cell centers...rather they are offset
c a little so that the cell edges are equidistant from the k and k-1
c grid points.  In extrapolating the coventient vectors to the 
c contravarient vector positions, this assymetry is accounted for
c using a linear interpolation of the k and k-1 values to the grid
c point location.
c----------------------------------------------------------------------
CVD$R VECTOR
      include 'incurv.h'

      real bt(nx,ny,nz,3),   !main cell covarient
     x     btmf(nx,ny,nz,3)  !main cell contravarient

      real bx1, bx2, by1, by2, bz1, bz2  !main cell center fields
      real zrat           !ratio for doing linear interpolation
                          !to grid point position.
      real zplus, zminus  !position of main cell edges up and down
      real b_j, b_jm, b_i, b_im !intermediate step in average process

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1
               im = i-1
               jm = j-1
               km = k-1

c The x component of B resides at the k and k-1 edges, so this
c requires the non-uniform grid interpolation

               zplus = (qz(k+1) + qz(k))/2.0
               zminus = (qz(k) + qz(k-1))/2.0
               zrat = (qz(k) - zminus)/(zplus - zminus)
   
               b_j = bt(i,j,km,1) 
     x               + zrat*(bt(i,j,k,1) - bt(i,j,km,1)) 
               b_jm = bt(i,jm,km,1)
     x                + zrat*(bt(i,jm,k,1) - bt(i,jm,km,1))
               bx1 = (b_j + b_jm)/2.0

               b_j = bt(ip,j,km,1) 
     x               + zrat*(bt(ip,j,k,1) - bt(ip,j,km,1)) 
               b_jm = bt(ip,jm,km,1)
     x                + zrat*(bt(ip,jm,k,1) - bt(ip,jm,km,1))
               bx2 = (b_j + b_jm)/2.0

               
               b_i = bt(i,j,km,2) 
     x               + zrat*(bt(i,j,k,2) - bt(i,j,km,2)) 
               b_im = bt(im,j,km,2)
     x                + zrat*(bt(im,j,k,2) - bt(im,j,km,2))           
               by1 = (b_i + b_im)/2.0

               b_i = bt(i,jp,km,2) 
     x               + zrat*(bt(i,jp,k,2) - bt(i,jp,km,2)) 
               b_im = bt(im,jp,km,2)
     x                + zrat*(bt(im,jp,k,2) - bt(im,jp,km,2))
               by2 = (b_i + b_im)/2.0


               bz1 = 0.25*(bt(i,j,k,3) + bt(i,jm,k,3) +
     x                     bt(im,jm,k,3) + bt(im,j,k,3))
               bz2 = 0.25*(bt(i,j,kp,3) + bt(i,jm,kp,3) +
     x                     bt(im,jm,kp,3) + bt(im,j,kp,3))

               btmf(i,j,k,1) = 0.5*(bx1+bx2)
               btmf(i,j,k,2) = 0.5*(by1+by2)
               btmf(i,j,k,3) = 0.5*(bz1+bz2)

 10            continue

c      call boundaries(btmf)
      call periodic(btmf)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE curlB(b1,np,aj)
c Calculates curl B / n*alpha.  The resulting "current" is called aj
c which is used in several other places in the code.  This curl is 
c performed on the main cell where B is covarient.  The resulting
c current is main cell contravarient.  Note that dz_cell is used for
c the cell dimensions since dz_grid is not equal to dz_cell on non-
c uniform grid.
c----------------------------------------------------------------------
CVD$R VECTOR
      include 'incurv.h'

      real b1(nx,ny,nz,3),
c     x     nf(nx,ny,nz),
     x     np(nx,ny,nz),
     x     aj(nx,ny,nz,3)

      real curl_B(3)      !dummy for holding curl vector
      real ntot(3)        !total density, np + nf

c      call periodic_scalar(np)
c      call periodic_scalar(nf)
      call periodic(b1)
cc     call fix_normal_b(b1)

      do 10 k=2,nz-1   
         do 10 j=2,ny-1
            do 10 i=2,nx-1

               ip = i+1
               jp = j+1
               kp = k+1

c               if (ip .gt. nx) then ip = nx
c               if (jp .gt. ny) then jp = ny
c               if (kp .gt. nz) then kp = nz

c               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
c     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
c               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
c     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
c               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
c     x                 + 0.5*(np(i,j,k)+np(i,j,kp))


               ntot(1) = 0.5*(np(i,j,k)+np(ip,j,k))
               ntot(2) = 0.5*(np(i,j,k)+np(i,jp,k))
               ntot(3) = 0.5*(np(i,j,k)+np(i,j,kp))

               curl_B(1) = (b1(i,j,k,3)/dy) - (b1(i,j-1,k,3)/dy) 
     x                    + (b1(i,j,k-1,2)/dz_cell(k))  
     x                    - (b1(i,j,k,2)/dz_cell(k))
               curl_B(2) = (b1(i,j,k,1)/dz_cell(k)) 
     x                     - (b1(i,j,k-1,1)/dz_cell(k))
     x                     - (b1(i,j,k,3)/dx) + (b1(i-1,j,k,3)/dx)
               curl_B(3) = (b1(i,j,k,2)/dx) - (b1(i-1,j,k,2)/dx) 
     x                     + (b1(i,j-1,k,1)/dy) - (b1(i,j,k,1)/dy)

               do 10 m=1,3
                  aj(i,j,k,m) = curl_B(m)/(ntot(m)*alpha)
 10            continue

c      call periodic(aj)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE curlE(E,curl_E)
c E is dual cell covarient, and curl_E will be returned as main
c cell covarient...as all magnetic fields are.  All i,j,k exclude
c boundaries.  Boundaries are taken care of in main fluid code.
c----------------------------------------------------------------------
      include 'incurv.h'

      real E(nx,ny,nz,3)      !E field, main cell contravarient
      real curl_E(nx,ny,nz,3) !curl of E, main cell covarient
      real lx, ly, lz         !lengths of dual cell edges

c      call periodic(E)

      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1

               lx = qx(i+1) - qx(i)
               ly = qy(j+1) - qy(j)
               lz = qz(k+1) - qz(k)

               curl_E(i,j,k,1) =  (E(i,j+1,k,3)/ly) - (E(i,j,k,3)/ly)
     x                       + (E(i,j,k,2)/lz) - (E(i,j,k+1,2)/lz)
               curl_E(i,j,k,2) =  (E(i,j,k,3)/lx) - (E(i+1,j,k,3)/lx)
     x                       + (E(i,j,k+1,1)/lz) - (E(i,j,k,1)/lz)
               curl_E(i,j,k,3) =  (E(i,j,k,1)/ly) - (E(i,j+1,k,1)/ly)
     x                       + (E(i+1,j,k,2)/lx) - (E(i,j,k,2)/lx)

 10          continue

c      call periodic(curl_E)

      return
      end
c----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_um_dot_BB(u,b,cc)
cc uf and btmf are gathered at main cell center and uf.B*B 
cc calculated.  Result returned to main cell contravarient
cc postion.
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real u(nx,ny,nz,3),  !main cell contravarient fluid velocity (uf)
c     x     b(nx,ny,nz,3),   !main cell contravarient bt (btmf)
c     x     cc(nx,ny,nz,3)   !(uf.B)*B

c      real ux,uy,uz,bx,by,bz   !temp vars for u and b at cell center
c      real temp                !used to vectorize loop
cc      real ct(nx,ny,nz,3)      !result are main cell center
c      real udotb               !u dot b
c      real zfrc(nz)            !0.5*dz_grid(k)/dz_cell(k)

c! first gather everything at center

c      call periodic(u)
c      call periodic(b)
cc      call fix_normal_b(b)

c      do 5 k=1,nz
c         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
c 5       continue

c      do 10 k=2,nz-1      
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1

c               im = i-1
c               jm = j-1     
c               km = k-1

c               ux = 0.5*(u(i,j,k,1) + u(im,j,k,1))
c               bx = 0.5*(b(i,j,k,1) + b(im,j,k,1))

c               uy = 0.5*(u(i,j,k,2) + u(i,jm,k,2))
c               by = 0.5*(b(i,j,k,2) + b(i,jm,k,2))

c               uz = zfrc(k)*(u(i,j,k,3) - u(i,j,km,3)) + u(i,j,km,3)
c               bz = zfrc(k)*(b(i,j,k,3) - b(i,j,km,3)) + b(i,j,km,3)

c               udotb = ux*bx + uy*by + uz*bz

c               ct(i,j,k,1) = udotb*bx
c               ct(i,j,k,2) = udotb*by
c               ct(i,j,k,3) = udotb*bz

c 10            continue

c      call periodic(ct)


cc extrapolate back to main cell contravarient positions.
cc ...just average across cells.

c      do 60 k=2,nz-1
c         do 60 j=2,ny-1
c            do 60 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

c               cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
c               cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
c               cc(i,j,k,3) = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))

c 60            continue

c      call periodic(cc)

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_um_dot_BB_old(u,b,cc)
cc uf and btmf are gathered at main cell center and uf.B*B 
cc calculated.  Result returned to main cell contravarient
cc postion.
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real u(nx,ny,nz,3),  !main cell contravarient fluid velocity (uf)
c     x     b(nx,ny,nz,3),   !main cell contravarient bt (btmf)
c     x     cc(nx,ny,nz,3)   !(uf.B)*B

c      real ux,uy,uz,bx,by,bz   !temp vars for u and b at cell center
c      real temp                !used to vectorize loop
c      real ct(nx,ny,nz,3)      !result are main cell center
c      real udotb               !u dot b
c      real zfrc(nz)            !0.5*dz_grid(k)/dz_cell(k)

c! first gather everything at center

c      call periodic(u)
c      call periodic(b)
cc      call fix_normal_b(b)

c      do 5 k=1,nz
c         zfrc(k) = 0.5*dz_grid(k)/dz_cell(k)
c 5       continue

c      do 10 i=2,nx-1      
c         do 10 j=2,ny-1
c            do 10 k=2,nz-1

c               im = i-1
c               jm = j-1     
c               km = k-1

c               if (i .eq. 2) then 
c                  ux = 2.0*u(2,j,k,1) - 
c     x                 0.5*(u(3,j,k,1) + u(2,j,k,1))
c                  bx = 2.0*b(2,j,k,1) - 
c     x                 0.5*(b(3,j,k,1) + b(2,j,k,1))
c               else
c                  ux = 0.5*(u(i,j,k,1) + u(im,j,k,1))
c                  bx = 0.5*(b(i,j,k,1) + b(im,j,k,1))
c               endif

c               if (j .eq. 2) then 
c                  uy = 2.0*u(i,2,k,2) - 
c     x                 0.5*(u(i,3,k,2) + u(i,2,k,2))
c                  by = 2.0*b(i,2,k,2) - 
c     x                 0.5*(b(i,3,k,2) + b(i,2,k,2))
c               else
c                  uy = 0.5*(u(i,j,k,2) + u(i,jm,k,2))
c                  by = 0.5*(b(i,j,k,2) + b(i,jm,k,2))
c               endif

c               if (k .eq. 2) then
c                  uz = 2.0*u(i,j,2,3) - 
c     x                 zfrc(k)*(u(i,j,3,3) - u(i,j,2,3)) + u(i,j,2,3) 
cc                  uz = 2.0*u(i,j,2,3) - 
cc     x                 0.5*(u(i,j,3,3) + u(i,j,2,3))
c                  bz = b(i,j,2,3)
cc                  bz = 2.0*b(i,j,2,3) - 
cc     x                 0.5*(b(i,j,3,3) + b(i,j,2,3))
c               else
c                  uz = zfrc(k)*(u(i,j,k,3) - u(i,j,km,3)) + u(i,j,km,3)
c                  bz = zfrc(k)*(b(i,j,k,3) - b(i,j,km,3)) + b(i,j,km,3)
cc                  uz = 0.5*(u(i,j,k,3) + u(i,j,km,3))
cc                  bz = 0.5*(b(i,j,k,3) + b(i,j,km,3))            
c               endif

c               udotb = ux*bx + uy*by + uz*bz

c               ct(i,j,k,1) = udotb*bx
c               ct(i,j,k,2) = udotb*by
c               ct(i,j,k,3) = udotb*bz

c 10            continue

c      call periodic(ct)

cc extrapolate back to main cell contravarient positions.
cc ...just average across cells.

c      do 60 i=2,nx
c         do 60 j=2,ny
c            do 60 k=2,nz

c               ip = i+1
c               jp = j+1
c               kp = k+1

c               if (i .eq. nx-1) then 
c                  cc(nx,j,k,1) = 2.0*ct(nx,j,k,1) - 
c     x                           0.5*(ct(nx,j,k,1) + ct(nx-1,j,k,1))
c               else
c                  cc(i,j,k,1) = 0.5*(ct(i,j,k,1) + ct(ip,j,k,1))
c               endif

c               if (j .eq. ny-1) then 
c                  cc(i,ny,k,2) = 2.0*ct(i,ny,k,2) - 
c     x                           0.5*(ct(i,ny,k,2) + ct(i,ny-1,k,2))
c               else
c                  cc(i,j,k,2) = 0.5*(ct(i,j,k,2) + ct(i,jp,k,2))
c               endif
                  
c               if (k .eq. nz-1) then
cc                  temp = 2.0*ct(i,j,nz,3) - 
cc     x                           0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3))
c                   temp = 0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)) +
c     x                   (2.0*dz_cell(nz)/dz_grid(nz))*(ct(i,j,nz,3) -
c     x                    0.5*(ct(i,j,nz,3) + ct(i,j,nz-1,3)))
c               else
c                  temp = 0.5*(ct(i,j,k,3) + ct(i,j,kp,3))
c               endif             
  
c               cc(i,j,k,3) = temp

c 60            continue

c      call periodic(cc)

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_ugradu_Lax(uf,ugradu,delta_t)
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real uf(nx,ny,nz,3),
c     x     ugradu(nx,ny,nz,3)

c      real ufc(nx,ny,nz,3)
c      real ax1,ax2,ay1,ay2,az1,az2       !misc const
c      real u1,u2,u3                      !temp vars

c      parameter(ad = 0.001)                 !coefficient to add extra
                                         !diffusion
c      call periodic(uf)

c      call face_to_center(uf,ufc)

c      do 10 k=2,nz-1
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1
c               ip=i+1
c               jp=j+1
c               kp=k+1
c               im=i-1
c               jm=j-1
c               km=k-1

cc xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c               ax1 = 0.5*ufc(i,j,k,1)/dx
c               ax2 = ad*abs(ufc(ip,j,k,1) - ufc(im,j,k,1))
c               ay1 = 0.5*ufc(i,j,k,2)/dy
c               ay2 = ad*abs(ufc(ip,j,k,2) - ufc(im,j,k,2))
c               u1 = ax1*(ufc(ip,j,k,1) - ufc(im,j,k,1)) - 
c     x              ax2*(ufc(im,j,k,1) - 2.0*ufc(i,j,k,1) +
c     x                      ufc(ip,j,k,1))
c               u2 = ay1*(ufc(i,jp,k,1) - ufc(i,jm,k,1)) - 
c     x              ay2*(ufc(i,jm,k,1) - 2.0*ufc(i,j,k,1) +
c     x                      ufc(i,jp,k,1)) 
c               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
c               az2 = ad*abs(ufc(ip,j,k,3) - ufc(im,j,k,3))
c               u3 = az1*(ufc(i,j,kp,1)-ufc(i,j,km,1)) -
c     x              az2*(ufc(i,j,km,1) - 2.0*ufc(i,j,k,1) +
c     x                   ufc(i,j,kp,1))
c               ct(i,j,k,1) = u1 + u2 + u3

cc yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

c               ax1 = 0.5*ufc(i,j,k,1)/dx
c               ax2 = ad*abs(ufc(i,jp,k,1) - ufc(i,jm,k,1))
c               ay1 = 0.5*ufc(i,j,k,2)/dy
c               ay2 = ad*abs(ufc(i,jp,k,2) - ufc(i,jm,k,2))
c               u1 = ax1*(ufc(ip,j,k,2) - ufc(im,j,k,2)) - 
c     x              ax2*(ufc(im,j,k,2) - 2.0*ufc(i,j,k,2) +
c     x                      ufc(ip,j,k,2))
c               u2 = ay1*(ufc(i,jp,k,2) - ufc(i,jm,k,2)) - 
c     x              ay2*(ufc(i,jm,k,2) - 2.0*ufc(i,j,k,2) +
c     x                      ufc(i,jp,k,2)) 
c               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
c               az2 = ad*abs(ufc(i,jp,k,3) - ufc(i,jm,k,3))
c               u3 = az1*(ufc(i,j,kp,2)-ufc(i,j,km,2)) -
c     x              az2*(ufc(i,j,km,2) - 2.0*ufc(i,j,k,2) +
c     x                   ufc(i,j,kp,2))
c               ct(i,j,k,2) = u1 + u2 + u3

cc zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
 
c               ax1 = 0.5*ufc(i,j,k,1)/dx
c               ax2 = ad*abs(ufc(i,j,kp,1) - ufc(i,j,km,1))
c               ay1 = 0.5*ufc(i,j,k,2)/dy
c               ay2 = ad*abs(ufc(i,j,kp,2) - ufc(i,j,km,2))
c               u1 = ax1*(ufc(ip,j,k,3) - ufc(im,j,k,3)) - 
c     x              ax2*(ufc(im,j,k,3) - 2.0*ufc(i,j,k,3) +
c     x                      ufc(ip,j,k,3))
c               u2 = ay1*(ufc(i,jp,k,3) - ufc(i,jm,k,3)) - 
c     x              ay2*(ufc(i,jm,k,3) - 2.0*ufc(i,j,k,3) +
c     x                      ufc(i,jp,k,3)) 
c               az1 = ufc(i,j,k,3)/(dz_grid(k)+dz_grid(kp))
c               az2 = ad*abs(ufc(i,j,kp,3) - ufc(i,j,km,3))
c               u3 = az1*(ufc(i,j,kp,3)-ufc(i,j,km,3)) -
c     x              az2*(ufc(i,j,km,3) - 2.0*ufc(i,j,k,3) +
c     x                   ufc(i,j,kp,3))
c               ct(i,j,k,3) = u1 + u2 + u3

c 10            continue

c      call periodic(ct)

cc interpolate back to contravarient positions.

c      do 20 k=2,nz-1
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1
c               ip = i+1
c               jp = j+1
c               kp = k+1
c               ugradu(i,j,k,1) = 0.5*(ct(i,j,k,1)+ct(ip,j,k,1))
c               ugradu(i,j,k,2) = 0.5*(ct(i,j,k,2)+ct(i,jp,k,2))
c               ugradu(i,j,k,3) = 0.5*(ct(i,j,k,3)+ct(i,j,kp,3))
c 20         continue

c       call periodic(ugradu)

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE get_Ef(Ef,aj,np,nf,up,uf,btmf,nu,ugradu,delta_t,
c     x                  gradPf)
cc Need to treat boundaries separately!!
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real Ef(nx,ny,nz,3),
c     x     aj(nx,ny,nz,3),
c     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
c     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     btmf(nx,ny,nz,3),
c     x     nu(nx,ny,nz),
c     x     ugradu(nx,ny,nz,3),
c     x     gradPf(nx,ny,nz,3)
 
c      real ntot(3)                !total plasma density
c      real fnp(3)                 !fraction, np/n
c      real aac(3),bbc(3),ccc(3)
c      real cc(nx,ny,nz,3)

c      call periodic_scalar(np)
c      call periodic_scalar(nf)

c      do 10 k=2,nz-1 
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

cc               if (ip .gt. nx) then ip = nx
cc               if (jp .gt. ny) then jp = ny
cc               if (kp .gt. nz) then kp = nz

c               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
c     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
c               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
c     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
c               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
c     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

c               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
c               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
c               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)
               
c               do 10 m=1,3
c                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m)
c 10               continue

c      call crossf(a,btmf,c)

c      call get_ugradu_Lax(uf,ugradu,delta_t)

c      do 20 k=2,nz-1 
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

cc               if (ip .gt. nx) then ip = nx
cc               if (jp .gt. ny) then jp = ny
cc               if (kp .gt. nz) then kp = nz

c               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))
c     x                 + 0.5*(np(i,j,k)+np(ip,j,k))
c               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))
c     x                 + 0.5*(np(i,j,k)+np(i,jp,k))
c               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))
c     x                 + 0.5*(np(i,j,k)+np(i,j,kp))

c               fnp(1) = 0.5*(np(i,j,k)+np(ip,j,k))/ntot(1)
c               fnp(2) = 0.5*(np(i,j,k)+np(i,jp,k))/ntot(2)
c               fnp(3) = 0.5*(np(i,j,k)+np(i,j,kp))/ntot(3)

c               do 20 m=1,3
c                  Ef(i,j,k,m) = c(i,j,k,m) - ugradu(i,j,k,m) 
c     x                          + nu(i,j,k)*fnp(m)*up(i,j,k,m)
c     x                          + nuei*aj(i,j,k,m) - gradPf(i,j,k,m)
cc     x                          + etar(i,j,k,m)*aj(i,j,k,m)
cc                  Ef(i,j,k,m) = c(i,j,k,m) + 
cc     x                          nu(i,j,k)*fnp(m)*up(i,j,k,m)
c 20            continue

c      call periodic(Ef)
cc      call fix_tangential_E(Ef)

c      Ef(nx-1:nx,:,:,3) = 0.0
c      Ef(nx-1:nx,:,:,2) = 0.0
c      Ef(nx-1:nx,:,:,1) = 0.0

c      return
c      end
cc----------------------------------------------------------------------



cc----------------------------------------------------------------------
c      SUBROUTINE get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
c     x                            delta_t)
cc This is the heart of the fluid velocity update.  It solves eqn. 18
cc (Dan's paper) for uf+
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real Ef(nx,ny,nz,3),
c     x     btmf(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3), !using particle velocities at t level n-1/2
c     x     nu(nx,ny,nz),
c     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
c     x     uplus(nx,ny,nz,3),
c     x     uminus(nx,ny,nz,3)

c      real a1,b1,b2,b3      !4 coefficients (see notebook solution)
c      real PP,QQ            !intermediate variables for a1,b1,b2,b3
c      real eta              !intermediate variable for PP,QQ
cc      real B(3)             !B for cross product call
cc      real Bsqrd                     !B*B
c      real um_x_B(nx,ny,nz,3)        !uf- X B
c      real um_dot_BB(nx,ny,nz,3)     !uf- . B
c      real ntot(3)                   !total density np + nf
c      real npave(3)
c      real btc(nx,ny,nz,3)
c      real bsqrd(nx,ny,nz),bsq(3)

c      do 10 i=2,nx-1    
c         do 10 j=2,ny-1
c            do 10 k=2,nz-1
c               do 10 m=1,3
c                  uminus(i,j,k,m) = uf(i,j,k,m) + 
c     x                              0.5*delta_t*Ef(i,j,k,m)
c 10            continue

c      call crossf(uminus, btmf, um_x_B)
c      call get_um_dot_BB(uminus , btmf, um_dot_BB)

c      call face_to_center(btmf,btc)

c      do 15 k=2,nz-1 
c         do 15 j=2,ny-1
c            do 15 i=2,nx-1
c               bsqrd(i,j,k) =  btc(i,j,k,1)**2 + btc(i,j,k,2)**2 + 
c     x               btc(i,j,k,3)**2
c 15            continue

c      call periodic_scalar(bsqrd)
              
c      do 20 k=2,nz-1
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

cc               if (ip .gt. nx) then ip = nx
cc               if (jp .gt. ny) then jp = ny
cc               if (kp .gt. nz) then kp = nz

c               bsq(1) = 0.5*(bsqrd(i,j,k) + bsqrd(ip,j,k))
c               bsq(2) = 0.5*(bsqrd(i,j,k) + bsqrd(i,jp,k))
c               bsq(3) = 0.5*(bsqrd(i,j,k) + bsqrd(i,j,kp))

c               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
c               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
c               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))

c               ntot(1) = 0.5*(nf(i,j,k)+nf(ip,j,k)) + npave(1)
c               ntot(2) = 0.5*(nf(i,j,k)+nf(i,jp,k)) + npave(2)
c               ntot(3) = 0.5*(nf(i,j,k)+nf(i,j,kp)) + npave(3)

c               do 20 m=1,3

c                  eta = npave(m)*delta_t/(2.0*ntot(m))
c                  QQ = eta/(1.0+nu(i,j,k)*eta)
c                  PP = (1.0-nu(i,j,k)*eta)/(1.0+nu(i,j,k)*eta)
c                  a1 = 1.0/(1.0 + QQ*QQ*bsq(m))
c                  b1 = PP - (QQ*QQ*bsq(m))
c                  b2 = (QQ*PP) + QQ
c                  b3 = (QQ*QQ*PP) + (QQ*QQ)

c                  uplus(i,j,k,m) = a1*(b1*uminus(i,j,k,m) + 
c     x                             b2*um_x_B(i,j,k,m) + 
c     x                             b3*um_dot_BB(i,j,k,m))

c 20            continue

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE predict_uf(Ef,b0,b1,b12,uf,uf2,ufp2,nu,np,nf,uplus, 
c     x                      uminus,ugradu,up,gradP,nuin,bdp,pf1)
cc Calculate the fluid velocity, uf,  at the new time step and replace
cc uf1 with the new value, uf, in preparation for the next time step.
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real Ef(nx,ny,nz,3),
c     x     b0(ny),
c     x     b1(nx,ny,nz,3),
c     x     b12(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     uf2(nx,ny,nz,3),
c     x     ufp2(nx,ny,nz,3),
c     x     nu(nx,ny,nz),
c     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
c     x     uplus(nx,ny,nz,3), 
c     x     uminus(nx,ny,nz,3),
c     x     ugradu(nx,ny,nz,3),
c     x     up(nx,ny,nz,3),
cc     x     gradP(nx,ny,nz,3),
c     x     nuin(nx,ny,nz),
c     x     bdp(nx,ny,nz,3),
c     x     pf1(nx,ny,nz)

c      real b1h(nx,ny,nz,3)
c      real bth(nx,ny,nz,3)
c      real btmfh(nx,ny,nz,3)
c      real ajh(nx,ny,nz,3)
c      real gradPf(nx,ny,nz,3)

c      real m_den
c      real delta_t

c      delta_t = 2.0*dtsub

c      do k=2,nz-1
c         do j=2,ny-1
c            do i=2,nx-1
c               m_den = mO*0.5*(nf(i+1,j,k)+nf(i,j,k))
c               gradPf(i,j,k,1) = (pf1(i+1,j,k)-pf1(i,j,k))/(dx*m_den)
c               m_den = mO*0.5*(nf(i,j+1,k)+nf(i,j,k))
c               gradPf(i,j,k,2) = (pf1(i,j+1,k)-pf1(i,j,k))/(dy*m_den)
c               m_den = mO*0.5*(nf(i,j,k+1)+nf(i,j,k))
c               gradPf(i,j,k,3) = (pf1(i,j,k+1)-
c     x                             pf1(i,j,k))/(dz_grid(k)*m_den)
c            enddo
c         enddo
c      enddo


c      do 10 k=1,nz
c         do 10 j=1,ny
c            do 10 i=1,nx
c               b1h(i,j,k,1) = 0.5*(b1(i,j,k,1) + b12(i,j,k,1))
c               b1h(i,j,k,2) = 0.5*(b1(i,j,k,2) + b12(i,j,k,2))
c               b1h(i,j,k,3) = 0.5*(b1(i,j,k,3) + b12(i,j,k,3))
c               bth(i,j,k,1) = b1h(i,j,k,1) + bdp(i,j,k,1)
c               bth(i,j,k,2) = b1h(i,j,k,2) + b0(j) + bdp(i,j,k,2)
c               bth(i,j,k,3) = b1h(i,j,k,3) + bdp(i,j,k,3)
c 10            continue

c      call cov_to_contra(bth,btmfh)
c      call curlB(b1h,nf,np,ajh)

c      call get_Ef(Ef,ajh,np,nf,up,uf,btmfh,nu,ugradu,delta_t,gradPf)
c      call get_uplus_uminus(Ef,btmfh,uf2,nu,np,nf,uplus,uminus,
c     x                      delta_t)

c      do 20 k=2,nz-1
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1
c               do 20 m=1,3
c                  ufp2(i,j,k,m) = uplus(i,j,k,m) + 
c     x                            0.5*delta_t*Ef(i,j,k,m) !-
cc     x                        0.5*delta_t*nuin(i,j,k)*uplus(i,j,k,m)
c 20            continue

cc      ufp2(nx-1:nx,:,:,1) = -vsw

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE correct_uf(Ef,btmf,uf,uf2,ufp2,nu,np,nf,uplus,uminus, 
c     x                      ugradu,aj,up,ufp1,gradP,nuin,pf)
cc Calculate the fluid velocity, uf,  at the new time step and replace
cc uf1 with the new value, uf, in preparation for the next time step.
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real Ef(nx,ny,nz,3),
c     x     btmf(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     uf2(nx,ny,nz,3),
c     x     ufp2(nx,ny,nz,3),
c     x     nu(nx,ny,nz),
c     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
c     x     uplus(nx,ny,nz,3), 
c     x     uminus(nx,ny,nz,3),
c     x     ugradu(nx,ny,nz,3),
c     x     aj(nx,ny,nz,3),
c     x     up(nx,ny,nz,3),
c     x     ufp1(nx,ny,nz,3),
cc     x     gradP(nx,ny,nz,3),
c     x     nuin(nx,ny,nz),
c     x     pf(nx,ny,nz)
 
c      real m_den
c      real delta_t
c      real gradPf(nx,ny,nz,3)

c      delta_t = dtsub

c      do k=2,nz-1
c         do j=2,ny-1
c            do i=2,nx-1
c               m_den = mO*0.5*(nf(i+1,j,k)+nf(i,j,k))
c               gradPf(i,j,k,1) = (pf(i+1,j,k)-pf(i,j,k))/(dx*m_den)
c               m_den = mO*0.5*(nf(i,j+1,k)+nf(i,j,k))
c               gradPf(i,j,k,2) = (pf(i,j+1,k)-pf(i,j,k))/(dy*m_den)
c               m_den = mO*0.5*(nf(i,j,k+1)+nf(i,j,k))
c               gradPf(i,j,k,3) = (pf(i,j,k+1)-
c     x                             pf(i,j,k))/(dz_grid(k)*m_den)
c            enddo
c         enddo
c      enddo

      
c      do 10 k=2,nz-1
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1
c               do 10 m=1,3
c                  ufp1(i,j,k,m) = 0.5*(uf(i,j,k,m) + ufp2(i,j,k,m))
c 10              continue

c      call get_Ef(Ef,aj,np,nf,up,ufp1,btmf,nu,ugradu,delta_t,gradPf)
c      call get_uplus_uminus(Ef,btmf,uf,nu,np,nf,uplus,uminus,
c     x                      delta_t)

c      do 20 k=2,nz-1
c         do 20 j=2,ny-1
c            do 20 i=2,nx-1
c               do 20 m=1,3
c                  uf2(i,j,k,m) = uf(i,j,k,m)
c                  uf(i,j,k,m) = uplus(i,j,k,m) + 0.5*dtsub*Ef(i,j,k,m)
cc     x                          - 0.5*dtsub*nuin(i,j,k)*uplus(i,j,k,m)
c 20            continue

cc      uf(nx-1:nx,:,:,1) = -vsw

c      return
c      end
cc----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_E(E,b0,bt,btmf,aj,up,np,nu)
c E must be at time level m. We have uf at levels m-1/2 and m+1/2, so
c the average value is used for uf in the calculation of ui.
c----------------------------------------------------------------------
CVD$R VECTOR
      include 'incurv.h'

      real E(nx,ny,nz,3),
     x     b0(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     btmf(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     uf2(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3)

      real ntot(3)         !total density np + nf
      real fnp(3),fnf(3)   !fraction np and nf of n
      real npave(3)

c      real a(nx,ny,nz,3), 
c     x     c(nx,ny,nz,3)  !dummy vars for doing cross product

      call periodic_scalar(np)
c      call periodic_scalar(nf)

      do 10 k=2,nz-1    
         do 10 j=2,ny-1
            do 10 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

c               if (ip .gt. nx) then ip = nx
c               if (jp .gt. ny) then jp = ny
c               if (kp .gt. nz) then kp = nz

c               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
c               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
c               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))

c               ntot(1) = npave(1) + 0.5*(nf(i,j,k)+nf(ip,j,k))
c               ntot(2) = npave(2) + 0.5*(nf(i,j,k)+nf(i,jp,k))
c               ntot(3) = npave(3) + 0.5*(nf(i,j,k)+nf(i,j,kp))
               
c               fnp(1) = npave(1)/ntot(1)
c               fnp(2) = npave(2)/ntot(2)
c               fnp(3) = npave(3)/ntot(3)

c               fnf(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))/ntot(1)
c               fnf(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))/ntot(2)
c               fnf(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))/ntot(3)

c               ntot = np(i,j,k) + nf(i,j,k)
c               fnp = np(i,j,k)/ntot
c               fnf = nf(i,j,k)/ntot

               do 10 m=1,3
                  a(i,j,k,m) = aj(i,j,k,m) - up(i,j,k,m)
c                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m) - 
c     x                         fnf(m)*0.5*(uf2(i,j,k,m)+uf(i,j,k,m))
c                  a(i,j,k,m) = - fnp(m)*up(i,j,k,m) - 
c     x                         fnf(m)*0.5*(uf2(i,j,k,m)+uf(i,j,k,m))
 10               continue

      call crossf(a,btmf,c)



               
      do 20 k=2,nz-1      
         do 20 j=2,ny-1   
            do 20 i=2,nx-1
               do 20 m=1,3 
                  E(i,j,k,m) = c(i,j,k,m) + nu(i,j,k)*aj(i,j,k,m)
c     x                         + nuei*aj(i,j,k,m) !- gradP(i,j,k,m)
c     x                         + etar(i,j,k,m)*aj(i,j,k,m)
 20               continue


c      call fix_tangential_E(E)
      call periodic(E)

c      call fix_tangential_E(E)

c      E(:,:,nz,3) = -vsw*q*b0(:,:,nz,2)/mO
c      E(:,:,nz,2) = -vsw*q*b0(:,:,nz,3)/mO
c      E(:,:,nz,2) = 0.0
c      E(:,:,nz,1) = 0.0

c      E(:,:,1,3) = -vsw*q*b0(:,:,1,2)/mO
c      E(:,:,1,2) = -vsw*q*b0(:,:,1,3)/mO
c      E(:,:,1:2,2) = 0.0
c      E(:,:,1:2,1) = 0.0


c      E(nx-1:nx,:,:,3) = vsw*q*b0_init/mO

c      E(nx-1:nx,:,:,2) = 0.0
cc      E(nx-1:nx,:,:,1) = 0.0

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE predict_B(b0,b1,b12,b1p2,bt,btmf,E,aj,up,np,nu)
c Predictor step in magnetic field update.
c----------------------------------------------------------------------
CVD$R VECTOR
      include 'incurv.h'

      real b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b12(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     bt(nx,ny,nz,3),
     x     btmf(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
c     x     uf2(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)   !curl of E

      call get_E(E,b0,bt,btmf,aj,up,np,nu)  !E at time level m 


      call curlE(E,curl_E)
c      call fix_tangential_E(E)
c      call periodic(E)
c      call fix_tangential_E(E)

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1
               do 10 m=1,3
                  b1p2(i,j,k,m)=(1.0/12.0)*(b12(i+1,j,k,m)+
     x                 b12(i-1,j,k,m)+
     x                 b12(i,j+1,k,m)+b12(i,j-1,k,m)+
     x                 b12(i,j,k+1,m)+b12(i,j,k-1,m)+
     x                 6.0*b12(i,j,k,m)) -
     x                 2.0*dtsub*curl_E(i,j,k,m)


c                  b1p2(i,j,k,m) = b12(i,j,k,m) - 
c     x                            2.0*dtsub*curl_E(i,j,k,m)
 10               continue

c      call boundaries(b1p2)
c      call damp(b1p2)
      call periodic(b1p2)
c      call fix_normal_b(b1p2)


      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE get_Ep1(E,b0,b1,b1p2,aj,up,np,nu)
c The main feature here is that E must be calculated at time level
c m + 1/2.  That means that we need B at m + 1/2.  So b1p1 is
c calculated as 0.5*(b1 + b1p2).  uf and np are already at time level
c m + 1/2, so they are used as is. 
c----------------------------------------------------------------------
CVD$R VECTOR
      include 'incurv.h'

      real E(nx,ny,nz,3),
     x     b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3),
c     x     bdp(nx,ny,nz,3)

      real b1p1(nx,ny,nz,3)   !b1 at time level m + 1/2
      real btp1(nx,ny,nz,3)   !bt at time level m + 1/2
      real btp1mf(nx,ny,nz,3) !btp1 at contravarient position
      real ntot(3)            !total density np + nf
      real fnp(3),fnf(3)      !fraction np and nf of n
      real npave(3)

c      real a(nx,ny,nz,3),
c     x     c(nx,ny,nz,3)    !dummy vars for doing cross product


      do 5 k=1,nz
         do 5 j=1,ny
            do 5 i=1,nx
               btp1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1))
               b1p1(i,j,k,1) = 0.5*(b1p2(i,j,k,1) + b1(i,j,k,1))
               btp1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2))
               b1p1(i,j,k,2) = 0.5*(b1p2(i,j,k,2) + b1(i,j,k,2))
               btp1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3))
               b1p1(i,j,k,3) = 0.5*(b1p2(i,j,k,3) + b1(i,j,k,3))
 5             continue

      call curlB(btp1,np,aj)

      call periodic_scalar(np)
c      call periodic_scalar(nf)

      do 10 k=2,nz-1       
         do 10 j=2,ny-1
            do 10 i=2,nx-1

c               ip = i+1
c               jp = j+1
c               kp = k+1

c               if (ip .gt. nx) then ip = nx
c               if (jp .gt. ny) then jp = ny
c               if (kp .gt. nz) then kp = nz

c               npave(1) = 0.5*(np(i,j,k)+np(ip,j,k))
c               npave(2) = 0.5*(np(i,j,k)+np(i,jp,k))
c               npave(3) = 0.5*(np(i,j,k)+np(i,j,kp))

c               ntot(1) = npave(1) + 0.5*(nf(i,j,k)+nf(ip,j,k))
c               ntot(2) = npave(2) + 0.5*(nf(i,j,k)+nf(i,jp,k))
c               ntot(3) = npave(3) + 0.5*(nf(i,j,k)+nf(i,j,kp))
               
c               fnp(1) = npave(1)/ntot(1)
c               fnp(2) = npave(2)/ntot(2)
c               fnp(3) = npave(3)/ntot(3)

c               fnf(1) = 0.5*(nf(i,j,k)+nf(ip,j,k))/ntot(1)
c               fnf(2) = 0.5*(nf(i,j,k)+nf(i,jp,k))/ntot(2)
c               fnf(3) = 0.5*(nf(i,j,k)+nf(i,j,kp))/ntot(3)

c               ntot = np(i,j,k) + nf(i,j,k)
c               fnp = np(i,j,k)/ntot
c               fnf = nf(i,j,k)/ntot

               do 10 m=1,3
c                  a(i,j,k,m) = aj(i,j,k,m) - fnp(m)*up(i,j,k,m) - 
c     x                                       fnf(m)*uf(i,j,k,m)
                  a(i,j,k,m) = aj(i,j,k,m) - up(i,j,k,m)
c                  a(i,j,k,m) = - fnp(m)*up(i,j,k,m) - 
c     x                           fnf(m)*uf(i,j,k,m)
 10               continue

      call cov_to_contra(btp1,btp1mf)

      call crossf(a,btp1mf,c)

               
      do 20 k=2,nz-1       
         do 20 j=2,ny-1     
            do 20 i=2,nx-1  
               do 20 m=1,3 
                  E(i,j,k,m) = c(i,j,k,m) + nu(i,j,k)*aj(i,j,k,m)
c     x                         + nuei*aj(i,j,k,m) !- gradP(i,j,k,m)
c     x                         + etar(i,j,k,m)*aj(i,j,k,m)
 20               continue

c      call fix_tangential_E(E)
      call periodic(E)
c      call fix_tangential_E(E)

c      E(:,:,nz,3) = -vsw*q*b0(:,:,nz,2)/mO
c      E(:,:,nz,2) = 0.0
c      E(:,:,nz,1) = 0.0

c      E(:,:,1,3) = -vsw*q*b0(:,:,1,2)/mO
c      E(:,:,1:2,2) = 0.0
c      E(:,:,1:2,1) = 0.0

c      E(nx-1:nx,:,:,3) = vsw*q*b0_init/mO
c      E(nx-1:nx,:,:,2) = 0.0
cc      E(nx-1:nx,:,:,1) = 0.0

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE correct_B(b0,b1,b1p2,E,aj,up,np,nu)
c Corrector step in magnetic field update.
c----------------------------------------------------------------------
CVD$R VECTOR
      include 'incurv.h'

      real b0(nx,ny,nz,3),
     x     b1(nx,ny,nz,3),
     x     b1p2(nx,ny,nz,3),
     x     E(nx,ny,nz,3),
     x     aj(nx,ny,nz,3),
     x     up(nx,ny,nz,3),
c     x     uf(nx,ny,nz,3),
     x     np(nx,ny,nz),
c     x     nf(nx,ny,nz),
     x     nu(nx,ny,nz)
c     x     gradP(nx,ny,nz,3),
c     x     bdp(nx,ny,nz,3)

      real curl_E(nx,ny,nz,3)            !curl of E

      call get_Ep1(E,b0,b1,b1p2,aj,up,np,nu)  
                                                   !E at time level m 

      call curlE(E,curl_E)
c      call fix_tangential_E(E)
c      call periodic(E)
c      call fix_tangential_E(E)

c      write(*,*) 'E cb...',E(23,8,14,1),E(23,8,14,2),E(23,8,14,3)

      do 10 k=2,nz-1
         do 10 j=2,ny-1
            do 10 i=2,nx-1
               do 10 m=1,3
                  b1p2(i,j,k,m)=(1.0/12.0)*(b1(i+1,j,k,m)+b1(i-1,j,k,m)+
     x                 b1(i,j+1,k,m)+b1(i,j-1,k,m)+
     x                 b1(i,j,k+1,m)+b1(i,j,k-1,m)+6.0*b1(i,j,k,m)) -
     x                          dtsub*curl_E(i,j,k,m)
c                  b1p2(i,j,k,m) = b1(i,j,k,m) - 
c     x                            dtsub*curl_E(i,j,k,m)
 10               continue

c      call boundaries(b1p2)
c      call damp(b1p2)
      call periodic(b1p2)
c      call fix_normal_b(b1p2)

      return
      end
c----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE predict_nf(nf,nf1,nf3,nfp1,uf,divu,b1)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real nf(nx,ny,nz),
c     x     nf1(nx,ny,nz),
c     x     nf3(nx,ny,nz),
c     x     nfp1(nx,ny,nz),
c     x     uf(nx,ny,nz,3),
c     x     divu(nx,ny,nz),
c     x     b1(nx,ny,nz,3)
 
c      real minnf,maxnf
c      real flx(nx,ny,nz,3)
c      real ufc(nx,ny,nz,3)
c      real b1c(nx,ny,nz,3)

cc      call face_to_center(uf,ufc)
cc      call face_to_center(b1,b1c)
 
c      minnf = 10.0e20
c      maxnf = 10.0e20

c      call periodic_scalar(nf1)

c      do 5 i=2,nx-1
c         do 5 j=2,ny-1
c            do 5 k=2,nz-1
c               flx(i,j,k,1) = 0.5*(nf1(i,j,k)+nf1(i+1,j,k))*uf(i,j,k,1)
c               flx(i,j,k,2) = 0.5*(nf1(i,j,k)+nf1(i,j+1,k))*uf(i,j,k,2)
c               flx(i,j,k,3) = 0.5*(nf1(i,j,k)+nf1(i,j,k+1))*uf(i,j,k,3)
c 5             continue

c      call periodic(flx)

c      do 10 i=2,nx-1
c         do 10 j=2,ny-1
c            do 10 k=2,nz-1
c               nfp1(i,j,k) = nf3(i,j,k) 
c     x         - (2.0*dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
c     x         - (2.0*dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
c     x         - (2.0*dtsub/(dz_cell(k)))*
c     x           (flx(i,j,k,3)-flx(i,j,k-1,3))
c               divu(i,j,k) = (uf(i,j,k,1) - uf(i-1,j,k,1))/dx +
c     x                       (uf(i,j,k,2) - uf(i,j-1,k,2))/dy +
c     x                       (uf(i,j,k,3) - uf(i,j,k-1,3))/dz_cell(k)
c 10            continue

c      do 50 i=2,nx-1
c         do 50 j=2,ny-1
c            do 50 k=2,nz-1
c               nf(i,j,k) = 0.5*(nfp1(i,j,k) + nf1(i,j,k))
c               nf3(i,j,k) = nf1(i,j,k)
c 50            continue

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE correct_nf(nf,nf1,ufp1)
cc----------------------------------------------------------------------

c      include 'incurv.h'

c      real nf(nx,ny,nz),
c     x     nf1(nx,ny,nz),
c     x     ufp1(nx,ny,nz,3) 

c      real flx(nx,ny,nz,3)
c      real minnf,maxnf
c      real ufp1c(nx,ny,nz,3)

c      minnf = 10.0000000000000000e20
c      maxnf = 10.0000000000000000e20

cc      call face_to_center(ufp1,ufp1c)

c      do 5 i=2,nx-1
c         do 5 j=2,ny-1
c            do 5 k=2,nz-1
c             flx(i,j,k,1) = 0.5*(nf(i,j,k)+nf(i+1,j,k))*ufp1(i,j,k,1)
c             flx(i,j,k,2) = 0.5*(nf(i,j,k)+nf(i,j+1,k))*ufp1(i,j,k,2)
c             flx(i,j,k,3) = 0.5*(nf(i,j,k)+nf(i,j,k+1))*ufp1(i,j,k,3)
c 5           continue

c      call periodic(flx)

c      do 10 i=2,nx-1
c         do 10 j=2,ny-1
c            do 10 k=2,nz-1
c               nf1(i,j,k) = nf1(i,j,k) 
c     x            - (dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
c     x            - (dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
c     x            - (dtsub/(dz_cell(k)))*
c     x              (flx(i,j,k,3)-flx(i,j,k-1,3))
c               if (nf(i,j,k) .lt. minnf) then minnf = nf(i,j,k)
c               if (nf(i,j,k) .ge. maxnf) then maxnf = nf(i,j,k)
c 10            continue


c      call periodic_scalar(nf1)

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE trans_nf_Lax(nf,nf1,nfp1,uf)
cc----------------------------------------------------------------------
c      include 'incurv.h'

c      real nf(nx,ny,nz),
c     x     nf1(nx,ny,nz),
c     x     nfp1(nx,ny,nz),
c     x     uf(nx,ny,nz,3)
 
c      real minnf,maxnf
c      real flx(nx,ny,nz,3)
c      real ufc(nx,ny,nz,3)

cc      call face_to_center(uf,ufc)
 
c      minnf = 10.0e20
c      maxnf = 10.0e20

c      do 3 i=2,nx-1
c         do 3 j=2,ny-1
c            do 3 k=2,nz-1
c               nf1(i,j,k) = nfp1(i,j,k)
c 3             continue      

c      call periodic_scalar(nf1)

c      do 5 i=2,nx-1
c         do 5 j=2,ny-1
c            do 5 k=2,nz-1
c               flx(i,j,k,1) = 0.5*(nf1(i,j,k)+nf1(i+1,j,k))*uf(i,j,k,1)
c               flx(i,j,k,2) = 0.5*(nf1(i,j,k)+nf1(i,j+1,k))*uf(i,j,k,2)
c               flx(i,j,k,3) = 0.5*(nf1(i,j,k)+nf1(i,j,k+1))*uf(i,j,k,3)
c 5             continue

c      call periodic(flx)

c      do 10 i=2,nx-1
c         do 10 j=2,ny-1
c            do 10 k=2,nz-1
c               nfp1(i,j,k) = (1.0/6.0)*(nf1(i+1,j,k)+nf1(i-1,j,k)+
c     x                               nf1(i,j+1,k)+nf1(i,j-1,k)+
c     x                               nf1(i,j,k+1)+nf1(i,j,k-1))  
c     x         - 0.5*(dtsub/dx)*(flx(i+1,j,k,1)-flx(i-1,j,k,1))
c     x         - 0.5*(dtsub/dy)*(flx(i,j+1,k,2)-flx(i,j-1,k,2))
c     x         - (dtsub/(dz_grid(k)+dz_grid(k+1)))*
c     x           (flx(i,j,k+1,3)-flx(i,j,k-1,3))
c 10            continue

c      do 50 i=2,nx-1
c         do 50 j=2,ny-1
c            do 50 k=2,nz-1
c               nf(i,j,k) = 0.5*(nf1(i,j,k) + nfp1(i,j,k))
c 50            continue

c      call periodic_scalar(nf)

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE trans_nf_LaxWend1(nf,nf1,nfp1,uf)
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real nf(nx,ny,nz),
c     x     nf1(nx,ny,nz),
c     x     nfp1(nx,ny,nz),
c     x     uf(nx,ny,nz,3)
 
c      real minnf,maxnf
c      real flx(nx,ny,nz,3)

c      do 3 k=2,nz-1
c         do 3 j=2,ny-1
c            do 3 i=2,nx-1
c               nf1(i,j,k) = nfp1(i,j,k)
c 3             continue   
   
c      call periodic_scalar(nf1)

c      do 5 k=2,nz-1
c         do 5 j=2,ny-1
c            do 5 i=2,nx-1
c               flx(i,j,k,1) = 0.5*(nf1(i,j,k)+nf1(i+1,j,k))*uf(i,j,k,1)
c               flx(i,j,k,2) = 0.5*(nf1(i,j,k)+nf1(i,j+1,k))*uf(i,j,k,2)
c               flx(i,j,k,3) = 0.5*(nf1(i,j,k)+nf1(i,j,k+1))*uf(i,j,k,3)
c 5             continue

c      call periodic(flx)
c      call periodic_scalar(nf1)

c      do 10 k=2,nz-1
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1
c               nf(i,j,k) = (1.0/12.0)*(nf1(i+1,j,k)+nf1(i-1,j,k)+
c     x                     nf1(i,j+1,k)+nf1(i,j-1,k)+
c     x                     nf1(i,j,k+1)+nf1(i,j,k-1)+6.0*nf1(i,j,k))  
c     x         - 0.5*(dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
c     x         - 0.5*(dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
c     x         - 0.5*(dtsub/dz_grid(k))*
c     x           (flx(i,j,k,3)-flx(i,j,k-1,3))
c 10            continue

c      call periodic_scalar(nf)
c      nf(nx-1:nx,:,:) = nf_init*0.01

c      return
c      end
cc----------------------------------------------------------------------



cc----------------------------------------------------------------------
c      SUBROUTINE trans_nf_LaxWend2(nf,nf1,nfp1,ufp1)
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real nf(nx,ny,nz),
c     x     nf1(nx,ny,nz),
c     x     nfp1(nx,ny,nz),
c     x     ufp1(nx,ny,nz,3)
 
c      real minnf,maxnf
c      real flx(nx,ny,nz,3)

c      call periodic_scalar(nf)

c      do 5 k=2,nz-1
c         do 5 j=2,ny-1
c            do 5 i=2,nx-1
c               flx(i,j,k,1) = 0.5*(nf(i,j,k)+nf(i+1,j,k))*ufp1(i,j,k,1)
c               flx(i,j,k,2) = 0.5*(nf(i,j,k)+nf(i,j+1,k))*ufp1(i,j,k,2)
c               flx(i,j,k,3) = 0.5*(nf(i,j,k)+nf(i,j,k+1))*ufp1(i,j,k,3)
c 5             continue

c      call periodic(flx)
c      call periodic_scalar(nf1)

c      do 10 k=2,nz-1
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1
c               nfp1(i,j,k) = nf1(i,j,k) +
c     x                     0.002*(nf1(i+1,j,k)+nf1(i-1,j,k)+
c     x                     nf1(i,j+1,k)+nf1(i,j-1,k)+
c     x                     nf1(i,j,k+1)+nf1(i,j,k-1)-6.0*nf1(i,j,k))
c     x         - (dtsub/dx)*(flx(i,j,k,1)-flx(i-1,j,k,1))
c     x         - (dtsub/dy)*(flx(i,j,k,2)-flx(i,j-1,k,2))
c     x         - (dtsub/dz_grid(k))*(flx(i,j,k,3)-flx(i,j,k-1,3))
c 10            continue

c      call periodic_scalar(nfp1)
c      nfp1(nx-1:nx,:,:) = nf_init*0.01

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE trans_pf_LaxWend1(pf,pf1,uf)
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real pf(nx,ny,nz),
c     x     pf1(nx,ny,nz),
c     x     uf(nx,ny,nz,3)
 
c      real minnf,maxnf
c      real t1(nx,ny,nz,3)
c      real t1c(nx,ny,nz)
c      real t2(nx,ny,nz)

c      parameter(gamma = 5./3.)

c      call periodic_scalar(pf1)

c      do 5 k=2,nz-1
c         do 5 j=2,ny-1
c            do 5 i=2,nx-1
c               t1(i,j,k,1) = uf(i,j,k,1)*(pf1(i+1,j,k)-pf1(i,j,k))/dx 
c               t1(i,j,k,2) = uf(i,j,k,2)*(pf1(i,j+1,k)-pf1(i,j,k))/dy
c               t1(i,j,k,3) = 
c     x               uf(i,j,k,3)*(pf1(i,j,k+1)-pf1(i,j,k))/dz_grid(k)
c               t2(i,j,k) = 
c     x            gamma*pf1(i,j,k)*((uf(i,j,k,1)-uf(i-1,j,k,1))/dx +
c     x                   (uf(i,j,k,2)-uf(i,j-1,k,2))/dy +
c     x                   (uf(i,j,k,3)-uf(i,j,k-1,3))/dz_cell(k))
c 5             continue

c      call periodic(t1)

c      do k = 2,nz-1
c         do j = 2,ny-1
c            do i = 2,nx-1
c               t1c(i,j,k) = (1./2.)*(t1(i,j,k,1)+t1(i-1,j,k,1)) +
c     x                      (1./2.)*(t1(i,j,k,2)+t1(i,j-1,k,2)) +
c     x                      (1./2.)*(t1(i,j,k,3)+t1(i,j,k-1,3))
c            enddo
c         enddo
c      enddo

cc      call periodic_scalar(t1c)
cc      call periodic_scalar(t2)

c      do 10 k=2,nz-1
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1
c               pf(i,j,k) = (1.0/12.0)*(pf1(i+1,j,k)+pf1(i-1,j,k)+
c     x                     pf1(i,j+1,k)+pf1(i,j-1,k)+
c     x                     pf1(i,j,k+1)+pf1(i,j,k-1)+6.0*pf1(i,j,k))  
c     x         - 0.5*dtsub*(t1c(i,j,k)+t2(i,j,k))
c               if (pf(i,j,k) .lt. 0.0) then
c                 ! write(*,*) 'error...',t1c(i,j,k),t2(i,j,k)
c               endif

c 10            continue

cc      pf = abs(pf)

c      call periodic_scalar(pf)
c      pf(nx-1:nx,:,:) = nf_init*0.01*kboltz*tempf0

c      return
c      end
cc----------------------------------------------------------------------


cc----------------------------------------------------------------------
c      SUBROUTINE trans_pf_LaxWend2(pf,pf1,ufp1)
cc----------------------------------------------------------------------
cCVD$R VECTOR
c      include 'incurv.h'

c      real pf(nx,ny,nz),
c     x     pf1(nx,ny,nz),
c     x     ufp1(nx,ny,nz,3)
 
c      real minnf,maxnf
c      real t1(nx,ny,nz,3)
c      real t1c(nx,ny,nz)
c      real t2(nx,ny,nz)
c      real pfp1(nx,ny,nz)

c      parameter(gamma = 5./3.)

cc      call periodic_scalar(pf1)

c      do 5 k=2,nz-1
c         do 5 j=2,ny-1
c            do 5 i=2,nx-1
c               t1(i,j,k,1) = ufp1(i,j,k,1)*(pf(i+1,j,k)-pf(i,j,k))/dx 
c               t1(i,j,k,2) = ufp1(i,j,k,2)*(pf(i,j+1,k)-pf(i,j,k))/dy
c               t1(i,j,k,3) = 
c     x               ufp1(i,j,k,3)*(pf(i,j,k+1)-pf(i,j,k))/dz_grid(k)
c               t2(i,j,k) = 
c     x           gamma*pf(i,j,k)*((ufp1(i,j,k,1)-ufp1(i-1,j,k,1))/dx +
c     x                   (ufp1(i,j,k,2)-ufp1(i,j-1,k,2))/dy +
c     x                   (ufp1(i,j,k,3)-ufp1(i,j,k-1,3))/dz_cell(k))
c 5             continue

c      call periodic(t1)

c      do k = 2,nz-1
c         do j = 2,ny-1
c            do i = 2,nx-1
c               t1c(i,j,k) = (1./2.)*(t1(i,j,k,1)+t1(i-1,j,k,1)) +
c     x                      (1./2.)*(t1(i,j,k,2)+t1(i,j-1,k,2)) +
c     x                      (1./2.)*(t1(i,j,k,3)+t1(i,j,k-1,3))
c            enddo
c         enddo
c      enddo

cc      call periodic_scalar(t1c)
cc      call periodic_scalar(t2)

c      do 10 k=2,nz-1
c         do 10 j=2,ny-1
c            do 10 i=2,nx-1
c               pfp1(i,j,k) =  
c     x                     0.002*(pf1(i+1,j,k)+pf1(i-1,j,k)+
c     x                     pf1(i,j+1,k)+pf1(i,j-1,k)+
c     x                     pf1(i,j,k+1)+pf1(i,j,k-1)-6.0*pf1(i,j,k))
c     x         +   pf1(i,j,k)
c     x         - dtsub*(t1c(i,j,k)+t2(i,j,k))
c 10            continue

c      pf1 = pfp1
cc      pf1 = abs(pf1)

c      call periodic_scalar(pf1)
c      pf1(nx-1:nx,:,:) = nf_init*0.01*kboltz*tempf0

c      return
c      end
cc----------------------------------------------------------------------














