c----------------------------------------------------------------------
      SUBROUTINE fgrd_reg()
c----------------------------------------------------------------------
      include 'incurv.h'

      rk=nz/2 + nz/4
      rj=ny/2
      ri=nx/2

      do 10 i=1,nx
         qx(i) = i*dx
 10      continue

      do 20 j=1,ny
         qy(j) = j*dy
 20      continue

      do 30 k=1,nz
         dz_grid(k) = delz
         dz_cell(k) = delz
         qz(k) = k*dz_grid(k)
 30      continue

      open(40,file='coord.dat',status='unknown',form='unformatted')
      write(40) nx
      write(40) ny
      write(40) nz
      write(40) qx
      write(40) qy
      write(40) qz
      write(40) dz_grid
      write(40) dz_cell
      close(40)

      

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE fgrd_irr()
c----------------------------------------------------------------------
      include 'incurv.h'

      rk=nz/2 + nz/4
      rj=ny/2
      ri=nx/2

      do 10 i=1,nx
         qx(i) = i*dx
 10      continue

      do 20 j=1,ny
         qy(j) = j*dy
 20      continue

      do 30 k=1,nz
c         dz_grid(k) =  delz+0.01*((nz/2.0)**2 - (k-nz/2.0)**2)
         dz_grid(k) = delz + nz/2 - abs(k-nz/2)
         qz(k) = qz(k-1)+dz_grid(k)
 30      continue

      dz_cell(1) = dz_grid(1)
      dz_cell(nz) = dz_grid(nz)
      do 40 k=2,nz-1 
         dz_cell(k) = ((qz(k+1) + qz(k))/2.0) - 
     x                ((qz(k) + qz(k-1))/2.0) 
 40      continue

      open(40,file='coord.dat',status='unknown',form='unformatted')
      write(40) nx
      write(40) ny
      write(40) nz
      write(40) qx
      write(40) qy
      write(40) qz
      write(40) dz_grid
      write(40) dz_cell
      close(40)

      

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE fsetup(nf,b0,bt)
c----------------------------------------------------------------------
      include 'incurv.h'

      real eoverm
      parameter(eoverm = q/mO)
 

      real nf(nx,ny,nz),
     x     b0(nz),
     x     bt(nx,ny,nz,3)

      do 10 k=1,nz
         b0(k) = 3e-5*eoverm   !/sec
10       continue

      do 20 i=1,nx
         do 20 j=1,ny
            do 20 k=1,nz
               nf(i,j,k) = 10e20 
               bt(i,j,k,3) = b0(k)                  
 20            continue

      open(30,file='nf.dat',status='unknown',form='unformatted')
      write(30) nz
      write(30) nf
      close(30)

      open(40,file='b0.dat',status='unknown',form='unformatted')
      write(40) nz
      write(40) b0
      close(40)

      return
      end
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      SUBROUTINE f_init(b0,b1,b12,uf,uf1,nf,up,np,input_p)
c----------------------------------------------------------------------
      include 'incurv.h'

      real b0(nz)
      real b1(nx,ny,nz,3)
      real b12(nx,ny,nz,3)
      real uf(nx,ny,nz,3)
      real uf1(nx,ny,nz,3)
      real nf(nx,ny,nz)
      real up(nx,ny,nz,3)
      real np(nx,ny,nz)
      real input_p(3)

      real mO_q
      parameter(mO_q = mO/q)
      real vol
      real denf

      input_E = 0.0

      do 10 i=2,nx-1
         do 10 j=2,ny-1

c            b1(i,ny/2,k,3) = b0(k)/2.0 
c            b12(i,ny/2,k,3) = b0(k)/2.0
c            b1(i,(1+ny/2),k,3) = b0(k)/4.0 
c            b12(i,(1+ny/2),k,3) = b0(k)/4.0
c            b1(i,(ny/2-1),k,3) = b0(k)/4.0 
c            b12(i,(ny/2-1),k,3) = b0(k)/4.0


c            b1(i,j,nz/2,2) = 1.0
c            b12(i,j,nz/2,2) = 1.0
c            b1(i,j,(1+nz/2),2) = 0.5
c            b12(i,j,(1+nz/2),2) = 0.5
c            b1(i,j,((nz/2)-1),2) = 0.5
c            b12(i,j,((nz/2)-1),2) = 0.5

            
c            uf1(i,rj,k,2) = 1.0
c            uf(i,rj,k,2) = 1.0

c            uf1(i,(1+rj),k,2) = 0.95
c            uf(i,(1+rj),k,2) = 0.95
c            uf1(i,(rj-1),k,2) = 0.95
c            uf(i,(rj-1),k,2) = 0.95

c            uf1(i,(2+rj),k,2) = 0.85
c            uf(i,(2+rj),k,2) = 0.85
c            uf1(i,(rj-2),k,2) = 0.85
c            uf(i,(rj-2),k,2) = 0.85
            
c            uf1(i,(3+rj),k,2) = 0.7
c            uf(i,(3+rj),k,2) = 0.7
c            uf1(i,(rj-3),k,2) = 0.7
c            uf(i,(rj-3),k,2) = 0.7

c            uf1(i,(4+rj),k,2) = 0.5
c            uf(i,(4+rj),k,2) = 0.5
c            uf1(i,(rj-4),k,2) = 0.5
c            uf(i,(rj-4),k,2) = 0.5

c            uf1(i,(5+rj),k,2) = 0.25
c            uf(i,(5+rj),k,2) = 0.25
c            uf1(i,(rj-5),k,2) = 0.25
c            uf(i,(rj-5),k,2) = 0.25

c            uf1(i,(6+rj),k,2) = 0.1
c            uf(i,(6+rj),k,2) = 0.1
c            uf1(i,(rj-6),k,2) = 0.1
c            uf(i,(rj-6),k,2) = 0.1

            uf1(i,j,rk,2) = 1.0
            uf(i,j,rk,2) = 1.0

            uf1(i,j,(1+rk),2) = 0.98
            uf(i,j,(1+rk),2) = 0.98
            uf1(i,j,((rk)-1),2) = 0.98
            uf(i,j,((rk)-1),2) = 0.98

            uf1(i,j,(2+rk),2) = 0.85
            uf(i,j,(2+rk),2) = 0.85
            uf1(i,j,((rk)-2),2) = 0.85
            uf(i,j,((rk)-2),2) = 0.85

            uf1(i,j,(3+rk),2) = 0.70
            uf(i,j,(3+rk),2) = 0.70
            uf1(i,j,((rk)-3),2) = 0.70
            uf(i,j,((rk)-3),2) = 0.70

            uf1(i,j,(4+rk),2) = 0.5
            uf(i,j,(4+rk),2) = 0.5
            uf1(i,j,((rk)-4),2) = 0.5
            uf(i,j,((rk)-4),2) = 0.5

            uf1(i,j,(5+rk),2) = 0.25
            uf(i,j,(5+rk),2) = 0.25
            uf1(i,j,((rk)-5),2) = 0.25
            uf(i,j,((rk)-5),2) = 0.25

            uf1(i,j,(6+rk),2) = 0.1
            uf(i,j,(6+rk),2) = 0.1
            uf1(i,j,((rk)-6),2) = 0.1
            uf(i,j,((rk)-6),2) = 0.1

c            up(i,j,(6+rk),1) = 10.0
c            up(i,j,(rk-6),1) = 10.0
c            np(i,j,(6+rk)) = 10e22
c            np(i,j,(rk-6)) = 10e22

c            up(i,j,(5+rk),1) = 10.0
c            up(i,j,(rk-5),1) = 10.0
c            np(i,j,(5+rk)) = 10e22
c            np(i,j,(rk-5)) = 10e22

c            up(i,j,(4+rk),1) = 10.0
c            up(i,j,(rk-4),1) = 10.0
c            np(i,j,(4+rk)) = 10e22
c            np(i,j,(rk-4)) = 10e22

c            up(i,j,(3+rk),1) = 10.0
c            up(i,j,(rk-3),1) = 10.0
c            np(i,j,(3+rk)) = 10e22
c            np(i,j,(rk-3)) = 10e22

c            up(i,j,(2+rk),1) = 10.0
c            up(i,j,(rk-2),1) = 10.0
c            np(i,j,(2+rk)) = 10e22
c            np(i,j,(rk-2)) = 10e22

c            up(i,j,(1+rk),1) = 10.0
c            up(i,j,(rk-1),1) = 10.0
c            np(i,j,(1+rk)) = 10e22
c            np(i,j,(rk-1)) = 10e22

c            up(i,j,(rk),1) = 10.0
c            up(i,j,(rk),1) = 10.0
c            np(i,j,(rk)) = 10e22
c            np(i,j,(rk)) = 10e22


 10         continue

      do 20 i=2,nx-1
         do 20 j=2,ny-1
            do 20 k=2,nz-1
               do 20 m=1,3
                  vol = dx*dy*dz_cell(k)*(km_to_m)**3
                  denf = nf(i,j,k)/(km_to_m**3)
                  input_E = input_E + 
     x                        0.5*mO*vol*denf*(uf(i,j,k,m)*km_to_m)**2
     x                     +  (vol/(2.0*mu0))*(mO_q*b1(i,j,k,m))**2
                  input_p(m) = input_p(m) + mO*vol*denf*uf(i,j,k,m)
 20               continue

      write(*,*) 'Input E...',input_E
      write(*,*) 'Input P...',input_p(1),input_p(2),input_p(3)

      open(10, file='ufinit.dat', status='unknown', form='unformatted')
      write(10) nx
      write(10) ny
      write(10) nz
      write(10) uf
      close(10)


      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE f_diag(uf,nf,b1,E,up,np)
c----------------------------------------------------------------------
      include 'incurv.h'

      real uf(nx,ny,nz,3)
      real nf(nx,ny,nz)
      real b1(nx,ny,nz,3)
      real E(nx,ny,nz,3)
      real up(nx,ny,nz,3)
      real np(nx,ny,nz)

      real mO_q
      parameter(mO_q = mO/q)

      real Eu               !kinetic energy of fluid flow
      real Eup
      real EB               !Magnetic field energy density
      real EE               !Electric field energy density
      real vol              !volume of cell
      real total_E          !total energy
      real denf             !fluid density
      real denp

      Eu = 0.0
      Eup = 0.0
      EB = 0.0
      EE = 0.0
      do 10 i=2,nx-1
         do 10 j=2,ny-1
            do 10 k=2,nz-1
               do 10 m=1,3
                  vol = dx*dy*dz_cell(k)*km_to_m**3
                  denf = nf(i,j,k)/(km_to_m**3)
                  denp = np(i,j,k)/(km_to_m**3)
                  Eu = Eu + 0.5*mO*denf*vol*(uf(i,j,k,m)*km_to_m)**2
                  Eup = Eup+0.5*mBa*denp*vol*(up(i,j,k,m)*km_to_m)**2
                  EB = EB + (vol/(2.0*mu0))*(mO_q*b1(i,j,k,m))**2
                  EE = EE + (epsilon*vol/2.0)*(mO_q*E(i,j,k,m))**2
 10               continue

      write(*,*) 'Initial energy........',input_E
      write(*,*) 'Total up ke...........',Eup
      write(*,*) 'Total uf ke...........',Eu
      write(*,*) 'Total B energy........',EB
      write(*,*) 'Total E energy........',EE
      write(*,*) 'Total energy..........',Eup+Eu+EB+EE
      write(*,*) 'Normalized energy.....',(Eup+Eu+EB+EE)/input_E
      write(*,*) ' '

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      SUBROUTINE foutput(uf,b1,bt,aj)
c----------------------------------------------------------------------
      include 'incurv.h'

      real uf(nx,ny,nz,3)
      real b1(nx,ny,nz,3)
      real bt(nx,ny,nz,3)
      real aj(nx,ny,nz,3)

      open(10, file='uf.dat', status='unknown', form='unformatted')
      write(10) nx
      write(10) ny
      write(10) nz
      write(10) uf
      close(10)

      open(20, file='aj.dat', status='unknown', form='unformatted')
      write(20) nx
      write(20) ny
      write(20) nz
      write(20) aj
      close(20)

      open(30, file='b1.dat', status='unknown', form='unformatted')
      write(30) nx
      write(30) ny
      write(30) nz
      write(30) b1
      close(30)

      open(40, file='bt.dat', status='unknown', form='unformatted')
      write(40) nx
      write(40) ny
      write(40) nz
      write(40) bt
      close(40)


      return
      end
c----------------------------------------------------------------------














