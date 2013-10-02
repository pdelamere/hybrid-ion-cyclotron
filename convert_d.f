c----------------------------------------------------------------------
c converts 3d scalar variables from cray to IDL f77_unformatted
c compatable arrays
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_np()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real np(nx,ny,nz)

      open(10,file='c.npall.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='npall.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'np frame #....',m
         write(20) m
         read(10) np
         write(20) np
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_ne()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real ne(nx,ny,nz)

      open(10,file='c.neall.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='neall.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'ne frame #....',m
         write(20) m
         read(10) ne
         write(20) ne
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_nn()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real nn(nx,ny,nz)

      open(10,file='c.nnall.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='nnall.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'nn frame #....',m
         write(20) m
         read(10) nn
         write(20) nn
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_nf()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real nf(nx,ny,nz)

      open(10,file='c.nfall.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='nfall.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'nf frame #....',m
         write(20) m
         read(10) nf
         write(20) nf
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_pf()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real pf(nx,ny,nz)

      open(10,file='c.pf.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='pf.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'pf frame #....',m
         write(20) m
         read(10) pf
         write(20) pf
 30      continue

      return
      end
c----------------------------------------------------------------------




c----------------------------------------------------------------------
      subroutine convert_b1()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real arr(nx,ny,nz,3)

      open(10,file='c.b1all.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='b1all.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'b1 frame #....',m
         write(20) m
         read(10) arr
         write(20) arr
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_uf()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real arr(nx,ny,nz,3)

      open(10,file='c.ufall.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='ufall.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'uf frame #....',m
         write(20) m
         read(10) arr
         write(20) arr
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_up()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real arr(nx,ny,nz,3)

      open(10,file='c.up.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='upall.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'up frame #....',m
         write(20) m
         read(10) arr
         write(20) arr
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_aj()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real arr(nx,ny,nz,3)

      open(10,file='c.ajall.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='ajall.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'aj frame #....',m
         write(20) m
         read(10) arr
         write(20) arr
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_E()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real arr(nx,ny,nz,3)

      open(10,file='c.Eall.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='Eall.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'E frame #....',m
         write(20) m
         read(10) arr
         write(20) arr
 30      continue

      return
      end
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      subroutine convert_energy()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real in_E
      real in_EeP
      real Evp
      real Euf
      real EB1
      real EB1x
      real EB1y
      real EB1z
      real EE
      real EeP
      real in_chex
      real in_bill

      open(10,file='c.energy.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='energy.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout

      do 30 i=1,(nt-mbegin)
         read(10) m
         write(*,*) 'energy time step #....',m
         write(20) m
         read(10) in_E,in_EeP,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,
     x            in_chex,in_bill
         write(20) in_E,in_EeP,Evp,Euf,EB1,EB1x,EB1y,EB1z,EE,EeP,
     x             in_chex,in_bill
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_mom()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real pup(3)
      real puf(3)
      real peb(3)
      real input_p(3)

      open(10,file='c.momentum.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='momentum.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout

      do 30 i=1,(nt-mbegin)
         read(10) m
         write(*,*) 'momentum time step #....',m
         write(20) m
         read(10) pup, puf, peb, input_p
         write(*,*) pup, puf, peb, input_p
         write(20) pup, puf, peb, input_p
 30      continue

      return
      end
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      subroutine convert_p_cons()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      integer*4 n
      real pup(3)
      real puf(3)
      real input_p(3)

      open(10,file='c.p_conserve.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='p_conserve.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout

      do 30 i=1,(nt-mbegin)*10
         read(10) m
         read(10) n
c         write(*,*) 'p_conserve time step #....',m,n
         write(20) m
         write(20) n
         read(10) pup, puf, input_p
         write(*,*) pup
         write(*,*) puf
         write(*,*) input_p
         write(20) pup, puf, input_p
 30      continue

      return
      end
c----------------------------------------------------------------------




c----------------------------------------------------------------------
      subroutine convert_chex()
c----------------------------------------------------------------------
      include 'incurv.h'

      real chex

      open(10,file='c.chex.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='chex.dat',status='unknown',
     x        form='unformatted')

      do 30 i=1,(nt-mbegin)
         write(*,*) 'chex.....i',i
         read(10) chex
         write(20) chex
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_bill()
c----------------------------------------------------------------------
      include 'incurv.h'

      real bill

      open(10,file='c.bill.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='bill.dat',status='unknown',
     x        form='unformatted')

      do 30 i=1,(nt-mbegin)
         write(*,*) 'bill....i',i
         read(10) bill
         write(20) bill
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_coord()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer nnx,nny,nnz
      real qqx(nx),qqy(ny),qqz(nz),ddz_grid(nz),ddz_cell(nz)

      open(30,file='c.coord.dat',
     x        form='unformatted',status='unknown')
     
      open(40,file='coord.dat',status='unknown',
     x        form='unformatted')

      write(*,*) 'Converting coords.dat....'

      read(30) nnx
      write(*,*) 'nnx...',nnx
      write(40) nnx 

      read(30) nny
      write(40) nny

      read(30) nnz
      write(40) nnz

      write(*,*) nnx,nny,nnz

      read(30) qqx
      write(40) qqx

      read(30) qqy
      write(40) qqy

      read(30) qqz
      write(40) qqz

      read(30) ddz_grid
      write(40) ddz_grid

      read(30) ddz_cell
      write(40) ddz_cell

      close(30)
      close(40)

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_ugradu()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real arr(nx,ny,nz,3)

      open(10,file='c.ugradu.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='ugraduall.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'ugradu frame #....',m
         write(20) m
         read(10) arr
         write(20) arr
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_Epar()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real arr(nx,ny,nz)
      real arr1(nx,ny,nz,3)

      open(10,file='c.Eparall.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='Eparall.dat',status='unknown',
     x        form='unformatted')

      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'Epar frame #....',m
         write(20) m
         read(10) arr
c         do 40 ii = 1,nx
c            do 40 jj = 1,ny
c               do 40 kk = 1,nz
c                  arr1(ii,jj,kk,3) = arr(ii,jj,kk)
c 40      continue
         write(20) arr
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_ddj()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real ddj(nx,ny,nz)

      open(10,file='c.ddj.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='ddj.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'ddj frame #....',m
         write(20) m
         read(10) ddj
         write(20) ddj
 30      continue

      return
      end
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      subroutine convert_t1()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real t1(nx,ny,nz)

      open(10,file='c.t1.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='t1.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 't1 frame #....',m
         write(20) m
         read(10) t1
         write(20) t1
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_t2()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real t2(nx,ny,nz)

      open(10,file='c.t2.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='t2.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 't2 frame #....',m
         write(20) m
         read(10) t2
         write(20) t2
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_t3()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real t3(nx,ny,nz)

      open(10,file='c.t3.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='t3.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 't3 frame #....',m
         write(20) m
         read(10) t3
         write(20) t3
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_ve()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real ve(Ni_max,3)

      open(10,file='c.ve.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='ve.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 've frame #....',m
         write(20) m
         read(10) ve
         write(20) ve
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_xe()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real xe(Ni_max,3)

      open(10,file='c.xe.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='xe.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) nx
      write(20) ny
      write(20) nz

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'xe frame #....',m
         write(20) m
         read(10) xe
         write(20) xe
 30      continue

      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_vp()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real vp(Ni_max,3)

      open(10,file='c.vp.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='vp.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) Ni_max

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'vp frame #....',m
         write(20) m
         read(10) vp
         write(20) vp
 30      continue

      return
      end
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      subroutine convert_xp()
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real xp(Ni_max,3)

      open(10,file='c.xp.dat',
     x        form='unformatted',status='unknown')
     
      open(20,file='xp.dat',status='unknown',
     x        form='unformatted')

c      read(10) nnt,nnout,nnx,nny,nnz
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) Ni_max

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'xp frame #....',m
         write(20) m
         read(10) xp
         write(20) xp
 30      continue

      return
      end
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      program main
c----------------------------------------------------------------------

      write(*,*) 'starting convert...'

c      call convert_Epar()
      call convert_np()
c      call convert_ne()
c      call convert_nn()
      call convert_nf()
      call convert_uf()
      call convert_ugradu()
      call convert_b1()
      call convert_aj()
      call convert_E()
c      call convert_energy()
c      call convert_mom()
      call convert_up()
c      call convert_chex()
c      call convert_bill
      call convert_pf()
      call convert_coord()
c      call convert_p_cons()
c      call convert_ddj()
c      call convert_t1()
c      call convert_t2()
c      call convert_t3()
      call convert_vp()
      call convert_xp()
      

      stop
      end
c----------------------------------------------------------------------






