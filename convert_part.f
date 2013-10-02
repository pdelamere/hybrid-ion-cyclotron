c----------------------------------------------------------------------
c converts 3d scalar variables from cray to IDL f77_unformatted
c compatable arrays
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine convert_vp(filenum)
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real vp(Ni_max,3)
      character(2) filenum

      open(10,file='c.vp_'//trim(filenum)//'.dat',
     x  status='unknown',form='unformatted',convert='little_endian')

       open(20,file='vp_'//trim(filenum)//'.dat',
     x      status='unknown',form='unformatted')

      write(*,*) 'Ni_max...',Ni_max
      write(*,*) nt,nout,nx,ny,nz,Ni_max
      write(20) nt
      write(20) nout
      write(20) Ni_max

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'vp frame #....',m,Ni_max
         write(20) m
         read(10) vp
         write(20) vp
 30      continue

      return
      end
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      subroutine convert_xp(filenum)
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real xp(Ni_max,3)
      character(2) filenum

      open(10,file='c.xp_'//trim(filenum)//'.dat',
     x  status='unknown',form='unformatted',convert='little_endian')

       open(20,file='xp_'//trim(filenum)//'.dat',
     x      status='unknown',form='unformatted')

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
      subroutine convert_mrat(filenum)
c----------------------------------------------------------------------
      include 'incurv.h'

      integer*4 m
      real mrat_arr(Ni_max)
      character(2) filenum

      open(10,file='c.mrat_'//trim(filenum)//'.dat',
     x     status='unknown',form='unformatted',convert='little_endian')

      open(20,file='mrat_'//trim(filenum)//'.dat',
     x     status='unknown',form='unformatted')
      
      write(*,*) nt,nout,nx,ny,nz
      write(20) nt
      write(20) nout
      write(20) Ni_max

      do 30 i=1,(nt-mbegin)/nout
         read(10) m
         write(*,*) 'mrat frame #....',m
         write(20) m
         read(10) mrat_arr
         write(20) mrat_arr
 30      continue

      return
      end
c----------------------------------------------------------------------



c----------------------------------------------------------------------
      program main
c----------------------------------------------------------------------

      character(2) filenum(16) !max 16 processors

      write(*,*) 'starting convert...'

      filenum = (/'0 ','1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ',
     x     '9 ','10','11','12','13','14','15'/)

      do i = 1,10

         write(*,*) 'converting file...',i

         call convert_vp(filenum(i))
         call convert_xp(filenum(i))
         call convert_mrat(filenum(i))
      
      enddo

      stop
      end
c----------------------------------------------------------------------






