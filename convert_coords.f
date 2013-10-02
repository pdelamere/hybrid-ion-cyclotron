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
      program main
c----------------------------------------------------------------------
      include 'incurv.h'

      call convert_coord()
      

      stop
      end
c----------------------------------------------------------------------

