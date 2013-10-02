pro b_time_series,nf


  xsz = 1000.
  ysz = 500.
  file = 'ion_cyclotron.mp4'
  width = xsz
  height = ysz
  frames = 180
  fps = 30
  speed = 2
  samplerate = 22050L
 
  ; Create object and initialize video/audio streams
  oVid = IDLffVideoWrite(file)
  vidStream = oVid.AddVideoStream(width, height, fps)
 ; audStream = oVid.AddAudioStream(samplerate)



;dir = '/Volumes/Macd97-2/hybrid/SWAP/run_test/'
dir = './tmp/'

nframe=nf
read_para,dir
restore,filename=dir+'para.sav'
read_coords,dir,x,y,z
@get_const

dt = 2.0*dt

;xsz = 1000.
;ysz = 1000.
;window,xsize=xsz,ysize=ysz
XINTERANIMATE, SET=[xsz,ysz, nframe], /SHOWLOAD 

by = 0.0
tm = 0.0
!p.multi=[0,1,2]
for nfrm = 500,nframe do begin

   c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1

   by = [by,b1(1,1,nz/2+10,1)]
   tm = [tm,dt*nfrm]
   

;   plot,z,b1(1,1,*,0),yrange=[-3,3]
   plot,z,(b1(1,1,*,1)*16*mp/q)/1e-9,yrange=[-200,200]
   w = window(dimensions=[xsz,ysz],/buffer)   
   p = plot((z(nz/2)-z)/1800.,(b1(1,1,*,1)*16*mp/q)/1e-9,$
       /current,/buffer,yrange=[-100,100],xrange=[-2,2],$
       xtitle='z ($R_{Io}$)',ytitle='$B_y$ (nT)',font_size=18,thick=2,$
           xstyle=1,xthick=2,ythick=2)
   time = oVid.Put(vidStream, w.CopyWindow())
   w.close

   plot,tm,(by*16*mp/q)/1e-9,xtitle='time (s)',xrange=[0,dt*nframe]

   img = tvrd(0,0,xsz,ysz)

;   tvscl,img
   xinteranimate, frame = nfrm-1, image = img

;   plot,tm,by,xtitle='time (s)'

endfor

oVid.Cleanup

by = by(0:*)
nstep = n_elements(by)

;plot,tm,by
;stop

Ni = nstep
Ti = dt
f21 = Ni/2 + 1
f = indgen(Ni)
f(f21) = f21 - Ni + findgen(f21 - 2)
f = f/(Ni*Ti)

ftarr = fft(by)
ftarr = abs(ftarr)^2

wp = shift(abs(ftarr),-f21)
f = shift(f(0:*),-f21)
window,1
!p.multi=[0,1,1]
plot,/xlog,/ylog,f,wp,xrange=[1./(nstep*dt),max(f)/2],charsize=1.0,/xsty,/ysty

omegai = (q*1700e-9/(64*mp))/(2*!pi)
;plot,/ylog,f,wp,xrange=[0,2],/xsty,/ysty
plots,[omegai,omegai],[0.00001,1.0],linestyle=1


;plot,/xlog,/ylog,f,wp,charsize=1.0,/xsty,/ysty

xinteranimate,/keep_pixmaps

w = window(dimensions=[800,600])
fnsz=24
p = plot(tm*omegai,(by*16*mp/q)/1e-9,xtitle='$\Omega_{SO_2^+} t$',$
         ytitle='$B_y$ (nT)',xrange=[0,dt*nframe]*omegai,font_size=fnsz,thick=2,$
         /current)

w1 = window(dimensions=[800,600])
p1 = plot(f/omegai,wp,/ylog,xrange=[0.1,4],$
          xtitle='Frequency ($\Omega_{SO_2^+}$)',$
          ytitle='Wave Power',font_size=fnsz,thick=2,$
          /current,'o',sym_filled=1)
p1.yrange=[1e-9,1e-2]
p2 = plot([1,1],p1.yrange,'--',/overplot,thick=2)

return
end
