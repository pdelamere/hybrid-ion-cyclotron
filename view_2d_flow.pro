;dir = '/Volumes/Macd97-2/hybrid/Io/local/run_1_40/'
;dir = '/Volumes/Macd97-2/hybrid/Io/local/'
dir = '/Volumes/Macd97-2/hybrid/SWAP/run_test/'
;dir = '/Volumes/Macd97-2/hybrid/3d_buf/janus/run_2/'
;dir = './'
read_para,dir

restore,filename=dir+'para.sav'

read_coords,dir,x,y,z

nframe=100
procnum = 5
zslc = 1
buffstat = 1

xrng=90
yrng=90
ftsz=14   
tk = 2
cbt =0.02
cbs =1.3
sns = 5
ssz = 0.5

xsz = 1500.
ysz = 1000.
;window,1,xsize=xsz,ysize=ysz
XINTERANIMATE, SET=[xsz,ysz, nframe], /SHOWLOAD 

for nfrm = 1,nframe do begin

;   c_read_2dxy_vec_m_32,dir,'c.up_'+strtrim(string(procnum),2),nfrm,up,upxz
   c_read_3d_m_32,dir,'c.np',nfrm,np
   c_read_3d_vec_m_32,dir,'c.up',nfrm,up
   c_read_3d_vec_m_32,dir,'c.b1',nfrm,b1
;   c_read_2dxy_m_32,dir,'c.temp_p_'+strtrim(string(procnum),2),nfrm,tp,tpxz
;   c_read_2dxy_vec_m_32,dir,'c.b1_'+strtrim(string(procnum),2),nfrm,b1,b1xz
   
   dx = x(1)-x(0)
   dy = y(1)-y(0)
   dt = 1
   
   minx = 0
   maxx = max(x)/dx -1
   
   minz = 0
   maxz = nz-1                  ;max(z)/dx
   
   Rio = 400/dx
   
;   up2d = reform(up(*,*,*))
;   np2d = reform(np(*,*))
;   tp2d = reform(tp(*,*))
;   b12dxy = reform(b1(*,*,*))
;   b12d = reform(b1(*,ny/2,*,*))
   
;   ntot2d = 1.0*np2d;<10000e15
   
;   gamma = 5./3.
;   mp = 1.67e-27
;   ent = pf/(22.*mp*nf2d)^gamma
   
   
;   y = reverse(y)
;   x = reverse(x)
   
;   y = y(ny/2) - y
;   x = x(nx/2) - x
;   z = z(zslc) - z
   
   minx = 0
   maxx = max(x)/dx -1
   
   minz = 0
   maxz = max(y)/dx - 2
   
   
;   up2dx=rotate(up2d(*,*,0),2)
;   up2dy=rotate(up2d(*,*,1),2)
   
;   b12dx=rotate(b12d(*,*,0),2)
   
;   np2d = rotate(np2d(*,*),2)
;   ntot2d = rotate(ntot2d(*,*),2)
;   tp2d = rotate(tp2d(*,*),2)
;;   ent = rotate(ent,2)
   
;   utot = fltarr(nx,ny,3)
;   utot(*,*,0) = rotate(up2d(*,*,0),2)
;   utot(*,*,1) = rotate(up2d(*,*,1),2)
;   utot(*,*,2) = rotate(up2d(*,*,2),2)
   
   @get_const
   b1 = b1*(mproton)/q

   bt = sqrt(b1(*,*,*,0)^2 + b1(*,*,*,1)^2 + b1(*,*,*,2)^2)

   w = window(window_title='2d_Io',dimensions=[xsz,ysz],margin=0,$
              buffer=buffstat)

   im1 = image(reform(np(*,1,*))/1e15,x/dx,z/dx,axis_style=2,xtickdir=1,$
               ytickdir=1,rgb_table=33,/current,layout=[2,2,1],font_size=ftsz,$
               buffer=buffstat,dimensions=[xsz,ysz])

;   im1.refresh,/disable
;   im1.xrange=[-xrng,xrng]
;   im1.yrange=[-yrng,yrng]
   im1.xtitle='x ('+strmid(strtrim(string(dx),2),0,3)+' km)'
   im1.ytitle='y ('+strmid(strtrim(string(dx),2),0,3)+' km)'
   

;   whx=where(x/dx gt -xrng and x/dx lt xrng)
;   why=where(y/dx gt -yrng and y/dx lt yrng)
   
   s1 = streamline(reform(up(*,1,*,0)),reform(up(*,1,*,2)),x/dx,z/dx,$
                   x_streamparticles=21,y_streamparticles=21,$
                   color='white',overplot=1,streamline_stepsize=ssz,$
                   streamline_nsteps=sns,axis_style=2,arrow_size=1.0,$
                   auto_range=1,$
                   xtickdir=1,ytickdir=1,$
;                         xrange=[-xrng,xrng],$
;                         yrange=[-yrng,yrng],$
                   thick=tk,font_size=ftsz,arrow_thick=tk,buffer=buffstat)
   
;   cbs1 = colorbar(target=s1)

   cb1 = colorbar(target=im1,title='Density',orientation=1,textpos=1,$
                  font_size=ftsz)
   cb1.translate,-cbt,0,/normal
   cb1.scale,0.8,cbs


   im1 = image(reform(b1(*,1,*,0))/1e-9,x/dx,z/dx,axis_style=2,xtickdir=1,$
               ytickdir=1,rgb_table=33,/current,layout=[2,2,2],font_size=ftsz,$
               buffer=buffstat,dimensions=[xsz,ysz])

   cb1 = colorbar(target=im1,title='Magnetic field',orientation=1,textpos=1,$
                  font_size=ftsz)
   cb1.translate,-cbt,0,/normal
   cb1.scale,0.8,cbs


   im1 = image(reform(up(*,1,*,0)),x/dx,z/dx,axis_style=2,xtickdir=1,$
               ytickdir=1,rgb_table=33,/current,layout=[2,2,3],font_size=ftsz,$
               buffer=buffstat,dimensions=[xsz,ysz])

   cb1 = colorbar(target=im1,title='Flow',orientation=1,textpos=1,$
                  font_size=ftsz)
   cb1.translate,-cbt,0,/normal
   cb1.scale,0.8,cbs



;   im1.refresh


   im1.refresh

   img = im1.CopyWindow()
   
   tvscl,img,true=1
   xinteranimate, frame = nfrm-1, image = img

;   im1.delete
;   im.delete
;   s1.delete
;   cb.delete

;   w.close
   
endfor

xinteranimate,/keep_pixmaps


end
