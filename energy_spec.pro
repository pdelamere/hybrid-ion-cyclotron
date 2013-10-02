;----------------------------------------------------------------
PRO read_part,file,nfrm,Ni_max,xp
;----------------------------------------------------------------

;Ni_max=long(0)
;nt=0l
;ntout=0l
frm=0l

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
;readu,1,nt
;readu,1,ntout
;readu,1,Ni_max
;print,nt,ntout,Ni_max




xp=fltarr(Ni_max,3,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end

;----------------------------------------------------------------


;----------------------------------------------------------------
PRO read_part_scalar,file,nfrm,Ni_max,xp
;----------------------------------------------------------------

;Ni_max=long(0)
;nt=0
;ntout=0
frm=0

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
;readu,1,nt
;readu,1,ntout
;readu,1,Ni_max
;print,nt,ntout,Ni_max

xp=fltarr(Ni_max,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end
;----------------------------------------------------------------


;----------------------------------------------------------------
pro get_dist,xyz_traj,lxyz,nlxyz,t
;----------------------------------------------------------------
close,1

mp=1.67e-27
nfile = '1'
nfrm = 4
;define look direction of instrument
;theta = 180.*!dtor
;phi = !pi/4
;dtheta= 20.*!dtor
;dphi = 20.*!dtor

dir = '/Volumes/MacD97-2/hybrid/SWAP/run_test/'

read_para,dir

restore,filename=dir+'para.sav'

read_coords,dir,x,y,z

part_dir = dir+'part_1/'
fluid_fields_dir = dir

print,part_dir

c_read_3d_m_32,fluid_fields_dir,'c.np',nfrm,npt
c_read_3d_m_32,fluid_fields_dir,'c.np_b',nfrm,npb
c_read_3d_vec_m_32,fluid_fields_dir,'c.b1',nfrm,b1
c_read_3d_vec_m_32,fluid_fields_dir,'c.up',nfrm,up
c_read_3d_m_32,fluid_fields_dir,'c.temp_p',nfrm,ti

up2d = fltarr(nx,nz,3)
up2d(*,*,*) = reform(up(*,1,*,*))

np = npt+npb

xx = x

z1=nz/2
z2 = nz/2

x1 = 0
x2 = nx-1

zz = z(nz/2-z1:nz/2+z2)

xx  = x(x1:x2)

np_img = reform(np(x1:x2,1,nz/2-z1:nz/2+z2))
b1_img = reform(sqrt(b1(*,1,*,0)^2 + b1(*,1,*,1)^2 + b1(*,1,*,2)^2))

;window,0
;!p.multi=[0,1,1]
;device,decompose=0
;loadct,39
xx = indgen(nx)
zz = indgen(nz)
;xx = xg
;zz = z
;window,0
x1 = 0
x2 = nx-1
z1 = 0 ;nz/2-10
z2 = nz-1 ;nz/2+10

sz = size(b1_img)
im = image(b1_img(x1:x2,z1:z2),xx(x1:x2),zz(z1:z2),axis_style=2,rgb_table=33,$
          xtickdir=1,ytickdir=1)

scl = 500/sz(1)
im.scale,scl,scl,0

consx = 1
consz = consx
sz = size(up2d(x1:x2,z1:z2))

xxx = congrid(xx,sz(1)/consz)
zzz = congrid(zz(z1:z2),sz(2)/consz)

vx = congrid(reform(up2d(x1:x2,z1:z2,0)),sz(1)/consz,sz(2)/consz)
vz = congrid(reform(up2d(x1:x2,z1:z2,2)),sz(1)/consz,sz(2)/consz)

scl=1
vec = vector(vx(0:*:scl,0:*:scl),vz(0:*:scl,0:*:scl),xxx(0:*:scl),zzz(0:*:scl),$
             /overplot,color="white",thick=1,$
             head_size=0.5,x_subsample=5,y_subsample=5,length_scale=0.9,$
            xrange=[xx(x1),xx(x2)],yrange=[zz(z1),zz(z2)])   
   


xcur = nx/2 ;fix(xcur)
zcur = nz/2 ;fix(zcur)
ycur=1

print,xcur,ycur,zcur

ztemp = fltarr(35)

dv = 0.5
nn = 2000.
alpha=1.9263418e-20
d3v = dv^3
vxyp = fltarr(nn/dv,nn/dv)
vxzp = fltarr(nn/dv,nn/dv)
vyzp = fltarr(nn/dv,nn/dv)

vx = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
vy = -(nn/(2*dv))*dv + findgen(nn/dv)*dv
vz = -(nn/(2*dv))*dv + findgen(nn/dv)*dv

dx = 200.0
ndx = 5.0

b = b1(xcur,ycur,zcur,*)
bx = b1(xcur,ycur,zcur,0)
by = b1(xcur,ycur,zcur,1)
bz = b1(xcur,ycur,zcur,2)


bxy = b(1)/b(0)
bzx = b(2)/b(0)
bzy = b(2)/b(1)


!p.multi=[0,1,1]
;f_read_coord,fluid_fields_dir+'coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfil = 0
xfile = dir+'c.xp_'+strtrim(string(nfil),2)
vfile = dir+'c.vp_'+strtrim(string(nfil),2)
mratfile = dir+'c.mrat_'+strtrim(string(nfil),2)

read_part,xfile,nfrm,Ni_max,xp
read_part,vfile,nfrm,Ni_max,vp
read_part_scalar,mratfile,nfrm,Ni_max,mrat


for nfil = 1,2 do begin
    xfile = dir+'c.xp_'+strtrim(string(nfil),2)
    vfile = dir+'c.vp_'+strtrim(string(nfil),2)
    mratfile = dir+'c.mrat_'+strtrim(string(nfil),2)

    read_part,xfile,nfrm,Ni_max,xpp
    read_part,vfile,nfrm,Ni_max,vpp
    read_part_scalar,mratfile,nfrm,Ni_max,mratt

    xp = [xp,xpp]
    vp = [vp,vpp]
    mrat = [mrat,mratt]

endfor

;save,filename='part_arr.sav',xp,vp,mrat
;restore,'part_arr.sav'

i = xcur
j = ycur
k = zcur

dE = 20.
hmin = 1.0
hmax = 10000.

nh = fix((hmax-hmin)/dE)+1

zcnt=10
nnz = (nz/2+zcnt) - (nz/2 - zcnt)
evst = fltarr(nnz+1,nh)

;plots,[xx(i),xx(i)],[zz(nz/2-zcnt),zz(nz/2+zcnt)],/data,linestyle=0,thick=4

cnt = 0
zind = nz/2
;for zind = nz/2-zcnt,nz/2+zcnt do begin

   i = xcur
   j = 1
   k = zind
   print,zind,np_img(i,k)

;   wset,0
   plots,xx(i),zz(k),/data,psym=6

   wh = where((xp(*,0) ge x(i)-ndx*dx) and (xp(*,0) le x(i)+ndx*dx) and $
              (xp(*,1) ge y(j)-ndx*dx) and (xp(*,1) le y(j)+ndx*dx) and $ 
              (xp(*,2) ge z(k)-ndx*dx) and (xp(*,2) le z(k)+ndx*dx))
   print,n_elements(wh),z(k)

   if (wh(0) gt -1) then begin
      e_arr = 0
      cnt_arr = 0
      for l = 0ll,n_elements(wh)-1 do begin
         ii = fix(vp(wh(l),0)/dv) + (nn/(2*dv))
         jj = fix(vp(wh(l),1)/dv) + (nn/(2*dv))
         kk = fix(vp(wh(l),2)/dv) + (nn/(2*dv))
         vxyp(ii,jj) = 1.0 + vxyp(ii,jj)
         vxzp(ii,kk) = 1.0 + vxzp(ii,kk)
         vyzp(jj,kk) = 1.0 + vyzp(jj,kk)
         vpp = reform(vp(wh(l),*))
;      vpp2 = sqrt(vpp(0)^2 + vpp(1)^2 + vpp(2)^2)
;      vpp = vpp/vpp2
         vdotl = transpose(vpp(*))#[1,0,0]

         if (vdotl gt cos(120*!dtor)) then begin
            e_arr = [e_arr,(vpp(0)^2 + vpp(1)^2 + vpp(2)^2)/mrat(wh(l))]
            cnt_arr = [cnt_arr,vdotl]
         endif
         
      endfor
   endif

   sz = size(vxzp)
   im1 = image(smooth(alog(vxzp>1),2),vx,vz,rgb_table=33,$
               axis_style=2,xtickdir=1,ytickdir=1,$
              xtitle='$v_x$',ytitle='$v_z$')
;   scl = 500./sz(1)
   scl=1
   im1.scale,scl,scl,1


   help,e_arr
   e_arr = 0.5*1.67e-27*e_arr*1e6/1.6e-19   ;convert to eV

   h = histogram(e_arr,binsize=dE,min = hmin, max = hmax,reverse_indices=ri)
   h1 = fltarr(n_elements(h))
   for i = 0,n_elements(h)-2 do begin
      if (ri[i] ne ri[i+1]) then begin
         h1(i) = total(cnt_arr(ri(ri(i):ri(i+1)-1)))
      endif
   endfor

   h = h1

   xE = dE*(findgen(n_elements(h)) + hmin)

   plot,xE,h,psym=10,/xlog,xrange=[10,10000],/ylog,yrange=[1,max(h)]
;   stop

;rebin to log scale
   s = {energy_bin, e_min: 0.0, e_max: 0.0}
   close,2
   openr,2,'caps_e_bin.dat'

   levst = 0
   lxE = 10
   while not(eof(2)) do begin
      readf,2,s
      emin = s.e_min
      emax = s.e_max      
      if ((emin ge 1.0) and (emax le 10000)) then begin
         wh = where((xE gt emin) and (xE le emax))
                                ; if (wh(0) ge 0) then levst = [levst,total(h(wh))]
         if (wh(0) ge 0) then levst = [total(h(wh)),levst]
         if (wh(0) eq -1) then levst = [0,levst]
         lxE = [(emax+emin)/2,lxE]
      endif
   endwhile

   emin = 1
   emax = 1

   lxE = lxE(0:n_elements(lxE)-2)
   levst = levst(0:n_elements(levst)-2)
   plot,lxE,levst,psym=10,/xlog,xrange=[10,10000],/xsty
;   b = barplot(lxE,levst,/xlog,xrange=[10,10000],/current)
   ;b = barplot(lxE,levst,xrange=[10,10000],nbars=1,width=2.0)
   
 

   if (cnt eq 0) then evst = fltarr(nnz+1,n_elements(lxE))

;   evst(cnt,*) = alog(levst>1)
   evst(cnt,*) = levst


   cnt = cnt+1
;endfor

;c = contour(evst,findgen(nnz+1),lxE,rgb_table=33,/ylog)


;im = image(evst,rgb_table=33)
;im.scale,10,10,1
;yax = axis('y',location=[0,0],title='Energy (eV)')


pg_plotimage,evst,findgen(nnz+1),lxE,/ylog
pg_plotimage_new,evst,findgen(nnz+1),lxE,/ylog,ytitle='Energy (eV)'




;c = image(evst,findgen(nnz+1),lxE,/ylog,/fill,rgb_table=33)
;c = image(newimage,rgb_table=33,axis_style=2,/ylog,/xlog);,axis_style=2,xtickdir=1,ytickdir=1)
;c.scale,10,10,1
;yaxis = axis('y',location=[0,0],title='Energy (eV)',/log,tickdir=1)

;c.axis_style=2

return
end
;----------------------------------------------------------------------


;----------------------------------------------------------------------
pro read_traj,t
;----------------------------------------------------------------------

t = {trajectory,Julian: 0.0, x: 0.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, $
     vz: 0.0,sx: 0.0, sy: 0.0, sz: 0.0, px: 0.0, py: 0.0, pz: 0.0} 

d = {trajectory,Julian: 0.0, x: 0.0, y: 0.0, z: 0.0, vx: 0.0, vy: 0.0, $
     vz: 0.0,sx: 0.0, sy: 0.0, sz: 0.0, px: 0.0, py: 0.0, pz: 0.0} 

close,1
openr, 1, 'trajectory1.csv'

junk = ' ' 
readf, 1, junk
print,junk


while not(eof(1)) do begin
   readf,1,d
   t = [t,d]
endwhile

t = t(1:*)

close,1

return
end
;----------------------------------------------------------------------


;----------------------------------------------------------------------
;main program
;----------------------------------------------------------------------

;set_plot,'ps'
;;device,/landscape
;!p.font=0
;device,filename='vdist.ps'
;!p.thick=2.0
;!x.thick=2.0
;!y.thick=2.0
!p.multi=[0,2,3]
;!p.charsize=1.5
;@x6x9
;device,/color


;read_traj,t
;rundir = '/Volumes/MacD97-2/hybrid/SWAP/run_test/'
;f_read_coord,rundir+'coord.dat',x,y,z,dzc,dzg,nx,ny,nz



;xyz_traj=[nx/2,ny/2,nz/2]
;lxyz = [1.0,0.0,0.0]
;nlxyz = [0.0,1.0,0.0]

get_dist,xyz_traj,lxyz,nlxyz

;stop

;t2 = -13.5e3/cos(15*!dtor) - atan(15*!dtor)*xx
;dx = 2200.

;xx = x - x(55)
;yy = y - y(ny/2)

;nstep = 10
;for i = 0,nstep-1 do begin 
;   isc = (nstep-i)*10
;   jsc = ny/2 - round((13.5e3/dx)/cos(15*!dtor)) - atan(15*!dtor)*(isc-(nx-55))
;   print,isc,jsc
;   xyz_traj=[isc,jsc,nz/2-35]
;   get_dist,xyz_traj,lxyz,nlxyz,t
;endfor

;evice,/close

stop
end
;----------------------------------------------------------------------
