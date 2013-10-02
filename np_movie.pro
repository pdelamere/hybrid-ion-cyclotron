close,1
rundir = '/Volumes/MacD97-2/hybrid/SWAP/run_pickup'
f_read_coord,rundir+'/coord.dat',x,y,z,dzc,dzg,nx,ny,nz

nfrm = 20

loadct,33
s_max = 255
s_min = 0
nnx = nx
nny = nz
stretch,s_min,s_max
device,decompose=0

zm = 4

img = fltarr(nnx,nny,nfrm)
 
; Initialize XINTERANIMATE: 
XINTERANIMATE, SET=[zm*nnx,zm*nny, nfrm], /SHOWLOAD 
 
lnvy = fltarr(nfrm)

for i = 1,nfrm do begin 
   f_read_3d_m_32,rundir+'/np_b_1',i,npb
   f_read_3d_m_32,rundir+'/npall_1',i,npt

   np = npt;+ npb
;   np = npt
   print,min(np)

;   f_read_3d_vec_m_32,'run2/b1all_2',i,b1

;   lnvy(i-1) = alog(max(abs(up(*,1,30:120,2))))
   h=bytscl(reform(np(*,1,*)))

;   h=rebin(reform(sqrt(b1(*,1,*,0)^2+b1(*,1,*,2)^2)),nnx,nny)
   img(0:nx-1,*,i-1) = h(*,*)
;   img(nx-1:2*nx-3,*,i-1) = h(1:*,*)
;   xinteranimate, frame = i-1, image = img<s_max
endfor

;img = bytscl(img)

img = bytscl(img)
;img(*,nz/2-96/2,*) = 255
;img(*,nz/2+96/2,*) = 255


for i = 0,nfrm-1 do begin
   xinteranimate, frame = i, image = rebin(img(*,*,i)<255,zm*nnx,zm*nny)
endfor

;plot,lnvy

xinteranimate,/keep_pixmaps



end
