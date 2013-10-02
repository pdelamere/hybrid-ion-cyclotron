;----------------------------------------------------------------------
;main program
;----------------------------------------------------------------------


dir = '/Volumes/MacD97-2/hybrid/SWAP/run_test/'

read_para,dir

restore,filename=dir+'para.sav'

read_coords,dir,x,y,z

xparr = fltarr(nt,100,3,/nozero)

file = dir+'c.test_part.dat'
print,' reading...',file
openr,1,file,/f77_unformatted

xp = fltarr(100,3)

readu,1,xp
print,xp
frmcnt = 1
xparr(0,*,*) = xp
i = 1

while not(eof(1)) do begin

   readu,1,xp
   xparr(i,*,*) = xp
   i = i+1

endwhile

close,1

w = window(dimensions=[1000,1000])

x = findgen(nx)
z = findgen(nz)

p = plot(/test,/nodata,xrange=[min(x),max(x)],zrange=[min(z),max(z)],$
        /current,aspect_ratio=1)

;sz=size(xp)
for i = 0,9 do begin
   p = plot(xparr(*,i,0)/dx,xparr(*,i,2)/dx,/current,'.',/overplot,$
            xrange=[min(x),max(x)],yrange=[min(z),max(z)],$
           aspect_ratio=1)
endfor



end
;----------------------------------------------------------------------
