;------------------------------------------------
;This script is used to calculate 3 jet indices like jet speed (zonal wind max), jet latitude (latitude of U max) and jet zonality. 
;references: Riboldi et al. (2023) GRL , Galfi and Messori (2023) ERL
;250 hPa zonal wind is used.
;Information
;author: cui xiangyang
;time: 2023-10-21
;release: v1
;-------------------------------------------------
;get variables from shell scripts

datadir = getenv("DATADIR")
outdir  = getenv("outdir")
var     = getenv("var")
;-------------------------------------------------

f = addfile(datadir,"r")
 
u = f->u(:,0,:,:)
time = u&time
lat  = u&lat
lon  = u&lon
;------------------------------------------------
;calculate jet speed and jet latitude. 
jets = new((/dimsizes(time)/),float)
jets = 0.
jets!0 = "time"
jets&time = time
jets@units = "m/s"
jets@long_name = "Jet Speed(Umax)"

jetl = jets
jetl@units = ""
jetl@long_name = "Jet latitude(Umax)"

Ut = u(:,{15:75},{-60:0})
lattmp = Ut&lat

do i = 0, dimsizes(time)-1
    umax = dim_avg_n(Ut(i,:,:),1)
    jets(i) = max(umax)
    mm  = maxind(umax)
    l1  = lattmp(mm)
    jetl(i) = l1
    delete(umax)
end do
delete(lattmp)
;-------------------------------------------------
;calculate jet zonality.
jetz = jets
jetz@units = ""
jetz@long_name = "Jet zonality"

us = u(:,{10:80},{-75:50})
lattmp = us&lat
lontmp = us&lon
do i = 0, dimsizes(time)-1 
    j = 0
    phi_delta = abs( lattmp( maxind(us(i,:,j+1)) ) - lattmp( maxind(us(i,:,j)) ) )
    do j  = 1, dimsizes(lontmp)-1
        m1 = dim_maxind(us(i,:,j), 0)
        m2 = dim_maxind(us(i,:,j-1), 0)
	    phi_deltatmp = abs( lattmp(m1) - lattmp(m2) )
	    if (phi_deltatmp.gt.phi_delta) then
	        phi_delta = phi_deltatmp
	    end if
        delete(phi_deltatmp)
    end do
    jetz(i) = phi_delta
    delete(phi_delta)
end do
;------------------------------------------------

system("/bin/rm -f "+outdir)
ncdf = addfile(outdir,"c")
ncdf->jets = jets
ncdf->jetl = jetl
ncdf->jetz = jetz

