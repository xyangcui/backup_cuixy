;-----------------------------------------------------------------
;; draw general budget terms. 
;; MEJS box (10-70E,10-40N)
;; EAJS box (110-180E,15-50N)
;; time: 2023-12-13
;-----------------------------------------------------------------
rad   = 4.*atan(1.)/180. 
tday  = 24.*60.*60.
mlatn = 30
mlats = 20
mlonw = 30
mlone = 70

elatn = 40
elats = 20
elonw = 110
elone = 180

ddir = "/home/sunming/data5/cuixy/global_jets/data/temp_budget"
tdir = ddir+"/"+"anom_dtdt.nc"
udir = ddir+"/"+"anom_udtdx.nc"
vdir = ddir+"/"+"anom_vdtdy.nc"
wdir = ddir+"/"+"anom_wdtdp.nc"
rdir = ddir+"/"+"anom_wrtcp.nc"
pdirmejs = "/home/sunming/data5/cuixy/global_jets/plot/MEJS_budget"
pdireajs = "/home/sunming/data5/cuixy/global_jets/plot/EAJS_budget"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------

ft  = addfile(tdir,"r")
fu  = addfile(udir,"r")
fv  = addfile(vdir,"r")
fw  = addfile(wdir,"r")
fr  = addfile(rdir,"r")

t = ft->dtdt
t = t
u = fu->udtdx
u = u*tday
v = fv->vdtdy
v = v*tday
w = fw->wdtdp
w = w*tday
r1 = fr->wrtcp
r1 = r1*tday
;residual
r = t + (u + v + w - r1)
copy_VarCoords(t, r)
u = -1*u
v = -1*v
w = -1*w+r1
lat = t&latitude
;---------------------------------------------------------------=
;; MJO phases.
time = ft->time
TIME = cd_calendar(time, 0)
year = toint(TIME(:,0))
mon  = toint(TIME(:,1))
day  = toint(TIME(:,2))

rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
delim  = ","
rmmf   = asciiread(rmmdir,-1,"string")
iyear  = str_get_field(rmmf,1,delim)
imon   = str_get_field(rmmf,2,delim)
iday   = str_get_field(rmmf,3,delim)
phaset = str_get_field(rmmf,6,delim)
ampt   = str_get_field(rmmf,7,delim)

iyear  := toint(iyear(1::))
imon   := toint(imon(1::))
iday   := toint(iday(1::))
phaset := toint(phaset(1::))
ampt   := tofloat(ampt(1::))

amp   = new((/dimsizes(time)/),float)
phase = new((/dimsizes(time)/),integer)
amp   = 0.
phase = 0

do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    amp(i)   = ampt(nt)
    phase(i) = phaset(nt)
    delete(nt)
end do

nn = ispan(0, dimsizes(time)-1, 1)
amp!0   = "timeind"
phase!0 = "timeind"
amp&timeind   = nn
phase&timeind = nn
delete([/nn,ampt,phaset,iyear,imon,iday/])

print("complete RMM selling "+ systemfunc("date"))

n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

amp   := amp(n1:n2)
phase := phase(n1:n2)
mon   := mon(n1:n2)
time  := time(n1:n2)

nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
amp   := amp(nt)
phase := phase(nt)

delete(nt)
ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.

;; store lead-lag values for plotting. 
llmejs = new((/8,41,5/),float)
llmejs = 0.
tm = t(:,{mlats:mlatn},{mlonw:mlone})
;mejs
clatjs = cos(tm&latitude*rad)

tjs = wgt_areaave_Wrap(t(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
ujs = wgt_areaave_Wrap(u(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
vjs = wgt_areaave_Wrap(v(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
wjs = wgt_areaave_Wrap(w(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
rjs = wgt_areaave_Wrap(r(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
do i = 0, 7
    ip = i + 1
    nt = ind(phase.eq.ip.and.amp.gt.1.)
    nind = ind_DJF(nt)
    do j = 0, 40
        ntt   = -15 + j
        nindd = nind + ntt
        llmejs(i,j,0) = dim_avg_n(tjs(nindd), 0)
        llmejs(i,j,1) = dim_avg_n(ujs(nindd), 0)
        llmejs(i,j,2) = dim_avg_n(vjs(nindd), 0)
        llmejs(i,j,3) = dim_avg_n(wjs(nindd), 0)
        llmejs(i,j,4) = dim_avg_n(rjs(nindd), 0)
        delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nind/])
end do

delete([/clatjs,tjs,ujs,vjs,wjs,rjs/])
print("calc MEJS data done. "+ systemfunc("date"))

;; EAJS.
lleajs = new((/8,41,5/),float)
lleajs = 0.
;mejs
te = t(:,{elats:elatn},{elonw:elone})
clatjs = cos(te&latitude*rad)

tjs = wgt_areaave_Wrap(t(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
ujs = wgt_areaave_Wrap(u(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
vjs = wgt_areaave_Wrap(v(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
wjs = wgt_areaave_Wrap(w(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
rjs = wgt_areaave_Wrap(r(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
do i = 0, 7
    ip = i + 1
    nt = ind(phase.eq.ip.and.amp.gt.1.)
    nind = ind_DJF(nt)
    do j = 0, 40
        ntt   = -15 + j
        nindd = nind + ntt
        lleajs(i,j,0) = dim_avg_n(tjs(nindd), 0)
        lleajs(i,j,1) = dim_avg_n(ujs(nindd), 0)
        lleajs(i,j,2) = dim_avg_n(vjs(nindd), 0)
        lleajs(i,j,3) = dim_avg_n(wjs(nindd), 0)
        lleajs(i,j,4) = dim_avg_n(rjs(nindd), 0)
        delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nind/])
end do

delete([/clatjs,tjs,ujs,vjs,wjs,rjs/])
print("calc EAJS data done. "+ systemfunc("date"))

;------------------------------------------------------
;; plot format: Line chart. about 5 lines. RGBA color.
;; dtdt in (0,0,0,1.). udtdx in (112,160,205,1.0)
;; vdtdy in (0,52,102,1.) r in (178,178,178)


res = True
res@gsnDraw  = False
res@gsnFrame = False
;res@gsnRightStringFontHeightF  = 0.012
;res@gsnLeftStringFontHeightF   = 0.015

res@tmXTOn = False
res@tmYROn = False
;res@tmYRBorderOn = False
;res@tmXTBorderOn = False

res@tiXAxisString  = "Lag days"
res@tmXBMode   = "Explicit"
res@tmXBValues = (/ 0. , 5., 10., 15., 20., 25., 30., 35., 40. /)
res@tmXBLabels = (/"-15","-10","-5","0","5","10","15","20","25" /)

res@tiYAxisFontHeightF = 0.012
res@tiYAxisString = "Temperature variation (K/day)"
res@tmYLMode        = "Manual"
res@tmYLMinorPerMajor = 0
res@tmYLTickSpacingF  = 0.3
res@tmYLLabelFontHeightF = 0.015
res@tmYLMajorLengthF   = -0.01
res@trYMaxF  = 1.5
res@trYMinF  = -1.5
; t u v w r
color = (/"black","red","blue","green","purple"/)
xlabel = ispan(0,40,1)
;----------------------------------------------------------------
;; MEJS
wks = gsn_open_wks("eps", pdirmejs)

dum = new(5,graphic)
res = True
res@xyLineThicknessF = 2.0
res@xyLineColor = color(0)
res@gsnLeftString = "MEJS"
dum(0) = gsn_csm_xy(wks, xlabel, llmejs(0,:,0), res)
res@gsnLeftString = ""
do i = 1,4
  res@xyLineColor = color(i)
  dum(i) = gsn_xy(wks, xlabel, llmejs(0,:,i), res)
end do

overlay(dum(0),dum(1))
overlay(dum(0),dum(2))
overlay(dum(0),dum(3))
overlay(dum(0),dum(4))

draw(dum(0))
delete([/wks,dum/])

wks = gsn_open_wks("eps", pdireajs)

dum = new(5,graphic)
res = True
res@xyLineThicknessF = 2.0
res@xyLineColor = color(0)
res@gsnLeftString = "EAJS"
dum(0) = gsn_csm_xy(wks, xlabel, lleajs(0,:,0), res)
res@gsnLeftString = ""
do i = 1,4
  res@xyLineColor = color(i)
  dum(i) = gsn_xy(wks, xlabel, lleajs(0,:,i), res)
end do

overlay(dum(0),dum(1))
overlay(dum(0),dum(2))
overlay(dum(0),dum(3))
overlay(dum(0),dum(4))

draw(dum(0))
delete(wks)
















