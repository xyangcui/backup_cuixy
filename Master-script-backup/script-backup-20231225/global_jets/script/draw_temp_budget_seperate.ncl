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
adir = ddir+"/anom_udtdx_seperate.nc"
pdirmejs = "/home/sunming/data5/cuixy/global_jets/plot/MEJS_budget"
pdireajs = "/home/sunming/data5/cuixy/global_jets/plot/EAJS_budget"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------

f  = addfile(adir,"r")

ultl = f->ultl
ulti = f->ulti
ulth = f->ulth
uitl = f->uitl
uiti = f->uiti
uith = f->uith
uhtl = f->uhtl
uhti = f->uhti
uhth = f->uhth

;---------------------------------------------------------------=
;; MJO phases.
time = f->time
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
;; EAJS.
lleajs = new((/9/),float)
lleajs = 0.
;mejs
te = ultl(:,{elats:elatn},{elonw:elone})
clatjs = cos(te&latitude*rad)

;ultljs = wgt_areaave_Wrap(ultl(:,{elats:elatn},{elonw:elone})\
;                         , clatjs, 1.0, 1)
;ultijs = wgt_areaave_Wrap(ulti(:,{elats:elatn},{elonw:elone})\
;                         , clatjs, 1.0, 1)
;ulthjs = wgt_areaave_Wrap(ulth(:,{elats:elatn},{elonw:elone})\
;                         , clatjs, 1.0, 1)

;uitljs = wgt_areaave_Wrap(uitl(:,{elats:elatn},{elonw:elone})\
;                         , clatjs, 1.0, 1)
;uitijs = wgt_areaave_Wrap(uiti(:,{elats:elatn},{elonw:elone})\
;                         , clatjs, 1.0, 1)
;uithjs = wgt_areaave_Wrap(uith(:,{elats:elatn},{elonw:elone})\
;                         , clatjs, 1.0, 1)

;uhtljs = wgt_areaave_Wrap(uhtl(:,{elats:elatn},{elonw:elone})\
;                         , clatjs, 1.0, 1)
;uhtijs = wgt_areaave_Wrap(uhti(:,{elats:elatn},{elonw:elone})\
;                         , clatjs, 1.0, 1)
;uhthjs = wgt_areaave_Wrap(uhth(:,{elats:elatn},{elonw:elone})\
;                         , clatjs, 1.0, 1)
delete([/te,clatjs/])
;;MEJS
llmejs = new((/9/),float)
llmejs = 0.
;mejs
te = ultl(:,{mlats:mlatn},{mlonw:mlone})
clatjs = cos(te&latitude*rad)

ultljs = wgt_areaave_Wrap(ultl(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
ultijs = wgt_areaave_Wrap(ulti(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
ulthjs = wgt_areaave_Wrap(ulth(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)

uitljs = wgt_areaave_Wrap(uitl(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
uitijs = wgt_areaave_Wrap(uiti(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
uithjs = wgt_areaave_Wrap(uith(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)

uhtljs = wgt_areaave_Wrap(uhtl(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
uhtijs = wgt_areaave_Wrap(uhti(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
uhthjs = wgt_areaave_Wrap(uhth(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)

;; -5 -- 0?
i = 0
  ip = i + 1
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind  = ind_DJF(nt)
  ntt   = -5
  nin = nind + ntt
  do j = 1, 5
    ntt   = -5 + j
    nindd = nind + ntt
    nin  := array_append_record(nin, nindd, 0)
  end do
  lleajs(0) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  lleajs(1) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  lleajs(2) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  lleajs(3) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  lleajs(4) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  lleajs(5) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  lleajs(6) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  lleajs(7) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  lleajs(8) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)
;print(lleajs)
print("calc EAJS data done. "+ systemfunc("date"))

i = 0
  ip = i + 1
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind  = ind_DJF(nt)
  ntt   = -5
  nin = nind + ntt
  do j = 1, 20
    ntt   =  ntt + j
    nindd = nind + ntt
    nin  := array_append_record(nin, nindd, 0)
  end do
  llmejs(0) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  llmejs(1) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  llmejs(2) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  llmejs(3) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  llmejs(4) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  llmejs(5) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  llmejs(6) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  llmejs(7) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  llmejs(8) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)
print(llmejs)
print("calc MEJS data done. "+ systemfunc("date"))

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

wks = gsn_open_wks("x11", pdireajs)

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
















