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
mlonw = 60
mlone = 80

elatn = 35
elats = 25
elonw = 170
elone = 190

ddir = "/home/sunming/data5/cuixy/global_jets/data/temp_budget"
tdir = ddir+"/"+"anom_dtdt_r2.5.nc"
udir = ddir+"/"+"anom_udtdx_r2.5.nc"
vdir = ddir+"/"+"anom_vdtdy_r2.5.nc"
wdir = ddir+"/"+"anom_wdtdp_r2.5.nc"
;rdir = ddir+"/"+"anom_wrtcp_r2.5.nc"
adir = ddir+"/anom_udtdx_seperate.nc"
pdireajs = "/home/sunming/data5/cuixy/global_jets/plot/temp_budget"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
;-----------------------------------------------------------------
;;seperate
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
;-----------------------------------------------------------------
ft  = addfile(tdir,"r")
fu  = addfile(udir,"r")
fv  = addfile(vdir,"r")
fw  = addfile(wdir,"r")
;fr  = addfile(rdir,"r")

t = ft->dtdt
t = t
u = fu->udtdx
u = u*tday
v = fv->vdtdy
v = v*tday
w = fw->wdtdp
w = w*tday
;r1 = fr->wrtcp
;r1 = r1*tday
;residual
r = t + (u + v + w) ;- r1)
copy_VarCoords(t, r)
u = -1*u
v = -1*v
w = -1*w;+r1
lat = t&lat
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
;;MEJS.
llmejs = new((/8,41,5/),float)
llmejs = 0.
tm = t(:,{mlats:mlatn},{mlonw:mlone})
clatjs = cos(tm&lat*rad)
print(clatjs)
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
;; seperate.
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

;seperate.
lsmejs = new((/10/),float)
lsmejs = 0.
i = 0
  ip = i + 1
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind  = ind_DJF(nt)
  ntt   = 2
  nin = nind + ntt
  do j = 1, 5
    ntt   = 2 + j
    nindd = nind + ntt
    nin  := array_append_record(nin, nindd, 0)
  end do
  lsmejs(0) = dim_avg_n(ujs(nin), 0)
  lsmejs(1) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  lsmejs(2) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  lsmejs(3) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  lsmejs(4) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  lsmejs(5) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  lsmejs(6) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  lsmejs(7) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  lsmejs(8) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  lsmejs(9) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)

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
te = t(:,{elats:elatn},{elonw:elone})
clatjs = cos(te&lat*rad)

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
;; seperate.
ultljs = wgt_areaave_Wrap(ultl(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
ultijs = wgt_areaave_Wrap(ulti(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
ulthjs = wgt_areaave_Wrap(ulth(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)

uitljs = wgt_areaave_Wrap(uitl(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
uitijs = wgt_areaave_Wrap(uiti(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
uithjs = wgt_areaave_Wrap(uith(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)

uhtljs = wgt_areaave_Wrap(uhtl(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
uhtijs = wgt_areaave_Wrap(uhti(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
uhthjs = wgt_areaave_Wrap(uhth(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
;seperate.
delete(nin)
lseajs = new((/10/),float)
lseajs = 0.
i = 0
  ip = i + 1
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind  = ind_DJF(nt)
  ntt   = -5
  nin   = nind + ntt
  do j = 1, 7
    ntt   = -5 + j
    nindd = nind + ntt
    nin  := array_append_record(nin, nindd, 0)
  end do
  lseajs(0) = dim_avg_n(ujs(nin), 0)
  lseajs(1) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  lseajs(2) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  lseajs(3) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  lseajs(4) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  lseajs(5) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  lseajs(6) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  lseajs(7) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  lseajs(8) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  lseajs(9) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)
print(lseajs)
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
res@gsnRightStringFontHeightF  = 0.012
res@gsnLeftStringFontHeightF   = 0.015

res@tmXTOn = False
res@tmYROn = False

res@tiXAxisString  = "Lag days"
res@tmXBMode   = "Explicit"
res@tmXBValues = (/ 0. , 5., 10., 15., 20., 25., 30., 35., 40. /)
res@tmXBLabels = (/"-15","-10","-5","0","5","10","15","20","25" /)
res@tmXBLabelFontHeightF = 0.012
res@tiXAxisFontHeightF = 0.012

res@tiYAxisFontHeightF = 0.012

res@tmYLMode        = "Manual"
res@tmYLMinorPerMajor = 0
res@tmYLTickSpacingF  = 0.3
res@tmYLLabelFontHeightF = 0.012

res@tmYLMajorLengthF   = -0.01
res@trYMaxF  = 1.5
res@trYMinF  = -1.5
res@gsnYRefLine = 0.
res@gsnYRefLineDashPattern = 1

res@xyLineThicknessF = 4.0
res@vpWidthF  = 0.38
res@vpHeightF = 0.4
; t u v w r
color = (/"black","royalblue4","skyblue3","darkgreen",\
         "darkgoldenrod4"/)  ; gray69
print(color) 
xlabel = ispan(0,40,1)

;; bar chart.
resb = True
resb@gsnDraw     = False
resb@gsnFrame    = False
resb@gsnMaximize = False

resb@gsnLeftStringFontHeightF   = 0.012

resb@tmXTOn = False
resb@tmYROn = False

resb@tiYAxisFontHeightF = 0.012
resb@tmYLMode           = "Manual"
resb@tmYLMinorPerMajor  = 0
resb@tmYLTickSpacingF   = 0.1
resb@tmYLLabelFontHeightF = 0.012

resb@tmXBMode    = "Explicit"
resb@tmXBValues  = (/ 0.5 , 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5,5. /)
resb@tmXBLabels  = (/"~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'","~F21~LL","~F21~LI",\
                    "~F21~LH","~F21~IL","~F21~II","~F21~IH","~F21~HL","~F21~HI","~F21~HH"/)
resb@tmXBLabelFontHeightF = 0.012

resb@gsnXYBarChart = True
resb@gsnXYBarChartBarWidth = 0.3

resb@tmYLMajorLengthF   = -0.01
resb@gsnYRefLine        = 0.
resb@trXMinF  = 0.2
resb@trXMaxF  = 5.3
;----------------------------------------------------------------
;; EAJS
wks = gsn_open_wks("x11", pdireajs)
;drawNDCGrid(wks)
dum = new(5,graphic)
res@vpXF = 0.59
res@vpYF = 0.95
res@trYMaxF  = 1.5
res@trYMinF  = -1.5
res@tmYLTickSpacingF  = 0.3
res@gsnLeftString = "~F21~(b) EAJS"
res@gsnRightString = "~F21~ave: "+elats+"-"+elatn+"N,"+elonw+"-"+elone+"E"
res@xyLineColor = color(0)
dumea = gsn_csm_xy(wks, xlabel, lleajs(0,:,0), res)
res@gsnLeftString = ""
res@gsnRightString = ""
do i = 1,4
  res@xyLineColor = color(i)
  dum(i) = gsn_xy(wks, xlabel, lleajs(0,:,i), res)
end do

overlay(dumea,dum(1))
overlay(dumea,dum(2))
overlay(dumea,dum(3))
overlay(dumea,dum(4))
;;draw two lines.
resl = True
resl@gsLineColor = "gray64"
resl@gsLineDashPattern = 2
resl@gsLineThicknessF  = 2.0   
dum11 = gsn_add_polyline(wks, dumea, (/27,27/), (/-10,10/), resl)
dum12 = gsn_add_polyline(wks, dumea, (/32,32/), (/-10,10/), resl)

;; text.
rest               = True            
rest@txFontHeightF = 0.012    
text1  = "~F21~(~F34~6~F21~T/~F34~6~F21~t~F21~)'"                 
text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
text3  = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
text4  = "~F21~(~F33~w~F33~s~F21~)'"
text5  = "~F21~Q'/c~B~p~N~"
rest@txFontColor = color(0)
gsn_text_ndc(wks, text1, 0.15, 0.72, rest)
rest@txFontColor = color(1)
gsn_text_ndc(wks, text2, 0.20, 0.89, rest)
rest@txFontColor = color(2)
gsn_text_ndc(wks, text3, 0.13, 0.82, rest)
rest@txFontColor = color(3)
gsn_text_ndc(wks, text4, 0.2, 0.64, rest)
rest@txFontColor = color(4)
gsn_text_ndc(wks, text5, 0.12, 0.77, rest)
draw(dumea)

dum = new(5,graphic)
res@vpXF = 0.09
res@vpYF = 0.95
res@trYMaxF  = 1.5
res@trYMinF  = -1.5
res@tmYLTickSpacingF  = 0.3
res@tiYAxisString = ""
res@gsnLeftString = "~F21~(a) MEJS"
res@tiYAxisString = "Temperature tendency (K/day)"
res@gsnRightString = "~F21~ave: "+mlats+"-"+mlatn+"N,"+mlonw+"-"+mlone+"E"
res@xyLineColor = color(0)

dumme = gsn_csm_xy(wks, xlabel, llmejs(0,:,0), res)
res@gsnLeftString = ""
res@gsnRightString = ""
do i = 1,4
  res@xyLineColor = color(i)
  dum(i) = gsn_xy(wks, xlabel, llmejs(0,:,i), res)
end do

overlay(dumme,dum(1))
overlay(dumme,dum(2))
overlay(dumme,dum(3))
overlay(dumme,dum(4))
dum21 = gsn_add_polyline(wks, dumme, (/17,17/), (/-10,10/), resl)
dum22 = gsn_add_polyline(wks, dumme, (/22,22/), (/-10,10/), resl)
rest@txFontColor = color(0)
gsn_text_ndc(wks, text1, 0.71, 0.7, rest)
rest@txFontColor = color(1)
gsn_text_ndc(wks, text2, 0.72, 0.58, rest)
rest@txFontColor = color(2)
gsn_text_ndc(wks, text3, 0.79, 0.85, rest)
rest@txFontColor = color(3)
gsn_text_ndc(wks, text4, 0.72, 0.87, rest)
rest@txFontColor = color(4)
gsn_text_ndc(wks, text5, 0.65, 0.76, rest)
draw(dumme)


;;draw Bar chart.
resb@trYMaxF  = 0.6
resb@trYMinF  = -2.0
resb@tmYLTickSpacingF   = 0.3
resb@vpXF = 0.1
resb@vpYF = 0.47
resb@vpWidthF  = 0.8
resb@vpHeightF = 0.18
resb@tmXBLabelsOn = True
resb@gsnLeftString = "~F21~(c) EAJS"
dumeabar = gsn_csm_xy(wks,fspan(0.5,5.,10),lseajs,resb)
;;polyline.
respoly = True
respoly@gsLineColor = "gray64"
respoly@gsLineDashPattern = 2
respoly@gsLineThicknessF  = 5.0
;dum3 = gsn_add_polyline(wks,dumeabar,(/0.7,0.7/),(/0.3,-0.5/),respoly)
;; shade color for lseajs(0) and lleajs(2).
resea = True
resea@gsFillColor = color(1);
dum4 = gsn_add_polygon(wks, dumeabar, (/0.35,0.65,0.65,0.35,0.35/), (/0.,0.,lseajs(0),lseajs(0),0./), resea)
dum5 = gsn_add_polygon(wks, dumeabar,(/1.35,1.65,1.65,1.35,1.35/), (/0.,0.,lseajs(2),lseajs(2),0./), resea)
;draw(dumeabar)
;;mejs
resb@trYMaxF  = 1.5
resb@trYMinF  = -0.1
resb@tmYLTickSpacingF   = 0.3
resb@vpXF = 0.1
resb@vpYF = 0.22
resb@vpWidthF  = 0.8
resb@vpHeightF = 0.18
resb@gsnLeftString = "~F21~(d) MEJS"
resb@gsnYRefLine        = 0.
resb@tmXBLabelsOn  = True
dummebar = gsn_csm_xy(wks,fspan(0.5,5.,10),lsmejs,resb)
;;polyline.
respoly = True
respoly@gsLineColor = "gray64"
respoly@gsLineDashPattern = 2
respoly@gsLineThicknessF  = 5.0
dum3 = gsn_add_polyline(wks,dummebar,(/0.7,0.7/),(/1.5,-0.1/),respoly)
;; shade color for lseajs(0) and lleajs(2).
resea = True
resea@gsFillColor = color(1)
dum4 = gsn_add_polygon(wks, dummebar, (/0.35,0.65,0.65,0.35,0.35/), (/0.,0.,lsmejs(0),lsmejs(0),0./), resea)
dum5 = gsn_add_polygon(wks, dummebar,(/1.35,1.65,1.65,1.35,1.35/), (/0.,0.,lsmejs(2),lsmejs(2),0./), resea)
;draw(dummebar)

i = 1
do while (i.eq.1)
    i = 1
end do
