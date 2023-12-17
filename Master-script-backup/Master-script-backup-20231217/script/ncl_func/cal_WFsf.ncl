;===================================================================================================
; File:         cal_WFsf.ncl
; Category:     NCAR command language script
; Version:      01.00.0000
; Author(s):    Ming Sun (gosun1994@gmail.com)
; Copyright:    (C) 2020 Ming Sun
; License:      GNU General Public License version 3 (GPLv3)
; Date Created: 2020-02-04 by Ming Sun
; Last Updated: 2020-02-04 by Ming Sun
;---------------------------------------------------------------------------------------------------
; Description:  read the stream function and mean u;v to calculate wave activity flux
;===================================================================================================
	
undef("cal_WFsf")

function cal_WFsf( psi[*][*]:float, um[*][*]:float, vm[*][*]:float )
local DIM,nlatitude,nlongitude,latitude,longitude,\\
      cumag,dpsidlongitude,ddpsidlongitudelongitude,dpsidlatitude,ddpsidlongitudelatitude,ddpsidlatitudelatitude,\\
      xuterm,xvterm,yvterm,dlongitude,dlatitude,nl,dX,ml,dY,px,py,out

begin
  DIM        = dimsizes(psi)
  nlatitude  = DIM(0)
  nlongitude = DIM(1)

  latitude   = psi&latitude
  longitude  = psi&longitude

  cumag = sqrt(um^2 + vm^2)
  copy_VarMeta(um, cumag)
  cumag@_FillValue = um@_FillValue

;--------Psi_x Psi_y Psi_xx Psi_yy Psi_xy------------------------------------
 dpsidlongitude           = new ((/nlatitude,nlongitude/), float, psi@_FillValue) ;ntime*nlatitude*nlongitude
 ddpsidlongitudelongitude = new ((/nlatitude,nlongitude/), float, psi@_FillValue)
 dpsidlatitude            = new ((/nlatitude,nlongitude/), float, psi@_FillValue)
 ddpsidlongitudelatitude  = new ((/nlatitude,nlongitude/), float, psi@_FillValue)
 ddpsidlatitudelatitude   = new ((/nlatitude,nlongitude/), float, psi@_FillValue)
;--------3 terms consists of Partial Psi------------------------------------
 xuterm          = new ((/nlatitude,nlongitude/), float, psi@_FillValue)
 xvterm          = new ((/nlatitude,nlongitude/), float, psi@_FillValue)
 yvterm          = new ((/nlatitude,nlongitude/), float, psi@_FillValue)
;------Here dx dy in rads-------------------- 
 dlongitude = (longitude(2)-longitude(1))*0.0174533 ;; one degree equals to how many rads
 dlatitude = (latitude(2)-latitude(1))*0.0174533

 do nl=0,nlatitude-1
   dX = 6371000.*cos(0.0174533*latitude(nl))*dlongitude ;; dx in distance  
   dpsidlongitude(nl:nl,:)           =  doubletofloat(center_finite_diff_n(psi(nl:nl,:),dX,True,0,1))
   ddpsidlongitudelongitude(nl:nl,:) =  doubletofloat(center_finite_diff_n(dpsidlongitude(nl:nl,:),dX,True,0,1))
 end do 

 do ml=0,nlongitude-1
   dY = 6371000.*dlatitude
   dpsidlatitude(:,ml:ml)            = doubletofloat(center_finite_diff_n(psi(:,ml:ml), dY, True,0,0))
   ddpsidlongitudelatitude(:,ml:ml)  = doubletofloat(center_finite_diff_n(dpsidlongitude(:,ml:ml),dY,True,0,0))
   ddpsidlatitudelatitude(:,ml:ml)   = doubletofloat(center_finite_diff_n(dpsidlatitude(:,ml:ml),dY,True,0,0))
 end do

xuterm = (dpsidlongitude*dpsidlongitude - psi*ddpsidlongitudelongitude)
xvterm = (dpsidlongitude*dpsidlatitude  - psi*ddpsidlongitudelatitude)
yvterm = (dpsidlatitude*dpsidlatitude   - psi*ddpsidlatitudelatitude)

px = (um*xuterm + vm*xvterm)/(2.*cumag)
py = (um*xvterm + vm*yvterm)/(2.*cumag)

copy_VarCoords(psi, px)
copy_VarCoords(psi, py)

  out = new((/2,nlatitude,nlongitude/),typeof(psi))
  out(0,:,:) = px
  out(1,:,:) = py

  return (out)

end