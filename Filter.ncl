;-------------------------------------------------------------------
; This function is used to filter.
; wts: weight for Lanczos filter.
; xClmDay: climatology
; xClmDay_sm: smoothed climatology, subtracting first 3 harmonics.
; xAnom: daily anomalies.
; x1: filtered result.
; npts: notes of lanczos filter. default:201
; nwa: larger value.  default: 100
; nwb: smaller value. default: 20
;-------------------------------------------------------------------
function iso_filter(x,npts,nwa,nwb)
local time,ymdhms,iyear,imon,iday,ddd,yyyyddd,wts,xClmDay,xClmDay_sm
local xAnom,x1,npts,nwa,nwb,x
begin
    npts = 201
    nwa  = 100
    nwb  = 20
    time    = x&time
    ymdhms  = cd_calendar(time, 0)
    iyear   = floattoint( ymdhms(:,0) )
    imon    = floattoint( ymdhms(:,1) )   
    iday    = floattoint( ymdhms(:,2) )
    ddd     = day_of_year(iyear, imon, iday) 
    yyyyddd = iyear*1000 + ddd 
    wts     = filwgts_lanczos (npts, 2, 1./nwa,1./nwb , 1.0)

    xClmDay    = clmDayTLL(x, yyyyddd)
    xClmDay_sm = smthClmDayTLL(xClmDay, 3)
    xAnom      = calcDayAnomTLL(x, yyyyddd, xClmDay_sm)
    x1         = wgt_runave_n_Wrap(xAnom, wts, 0, 0)

    return(x1)
end
