;--------------------------------------------------------
;; spectral analysis.
;; original data remove series and dtrend.
;; by Backward self-correlation.
;; smoothed by Hamming window.
;; w(n) = a - (1-a)*cos(2pi*n/(N-1)) a=0.54
;--------------------------------------------------------
pi = atan(1.0) * 4.
function spec_anal(x[*]:float)
local x
begin
;; some constants.
  a  = 0.54
  x := dim_standardize(x, 0) ; standard.
  n  = dimsizes(x)           ; time length.
  M  = toint(n/3.)           ; max backward step.
  sigma = stddev(x)          ; standard deviation.
  R = new((/M+1/),float)     ; backward correlation.
  R = 0.
  xmean = dim_avg(x)         ; averaged x.
  y  = R                     ; raw spectral.
  ys = new((/2,dimsizes(R)/),float)
  ys = 0.
;; backward correlation.
  do tao = 0, M
    tmax  = n-M-1
    rtemp = 0.
    do j = 0, tmax
      rtemp = ((x(j)-xmean)/sigma)*((x(j+tao)-xmean)/sigma) + rtemp
    end do
    R(tao) = rtemp/(n-tao)
  end do

;; calculate raw spectral.
  do k = 1, M-1
    ytemp = 0.
    do j = 1, M-1
      ytemp = R(j)*cos(pi*k*j/M) + ytemp
    end do
    y(k) = (R(0)+2*ytemp+R(M)*cos(pi*k))/M
  end do

  ytmp = 0.
  k = 0
  do j = 1, M-1
    ytemp = R(j)*cos(pi*k*j/M) + ytemp
  end do
  y(k) = (R(0)+2*ytemp+R(M)*cos(pi*k))/(2*M)

  ytmp = 0.
  k = M
  do j = 1, M-1
    ytemp = R(j)*cos(pi*k*j/M) + ytemp
  end do
  y(k) = (R(0)+2*ytemp+R(M)*cos(pi*k))/(2*M)

;; smooth raw spectral.
  do k = 0, M
    w = a - (1.-a)*cos(2*pi*k/M)
    ys(0,k) = w*y(k)
  end do

;; calc noisy spectral.
  talfa = 1.645
  kar2  = 11.07
  Ralpa = ( -1.+ talfa*sqrt(n-2) )/(n-1)

  if (R(1).gt.Ralpa) then
    do k = 0, M
      ys(1,k) = ys(0,k)*(1-R(1)*R(1))/(1+R(1)*R(1)-2*R(1)*cos(pi*k/M))\
      *kar2*M/(2.*n-M/2.)
    end do
  else
    do k = 0, M
      ys(1,k) = ys(0,k)*kar2*M/(2.*n-M/2.)
    end do
  end if

  return (ys)

end 