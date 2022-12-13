undef("diffx")
function diffx(x[*][*][*][*]:numeric)
local df,dim,lat,lon,NLATS,NLONS,i,j,dx,rd
begin
  rd = 6378140
  dim = dimsizes(x)
  lat =  (x&latitude)
  lon =  (x&longitude)
  NLATS = dim(2)
  NLONS = dim(3)
  df = x
  do i = 0, NLATS - 1
    j = 0
    dx = (lon(j+1) - lon(j) + 360 - lon(NLONS - 1))*rd*2*3.14159/360
    dx = dx*cos(lat(i)/180.*3.14159)
    df(:,:,i,j) = (x(:,:,i,j+1) - x(:,:,i,NLONS - 1)) / dx
    do j = 1, NLONS - 2
      dx = (lon(j+1) - lon(j-1))*rd*2*3.14159/360
      dx = dx*cos(lat(i)/180.*3.14159)
      df(:,:,i,j) = (x(:,:,i,j+1) - x(:,:,i,j-1)) / dx
    end do
    j = NLONS - 1
    dx = (360 - lon(j-1))*rd*2*3.14159/360
    dx = dx*cos(lat(i)/180.*3.14159)
    df(:,:,i,j) = (x(:,:,i,0) - x(:,:,i,j-1)) / dx
  end do
  return df
end

undef("diffy")
function diffy(x[*][*][*][*]:numeric)
local df,dim,lat,lon,NLATS,i,dx,rd
begin
  rd = 6378140
  dim = dimsizes(x)
  lat =  (x&latitude)
  NLATS = dim(2)
  df = x
  i = 0
  ;lattitude N-->S; dy ahead north is positive
  dx = (lat(i) - lat(i+1))*rd*2*3.14159/360
  df(:,:,i,:) = (x(:,:,i,:) - x(:,:,i+1,:)) / dx
  do i = 1, NLATS - 2
    dx = (lat(i-1) - lat(i+1))*rd*2*3.14159/360
    df(:,:,i,:) = (x(:,:,i-1,:) - x(:,:,i+1,:)) / dx
  end do
  i = NLATS - 1
  dx = (lat(i-1) - lat(i))*rd*2*3.14159/360
  df(:,:,i,:) = (x(:,:,i-1,:) - x(:,:,i,:)) / dx
  return df
end

undef("diffp")
function diffp(x[*][*][*][*]:numeric)
local df,dim,lev,llev,NLEVS,i,dx,dx1,dx2
begin
  dim = dimsizes(x)
  lev =  (x&level*100.)
  llev = log(lev)
  NLEVS = dim(1)
  df = x

  i = 0
  dx = llev(i+1) - llev(i)
  df(:,i,:,:) = (x(:,i+1,:,:) - x(:,i,:,:)) / dx / lev(i)

  do i = 1, NLEVS - 2
    dx1 = llev(i) - llev(i-1)
    dx2 = llev(i+1) - llev(i-1)
    df(:,i,:,:) = ((x(:,i+1,:,:) - x(:,i,:,:))*dx1/dx2 + \
                   (x(:,i,:,:) - x(:,i-1,:,:))*dx2/dx1 ) / (dx1+dx2) / lev(i)
  end do

  i = NLEVS - 1
  dx = llev(i) - llev(i-1)
  df(:,i,:,:) = (x(:,i,:,:) - x(:,i-1,:,:)) / dx / lev(i)
  return df
end

undef("difft")
function difft(x[*][*][*][*]:numeric,y[*][*][*]:numeric)
local df,dim,lat,lon,NTIME,i,dy
begin
  dim = dimsizes(x)
  NTIME = dim(0)
  df = x
  dy = 1
  i = 0
  df(i,:,:,:) = ( x(i,:,:,:) - y) /24/3600
  do i = 1, NTIME - 1
    df(i,:,:,:) = (x(i,:,:,:) - x(i-1,:,:,:)) /24/3600/dy
  end do
  return df
end


begin
;************************************************
;define the parameter in this script
;************************************************
  SYEAR = 2014         ;start year
  EYEAR = 2014       ;end year
  FPATH = "./"
  NTIME = (EYEAR-SYEAR+1)*2
  R = 287.05
  Cp = 1005
;******************************************************
; define the Matrix
;******************************************************
  latS = 18
  latN =  22
  do  yr = SYEAR , EYEAR
    name = FPATH + "2014_JA_input.nc" ;"temp/air."+yr+".nc"

;     if ( yr .gt. SYEAR ) then
;       xyr = yr - 1
;       name1 = FPATH+"temp/air."+xyr+".nc"
;       lf = addfile(name1,"r")
;       ltmp = short2flt (lf->air(364,{1000:150},{latS:latN},:)) ; This part changed from MERRA Q1
;       LT = ltmp
;    ;  copyatt(LT,ltmp)
;     end if

    f = addfile(name,"r")
    list_filevars(f)
    tmp = short2flt (f->temp(:,{1000:100},{latS:latN},:))
    T = tmp
    copyatt(T,tmp)

    printVarSummary (T)

    delete(tmp)
    ; delete(f)
    PT = T
    printVarSummary (PT)
    dim = dimsizes(T)
    tt=dim(0)

    NLEVS = 18; to be counted from data from 1000hpa to 100hpa
    do i = 0, NLEVS - 1
      zl = ((/ T&level(i) /))
; print (zl)
      PT(:,i,:,:) = T(:,i,:,:)*((1000./zl)^(R/Cp))
    end do

    if ( yr .eq. SYEAR ) then
      y = PT(0,:,:,:) ;for first  year tend is forward differential
      ;different from Jian' program  see difft function line 93
      else              ; for other year last date of previous year is taken
      y = LT(:,:,:)   ; This part changed from MERRA Q1
    end if


    Q1 = difft(PT,y)
;      printVarSummary(Q1)

    ;for test
;      X = T(0,{1000},:,:)
;      printVarSummary(X)
;      N = num(ismissing(X))
;      print(N+"")
;      ;exit

    ; name = FPATH+"uwnd/uwnd."+yr+".nc"
    ; f= addfile(name,"r")
    tmp = short2flt (f->u(:,{1000:150},{latS:latN},:))
    U = tmp
;      U = tmp*tmp@scale_factor+tmp@add_offset
    copyatt(U,tmp)
    delete(tmp)
    ; delete(f)
    printVarSummary (U)
    Q1 = Q1 + diffx(PT)*U
    delete(U)
    ; name = FPATH+"vwnd/vwnd."+yr+".nc"
    ; f= addfile(name,"r")
    tmp = short2flt (f->v(:,{1000:150},{latS:latN},:))
    V = tmp
;      V = tmp*tmp@scale_factor+tmp@add_offset
    copyatt(V,tmp)
    delete(tmp)
    ; delete(f)
    printVarSummary (V)

    Q1 = Q1 + diffy(PT)*V
    delete(V)
    ; name = FPATH+"omega/omega."+yr+".nc"
    ; f= addfile(name,"r")
    tmp = short2flt (f->w(:,{1000:150},{latS:latN},:))
    W = tmp
;      W = tmp*tmp@scale_factor+tmp@add_offset
    copyatt(W,tmp)
    delete(tmp)
    delete(f)
    printVarSummary (W)

    Q1 = Q1 + diffp(PT)*W
    delete(W)

    Q1 = Q1*T/PT*3600
    dsizeQ = dimsizes(Q1)
    print(dsizeQ)

    copyatt(Q1,T)
    delete(T)
    delete(PT)

    printVarSummary(Q1)
    fQ1 = Q1
    lat = Q1&latitude
    lon = Q1&longitude
    lev = Q1&level*1.
    time = ispan(1,tt,1)*1.0
    print(lev)
    delete(Q1)
    fQ1!0 = "time"
    fQ1!1 = "lev"
    fQ1!2 = "lat"
    fQ1!3 = "lon"
;
    ; printVarSummary(fQ1)
    ; exit
    ; fQ1&lev = lev
    ; fQ1&time = time
    ; fQ1&lon = lon
    ; fQ1&lat = lat


    name = "Q1."+yr+".nc"
    system("rm -rf "+name)
    f = addfile(name,"c")
    f->Q1 = fQ1(:,:,{latS:latN},:)

    delete(f)
    delete(fQ1)
    delete(time)
    print("Year "+yr+" done!")
  end do
end
