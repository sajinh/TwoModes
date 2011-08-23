load "$SysE/lib/ncl/helper_libs.ncl"
load "$SysE/lib/ncl/data_analysis_helpers.ncl"

; The script expects to be provided with a file containing an input
; variable containing coordinates time, lat and lon in the order
; mentioned here
; A climatology is constructed for each year of the data based on
; annual harmonics
; This climatology is removed from the data to yield daily anomalies
begin
  var="olr"
  infile=var+"_fill.nc"
  outfile=var+"_hig.nc"

  nhar = 4 ; number of harmonics to be used to construct annual cycle

  fin=addfile(infile,"r")
  varx=fin->$var$(:,{-30:30},{40:120})
  time=varx&time
  ntim=dimsizes(time)

  ssta=varx(lat|:,lon|:,time|:)
  ; Remove any linear trend present in data
  ssta=dtrend_msg(time,varx(lat|:,lon|:,time|:),True,False)
  varx=ssta(time|:,lat|:,lon|:)
  var_raw=varx

  year1=ut_calendar(time(0),0)
  year2=ut_calendar(time(ntim-1),0)
  yr1=floattointeger(year1(:,0))
  yr2=floattointeger(year2(:,0))
  mo1=floattointeger(year1(:,1))
  mo2=floattointeger(year2(:,1))
  dy1=floattointeger(year1(:,2))
  dy2=floattointeger(year2(:,2))
  t1=0
  t2=0
  print(yr1+":"+yr2)

  ; this complexity is to add the ability to provide data
  ; that need not have complete data for one year at the end-points
  if (yr1.ne.yr2) then
    ; For the first year, time goes from tim1 to the end of the year
    find_t1t2(yr1,mo1,dy1,yr1,12,31,time,t1,t2)
    varx(t1:t2,:,:)=constructed_anncycle(varx(t1:t2,:,:), nhar)

    do iyears=yr1+1,yr2-1
      find_t1t2(iyears,01,01,iyears,12,31,time,t1,t2)
      varx(t1:t2,:,:)=constructed_anncycle(varx(t1:t2,:,:), nhar)
    end do

    ; For the last year, time goes from beginning of year to tim2
    find_t1t2(yr2,01,01,yr2,mo2,dy2,time,t1,t2)
    varx(t1:t2,:,:)=constructed_anncycle(varx(t1:t2,:,:), nhar)
  else
	; If year1 and year2 are same (there is only one year of data
    find_t1t2(yr1,mo1,dy1,yr2,mo2,dy2,time,t1,t2)
    varx(t1:t2,:,:)=constructed_anncycle(varx(t1:t2,:,:), nhar)
  end if

  ; varx contains the annual cycle of the data
  ; we remove this to get the higher frequency daily anomalies in
  ; the data set
  var_raw=var_raw-varx
  fo=addfile(outfile,"c")
  fo->$var$=var_raw
end
