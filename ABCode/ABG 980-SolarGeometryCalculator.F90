! TRNSYS Type 980: Solar Geometry Calculator (Incidence Angle, Tracking Angle and Sunrise, Sunset)
! ----------------------------------------------------------------------------------------------------------------------
!
! Try GitHub Try Try
! This routine calculates the sun position and the main solar angles as well as the sunrise and sunset. 
! The Excelergy solar geometry module has been converted to TRNSYS
!
! 
! ----------------------------------------------------------------------------------------------------------------------
! Nb | Parameters    | Description                                                    | Input  Units   | Internal Units 
! ---|---------------|----------------------------------------------------------------|----------------|----------------
!  1 | TimeZone      | Julian Day                                                     | 
!  2 | LongD         | Longitude
!  3 | LatD          | Latitude
!  4 | Altitude      | Altitude
!  5 | ColTilt       | Collector tilt or inclination
!  6 | ColAz         | Collector rotation

! ----------------------------------------------------------------------------------------------------------------------
! Nb | Inputs        | Description                                                    | Input  Units   | Internal Units 
! ---|---------------|----------------------------------------------------------------|----------------|----------------
!  1 | n             | Julian Day                                                     | 
!  2 | StdTime       | StandardTime                                                   | hr             | hr
!  3 | DNI           | Direct Normal Irradiation                                      | kJ/hr-m^2      | kJ/hr-m^2


! Outputs 
! ----------------------------------------------------------------------------------------------------------------------
! Nb | Variable      | Description                                                    | Output  Units 
! ---|---------------|----------------------------------------------------------------|---------------------------------
!  1 | ANI        | Product of cosine of theta and DNI                             | kJ/hr-m^2
!  2 | theta      | Incidence Angle
!  3 | beta       | Tracking angle
!  4 | SolarTime  | Solar Time
!  5 | Sunrise    | Sunrise
!  6 | Sunset     | Sunset
!  7 | SolarZenith | Solar zenith angle (converting SolarZenith_rad to degrees in output line)
!  8 | SolarAz     | Solar azimuth angle (converting SolarAz_rad to degrees in output line)

!
! Author:   Jenny L. Lohr & Antonio Gavilan
!   
!
!
! Revision history
! ---------------------------------------------------------------------------------------------------------------------
! AGM, March 14, 2012: solar geometry calculations added to the original component DNIxcos(theta) to actually calculate incidence angle, tracking angle, sunrise and sunset
!
! ----------------------------------------------------------------------------------------------------------------------
! Copyright ?2010 Abengoa Solar, Inc. All rights reserved.

subroutine type980(time,xin,out,t,dtdt,par,info,iCntrl,*)
!dec$attributes dllexport :: type980

use TrnsysConstants
use TrnsysFunctions
implicit none

! --- TRNSYS declarations ----------------------------------------------------------------------------------------------
integer :: NI,NP,ND,NO,NS    ! Nb of inputs, parameters, derivatives, outputs, storage locations
parameter (NI=3, NP =6, ND=0, NO=9, NS=0)
real(8), intent(in)    :: time, xin(NI), par(NP), t(*)
real(8), intent(inout) :: dtdt(*)
real(8), intent(out)   :: out(NO)
integer, intent(inout) :: info(15), iCntrl
character*3 YCHECK(NI),OCHECK(NO)

! --- Local variables --------------------------------------------------------------------------------------------------

integer :: mode,m,month,n, julian_day
real(8) :: stored(NS),time0,tfinal,delt

!Parameters
real(8) :: TimeZone,LongD,LatD,Altitude,ColTilt,ColAz

!Inputs
real(8) :: StdTime,DNI

!Outputs
real(8) :: SolarTime,theta,beta,ANI

!intermediate
real(8) :: B,B_rad,pi,EOT,dec ,dec_rad,StdLongD,SolarNoon,Lat_rad,DayLightHours
real(8) :: SunRise,SunSet,HourAngle,HourAngle_rad ,SolarAlt_rad,SolarAz_rad,SolarAz,costheta
real(8) :: SolarZenith_rad,SolarZenith,ColTilt_rad ,ColAz_rad,CosTh,Theta_rad,beta_rad,theta_rad_2,theta_2

! Get global TRNSYS simulation variables
time0  =getSimulationStartTime()
tfinal =getSimulationStopTime()
delt   =getSimulationTimeStep()

! --- Initial call to detect the TRNSYS version for which this Type is written -----------------------------------------
if (info(7) .eq. -2) then
    info(12) = 16
    return 1
endif

! --- All other calls --------------------------------------------------------------------------------------------------
! Always Read parameters
TimeZone =  par(1)
LongD =     par(2) 
LatD =      par(3) 
Altitude =  par(4)
ColTilt =   par(5) 
ColAz =     par(6) 

!Inputs

julian_day           = xin(1)
StdTime     = xin(2)
DNI         = xin(3)       ![kJ/hr-m^2]



! --- Second call in simulation: initialization call (not a simulation call) -------------------------------------------
if (info(7) .eq. -1) then
    info(9)  = 1     ! This type's outputs depend upon the passage of time
    info(6)  = NO    

    call typeck(1,info,NI,NP,ND)    ! Check the nb of inputs/parameters/derivatives
    call RCHECK(info,YCHECK,OCHECK) ! Check units

    return 1    ! Exit - End of the routine for very first call of the simulation
endif

! --- Third call in simulation: initial timestep manipulations, no iterations------------------------------------
if (time < (time0 + delt/2.) ) then

  !Verify values of parameters
  if (theta        < -180.)     call TYPECK(-4,info,0,2,0)
  if (theta        >  180.)     call TYPECK(-4,info,0,2,0)
  if (DNI    < 0.)     call TYPECK(-4,info,0,3,0)  
  
  !   Initialize outputs
    out(1)  = 0.
    out(2) = 0.
    out(3) = 0.
    out(4) = 0.
    out(5) = 0.
    out(6) = 0.
        
  return 1
endif
 
 ! --- Iterative calls -----------------------------------------------------------------------------------------
! Calculations

! Incidence angle and tracking angle calculations
pi = 3.1416

!2- Array of accumulated days for each month"

	
if(julian_day<=31)then
    month = 1
elseif(julian_day<=59)then
    month = 2
elseif(julian_day<=90)then
    month = 3    
elseif(julian_day<=120)then
    month = 4
elseif(julian_day<=151)then
    month = 5
elseif(julian_day<=181)then
    month = 6
elseif(julian_day<=212)then
    month = 7
elseif(julian_day<=243)then
    month = 8
elseif(julian_day<=273)then
    month = 9
elseif(julian_day<=304)then
    month = 10
elseif(julian_day<=334)then
    month = 11
else
    month = 12
endif
	
!2- Day of the year"
	n  = julian_day

!3- New B per Duffie & Beckman 1.4.2
	B = (n - 1) * 360. / 365. 
	B_rad = B * pi/180.	


!4- Equation of Time in minutes"
	EOT = 229.2 * (0.000075 + 0.001868 * Cos(B_rad) - 0.032077 * Sin(B_rad) - 0.014615 * Cos(B_rad * 2) - 0.04089 * Sin(B_rad * 2.))

!5- Declination   (per Duffie & Beckman 1.6.1a)"
	dec = 23.45 * Sin(360. * (284. + n) / 365. * Pi / 180.)
	dec_rad = dec * pi / 180.

!6- Standard longitude
	StdLongD = TimeZone * (15.) !Now works with eastern hemisphere locations
	!StdLongD = TimeZone * (-15.)

!7- Calculation of the Solar Noon in hours
	SolarNoon = 12. + (StdLongD - LongD) / 15. - EOT / 60. !Now works with eastern hemisphere locations
	!SolarNoon = 12. - (StdLongD - LongD) / 15. - EOT / 60.

!8- Latitude in radians
	Lat_rad = LatD * pi / 180.

!9- Calculation of the number of daylight hours"
	DayLightHours = 2. / 15. * ACOS(-TAN(Lat_rad) * TAN(Dec_rad)) * 180. / Pi
  
!10- Calculation of Sun Rise and Sun Set in hours"
	SunRise = SolarNoon - DayLightHours / 2.
	SunSet = SolarNoon + DayLightHours / 2.

!11- Solar time - in hours
	!SolarTime = StdTime + (StdLongD - LongD) / 15. + EOT / 60.
	SolarTime = StdTime - (StdLongD - LongD) / 15. + EOT / 60. !This should be modified KQW
!12- Calculation of Hour Angle in radians
	HourAngle = (SolarTime - 12.) * 15. 
	HourAngle_rad = HourAngle * pi / 180.

!13- Solar Altitude   (Radians)
	SolarAlt_rad = ASIN( Sin(dec_rad) * Sin(Lat_rad) + Cos(Lat_rad) * Cos(Dec_rad) * Cos(HourAngle_rad))

!14- Solar azimuth"
	SolarAz_rad = ACOS((Sin(Dec_rad) * Cos(Lat_rad) - Cos(Dec_rad) * Cos(HourAngle_rad) * Sin(Lat_rad)) / Cos(SolarAlt_rad))	
	SolarAz = SolarAz_rad * 180./pi
    !if (HourAngle > 0.) then
    !    SolarAz_rad = 2.*pi - SolarAz_rad
    !else
    !    SolarAz_rad = SolarAz_rad
    !endif

!15- Solar Zenith
	SolarZenith_rad = Pi / 2. - SolarAlt_rad
    SolarZenith = SolarZenith_rad * 180./pi
!16- Calculation of Solar Incidence Angle for Trough"
	!Stine Reference
	ColTilt_rad = ColTilt * pi / 180.
	ColAz_rad = ColAz * pi / 180.

	CosTh = sqrt(1. - (Cos(SolarAlt_rad - ColTilt_rad) - Cos(ColTilt_rad) * Cos(SolarAlt_rad) * (1 - Cos(SolarAz_rad - ColAz_rad))) ** 2.)
      
	Theta_rad = ACOS(CosTh)
	Theta = Theta_rad * 180./pi
     
!17- Calculation of Tracking Angle for Trough"
	!Stine Reference

	beta_rad = ATAN(Cos(SolarAlt_rad) * Sin(SolarAz_rad - ColAz_rad) / (Sin(SolarAlt_rad - ColTilt_rad) + Sin(ColTilt_rad) * Cos(SolarAlt_rad) * (1. - Cos(SolarAz_rad - ColAz_rad))))
	beta = beta_rad * 180. / pi 

!18- Calculation of Tracking Angle for Trough"    
    ANI= DNI * COSD(theta)

!19 - Convert Solar azimuth angle into the following convention: Referenced from due-south, (+) west of south, i.e. CW  & (-) east of south, i.e. CCW
	SolarAz_rad = (Cos(SolarZenith_rad) * Sin(Lat_rad) - Sin(Dec_rad)) / (Sin(SolarZenith_rad) * Cos(Lat_rad))
	if (ABS(SolarAz_rad) > 1 ) then !this value can not be beyond range of -1 to 1, becuase the arc_cosine will be taken
	    SolarAz_rad = 1
	endif
	SolarAz_rad = ACOS( SolarAz_rad )	 !Equation 1.6.6 Duffie and Beckman
    
    
    If (HourAngle  < 0 ) then
        SolarAz_rad =  - SolarAz_rad
    EndIf
    
!Outputs
out(1) = ANI       ![kJ/hr-m^2]    
out(2) = theta
out(3) = beta
out(4) = SolarTime
out(5) = SunRise 
out(6) = SunSet
out(7) = SolarZenith_rad * 180./pi  !Solar zenith angle in degrees
out(8) = SolarAz_rad * 180./pi  !Solar azimuth angle in degrees
out(9) = month


	
return 1

end subroutine type980