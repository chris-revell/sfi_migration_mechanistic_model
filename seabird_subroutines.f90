! Christopher Revell, University of Cambridge, 2016
MODULE SEABIRD_SUBROUTINES

  IMPLICIT NONE

CONTAINS

  subroutine forcecalc(force,currentlat,currentlon,resources_filtered)!wind_zonal_current,wind_merid_current,resources_filtered,a)
    integer*4                             :: i,j
    real*4,intent(in)                     :: currentlat,currentlon
  !  real*4,intent(in)                     :: a,wind_zonal_current,wind_merid_current
    real*4,intent(in),dimension(720,1440) :: resources_filtered
    real*4,intent(inout),dimension(2)     :: force
    integer*4,dimension(2)                :: resources_shape
    real*4                                :: resourcelat,resourcelon,bearing,d_latlong,forcemagnitude,dist

    resources_shape = (/720,1440/)
    d_latlong = 180.0/720.0
    force = (/0.0,0.0/)

    do i=1,resources_shape(1)
      do j=1,resources_shape(2)
        resourcelat = (resources_shape(1)/2.0-i-0.5)*d_latlong
        resourcelon = (j+0.5-resources_shape(2)/2.0)*d_latlong
        bearing = initialbearing(currentlat,currentlon,resourcelat,resourcelon)
        dist = realdistance(resourcelat,resourcelon,currentlat,currentlon)
        if (dist.LT.0.00000001) then
          CYCLE
        else
          forcemagnitude = resources_filtered(i,j)/dist
          force(1) = force(1)+forcemagnitude*cos(bearing)
          force(2) = force(2)+forcemagnitude*sin(bearing)
        endif
      enddo
    enddo
  end subroutine forcecalc

  real*4 function realdistance(lata,lona,latb,lonb)
    real*4,intent(in)    :: lata,lona,latb,lonb
    real*4               :: delta_long,term1,term2
    real*4               :: pi  = 4*atan(1.0)

    delta_long = lonb-lona
    term1      = sin(lata*pi/180.0)*sin(latb*pi/180.0)
    term2      = cos(lata*pi/180.0)*cos(latb*pi/180.0)*cos(delta_long*pi/180.0)
    if ((term1+term2).GT.1.0) then
        realdistance = 6371.0*acos(1.0)
    elseif ((term1+term2).LT.-1.0) then
        realdistance = 6371.0*acos(-1.0)
    else
        realdistance = 6371.0*acos(term1+term2) ! 6371 is the radius of the earth in km (assuming spherical)
    endif

  end function

  real*4 function initialbearing(lat1,lon1,lat2,lon2)
    real*4,intent(in) :: lat1,lon1,lat2,lon2
    real*4            :: a,b,deltalon
    real*4            :: pi = 4*atan(1.0)

    deltalon = (lon2-lon1)*pi/180.0
    a = sin(deltalon)*cos(lat2*pi/180.0)
    b = cos(lat1*pi/180.0)*sin(lat2*pi/180.0) - cos(lat2*pi/180.0)*sin(lat1*pi/180.0)*cos(deltalon)
    initialbearing = atan2(a,b)*180.0/pi

  end function

END MODULE
