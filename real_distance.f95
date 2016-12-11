!Assumes 720x1440 environment arrays
function realdistance(a1,a2,b1,b2)

  integer,intent(in) :: a1,a2,b1,b2
!  real, intent(out)  :: realdistance
  real               :: d_latlong,delta_long,term1,term2
  real,dimension(2)  :: latlonga,latlongb
  real,parameter     :: pi  = 4*atan (1.0)

  d_latlong = 180.0/720.0

  latlonga   = (/(720.0/2.0-a1-0.5)*d_latlong, (a2+0.5-1440.0/2.0)*d_latlong/)
  latlongb   = (/(720.0/2.0-b1-0.5)*d_latlong, (b2+0.5-1440.0/2.0)*d_latlong/)
  latlonga   = latlonga*pi/180.0
  latlongb   = latlongb*pi/180.0
  delta_long = latlongb(2)-latlonga(2)
  term1      = sin(latlonga(1))*sin(latlongb(1))
  term2      = cos(latlonga(1))*cos(latlongb(1))*cos(delta_long)
  !When the bird doesn't move, and latlonga = latlongb, small rounding errors can lead to taking the arccos of a number a tiny bit higher than 1, eg 1.0000000000000002, so it's safer to set the realdistance equal to 0 in this case rather than doing the full calculation
  if (a1.EQ.b1.AND.a2.EQ.b2) then
      realdistance = 0.0
  elseif ((term1+term2).GT.1.0) then
      realdistance = 6371.0*acos(1.0)
  elseif ((term1+term2).LT.-1.0) then
      realdistance = 6371.0*acos(-1.0)
  else
      realdistance = 6371.0*acos(term1+term2) ! 6371 is the radius of the earth in km (assuming spherical)
  endif

end function realdistance
