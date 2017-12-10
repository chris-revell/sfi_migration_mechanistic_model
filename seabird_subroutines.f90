! Christopher Revell, University of Cambridge, 2016

subroutine boltzmanncalc(currentlat,currentlon,resources_filtered,a,kT)

  integer*8                             :: i,j,k,l
  real*8,intent(inout),dimension(3,3)   :: boltzmann_factors
  integer*8,intent(in)                  :: currentlat,currentlon,initiala,initialb
  real*8,intent(in)                     :: a,b,c,kT,wind_m,wind_z,t
  real*8,intent(in),dimension(720,1440) :: resources_filtered
  real*8,intent(in),dimension(720,1440) :: earth
  integer*8,dimension(2)                :: resources_shape
  integer*8,dimension(2)                :: state_index
  real*8,dimension(2)                   :: wind_vector
  integer*8,dimension(2)                :: displacement_vector
  real*8                                :: state_potential,wind_magnitude,displacement_vector_magnitude,wind_mag_sq,disp_mag_sq
  real*8                                :: min_potential

  integer*8,dimension(2)                :: force

  resources_shape = (/720,1440/)

!Calculate potentials in new possible states and convert to Boltzmann factors

  force = (/0,0/)

  do i=1,resources_shape(1)
    do j=1,resources_shape(2)
      lat = (resources_shape[0]/2.0-i-0.5)*d_latlong
      lon = (j+0.5-resources_shape[1]/2.0)*d_latlong
      
      forcemagnitude = resources_filtered(i,j)/realdistance(lat,lon,currentlat,currentlon)




  !do i=1,3
  !  do j=1,3
  !    state_index = (/(currentlat+i-2),MOD((currentlon+j-2),resources_shape(2))/)
  !    if (i.EQ.2.AND.j.EQ.2) then
  !        boltzmann_mask(i,j) = .FALSE.
  !    elseif (earth(state_index(1),state_index(2)).EQ.1) then
  !        boltzmann_mask(i,j) = .FALSE.
  !    else
  !      state_potential = 0.0
  !      do k=1,resources_shape(1)
  !        do l=1,resources_shape(2)
  !          if (earth(k,l).EQ.0.AND.resources_filtered(k,l).GT.0.AND.k.NE.state_index(1).AND.l.NE.state_index(2)) then
  !            state_potential = state_potential + resources_filtered(k,l)/(realdistance(k,l,state_index(1),state_index(2)))
  !          endif
  !        enddo
  !      enddo
  !      wind_vector = (/wind_m,wind_z/) !In form [y,x] for ease of translation to np arrays.
  !      wind_mag_sq = DOT_PRODUCT(wind_vector,wind_vector)
  !      wind_magnitude = SQRT(wind_mag_sq)
  !      displacement_vector = (/i,j/)
  !      disp_mag_sq = DOT_PRODUCT(displacement_vector,displacement_vector)
  !      displacement_vector_magnitude = sqrt(disp_mag_sq)
  !      state_potential=state_potential+a*wind_magnitude*DOT_PRODUCT(wind_vector,displacement_vector)/displacement_vector_magnitude
!
  !      breeding_dist_dif = realdistance(initiala,initialb,state_index(1),state_index(2)) &
  !                          - realdistance(initiala,initialb,currentlat,currentlon)
  !      state_potential = state_potential - SIGN(1.0,breeding_dist_dif)*b*(ABS(breeding_dist_dif))**(t*c/8760.0)
  !      boltzmann_factors(i,j) = state_potential/kT
  !    endif
  !  enddo
  !enddo
!
  !min_potential = MINVAL(boltzmann_factors,MASK=boltzmann_mask)
!
  !do i=1,3
  !  do j=1,3
  !    if (boltzmann_mask(i,j)) then
  !      boltzmann_factors(i,j) = EXP(boltzmann_factors(i,j)-min_potential)
  !    else
  !      boltzmann_factors(i,j) = 0
  !    endif
  !  enddo
  !enddo

end subroutine boltzmanncalc

!Assumes 720x1440 environment arrays
function realdistance(lata,lona,latb,lonb)

  real*8,intent(in) :: lata,lona,latb,lonb
  real*8               :: delta_long,term1,term2
  real*8,parameter     :: pi  = 4*atan (1.0)

  delta_long = lonb-lona
  term1      = sin(lata*pi/180.0)*sin(latb*pi/180.0)
  term2      = cos(lata*pi/180.0)*cos(latb*pi/180.0)*cos(delta_long*pi/180.0)
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
