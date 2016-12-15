
subroutine boltzmanncalc(boltzmann_factors,currenta,currentb,initiala,initialb,earth,wind_m,wind_z,resources_filtered,a,b,c,kT,t)

  integer                             :: i,j,k,l
  real,intent(inout),dimension(3,3)   :: boltzmann_factors
  integer,intent(in)                  :: currenta,currentb,initiala,initialb
  real,intent(in)                     :: a,b,c,kT,wind_m,wind_z,t
  real,intent(in),dimension(720,1440) :: resources_filtered
  real,intent(in),dimension(720,1440) :: earth
  integer,dimension(2)                :: resources_shape
  integer,dimension(2)                :: state_index
  real,dimension(2)                   :: wind_vector
  integer,dimension(2)                :: displacement_vector
  real                                :: state_potential,wind_magnitude,displacement_vector_magnitude,wind_mag_sq,disp_mag_sq

  resources_shape = (/720,1440/)



!Calculate potentials in new possible states and convert to Boltzmann factors

  do i=1,3
    do j=1,3
      state_index = (/(currenta+i-2),MOD((currentb+j-2),resources_shape(2))/)
      if (i.EQ.2.AND.j.EQ.2) then
          CYCLE
      elseif (earth(state_index(1),state_index(2)).EQ.1) then
          CYCLE
      else
        state_potential = 0.0
        do k=1,resources_shape(1)
          do l=1,resources_shape(2)
            if (earth(k,l).EQ.0.AND.resources_filtered(k,l).GT.0.AND.k.NE.state_index(1).AND.l.NE.state_index(2)) then
              state_potential = state_potential + resources_filtered(k,l)/realdistance(k,l,state_index(1),state_index(2))
            endif
          enddo
        enddo
        wind_vector = (/wind_m,wind_z/) !In form [y,x] for ease of translation to np arrays.
        wind_mag_sq = DOT_PRODUCT(wind_vector,wind_vector)
        wind_magnitude = SQRT(wind_mag_sq)
        displacement_vector = (/i,j/)
        disp_mag_sq = DOT_PRODUCT(displacement_vector,displacement_vector)
        displacement_vector_magnitude = sqrt(disp_mag_sq)
        state_potential=state_potential+a*wind_magnitude*DOT_PRODUCT(wind_vector,displacement_vector)/displacement_vector_magnitude

        breeding_dist_dif = realdistance(initiala,initialb,state_index(1),state_index(2)) &
                            - realdistance(initiala,initialb,currenta,currentb)
        state_potential = state_potential - SIGN(1.0,breeding_dist_dif)*b*(ABS(breeding_dist_dif))**(t*c/8760.0)

        boltzmann_factors(i,j) = EXP(state_potential/kT)
      endif
    enddo
  enddo

end subroutine boltzmanncalc

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
