! Christopher Revell, University of Cambridge, 2016

subroutine boltzmanncalc(boltzmann_factors,currenta,currentb,destinationa,destinationb,elevation,earth,a,b,kT)

  integer*8                             :: i,j,k,l
  real*8,intent(inout),dimension(3,3)   :: boltzmann_factors
  real*8,intent(in),dimension(3,3)      :: elevation
  integer*8,intent(in)                  :: currenta,currentb,destinationa,destinationb
  real*8,intent(in)                     :: a,b,kT
  real*8,intent(in),dimension(720,1440) :: earth
  integer*8,dimension(2)                :: earthshape
  integer*8,dimension(2)                :: state_index
  real*8                                :: state_potential
  real*8                                :: min_potential
  logical,dimension(3,3)                :: boltzmann_mask

  earthshape = (/720,1440/)

  boltzmann_mask(:,:) = .TRUE.

!Calculate potentials in new possible states and convert to Boltzmann factors

  do i=1,3
    do j=1,3
      state_index = (/(currenta+i-2),MOD((currentb+j-2),earthshape(2))/)
      if (i.EQ.2.AND.j.EQ.2) then
          boltzmann_mask(i,j) = .FALSE.
      else
        state_potential = elevation(state_index(1),state_index(2))
        state_potential = state_potential + a*ABS(realdistance(currenta,currentb,destinationa,destinationb))
        do k=1,earthshape(1)
          do l=1,earthshape(2)
            if (k.NE.state_index(1).AND.l.NE.state_index(2)) then
              state_potential = state_potential + b*earth(k,l)/(realdistance(k,l,state_index(1),state_index(2)))
            endif
          enddo
        enddo
        boltzmann_factors(i,j) = state_potential/kT
      endif
    enddo
  enddo

  min_potential = MINVAL(boltzmann_factors,MASK=boltzmann_mask)

  do i=1,3
    do j=1,3
      if (boltzmann_mask(i,j)) then
        boltzmann_factors(i,j) = EXP(boltzmann_factors(i,j)-min_potential)
      else
        boltzmann_factors(i,j) = 0
      endif
    enddo
  enddo

end subroutine boltzmanncalc

!Assumes 720x1440 environment arrays
function realdistance(a1,a2,b1,b2)

  integer*8,intent(in) :: a1,a2,b1,b2
  real*8               :: d_latlong,delta_long,term1,term2
  real*8,dimension(2)  :: latlonga,latlongb
  real*8,parameter     :: pi  = 4*atan (1.0)

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
