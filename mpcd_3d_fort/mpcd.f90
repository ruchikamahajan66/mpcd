# This code is based on Multi-particle collision dynamics simulation technique that is used to model fluids 

# global simulation parameters defined

module mpcd_params
implicit none

integer, parameter:: tot_part=80000, max_part = 40
integer, parameter :: lx=20, ly=20, lz=20
real*8, parameter :: llx=dfloat(lx), lly=dfloat(ly), llz=dfloat(lz)
real*8 :: pos(3*tot_part), vel(3*tot_part)
real*8, parameter :: kb_temp=1.0d0, mass=1.0d0, dt=0.1d0, force=0.00348d0
real*8, parameter :: mass_inv = 1.0d0/mass
real*8, parameter :: vconst=dsqrt(12.0d0*kb_temp/mass), tot_partp = dfloat(tot_part)
real*8, parameter :: llxby2=llx/2.0d0
real*8, parameter :: llyby2=lly/2.0d0
real*8, parameter :: llzby2=llz/2.0d0
integer :: zzz=-4280145, zzzz=87383143, nofiter=5000, nofiter1=2000
integer :: part_no(lx*ly*lz), cell_part(max_part,lx*ly*lz)


end module mpcd_params




!#########################################################################

program mpcd
use mpcd_params
implicit none

integer:: i,j,step
real*8 :: random,zz1
real*8 :: avg_vel(20),avg_vel1(20),no_part(20)

common /comm1/ step

zz1=random(zzz)

call posinit
call velinit

  open(085,file='mpcenergy.dat',status='unknown',form='formatted')

do step=1,nofiter

  call streaming
  call celllist
  call collision

  write(*,*) step

end do


end program mpcd




!#########################################################################
subroutine posinit
use mpcd_params
implicit none

real*8  :: random
integer :: i

do i=1,tot_part

  pos(3*i-2) = random(zzzz)*llx
  pos(3*i-1) = random(zzzz)*lly
  pos(3*i)   = random(zzzz)*llz
    
end do

end subroutine posinit




!###########################################################################
subroutine velinit
use mpcd_params
implicit none

integer :: i,j
real*8 :: random
real*8 :: av_velx,av_vely,av_velz

av_velx=0.0d0
av_vely=0.0d0
av_velz=0.0d0

do i=1,tot_part

  vel(3*i-2) = vconst*(random(zzzz)-0.5)
  av_velx = av_velx + vel(3*i-2)

  vel(3*i-1) = vconst*(random(zzzz)-0.5)
  av_vely = av_vely + vel(3*i-1)

  vel(3*i) = vconst*(random(zzzz)-0.5)
  av_velz = av_velz + vel(3*i)

end do

av_velx = av_velx/tot_partp
av_vely = av_vely/tot_partp
av_velz = av_velz/tot_partp

do i=1,tot_part

  vel(3*i-2) = vel(3*i-2) - av_velx
  vel(3*i-1) = vel(3*i-1) - av_vely
  vel(3*i)   = vel(3*i)   - av_velz

end do


end subroutine velinit





!#########################################################################
subroutine celllist
use mpcd_params
implicit none

integer:: i,j,cell_no
real*8 :: r1,r2,r3,random
real*8 :: temp_pos(3*tot_part)

part_no = 0

r1 = random(zzzz) - 0.5d0
r2 = random(zzzz) - 0.5d0
r3 = random(zzzz) - 0.5d0

do i=1,tot_part

  temp_pos(3*i-2) = pos(3*i-2) - r1
  temp_pos(3*i-1) = pos(3*i-1) - r2
  temp_pos(3*i)   = pos(3*i)   - r3

  if (temp_pos(3*i-2).lt.0.0d0) then
    temp_pos(3*i-2) = temp_pos(3*i-2) + llx
  else if (temp_pos(3*i-2).gt.llx) then
    temp_pos(3*i-2) = temp_pos(3*i-2) - llx
  end if

  if (temp_pos(3*i-1).lt.0.0d0) then
    temp_pos(3*i-1) = temp_pos(3*i-1) + lly
  else if (temp_pos(3*i-1).gt.lly) then
    temp_pos(3*i-1) = temp_pos(3*i-1) - lly
  end if

  if (temp_pos(3*i).lt.0.0d0) then
    temp_pos(3*i) = temp_pos(3*i) + llz
  else if (temp_pos(3*i).gt.llz) then
    temp_pos(3*i) = temp_pos(3*i) - llz
  end if

end do 


do i=1,tot_part

  cell_no = 1 + int(temp_pos(3*i-2)) + lx*int(temp_pos(3*i-1)) + lx*ly*int(temp_pos(3*i))
  part_no(cell_no) = part_no(cell_no) + 1

  j = part_no(cell_no)
  cell_part(j,cell_no) = i

end do


end subroutine celllist






!##########################################################################
subroutine collision
use mpcd_params
implicit none

integer :: i,j,k,step
real*8 :: cell_vel1(lx*ly*lz), cell_vel2(lx*ly*lz), cell_vel3(lx*ly*lz)
real*8 :: rx,ry,rz,theta,random,phi,rho
real*8 :: rot11,rot12,rot13,rot21,rot22,rot23,rot31,rot32,rot33
real*8 :: del_vx,del_vy,del_vz,var,scale_fac
real*8 :: del_vx1(tot_part),del_vy1(tot_part),del_vz1(tot_part)

common /comm1/ step

cell_vel1 = 0.0d0
cell_vel2 = 0.0d0
cell_vel3 = 0.0d0

!SRD

do i=1,lx*ly*lz

  if (part_no(i).gt.1) then

    do j=1,part_no(i)
      k = cell_part(j,i)

      cell_vel1(i) = cell_vel1(i) + vel(3*k-2)
      cell_vel2(i) = cell_vel2(i) + vel(3*k-1)
      cell_vel3(i) = cell_vel3(i) + vel(3*k)

    end do

    cell_vel1(i) = cell_vel1(i)/dfloat(part_no(i))
    cell_vel2(i) = cell_vel2(i)/dfloat(part_no(i))
    cell_vel3(i) = cell_vel3(i)/dfloat(part_no(i))

    rho = 2.0d0*random(zzzz)-1.0d0
    phi = 4.0d0*asin(1.0d0)*random(zzzz)
    rx = cos(phi)*dsqrt(1-rho*rho)
    ry = sin(phi)*dsqrt(1-rho*rho)
    rz = rho
    theta = 2.0d0*asin(1.0d0)*130.0d0/180.0d0

    rot11 = (1.0d0-cos(theta))*rx*rx + cos(theta)
    rot12 = (1.0d0-cos(theta))*rx*ry - sin(theta)*rz
    rot13 = (1.0d0-cos(theta))*rx*rz + sin(theta)*ry
    rot21 = (1.0d0-cos(theta))*ry*rx + sin(theta)*rz
    rot22 = (1.0d0-cos(theta))*ry*ry + cos(theta)
    rot23 = (1.0d0-cos(theta))*ry*rz - sin(theta)*rx
    rot31 = (1.0d0-cos(theta))*rz*rx - sin(theta)*ry
    rot32 = (1.0d0-cos(theta))*rz*ry + sin(theta)*rx
    rot33 = (1.0d0-cos(theta))*rz*rz + cos(theta)


    do j=1,part_no(i)
      k = cell_part(j,i)

      del_vx = vel(3*k-2) - cell_vel1(i)
      del_vy = vel(3*k-1) - cell_vel2(i)
      del_vz = vel(3*k)   - cell_vel3(i)

      vel(3*k-2) = cell_vel1(i) + rot11*del_vx + rot12*del_vy + rot13*del_vz
      vel(3*k-1) = cell_vel2(i) + rot21*del_vx + rot22*del_vy + rot23*del_vz
      vel(3*k)   = cell_vel3(i) + rot31*del_vx + rot32*del_vy + rot33*del_vz

    end do
  end if
end do

!Thermostat

if (mod(step,25)==0) then

  do i=1,lx*ly*lz

    var = 0.0d0

    do j=1,part_no(i)
      k=cell_part(j,i)

      del_vx1(k) = vel(3*k-2) - cell_vel1(i)
      del_vy1(k) = vel(3*k-1) - cell_vel2(i)
      del_vz1(k) = vel(3*k)   - cell_vel3(i)

      var = var + del_vx1(k)*del_vx1(k) + del_vy1(k)*del_vy1(k) + del_vz1(k)*del_vz1(k)

    end do

    scale_fac = dsqrt(3.0d0*dfloat(part_no(i)-1)*kb_temp*mass_inv/var)

    do j=1,part_no(i)
      k=cell_part(j,i)

      vel(3*k-2) = cell_vel1(i) + del_vx1(k)*scale_fac
      vel(3*k-1) = cell_vel2(i) + del_vy1(k)*scale_fac
      vel(3*k)   = cell_vel3(i) + del_vz1(k)*scale_fac

    end do
  end do
end if


end subroutine collision






!###################################################################################
subroutine streaming
use mpcd_params
implicit none

integer:: i,j,k,step
real*8 :: energy,dt1,momentum
real*8 :: posx,posz,velx

common /comm1/ step

energy = 0.0d0
momentum = 0.0d0

do i=1,tot_part

  pos(3*i-2) = pos(3*i-2) + vel(3*i-2)*dt
  pos(3*i-1) = pos(3*i-1) + vel(3*i-1)*dt
  pos(3*i) = pos(3*i) + vel(3*i)*dt

  if (pos(3*i-2).lt.0.0d0) then
    pos(3*i-2) = pos(3*i-2) + llx
  else if (pos(3*i-2).gt.llx) then
    pos(3*i-2) = pos(3*i-2) - llx
  end if

  if (pos(3*i-1).lt.0.0d0) then
    pos(3*i-1) = pos(3*i-1) + lly
  else if (pos(3*i-1).gt.lly) then
    pos(3*i-1) = pos(3*i-1) - lly
  end if

  if (pos(3*i).lt.0.0d0) then
    pos(3*i) = pos(3*i) + llz
  else if (pos(3*i).gt.lly) then
    pos(3*i) = pos(3*i) - llz
  end if

  energy = energy + 0.5d0*mass*(vel(3*i-2)**2 + vel(3*i-1)**2 + vel(3*i)**2)/tot_partp
  momentum = momentum + mass*sqrt(vel(3*i-2)**2 + vel(3*i-1)**2 + vel(3*i)**2)/tot_partp
end do
  write(*,*) energy
  write(85,21) step,energy
  21 format(2g18.10)


end subroutine streaming







!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION random(IDUM)
implicit none

!  random returns a unifom random deviate on the interval [0,1]
! __________________________________________________________________
!
INTEGER :: IDUM
REAL*8 :: RAN2,random
integer,parameter :: IM1=2147483563,IM2=2147483399
integer,parameter :: IMM1=IM1-1,                                   &
IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
           NTAB=32
        integer,parameter :: NDIV=1+IMM1/NTAB
         real*8,parameter :: EPS=1.2e-7,RNMX=1.-EPS,AM=1./IM1
        INTEGER :: IDUM2,J,K,IV(NTAB),IY
        DATA IDUM2/123456789/, iv/NTAB*0/, iy/0/
        IF (IDUM.LE.0) THEN
          IDUM=MAX(-IDUM,1)
          IDUM2=IDUM
          DO  J=NTAB+8,1,-1
             K=IDUM/IQ1
             IDUM=IA1*(IDUM-K*IQ1)-K*IR1
             IF (IDUM.LT.0) IDUM=IDUM+IM1
             IF (J.LE.NTAB) IV(J)=IDUM
          end do
          IY=IV(1)
        ENDIF
        K=IDUM/IQ1
        IDUM=IA1*(IDUM-K*IQ1)-K*IR1
        IF (IDUM.LT.0) IDUM=IDUM+IM1
        K=IDUM2/IQ2
        IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
        IF (IDUM2.LT.0) IDUM2=IDUM2+IM2
        J=1+IY/NDIV
        IY=IV(J)-IDUM2
        IV(J)=IDUM
        IF(IY.LT.1)IY=IY+IMM1
        random=MIN(AM*IY,RNMX)

        END function random
