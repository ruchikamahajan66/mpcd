# This code is based on Multi-particle collision dynamics simulation technique in 2D that is used to model fluids 


module mpcd_params
implicit none

integer, parameter:: tot_part=8000, max_part = 100
integer, parameter :: lx=20, ly=20
real*8, parameter :: llx=dfloat(lx), lly=dfloat(ly)
real*8 :: pos(2*tot_part), vel(2*tot_part)
real*8, parameter :: kb_temp=1.0d0, mass=1.0d0, dt=0.1d0, force=0.1d0
real*8, parameter :: mass_inv = 1.0d0/mass
real*8, parameter :: vconst=dsqrt(12.0d0*kb_temp/mass), tot_partp = dfloat(tot_part)
integer :: zzz=-4280145, zzzz=87383143, nofiter=50000, nofiter1=2000
integer :: part_no(lx*ly), cell_part(max_part,lx*ly)


end module mpcd_params




!#########################################################################
program mpcd
use mpcd_params
implicit none

integer:: i,j,step,step1
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

  !write(*,*) step

end do


end program mpcd




!#########################################################################
subroutine posinit
use mpcd_params
implicit none

real*8  :: random
integer :: i

do i=1,tot_part

  pos(2*i-2) = random(zzzz)*llx
  pos(2*i-1) = random(zzzz)*lly
      
end do

end subroutine posinit




!###########################################################################
subroutine velinit
use mpcd_params
implicit none

integer :: i,j
real*8 :: random
real*8 :: av_velx,av_vely

av_velx=0.0d0
av_vely=0.0d0

do i=1,tot_part

  vel(2*i-2) = vconst*(random(zzzz)-0.5)
  av_velx = av_velx + vel(2*i-2)

  vel(2*i-1) = vconst*(random(zzzz)-0.5)
  av_vely = av_vely + vel(2*i-1)


end do

av_velx = av_velx/tot_partp
av_vely = av_vely/tot_partp

do i=1,tot_part

  vel(2*i-2) = vel(2*i-2) - av_velx
  vel(2*i-1) = vel(2*i-1) - av_vely

end do


end subroutine velinit





!#########################################################################
subroutine celllist
use mpcd_params
implicit none

integer:: i,j,cell_no
real*8 :: r1,r2,random
real*8 :: temp_pos(2*tot_part)

part_no = 0

r1 = random(zzzz) - 0.5d0
r2 = random(zzzz) - 0.5d0

do i=1,tot_part

  temp_pos(2*i-2) = pos(2*i-2) - r1
  temp_pos(2*i-1) = pos(2*i-1) - r2

  if (temp_pos(2*i-2).lt.0.0d0) then
    temp_pos(2*i-2) = temp_pos(2*i-2) + llx
  else if (temp_pos(2*i-2).gt.llx) then
    temp_pos(2*i-2) = temp_pos(2*i-2) - llx
  end if

  if (temp_pos(2*i-1).lt.0.0d0) then
    temp_pos(2*i-1) = temp_pos(2*i-1) + lly
  else if (temp_pos(2*i-1).gt.lly) then
    temp_pos(2*i-1) = temp_pos(2*i-1) - lly
  end if

end do 


do i=1,tot_part

  cell_no = 1 + int(temp_pos(2*i-2)) + lx*int(temp_pos(2*i-1))
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
real*8 :: cell_vel1(lx*ly), cell_vel2(lx*ly)
real*8 :: rot11,rot12,rot21,rot22,theta
real*8 :: del_vx,del_vy,random
real*8 :: del_vx1(tot_part),del_vy1(tot_part)

common /comm1/ step

   cell_vel1 = 0.0d0
  cell_vel2 = 0.0d0

!SRD

do i=1,lx*ly

  if (part_no(i).gt.1) then

    do j=1,part_no(i)
      k = cell_part(j,i)

      cell_vel1(i) = cell_vel1(i) + vel(2*k-2)
      cell_vel2(i) = cell_vel2(i) + vel(2*k-1)

    end do

    cell_vel1(i) = cell_vel1(i)/dfloat(part_no(i))
    cell_vel2(i) = cell_vel2(i)/dfloat(part_no(i))

    ! write(*,*) step,cell_vel1(i),cell_vel2(i)

    theta=2.0*3.14*random(zzzz);

           rot11 =  cos(theta);
           rot12 = -sin(theta);
           rot21 =  sin(theta);
           rot22 =  cos(theta);


    do j=1,part_no(i)
      k = cell_part(j,i)

      del_vx = vel(2*k-2) - cell_vel1(i)
      del_vy = vel(2*k-1) - cell_vel2(i)

      vel(2*k-2) = cell_vel1(i) + rot11*del_vx + rot12*del_vy 
      vel(2*k-1) = cell_vel2(i) + rot21*del_vx + rot22*del_vy 

    end do
  end if
end do

end subroutine collision






!###################################################################################
subroutine streaming
use mpcd_params
implicit none

integer:: i,j,k,step
real*8 :: energy,dt1,momentum
real*8 :: posx,posz,velx
real*8 :: dtafter

common /comm1/ step

energy = 0.0d0
momentum = 0.0d0

do i=1,tot_part

  pos(2*i-2) = pos(2*i-2) + vel(2*i-2)*dt!+(dt**2*force)/2.0*mass
  pos(2*i-1) = pos(2*i-1) + vel(2*i-1)*dt

 ! vel(2*i-2) = vel(2*i-2)  +(force*dt*mass_inv)

  if (pos(2*i-2).lt.0.0d0) then
    pos(2*i-2) = pos(2*i-2) + llx
  else if (pos(2*i-2).gt.llx) then
    pos(2*i-2) = pos(2*i-2) - llx
  end if


 if (pos(2*i-1).lt.0.0d0) then

		       dtafter = abs(pos(2*i-1))/abs(vel(2*i-1))

			pos(2*i-2)=pos(2*i-2)-2*dtafter*vel(2*i-2)
			pos(2*i-1)=pos(2*i-1)-2*dtafter*vel(2*i-1)
			
                        vel(2*i-2) = -1*vel(2*i-2)
			vel(2*i-1)=-1*vel(2*i-1)		

 else if (pos(2*i-1).gt.lly) then
		
		  dtafter=abs(pos(2*i-1)-lly)/abs(vel(2*i-1))
                     
			pos(2*i-2)=pos(2*i-2)-2*dtafter*vel(2*i-2)
			pos(2*i-1)=pos(2*i-1)-2*dtafter*vel(2*i-1)

                        vel(2*i-2) = -1*vel(2*i-2)
			vel(2*i-1)=-1*vel(2*i-1)
endif

  energy = energy + 0.5d0*mass*(vel(2*i-2)**2 + vel(2*i-1)**2) /tot_partp
  momentum = momentum + mass*sqrt(vel(2*i-2)**2 + vel(2*i-1)**2) /tot_partp
end do
  write(*,*) step,energy,momentum
  write(85,21) step,energy,momentum
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
