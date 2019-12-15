module global
implicit none
real*8,save :: pi,A_egammao,sigma_ang,sigma_egammao
real*8,save :: A_ang,log_angc,log_egammaoc,theta_cut,s,sigma_gaussmd,A_obs
real*8,save :: A_t90,sigma_t90,log_t90c,Ep
real*8,parameter :: omegam=0.3,omegal=0.7,h=0.71,a=-1.,b=-2.25
end module global


program main
use global, only : Ep,pi,h
implicit none
integer*4,parameter :: n=167,ns=367,num=14,nums=18,numz=100
integer*4 :: i,k
real*8,dimension(num) :: cnte,cntp,cntt90
real*8,dimension(numz) :: cntz,cntl
real*8,dimension(nums) :: cnts
real*8,dimension(n) :: z,t90,expoure,fluence_z,fluxph,fluxerg,eiso,Lp_saved,Epeak
real*8,dimension(ns) :: s,zs
character(len=10),dimension(n) :: grbname
real*8 :: z_1,z_2,Lp1,Lp2,logfluence1,logfluence2,eiso1,eiso2,pea1,pea2,t901,t902,threshold
real*8 :: z_step,z_step1,z_step2,Lpstep,fluencestep,isostep,peastep,t90step,zsum
real*8 :: a1,a2,dl_cm,ss,kc
character(len=100) :: head
external ez,band,bandE

!set threshold
! threshold=2.6
pi=3.141592653


open(14,file='data/fluence_z.txt')
read(14,*)head
do i=1,ns
	read(14,*)zs(i),s(i)
enddo

!open(15,file='data/swift_lgrb.txt')
open(15,file='data/cpl_pl.txt')
read(15,*)head
do i=1,n
  read(15,*)grbname(i),z(i),t90(i),expoure(i),fluence_z(i),fluxph(i),fluxerg(i),Epeak(i)
enddo

! do i=1,n
!   write(*,'(a10,f6.3,2f10.3,e12.3,f8.3,e12.3)')grbname(i),z(i),t90(i),expoure(i),s(i),fluxph(i),fluxerg(i)
! enddo

!calculate Eiso & Lp
do i=1,n
    call qromb(ez,0.0d0,z(i),ss)
    dl_cm=9.26d27/h*(1.0+z(i))*ss
		Ep=Epeak(i)
    call qromb(bandE,1.5d1,1.5d2,a1)
    call qromb(bandE,1.0d0/(1.0+z(i)),1.0d4/(1.0+z(i)),a2)
    kc=a2/a1
    eiso(i)=4*pi*dl_cm**2*fluence_z(i)*kc/(1+z(i))
		!Lp_saved(i)=4*pi*dl_cm**2*fluxerg(i)*kc
		Lp_saved(i)=2*eiso(i)/t90(i)*(1+z(i))
enddo

! do i=1,n
!   write(*,*)eiso(i)
! enddo

open(20,file='obsresult/z.txt')
open(21,file='obsresult/lp.txt')
open(22,file='obsresult/fluence_pdf.txt')
open(30,file='obsresult/fluence.txt')
open(23,file='obsresult/eiso.txt')
open(24,file='obsresult/peakflux.txt')
open(25,file='obsresult/t90.txt')

! write(20,'(a10)')'index,z'
! write(21,'(a10)')'index,lp'
! write(22,'(a10)')'index,fluence'
! write(23,'(a10)')'index,eiso'
! write(24,'(a10)')'index,peakflux'


cntz=0.;cntl=0.;cnts=0.;cntp=0.;cnte=0.;cntt90=0.
!calculate z_cdf ..
z_1=0.001
z_2=10.
z_step=(z_2-z_1)/numz
cntz(1)=n
do k=1,numz
	do i=1,n-1
		if (z(i)>=(z_1+(k-1)*z_step) .and. z(i)<(z_1+k*z_step)) then
			cntz(k)=cntz(k)-1.
		endif
	enddo
	cntz(k+1)=cntz(k)
enddo
cntz=cntz/real(n)
open(37,file='obsresult/z_cdf.txt')
write(37,*)'index,z'
do k=1,numz	
	write(37,*)z_1+(k-0.5)*z_step,',',cntz(k)
enddo

cntz=0
cntz(1)=ns
do k=1,numz
	do i=1,ns-1
		if (zs(i)>=(z_1+(k-1)*z_step) .and. zs(i)<(z_1+k*z_step)) then
			cntz(k)=cntz(k)-1.
		endif
	enddo
	cntz(k+1)=cntz(k)
enddo
cntz=cntz/real(ns)
open(38,file='obsresult/sz.txt')
write(38,*)'index,z'
do k=1,numz	
	write(38,*)z_1+(k-0.5)*z_step,',',cntz(k)
enddo

Lp1=46.
Lp2=54.
Lpstep=(Lp2-Lp1)/numz
do k=1,numz
	cntl(1)=n
  do i=1,n-1
    if ((log10(Lp_saved(i)))>=(Lp1+(k-1)*Lpstep) .and. log10(Lp_saved(i))<(Lp1+k*Lpstep)) then
      cntl(k)=cntl(k)-1.
    endif
  enddo
  cntl(k+1)=cntl(k)
enddo
cntl=cntl/n
write(21,*)'asd,asd'
do k=1,numz
	write(21,*)Lp1+(k-0.5)*Lpstep,',',cntl(k)
enddo

!calculate fluence_pdf
logfluence1=3.
logfluence2=8.
fluencestep=(logfluence2-logfluence1)/nums
do k=1,nums
  do i=1,ns
    if (log10(s(i))<=(-(logfluence1+(k-1)*fluencestep)) .and. log10(s(i))>(-(logfluence1+k*fluencestep))) then
      cnts(k)=cnts(k)+1.
    endif
  enddo
enddo
cnts=cnts/sum(cnts)
write(22,*)'index,s,bin'
write(30,*)'index,s'
do k=1,nums
	write(22,*)logfluence1+(k-0.5)*fluencestep,',',cnts(k)/fluencestep,',',fluencestep/2
	write(30,*)logfluence1+(k-0.5)*fluencestep,cnts(k)
enddo

!calculate all kinds of pdf ..
do k=1,num
  !z cdf
	!z_1=0.001
	!z_2=10.
	!z_step=(z_2-z_1)/num
	!cntz(1)=n
	!do i=1,n-1
		!if (z(i)>=(z_1+(k-1)*z_step) .and. z(i)<(z_1+k*z_step)) then
			!cntz(k)=cntz(k)-1.
		!endif
	!enddo
  ! cntz(k+1)=cntz(k)
  ! z_1=0.0001
  ! z_2=10.
  ! z_step=(z_2-z_1)/num
  ! do i=1,n
  !   if (z(i)>=(z_1+(k-1)*z_step) .and. z(i)<(z_1+k*z_step)) then
  !     cntz(k)=cntz(k)+1
  !   endif
  ! enddo

  !lp cdf
  !Lp1=46.
  !Lp2=54.
  !Lpstep=(Lp2-Lp1)/num
  !cntl(1)=n
  !do i=1,n-1
    !if ((log10(Lp_saved(i)))>=(Lp1+(k-1)*Lpstep) .and. log10(Lp_saved(i))<(Lp1+k*Lpstep)) then
      !cntl(k)=cntl(k)-1.
    !endif
  !enddo
  !cntl(k+1)=cntl(k)

  !fluence
  !logfluence1=3.
  !logfluence2=8.
  !fluencestep=(logfluence2-logfluence1)/nums
  !do i=1,n
    !if (log10(s(i))<=(-(logfluence1+(k-1)*fluencestep)) .and. log10(s(i))>(-(logfluence1+k*fluencestep))) then
      !cnts(k)=cnts(k)+1.
    !endif
  !enddo
  !t90
  t901=0.3
  t902=3.0
  t90step=(t902-t901)/num
  do i=1,n
    if ((log10(t90(i)))>=(t901+(k-1)*t90step) .and. log10(t90(i))<(t901+k*t90step)) then
      cntt90(k)=cntt90(k)+1.
    endif
  enddo
  !Eiso
  ! eiso1=49. 
  ! eiso2=55.
  ! isostep=(eiso2-eiso1)/num
  ! do i=1,n
  !   if (log10(eiso(i))>=(eiso1+(k-1)*isostep) .and. log10(eiso(i))<(eiso1+k*isostep)) then
  !     cnte(k)=cnte(k)+1.
  !   endif
  ! enddo
  !cdf of peak
  pea1=log10(0.2)
  pea2=2.17609
  peastep=(pea2-pea1)/num
  cntp(1)=n
  do i=1,n-1
    if (log10(fluxph(i))>=(pea1+(k-1)*peastep) .and. log10(fluxph(i))<(pea1+k*peastep)) then
      cntp(k)=cntp(k)-1
    endif
  enddo
  cntp(k+1)=cntp(k)
enddo

do k=1,num
  eiso1=48. 
  eiso2=56.
  isostep=(eiso2-eiso1)/num
  do i=1,n
    if (log10(eiso(i))>=(eiso1+(k-1)*isostep) .and. log10(eiso(i))<(eiso1+k*isostep)) then
      cnte(k)=cnte(k)+1.
    endif
  enddo
enddo

! 0.01-4. of z
do k=1,10
 z_1=0.0001
  z_2=4.
  z_step=0.4
  do i=1,n
    if (z(i)>=(z_1+(k-1)*z_step) .and. z(i)<(z_1+k*z_step)) then
      cntz(k)=cntz(k)+1
    endif
  enddo
enddo
!4.-10. of z
z_step1=0.4+0.4/2
z_step2=0.4*2+0.4/2
do i=1,n
  if (z(i)>=(z_2) .and. z(i)<(z_2+z_step1)) then
    cntz(11)=cntz(11)+1
  endif
enddo
do i=1,n
  if (z(i)>=(z_2+z_step1) .and. z(i)<(z_2+z_step1+z_step2)) then
      cntz(12)=cntz(12)+1
  endif
enddo


! export pdf result:
! cntz=cntz/sum(cntz)/z_step
zsum=sum(cntz)
do k=1,10
  cntz(k)=cntz(k)/zsum/z_step
enddo
cntz(11)=cntz(11)/zsum/z_step1
cntz(12)=cntz(12)/zsum/z_step2
!cntl=cntl/n
!cnts=cnts/sum(cnts)
cnte=cnte/sum(cnte)/isostep
cntp=cntp/n
cntt90=cntt90/sum(cntt90)/t90step

! do i=1,num
!   write(*,*)cnts(i)
! enddo

write(20,*)'index,z,bin'
!write(21,*)'asd,asd'
!write(22,*)'index,s,bin'
!write(30,*)'index,s,bin'
write(23,*)'asd,asd'
write(24,*)'asd,asd'
write(25,*)'asd,asd'

do k=1,num
    ! write(20,*)z_1+(k-0.5)*z_step,',',cntz(k)
    !write(21,*)Lp1+(k-0.5)*Lpstep,',',cntl(k)
    !write(22,*)logfluence1+(k-0.5)*fluencestep,',',cnts(k)/fluencestep,',',fluencestep/2
    !write(30,*)logfluence1+(k-0.5)*fluencestep,cnts(k)
		write(23,*)eiso1+(k-0.5)*isostep,',',cnte(k)
    write(24,*)pea1+(k-0.5)*peastep,',',cntp(k)
    write(25,*)t901+(k-0.5)*t90step,',',cntt90(k)
enddo

!do k=1,num
  !write(23,*)eiso1+(k-0.5)*isostep,',',cnte(k)
!enddo
do k=1,10
    write(20,*)z_1+(k-0.5)*z_step,log10(cntz(k))
enddo
write(20,*)4.4,log10(cntz(11))
write(20,*)5.4,log10(cntz(12))


open(31,file='obsresult/lp_repo.txt')
open(32,file='obsresult/z_repo.txt')
open(33,file='obsresult/pea_repo.txt')
open(34,file='obsresult/eiso_repo.txt')
open(35,file='obsresult/t90_repo.txt')
open(36,file='obsresult/fluence_repo.txt')
open(39,file='obsresult/fluence_zrepo.txt')
write(31,*)'lp'
write(32,*)'z'
write(33,*)'peak_flux_erg'
write(34,*)'eiso'
write(35,*)'t90'
write(36,*)'fluence'
write(39,*)'fluence'
do i=1,n
  write(31,*)Lp_saved(i)
  write(32,*)z(i)
  write(33,*)fluxerg(i)*1.0d8
  write(34,*)eiso(i)
  write(35,*)t90(i)
	write(39,*)s(i)
enddo
do i=1,ns
	write(36,*)s(i)
enddo


end program main

function ez(zprime)
use global, only : omegal, omegam
implicit none
real*8 zprime,ez
ez=1.0/sqrt(omegal+omegam*(1+zprime)**3)
end function



function bandE(E)
use global, only : Ep,a,b
implicit none
real*8 bandE,E
if (E<=(a-b)*Ep/(2.0+a)) then
  bandE=(E/100.0)**a*exp(-E*(2.0+a)/Ep)*E
else
  bandE=((a-b)*Ep/(100.0*(2.0+a)))**(a-b)*exp(b-a)*(E/100.0)**b*E
endif
end function



include 'mathlib/trapzd.f90'
include 'mathlib/trapz.f90'
include 'mathlib/qromb.f90'
include 'mathlib/gammp.f90'
include 'mathlib/gammln.f90'
include 'mathlib/erf.f90'
include 'lib/tzt.f90'
include 'mathlib/polint.f90'
include 'mathlib/gcf.f90'
include 'mathlib/gser.f90'
include 'lib/locate.f90'
include 'mathlib/ran1.f90'
