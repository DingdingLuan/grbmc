module global
implicit none
real*8,save :: pi,A_egammao,sigma_ang,sigma_egammao,Ep,theta_jet
real*8,save :: A_ang,log_angc,log_egammaoc,theta_cut,s,sigma_gaussmd,A_obs
real*8,save :: A_t90,sigma_t90,log_t90c
real*8,parameter :: omegam=0.3,omegal=0.7,h=0.71
end module global

program main
use global
implicit none
integer*4 :: i,j,LM,k,ssb
integer*4,parameter :: nsi=500,negamma=500,num=100000,n=50,goalnum=50000
integer*4 :: nx,ny,selectnum
real*8 :: cos_arpha,doplr,Gamma,Gamma_c,b_slope,step_theta_rad,step_fai,betta,s_1st,s_2nd,eps0
real*8,dimension(nsi) :: Faiinte,GReisop,fai,theta_rad
real*8 :: cnt(n),cnt1(50),cnt2(20),cntt90(50)
real*8,dimension(nsi) :: log_ang,z,log_egammao,cal_obss,log_t90
real*8,dimension(num) :: z0,ang,logegammao,obs_ang,t90
real*8 :: ang1,ang2,tot,egammao1,egammao2,eiso1,eiso2,k1,k2,logfluence1,logfluence2,isostep,fluencestep,z_bin
real*8,dimension(nsi) :: tz,dndz,bin,f_metal
real*8,dimension(num) :: z_saved,pea_saved,Lp_saved,Eiso_saved,t90_saved
real*8 :: x,y,obs,psi,pEgammao,alpha,beta,pea,ss,dl_cm,dtdz,dvdz,gammln,gammp,z1,z2,dl_m,dl_mpc,sfr,obss
real*8,dimension(num) :: Eiso_const,Eiso_gauss,Eiso_pl,peak_const,peak_gauss,peak_pl
real*8,dimension(num) :: fluence_const,fluence_gauss,fluence_pl
real*8 :: pea1,pea2,peastep,Lp1,Lp2,Lpstep,z_1,z_2,z_step,t90step,Ep_prime
real*8 :: t901,t902,a1,a2,pea_erg,tn,threshold,eps,repo,Lp,fluence_sub
character(len=10) :: head
real*8 :: obsfluence(18),indexs(18),chsq,df,probb
external psi,pEgammao,band,bandE,ez,obss,tn,bplterm

!set parameters:
pi=3.14159265
!threshold=0.3
!threshold=1.27d-8
!threshold=5.89d-8   !flunce limited by obs data
log_egammaoc=52.62
sigma_egammao=0.64


  !for gauss model:
    !sigma_gaussmd=666
  !for broken pl model:
    s=2.


!read obsEiso pdf data:
open(77,file='obsresult/fluence.txt')
read(77,'(a10)') head
do ssb=1,18
	read(77,*)indexs(ssb),obsfluence(ssb)
	! write(*,*)obsfluence(ssb)
enddo


!Calculate the pdf of redshift:
z1=0.001;z2=10.0
call tzt(z1,z2,nsi,z,tz)
dndz=0.0
do i=1,nsi
  call qromb(ez,0.0d0,z(i),ss)
  dl_m=9.26d25/h*(1.0+z(i))*ss  ! in units M
  dl_mpc=dl_m*3.24d-23 ! in Mpc
  dtdz=9.78d9/h/((1.0+z(i))*sqrt(omegal+omegam*(1.0+z(i))**3)) ! in yr
  dvdz=4.0*pi*3.065948d-7*dl_mpc**2/(1.0+z(i))*dtdz ! in Mpc^3
  sfr=0.01*(1.0+z(i))**2.6/(1.0+((1.0+z(i))/3.2)**6.2) !SFR of Madau 2017
  ! if (z(i)<=1.0) then
  !   sfr=(1.0+z(i))**3.44
  ! else
  !   sfr=2.0**3.44
  !end if
	f_metal(i)=gammp(0.84d0,0.1**2*10.0**(0.3*z(i)))
	dndz(i)=f_metal(i)*sfr*dvdz/(1+z(i))
	!dndz(i)=sfr*dvdz/(1+z(i))
enddo
call trapz(z,dndz,nsi,tot)
! open(11,file="pdf_redshift.txt")
! do i=1,nsi
!   write(11,"(f7.3,2e11.3)") z(i),dndz(i)/tot,dndz(i)
! end do

! Generate redshift z0:
call random_seed()
!idum=-1257
i=1;cnt=0
do while (i<=num)
  call random_number(x)
  call random_number(y)
  x=(z2-z1)*x
  y=y*maxval(dndz/tot)
    call locate(z,nsi,x,j)
    if (y<=dndz(j)/tot) then
      z0(i)=x
      do k=1,50
        if (x>=(k-1)*0.2 .and. x<k*0.2) then
          cnt(k)=cnt(k)+1
        end if
      end do
      i=i+1
      !write(*,*)i
    end if
end do
! open(21,file='result/input_z_dis.txt')
! do k=1,50
!   write(21,*) (k-0.5)*0.2,cnt(k)/0.2/real(num)
! end do




!normalize log-normal functions:

!for jet-ang:
A_ang=1.0
sigma_ang=0.6
log_angc=-1.27
ang1=-2.5;ang2=-0.1
call qromb(psi,ang1,ang2,tot)
A_ang=1.0/tot
!for energy per solid angle:
A_egammao=1.0
egammao1=log_egammaoc-5.5;egammao2=log_egammaoc+5.5
call qromb(pEgammao,egammao1,egammao2,tot)
A_egammao=1./tot
!for t90:
A_t90=1.0
sigma_t90=0.526
log_t90c=1.69
t901=0.3;t902=3.
call qromb(tn,t901,t902,tot)
A_t90=1.0/tot
!for obs angle:
A_obs=1.0



!calculate pdf of jet_theta :
log_ang=(/((ang1+(ang2-ang1)/(nsi-1)*LM),LM=0,nsi-1)/)
!calculate pdf of Egamma0:
log_egammao=(/((egammao1+(egammao2-egammao1)/(negamma-1)*LM),LM=0,negamma-1)/)
!calculate pdf of t90:
log_t90=(/((t901+(t902-t901)/(nsi-1)*LM),LM=0,nsi-1)/)
!calculate pdf of obsangle:
cal_obss=(/((0+pi/2/(nsi-1)*LM),LM=0,nsi-1)/)


!generate jet-theta 
i=1;cnt=0
call random_seed()
do while (i<=num)
  call random_number(x)
  call random_number(y)
  x=(ang2-ang1)*x+ang1
  y=y*psi(log_angc)
  call locate(log_ang,nsi,x,j)
  if (y<=psi(log_ang(j))) then
    ang(i)=10**x
    do k=1,12
      if (x>=(ang1+(k-1)*0.2) .and. x<(ang1+k*0.2)) then
        cnt(k)=cnt(k)+1
      endif
    enddo
    i=i+1
    !write(*,*)i
  endif
enddo
open(33,file='result/input_angle.txt')
do k=1,12
  write(33,*) ang1+(k-0.5)*0.2,cnt(k)/0.2/real(num)
end do
ang=10**(-1.27)

!generate obs angle:
i=1;cnt=0
call random_seed()
do while (i<=num)
  call random_number(x)
  call random_number(y)
  x=pi/2*x
  call locate(cal_obss,nsi,x,j)
  if (y<=obss(cal_obss(j))) then
    obs_ang(i)=x
    do k=1,12
      if (x>=(0.+(k-1)*pi/2/12) .and. x<(0.+k*pi/2/12)) then
        cnt(k)=cnt(k)+1
      endif
    enddo
    i=i+1
    !write(*,*)i
  endif
enddo
! open(13,file='result/input_obsangle.txt')
! do k=1,12
!   write(13,*) 0.+(k-0.5)*pi/2/12,cnt(k)/0.2/real(num)
! enddo

!generate t90 :
i=1;cnt=0
call random_seed()
do while (i<=num)
  call random_number(x)
  call random_number(y)
  x=(t902-t901)*x+t901
	y=y*tn(log_t90c)
  call locate(log_t90,nsi,x,j)
  if (y<=tn(log_t90(j))) then
    t90(i)=10**x
    do k=1,12
      if (x>=(t901+(k-1)*0.225) .and. x<(t901+k*0.225)) then
        cnt(k)=cnt(k)+1
      endif
    enddo
    i=i+1
    !write(*,*)i
  endif
enddo
! open(32,file='result/input_t90.txt')
! do k=1,12
!   write(32,*) t901+(k-0.5)*0.225,cnt(k)/0.225/real(num)
! end do


!generate Egamma0:
i=1;cnt=0
call random_seed()
do while (i<=num)
  call random_number(x)
  call random_number(y)
  x=(egammao2-egammao1)*x+egammao1
  y=y*pEgammao(log_egammaoc)
	call locate(log_egammao,negamma,x,j)
  if (y<=pEgammao(log_egammao(j))) then
    logegammao(i)=x
    do k=1,10
      if (x>=(egammao1+(k-1)*0.5) .and. x<(egammao1+k*0.5)) then
        cnt(k)=cnt(k)+1
      endif
    enddo
    i=i+1
    !write(*,*)i
  endif
enddo
! open(34,file='result/input_egamma0.txt')
! do k=1,10
!   write(34,*) egammao1+(k-0.5)*0.5,cnt(k)/0.5/real(num)
! end do



Gamma_c=400.
b_slope=4.
!Calculations:
j=1;i=1
z_saved=0.;Lp_saved=0.;pea_saved=0.;Eiso_saved=0.
fluence_pl=0.
!do while (j.le.goalnum)
do while (i.le.num)
!   call random_seed()
!   call random_number(x)
!   x=num*x
!   i=int(x)
  !obs=obs_ang(i) !observation angle in rad
  !call qromb(ez,0.0d0,z0(i),ss)
  !dl_cm=9.26d27/h*(1.0+z0(i))*ss
  !theta_jet=ang(i)
  !call qromb(bplterm,theta_jet,pi/2,repo)      
  !eps0=10**logegammao(i)/(4*pi*((1-cos(theta_jet))+repo))
  !!Pl Model:
!!   if (obs<=ang(i)) then
!!     Eiso_pl(i)=eps*4*pi
!!   else
!!     Eiso_pl(i)=eps*(obs/ang(j))**(-s)*4*pi
!!   endif								!classical method
	!!Calculate Eiso in relative:
	!theta_rad=(/((0.+(2*pi)/(nsi-1)*LM),LM=0,nsi-1)/)
	!fai=(/((0.+(2*pi)/(nsi-1)*LM),LM=0,nsi-1)/)
	!step_theta_rad=theta_rad(2)-theta_rad(1)
	!step_fai=fai(2)-fai(1)
	!do nx=1,nsi-1
    !if (theta_rad(nx) .le. ang(i)) then
      !eps=eps0
      !Gamma=Gamma_c 
		!else
      !eps=eps0*(theta_rad(nx)/ang(i))**(-s)
      !Gamma=1+(Gamma_c-1)*(theta_rad(nx)/ang(i))**(-b_slope)
    !endif
    !betta=sqrt(1-1/Gamma**2)
		!!set fai repo
		!do ny=1,nsi-1
			!cos_arpha=cos(theta_rad(nx))*cos(obs)+sin(theta_rad(nx))*sin(fai(ny))*sin(obs)
			!doplr=1/(Gamma*(1-betta*cos_arpha))
			!Faiinte(ny)=doplr**3/Gamma*eps*sin(theta_rad(nx))
		!enddo
		!!cal 1_step fai integr
		!s_1st=0.
		!do ny=1,nsi-2
			!s_1st=s_1st+0.5*(Faiinte(ny)+Faiinte(ny+1))*step_fai
		!enddo
		!GReisop(nx)=s_1st
	!enddo
	!!cal 2_step theta_rad integr
	!s_2nd=0.
	!do nx=1,nsi-2
		!s_2nd=s_2nd+0.5*(GReisop(nx)+GReisop(nx+1))*step_theta_rad
	!enddo
	Eiso_pl(i)=10.0d51  ! Eiso in rest frame
	!!Ep=10**(-23.04+0.48*log10(Eiso_pl(i)))
	!Ep=100.
	!!call qromb(band,1.5d1,1.5d2,k1)
	!call qromb(band,1.0d0,1.0d3,k1)    ! from 1kev-1000kev
	!!call qromb(bandE,1.5d1*(1+z0(i)),1.5d2*(1+z0(i)),k1)
  !!call qromb(band,1.5d1,1.5d2,k1)
	!call qromb(bandE,1.0d0/(1.0+z0(i)),1.0d4/(1.0+z0(i)),k2)
	!!call qromb(bandE,1.0d0,1.0d4,k2)
	!k2=k2*1.60219d-9
	!Lp=2*Eiso_pl(i)/t90(i)*(1+z0(i))   !light_curve is triangle & t90 in observor frame
	!!Lp=2*Eiso_pl(i)/t90(i)   !light_curve is triangle & t90 in observor frame
	!!pea_erg=Lp*k1/(4*pi*dl_cm**2*k2)
	!!pea=2*Eiso_pl(i)*k1/(4*pi*dl_cm**2*t90(i)*k2)*(1+z0(i))  !light_curve is triangle
	!pea=Lp*k1/(4*pi*dl_cm**2*k2)  
	!!pea=2*Eiso_pl(i)*k1/(4*pi*dl_cm**2*t90(i)*k2)  !light_curve is triangle
	!Ep_prime=Ep/(1+z0(i))
	!threshold=1.514-0.8801*log10(Ep_prime)-0.06578*log10(Ep_prime)**2+0.01319*log10(Ep_prime)**3-0.02045*log10(Ep_prime)**4
	!threshold=10**threshold
	!call qromb(band,1.5d1,1.5d2,a1)
	!call qromb(bandE,1.5d1,1.5d2,a2)
  !!fluence_sub=Eiso_pl(i)*(1+z0(i))/(4*pi*dl_cm**2*k2/k1/1.60219d-9) !Regard Eiso as Earth_observec Nov.3
  !!call qromb(band,1.5d1/(1+z0(i)),1.5d2/(1+z0(i)),a1)
  !!call qromb(bandE,1.5d1/(1+z0(i)),1.5d2/(1+z0(i)),a2)
	!pea_erg=pea*a2/a1*1.60219d-9
         !!write(*,*)pea_erg/threshold
  !!if (fluence_sub .ge. threshold) then
	!!if (pea_erg .ge. threshold) then
	!if (pea .ge. threshold) then
    !!peak_pl(i)=pea
		!k1=a2
    !call qromb(bandE,1.5d1/(1+z0(i)),1.5d2/(1+z0(i)),k1)
    !call qromb(bandE,1.5d1,1.5d2,k1)
    !fluence_pl(j)=Eiso_pl(i)*(1+z0(i))/(4*pi*dl_cm**2*k2/k1/1.60219d-9)
		!fluence_pl(j)=Eiso_pl(i)*(1+z0(i))/(4*pi*dl_cm**2*k2/k1/1.60219d-9) 
		fluence_pl(j)=10.0d-8
		!call qromb(band,1.5d1,1.5d2,a1)
		!call qromb(bandE,1.5d1,1.5d2,a2)
		!pea_erg=a2/a1*pea
    z_saved(j)=z0(i)
    Eiso_saved(j)=Eiso_pl(i)
    !write(*,*)Eiso_saved(j)
		!Lp_saved(j)=pea_erg*4*pi*dl_cm**2*k2/k1/1.60219d-9
		Lp_saved(j)=Lp
		!write(*,*)Lp_saved(j)
    pea_saved(j)=pea_erg
    !pea_saved(j)=pea
    t90_saved(j)=t90(i)
    j=j+1    
  !else
  !   peak_pl(i)=0.
  !endif
  !write(*,*)real(j)/real(goalnum)*100,'%'
	write(*,*)real(i)/real(num)*100,'%'
  i=i+1
enddo
open(29,file='result/bpl/output_number.txt')
write(29,*)j-1
selectnum=j-1

!write(*,*)j-1
!calculate Eiso probability distribution
eiso1=49.
eiso2=55.
isostep=(eiso2-eiso1)/n
cnt=0;cnt1=0;cnt2=0;
do i=1,selectnum
  !if (peak_const(i)/=0.) then
    do k=1,n
        if (log10(Eiso_saved(i))>=(eiso1+(k-1)*isostep) .and. log10(Eiso_saved(i))<(eiso1+k*isostep)) then
        cnt(k)=cnt(k)+1
      endif
    enddo
  !endif
enddo
open(35,file='result/bpl/output_Eiso.txt')
cnt=cnt/sum(cnt)/isostep
write(35,*)'Eiso'
do k=1,n
  write(35,*)cnt(k)
enddo
open(38,file='result/bpl/Eiso_index.txt')
write(38,*)'index'
do k=1,n
  write(38,*) eiso1+(k-0.5)*isostep
enddo

!calculate t90 probability distribution
cntt90=0.
t90step=(t902-t901)/20
do i=1,selectnum
	do k=1,20
		if ((log10(t90_saved(i)))>=(t901+(k-1)*t90step) .and. log10(t90_saved(i))<(t901+k*t90step)) then
			cntt90(k)=cntt90(k)+1.
		endif
	enddo
enddo
cntt90=cntt90/sum(cntt90)/t90step
open(51,file='result/bpl/t90_index.txt')
open(52,file='result/bpl/output_t90.txt')
write(51,*)'index'
write(52,*)'t90'
do k=1,20
	write(51,*)t901+(k-0.5)*t90step
enddo
do k=1,20
	write(52,*)cntt90(k)
enddo

!calculate fluence probability distribution
logfluence1=-3.
logfluence2=-8.
fluencestep=(logfluence1-logfluence2)/25 !!! gai cheng 18 le 
cnt=0
do i=1,selectnum
  !if (fluence_pl(i)/=0.) then
		!do k=1,n
		do k=1,25    !!! gai cheng 25 le 
			if (log10(fluence_pl(i))<=(logfluence1-(k-1)*fluencestep) .and. log10(fluence_pl(i))>(logfluence1-k*fluencestep)) then
        cnt(k)=cnt(k)+1
			endif
		enddo
  !endif
enddo
cnt=cnt/sum(cnt)
open(39,file='result/bpl/output_fluence.txt')
! calculate the chi-square value for pdf_s
call chstwo(cnt,obsfluence,25,0,df,chsq,probb)
write(*,*)"Chi_Square Value is:",chsq
write(*,*)"q Value is:",probb
write(*,*)"Selected number is:",selectnum
fluencestep=(logfluence1-logfluence2)/n !!! gai cheng 18 le 
do i=1,selectnum
  !if (fluence_pl(i)/=0.) then
		!do k=1,n
		do k=1,n     
			if (log10(fluence_pl(i))<=(logfluence1-(k-1)*fluencestep) .and. log10(fluence_pl(i))>(logfluence1-k*fluencestep)) then
        cnt(k)=cnt(k)+1
			endif
		enddo
  !endif
enddo
cnt=cnt/sum(cnt)
cnt=cnt/fluencestep
write(39,*)'fluence'
do k=1,n
  write(39,*)cnt(k)
  !write(*,*) logfluence1-(k-0.5)*fluencestep,cnt(k)
enddo
open(42,file='result/bpl/fluence_index.txt')
write(42,*)'index'
do k=1,n
  write(42,*) logfluence1-(k-0.5)*fluencestep
enddo


!calculate Lp probability distribution:
Lp1=46.
Lp2=54.
Lpstep=(Lp2-Lp1)/n
cnt=0.;cnt(1)=selectnum
do k=1,n
    do i=1,selectnum-1
        if ((log10(Lp_saved(i)))>=(Lp1+(k-1)*Lpstep) .and. (log10(Lp_saved(i)))<(Lp1+k*Lpstep)) then
          cnt(k)=cnt(k)-1
        endif
    enddo
    cnt(k+1)=cnt(k)
enddo
open(36,file='result/bpl/output_Lp.txt')
cnt=cnt/real(selectnum)
write(36,*)'Lp'
do k=1,n
  write(36,*)cnt(k)
enddo
open(40,file='result/bpl/Lp_index.txt')
write(40,*)'index'
do k=1,n
  write(40,*)Lp1+(k-0.5)*Lpstep
  !write(*,*)Lp1+(k-0.5)*Lpstep,cnt(k)
enddo


!calculate pea_saved CDF probability distribution:
pea1=log10(threshold)
pea2=2.17609
peastep=(pea2-pea1)/n
cnt=0.;cnt(1)=selectnum
do k=1,n
    do i=1,selectnum-1
        if (log10(pea_saved(i))>=(pea1+(k-1)*peastep) .and. log10(pea_saved(i))<(pea1+k*peastep)) then
          cnt(k)=cnt(k)-1
        endif
    cnt(k+1)=cnt(k)
    enddo
enddo
open(41,file='result/bpl/output_pea.txt')
!write(*,*)'test'
!write(*,*)cnt(1),selectnum
cnt=cnt/real(selectnum)

write(41,*)'pea'
do k=1,n
  write(41,*)cnt(k)
enddo
open(43,file='result/bpl/pea_index.txt')
write(43,*)'index'
do k=1,n
  write(43,*)pea1+(k-0.5)*peastep
enddo



!calculate z_saved probability distribution:
!z_1=0.01
!z_2=7.
z_step=(z2-z1)/n
cnt=0;cnt(1)=selectnum
do k=1,n
    do i=1,selectnum-1
        if (z_saved(i)>=(z1+(k-1)*z_step) .and. z_saved(i)<(z1+k*z_step)) then
          cnt(k)=cnt(k)-1
        endif
    enddo
    cnt(k+1)=cnt(k)
enddo
open(44,file='result/bpl/output_z_cdf.txt')
cnt=cnt/real(selectnum)
write(44,*)'z'
do k=1,n
  write(44,*)cnt(k)
enddo
open(45,file='result/bpl/z_cdf_index.txt')
write(45,*)'index'
do k=1,n
  write(45,*)z1+(k-0.5)*z_step
enddo


!cnt=0
!do i=1,selectnum
    !do k=1,n
        !if (z_saved(i)>=(z1+(k-1)*z_step) .and. z_saved(i)<(z1+k*z_step)) then
          !cnt(k)=cnt(k)+1
        !endif
    !enddo
!enddo
!open(44,file='result/bpl/output_z.txt')
!cnt=cnt/sum(cnt)/z_step
!write(44,*)'z'
!do k=1,n
  !write(44,*)log10(cnt(k))
!enddo
!open(45,file='result/bpl/z_index.txt')
!write(45,*)'index'
!do k=1,n
  !write(45,*)z1+(k-0.5)*z_step
!enddo


! output the all kinds of variables_repo
open(48,file='result/bpl/z_repo.txt')
open(49,file='result/bpl/lp_repo.txt')
open(50,file='result/bpl/pea_repo.txt')
open(51,file='result/bpl/eiso_repo.txt')
open(52,file='result/bpl/t90_repo.txt')
open(53,file='result/bpl/fluence_repo.txt')
do i=1,selectnum
		write(48,*)z_saved(i)
		write(49,*)Lp_saved(i)
		write(50,*)pea_saved(i)*1.0d8
		write(51,*)Eiso_saved(i)
		write(52,*)t90_saved(i)
		write(53,*)fluence_pl(i)
enddo
open(54,file='result/bpl/output_parameters.txt')
write(54,*)'parameters'
write(54,*)chsq
write(54,*)log_egammaoc
write(54,*)sigma_egammao
write(54,*)Gamma_c
write(54,*)b_slope



end program main







function psi(log_ang)
! the intrinsic distribution of angle
use global, only : A_ang,sigma_ang,log_angc,pi
implicit none
real*8 psi,log_ang
psi=A_ang/sqrt(2.0*pi)/sigma_ang*exp(-(log_ang-log_angc)**2/(2.0*sigma_ang**2))
end function


function obss(ang_obss)
use global,only : pi,A_obs
implicit none
real*8 ang_obss,obss
obss=A_obs*sin(ang_obss)
end function

function bplterm(thetao)
use global, only : s,theta_jet
implicit none
real*8 bplterm,thetao
bplterm=(thetao/theta_jet)**(-s)*sin(thetao)
end function

function pEgammao(log_egammao)
use global, only:A_egammao,sigma_egammao,log_egammaoc,pi
implicit none
real*8 pEgammao,log_egammao
!pEgammao=A_egammao/(sqrt(2.0*pi)*sigma_egammao)*exp(-(log_egammao-log_egammaoc)**2/(2.0*sigma_egammao**2))
pEgammao=A_egammao*exp(-(log_egammao-log_egammaoc)**2/(2.0*sigma_egammao**2))
end function

function ez(zprime)
use global, only : omegal, omegam
implicit none
real*8 zprime,ez
ez=1.0/sqrt(omegal+omegam*(1+zprime)**3)
end function


function band(E)
!band funciton for lgrb:
use global, only : Ep
implicit none
real*8 band,E,a,b
a=-1.0
b=-2.25
if (E<=(a-b)*Ep/(2.0+a)) then
  band=(E/100.0)**a*exp(-E*(2.0+a)/Ep)
else
  band=((a-b)*Ep/(100.0*(2.0+a)))**(a-b)*exp(b-a)*(E/100.0)**b
endif
end function


function bandE(E)
!band funciton for lgrb:
use global, only : Ep
implicit none
real*8 bandE,E,a,b
a=-1.0
b=-2.25
if (E<=(a-b)*Ep/(2.0+a)) then
  bandE=(E/100.0)**a*exp(-E*(2.0+a)/Ep)*E
else
  bandE=((a-b)*Ep/(100.0*(2.0+a)))**(a-b)*exp(b-a)*(E/100.0)**b*E
endif
end function

function tn(log_t90)
use global, only : A_t90,sigma_t90,log_t90c,pi
implicit none
real*8 tn,log_t90
tn=A_t90/sqrt(2.0*pi)/sigma_t90*exp(-(log_t90-log_t90c)**2/(2.0*sigma_t90**2))
end function


include 'mathlib/sort.f90'
include 'mathlib/probks.f90'
include 'mathlib/kstwo.f90'
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
include 'mathlib/chstwo.f90'
include 'mathlib/gammq.f90'
