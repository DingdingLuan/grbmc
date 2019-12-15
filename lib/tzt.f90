subroutine tzt(z1,z2,nz,z,t)
   implicit none;
 integer*4::s,nz,Lz;
 integer*4,parameter :: nz1=50001
 real*8,dimension(nz)::z,t;
 real*8,dimension(nz1)::zn,tz;
 real*8::z0,z1,z2,zz,omegam,omegal;
omegam=0.266
omegal=0.734
z0=1000.0
z=(/((z1+(z2-z1)/(nz-1)*Lz),Lz=0,nz-1)/);
!t=-1.240191545d10*z+4.703083872d7+68927.1838&
!&*sqrt(3.240302366d10*z**2-1.868729701d9*z+3.940756737d10);
do s=1,nz
	zz=z(s)
	zn=(/((zz+(z0-zz)/(nz1-1)*Lz),Lz=0,nz1-1)/);
	tz=9.78d9/0.71/(1+zn)/sqrt(omegal+omegam*(1+zn)**3)
	call trapz(zn,tz,nz1,t(s))
enddo
   return;
end subroutine
