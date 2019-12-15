SUBROUTINE qromb(func,a,b,ss)
INTEGER*4 JMAX,JMAXP,K,KM
REAL*8 a,b,func,ss,EPS
EXTERNAL func
PARAMETER (EPS=1.0d-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
!CU    USES polint,trapzd
INTEGER*4 j
REAL*8 dss,h(JMAXP),s(JMAXP)
h(1)=1.
do j=1,JMAX
	call trapzd(func,a,b,s(j),j)
	if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.0d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
	endif
	s(j+1)=s(j)
	h(j+1)=0.25*h(j)
end do
pause 'too many steps in qromb'
END
