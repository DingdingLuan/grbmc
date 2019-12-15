SUBROUTINE gser(gamser,a,x,gln)
INTEGER*4 ITMAX
REAL*8 a,gamser,gln,x,EPS
PARAMETER (ITMAX=100,EPS=3.d-7)
!CU    USES gammln
INTEGER*4 n
REAL*8 ap,del,sum,gammln
EXTERNAL gammln
gln=gammln(a)
if(x.le.0.)then
	if(x.lt.0.)pause 'x < 0 in gser'
	gamser=0.
	return
endif
ap=a
sum=1./a
del=sum
do n=1,ITMAX
	ap=ap+1.
	del=del*x/ap
	sum=sum+del
	if(abs(del).lt.abs(sum)*EPS)goto 1
end do
pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
return
END
