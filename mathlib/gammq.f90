FUNCTION gammq(a,x)
REAL*8 a,gammq,x
!CU    USES gcf,gser
REAL*8 gammcf,gamser,gln
if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
if(x.lt.a+1.)then
  call gser(gamser,a,x,gln)
  gammq=1.-gamser
else
  call gcf(gammcf,a,x,gln)
  gammq=gammcf
endif
return
END
