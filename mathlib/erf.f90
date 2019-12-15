FUNCTION erf(x)
REAL*8 erf,x
! USES gammp
REAL*8,EXTERNAL :: gammp
if(x.lt.0.)then
	erf=-gammp(0.5d0,x**2)
else
	erf=gammp(0.5d0,x**2)
endif
return
END
