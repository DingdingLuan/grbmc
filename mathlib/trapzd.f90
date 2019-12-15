SUBROUTINE trapzd(func,a,b,s,n)
INTEGER*4 n
REAL*8 a,b,s,func
EXTERNAL func
INTEGER*4 it,j
REAL*8 del,sum,tnm,x
if (n.eq.1) then
	s=0.5*(b-a)*(func(a)+func(b))
else
	it=2**(n-2)
	tnm=it
	del=(b-a)/tnm
	x=a+0.5*del
	sum=0.
	do j=1,it
		sum=sum+func(x)
		x=x+del
    enddo
	s=0.5*(s+(b-a)*sum/tnm)
endif
return
END
