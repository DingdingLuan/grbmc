FUNCTION probks(alam)
REAL*8 probks,alam,EPS1,EPS2
PARAMETER (EPS1=0.001, EPS2=1.e-8)
INTEGER*4 j
REAL*8 a2,fac,term,termbf
a2=-2.*alam**2
fac=2.
probks=0.
termbf=0.
do j=1,100
	term=fac*exp(a2*j**2)
	probks=probks+term
	if(abs(term).le.EPS1*termbf.or.abs(term).le.EPS2*probks)return
	fac=-fac
	termbf=abs(term)
end do
probks=1.
return
END
