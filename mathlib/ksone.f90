SUBROUTINE ksone(data,n,func,d,prob)
INTEGER*4 n
REAL*8 d,data(n),func,prob
EXTERNAL func
!USES probks,sort
INTEGER*4 j
REAL*8 dt,en,ff,fn,fo,probks
call sort(n,data)
en=n
d=0.
fo=0.
do j=1,n
	fn=j/en
	ff=func(data(j))
	dt=max(abs(fo-ff),abs(fn-ff))
	if(dt.gt.d)d=dt
	fo=fn
enddo
en=sqrt(en)
prob=probks((en+0.12+0.11/en)*d)
return
END
