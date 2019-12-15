SUBROUTINE kstwo(data1,n1,data2,n2,d,prob)
INTEGER*4 n1,n2
REAL*8 d,prob,data1(n1),data2(n2)
! USES probks,sort
INTEGER*4 j1,j2
REAL*8 d1,d2,dt,en1,en2,en,fn1,fn2,probks
call sort(n1,data1)
call sort(n2,data2)
en1=n1
en2=n2
j1=1
j2=1
fn1=0.
fn2=0.
d=0.
do while(j1.le.n1.and.j2.le.n2)
	d1=data1(j1)
	d2=data2(j2)
	if(d1.le.d2)then
		fn1=j1/en1
          j1=j1+1
  	endif
	if(d2.le.d1)then
          fn2=j2/en2
          j2=j2+1
	endif
	dt=abs(fn2-fn1)
	if(dt.gt.d)d=dt
enddo
en=sqrt(en1*en2/(en1+en2))
prob=probks((en+0.12+0.11/en)*d)
return
END
