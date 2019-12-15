subroutine line(a,b,n,v)
implicit none;
 integer ::i,n ;
 real*8::a,b;
 real*8,dimension(n)::v;
v=(/((a+(b-a)/(n-1)*i),i=0,n-1)/);
return;
end subroutine
