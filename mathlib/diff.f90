subroutine diff(x,n,v)
implicit none;
 integer ::i,n ;
 double precision,dimension(n)::x,v;
 do i=1,(n+2)/2
v(i)=x(i+1)-x(i);
v(n-i+1)=x(n-i+1)-x(n-i);
end do;
v(n/2+1)=0.5*(v(n/2)+v(n/2+2));
return;
end subroutine
