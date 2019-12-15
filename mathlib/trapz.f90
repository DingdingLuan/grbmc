subroutine trapz(x,v,n,F)
    implicit none;
 integer*4 i,n;
 real*8,dimension(n):: x,v;
 real*8:: F;
F=.5*(v(1)*(x(2)-x(1))+v(n)*(x(n)-x(n-1))&
&+DOT_PRODUCT(v(2:n-1),x(3:n)-x(1:n-2)));
   return;
end subroutine 
