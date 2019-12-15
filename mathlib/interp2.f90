subroutine interp2(x,y,v,n1,n2,a,b,NZ)
   implicit none
   integer*4::n1,n2;
   integer*4,dimension(1)::sx,sy;
   real*8(n2,n1)::v;
   real*8(n1):: x;
   real*8(n2):: y;
   real*8(1,1):: NZ1,v1,v2,v3,v4,v5;
   real*8(1):: x1,x2,x3,y1,y2,y3
   real*8::a,b,NZ
   sx=minloc(abs(x-a));sy=minloc(abs(y-b));
   y1=y(sy);x1=x(sx);
   y2=y(sy-1);x2=x(sx-1);y3=y(sy+1);x3=x(sx+1);
   v1=v(sy,sx);v2=v(sy,sx-1);v3=v(sy-1,sx);v4=v(sy+1,sx);
   v5=v(sy,sx+1);
   if((y1(1)>b).AND.(x1(1)>a))then
      NZ1=v1(1,1)-(v1(1,1)-v2(1,1))/(x1(1)-x2(1))*(x1(1)-a)&
      &-(v1(1,1)-v3(1,1))/(y1(1)-y2(1))*(y1(1)-b);
      else if((y1(1)<=b).AND.(x1(1)>a))then
      NZ1=v1(1,1)-(v1(1,1)-v2(1,1))/(x1(1)-x2(1))*(x1(1)-a)&
      &-(v4(1,1)-v1(1,1))/(y3(1)-y1(1))*(y1(1)-b);
      else if(y1(1)>b.AND.x1(1)<=a)then
      NZ1=v1(1,1)-(v5(1,1)-v1(1,1))/(x3(1)-x1(1))*(x1(1)-a)&
      &-(v1(1,1)-v3(1,1))/(y1(1)-y2(1))*(y1(1)-b);
   else
      NZ1=v1(1,1)-(v5(1,1)-v1(1,1))/(x3(1)-x1(1))*(x1(1)-a)&
      &-(v4(1,1)-v1(1,1))/(y3(1)-y1(1))*(y1(1)-b);
    end if
    NZ=NZ1(1,1);
   return;
end subroutine


