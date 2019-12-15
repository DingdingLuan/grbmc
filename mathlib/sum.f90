subroutine sum(v,n,SU)
    implicit none ;
    integer*4 i,n
    double precision,dimension(n):: v
    double precision:: SU
SU=0.0d0;
   do i=1,n
 SU=SU+v(i);
    enddo ;
   return;
end subroutine 
