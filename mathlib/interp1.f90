subroutine interp1(x,v,n1,n2,C,NZ)
   implicit none;
 integer*4 i,n1,n2;
   double precision,dimension(n1):: x,v;
 double precision,dimension(n2):: C,NZ;
 double precision,dimension(1)::NZ1,P1,P2,P3,V1,V2,V3;
 do i=1,n2
 P1=x(minloc(abs(x-C(i)))-1);P2=x(minloc(abs(x-C(i))));
 P3=x(minloc(abs(x-C(i)))+1);
 V1=v(minloc(abs(x-C(i)))-1);V2=v(minloc(abs(x-C(i))));
 V3=v(minloc(abs(x-C(i)))+1);
if (C(i)==x(i)) then 
NZ1(1)=v(i);
else
 if (C(i)>P2(1)) then
  NZ1(1)=V2(1)+(V3(1)-V2(1))/(P3(1)-P2(1))*(C(i)-P2(1));
  else 
  NZ1(1)=V2(1)+(V2(1)-V1(1))/(P2(1)-P1(1))*(C(i)-P2(1));
  end if;
end if;
  NZ(i)=NZ1(1);
  end do ;
  
   return;
end subroutine
