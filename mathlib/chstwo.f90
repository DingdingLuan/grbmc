SUBROUTINE chstwo(bins1,bins2,nbins,knstrn,df,chsq,prob)
INTEGER*4 knstrn,nbins
REAL*8 chsq,df,prob,bins1(nbins),bins2(nbins)
!CU    USES gammq
INTEGER*4 j
REAL*8 gammq
df=nbins-knstrn
 chsq=0.
do j=1,nbins
  if(bins1(j).eq.0..and.bins2(j).eq.0.)then
    df=df-1.
  else
    chsq=chsq+(bins1(j)-bins2(j))**2/(bins1(j)+bins2(j))
  endif
enddo
prob=gammq(0.5*df,0.5*chsq)
return
END
