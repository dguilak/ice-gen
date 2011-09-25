	subroutine dist(nm,k1,nmm,k2,d)
c min. image distance between atoms
	implicit real*8(a-h,o-z)
	parameter (nmoll=4000,natmm=6000,ndimm=18000)
	common /xlattt/xlat(nmoll,9,2)
	common/parambo/side(3)
	dimension ve(3)
	i1=(k1-1)*3
	i2=(k2-1)*3
	do 10 k=1,3
	   ve(k)=xlat(nm,i1+k,1)-xlat(nmm,i2+k,1)
	   if(ve(k).lt.-side(k)*0.5d0)ve(k)=ve(k)+side(k)
	   if(ve(k).gt.side(k)*0.5d0)ve(k)=ve(k)-side(k)
10	continue
	d=sqrt(ve(1)**2+ve(2)**2+ve(3)**2)
	return
	end

