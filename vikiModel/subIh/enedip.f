	subroutine enedip(ep)
c calculate interaction energy for collection of dipoles
c each molecule represented as dipole
	implicit real*8 (a-h,o-z)
	parameter (nmoll=4000,natmm=6000,ndimm=18000)
	common/parambo/side(3)
	common/numbs/nmol,natm,ndim
	common/xladip/xld(nmoll*2,6),dipm(nmoll,6)
	common/param2/chap(3,3),a12,c6,fam(9),rhh,rsu,pi,cut2,cutt,dipli
	dimension y1(6),y2(6),y11(6),lmax(3),ve(3)
	ndip=nmol*2
	do k=1,3
	   xl=0.5d0+cutt/side(k)
	   lmax(k)=xl
	end do
cc	print*,'got to epot'
	ep=0.0d0
	do 201 km=1,nmol
	 kmm=km
c	 kmm=km+1
	 nd=0
	 do 202 km2=1,nmol
c	  do 2202 nbo=1,2
	    nd=nd+1
	    do 1 i=1,6
1	       y2(i)=dipm(km,i)
c1	       y11(i)=dipm(km,i)
	    do 402 i=1,6
402	       y11(i)=dipm(km2,i)
c402	       y2(i)=xld(nd,i)
c find y11 - nearest image km location
	    do 4302 k=1,3
		ve(k)=y11(k)-y2(k)
		if(ve(k).gt.side(k)*0.5d0)y11(k)=y11(k)-side(k)
		if(ve(k).lt.-side(k)*0.5d0)y11(k)=y11(k)+side(k)
4302	    continue
c   go over all unit cells
	    kd=iabs(km2-km)
	    do l1=-lmax(1),lmax(1)
	    do l2=-lmax(2),lmax(2)
	    do l3=-lmax(3),lmax(3)
	      kl=iabs(l1)+iabs(l2)+iabs(l3)
c skip interaction with itself 
	      if((kl+kd).ne.0)then
	        y1(1)=y11(1)+l1*side(1)
		y1(2)=y11(2)+l2*side(2)
		y1(3)=y11(3)+l3*side(3)
	        do k=1,3
		   ve(k)=y1(k)-y2(k)
		   y1(3+k)=y11(k+3)
		end do
	     	d2=ve(1)**2+ve(2)**2+ve(3)**2
	     	if(d2.gt.cut2)go to 1202
		di=sqrt(d2)
		do k=1,3
		   ve(k)=ve(k)/di
		end do
		pr1=ve(1)*y1(4)+ve(2)*y1(5)+ve(3)*y1(6)
		pr2=ve(1)*y2(4)+ve(2)*y2(5)+ve(3)*y2(6)
		pr3=y1(4)*y2(4)+y1(5)*y2(5)+y1(6)*y2(6)
		addd=(pr3-3.0*pr1*pr2)/(di*d2)
	     	ep=ep+addd*0.5
c		if(kd.eq.0)then
c	     	   ep=ep+addd*0.5d0
c		else
c	     	   ep=ep+addd
c		end if
1202	        continue
	      end if
c	      if(km.eq.1.and.km2.eq.2.and.kl.eq.0)xeli=addd
	    end do
	    end do
	    end do
2202	  continue
202	 continue
201	continue
c remove inter between dipoles of same molecule
c	ep=ep-xeli*nmol
c rescale to represent correct dipole-dipole interaction
	ep=ep*dipli
c	ep=ep/nmol
	return
	end
