	subroutine epot
C 	POTENTIAL ENERGY TIPS2 BETWEEN PAIRS OF MOLECULES
c also, difference between perfect crystal, and shifted
C CAREFUL USING xlat not coo COORDINATES!
c including possibility of potl cutoff longer than 1/2 unit cell
c using periodic boundary conditions
c distance checkup on O only no H
	implicit real*8(a-h,o-z)
	parameter (nmoll=4000,natmm=6000,ndimm=18000)
	common/numbs/nmol,natm,ndim
c CAREFUL COMPLICATION OF INTERACTION WITH SELF IMAGE
c factor 0.5 (?!)
	common /xlattt/xlat(nmoll,9,2)
c	common /coovel/coo(ndimm),cob(ndimm),coa(ndimm),vel(ndimm)
	dimension y11(12),y1(12),y2(12),ve(3)
	dimension z11(12),z1(12),z2(12),we(3)
	common /unitc/ aerg,aev,acm1,atem,akm,amev,amu,ams,atim,aps
	common/pottts/eppo(nmoll),eppp,ehh,eldd
	common/param1/roh,ang,qh,ao2,co2,rom,qmm,qmh,qhh,ff,ff12,ddpl
	common/param2/chap(3,3),a12,c6,fam(9),rhh,rsu,pi,cut2,cutt,dipli
	common/parambo/side(3)
	common/param5/roo,rfcc,cutne,ctn2
	common/param7/cchap(3,3),rroh,rrsu,rrhh
	dimension lmax(3)
	common/catttts/icat(4)
	cu=sqrt(cut2)
c icat includes categories of pair interactions
	do k=1,4
	   icat(k)=0
	end do
	do k=1,3
	   xl=0.5d0+cu/side(k)
	   lmax(k)=xl
	end do
cc	print*,'got to epot'
	ep=0.0d0
	ehh=0.0d0
	eldd=0.0d0
	do 77 i=1,nmol
77	   eppo(i)=0.0d0
c	nm1=nmol-1
c	do 201 km=1,nm1

	do 201 km=1,nmol
	 kmm=km
c	 kmm=km+1
	 do 202 km2=kmm,nmol
	    do 1 i=1,9
	       z11(i)=xlat(km,i,2)
1	       y11(i)=xlat(km,i,1)
	    do 402 i=1,9
	       z2(i)=xlat(km2,i,2)
402	       y2(i)=xlat(km2,i,1)
c find y11 - nearest image km location
	    do 4302 k=1,3
		ve(k)=y11(k)-y2(k)
		if(ve(k).gt.side(k)*0.5d0)then
		    ve(k)=ve(k)-side(k)
		    do 4303 j=k,9,3
			z11(j)=z11(j)-side(k)
4303			y11(j)=y11(j)-side(k)
		end if
		if(ve(k).lt.-side(k)*0.5d0)then
		    ve(k)=ve(k)+side(k)
		    do 4304 j=k,9,3
			z11(j)=z11(j)+side(k)
4304			y11(j)=y11(j)+side(k)
		end if
4302	    continue
c	    d2=ve(1)**2+ve(2)**2+ve(3)**2
c	    if(d2.gt.cut2)go to 202
c   go over all unit cells
	    kd=iabs(km2-km)
	    do l1=-lmax(1),lmax(1)
	    do l2=-lmax(2),lmax(2)
	    do l3=-lmax(3),lmax(3)
	      kl=iabs(l1)+iabs(l2)+iabs(l3)
c skip interaction with itself 
	      if((kl+kd).ne.0)then
	       	do i=0,6,3
	           y1(i+1)=y11(i+1)+l1*side(1)
		   y1(i+2)=y11(i+2)+l2*side(2)
		   y1(i+3)=y11(i+3)+l3*side(3)
	           z1(i+1)=z11(i+1)+l1*side(1)
		   z1(i+2)=z11(i+2)+l2*side(2)
		   z1(i+3)=z11(i+3)+l3*side(3)
		end do
	        do k=1,3
		   ve(k)=y1(k)-y2(k)
		end do
c   LENNARD JONES
	     	d2=ve(1)**2+ve(2)**2+ve(3)**2
	     	if(d2.gt.cut2)go to 1202
		ine=0
		if(d2.lt.ctn2)ine=1
	     	d6=d2*d2*d2
	     	addd=(ao2/d6-co2)/d6
c   COULOMBIC
c   find m poxints
	     	do 3 k=1,3
	           z1(9+k)=z1(k)*ff12+(z1(k+3)+z1(k+6))*ff
3	           y1(9+k)=y1(k)*ff12+(y1(k+3)+y1(k+6))*ff
	     	do 4 k=1,3
	           z2(9+k)=z2(k)*ff12+(z2(k+3)+z2(k+6))*ff
4	           y2(9+k)=y2(k)*ff12+(y2(k+3)+y2(k+6))*ff
		eh=0.0d0
		ehd=0.0d0
	     	do 5 k1=3,9,3
	     	do 5 k2=3,9,3
	          di=sqrt((y1(k1+1)-y2(k2+1))**2+(y1(k1+2)-y2(k2+2))**2+
     x	          (y1(k1+3)-y2(k2+3))**2)
	          ad=chap(k1/3,k2/3)/di
	         ddi=sqrt((z1(k1+1)-z2(k2+1))**2+(z1(k1+2)-z2(k2+2))**2+
     x	          (z1(k1+3)-z2(k2+3))**2)
	          aad=cchap(k1/3,k2/3)/ddi
	          addd=addd+ad
c		  if(k1.ne.9.and.k2.ne.9)eh=eh+ad
		  eh=eh+ad
		  ehd=ehd+ad-aad
5	     	continue
	     	eppo(km)=eppo(km)+addd*0.5d0
		if(kd.eq.0)then
	     	   ep=ep+addd*0.5d0
		   ehh=ehh+eh*0.5d0
		   eldd=eldd+ehd*0.5d0
		else
	     	   ep=ep+addd
		   ehh=ehh+eh
		   eldd=eldd+ehd
	     	   eppo(km2)=eppo(km2)+addd*0.5d0
		end if
		if(ine.eq.1)then
c nearest neighbors
		   ae=addd/akm
		   if(ae.lt.-6.0)icat(1)=icat(1)+1
		   if(ae.lt.-5.5.and.ae.ge.-6.0)icat(2)=icat(2)+1
		   if(ae.lt.-4.45.and.ae.ge.-5.5)icat(3)=icat(3)+1
		   if(ae.ge.-4.45)icat(4)=icat(4)+1
		end if
1202	        continue
	      end if
	    end do
	    end do
	    end do
202	 continue
201	continue
c	print*,'finished epot'
	eppp=ep
c	do 111 nm=1,nmol
c111	write(199,*)nm,eppo(nm)/akm
	return
	end

