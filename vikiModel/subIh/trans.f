	subroutine trans(iho)
	implicit real*8(a-h,o-z)
c write graphic file for oxygen lattice
c if iho=1 oxygens only
c =2 include hydrogens
	parameter (nmoll=4000,natmm=6000,ndimm=18000)
	common /unitc/ aerg,aev,acm1,atem,akm,amev,amu,ams,atim,aps
	common/numbs/nmol,natm,ndim
	common /xlattt/xlat(nmoll,9,2)
	dimension ve(3)
	character*5 atoc(3),wat,car
	data atoc,wat,car/'    O','    H','    H','  H2O','    C'/
	write(99,777)
777	format('* title')
	write(99,778)
778	format('*')
c CHARMM COORDINATES
	na=nmol
	iaa=1
	if(iho.eq.2)then
	   na=natm
           iaa=3
	end if
305	format(2i5,2a5,3f10.5,2a5,f10.5)
	write(99,305) na
	kaa=0
	do 303 nm=1,nmol
	   do 304 ia=1,iaa
	      kaa=kaa+1
	      io=(ia-1)*3
	   if(iho.eq.1)then
	      do 3306 j=1,3
3306	         ve(j)=xlat(nm,io+j,1)*ams*0.5
	      write(99,305)kaa,nm,wat,car,ve,wat,car,dble(ia)
	   else
	      do 306 j=1,3
306	         ve(j)=xlat(nm,io+j,1)*ams
	      write(99,305)kaa,nm,wat,atoc(ia),ve,wat,atoc(ia),dble(ia)
	   end if
304	continue
303	continue
	return
	end
