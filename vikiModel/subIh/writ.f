	subroutine writ(ep1,dip,n4)
c write output config - everything in au except dipole
c see subroutin hydro
	implicit real*8(a-h,o-z)
	parameter (nmoll=4000,natmm=6000,ndimm=18000)
	common /xlattt/xlat(nmoll,9,2)
	common/numbs/nmol,natm,ndim
	common/parambo/side(3)
	dimension dip(n4)
	write(9,*)nmol,ep1,(dip(k),k=1,4),(side(k),k=1,3)
	do 1 nm=1,nmol
	   write(9,*)(xlat(nm,k,1),k=1,9)
1	continue
	close (9)
	open(file='fort.9',form='formatted',status='old',
     x	unit=9,access='append')
	return
	end

