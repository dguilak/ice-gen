	subroutine build
c building a piece bound by hexagonal surfaces
	implicit real*8(a-h,o-z)
c nmoll max num of water molecules
	parameter (nmoll=4000,natmm=6000,ndimm=18000)
	common /unitc/ aerg,aev,acm1,atem,akm,amev,amu,ams,atim,aps
	common/xmmmmm/ xmas(3),xmass(9),xiner(3),xmol
	common/param1/roh,ang,qh,ao2,co2,rom,qmm,qmh,qhh,ff,ff12,ddpl
	common/param2/chap(3,3),a12,c6,fam(9),rhh,rsu,pi,cut2,cutt,dipli
	common/param6/ksha(3,2)
	common/param3/xredi(3),xdis2(3),xdee(3)
	common/param5/roo,rfcc,cutne,ctn2
	common/parambo/side(3)
	common/numbs/nmol,natm,ndim
	dimension dell(3),vv(3),unc(3,4,3)
	dimension ixyz(3)
c coordinates of water
	common /xlattt/xlat(nmoll,9,2)
c	print*,'shift of lattice origin from 0,0,0 in units of'
c	print*,'nearest neigh oo distance?'
c	read*,dell(1),dell(2),dell(3)
	dell(1)=0.01
	dell(2)=0.01
	dell(3)=0.01
c see fig. 28 Kittel p.31
c and fig 31 p. 32
c side of hexagon seen from top
	a=roo*4.0d0/(3.0d0*sqrt(2.0d0))
c lattice vectors
	vv(1)=roo*4.0d0/sqrt(6.0d0)
	vv(2)=a*3.0d0
	vv(3)=roo*(1.0d0+1.0d0/3.0d0)
	print*,'8 atoms in a cell'
	print*,'unit vectors in ams:'
	print*,(vv(k)*ams,k=1,3)
	write(3,*)'a=',a*ams,' vv=',(vv(k)*ams,k=1,3)
700	print*,'nx,ny,nz?'
	print*,'nz=#hex layers must be multiple of 2'
c	print*,'nx probably should be > 2'
c number of repeating units in each direction
	read*,nx,ny,nz
c	if(amin0(nx,ny,nz).lt.2)then
c	   print*,'too small nx or ny or nz'
c	   stop
c	end if
	side(1)=nx*vv(1)
	side(2)=ny*vv(2)
	side(3)=nz*vv(3)
	print*,'sides of box A',side(1)*ams,side(2)*ams,side(3)*ams
	print*,'satisfied? enter 0'
	read*,isat
	if(isat.ne.0)go to 700
	write(3,*)'nx,ny,nz',nx,ny,nz
	write(3,*)'sides of box A',side(1)*ams,side(2)*ams,side(3)*ams
c location of atoms in surface unit cell
	do 100 i=1,3
	  do 100 j=1,4
	   do 100 k=1,3
100		unc(i,j,k)=0.0d0
c unit cell of bottom layer
	unc(1,1,1)=a*sqrt(3.0d0)/2.0d0
	unc(1,2,2)=a*0.5d0
	unc(1,2,3)=roo/3.0d0
	unc(1,3,2)=a*1.5d0
	unc(1,4,1)=a*sqrt(3.0d0)/2.0d0
	unc(1,4,2)=a*2.0d0
	unc(1,4,3)=roo/3.0d0
c unit cell of 2-nd layer
	unc(2,1,1)=a*sqrt(3.0d0)/2.0d0
	unc(2,1,3)=roo/3.0d0
	unc(2,2,2)=a*0.5d0
	unc(2,3,2)=a*1.5d0
	unc(2,3,3)=roo/3.0d0
	unc(2,4,1)=a*sqrt(3.0d0)/2.0d0
	unc(2,4,2)=a*2.0d0
c generate oxygen lattice
	nmol=0
	do 1 iz=0,nz-1
	   ila=mod(iz,2)+1
	   ixyz(3)=iz
	   do 2 ix=0,nx-1
	        ixyz(1)=ix
		do 3 iy=0,ny-1
	            ixyz(2)=iy
		    do 4 iat=1,4
		   	nmol=nmol+1
		   	if(nmol.gt.nmoll)then
		      	  print*,'dimension overrun,increase nmoll',nmol
		      	  stop
		        end if
			do 5 k=1,3
5		           xlat(nmol,k,1)=unc(ila,iat,k)+vv(k)*ixyz(k)
c		   	write(97,*)'iz,ix,iy,iat',iz,ix,iy,iat
c		   	write(97,*)'nmol',nmol
c		   	write(97,*)(xlat(nmol,k,1)*ams,k=1,3)
4		    continue
3	        continue
2	    continue
1	continue
	print*,'nmol=',nmol
	write(3,*)'nmol=',nmol
	natm=nmol*3
	ndim=natm*3
c rescale to the right lattice constant
c shift by vector dell
	do 6 nm=1,nmol
	   do 7 k=1,3
	  xlat(nm,k,1)=xlat(nm,k,1)+dell(k)*roo
	  if(xlat(nm,k,1).gt.side(k))xlat(nm,k,1)=xlat(nm,k,1)-side(k)
	  if(xlat(nm,k,1).lt.0.0d0)xlat(nm,k,1)=xlat(nm,k,1)+side(k)
	   xlat(nm,k,2)=xlat(nm,k,1)
7	   continue
c	   if(nm.lt.6) write(99,*)'nm=',nm,(xlat(nm,k)*ams,k=1,3)
6	continue
	return
	end
