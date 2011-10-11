c looking for structure with zero dipole
c here evaluating E-difference with respect to 
c periodic ice crystal, with H in the middle
c a version employing hexagonal cell of ice IH
c version for spc potl, takefrom jpc 79,926 (83?)
c tetrahedral hoh, eliminates O-contribution to E-differences
c remake manc2.f
c adding option for replacing system by a set of dipoles
c remake surh.f, many.f, man1.f
c building models for ice IH bounded by hexagonal surfaces
c CUBIC ice version
c looking for possible structures / energies for given size
c of unit cell
c now writing structures on disk in integer format, 
c to avoid memory problems,
c trying to get all DIFFERENT structures
c Structures with same energy, dipole size
c same num of pair interactions of 4 kinds
c assumed identical up to some symmetry operation
c CAREFUL CUBIC UNIT CELL ASSUMED
C MAY AFFECT RESULTS
c CHECK DIFFERENT UNIT CELL ARRANGEMENT
c periodic boundaries
c using Kitell p.32 fig. 31
c see also G. Kroes Sur Sci 275, 365 (1992)
c TIPS2 potential
c recording structures on fort.9 formatted file
c fort.3 long log, fort.33 short log
	program findstru
	implicit real*8(a-h,o-z)
	write(3,*)'program to calculate different periodic structures'
	write(3,*)'possible for ice Ih, given unit cell'
	write(3,*)'E-dif with respect to periodic'
c set constants relevant to water
	call cons
c build piece of O-structure
	call build
c analysis of structure
	call anal
c generate possible structures by
c distributing  hydrogens in between OO
	call hydro
cc write graphic file for piece of H2O structure
c	call trans(2)
cc write restart file for trajectory
c	call writ
	stop
	end
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
	subroutine cons
	implicit real*8(a-h,o-z)
c prepare data in a.u. for entire program
	common /unitc/ aerg,aev,acm1,atem,akm,amev,amu,ams,atim,aps
	common/xmmmmm/ xmas(3),xmass(9),xiner(3),xmol
	common/param1/roh,ang,qh,ao2,co2,rom,qmm,qmh,qhh,ff,ff12,ddpl
	common/param2/chap(3,3),a12,c6,fam(9),rhh,rsu,pi,cut2,cutt,dipli
	common/param6/ksha(3,2)
	common/param3/xredi(3),xdis2(3),xdee(3)
	common/param4/xinv(3,2),corot(9)
	common/param5/roo,rfcc,cutne,ctn2
	common/param7/cchap(3,3),rroh,rrsu,rrhh
c conversion factors a.u. to others
c erg - au
	aerg=2.2940e10
c ev - au
	aev=3.6752e-2
c cm-1 au
	acm1=4.5563e-6
c Kelvin - au
	atem=3.163e-6
c kcal/mol - au
	akm=1.5931e-3
c milli ev - au
	amev=3.6752e-5
c amu - a.u. mass
	amu=1822.84
c angstrom - au
	ams=0.52917
c au - second
	atim=2.4189e-17
c au - picosecond
	aps=2.4189e-5
        pi=4.0d0*atan(1.0d0)
c atomic masses
	xmas(1)=15.995d0*amu
	xmas(2)=1.008d0*amu
	xmas(3)=1.008d0*amu
	xmol=xmas(1)+xmas(2)+xmas(3)
c parameters of TIPS2 in a.u.
c TIPS2 h charge
c changing parameters to SPC
	qh=0.535d0
	qh=0.41
c "dipole moment" at center of bond
	rom=0.15d0/ams
	rom=0.0d0
	ao2=695000.0d0*akm/ams**12
	ao2=629400.0d0*akm/ams**12
	co2=600.0d0*akm/ams**6
	co2=625.5d0*akm/ams**6
	a12=ao2*12.0d0
	c6=co2*6.0d0
	qm=-2.0*qh
	qmm=qm**2
	qmh=qm*qh
	qhh=qh**2
	chap(1,1)=qhh
	chap(2,2)=qhh
	chap(1,2)=qhh
	chap(2,1)=qhh
	chap(3,3)=qmm
	chap(3,1)=qmh
	chap(1,3)=qmh
	chap(3,2)=qmh
	chap(2,3)=qmh
	roo=2.75/ams
	print*,'roo?'
	read*,roo
	write(3,*)'roo=',roo
	roo=roo/ams
c O-H distance
	roh=0.9572d0/ams
	roh=1.0000d0/ams
	rroh=roo*0.5d0
	rati=(roh/rroh)**2
	do k=1,3
	  do kk=1,3
		cchap(k,kk)=chap(k,kk)*rati
	  end do
	end do
c HOH angle
	ang=104.52 *pi/180.0d0
	ang=109.47 *pi/180.0d0
	rhh=roh*sin(ang/2.0d0)*2.0d0
	rrhh=rroh*sin(ang/2.0d0)*2.0d0
	print*,'rhh=',rhh*ams
	rsu=roh*cos(ang/2.0d0)*2.0d0
	rrsu=rroh*cos(ang/2.0d0)*2.0d0
	rsu2=rsu*0.5d0
c maybe define rsu2 with respect to 109 angle??
	ff=rom/rsu
	ff12=1.0d0-2.0d0*ff
c corot - coordinates of atoms in body fixed rotating system
	do 320 i=1,9
320	corot(i)=0.0d0
	aa=2.0d0*xmas(2)*rsu2/xmol
	corot(3)=aa
	corot(5)=-rhh/2.0d0
	corot(6)=aa-rsu2
	corot(8)=rhh/2.0d0
	corot(9)=aa-rsu2
c XINER MOMENTS OF INERTIA X Y Z = 1 2 3
	xiner(3)=xmas(2)*rhh**2/2.0d0
	xiner(2)=2.0d0*xmas(1)*xmas(2)*rsu2**2/xmol
	xiner(1)=xiner(2)+xiner(3)
c QUANTITIES FOR SHAKE
	xredi(1)=2.0d0/xmas(1)+2.0d0/xmas(2)
	xredi(2)=2.0d0/xmas(2)+2.0d0/xmas(3)
	xredi(3)=2.0d0/xmas(1)+2.0d0/xmas(3)
	xdis2(1)=roh**2
	xdis2(2)=rhh**2
	xdis2(3)=roh**2
	do 26 ic=1,3
	ksha(ic,1)=-1+(ic-1)*3
	ic1=mod(ic,3)
	ksha(ic,2)=-1+ic1*3
	xinv(ic,1)=1.0d0/xmas(ic)
	xinv(ic,2)=1.0d0/xmas(ic1+1)
26	continue
	do 31 i=1,3
31	xdee(i)=1.0/xdis2(i)
c OO distance in crystal cubic, eisenberg+kauzmann
c	dipli=qh*(0.5*roo-roh)*qh
c	ddpl=(0.5*roo-roh)
c	dipli=qh*(0.5*roo-roh)*qh*2*rsu2
	dipli=(qh*2*rsu2)**2
	dne=sqrt(3.0d0)*0.25d0
	rfcc=roo/dne
c cutne is cutoff for nearest neighbor distance, ctn2 - squared
	cutne=(roo+0.2)
	ctn2=cutne**2
	print*,'cutoff for potential [A]?'
	read*,cutt
	write(3,99)cutt
99	format('cutoff dist.[A]:',f15.3)
	cutt=cutt/ams
	cut2=cutt**2
	return
	end
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
	subroutine anal
c analysis of O-structure
	implicit real*8(a-h,o-z)
	parameter (nmoll=4000,natmm=6000,ndimm=18000)
	parameter(nemax=4)
c nemax is maximal number of nearest neighbors allowed
	common /unitc/ aerg,aev,acm1,atem,akm,amev,amu,ams,atim,aps
	common/xmmmmm/ xmas(3),xmass(9),xiner(3),xmol
	common/param1/roh,ang,qh,ao2,co2,rom,qmm,qmh,qhh,ff,ff12,ddpl
	common/param2/chap(3,3),a12,c6,fam(9),rhh,rsu,pi,cut2,cutt,dipli
	common/param6/ksha(3,2)
	common/param3/xredi(3),xdis2(3),xdee(3)
	common/param4/xinv(3,2),corot(9)
	common/param5/roo,rfcc,cutne,ctn2
	common /xlattt/xlat(nmoll,9,2)
	common/xnnneigs/ ve(3),xnei(nmoll,nemax,3)
	common/nnneigs/ nei(nmoll,nemax),nne(nmoll)
	dimension hisn(0:nemax),box(3,100)
	common/numbs/nmol,natm,ndim
	common/parambo/side(3)
c nne counts number of nearest neighbors
c nei contains matrix of nearest neighbors
c hisn - histogram of neighbor numbers
	do 3 nm=1,nmol
3	   nne(nm)=0
	do 4 ne=0,nemax
4	   hisn(ne)=0
c locate nearest neighbors
	do 1 nm=1,nmol-1
	   do 2 nmm=nm+1,nmol
	    	do 10 k=1,3
		   ve(k)=xlat(nm,k,1)-xlat(nmm,k,1)
c		   if(k.eq.3) go to 10
		   if(ve(k).lt.-side(k)*0.5d0)ve(k)=ve(k)+side(k)
		   if(ve(k).gt.side(k)*0.5d0)ve(k)=ve(k)-side(k)
10	  	continue
	   	d2=ve(1)**2+ve(2)**2+ve(3)**2
	        if(d2.lt.ctn2)then
c		   if(nm.eq.1)then
c			d22=sqrt(d2)*ams
c		 	ctn22=sqrt(ctn2)*ams
c			print*,d22,ctn22,nm,nmm
c			write(98,*)nm,nmm,d22,ctn22
c			write(98,*)(xlat(nm,k,1)*ams,k=1,3)
c			write(98,*)(xlat(nmm,k,1)*ams,k=1,3)
c		   end if
		   nne(nm)=nne(nm)+1
		   if(nne(nm).gt.nemax)then
		      print*,'near neigh dim. overrun'
		      print*,nm,nne(nm)
		      stop
		   end if
		   nei(nm,nne(nm))=nmm
		   nne(nmm)=nne(nmm)+1
		   if(nne(nmm).gt.nemax)then
		      print*,'near neigh dim. overrun'
		      print*,nmm,nne(nmm)
		      stop
		   end if
		   nei(nmm,nne(nmm))=nm
		   do 11 k=1,3
		 	xnei(nm,nne(nm),k)=-ve(k)
		 	xnei(nmm,nne(nmm),k)=ve(k)
11		   continue
	        end if
2	   continue
1	continue
	do 5 nm=1,nmol
	   n=nne(nm)
	   hisn(n)=hisn(n)+1
5	continue
	su=0.0d0
	do 6 ne=0,nemax
	   su=su+hisn(ne)
6	   print*,'coor:',ne,hisn(ne)
c cosine of tetrahedral angle from Kittel p.31
        cost=-1.0d0/3.0d0
c check tetrahedrality of bonds
	do 12 nm=1,nmol
	   do 13 k=1,nne(nm)-1
 	    do 14 kk=k+1,nne(nm)
		co=xnei(nm,k,1)*xnei(nm,kk,1)+
     x		xnei(nm,k,2)*xnei(nm,kk,2)+
     x		xnei(nm,k,3)*xnei(nm,kk,3)
		co=co/roo**2
		if(abs(co-cost).gt.1.d-5)then
		   print*,'no good cos',co,cost
		   stop
		end if
14	    continue
13	   continue
12	continue
	return
	end
	subroutine hydro
c performs searches for different structures obeying ice rules
c spread hydrogens at random in between pairs of OO
c according to ice rules
c need to supply nearest neighbor list nei(nmoll,nemax)
c and coordinates of nearest neighbors xnei with respect to O-atom
c H spread initially at random between OO
c then random walk of H from higher coor O to lower coor O
c probablity 1 of step if coor diffr bigger than 1
c probablity 0.5 of step if coor diffr eq 1
c finally feeds H-coors into xlat
c water placed in plane defined by two OO vectors
c HOH angle set to "true" water value ang, not tetrahedral 109
C MAYBE BETTER TO USE TRUE ICE VALUE! BIGGER ANgle
C AND LONGER BOND - LIKE IN ICE NOT TIPS2 GAS
c bisector of H2O alligned along bisector of two OO vectors
	parameter (nmoll=4000,natmm=6000,ndimm=18000)
	parameter(nemax=4,nsmax=10000)
c nsmax is max num of distinct structures expected to be found
	implicit real*8(a-h,o-z)
	common /unitc/ aerg,aev,acm1,atem,akm,amev,amu,ams,atim,aps
	common /xlattt/xlat(nmoll,9,2)
	common/param1/roh,ang,qh,ao2,co2,rom,qmm,qmh,qhh,ff,ff12,ddpl
	common/param2/chap(3,3),a12,c6,fam(9),rhh,rsu,pi,cut2,cutt,dipli
	common/param5/roo,rfcc,cutne,ctn2
	common/param6/ksha(3,2)
	common/param7/cchap(3,3),rroh,rrsu,rrhh
	common/pottts/eppo(nmoll),eppp,ehh,eldd
	common/catttts/icat(4)
	common/xnnneigs/ ve(3),xnei(nmoll,nemax,3)
	common/nnneigs/ nei(nmoll,nemax),nne(nmoll)
	common/numbs/nmol,natm,ndim
	common/parambo/side(3)
	common/xladip/xld(nmoll*2,6),dipm(nmoll,6)
c xld 1-3 bond center 4-6 bond dipole
c dipm 1-3 O-location 4-6 mol dipole
c for dipolar potl. approximation
	dimension rh(10),vv(2,3),bis(3),vhh(3),box(3,100),bmax(3)
	integer*2 ne2(nmoll,nemax)
	dimension nih(0:nemax),nnh(nmoll)
c nnh(nm) - num of OH-bonds
c ne2(nm,k)=1 H belongs to nm-th O
c ne2(nm,k)=-1 H belongs to neighbor of nm-th O 
c (neighbor index nmm=nei(nm,k))
c	dimension neli(nsmax,nmoll,nemax)
	integer*2 neli(nmoll,nemax)
c neli is listing of H-atom structures, 
c in ne2 format
c to compare to new structures & verify if new struct. obtained
	dimension dip(4),dipo(nsmax,4),epo(nsmax,4),epdo(nsmax)
	dimension inews(nsmax),icatt(nsmax,4)
c dip is dipole vector, last component 4 is abs size
c dipole calculated as vectorial sum in units of molecular
c dipole, set to unity along the bisector
	print*,'dseed?'
	read*,dseed
	write(3,*)'dseed=',dseed
	print*,'# improvement rounds, #steps in round? (eg~1000,10000)'
	read*,nimp,nmcma
	print*,'# max tries?'
	read*,maxt
	print*,'after how many futile attempts to stop?'
	read*,isto
	print*,'e-,dip. identity criterions?'
	read*,cri1,cri2
	write(3,*)'# improvement rounds:',nimp
	write(3,*)'# mc steps in a round',nmcma
	write(3,*)'stopping criterion: futile tries:',isto
	write(3,*)'max tries:',maxt
	write(3,*)'criterions:',cri1,cri2
	close (3)
	open(file='fort.3',form='formatted',status='old',
     x	unit=3,access='append')
c nst # distinct structs found
c ntry # of tries
c nfail # consecutive failures
	nst=0
	ntry=0
	nfail=0
	nfa1=0
	nfa2=0
5555	ntry=ntry+1
c	print*,'dseed?'
c	read*,dseed
c	print*,'# try:',ntry
c 	write(3,*)'# try:',ntry
	do 1 nm=1,nmol
	  nnh(nm)=0
	  do 11 ne=1,4
11	     ne2(nm,ne)=0
1	continue
c distribute H-atoms at random between OO pairs
	do 2 nm=1,nmol
	  nr=4
  	  call ggubs(dseed,nr,rh)
	  do 3 k=1,4
		if(ne2(nm,k).ne.0)go to 3
		if(rh(k).lt.0.5d0)then
		   ne2(nm,k)=1
c H belongs to oxygen nm
		   nnh(nm)=nnh(nm)+1
		   nmm=nei(nm,k)
		   do 4 kk=1,4
			if(nei(nmm,kk).ne.nm)go to 4
			ne2(nmm,kk)=-1
			go to 3
4		   continue
		else
		   ne2(nm,k)=-1
c H belongs to neighbor nmm
		   nmm=nei(nm,k)
		   do 5 kk=1,4
			if(nei(nmm,kk).ne.nm)go to 5
			ne2(nmm,kk)=1
			nnh(nmm)=nnh(nmm)+1
			go to 3
5		   continue
		end if
3	  continue
2	continue
c random walk: walk H-atoms from higher coordinated O to lower coord. O
c until all O have coordination 2
	do 6 kmc=1,nimp
	  do 8 ne=0,4
8	      nih(ne)=0
	  do 9 nm=1,nmol
	      n=nnh(nm)
	      nih(n)=nih(n)+1
9	  continue
c	  print*,'O-coordination table:'
c	  do 12 ne=0,4
c12	      print*,ne,nih(ne),nmol
	  if(nih(2).eq.nmol)then
c finished - all O 2-coordinated
	     go to 16
	  end if
	  do 7 nmc=1,nmcma
c pick at random a molecule and one of its 4 neighbors
	     nr=3
  	     call ggubs(dseed,nr,rh)
	     nm=(dble(nmol)*rh(1)+1.0d0)
	     km=(4.0d0*rh(2)+1.0d0)
	     nmm=nei(nm,km)
	     i0=0
c what would be new coord after exchange of H?
	     k1=nnh(nm)-ne2(nm,km)
	     k2=nnh(nmm)+ne2(nm,km)
	     kd=iabs(k2-k1)
	     nd=iabs(nnh(nm)-nnh(nmm))
c if diffr. of coordinations would decrease , carry out the step
	     if(kd.lt.nd)i0=1
c if diffr. of coordinations would remain same, carry out the step
c with prob. 0.5
	     if(kd.eq.nd.and.rh(3).gt.0.5d0)i0=1
	     if(i0.eq.0)go to 7
	     nnh(nm)=k1
	     nnh(nmm)=k2
 	     ne2(nm,km)=-ne2(nm,km)
	     do 13 kk=1,4
		if(nei(nmm,kk).eq.nm)then
		   ne2(nmm,kk)=-ne2(nmm,kk)
		   go to 7
	        end if
13	     continue
7	  continue
6	continue
c failed attempt
	nfail=nfail+1
	nfa1=nfa1+1
c	write(3,*)'failed to converge,nfail=',nfail
c	print*,'failed to converge'
	if(nfail.lt.isto.and.ntry.lt.maxt)go to 5555
	write(3,*)'program finished,nfail=',nfail
	print*,'program finished,nfail=',nfail
	go to 2666
16	continue
c calculate configuration, energy,dipole
c check if new structure
c using energy, dipole equality criterion
	do k=1,4
	   dip(k)=0.0d0
	end do
	do 14 nm=1,nmol
	   kk=0
	   do 15 k=1,4
		if(ne2(nm,k).ne.1)go to 15
		kk=kk+1
		do 17 j=1,3
		   vv(kk,j)=xnei(nm,k,j)
17		continue
15	   continue
	   if(kk.ne.2)then
	        write(3,*)'something sick with num. of OH bonds',kk
	        print*,'something sick with num. of OH bonds',kk
		stop
	   end if
	   sb=0.0d0
	   sh=0.0d0
	   do 18 j=1,3
		bis(j)=vv(1,j)+vv(2,j)
		sb=sb+bis(j)**2
		vhh(j)=vv(1,j)-vv(2,j)
		sh=sh+vhh(j)**2
18	   continue
	   sb1=1.0/sqrt(sb)
	   sh1=1.0/sqrt(sh)
	   sb=0.5d0*rsu*sb1
	   sh=0.5d0*rhh*sh1
	   ssb=0.5d0*rrsu*sb1
	   ssh=0.5d0*rrhh*sh1
	   do 19 j=1,3
		xlat(nm,3+j,1)=xlat(nm,j,1)+bis(j)*sb+vhh(j)*sh
		xlat(nm,6+j,1)=xlat(nm,j,1)+bis(j)*sb-vhh(j)*sh
		xlat(nm,3+j,2)=xlat(nm,j,1)+bis(j)*ssb+vhh(j)*ssh
		xlat(nm,6+j,2)=xlat(nm,j,1)+bis(j)*ssb-vhh(j)*ssh
		dip(j)=dip(j)+bis(j)*sb1
	 	dipm(nm,j)=xlat(nm,j,1)
		dipm(nm,j+3)=bis(j)*sb1
19	   continue
14	continue
	do j=1,3
	   dip(j)=dip(j)/nmol
	end do
	dip(4)=sqrt(dip(1)**2+dip(2)**2+dip(3)**2)
c check potential energy
	call epot
	ep1=eppp/nmol/akm
	eh1=ehh/nmol/akm
	eld1=eldd/nmol/akm
c check if this structure is new
c energy criterion+dipole equality+different types of
c pair interaction number equality
	do 5001 nnn=1,nst
	   dd1=abs(ep1-epo(nnn,1))
	   dd2=abs(dip(4)-dipo(nnn,4))
	   idi=iabs(icat(1)-icatt(nnn,1))+iabs(icat(2)-icatt(nnn,2))+
     x	   iabs(icat(3)-icatt(nnn,3))+iabs(icat(4)-icatt(nnn,4))
	   if(dd1.lt.cri1.and.dd2.lt.cri2.and.idi.eq.0)then
	      inews(nnn)=inews(nnn)+1
	      nfail=nfail+1
	      nfa2=nfa2+1
	      if(nfail.lt.isto.and.ntry.lt.maxt)go to 5555
	      write(3,*)'program finished,nfail=',nfail
	      print*,'program finished,nfail=',nfail
	      go to 2666
	   end if
5001	continue
c new structure found
c write xlat epot dipole
c also, calculate dipole per molecule
	nfail=0
	nst=nst+1
	do k=1,4
	   dipo(nst,k)=dip(k)
	   icatt(nst,k)=icat(k)
	end do
	inews(nst)=1
	epo(nst,1)=ep1
c translate to system of dipoles
	ndip=0
	do nm=1,nmol
	   kk=0
	   do 155 k=1,4
		if(ne2(nm,k).ne.1)go to 155
		kk=kk+1
		ndip=ndip+1
		do j=1,3
		   xld(ndip,j)=xlat(nm,j,1)+0.5*xnei(nm,k,j)
	if(xld(ndip,j).gt.side(j))xld(ndip,j)=xld(ndip,j)-side(j)
	if(xld(ndip,j).lt.0.0d0)xld(ndip,j)=xld(ndip,j)+side(j)
		   xld(ndip,3+j)=-xnei(nm,k,j)/roo
c the size of dipole is set here to unity
		end do
155	   continue
	end do
	call enedip(epd)
	ed1=epd/nmol/akm
	epdo(nst)=ed1
c	epo(nst,2)=eh1
	epo(nst,2)=eld1
	epo(nst,3)=ed1
	epo(nst,4)=dip(4)
	write(33,3334)ep1,eld1,ed1,dip(4),icat,nst
	close (33)
	open(file='fort.33',form='formatted',status='old',
     x	unit=33,access='append')
	if(dip(4).lt.1.d-6)then
	   call writ(ep1,dip,4)
	   print*,'found zero dipole struct'
	   write(3,*)'found zero dipole struct'
	   go to 2666
	end if
	if(nst.lt.nsmax.and.ntry.lt.maxt)go to 5555
2666	continue
	write(3,*)'finished,# structs found,requested',nst,nsmax
	write(3,*)'tries made, max:',ntry,maxt
	print*,'finished,# structs found,requested',nst,nsmax
	print*,'tries made, max:',ntry,maxt
	print*,'2 kinds of failures:',nfa1,nfa2
	write(3,*)'non-convergence failures:',nfa1
	write(3,*)'repeat stucture failures:',nfa2
c write final listing of structures
	rewind 33
	do nn=1,nst
	 write(33,3334)epo(nn,1),epo(nn,2),epo(nn,3),(dipo(nn,k),k=4,4),
     x	   (icatt(nn,k),k=1,4),inews(nn)
3334	   format(3f9.4,e12.5,4i4,i6)
	end do
c order final structures by energy
c write ordered list
1002	continue
	i0=0
	do 1003 i=1,nst-1
	   if(epo(i,1).le.epo(i+1,1))go to 1003
	   do k=1,4
	      xx=epo(i,k)
	      epo(i,k)=epo(i+1,k)
	      epo(i+1,k)=xx
	   end do
	   i0=1
1003	continue
	if(i0.eq.1)go to 1002
	write(3,*)'lowest energy tot/eld/edip/dipole'
	write(3,*)(epo(1,k),k=1,4)
	print*,'lowest energy tot/ehh/edip/dipole',(epo(1,k),k=1,4)
	do i=1,nst
	   dii=epo(i,1)-epo(1,1)
	   ra=1.0
	   if(dii.ne.0.0d0)ra=(epo(i,3)-epo(1,3))/dii
	   write(34,666)epo(i,1)-epo(1,1),epo(i,4)
	   write(35,666)epo(i,2)-epo(1,2),epo(i,4)
	   write(36,666)epo(i,3)-epo(1,3),epo(i,4),ra
666	format(f10.4,f12.5,f10.5)
	end do
	isu=icat(1)+icat(2)+icat(3)+icat(4)
	write(3,*)'total bonds',isu
	print*,'total bonds',isu
cc check distribution of distances
c	do 26 k=1,2
c	  do 24 n=1,100
c24	   box(k,n)=0.0d0
c26	continue
c	do 20 nm=1,nmol
c	  call dist(nm,2,nm,3,d)
c	  n=(d*ams/0.1d0+1.0d0)
c	  if(n.lt.100)box(2,n)=box(2,n)+1.0d0
c	  do 22 nmm=nm+1,nmol
c	     call dist(nm,1,nmm,1,d)
c	     n=(d*ams/0.1d0+1.0d0)
c	     if(n.lt.100)box(1,n)=box(1,n)+1.0d0
c	     do 21 k=2,3
c	        do 23 kk=2,3
c	           call dist(nm,k,nmm,kk,d)
c	           n=(d*ams/0.1d0+1.0d0)
c	           if(n.lt.100)box(2,n)=box(2,n)+1.0d0
c23		continue
c21	     continue
c22	  continue
c20	continue
c	do 27 k=1,2
c	   bmax(k)=0.0d0
c	   do 28 n=1,100
c28		if(box(k,n).gt.bmax(k))bmax(k)=box(k,n)
c	   do 25 n=1,100
c	        r=(n-0.5d0)*0.1
c	        write(90+k,*)r,box(k,n)/bmax(k)
c25	   continue
c27	continue
	return
	end
      SUBROUTINE GGUBS (DSEED,NR,R)
	implicit real*8 (a-h,o-z)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL*8           R(NR)
      DOUBLE PRECISION   DSEED
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL*8           S2P31,S2P31M,SEED
C                                  S2P31M = (2**31) - 1
C                                  S2P31 = (2**31)
      DATA               S2P31M/2147483647.E0/,S2P31/2147483648.E0/
C                                  FIRST EXECUTABLE STATEMENT
      SEED = DSEED
      DO 5 I=1,NR
         SEED = dMOD(16807.E0*SEED,S2P31M)
    5 R(I) = SEED / S2P31
      DSEED = SEED
      RETURN
      END
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
