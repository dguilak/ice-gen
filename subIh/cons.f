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
