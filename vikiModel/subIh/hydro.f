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
