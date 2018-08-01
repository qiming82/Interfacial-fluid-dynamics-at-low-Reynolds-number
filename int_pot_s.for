      subroutine int_pot_s
     +
     +   (z1,r1
     +   ,s1,s2
     +   ,z0,r0,s0
     +   ,n,d
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,y1
     +   )
      Implicit Double Precision (a-h,o-z)

	dimension x1(10),w1(10),x2(10),w2(10),amg(10),ams(10)
	dimension anzg(10),anrg(10),dzg(10),drg(10),dsg(10)
      dimension Bg(10),Bgi(10),Bs(10)   
      dimension sg(10),ss(10),rg(10),rs(10),zg(10),zs(10)
      dimension dss(10),drs(10),dzs(10)
	dimension xh2(10),rh2(10),akg(10),akgi(10),aks(10)
	dimension pf1(10),pf2(10),pe1(10),pe2(10)
      dimension Gi(10),G(10),Fg(10),Fgi(10),Eg(10),Egi(10)
	dimension Gz(10),Gr(10),Gzi(10),Gri(10)
	dimension GS(10),GSi(10),GSp(10),Gp(10),Gzp(10),Grp(10)
	dimension GSz(10),GSr(10),GSzi(10),GSri(10),GSzp(10),GSrp(10)

	eps = 0.0000000001D0
      
	pi  = 3.1415 92653 58979 32384 D0
      RL=1.D0
	d2=2.0D0*d
      
	call Gauss_Legendre (8,x1,w1)
	call LogGL (8,x2,w2)

      if (n.eq.1) then
	  do i=1,8
	     amg(i)=(1.0D0-x1(i))/2.0D0
	   end do
	else if (n.eq.2) then
	   do i=1,8
	      amg(i)=(1.0D0+x1(i))/2.0D0
	    end do
	end if

	a = (s2-s1)/2.0D0
	b = (s2+s1)/2.0D0

	do j=1,8
	   sg(j) = a*x1(j) + b
c	   anzg(j) = (anx2-anx1)/2.0D0*x1(j) + (anx2+anx1)/2.0D0
c	   anrg(j) = (any2-any1)/2.0D0*x1(j) + (any2+any1)/2.0D0
	   xd = sg(j) - s1
	   rg(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	   drg(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	   zg(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	   dzg(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	   dsg(j) = dsqrt(drg(j)**2+dzg(j)**2)
	   anzg(j) = -drg(j)/dsg(j)
	   anrg(j) = dzg(j)/dsg(j)
	   akg(j) = 4.0D0*r0*rg(j)/((zg(j)-z0)**2+(rg(j)+r0)**2)
c	   akgi(j) = 4.0D0*r0*rg(j)/((zg(j)+z0)**2+(rg(j)+r0)**2)
	end do
 
      Iopt = 1  
c ipot  1    both green's function and its derivative

	do j=1,8
      	ssss = d2 - rg(j) - r0
	    sss2 = ssss**2
	    call B0(rg(j),zg(j),r0,z0,Bg(j))
	    call B0(rg(j),zg(j),r0,-z0,Bgi(j))
	    call ell_int2(akg(j),Fg(j),Eg(j))
c	    call ell_int(akgi(j),Fgi(j),Egi(j))
	    call polyn (akg(j),pf1(j),pf2(j),pe1(j),pe2(j))
          call lgf_ax_fs(Iopt,zg(j),rg(j),z0,r0,G(j),Gz(j),Gr(j))
          call lgf_ax_fs(Iopt,zg(j),rg(j),-z0,r0,Gi(j),Gzi(j),Gri(j))
          call lgf_ax_fs(Iopt,zg(j),rg(j),RL-z0,r0,Gp(j),Gzp(j),Grp(j))
        dx = z0-zg(j)
	    G(j)=G(j)-0.5d0/dsqrt(dx*dx+sss2)
c	    G(j)=G(j)-0.5D0/dabs(zg(j)+RL/2.D0)
c     +	   -0.25d0/dsqrt(dx*dx+sss2)
        dx = -z0-zg(j)
	    Gi(j)=Gi(j)-0.5d0/dsqrt(dx*dx+sss2)
c	    Gi(j)=Gi(j)-0.5D0/dabs(zg(j)+RL/2.D0)
c     +	   -0.25d0/dsqrt(dx*dx+sss2)
        dx = z0+zg(j)-RL
	    Gp(j)=Gp(j)-0.5d0/dsqrt(dx*dx+sss2)
c	    Gp(j)=Gp(j)-0.5D0/dabs(zg(j)-RL+RL/2.D0)
c     +	   -0.25d0/dsqrt(dx*dx+sss2)
      end do

c********
c first   single layer integrals
c -------
      if (dabs(s1-s0).lt.eps .and. r0.gt.eps) then 
	   xx_tmp1 = 0.0 D0
	   xx_tmp2 = 0.0 D0
	   xx_tmp3 = 0.0 D0
	   xx2 = 0.0D0
	   xx3 = 0.0D0
         
	  do j=1,8
	     Fg(j) = pf1(j) - dlog((1.D0-akg(j))/(sg(j)-s0)**2)*pf2(j)
           ssss = d2 - rg(j) - r0
	     sss2 = ssss**2
         	 dx = z0 - zg(j)
	     dx2 = dx*dx
		 g(j)=fg(j)*Bg(j)- 0.0D0/dabs(zg(j)+RL/2.0D0)
     +        -0.5d0/dsqrt(dx2+sss2)
     	  end do	   

c  regular part        
	  do j=1,8
	     xx_tmp1 = xx_tmp1 + a*amg(j)*g(j)*dsg(j)*rg(j)*w1(j)
	     xx_tmp2 = xx_tmp2 + a*amg(j)*Gi(j)*dsg(j)*rg(j)*w1(j)
	     xx_tmp3 = xx_tmp3 + a*amg(j)*Gp(j)*dsg(j)*rg(j)*w1(j)
	     xx2 = xx2 + amg(j)*Bg(j)*dsg(j)*rg(j)*pf2(j)*w1(j)
	  end do

        xx2 = -xx2*dlog(s2-s1)*a

c  singular part
        do j=1,8
	     ss(j) = (s2-s1)*x2(j) + s1
	     if(n.eq.1)then
	     ams(j) = 1.D0 - x2(j)
	     else if (n.eq.2)then
	     ams(j) = x2(j)
	     end if
	     xd = ss(j)-s1
     	     rs(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	     zs(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	     drs(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	     dzs(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	     dss(j) = dsqrt(drs(j)**2+dzs(j)**2)
	     xh2(j) = (zs(j)-z0)**2
	     rh2(j) = xh2(j) + (rs(j)-r0)**2
	     Rd = dsqrt(xh2(j) + (rs(j)+r0)**2)
	     aks(j) = 4.0D0*r0*rs(j)/Rd**2
           call B0(rs(j),zs(j),r0,z0,Bs(j))
c	  end do
c	  do j=1,8
           call polyn (aks(j),Pff1,Pff2,Pee1,Pee2)
c	     R = sqrt(xh2(j) + (rs(j) + r0)**2)

	     xx3 = xx3 + ams(j)*Bs(j)*dss(j)*rs(j)*Pff2*w2(j)

	  end do
        xx3 = xx3*(s2-s1)

        if (dabs(z0).lt.eps) then
	     y1 = 2.D0*(xx_tmp1+2.0D0*(xx2 + xx3))+xx_tmp3
	   else
	     y1 = xx_tmp1+xx_tmp2+xx_tmp3+2.0D0*(xx2 + xx3)
	  end if

	else if (dabs(s2-s0).lt.eps .and. r0.gt.eps) then
	   xx_tmp1 = 0.0 D0
	   xx_tmp2 = 0.0D0
	   xx_tmp3 = 0.0D0
	   xx2 = 0.0D0
	   xx3 = 0.0D0
         
	  do j=1,8
	     Fg(j) = pf1(j) - dlog((1.D0-akg(j))/(sg(j)-s0)**2)*pf2(j)
           ssss = d2 - rg(j) - r0
	     sss2 = ssss**2
         	 dx = z0 - zg(j)
	     dx2 = dx*dx
		 g(j)=fg(j)*Bg(j)- 0.0D0/dabs(zg(j)+RL/2.0D0)
     +        -0.5d0/dsqrt(dx2+sss2)
        end do	   

c  regular part        
	  do j=1,8
	     xx_tmp1 = xx_tmp1 + a*amg(j)*g(j)*dsg(j)*rg(j)*w1(j)
	     xx_tmp2 = xx_tmp2 + a*amg(j)*Gi(j)*dsg(j)*rg(j)*w1(j)
	     xx_tmp3 = xx_tmp3 + a*amg(j)*Gp(j)*dsg(j)*rg(j)*w1(j)	     
	     xx2 = xx2 + amg(j)*Bg(j)*dsg(j)*rg(j)*pf2(j)*w1(j)
	  end do

        xx2 = -xx2*dlog(s2-s1)*a

c  singular part
        do j=1,8
	     ss(j) = (s1-s2)*x2(j) + s2
	     if(n.eq.1)then
	     ams(j) = x2(j)
	     else if (n.eq.2)then
	     ams(j) = 1.D0-x2(j)
	     end if
	     xd = ss(j)-s1
     	     rs(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	     zs(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	     drs(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	     dzs(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	     dss(j) = dsqrt(drs(j)**2+dzs(j)**2)
	     xh2(j) = (zs(j)-z0)**2
	     rh2(j) = xh2(j) + (rs(j)-r0)**2
	     Rd = dsqrt(xh2(j) + (rs(j)+r0)**2)
	     aks(j) = 4.0D0*r0*rs(j)/Rd**2
           call B0(rs(j),zs(j),r0,z0,Bs(j))
c	  end do
c	  do j=1,8
           call polyn (aks(j),Pff1,Pff2,Pee1,Pee2)
c	     R = sqrt(xh2(j) + (rs(j) + r0)**2)

	     xx3 = xx3 + ams(j)*Bs(j)*dss(j)*rs(j)*Pff2*w2(j)

	  end do
        xx3 = xx3*(s2-s1)

         if(dabs(z0-RL/2.D0).lt.eps) then
	     y1 = 2.D0*(xx_tmp1 + 2.0D0*(xx2 + xx3)) + xx_tmp2
	   else
	     y1 = xx_tmp1 + 2.0D0*(xx2 + xx3) + xx_tmp2 + xx_tmp3
	   end if
	else
	   xx_tmp1 = 0.0D0

	  do j=1,8
	     xx_tmp1=xx_tmp1+
     +       a*amg(j)*(G(j)+Gi(j)+Gp(j))*dsg(j)*rg(j)*w1(j)
c	     xx_tmp1=xx_tmp1+a*amg(j)*(G(j)-Gi(j))*dsg(j)*rg(j)*w1(j)
	  end do

	  y1 = xx_tmp1

	end if
c
c ---- done
	     
	return
	end
c-----------------------------
	subroutine B0(r,z,r0,z0,y)
      Implicit Double Precision (a-h,o-z)
      pi  = 3.1415 92653 58979 32384 D0

	rh = (z-z0)**2 + (r+r0)**2
	y = 1.D0/(pi*dsqrt(rh))

	return 
	end



