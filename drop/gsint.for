      subroutine gsint
     +
     +   (z1,r1
     +   ,s1,s2
     +   ,qx1,qx2,qy1,qy2
     +   ,s0,z0,r0
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,y11,y12
     +   ,y21,y22
     +   )

      Implicit Double Precision (a-h,o-z)
	
	dimension x1(10),x2(10),w1(10),w2(10)
      dimension sg(10),ss(10),qxg(10),qyg(10),rg(10),rs(10),zg(10)
	dimension zs(10),xh2(10),rh2(10),akg(10),aks(10),pf1(10),pf2(10)
	dimension qxs(10),qys(10),pe1(10),pe2(10)
      dimension E11(10), E12(10), E21(10), E22(10)
      dimension G11(10), G12(10), G21(10), G22(10)
      

	eps = 0.000000001D0
      pi  = 3.1415 92653 58979 32384 D0
        m = 4
	call Gauss_Legendre (m,x1,w1)
	call LogGL (m,x2,w2)

	a = (s2-s1)/2.0D0
	b = (s2+s1)/2.0D0
      r02 = r0*r0
c	write (*,*) 'r1=', r1
	do j=1,m
	   sg(j) = a*x1(j) + b
	   qxg(j) = (qx2-qx1)/2.0D0*x1(j) + (qx2+qx1)/2.0D0
	   qyg(j) = (qy2-qy1)/2.0D0*x1(j) + (qy2+qy1)/2.0D0
	   xd = sg(j) - s1
	   rg(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	   zg(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	   xh2(j) = (zg(j)-z0)**2
	   rh2(j) = xh2(j) + (rg(j)-r0)**2
	   akg(j) = 4.0D0*r0*rg(j)/(xh2(j)+(rg(j)+r0)**2)
	end do
      
	iopt = 1
	do j=1,m
	    call polyn (akg(j),pf1(j),pf2(j),pe1(j),pe2(j))
          call sgf_ax_g(iopt,zg(j),rg(j),z0,r0,E11(j),E12(j)
     +		     ,E21(j),E22(j),iaxis)
          call sgf_ax_g(iopt,zg(j),rg(j),-z0,r0,G11(j),G12(j)
     +	         ,G21(j),G22(j),iaxis)
      end do
c -------
      if (dabs(s1-s0).lt.eps .and. r0.gt.eps) then 
	   xx_tmp1 = 0.0 D0
	   xx_tmp2 = 0.0 D0     
	   yy_tmp1 = 0.0 D0
	   yy_tmp2 = 0.0 D0
	   xx2 = 0.0D0
	   yy2 = 0.0D0
	   xx3 = 0.0D0
	   yy3 = 0.0D0
         
	  do j=1,m
           call G_singlr (rg(j),zg(j),r0,z0,sg(j),s0,E11(j),E22(j))
        end do	   


c  regular part        
	  do j=1,m
	     xx_tmp1 = xx_tmp1 + a*E11(j)*qxg(j)*w1(j)
	     xx_tmp2 = xx_tmp2 - a*G11(j)*qxg(j)*w1(j)
	     
	     yy_tmp1 = yy_tmp1 + a*E22(j)*qyg(j)*w1(j)
	     yy_tmp2 = yy_tmp2 + a*G22(j)*qyg(j)*w1(j)
  
           Rd = dsqrt(xh2(j) + (rg(j)+r0)**2)
           r2 = rg(j)*rg(j)
	     xx2 = xx2 - rg(j)/Rd*(pf2(j) +
     +               xh2(j)/rh2(j)*pe2(j))*qxg(j)*w1(j)
	     yy2 = yy2 - 1.0D0/(r0*Rd)*((r02 + 
     +	            r2 + 2.0D0*xh2(j))*pf2(j) -
     +                (2.0D0*xh2(j)**2 + 3.0D0*xh2(j)*(r2+r02) +
     +                 (r02-r2)**2)/rh2(j)*pe2(j))*qyg(j)*w1(j)
	  end do
        xx2 = xx2*dlog(s2-s1)*a*4.0D0
	  yy2 = yy2*dlog(s2-s1)*a*2.0D0

c  singular part
        do j=1,m
	     ss(j) = (s2-s1)*x2(j) + s1
	     qxs(j) = (qx2-qx1)*x2(j) + qx1
	     qys(j) = (qy2-qy1)*x2(j) + qy1
	     xd = ss(j)-s1
     	     rs(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	     zs(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	     xh2(j) = (zs(j)-z0)**2
	     rh2(j) = xh2(j) + (rs(j)-r0)**2
	     Rd = dsqrt(xh2(j) + (rs(j)+r0)**2)
	     aks(j) = 4.0D0*r0*rs(j)/Rd**2
c	  end do
c	  do j=1,m
           call polyn (aks(j),Pff1,Pff2,Pee1,Pee2)
c	     R = sqrt(xh2(j) + (rs(j) + r0)**2)
	     xx3 = xx3 + rs(j)/Rd*(Pff2+
     +	         xh2(j)/rh2(j)*Pee2)*qxs(j)*w2(j)
           yy3 = yy3 + 1.0D0/(r0*Rd)*((r02+rs(j)**2+2.0D0*xh2(j))*Pff2
     +            -(2.0D0*xh2(j)**2 + 3.0D0*xh2(j)*(rs(j)**2+r02)+
     +              (r02-rs(j)**2)**2)/rh2(j)*Pee2)*qys(j)*w2(j)
	  end do
        xx3 = xx3*(s2-s1)*4.0D0
        yy3 = yy3*(s2-s1)*2.0D0

        if (dabs(z0).lt.eps) then
	     y11 = 0.0D0
	     y22 = 2.0D0*(yy_tmp1+2.0D0*(yy2+yy3))
	   else
	     y11 = xx_tmp1 + 2.0D0*(xx2 + xx3) + xx_tmp2
	     y22 = yy_tmp1 + yy_tmp2 + 2.0D0*(yy2 + yy3)
	  end if

	else if (dabs(s2-s0).lt.eps .and. r0.gt.eps) then
	   xx_tmp1 = 0.0 D0
	   xx_tmp2 = 0.0 D0     
	   yy_tmp1 = 0.0 D0
	   yy_tmp2 = 0.0 D0
	   xx2 = 0.0 D0
	   yy2 = 0.0 D0
	   xx3 = 0.0 D0
	   yy3 = 0.0 D0
         
	  do j=1,m
           call G_singlr (rg(j),zg(j),r0,z0,sg(j),s0,E11(j),E22(j))
        end do	   

c  regular part        
	  do j=1,m
	     xx_tmp1 = xx_tmp1 + a*E11(j)*qxg(j)*w1(j)
	     xx_tmp2 = xx_tmp2 - a*G11(j)*qxg(j)*w1(j)
	     
	     yy_tmp1 = yy_tmp1 + a*E22(j)*qyg(j)*w1(j)
	     yy_tmp2 = yy_tmp2 + a*G22(j)*qyg(j)*w1(j)
  
           Rd = dsqrt(xh2(j) + (rg(j)+r0)**2)
           r2 = rg(j)*rg(j)
	     xx2 = xx2 - rg(j)/Rd*(pf2(j) +
     +               xh2(j)/rh2(j)*pe2(j))*qxg(j)*w1(j)
	     yy2 = yy2 - 1.0D0/(r0*Rd)*((r02 + 
     +	            r2 + 2.0D0*xh2(j))*pf2(j) -
     +                (2.0D0*xh2(j)**2 + 3.0D0*xh2(j)*(r2+r02) +
     +                 (r02-r2)**2)/rh2(j)*pe2(j))*qyg(j)*w1(j)
	  end do
        xx2 = xx2*dlog(s2-s1)*a*4.0D0
	  yy2 = yy2*dlog(s2-s1)*a*2.0D0

c  singular part
        do j=1,m
	     ss(j) = (s1-s2)*x2(j) + s2
	     qxs(j) = (qx1-qx2)*x2(j) + qx2
	     qys(j) = (qy1-qy2)*x2(j) + qy2
	     xd = ss(j)-s1
     	     rs(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	     zs(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	     xh2(j) = (zs(j)-z0)**2
	     rh2(j) = xh2(j) + (rs(j)-r0)**2
	     Rd = dsqrt(xh2(j) + (rs(j) + r0)**2)
	     aks(j) = 4.0D0*r0*rs(j)/Rd**2
c	  end do
c	  do j=1,m
           call polyn (aks(j),Pff1,Pff2,Pee1,Pee2)
c	     R = sqrt(xh2(j) + (rs(j) + r0)**2)
	     xx3 = xx3 + rs(j)/Rd*(Pff2+
     +	         xh2(j)/rh2(j)*Pee2)*qxs(j)*w2(j)
           yy3 = yy3 + 1.0D0/(r0*Rd)*((r02+rs(j)**2+2.0D0*xh2(j))*Pff2 
     +            -(2.0D0*xh2(j)**2 + 3.0D0*xh2(j)*(rs(j)**2+r02)+
     +              (r02-rs(j)**2)**2)/rh2(j)*Pee2)*qys(j)*w2(j)
	  end do
        xx3 = xx3*(s2-s1)*4.0D0
        yy3 = yy3*(s2-s1)*2.0D0

	     y11 = xx_tmp1 + 2.0D0*(xx2 + xx3) + xx_tmp2
	     y22 = yy_tmp1 + yy_tmp2 + 2.0D0*(yy2 + yy3)

	else
	   xx_tmp1 = 0.0D0
	   yy_tmp1 = 0.0D0     

	  do j=1,m
	    xx_tmp1 = xx_tmp1 + a*(E11(j)-G11(j))*qxg(j)*w1(j)
	    yy_tmp1 = yy_tmp1 + a*(E22(j)+G22(j))*qyg(j)*w1(j)
	  end do

	  y11 = xx_tmp1
	  y22 = yy_tmp1

	end if

c	sum1 = 0.0D0
	sum2 = 0.0D0
	sum3 = 0.0D0
c	sum4 = 0.0D0

      do j=1,m
c	  sum1 = sum1 + a*(Gxx1-Gxx2)*qxg(j)*w1(j)
	  sum2 = sum2 + a*(E12(j)-G12(j))*qyg(j)*w1(j)
	  sum3 = sum3 + a*(E21(j)+G21(j))*qxg(j)*w1(j)
c	  sum4 = sum4 + a*(Gyy1+Gyy2)*qyg(j)*w1(j)
      end do

c	y11 = y11 + sum1
	y12 = sum2
	y21 = sum3
c	y22 = y22 + sum4

      if (r0.lt.eps) then
	  y22=0.0D0
	end if
      return
	end
c---------------------------------------------------
      subroutine G_singlr(r,z,r0,z0,s,s0,E11,E22)

      Implicit Double Precision (a-h,o-z)
	eps = 0.000000001 D0

	xh = z - z0
	xh2 = xh*xh
	rh2 = xh2 + (r-r0)**2
	ak2 = 4.0D0*r*r0/(xh2+(r+r0)**2)
	r2 = r*r
	r02 = r0*r0

	call polyn(ak2,pf1,pf2,pe1,pe2)
      
c	if (s.eq.s0) then
	  F = pf1 - dlog((1.0D0-ak2)/(s-s0)**2)*pf2
	  E = pe1 - dlog((1.0D0-ak2)/(s-s0)**2)*pe2
c	else
c	  call ell_int(k2,F,E)
c	end if
      
	Rd = dsqrt(xh2+(r+r0)**2)
	E11 = 4.0D0*r/Rd*(F+xh2/rh2*E)
	if (r0.lt.eps) then
      	E22 = 0.0 D0
      else
	    E22 = 2.0D0/(r0*Rd)*((r02+r2+2.0D0*xh2)*F - (2.0D0*xh2**2
     +	     + 3.0D0*xh2*(r2+r02)+(r02-r2)**2)/rh2*E) 
      end if

      return

	end







