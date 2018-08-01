      subroutine gsint
     +
     +   (z1,r1
     +   ,s1,s2
     +   ,f1,f2
     +   ,z0,r0,f0
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,y11,y12
     +   ,y21,y22
     +   )

      Implicit Double Precision (a-h,o-z)
	
	dimension x1(10),x2(10),w1(10),w2(10)
	dimension fg(10),drg(10),dzg(10)
      dimension sg(10),ss(10),rg(10),rs(10),zg(10)
	dimension zs(10),xh2(10),rh2(10),anxg(10),anyg(10)
c	dimension qxs(10),qys(10),pe1(10),pe2(10)
      dimension E11(10), E12(10), E21(10), E22(10)
      dimension G11(10), G12(10), G21(10), G22(10)
      dimension F11(10), F12(10), F21(10), F22(10)
      

c period one
	RL = 1.0D0
	RLH = RL/2.0D0
	eps = 0.000000001D0
      pi  = 3.1415 92653 58979 32384 D0
	subtr = 4.0D0*pi

	call Gauss_Legendre (8,x1,w1)

	a = (s2-s1)/2.0D0
	b = (s2+s1)/2.0D0
      r02 = r0*r0
c	write (*,*) 'r1=', r1
	do j=1,8
	   sg(j) = a*x1(j) + b
	   fg(j) = (f2-f1)/2.0D0*x1(j) + (f2+f1)/2.0D0
c	   qxg(j) = (qx2-qx1)/2.0D0*x1(j) + (qx2+qx1)/2.0D0
c	   qyg(j) = (qy2-qy1)/2.0D0*x1(j) + (qy2+qy1)/2.0D0
	   xd = sg(j) - s1
	   rg(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	   zg(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	   drg(j) = 3.0*Ar*xd**2 + 2.0*Br*xd + Cr
	   dzg(j) = 3.0*Az*xd**2 + 2.0*Bz*xd + Cz
	   anxg(j) = -drg(j)
	   anyg(j) = dzg(j)
c	   xh2(j) = (zg(j)-z0)**2
c	   rh2(j) = xh2(j) + (rg(j)-r0)**2
c	   akg(j) = 4.0D0*r0*rg(j)/(xh2(j)+(rg(j)+r0)**2)
	end do
	iopt = 1

	sum1 = 0.0D0
	sum2 = 0.0D0
	sum3 = 0.0D0
	sum4 = 0.0D0

      do j=1,8
	  call Sgf_axp (rg(j),zg(j),r0,z0,Gxx1,Gxy1,Gyx1,Gyy1)
	  call Sgf_axp (rg(j),zg(j),r0,-z0,Gxx2,Gxy2,Gyx2,Gyy2)
	  sum1 = sum1 + a*(Gxx1-Gxx2)*anxg(j)*(fg(j)-f0)*w1(j)
	  sum2 = sum2 + a*(Gxy1-Gxy2)*anyg(j)*(fg(j)-f0)*w1(j)
	  sum3 = sum3 + a*(Gyx1+Gyx2)*anxg(j)*(fg(j)-f0)*w1(j)
	  sum4 = sum4 + a*(Gyy1+Gyy2)*anyg(j)*(fg(j)-f0)*w1(j)
      end do

	y11 = sum1
	y12 = sum2
	y21 = sum3
	y22 = sum4

      if (z0.lt.eps .or. dabs(z0-RLH).lt.eps) then
	   y11 = 0.0D0
	   y12 = 0.0D0
	end if

      return 
	end
c-------------------------------------------------
      subroutine Sgf_axp (r,z,r0,z0,Gxx,Gxy,Gyx,Gyy)

      Implicit Double Precision (a-h,o-z)
      Iopt = 1
	RL = 1.0D0
	Np=2
	Nsum=5
      call sgf_ax_1p
     +
     +   (Iopt
     +   ,z,r
     +   ,z0,r0
     +   ,RL
     +   ,Nsum,Np
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,QXXX,QXXY,QXYX,QXYY
     +   ,QYXX,QYXY,QYYX,QYYY
     +   ,PXX,PXY,PYX,PYY
     +   ,Iaxis
     +   )
      return
	end


c---------------------------------------------------
