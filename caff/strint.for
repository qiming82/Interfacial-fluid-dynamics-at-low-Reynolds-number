      subroutine strint
     +
     +   (z1,r1
     +   ,s1,s2
     +   ,anx1,anx2,any1,any2
     +   ,z0,r0
     +   ,n
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,y1,y2
     +   ,y3,y4
     +   )

      Implicit Double Precision (a-h,o-z)

	dimension x1(10),w1(10)
      dimension sg(10),rg(10),zg(10),amg(10)
	dimension anzg(10),anrg(10)
	dimension dzg(10),drg(10),dsg(10)
      dimension Qxxx1(10),Qxxy1(10),Qxyx1(10),Qxyy1(10)
      dimension Qyxx1(10),Qyxy1(10),Qyyx1(10),Qyyy1(10)
      dimension Qxxx2(10),Qxxy2(10),Qxyx2(10),Qxyy2(10)
      dimension Qyxx2(10),Qyxy2(10),Qyyx2(10),Qyyy2(10)
      

c period one
	RL = 1.0D0
	eps = 0.000000001D0
      pi  = 3.1415 92653 58979 32384 D0
      m=8
	call Gauss_Legendre (m,x1,w1)

      if (n.eq.1) then
	  do i=1,m
	     amg(i)=(1.0D0-x1(i))/2.0D0
	   end do
	else if (n.eq.2) then
	   do i=1,m
	      amg(i)=(1.0D0+x1(i))/2.0D0
	    end do
	end if

	a = (s2-s1)/2.0D0
	b = (s2+s1)/2.0D0

	do j=1,m
	   sg(j) = a*x1(j) + b
	   anzg(j) = (anx2-anx1)/2.0D0*x1(j) + (anx2+anx1)/2.0D0
	   anrg(j) = (any2-any1)/2.0D0*x1(j) + (any2+any1)/2.0D0
	   xd = sg(j) - s1
	   rg(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	   zg(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
c	   drg(j) = 3.0*Ar*xd**2 + 2.0*Br*xd + Cr
c	   dzg(j) = 3.0*Az*xd**2 + 2.0*Bz*xd + Cz
c	   dsg(j) = dsqrt(dzg(j)**2+drg(j)**2)
	end do
c---   calculate the integral
      yxxtmp1=0.0D0
	yxxtmp2=0.0D0
      yxrtmp1=0.0D0
	yxrtmp2=0.0D0
      yrxtmp1=0.0D0
	yrxtmp2=0.0D0
      yrrtmp1=0.0D0
	yrrtmp2=0.0D0

      do j=1,m
	  call Sgfs_axp(zg(j),rg(j),z0,r0
     +	  ,Qxxx1(j),Qxxy1(j),Qxyx1(j),Qxyy1(j)
     +	  ,Qyxx1(j),Qyxy1(j),Qyyx1(j),Qyyy1(j)
     +	  )
	  call Sgfs_axp(zg(j),rg(j),-z0,r0
     +	  ,Qxxx2(j),Qxxy2(j),Qxyx2(j),Qxyy2(j)
     +	  ,Qyxx2(j),Qyxy2(j),Qyyx2(j),Qyyy2(j)
     +	  )

	  yxxtmp1=yxxtmp1+a*amg(j)*(Qxxx1(j)*anzg(j)+
     +	  Qxxy1(j)*anrg(j))*w1(j)
	  yxxtmp2=yxxtmp2+a*amg(j)*(Qxxx2(j)*anzg(j)+
     +	  Qxxy2(j)*anrg(j))*w1(j)

	  yxrtmp1=yxrtmp1+a*amg(j)*(Qxyx1(j)*anzg(j)+
     +	  Qxyy1(j)*anrg(j))*w1(j)
	  yxrtmp2=yxrtmp2+a*amg(j)*(Qxyx2(j)*anzg(j)+
     +	  Qxyy2(j)*anrg(j))*w1(j)

	  yrxtmp1=yrxtmp1+a*amg(j)*(Qyxx1(j)*anzg(j)+
     +	  Qyxy1(j)*anrg(j))*w1(j)
	  yrxtmp2=yrxtmp2+a*amg(j)*(Qyxx2(j)*anzg(j)+
     +	  Qyxy2(j)*anrg(j))*w1(j)

	  yrrtmp1=yrrtmp1+a*amg(j)*(Qyyx1(j)*anzg(j)+
     +	  Qyyy1(j)*anrg(j))*w1(j)
	  yrrtmp2=yrrtmp2+a*amg(j)*(Qyyx2(j)*anzg(j)+
     +	  Qyyy2(j)*anrg(j))*w1(j)
      end do

	y1 = yxxtmp1-yxxtmp2
	y2 = yxrtmp1-yxrtmp2
	y3 = yrxtmp1+yrxtmp2
	y4 = yrrtmp1+yrrtmp2

c	if (z0.lt.eps .or. dabs(z0-RL/2.0D0).lt.eps) then
	if (z0.lt.eps) then
	  y1 = 0.0D0
	  y2 = 0.0D0
	end if


      return
	end

c-----
      subroutine Sgfs_axp
     +
     +	(X,Y
     +    ,X0,Y0
     +   ,QXXX,QXXY,QXYX,QXYY
     +   ,QYXX,QYXY,QYYX,QYYY
     +    )

      Implicit Double Precision (a-h,o-z)

      Iopt = 2
	RL = 1.0D0
	Np=2
	Nsum=5

      call sgf_ax_1p
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,RL
     +   ,Nsum,Np
     +   ,GXX,GXY
     +   ,GYX,GYY
     +   ,QXXX,QXXY,QXYX,QXYY
     +   ,QYXX,QYXY,QYYX,QYYY
     +   ,PXX,PXY,PYX,PYY
     +   ,Iaxis
     +   )

      return
	end


