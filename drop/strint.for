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

	dimension x1(20),w1(20)
      dimension sg(20),rg(20),zg(20),amg(20)
	dimension anzg(20),anrg(20)
      dimension Qxxx1(20),Qxxy1(20),Qxyx1(20),Qxyy1(20)
      dimension Qyxx1(20),Qyxy1(20),Qyyx1(20),Qyyy1(20)
      dimension Qxxx2(20),Qxxy2(20),Qxyx2(20),Qxyy2(20)
      dimension Qyxx2(20),Qyxy2(20),Qyyx2(20),Qyyy2(20)
      
	eps = 0.000000001D0
      pi  = 3.1415 92653 58979 32384 D0
         m = 4
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

	Iopt=2

      do j=1,m
	  call sgf_ax_str     
     +
     + (Iopt
     + ,zg(j),rg(j)
     + ,z0,r0
     + ,Qxxx1(j),Qxxy1(j),Qxyx1(j),Qxyy1(j)
     + ,Qyxx1(j),Qyxy1(j),Qyyx1(j),Qyyy1(j)
     + ,Iaxis
     + )
	  call sgf_ax_str     
     +
     + (Iopt
     + ,zg(j),rg(j)
     + ,-z0,r0
     + ,Qxxx2(j),Qxxy2(j),Qxyx2(j),Qxyy2(j)
     + ,Qyxx2(j),Qyxy2(j),Qyyx2(j),Qyyy2(j)
     + ,Iaxis
     + )


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

	if (dabs(z0).lt.eps) then
	  y1 = 0.0D0
	  y2 = 0.0D0
	end if

      return
	end

