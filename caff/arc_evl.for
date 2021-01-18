      subroutine arc_evl(z,r,s,N
     +  ,ar,br,cr
     +  ,az,bz,cz)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),s(neq),snew(neq)
	dimension dz(neq),dr(neq),ds(neq),dzg(neq),drg(neq),dsg(neq)
      dimension az(neq),bz(neq),cz(neq),ar(neq),br(neq),cr(neq)
	dimension T(neq),xg(10),wg(10),sg(10),err(neq)

      pi=4*datan(1.d0)
	eps=0.0000001D0
      m=6
	RL=1.D0

      call Gauss_Legendre (m,xg,wg)
      
	errorflag=1.D0

      s(1)=0.0D0
      
	do j=1,N
	   T(j)=dsqrt((z(j+1)-z(j))**2+(r(j+1)-r(j))**2)
         s(j+1)=s(j)+T(j)*2.D0/RL
      end do
      
	kk=1
	do while (errorflag.gt.eps)

          call splc_clm
     +  (N
     +  ,s,r
     +  ,0.0D0
     +  ,0.0D0
     +  ,ar,br,cr
     +  )
          call splc_clm
     +  (N
     +  ,s,z
     +  ,RL/2.D0
     +  ,RL/2.D0
     +  ,az,bz,cz
     +  )

	    cz(N+1)=RL/2.D0
	    cr(N+1)=0.0D0

c	    do j=1,N+1
c	       dr(j)=cr(j)
c	       ddr(j)=2.0D0*br(j)
c	       dz(j)=cz(j)
c	       ddz(j)=2.0D0*bz(j)
c	       ds(j)=dsqrt(dz(j)**2+dr(j)**2)
c	    end do
c --- compute the arclength using the new spline
          snew(1)=0.0D0
	    do j=1,N
	       sum=0.0D0
	       a = (s(j+1)-s(j))/2.D0
	       b = (s(j+1)+s(j))/2.D0
		   do k=1,m
		   	   sg(k) = a*xg(k) + b
			   xd = sg(k) - s(j)
	           drg(k) = 3.D0*ar(j)*xd**2 + 2.D0*br(j)*xd + cr(j)
	           dzg(k) = 3.D0*az(j)*xd**2 + 2.D0*bz(j)*xd + cz(j)
	           dsg(k) = dsqrt(drg(k)**2+dzg(k)**2)
	           sum = sum + a*dsg(k)*wg(k)
	       end do
	       snew(j+1)=snew(j)+sum*2.D0/RL
	     end do
c----  check convergence	
	     do j=1,N+1
	        err(j)=dabs(snew(j)-s(j))
	     end do
	     errorflag=maxval(err)
	     s(1:N+1)=snew(1:N+1)

	    kk = kk +1
	     if(kk.gt.50) then
	         write(*,*)'iteration fails...'
	          exit
	      end if
      end do

          call splc_clm
     +  (N
     +  ,s,r
     +  ,0.0D0
     +  ,0.0D0
     +  ,ar,br,cr
     +  )
          call splc_clm
     +  (N
     +  ,s,z
     +  ,RL/2.D0
     +  ,RL/2.D0
     +  ,az,bz,cz
     +  )

	    cr(N+1)=0.0D0
	    cz(N+1)=RL/2.D0
c---- done

      return
	end