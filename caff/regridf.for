      subroutine regridf(znew,rnew,z,r,s,N,umx,dt)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)

      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
	RL=snew(N+1)
c-----	
c      call spline(N+1,snew,znew,b1,c1,d1)
c	    call spline(N+1,snew,rnew,b2,c2,d2)
c      if(time.gt.0.1 .and. N.lt.100)then
c	  N=100
c	  dt=0.20D0*(1.0D0/N)**(1.5D0)
c	end if
         h=RL/N

      	do j=1,N+1
	      s(j)=(j-1)*h
	    end do      

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         end do
        z(1)=0.d0
         z(N+1)=0.5d0
        dt1=0.003d0/(dabs(umx)+20.d0)
        if(dt1.lt.dt) dt=dt1
      return
	end