      subroutine regridf_sft2(znew,rnew,snew,z,r,s
     +  ,d2,c2,b2
     +  ,d1,c1,b1
     +  ,N,gma,time,dt)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024*2)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension coefz(neq,neq),coefr(neq,neq)
	dimension gma(neq),gtmp(neq)

	RL=snew(N+1)

      gtmp=gma

      h=RL/N

      do j=1,N+1
	   s(j)=(j-1)*h
	end do      

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         end do
c----
      call spline(n+1,snew,gtmp,b1,c1,d1)

         do k=1,N+1
      	  gma(k)=seval(N+1,s(k),snew,gtmp,b1,c1,d1)
         end do

      return
	end
