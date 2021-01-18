      subroutine regridf2(znew,rnew,z,r,s,N,umx,dt)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
      dimension stmp(neq),sd(neq)
      dimension dr(neq),dz(neq),ddr(neq),ddz(neq),ds(neq),f(neq)

      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
	RL=snew(N+1)
      rmin = minval(rnew(1:N+1))
c-----	
       if(rmin.gt.0.012d0)then
        h=RL/N
      	do j=1,N+1
	      s(j)=(j-1)*h
	 enddo      
         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
       z(1)=0.0d0
       z(N+1)=0.5d0
      dtnew = 0.01d0*rmin
      if(dt.gt.dtnew) dt=dtnew
c--
      else
	    do j=1,N+1
	       dr(j)=b2(j)
	       ddr(j)=2.0D0*c2(j)
	       dz(j)=b1(j)
	       ddz(j)=2.0D0*c1(j)
	       ds(j)=dsqrt(dz(j)**2+dr(j)**2)
	       call df(r(j),dz(j),ddz(j),dr(j),ddr(j),ds(j),f(j))
           enddo       
       fmax = maxval(f(1:N+1))  ! max of curvature -- before regridding
c--
       Ir=0  ! index of minimum
	 do 1 k=1,N+1
         if(dabs(rmin-rnew(k)).lt.0.0000001D0)then
	      Ir=k
	   endif
 1         continue
        RLH = snew(Ir)

         do j=1,Ir-1
            sd(j) = dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)
        enddo
         do j=Ir+1,N+1
            sd(j-1) = dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)
        enddo
      suml=0.0d0
       do j=1,Ir-1
         suml=suml+sd(j)
      enddo
       hs1=RLH/suml
      suml=0.0d0
       do j=Ir,N
         suml=suml+sd(j)
      enddo
       hs2 = (RL-RLH)/suml
      stmp(1)=0.0d0
       stmp(Ir)=snew(Ir)
       do k=Ir-1,2,-1
           stmp(k)=stmp(k+1)-hs1*sd(k)
        enddo
       do j=Ir+1,N
         stmp(j)=stmp(j-1)+hs2*sd(j-1)
       enddo
      stmp(N+1)=RL
       s=stmp

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
       z(1)=0.0
       z(N+1)=0.5
      call arc_evl(z,r,s,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
      dtnew = 0.005d0*rmin/(fmax+umx)
      if(dt.gt.dtnew) dt=dtnew
      endif
c done
      return
	end
c---------------------------------------------------------------------
      subroutine initl(z,r,N1)
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)
       dimension z(neq),r(neq),q(neq)
       open(1,file='xy0.dat',status='old')
c       open(2,file='charge0.dat',status='old')
       do k=1,N1
          read(1,*) z(k), r(k)
c          read(2,*) z(k), q(k)
        enddo
      return
       end
