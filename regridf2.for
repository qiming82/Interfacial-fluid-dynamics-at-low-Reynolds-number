      subroutine regridf2(znew,rnew,z,r,s,N,time,dt)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
	dimension sold(neq),zold(neq),rold(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension xtmp(neq),ytmp(neq),stmp(neq)
	dimension xtmp0(neq),ytmp0(neq),stmp0(neq)
      dimension dd1(neq),cc1(neq),bb1(neq)
      dimension dd2(neq),cc2(neq),bb2(neq)
	dimension Th(neq)

      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
      RL=snew(N+1) ! total lenght of the arclength
	do k=1,N
       Th(k)=s(k+1)-s(k)       
	end do

	hmin=minval(Th(1:N))
	hmax=maxval(Th(1:N))
       rmin=minval(rnew(1:N+1))

c-----	
c      if(time.gt.0.1.and.N.lt.80)then
c	   N=80
c	   dt=0.6D0*(1.0D0/N)**(1.5D0)
c	end if
c	if(rmin.gt.0.010.and.rmin.lt.0.02.and.N.lt.200)then
c	   N=200
c	   dt=0.4D0*(1.0D0/N)**(1.5D0)
c	endif

c      if(rmin.gt.0.04.or.hmax.gt.1.0D0/N)then
      if(rmin.gt.0.0065)then

         h=RL/N

      	do j=1,N+1
	      s(j)=(j-1)*h
	    end do      

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         end do
      else
         do k=1,N+1
	  s(k)=snew(k)
	  z(k)=znew(k)
	  r(k)=rnew(k)
	  end do
	 
      Ir=0
	 do k=1,N+1
         if(dabs(rmin-r(k)).lt.0.0000000001D0)then
	      Ir=k
	   end if
	 end do
cc--- test value
c	Ir=21
      write(*,*) 'Ir=', Ir
c	 rmin=r(Ir)
	 zmin=z(Ir)
	 smin=s(Ir)

c	do k=1,N
c       Th(k)=s(k+1)-s(k)       
c	end do

c	h=minval(Th(1:N))
	
	if((hmin/rmin).gt.0.1D0.and.N.lt.128)then
c	write(*,*) 'h/rmin=', h/rmin
c	write(*,*) 'add points...'
c-----copy
	  do k=1,N+1
	   sold(k)=s(k)
	   zold(k)=z(k)
	   rold(k)=r(k)
	  end do
       nn=4  ! N=N+2*nn
	
	 do k=1,2*nn+1
         stmp0(k)=s(Ir-nn+k-1)
         xtmp0(k)=z(Ir-nn+k-1)
         ytmp0(k)=r(Ir-nn+k-1)
         dd1(k)=d1(Ir-nn+k-1)
         bb1(k)=b1(Ir-nn+k-1)
         cc1(k)=c1(Ir-nn+k-1)
         dd2(k)=d2(Ir-nn+k-1)
         bb2(k)=b2(Ir-nn+k-1)
         cc2(k)=c2(Ir-nn+k-1)
	end do
	
	RLtmp=dabs(s(Ir+nn)-s(Ir-nn))

      Nm=4*nn
	htmp=RLtmp/Nm

       do k=1,Nm+1
	    stmp(k)=stmp0(1)+(k-1)*htmp
c	    write(*,*) stmp(k)
      	xtmp(k)=seval(Nm,stmp(k),stmp0,xtmp0,bb1,cc1,dd1)
	    ytmp(k)=seval(Nm,stmp(k),stmp0,ytmp0,bb2,cc2,dd2)
       end do
	Nold=N
	N=Nold+2*nn
	write(*,*) 'Nold=', Nold
	write(*,*) 'N=', N
c      end if
	!  from j=1, Ir-nn-1  unchanged
      do j=Ir-nn,Ir-nn+Nm
	   s(j)=stmp(j-Ir+nn+1)
	   z(j)=xtmp(j-Ir+nn+1)
	   r(j)=ytmp(j-Ir+nn+1)
c	   write(*,*) stmp(j-Ir+nn+1)
	end do
      do j=Ir-nn+Nm,N+1
	   s(j)=sold(Ir+nn+(j-Ir+nn-Nm))
	   z(j)=zold(Ir+nn+(j-Ir+nn-Nm))
	   r(j)=rold(Ir+nn+(j-Ir+nn-Nm))
	end do
	dt=0.1D0*(1.0D0/N)**(1.5D0)
c	dt=0.5D0*h**(1.5D0)
	else
c	write(*,*) 'no points added...'
c	write(*,*)  'N=', N

	  do k=1,N+1
	   snew(k)=s(k)
	   znew(k)=z(k)
	   rnew(k)=r(k)
	  end do
         h=RL/N

      	do j=1,N+1
	      s(j)=(j-1)*h
	    end do      

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         end do

      end if
c	dt=0.5D0*hmin**(1.5D0)
	dt=0.1D0*(1.0D0/N)**(1.5D0)
      end if
c	dt=0.5D0*h**(1.5D0)

c	end if
c-------------------------------
       

      return
	end