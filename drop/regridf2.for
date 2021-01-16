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

      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
      RL=snew(N+1) ! total lenght of the arclength
c-----	
c      if(time.gt.0.1.and.N.lt.80)then
c	   N=80
c	   dt=0.6D0*(1.0D0/N)**(1.5D0)
c	end if
c	if(time.gt.0.20.and.N.lt.160)then
c	   N=160
c	   dt=0.4D0*(1.0D0/N)**(1.5D0)
c	endif
c      if(rmin.lt.0.01 .and. N.lt.200)then
c	   N=160
c	   dt=0.1D0*(1.0D0/N)**(1.5D0)
c	end if


         h=RL/N

      	do j=1,N+1
	      s(j)=(j-1)*h
	    end do      

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         end do

       rmin=minval(r(1:N-20))

c------- adaptive regridding
c-- modification: s in [Ir-nn, Ir+nn]
c      double 2*nn+1 -> 4*nn+1
      if(rmin.lt.0.008)then
c      if(time.gt.0.24)then
c       rmin=minval(r(N/4:3*N/4+10))
	 
      Ir=0
	 do k=1,N-20
         if(dabs(rmin-r(k)).lt.0.0000000001D0)then
	      Ir=k
	   end if
	 end do
cc--- test value
c	Ir=21
c      write(*,*) 'Ir=', Ir
c	 rmin=r(Ir)
	 zmin=z(Ir)
	 smin=s(Ir)
	if((h/rmin).gt.0.1D0.and.N.lt.320)then
c	write(*,*) 'h/rmin=', h/rmin
c	write(*,*) 'add points...'
          call splc_clm
     +  (N
     +  ,s,r
     +  ,0.0D0
     +  ,-pi
     +  ,d2,c2,b2
     +  )
          call splc_clm
     +  (N
     +  ,s,z
     +  ,pi
     +  ,0.0D0
     +  ,d1,c1,b1
     +  )

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
	write(9,*) 'Nold=', Nold
	write(9,*) 'N=', N
c      end if
	!  from j=1, Ir-nn-1  unchanged
      do j=Ir-nn,Ir-nn+Nm
	   s(j)=stmp(j-Ir+nn+1)
	   z(j)=xtmp(j-Ir+nn+1)
	   r(j)=ytmp(j-Ir+nn+1)
c	   write(*,*) stmp(j-Ir+nn+1)
	end do
      do j=Ir-nn+Nm+1,N+1
	   s(j)=sold(Ir+nn+(j-Ir+nn-Nm))
	   z(j)=zold(Ir+nn+(j-Ir+nn-Nm))
	   r(j)=rold(Ir+nn+(j-Ir+nn-Nm))
	end do
	dt=0.1D0*(1.0D0/N)**(1.5D0)
	else
	write(9,*) 'no points added...'
	write(9,*)  'N=', N
      end if
      end if
c-------------------------------
       

      return
	end
