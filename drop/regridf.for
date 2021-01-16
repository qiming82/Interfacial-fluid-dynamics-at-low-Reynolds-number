      subroutine regridf(znew,rnew,z,r,s,N,time,dt)
c
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension coefz(neq,neq),coefr(neq,neq)

      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
	RL=snew(N+1)
c-----	
         h=RL/N

      	do j=1,N+1
	      s(j)=(j-1)*h
	    end do      

         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         end do

	z(1)=0.d0
	r(N+1)=0.d0

      return
	end
c--------------------------------------------------------------------------------
      subroutine regridfadp1(znew,rnew,z,r,s,N,dt,unm)
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
	dimension sold(neq),zold(neq),rold(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension xtmp(neq),ytmp(neq),stmp(neq),sd(neq)
	dimension xtmp0(neq),ytmp0(neq),stmp0(neq)
      dimension dd1(neq),cc1(neq),bb1(neq)
      dimension dd2(neq),cc2(neq),bb2(neq)
      dimension dz(neq),dr(neq),ds(neq),f(neq),ddz(neq),ddr(neq)

      pi=4*datan(1.d0)

c first get the new arclength
      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
c   snew is the current arclength-spacing after advancing pts
      RL=snew(N+1) ! total lenght of the arclength
       rmin=minval(rnew(5:N-10)) ! rmin - before regridding
         h=RL/N
	 smx = 2.0*snew(N+1)
c-- calculate curvature first
	    do j=1,N+1
	       dr(j)=b2(j)
	       ddr(j)=2.0D0*c2(j)
	       dz(j)=b1(j)
	       ddz(j)=2.0D0*c1(j)
	       ds(j)=dsqrt(dz(j)**2+dr(j)**2)
	       call df(rnew(j),dz(j),ddz(j),dr(j),ddr(j),ds(j),f(j))
	    enddo
         fmax = maxval(f(5:N-10))
c IF *****************	
      IF(rmin.gt.0.02d0)then ! equal-spacing regridding
      	 do j=1,N+1
	    s(j)=(j-1)*h
	  enddo      
         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo

	z(1)=0.d0
	r(N+1)=0.d0

       dtnew=0.1d0/fmax/(unm+5.d0)
        if(dt.gt.dtnew) dt=dtnew

c ELSE ***************
       Else
c       write(7,*) 'adaptive'
       Ir=0  ! index of minimum
	 do 1 k=5,N-10
         if(dabs(rmin-rnew(k)).lt.0.0000001D0)then
	      Ir=k
	   endif
 1         continue
c ——————— find local max to define a cap which maintain the grids ———————
	rlmax = maxval(rnew(Ir:N+1))
		Ir2=0
		do k=Ir,N+1
		  if(dabs(rlmax-rnew(k)).lt.0.000001d0)then
		    Ir2=k
		  endif
		enddo
c—————————————————————————————————————
	 do j = Ir2, N
	    sd(j) = snew(j)
	enddo
         do 21 j=1,Ir-1
c            sd(j) = atan(dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2))  
            sd(j) = dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)  
            sdtmp = 0.5d0/dabs(f(j))
            if(sd(j).lt.sdtmp) sd(j)=sdtmp
 21       continue
         do 22 j=Ir+1,Ir2
c           sd(j-1)=dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)
c           sdtmp=0.5d0/dabs(f(j))
c             if(sd(j-1).lt.sdtmp) sd(j-1)=sdtmp
           sd(j-1)=atan(dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2))
            sdtmp = atan(1.d0/(1.0+dabs(f(j-1))))
            if(sd(j-1).lt.sdtmp) sd(j-1)=sdtmp
 22         continue

c 	   do j=N-4,N+1
c	      sd(j-1)= sd(j-2)
c	   enddo

c-- change !
        RLH = snew(Ir)
        suml=0.d0
        do j=1,Ir-1
            suml = suml + sd(j)
        enddo
         hs1 = RLH/suml
 
        suml=0.d0
        do j=Ir,N
            suml = suml + sd(j)
        enddo
         hs2 = (RL-RLH)/suml

         stmp(1)=0.0
          stmp(Ir)=snew(Ir)
         stmp(N+1)=RL
         do j=Ir-1,2,-1
              stmp(j)=stmp(j+1)-hs1*sd(j)
          enddo
         do j=Ir+1,N
             stmp(j)=stmp(j-1)+hs2*sd(j-1)
         enddo
	      hmin=hs1*hmin

        s = stmp   
         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
        end do
      call arc_evl(z,r,s,N

     +  ,d2,c2,b2

     +  ,d1,c1,b1)

	z(1)=0.d0
	r(N+1)=0.d0

       dtnew=0.01d0/(10.0+fmax)/(unm+10.d0)

        if(dt.gt.dtnew) dt=dtnew
      ENDIF
c***************************
c DONE !
      return
	end
c-------------------------------------------------------------------------
	subroutine regridfadp2(znew,rnew,z,r,s,N,dt,unm)
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
	dimension sold(neq),zold(neq),rold(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension xtmp(neq),ytmp(neq),stmp(neq),sd(neq)
	dimension xtmp0(neq),ytmp0(neq),stmp0(neq)
      dimension dd1(neq),cc1(neq),bb1(neq)
      dimension dd2(neq),cc2(neq),bb2(neq)
      dimension dz(neq),dr(neq),ds(neq),f(neq),ddz(neq),ddr(neq)

      pi=4*datan(1.d0)

c first get the new arclength
      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
c   snew is the current arclength-spacing after advancing pts
      RL=snew(N+1) ! total lenght of the arclength
       rmin=minval(rnew(5:N-10)) ! rmin - before regridding
         h=RL/N
	 smx = 2.0*snew(N+1)
c-- calculate curvature first
	    do j=1,N+1
	       dr(j)=b2(j)
	       ddr(j)=2.0D0*c2(j)
	       dz(j)=b1(j)
	       ddz(j)=2.0D0*c1(j)
	       ds(j)=dsqrt(dz(j)**2+dr(j)**2)
	       call df(rnew(j),dz(j),ddz(j),dr(j),ddr(j),ds(j),f(j))
	    enddo
         fmax = maxval(f(5:N-10))
c IF *****************	
      IF(rmin.gt.0.025d0)then ! equal-spacing regridding
      	 do j=1,N+1
	    s(j)=(j-1)*h
	  enddo      
         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
	z(1)=0.d0
	r(N+1)=0.d0
       dtnew=0.1d0/fmax/(unm+2.d0)
        if(dt.gt.dtnew) dt=dtnew

c ELSE ***************
       Else
c       write(7,*) 'adaptive'
       Ir=0  ! index of minimum
	 do 1 k=5,N-10
         if(dabs(rmin-rnew(k)).lt.0.0000001D0)then
	      Ir=k
	   endif
 1         continue
         do 21 j=1,Ir-1
            sd(j) = atan(dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2))  
c            sdtmp = atan(1.d0/(1.0+dabs(f(j))))
c            if(sd(j).lt.sdtmp) sd(j)=sdtmp
c            sd(j)=1.d0/dabs(f(j))
 21       continue
         do 22 j=Ir+1,N+1
           sd(j-1)=atan(dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2))
           sdtmp=0.5d0/dabs(f(j-1))
c           sd(j-1)=atan((znew(j)-znew(Ir))**2+rnew(j)**2)
             if(sd(j-1).lt.sdtmp) sd(j-1)=sdtmp
 22         continue
           hmin=maxval(sd(Ir-3:Ir+3))
	   hmaxx = maxval(sd(Ir:N))

c-- change !
        RLH = snew(Ir)
        suml=0.d0
        do j=1,Ir-1
            suml = suml + sd(j)
        enddo
         hs1 = RLH/suml
 
        suml=0.d0
        do j=Ir,N
            suml = suml + sd(j)
        enddo
         hs2 = (RL-RLH)/suml

         stmp(1)=0.0
          stmp(Ir)=snew(Ir)
         stmp(N+1)=RL
         do j=Ir-1,2,-1
              stmp(j)=stmp(j+1)-hs1*sd(j)
          enddo
         do j=Ir+1,N
             stmp(j)=stmp(j-1)+hs2*sd(j-1)
         enddo
	      hmin=hs1*hmin

        s = stmp   
         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
        end do
      call arc_evl(z,r,s,N

     +  ,d2,c2,b2

     +  ,d1,c1,b1)

	z(1)=0.d0
	r(N+1)=0.d0
       dtnew=0.1d0/(5.0+fmax)/(unm+5.d0)
        if(dt.gt.dtnew) dt=dtnew
cc=======================
	if((hmin/rmin).gt.0.1D0/pi.and.N.lt.350.and.hmaxx.lt.smx)then
c-- modification: s in [Ir-nn, Ir+nn]
c   double 2*nn+1-->4*nn+1
c-----copy
	  do 5 k=1,N+1
	   sold(k)=s(k)
	   zold(k)=z(k)
	   rold(k)=r(k)
 5        continue
       nn=4  ! N=N+2*nn
	
	 do 6 k=1,2*nn+1
         stmp0(k)=s(Ir-nn+k-1)
         xtmp0(k)=z(Ir-nn+k-1)
         ytmp0(k)=r(Ir-nn+k-1)
         dd1(k)=d1(Ir-nn+k-1)
         bb1(k)=b1(Ir-nn+k-1)
         cc1(k)=c1(Ir-nn+k-1)
         dd2(k)=d2(Ir-nn+k-1)
         bb2(k)=b2(Ir-nn+k-1)
         cc2(k)=c2(Ir-nn+k-1)
 6       continue
	
	RLtmp=dabs(s(Ir+nn)-s(Ir-nn))

      Nm=4*nn
	htmp=RLtmp/Nm

       do k=1,Nm+1
	    stmp(k)=stmp0(1)+(k-1)*htmp
      	xtmp(k)=seval(Nm+1,stmp(k),stmp0,xtmp0,bb1,cc1,dd1)
	    ytmp(k)=seval(Nm+1,stmp(k),stmp0,ytmp0,bb2,cc2,dd2)
       end do
	Nold=N
	N=Nold+2*nn
c	write(9,*) 'Nold=', Nold
c	write(9,*) 'N=', N, 'rmin=', rmin, '+'
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
      call arc_evl(z,r,s,N

     +  ,d2,c2,b2

     +  ,d1,c1,b1)
      
	z(1)=0.d0
	r(N+1)=0.d0

c	dtnew=0.02d0*hmin/fmax
c       if(dt.gt.dtnew) dt=dtnew
cc==========================
	else
	write(9,*) 'no points added...'
	write(9,*)  'N=', N, 'hovr=', hmin/rmin, '-'
cc=============================
       dtnew=0.1d0/(5.0+fmax)/(unm+5.d0)
        if(dt.gt.dtnew) dt=dtnew
       endif
      ENDIF
c*************************** check max separation then
c DONE !
      return
	end
c--------------------------------------------------------------------
      subroutine regridfadp3(znew,rnew,z,r,s,N,dt,unm)
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
	dimension sold(neq),zold(neq),rold(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
	dimension xtmp(neq),ytmp(neq),stmp(neq),sd(neq)
	dimension xtmp0(neq),ytmp0(neq),stmp0(neq)
      dimension dd1(neq),cc1(neq),bb1(neq)
      dimension dd2(neq),cc2(neq),bb2(neq)
      dimension dz(neq),dr(neq),ds(neq),f(neq),ddz(neq),ddr(neq)

      pi=4*datan(1.d0)

c first get the new arclength
      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
c   snew is the current arclength-spacing after advancing pts
      RL=snew(N+1) ! total lenght of the arclength
       rmin=minval(rnew(1:N-10)) ! rmin - before regridding
         h=RL/N
c-- calculate curvature first
	    do j=1,N+1
	       dr(j)=b2(j)
	       ddr(j)=2.0D0*c2(j)
	       dz(j)=b1(j)
	       ddz(j)=2.0D0*c1(j)
	       ds(j)=dsqrt(dz(j)**2+dr(j)**2)
	       call df(rnew(j),dz(j),ddz(j),dr(j),ddr(j),ds(j),f(j))
	    end do
         fmax = maxval(f(1:N+1))
c IF *****************	
      if(rmin.ge.0.02d0)then ! equal-spacing regridding
      	 do j=1,N+1
	    s(j)=(j-1)*h
	  enddo      
         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
c ELSE ***************
       else !------- adaptive regridding
c local spacing propotional to the dist from R_min
       Ir=0  ! index of minimum
	 do 1 k=1,N-10
         if(dabs(rmin-rnew(k)).lt.0.0000001D0)then
	      Ir=k
	   endif
 1         continue
        RLH = snew(Ir)
         
c distance from the minimum location
         do 21 j=1,Ir-1
            sd(j) = atan(dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2))  
            sdtmp = atan(1.d0/(1.0+dabs(f(j))))
            if(sd(j).lt.sdtmp) sd(j)=sdtmp
 21       continue
         do 22 j=Ir+1,N+1
           sd(j-1)=atan(dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2))
           sdtmp= atan(0.1d0/(1.0+dabs(f(j))))
             if(sd(j-1).lt.sdtmp) sd(j-1)=sdtmp
 22         continue
c-- change!
c         hmin=minval(sd(1:N))
           hmin=maxval(sd(Ir-1:Ir+6))
c           hmin=sd(Ir)
c         if(sd(Ir-1).le.sd(Ir)) hmin=sd(Ir-1)
c step constant
            suml=0.0d0
          do 31 j=1,Ir-1
              suml=suml+sd(j)
 31        continue
           hs1 = RLH/suml
            suml = 0.0d0
          do 32 j=Ir,N
              suml=suml+sd(j)
 32        continue
           hs2 = (RL-RLH)/suml
          stmp(1)=0.0d0
c          do 4 k=2,N
cc             stmp(k)=stmp(k-1)+hs*atan(0.7d0*sd(k))
c             stmp(k) = stmp(k-1) + hs*sd(k-1)
c 4       continue
          stmp(Ir)=snew(Ir)
           do k=Ir-1,2,-1
              stmp(k)=stmp(k+1)-hs1*sd(k)
           enddo
           do k=Ir+1,N
              stmp(k)=stmp(k-1) + hs2*sd(k-1)
           enddo
          stmp(N+1)=RL
	      hmin=hs1*hmin
c-- evaluate output z and r
          s=stmp
         do k=1,N+1
      	  z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	      r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo
      call arc_evl(z,r,s,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
       dtnew=0.1d0/(5.0+fmax)/(unm+5.d0)
        if(dt.gt.dtnew) dt=dtnew
cc=======================
	if((hmin/rmin).gt.0.25D0/pi.and.N.lt.400)then
c-- modification: s in [Ir-nn, Ir+nn]
c   double 2*nn+1-->4*nn+1
c-----copy
	  do 5 k=1,N+1
	   sold(k)=s(k)
	   zold(k)=z(k)
	   rold(k)=r(k)
 5        continue
       nn=2  ! N=N+2*nn
	
	 do 6 k=1,2*nn+1
         stmp0(k)=s(Ir-nn+k-1)
         xtmp0(k)=z(Ir-nn+k-1)
         ytmp0(k)=r(Ir-nn+k-1)
         dd1(k)=d1(Ir-nn+k-1)
         bb1(k)=b1(Ir-nn+k-1)
         cc1(k)=c1(Ir-nn+k-1)
         dd2(k)=d2(Ir-nn+k-1)
         bb2(k)=b2(Ir-nn+k-1)
         cc2(k)=c2(Ir-nn+k-1)
 6       continue
	
	RLtmp=dabs(s(Ir+nn)-s(Ir-nn))

      Nm=4*nn
	htmp=RLtmp/Nm

       do k=1,Nm+1
	    stmp(k)=stmp0(1)+(k-1)*htmp
      	xtmp(k)=seval(Nm,stmp(k),stmp0,xtmp0,bb1,cc1,dd1)
	    ytmp(k)=seval(Nm,stmp(k),stmp0,ytmp0,bb2,cc2,dd2)
       end do
	Nold=N
	N=Nold+2*nn
	write(9,*) 'Nold=', Nold
	write(9,*) 'N=', N, 'rmin=', rmin, '+'
        write(9,*) '--------------------------'
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
	dtnew=0.1d0*hmin/dabs(fmax)
       if(dt.gt.dtnew) dt=dtnew
      call arc_evl(z,r,s,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
cc==========================
	else
	write(9,*) 'no points added...'
	write(9,*)  'N=', N, 'rmin=', rmin, '-'
        write(9,*) '---------------------------'
cc=============================
c	dtnew=0.05d0*rmin/fmax
c       dtnew = 2.d0*rmin/(dabs(fmax)+1.0)
       dtnew=0.1d0/(5.0+fmax)/(unm+5.d0)
       if(dt.gt.dtnew) dt=dtnew
      end if
c ENDIF *****************
      end if
c-------------------------------       
      return
	end
