	subroutine regridfadp(znew,rnew,z,r,s,N,wall,dt,unm)
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
       r(1)=0.0d0
       z(N+1)=wall
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
            sdtmp = atan(1.d0/(1.0+dabs(f(j))))
            if(sd(j).lt.sdtmp) sd(j)=sdtmp
c            sd(j)=1.d0/dabs(f(j))
 21       continue
         do 22 j=Ir+1,N+1
           sd(j-1)=atan(dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2))
c           sdtmp=0.5d0/dabs(f(j-1))
c           sd(j-1)=atan((znew(j)-znew(Ir))**2+rnew(j)**2)
c             if(sd(j-1).lt.sdtmp) sd(j-1)=sdtmp
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
       r(1)=0.0d0
       z(N+1)=wall
       dtnew=0.1d0/(5.0+fmax)/(unm+5.d0)
        if(dt.gt.dtnew) dt=dtnew
cc=======================
	if((hmin/rmin).gt.0.1D0/pi.and.N.lt.300.and.hmaxx.lt.smx)then
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
	dtnew=0.02d0*hmin/fmax
       if(dt.gt.dtnew) dt=dtnew
cc==========================
	else
c	write(9,*) 'no points added...'
c	write(9,*)  'N=', N, 'rmin=', rmin, '-'
cc=============================
       dtnew=0.1d0/(5.0+fmax)/(unm+5.d0)
        if(dt.gt.dtnew) dt=dtnew
       endif
      ENDIF
c*************************** check max separation then
c DONE !
      return
	end
c--------------------------------------------------------------------------------
      subroutine regridfadp1(znew,rnew,z,r,s,N,wall,dt,unm)
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
       r(1)=0.0d0
       z(N+1)=wall
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
c           sdtmp=0.5d0/dabs(f(j-1))
c           sd(j-1)=atan((znew(j)-znew(Ir))**2+rnew(j)**2)
c             if(sd(j-1).lt.sdtmp) sd(j-1)=sdtmp
 22         continue
c           hmin=maxval(sd(Ir-3:Ir+3))
c	   hmaxx = maxval(sd(Ir:N))

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
       r(1)=0.0d0
       z(N+1)=wall

       dtnew=0.06d0/(5.0+fmax)/(unm+5.d0)

        if(dt.gt.dtnew) dt=dtnew
      ENDIF
c***************************
c DONE !
      return
	end
c--------------------------------------------------------------------------------      
      subroutine regridf_adp2(znew,rnew,z,r,s,f,Nx,wall,dt,unm)
c-- f curvature to apply adaptivity
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),znew(neq),rnew(neq),snew(neq),s(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
        dimension dg(neq),bg(neq),cg(neq)
        dimension dc(neq),bc(neq),cc(neq)
	dimension gma(neq),gtmp(neq),stmp(neq),sd(neq),f(neq)
c	dimension ctmp(neq),cs(neq),Cb(neq,neq)

      call arc_evl(znew,rnew,snew,Nx
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
	RL=snew(Nx+1)

        pi = 4.0*datan(1.0d0)
c	RL=snew(Nx+1)

       h=RL/Nx
       
       fmax = maxval(abs(f(1:Nx-15)))  ! max of curvature
       Ir=0  ! index of minimum
	 do 1 k=1,Nx-15
         if(abs(fmax-dabs(f(k))).lt.0.0000001D0)then
	     Ir=k
	  endif
 1         continue

c--------------------------------------------------------
c IF
      IF (fmax.le.7.0d0) then ! equal-spacing regridding
       do j=1,Nx+1
	   s(j)=(j-1)*h
	end do      

         do k=1,Nx+1
      	    z(k)=seval(Nx+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(Nx+1,s(k),snew,rnew,b2,c2,d2)
         end do
       r(1)=0.0d0
       z(Nx+1)=wall
	r(Nx+1)=1.0
c----------------------------------------------------------
c ELSE --------------
      ELSEIF (fmax.gt.7.0d0.and.Ir.eq.1)then ! adaptive
	 write(*,*) 'here', Ir, f(1),f(2)
         do 21 j=1,Nx-15
            fk=dabs(f(j+1))
c               sd(j) = 1./(fk+1.0)**1.5
		sd(j) = atan(1.0/(fk+1.0))
c		sd(j) = sqrt((znew(j+1)-znew(1))**2+1.5*rnew(j+1)**2)
 21         continue
	 do j=Nx-14,Nx
	     sd(j) = sd(j-1)
	 enddo
          suml = 0.0d0
          do 3 j=1,Nx
               suml = suml + sd(j)
 3           continue
           hs = RL/suml
           stmp(Nx+1)=RL

           stmp(1)=0.0d0
	  do j=2,Nx
	    stmp(j) = stmp(j-1) + hs*sd(j)
 	  enddo
c-- evaluate output z and r
         s = stmp
         do k=1,Nx+1
      	  z(k)=seval(Nx+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(Nx+1,s(k),snew,rnew,b2,c2,d2)
         end do

       r(1)=0.0d0
       z(Nx+1)=wall
	r(Nx+1)=1.0

c         dtnew=0.001d0/(fmax+1.0)
       dtnew=0.004d0/(5.0+fmax)/(unm+10.d0)
        if(dt.gt.dtnew) dt=dtnew

c-------------------------------------------------------
      ELSE
	write(*,*) 'else here', Ir, f(1), f(2)
        do 31 j=1,Nx-15
           fk=dabs(f(j+1))
c            sd(j) = 1.0/atan(fk+1.0)
            sd(j) = atan(1.0/(fk+1.0))
 31      continue
	 do j=Nx-14,Nx
	     sd(j) = sd(j-1)
	 enddo
c-- change !
        do 32 j=1,Nx
            suml = suml + sd(j)
 32      continue
         hs = RL/suml
         stmp(1)=0.0
         do 34 j=2,Nx
              stmp(j)=stmp(j-1)+hs*sd(j-1)
 34          continue
          stmp(Nx+1)=RL
        s = stmp   
         do k=1,Nx+1
      	  z(k)=seval(Nx+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(Nx+1,s(k),snew,rnew,b2,c2,d2)
         end do

       r(1)=0.0d0
       z(Nx+1)=wall
	r(Nx+1)=1.0

       dtnew=0.002d0/(1.0+fmax)/(unm+10.d0)
c       dtnew=0.0012d0/fmax
      if(dt.gt.dtnew) dt=dtnew
c-- gma and Cb unchanged, so nothing needs to do about them
      Endif
cc-- DONE !!

      Return
      End
