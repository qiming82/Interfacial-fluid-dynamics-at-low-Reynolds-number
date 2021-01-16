      subroutine regridf_adp(znew,rnew,z,r,s,N,time,dt)
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

c first get the new arclength
      call arc_evl(znew,rnew,snew,N
     +  ,d2,c2,b2
     +  ,d1,c1,b1)
c   snew is the current arclength-spacing after advancing pts
      RL=snew(N+1) ! total lenght of the arclength
       rmin=minval(rnew(1:N-10)) ! rmin - before regridding
         h=RL/N
c IF *****************	
      if(rmin.ge.0.01d0)then ! equal-spacing regridding
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
         
c distance from the minimum location
         do 21 j=1,Ir-1
            sd(j) = dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)  
 21       continue
         do 22 j=Ir+1,N+1
           sd(j-1)=dsqrt((znew(j)-znew(Ir))**2+rnew(j)**2)
 22         continue
         hmin=minval(sd(1:N))
c step constant
            suml=0.0d0
          do 3 j=1,N
             suml=suml+atan(0.5d0*sd(j))
 3          continue
           hs = RL/suml
          stmp(1)=0.0
          do 4 k=2,N
             stmp(k)=stmp(k-1)+hs*atan(0.5d0*sd(k))
 4       continue
          stmp(N+1)=RL
c-- evaluate output z and r
        s=stmp
         do k=1,N+1
      	   z(k)=seval(N+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(N+1,s(k),snew,rnew,b2,c2,d2)
         enddo

cc=======================
	if((hmin/rmin).gt.0.1D0.and.N.lt.400)then
c-- modification: s in [Ir-nn, Ir+nn]
c   double 2*nn+1-->4*nn+1
c-----copy
	  do 5 k=1,N+1
	   sold(k)=s(k)
	   zold(k)=z(k)
	   rold(k)=r(k)
 5        continue

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

	    b1(N+1)=0.0D0
	    b2(N+1)=-pi

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
	dt=0.2D0*hmin/2.0
cc==========================
	else
	write(9,*) 'no points added...'
	write(9,*)  'N=', N
cc=============================
      end if
c ENDIF *****************
      end if
c-------------------------------
       

      return
	end
