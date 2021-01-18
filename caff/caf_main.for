      program vjet_axisym

c      implicit real*8 (a-h,o-y)
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),z0(neq),r0(neq)
	dimension s(neq),f(neq),fp(neq),hh(neq)
	dimension zp(neq),rp(neq),znew(neq),rnew(neq),snew(neq),v(neq)
	dimension dz(neq),ddz(neq),dr(neq),ddr(neq),ds(neq),indx(neq)
	dimension dzp(neq),ddzp(neq),drp(neq),ddrp(neq),dsp(neq)
	dimension u1(neq),u2(neq),u1p(neq),u2p(neq),un(neq),unp(neq)
      dimension an1(neq),an2(neq),an1p(neq),an2p(neq),ut(neq),utp(neq)
      dimension anx(neq),anr(neq),anxp(neq),anrp(neq)
      dimension q1(neq),q2(neq),q1p(neq),q2p(neq),bb(neq),bbp(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
      dimension az(neq),bz(neq),cz(neq),ar(neq),br(neq),cr(neq)
      dimension ml(neq,neq), mu(neq,neq),Am(neq,neq)
      dimension bbz(neq),bbr(neq),bbzp(neq),bbrp(neq),w_trp(neq)
	dimension Sxx(neq,neq),Sxr(neq,neq),Srx(neq,neq),Srr(neq,neq)
	dimension Sxxp(neq,neq),Sxrp(neq,neq),Srxp(neq,neq),Srrp(neq,neq)
      dimension u_sxxtmp1(neq),u_sxrtmp1(neq),u_srxtmp1(neq)
	dimension u_srrtmp1(neq)
      dimension u_sxxtmp2(neq),u_sxrtmp2(neq),u_srxtmp2(neq)
	dimension u_srrtmp2(neq)
	dimension en(neq),enp(neq),ipiv(neq),fc(neq)

      pi=4*datan(1.d0)
	eps=0.000000001D0
	RL=1.0D0
	RLH=RL/2.0d0
      wk=2.0D0*pi
	a=(0.5d0)/wk
c      a=0.05d0/RL  !eps=1/10
	a1=a*0.2D0
	d=a*1.5D0  ! finite tube radius
c      d=0.0d0   ! infinite thread
c	a2=a*0.2d0

c number of segments of the interface
      N=128
	Nm=2*N+2

c    weights of trap rule
      do i=1,N+1
	   w_trp(i)=1.0D0
	end do
	w_trp(1)=0.5D0
	w_trp(N+1)=0.5D0

      h=1.0D0/N
      
	nplt = 160
c   prepare to store the data
      open (22,file='axialv.dat')
      open (23,file='xyval.dat')
      open (24,file='normalforce.dat')
      open (11,file='hmin.dat')
      open (12,file='wmax.dat')
      open (13,file='hmax.dat')
      open (14,file='en.dat')
      open (7,file='ampltest.dat')
      open (8,file='plottime.dat')
      open (30,file='info.txt')


      sum11=0.0D0
	sum12=0.0D0
	sum21=0.0D0
	sum22=0.0D0
      y11 = 0.0D0
	y12 = 0.0D0
	y21 = 0.0D0
	y22 = 0.0D0

	dt=0.0008D0
      dt=0.1d0*(1.D0/N)**(1.5D0)
	dt2=dt/2.0d0
	dt3=dt2/2.0d0
	dt4=dt3/2.0d0
	tf=3.0D0
	nplt=tf/dt/50
      nplt2=nplt/5

c     nplt=10
c initial settings...
        iopt = 0
        if(iopt.eq.0)then
	do j=1,N+1
	   s(j)=(j-1)*h
	   z0(j)=RLH*s(j)
	   z(j)=z0(j)
	   r0(j)=a-a1*dcos(wk/RL*z(j)) 
	   r(j)=r0(j)
	   write(23,*) z0(j), r0(j)
	end do      
        else
           N=80
            call initl(z0,r0,N+1)
             do j=1,N+1
              z(j) = z0(j)
              r(j) = r0(j)
              write(23,*) z0(j),r0(j)
             enddo
           call arc_evl(z,r,s,N,d2,c2,b2,d1,c1,d1)
            dt = 0.0001d0
        endif
c
         cvol0=0.0D0
c       get the volume of the jet of one period,
          do j=1,N
       	   hh(j)=z(j+1)-z(j)
	       cvol0=cvol0+hh(j)*(r(j)**2+r(j+1)**2)
	     end do
	     cvol0=cvol0*pi/2.0d0

      time=0.0D0
	write(8,*) time, vcol0

      ram=1.D0
c--- viscous Atwood #
      beta1=(ram-1.0D0)/(ram+1.0D0)
	beta1=beta1/(4.0D0*pi) ! at the RHS

c      beta2=-1.D0/(4.D0*pi*(1.D0+ram))
      beta2=-1.D0/2.D0

c	eb=0.d0*a*dlog(d/a)*dlog(d/a)
       eb = 0.d0*a

      jj=1
      jjj=1
c	write(*,*) h
c start evolution
      do while (time.lt.tf)

	   time = time + dt
c--------

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
     +  ,RLH
     +  ,RLH
     +  ,az,bz,cz
     +  )

	    cz(N+1)=RLH
	    cr(N+1)=0.0D0
c--- calculate the electric force
      if(eb.ne.0.0D0)then
           call geten
     +
     +   (z,r
     +   ,az,bz,cz
     +   ,ar,br,cr
     +   ,N,d
     +   ,s	
     +   ,en
     +   )
      else
	  do kk=1,N+1
	    en(kk)=0.0D0
	  end do
      end if

	    do j=1,N+1
	       dr(j)=cr(j)
	       ddr(j)=2.0D0*br(j)
	       dz(j)=cz(j)
	       ddz(j)=2.0D0*bz(j)

	       ds(j)=dsqrt(dz(j)**2+dr(j)**2)
             an1(j)=-dr(j)/ds(j)
	       an2(j)=dz(j)/ds(j)
	       anx(j)=-dr(j)
	       anr(j)=dz(j)
	       call df(r(j),dz(j),ddz(j),dr(j),ddr(j),ds(j),fc(j))
	       f(j) = fc(j) - 0.5D0*eb*(en(j)**2)
	    end do
	    
c      prepare the integration of the kernels...

c      step 1, single layers

	    do i=1,N+1
	        sum11=0.0D0
			sum12=0.0D0
			sum21=0.0D0
		    sum22=0.0D0
	       do k=1,N
c	    write(*,*) 'd=', d
	    if(d.lt.a)then
c	    write(*,*) 'infinite thread'
	         call gsint
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,f(k),f(k+1)
     +   ,z(i),r(i),f(i)
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,y11,y12
     +   ,y21,y22
     +   )
	     else
c	   write(*,*) 'finite tube!'
	         call gsint_ct
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,f(k),f(k+1)
     +   ,z(i),r(i),f(i),d
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,y11,y12
     +   ,y21,y22
     +   )
	    endif
	       sum11 = sum11 + y11
	       sum12 = sum12 + y12
	       sum21 = sum21 + y21
	       sum22 = sum22 + y22
	      end do
c  u1 = ux,  u2 = ur
	      bbz(i) = sum11 + sum12
	      bbr(i) = sum21 + sum22
	    end do

	if(dabs(ram-1.0D0).lt.eps)then
         do k=1,N+1
	     u1(k)=bbz(k)*beta2
	     u2(k)=bbr(k)*beta2
	  end do
	else
c bb is the right hand side of the system
		do k=1,N+1
            bb(k)=bbz(k)*beta2
	      bb(k+N+1)=bbr(k)*beta2  
		end do

c       step 2, double layer potentials...

         do i=1,N+1
	   
	       do k=1,N
          call strint
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,anx(k),anx(k+1),anr(k),anr(k+1)
     +   ,z(i),r(i)
     +   ,1
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,anx(k),anx(k+1),anr(k),anr(k+1)
     +   ,z(i),r(i)
     +   ,2
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )

	   if (k.eq.1) then
	      Sxx(i,k)=u_sxxtmp1(k)
	      Sxr(i,k)=u_sxrtmp1(k)
	      Srx(i,k)=u_srxtmp1(k)
	      Srr(i,k)=u_srrtmp1(k)
c	   endif
c	   if (k.gt.1) then 
         else
	      Sxx(i,k)=u_sxxtmp1(k)+u_sxxtmp2(k-1)
	      Sxr(i,k)=u_sxrtmp1(k)+u_sxrtmp2(k-1)
	      Srx(i,k)=u_srxtmp1(k)+u_srxtmp2(k-1)
	      Srr(i,k)=u_srrtmp1(k)+u_srrtmp2(k-1)
	   end if
	   end do

	    Sxx(i,N+1)=u_sxxtmp2(N)
	    Sxr(i,N+1)=u_sxrtmp2(N)
	    Srx(i,N+1)=u_srxtmp2(N)
	    Srr(i,N+1)=u_srrtmp2(N)
	   
	  end do
c        do i=1,N
c	    write (*,*) u_srxtmp2(i)
c	   end do
c	  write (*,*) Srx(N+1,N+1)
        Sxx=-beta1*Sxx
	  Sxr=-beta1*Sxr
	  Srx=-beta1*Srx
	  Srr=-beta1*Srr

	  do i=1,N+1
	    Sxx(i,i)=1.0D0+Sxx(i,i)
	    Srr(i,i)=1.0D0+Srr(i,i)
c	    write (*,*) Sxx(i,i)
	  end do
c	 exit

c  Am=[Sxx -Sxr; -Srx Srr]
        do i=1,N+1
	    do k=1,N+1
	      Am(i,k)=Sxx(i,k)
	      Am(i,k+N+1)=Sxr(i,k)
	      Am(i+N+1,k)=Srx(i,k)
	      Am(i+N+1,k+N+1)=Srr(i,k)
	    end do
	  end do

         call dgesv(Nm, 1, Am, neq, ipiv, bb, neq, info)
	   v(1:Nm)=bb(1:Nm)       
	   u1(1:N+1)=v(1:N+1)
	   u2(1:N+1)=v(N+2:2*N+2)


	end if


	    do k=1,N+1
	       un(k)=u1(k)*an1(k)+u2(k)*an2(k)
	       ut(k)=u1(k)*an2(k)-u2(k)*an1(k)
c	       write (*,*) un(k)
c	       zp(k)=z(k)+dt*un(k)*an1(k)
c	       rp(k)=r(k)+dt*un(k)*an2(k)
	       znew(k)=z(k)+dt*(un(k)*an1(k)+ut(k)*an2(k))
	       rnew(k)=r(k)+dt*(un(k)*an2(k)-ut(k)*an1(k))
c	       write(23,*) znew(k),rnew(k)
	     end do

c---------------------------

       unm=maxval(un)


c       call regridf2(znew,rnew,z,r,s,N,unm,dt)
       call regridf(znew,rnew,z,r,s,N,unm,dt)
c      z(1:n+1)=znew(1:n+1)
c	r(1:n+1)=rnew(1:n+1)
c	s(1:n+1)=2.0/RL*z(1:n+1)

         cvol=0.0D0
c       get the volume of the jet of one period,
          do j=1,N
       	   hh(j)=z(j+1)-z(j)
	       cvol=cvol+hh(j)*(r(j)**2+r(j+1)**2)
	     end do
	     cvol=cvol*pi/2.0
c	   write (*,*) 'ratio=', cvol/cvol0

c scale the jet
c         vscale=(cvol/cvol0)**(1.D0/3.D0)
c          do j=1,N+1
cc	       z(j)=z(j)/vscale
c	       r(j)=r(j)/vscale
c	    end do


c----------
         vnmax = maxval(un)
         v1max = maxval(dabs(u1))

         rmin = minval(r(1:N+1))
	   rmax = maxval(r(1:N+1))

	   if (mod(jj,nplt2).eq.0) then
	      write (11,*) time, rmin
	      write (12,*) time, v1max
	      write (13,*) time, rmax
	   amp=(rmax-rmin)/2.0
	      write(7,*) time, dlog(amp/a1)
         endif

	   if (mod(jj,nplt).eq.0) then
	      jjj=jjj+1
	      write(8,*) time, cvol
c	      do j=1,N+1
c	          write(23,*) z(j), r(j)
c	       enddo
	     if(mod(jjj,2).eq.0)then
	        do j=N+1,1,-1
		      write (14,*) z(j), en(j)**2	 
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (24,*) z(j), f(j)
	        end do
	     else
	        do j=1,N+1
		      write (14,*) z(j), en(j)**2	 
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (24,*) z(j), f(j)
	        end do
	     end if
          end if	      

c	   write(*,*) 'rmin=', rmin
c	   write(*,*) 'd-rmax', d-rmax
c	   write(*,*) 'time=', time
	   if(rmin.lt.0.015d0.and.dt.gt.dt2)then
	    dt=dt2
	   endif
	   if(rmin.lt.0.01d0.and.dt.gt.dt3)then
	    dt=dt3
	   endif
	   if(rmin.lt.0.005d0.and.dt.gt.dt4)then
	    dt=dt4
	   endif
	   if((d-rmax).lt.0.1d0*a.and.dt.gt.dt2)then
	    dt=dt2
	   endif
	   if((d-rmax).lt.0.05d0*a.and.dt.gt.dt3)then
	    dt=dt3
	   endif
	   if((d-rmax).lt.0.02d0*a.and.dt.gt.dt4)then
	    dt=dt4
	   endif

         if (dabs(vnmax).gt.80.0D0) then
	       write(30,*) "velocity too big"
	       exit
	    end if
         if (rmin.lt.0.001D0*a) then
	       write(30,*) "neck too small, pinching"
	       write(30,*) 'rmin=', rmin
	       exit
	    end if
         if (dabs(d-rmax).lt.0.005D0*a) then
	       write(30,*) "amplitude too big, touching"
	       write(30,*) 'd-rmax=', d-rmax
	       exit
	    end if
         jj = jj + 1
c	 write(*,*) jj
	end do

      write(30,*) "stop time is", time
      write(30,*) 'dt=', dt, 'rmin', rmin
	write(30,*) 'eb=', eb/a, 'd=', d/a
        write(30,*) 'ka=', wk*a
	write(30,*) '\lambda=', ram

c output data
      open (21,file='initialc.dat')
c      open (22,file='axialv.dat')
c      open (23,file='xyval.dat')
c      open (24,file='radialv.dat')
c	open (25,file='finalxy.txt')
	open (30,file='xy_final.dat')
	 open (41,file='finalw.dat')
c	 do j=1,N+1
c	 write(23,*) z(j), r(j)
c	 enddo

	 write(8,*)  time, cvol
	 
	   jjj=jjj+1
	     if(mod(jjj,2).eq.0)then
	        do j=N+1,1,-1
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (24,*) z(j), f(j)
	        end do
	     else
	        do j=1,N+1
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (24,*) z(j), f(j)
	        end do
	     end if

	 do j=1,N+1
	   write (21,*) z0(j), r0(j)
c	   write (22,*) z(j), (u1(j)+u1p(j))/2.D0
c	   write (23,*) z(j), r(j)
c	   write (24,*) z(j), (u2(j)+u2p(j))/2.D0
	   write (30,*) z(j), r(j)	   
	   write (41,*) z(j), u1(j)
	 end do

c close 23
      end program
