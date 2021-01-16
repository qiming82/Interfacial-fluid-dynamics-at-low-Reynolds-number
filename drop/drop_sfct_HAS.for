      program vdrop_axisym
c------------------------------------------
c  viscous drop with insoluble surfactants
c subject to nonlinear straining flow Sherwood 1984
c------------------------------------------

c      implicit real*8 (a-h,o-y)
      implicit double precision (a-h,o-z)
	parameter(neq=1024)

	dimension z(neq),r(neq),z0(neq),r0(neq),s(neq),f(neq)
	dimension zold(neq),rold(neq),sold(neq),dsold(neq)
	dimension znew(neq),rnew(neq),snew(neq),v(neq)
	dimension dz(neq),ddz(neq),dr(neq),ddr(neq),ds(neq),indx(neq)
	dimension u1(neq),u2(neq),un(neq),fold(neq)
      dimension an1(neq),an2(neq),ut(neq)
      dimension anx(neq),anr(neq)
      dimension q1(neq),q2(neq),q1p(neq),q2p(neq),bb(neq),bbp(neq)
      dimension az(neq),bz(neq),cz(neq),ar(neq),br(neq),cr(neq)
      dimension d1(neq),b1(neq),c1(neq),d2(neq),b2(neq),c2(neq)
      dimension ag(neq),bg(neq),cg(neq)
c      dimension az(513),bz(513),cz(513),ar(513),br(513),cr(513)
c      dimension sumxx(neq,neq), sumxr(neq,neq), sumrx(neq,neq)
      dimension ml(neq,neq), mu(neq,neq),Am(neq,neq),hh(neq)
      dimension bbz(neq),bbr(neq),bbzp(neq),bbrp(neq),w_trp(neq)
	dimension Sxx(neq,neq),Sxr(neq,neq),Srx(neq,neq),Srr(neq,neq)
      dimension u_sxxtmp1(neq),u_sxrtmp1(neq),u_srxtmp1(neq)
	dimension u_srrtmp1(neq)
      dimension u_sxxtmp2(neq),u_sxrtmp2(neq),u_srxtmp2(neq)
	dimension u_srrtmp2(neq)
      Dimension gma(neq),sfct(neq),dgm(neq),sftmp(neq),dsfct(neq)
      dimension v0(neq)
	integer ipiv(neq)


      pi=4*datan(1.d0)
	eps = 0.0000000001D0
	      
      open (22,file='axialv.dat')
      open (23,file='xyval.dat')
      open (24,file='radialv.dat')
c      open (25,file='normalv.dat')
c      open (34,file='maragoni.dat')
c      open (26,file='sfmass.dat')
      open (27,file='tangv.dat')
      open (28,file='curvkm.dat')
      open (29,file='curv.dat')
	open (30,file='info.txt')

      open (11,file='hmin.dat')
      open (12,file='wmax.dat')
c	open (13,file='gma.dat')
c	open (14,file='tension.dat')
      open (10,file='unmax.dat')
	open (17,file='tipu.dat')

c number of segments of the interface
      N=200
	Nm=2*N+2
      h=0.5D0/N
      sp0 = h

      ram=0.80D0     ! ram = 0 is bubble
      beta=-(ram-1.0D0)/(ram+1.0D0)
	beta=beta/(4.0D0*pi)
	idef=1
	Ca = 0.d0  ! capillary # ! following Eggleton
        write(30,*) 'Ca0=', Ca
        cst2 = 0.000d0  ! nonlinear extensional flow
c       const = 2.0*pi*(ram+1.0d0)
c	sgc=1.0d0+0.2*dlog(1.0d0-xx)
     	C = Ca*4.0D0*pi 
c	coef1=-0.50d0/(2.0*pi*(ram+1.0d0))
      coef1=-0.5d0   ! after rescale the velocity

      aa1=2.5      !1.8   # 100
	bb1=1.9      !1.2   # 100
	cc1=0.3      !0.5   # 100
      istart=0
c initial settings...
       if (istart.eq.0)then
	do j=1,N+1
	   s(j)=(j-1)*h
c	   z0(j)=dsin(pi*s(j))  !!
          z0(j)=2.0D0*(aa1+bb1)*dcos(pi*s(j)) + cc1*dcos(3.0D0*pi*s(j))
         z(j)=z0(j)
c	   r0(j)=dcos(pi*s(j))
          r0(j)=(aa1-bb1)*dsin(pi*s(j)) + cc1*dsin(3.0D0*pi*s(j))
	   r(j)=r0(j)
c	   write(23,*) z0(j), r0(j)
	   sfct(j)=1.0d0
	end do    
        else
          Ca=0.13d0
          C=4.d0*pi*Ca
            call initl(z0,r0,s,sfct,N+1)
         do j=1,N+1
           z(j)=z0(j)
	   r(j)=r0(j)
	   write(23,*) z0(j), r0(j)
         enddo
c initial interface
          call arc_evl(z,r,s,N
     +     ,ar,br,cr
     +     ,az,bz,cz)
        endif

       cst2=0.0

c        cvol0=4.D0/3.D0*pi
         cvol0=0.0D0
c--------------- get the volume----------------------------------------
          do j=1,N
       	   hh(j)=z0(j+1)-z0(j)
	       cvol0=cvol0+hh(j)*(r0(j)**2+r0(j+1)**2)
	     end do
	     cvol0=dabs(cvol0)*pi
c            cvol0=4.d0*pi/3.d0
c----------------------------------------------------------------------
       if(dabs(cvol0-4.0d0*pi/3.0d0).gt.0.00001d0)then
          vscale=(cvol0/4.0d0/pi*3.0d0)**(1.D0/2.D0)
           do j=1,N+1
 	       r0(j)=r0(j)/vscale
 	    end do	   
	 endif
	do j=1,N+1
	   z(j)=z0(N+1-(j-1))
	   r(j)=r0(N+1-(j-1))
	   write(23,*) z(j), r(j)
	end do
c-----------------------------------------------------------------------
	a0 = r(1)
c   prepare to store the data


      sum11=0.0D0
	sum12=0.0D0
	sum21=0.0D0
	sum22=0.0D0
      y11 = 0.0D0
	y12 = 0.0D0
	y21 = 0.0D0
	y22 = 0.0D0
c      nplt=1

      itmax=100
	reit=10
      tol=0.000000001D0
      time=0.0D0
c	dt=0.001D0
      dt=0.2d0*(1.0D0/N)**(1.50D0)

	dt2=dt/2.0d0
	dt3=dt2/2.0d0
	dt4=dt3/2.0d0
	dt5=dt4/2.0d0
        dt6=dt5/2.0d0
        dt7=dt6/2.0d0

	!
	tf=5.D0
	nplt=tf/dt/20
	nplt2=nplt/10

       jj=1
	jjj=1
	kkk=1

c start evolution
      do while (time.lt.tf)
	   time = time + dt
c--------------------------------

	 call splc_clm
     +  (N
     +  ,s,sfct
     +  ,0.0d0
     +  ,0.0d0
     +  ,ag,bg,cg
     +  )
	  cg(N+1)=0.0d0
c---- get derivatives first
c--
          call splc_clm
     +  (N
     +  ,s,r
     +  ,0.0D0
     +  ,-pi
     +  ,ar,br,cr
     +  )
          call splc_clm
     +  (N
     +  ,s,z
     +  ,pi
     +  ,0.0D0
     +  ,az,bz,cz
     +  )

	    cz(N+1)=0.0D0
	    cr(N+1)=-pi

	    do j=1,N+1
               dsfct(j)=cg(j)
c	       dgm(j)=cg(j)
	       dr(j)=cr(j)
	       ddr(j)=2.0D0*br(j)
	       dz(j)=cz(j)
	       ddz(j)=2.0D0*bz(j)
	       
		   ds(j)=dsqrt(dz(j)**2+dr(j)**2)
             an1(j)=-dr(j)/ds(j)
	       an2(j)=dz(j)/ds(j)
	       anx(j)=-dr(j)
	       anr(j)=dz(j)

	       call df(r(j),dz(j),ddz(j),dr(j),ddr(j),ds(j),f(j))
		   q1(j)= coef1*(-f(j)*dr(j))
	       q2(j)= coef1*(f(j)*dz(j))
	    end do
           
c--------
c instantaneous arclength 
c--------	    
       arcl=0.0
	   do k=1,N
	     arcl=arcl+(s(k+1)-s(k))*(r(k)*ds(k)+r(k+1)*ds(k+1))
	   enddo
c---------------------------------------------------  
c prepare the integration of the kernels...
c--------------------------------------------------

c      step 1, single layers

	    do i=1,N+1
	        sum11=0.0D0
			sum12=0.0D0
			sum21=0.0D0
		    sum22=0.0D0
	       do k=1,N
	         call gsint
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,q1(k),q1(k+1),q2(k),q2(k+1)
     +   ,s(i),z(i),r(i)
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,y11,y12
     +   ,y21,y22
     +   )
	       sum11 = sum11 + y11
	       sum12 = sum12 + y12
	       sum21 = sum21 + y21
	       sum22 = sum22 + y22
	      end do
c  u1 = ux,  u2 = ur
	      bbz(i) = C*z(i)*(1.0+cst2*z(i)**2) + sum11 + sum12
	      bbr(i) = -C/2.0*r(i)*(1.0+3.0*cst2*z(i)**2) + sum21 + sum22
	    end do
          
	if(dabs(ram-1.0D0).lt.eps)then
         do k=1,N+1
	     u1(k)=bbz(k)
	     u2(k)=bbr(k)
	  end do
	else
c bb is the right hand side of the system
		do k=1,N+1
            bb(k)=bbz(k)
	      bb(k+N+1)=bbr(k)  
		end do
c       step 2, double layer potentials...

         do i=1,N+1
	   
	       do k=1,N
	if (idef.eq.0)then
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
	else
          call strint2
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,anx(k),anx(k+1),anr(k),anr(k+1)
     +   ,z(i),r(i),an1(i),an2(i)
     +   ,1,arcL
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,u_sxxtmp1(k),u_sxrtmp1(k)
     +   ,u_srxtmp1(k),u_srrtmp1(k)
     +   )

          call strint2
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,anx(k),anx(k+1),anr(k),anr(k+1)
     +   ,z(i),r(i),an1(i),an2(i)
     +   ,2,arcL
     +   ,az(k),bz(k),cz(k)
     +   ,ar(k),br(k),cr(k)
     +   ,u_sxxtmp2(k),u_sxrtmp2(k)
     +   ,u_srxtmp2(k),u_srrtmp2(k)
     +   )
      endif
	   if (k.eq.1) then
	      Sxx(i,k)=u_sxxtmp1(k)
	      Sxr(i,k)=u_sxrtmp1(k)
	      Srx(i,k)=u_srxtmp1(k)
	      Srr(i,k)=u_srrtmp1(k)
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
        Sxx=-beta*Sxx
	  Sxr=-beta*Sxr
	  Srx=-beta*Srx
	  Srr=-beta*Srr

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

c-----------
c  update marker points
c-----------
	    do k=1,N+1
	       un(k)=u1(k)*an1(k)+u2(k)*an2(k)
	       ut(k)=u1(k)*an2(k)-u2(k)*an1(k)
c	       znew(k)=z(k)+dt*un(k)*an1(k)
c	       rnew(k)=r(k)+dt*un(k)*an2(k)
	       znew(k)=z(k)+dt*(un(k)*an1(k)+ut(k)*an2(k))
   	       rnew(k)=r(k)+dt*(un(k)*an2(k)-ut(k)*an1(k))
c	          write (29,*) z(k), u1(k)
	     end do
c          exit
           vnmax1 = maxval(abs(un))

c+++++++++++++++++++++++++++++++++++++++++++++

c--- calculate the new arclength
c      call arc_evl(znew,rnew,snew,N
c     +  ,d2,c2,b2
c     +  ,d1,c1,b1)  ! for implicit method
c          z = znew
c          r = rnew
c          s = snew

c           call regridf(znew,rnew,z,r,s,N,time,dt)
           call regridfadp2(znew,rnew,z,r,s,N,dt,vnmax1)
c------------------
         cvol=0.0D0
c       get the volume of the jet of one period,
c          do j=1,N
c       	   hh(j)=z(j+1)-z(j)
c	       cvol=cvol+hh(j)*(r(j)**2+r(j+1)**2)
c	     end do
c	     cvol=cvol*pi
c	   write (*,*) 'ratio=', cvol/cvol0

c scale the jet
          vscale=(cvol/cvol0)**(1.D0/3.D0)
c          do j=1,N+1
cc	       z(j)=z(j)/vscale
c	       r(j)=r(j)/vscale
c	    end do
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c         rmin = minval(r(N/4:7*N/8))
         rmin = minval(r(1:N-20))
         vnmax = maxval(abs(un))
         v1max = maxval(abs(u1))
          fmax = maxval(abs(f(1:N+1)))
          if(vnmax.ne.vnmax1)then
            write(30,*) 'something is wrong...'
            exit
          endif

       dtnew=0.02d0/(5.0+fmax)/(vnmax+5.d0)

        if(dt.gt.dtnew) dt=dtnew

c	if(rmin.lt.0.025) then
c	 if(nplt.gt.100) then
c	    nplt = 100
c	  endif
c	  nplt2 = 50
c	endif

	   if (mod(jj,nplt2).eq.0) then        
	      write (10,*) time, vnmax
	      write (11,*) time, rmin
	      write (12,*) time, v1max
		write (17,*) time, u1(N+1)
            write(28,*) time, fmax
c		write(*,*) jj, dt, fmax, vnmax
	    endif

	   if (mod(jj,nplt).eq.0) then
	        do j=1,N+1
		      write (22,*) z(j), u1(j)
		      write (24,*) z(j), u2(j)
	          write(23,*) z(j), r(j)
	          write (27,*) z(j), ut(j)
	          write (29,*) z(j), f(j)
	        end do
          end if
c------
         if (vnmax.gt.50.0D0) then
	       write(30,*) "velocity too big"
	       write(30,*) 'vmax=', vnmax
	       write(30,*) 'N=', N
	       exit
	    end if
         if (rmin.lt.0.00001d0) then
	       write(30,*) "neck too small"
	       write(30,*) 'rmin=', rmin
	       write(*,*) "neck too small"
	       write(*,*) 'rmin=', rmin
	       exit
	    end if
c         if (dabs(vnmax).lt.0.0005D0) then
c	       write(30,*) "velocity too small"
c             write(30,*) 'steady drop! vmax_min=', vnmax
c             write(30,*) 'N=', N
c	       exit
c	    end if
          if(dabs(vnmax).lt.0.0002d0)then
              RLL=z(N+1)
              B  =r(1)
              D=(RLL-B)/(RLL+B)
	      write(30,*) 'Ca=', Ca, 'with Df=', D
	        do j=1,N+1
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (27,*) z(j), ut(j)
	          write (29,*) z(j), f(j)
	        end do
             if(Ca.lt.0.03d0)then
	      Ca=Ca+0.005d0
	      C = Ca/sgc*4.0*pi
             else
               write(30,*) 'steady state!'
                exit
             endif
	    endif

         jj = jj + 1
	end do

      write(30,*) "stop time is"
	write(30,*) time
	RL=z(N+1)
	B=r(1)
	D=(RL-B)/(RL+B)
	write(30,*) 'ratio D=', D
	write(30,*) 'Ca=', Ca
      write(30,*) 'dt=', dt
      write(30,*) 'N=', N
c      write(30,*) 'lambda=', ram, 'Biot=', Bi/const, 'K=', PK
c	write(30,*) 'following Eggleton'
!
c output data
       open (21,file='initialc.dat')
       open (37,file='finalxy.dat')
c       open (38,file='finalgma.dat')
	 open (41,file='finalw.dat')

	        do j=1,N+1
c	          write(13,*) z(j),sfct(j)
c	          write(14,*) z(j),gma(j)
		      write (22,*) z(j), u1(j)
	          write(23,*) z(j), r(j)
	          write (24,*) z(j), u2(j)
c	          write (25,*) z(j), un(j)
	          write (27,*) z(j), ut(j)
	          write (29,*) z(j), f(j)
	        end do


	 do j=1,N+1
	   write (21,*) z0(j), r0(j)
	   write (37,*) z(j), r(j)
c	   write (38,*) s(j), sfct(j)
	   write (41,*) z(j), u1(j)
	 end do

c close 23
      end program
