      subroutine regridf_sft_adp(znew,rnew,snew,z,r,s,f
     +  ,d2,c2,b2
     +  ,d1,c1,b1
     +  ,Nx,gma,dt)
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

        pi = 4.0*datan(1.0d0)
	RL=snew(Nx+1)

       gtmp=gma

       h=RL/Nx
       
       fmax = maxval(f(1:Nx+1))  ! max of curvature
       Ir=0  ! index of minimum
	 do 1 k=1,Nx+1
         if(dabs(fmax-f(k)).lt.0.0000001D0)then
	     Ir=k
	  endif
 1         continue

c clamped cubic spline
          call splc_clm
     +  (Nx
     +  ,snew,gtmp
     +  ,0.0d0
     +  ,0.0d0
     +  ,dg,cg,bg
     +  )
	 bg(Nx+1)=0.0d0
c--------------------------------------------------------
c IF
      IF (fmax.le.10.0d0) then ! equal-spacing regridding
       do j=1,Nx+1
	   s(j)=(j-1)*h
	end do      

         do k=1,Nx+1
      	    z(k)=seval(Nx+1,s(k),snew,znew,b1,c1,d1)
	    r(k)=seval(Nx+1,s(k),snew,rnew,b2,c2,d2)
      	    gma(k)=seval(Nx+1,s(k),snew,gtmp,bg,cg,dg)
         end do
       gma(1)=gtmp(1)
       gma(Nx+1)=gtmp(Nx+1)
        z(1)=0.0d0
        r(Nx+1)=0.0d0

c-- interpolate for bulk concentration !
c      do 91 j=1,Ny+1
c	 do i=1,Nx+1
c	  ctmp(i)=cb(i,j)
c	 enddo
c          call splc_clm
c     +  (Nx
c     +  ,snew,ctmp
c     +  ,0.0d0
c     +  ,0.0d0
c     +  ,dc,cc,bc
c     +  )
c	 bc(Nx+1)=0.0d0
c       cb(1,j)=ctmp(1)  ! end pts unchanged
c       cb(Nx+1,j)=ctmp(Nx+1)

c         do k=2,Nx
c           cb(k,j)=seval(Nx+1,s(k),snew,ctmp,bc,cc,dc)
c	   enddo
c 91     continue
c----------------------------------------------------------
c ELSE --------------
      ELSEIF (fmax.gt.10.0d0.and.Ir.eq.Nx+1)then ! adaptive
         do 21 j=1,Ir-1
            fk=dabs(f(j))
             if(fk.gt.1.0d0)then
               sd(j) = 1.d0/fk
c             elseif(fk.ge.100.0)then
c                sd(j)=1.2d0/fk
              else
               sd(j) = 0.25d0*f(j)**2
              endif
 21         continue
          suml = 0.0d0
          do 3 j=1,Nx
               suml = suml + sd(j)
 3           continue
           hs = RL/suml
           stmp(Nx+1)=RL
           do 4 j=Nx,2,-1
               stmp(j)=stmp(j+1)-hs*sd(j)
 4            continue
           stmp(1)=0.0d0
c-- evaluate output z and r
         s = stmp
         do k=1,Nx+1
      	  z(k)=seval(Nx+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(Nx+1,s(k),snew,rnew,b2,c2,d2)
      	  gma(k)=seval(Nx+1,s(k),snew,gtmp,bg,cg,dg)
c          write(30,*) 'gmak', gma(k)
         end do
       gma(1)=gtmp(1)
       gma(Nx+1)=gtmp(Nx+1)
       z(1)=0.0d0
       r(Nx+1)=0.0d0
         dtnew=0.003d0/(fmax+10.0)
        if(dt.gt.dtnew) dt=dtnew
c-------------------------------------------------------
      ELSE
        do 31 j=1,Nx
          fk=dabs(f(j))
           if(fk.gt.1.d0)then
            sd(j) = 2.5d0/fk
c           elseif(fk.ge.100.d0)then
c             sd(j)=1.2d0/fk
           else
            sd(j) = 0.5d0*f(j)**2
           endif
 31      continue
c-- change !
        do 32 j=1,Nx
            suml = suml + sd(j)
 32      continue
         hs = RL/suml
         stmp(1)=0.d0
         do 34 j=2,Nx
              stmp(j)=stmp(j-1)+hs*sd(j-1)
 34          continue
          stmp(Nx+1)=RL
        s = stmp   
         do k=1,Nx+1
      	  z(k)=seval(Nx+1,s(k),snew,znew,b1,c1,d1)
	  r(k)=seval(Nx+1,s(k),snew,rnew,b2,c2,d2)
      	  gma(k)=seval(Nx+1,s(k),snew,gtmp,bg,cg,dg)
         end do
       gma(1)=gtmp(1)
       gma(Nx+1)=gtmp(Nx+1)
       z(1)=0.0d0
       r(Nx+1)=0.0d0

       dtnew=0.0015d0/(fmax+10.0)
      if(dt.gt.dtnew) dt=dtnew
c-- gma and Cb unchanged, so nothing needs to do about them
      Endif
cc-- DONE !!

      Return
      End
c-------------
      subroutine initld (z,r,sfct,s,N1)   
c
      implicit double precision (a-h,o-z)
	parameter(neq=1024)
       dimension z(neq),r(neq),sfct(neq)
       dimension s(neq),Cb(neq,neq),y(neq)
c
       open(1,file='xy0.dat',status='old')
       open(2,file='sfct0.dat',status='old')
       do k=1,N1
          read(1,*) z(k), r(k)
          read(2,*) s(k), sfct(k)
        enddo
      return
       end

