      subroutine int_pot_sc
     +
     +   (z1,r1
     +   ,s1,s2
     +   ,z0,r0
     +   ,n,d
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,y1
     +   )
      Implicit Double Precision (a-h,o-z)

	dimension x1(10),w1(10),amg(10),dsg(10)  
      dimension sg(10),rg(10),drg(10),zg(10),dzg(10)
	   
	eps = 0.000000001D0
      
	pi  = 3.1415 92653 58979 32384 D0
      RL=1.D0
	Nsum=5
	Np=2
	d2=2.0D0*d
      
	call Gauss_Legendre (8,x1,w1)

      if (n.eq.1) then
	  do i=1,8
	     amg(i)=(1.0D0-x1(i))/2.0D0
	   end do
	else if (n.eq.2) then
	   do i=1,8
	      amg(i)=(1.0D0+x1(i))/2.0D0
	    end do
	end if

	a = (s2-s1)/2.0D0
	b = (s2+s1)/2.0D0

	do j=1,8
	   sg(j) = a*x1(j) + b
c	   anzg(j) = (anx2-anx1)/2.0D0*x1(j) + (anx2+anx1)/2.0D0
c	   anrg(j) = (any2-any1)/2.0D0*x1(j) + (any2+any1)/2.0D0
	   xd = sg(j) - s1
	   rg(j) = Ar*xd**3 + Br*xd**2 + Cr*xd + r1
	   drg(j) = 3.D0*Ar*xd**2 + 2.D0*Br*xd + Cr
	   zg(j) = Az*xd**3 + Bz*xd**2 + Cz*xd + z1
	   dzg(j) = 3.D0*Az*xd**2 + 2.D0*Bz*xd + Cz
	   dsg(j) = dsqrt(drg(j)**2+dzg(j)**2)
c	   akg(j) = 4.0D0*r0*rg(j)/((zg(j)-z0)**2+(rg(j)+r0)**2)
c	   akgi(j) = 4.0D0*r0*rg(j)/((zg(j)+z0)**2+(rg(j)+r0)**2)
	end do
 
      Iopt = 1  
c ipot ne 1    both green's function and its derivative
      tmp = 0.0D0

	do j=1,8
          call lgf_ax_gp_ct
     +
     +   (Iopt,1
     +   ,zg(j),rg(j)
     +   ,z0,r0
     +   ,RL,d
     +   ,Nsum,Np
     +   ,gc1
     +   )
          call lgf_ax_gp_ct
     +
     +   (Iopt,0  ! 0 for one subtraction
     +   ,zg(j),rg(j)
     +   ,-z0,r0
     +   ,RL,d
     +   ,Nsum,Np
     +   ,gc2
     +   )
c        call Gc(rg(j),zg(j),r0,z0,gc1,d,1)
c	  call Gc(rg(j),zg(j),r0,-z0,gc2,d,0)

          call lgf_ax_gpc_ct
     +
     +   (zg(j),rg(j),z0,r0
c     +   (z0,r0,zg(j),rg(j)
     +   ,d,RL
     +   ,fc1
     +   )
          call lgf_ax_gpc_ct
     +
     +   (zg(j),rg(j),-z0,r0
c     +   (-z0,r0,zg(j),rg(j)
     +   ,d,RL
     +   ,fc2
     +   )
c        fc1=0.0d0
c	  fc2=0.0d0 
		tmp=tmp+a*amg(j)*(gc1+gc2+fc1+fc2)*dsg(j)*rg(j)*w1(j)
c		tmp=tmp+a*amg(j)*(gc1+gc2)*dsg(j)*rg(j)*w1(j)
	 end do

	  y1 = tmp

c
c ---- done
	     
	return
	end
c-----------------------------
	subroutine Gc(r,z,r0,z0,y,d,iop)
      Implicit Double Precision (a-h,o-z)
      dimension c(200)

      Iopt=1
	RL=1.D0
	m=50
      
c	do j=1,2*m+1
c	   c(j)=1.D0
c	end do
c	c(1)=0.5D0
c	c(2*m+1)=0.5D0
      sum=0.0D0
	
	do i=-m,m
	   if(i.ne.0)then
	       call lgf_ax_fs(Iopt,z,r,z0+i*RL,r0,G,Gz,Gr)
	       dx = z - z0 - i*RL
	       sss = 2.0D0*d - r - r0
	       G = G - 0.5d0/dsqrt(dx*dx+sss*sss) 
c	       call lgf_ax_fs(Iopt,z,d2-r,z0+i*RL,r0,Gi,Gzi,Gri)
	      sum =sum + G
          end if
	end do
	
	call lgf_ax_fs(Iopt,z,r,z0+RL,r0,G1,Gz,Gr)
      dx = z -z0-RL
	sss= 2.0d0*d-r-r0
	G1 = G1 - 0.5d0/dsqrt(dx*dx+sss*sss) 
      
	if(iop.eq.0)then
	   y = sum-G1
	else
	   y = sum
	end if
	return 
	end
