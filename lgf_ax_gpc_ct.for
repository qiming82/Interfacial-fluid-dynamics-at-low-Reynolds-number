      subroutine lgf_ax_gpc_ct 
     +
     +   (x,s       ! evaluation point
     +   ,x0,s0     ! singular point
     +   ,sc        ! tube radius
     +   ,RL        ! period
     +   ,gpc
     +   )


c-----------------------------------------------------
c  Flow due to a periodic array of point-force rings
c  inside a cylinder of radius sc
c
c  Singularity is located at (x0, s0)
c
c  Green's function is evaluated at (x, s)
c
c  SYMBOLS:
c  --------
c
c  RL: period along the x axis (cylinder axis)
c  sc: radius of cylinder
c
c                                   
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (tol=0.0000001)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c----------------------------------
c If the poionts are on the axis,
c move them slightly off the axis
c to facilitate the computation
c----------------------------------

      If( s.lt.tol) s  = tol
      If(s0.lt.tol) s0 = tol

c--------
c prepare
c--------

      dx = x-x0
      wn = pi2/RL       ! wave number

      Iswitch = 2       ! for the Bessel functions

c-------------------------------
c Compute the complementary part
c-------------------------------

      ssss = s+s0-2.0D0*sc

c---
c set truncation limit
c---

c     M = - int(25.0D0/ssss)
c     M = - int(15.0D0/ssss)
      M = 12
c      M = 10

      gc = 0.0D0

      Do 1 i=0,M

        t = i*wn

        If(i.eq.0) t = 0.0001D0
        oc = sc*t
        o0 = s0*t
        o  = s *t
        osn  = - ssss *t

        call bess_I01K01 (oc,Iswitch,BI0c,BI1c,BK0c,BK1c)
        call bess_I01K01 (o0,Iswitch,BI00,BI10,BK00,BK10)
        call bess_I01K01 (o,Iswitch,BI0,BI1,BK0,BK1)

        f = -1.0d0/pi*BI0*BI00/BI0c*BK0c 

        call bess_I01K01 (osn,Iswitch,BI0sn,BI1sn,BK0sn,BK1sn)

        f = f + 1.0d0/pi*BK0sn

        dxi = dx*t
        dc  = Dcos(dxi)
        ds  = Dsin(dxi)

        fc = 1.0D0
        If(i.eq.0) fc = 0.50D0

        gc = gc + f * dc * fc

c        write (6,100) i,fxx,fxs,fsx,fss
c        if(i.eq.0) write (6,*) fxx

  1   Continue

      ccff  = wn
c      hh  = Dsqrt(dx*dx+ssss*ssss)
      gc = ccff*gc


      gpc = gc

c-----
c Done
c-----

  100 Format (1x,i3,8(1x,f15.10))
  101 Format (9(1x,f15.10))

      Return
      End
