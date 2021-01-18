      subroutine sgf_ax_1p_ct 
     +
     +   (x,s       ! evaluation point
     +   ,x0,s0     ! singular point
     +   ,sc        ! tube radius
     +   ,RL        ! period
     +   ,Nsum,Np
     +   ,gxx,gxs
     +   ,gsx,gss
     +   )

c-----------------------------------------
c FDLIB and BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

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
c  N0, N1, N2: 
c
c    Number of terms retained in the infinite sum
c    of sum in real space
c
c    Computations are performed for 
c
c    N0=Nsum, N1=N0*Np, N2=N1*Np
c
c    and the result is extapolated using Aitken's method.
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
c     M = 20
      M = 10

c     write (6,*) "sgf_ax_ct: comp sum truncated at M=",M

      rmcxx = 0.0D0
      rmcxs = 0.0D0
      rmcsx = 0.0D0
      rmcss = 0.0D0

      Do 1 i=0,M

        t = i*wn

        If(i.eq.0) t = 0.0001D0
        oc = sc*t
        o0 = s0*t
        o  = s *t
        osn  = - ssss *t

        call bess_I01K01 (oc,Iswitch,BI0c,BI1c,BK0c,BK1c)
        call bess_I01K01 (o0,Iswitch,BI00,BI10,BK00,BK10)

        bxx =  (-2.0D0*BK0c+oc*BK1c)*BI00 - BK0c*o0*BI10 
        bxs =              -oc*BK0c *BI00 + BK1c*o0*BI10
        bsx =              -oc*BK1c *BI10 + BK0c*o0*BI00
        bss =  ( 2.0D0*BK1c+oc*BK0c)*BI10 - BK1c*o0*BI00 

        temp = oc*BI1c+2.0D0*BI0c
        Det  = t*( oc*BI0c*BI0c-BI1c*temp )
        fc   = 4.0D0/det

        axx = fc*  ( oc*BI0c*bxx - bxs*temp )
        asx = fc*  ( oc*BI0c*bsx - bss*temp )
        axs = fc*t*(    BI0c*bxs - bxx*BI1c ) 
        ass = fc*t*(    BI0c*bss - bsx*BI1c )

        call bess_I01K01 (o,Iswitch,BI0,BI1,BK0,BK1)

        fxx = t*BI0*axx + (o*BI1+2.0*BI0)*axs
        fxs = t*BI0*asx + (o*BI1+2.0*BI0)*ass
        fsx = t*BI1*axx +  o*BI0         *axs
        fss = t*BI1*asx +  o*BI0         *ass

        call bess_I01K01 (osn,Iswitch,BI0sn,BI1sn,BK0sn,BK1sn)

        fxx = fxx + 8.0D0*BK0sn

        dxi = dx*t
        dc  = Dcos(dxi)
        ds  = Dsin(dxi)

        fc = 1.0D0
        If(i.eq.0) fc = 0.50D0

        rmcxx = rmcxx + fxx * dc * fc
        rmcxs = rmcxs + fxs * ds * fc
        rmcsx = rmcsx + fsx * ds * fc
        rmcss = rmcss - fss * dc * fc

c        write (6,100) i,fxx,fxs,fsx,fss
c        if(i.eq.0) write (6,*) fxx

  1   Continue

      ccff  = s0*wn

      rmcxx = ccff*rmcxx
      rmcxs = ccff*rmcxs
      rmcsx = ccff*rmcsx
      rmcss = ccff*rmcss

c---------------------------------
c Sum in real space over the rings
c---------------------------------

      dx0 = dx

      ss0   = s*s0
      ss04  = 4.0D0*ss0
      rss0  = Dsqrt(ss0)
      s2    = s*s
      s02   = s0*s0
      sc2   = sc*sc
      rs0   = Dsqrt(s0)
      rs    = Dsqrt(s)
      sss02 = (s+s0)**2
      ssd02 = (s-s0)**2
      ssss2 = ssss**2

c---
c prepare for fast summation
c with Aitken extrapolation
c---

      N0 = Nsum
      N1 = N0*Np
      N2 = N1*Np

c---
c primary ring
c---

      dx  = dx0
      dx2 = dx*dx
      r2  = dx2+ssd02
      rk2 = ss04 / (dx2+sss02)
      rk  = Dsqrt(rk2)

      call ell_int (rk2,f,e)

      hh  = Dsqrt(dx2+ssss2)

      rmrxx = 2.0D0*rk*(f+dx2*e/r2)*rs0/rs  - pi4*s0/hh
      rmrxs = - rk*dx*(f-(s2-s02+dx2)*e/r2)/rss0
      rmrsx = + rk*dx*(f+(s2-s02-dx2)*e/r2)*rs0/(s*rs)
      rmrss = + rk*((s02+s2+2.0*dx2)*f
     +        -( 2.0D0*dx2*dx2 +  3.0D0*dx2*(s02+s2)
     +        +  (s2-s02)*(s2-s02) )*e/r2 )/(rs0*s*rs)

c---
c sum over periodic rings
c---

      Do 2 Ir=1,N2

        dx = dx0+Ir*RL     ! ring to the left

        dx2 = dx*dx
        r2  = dx2+ssd02
        rk2 = ss04 / (dx2+sss02)
        rk  = Dsqrt(rk2)

        call ell_int (rk2,f,e)

        hh  = Dsqrt(dx2+ssss2)
        xx  = 2.0*rk*(f+dx2*e/r2)*rs0/rs  - pi4*s0/hh
        xs  = - rk*dx*(f-(s2-s02+dx2)*e/r2)/rss0
        sx  = + rk*dx*(f+(s2-s02-dx2)*e/r2)*rs0/(s*rs)
        ss  = + rk*((s02+s2+2.0*dx2)*f
     +           -( 2.0*dx2*dx2 + 3.0*dx2*(s02+s2)
     +           +  (s2-s02)*(s2-s02) )*e/r2 )/(rs0*s*rs)

        If   (Ir.eq.N0         ! hold for extrapolation
     +    .or.Ir.eq.N1
     +    .or.Ir.eq.N2) then
         fxx  = xx
         fxs  = xs
         fsx  = sx
         fss  = ss
        End If

        rmrxx = rmrxx + xx
        rmrxs = rmrxs + xs
        rmrsx = rmrsx + sx
        rmrss = rmrss + ss

        dx  = dx0-Ir*RL    ! ring to the right

        dx2 = dx*dx
        r2  = dx2+ssd02
        rk2 = ss04 / (dx2+sss02)
        rk  = Dsqrt(rk2)

        call ell_int (rk2,f,e)

        hh = sqrt(dx2+ssss2)
        xx = 2.0D0 * rk * (f+dx2*e/r2)*rs0/rs  - pi4*s0/hh
        xs = - rk*dx*(f-(s2-s02+dx2)*e/r2)/rss0
        sx = + rk*dx*(f+(s2-s02-dx2)*e/r2)*rs0/(s*rs)
        ss = + rk*((s02+s2+2.0*dx2)*f
     +           -( 2.0D0 * dx2*dx2 +  3.0D0 * dx2*(s02+s2)
     +           +  (s2-s02)*(s2-s02) )*e/r2 )/(rs0*s*rs)

        If   (Ir.eq.N0           ! hold for extrapolation
     +    .or.Ir.eq.N1
     +    .or.Ir.eq.N2) then
         fxx = fxx + xx
         fxs = fxs + xs
         fsx = fsx + sx
         fss = fss + ss
        End If

        rmrxx = rmrxx + xx
        rmrxs = rmrxs + xs
        rmrsx = rmrsx + sx
        rmrss = rmrss + ss

      If(Ir.eq.N0) then
        Gxx0 = rmrxx - 0.5D0*fxx
        Gxs0 = rmrxs - 0.5D0*fxs
        Gsx0 = rmrsx - 0.5D0*fsx
        Gss0 = rmrss - 0.5D0*fss
      Else If(Ir.eq.N1) then
        Gxx1 = rmrxx - 0.5D0*fxx
        Gxs1 = rmrxs - 0.5D0*fxs
        Gsx1 = rmrsx - 0.5D0*fsx
        Gss1 = rmrss - 0.5D0*fss
      Else If(Ir.eq.N2) then
        Gxx2 = rmrxx - 0.5D0*fxx
        Gxs2 = rmrxs - 0.5D0*fxs
        Gsx2 = rmrsx - 0.5D0*fsx
        Gss2 = rmrss - 0.5D0*fss
      End If

c        write (6,100) i,xx,xs,sx,ss

   2  Continue

c---------------------
c Aitken extrapolation
c---------------------

      rmrxx = (Gxx0*Gxx2-Gxx1**2)/(Gxx2-2.0D0*Gxx1+Gxx0)
      rmrxs = (Gxs0*Gxs2-Gxs1**2)/(Gxs2-2.0D0*Gxs1+Gxs0)
      rmrsx = (Gsx0*Gsx2-Gsx1**2)/(Gsx2-2.0D0*Gsx1+Gsx0)
      rmrss = (Gss0*Gss2-Gss1**2)/(Gss2-2.0D0*Gss1+Gss0)

c     write (6,101) Gxx0,Gxx1,Gxx2,rmrxx,rmcxx
c     write (6,101) Gxs0,Gxs1,Gxs2,rmrxs,rmcxs
c     write (6,101) Gsx0,Gsx1,Gsx2,rmrsx,rmcsx
c     write (6,101) Gss0,Gss1,Gss2,rmrss,rmcss

c     rmrxx = Gxx2
c     rmrxs = Gxs2
c     rmrsx = Gsx2
c     rmrss = Gss2

c-----------------------------------------------
c add the primary and complementary componenents
c-----------------------------------------------

      gxx = rmrxx + rmcxx
      gxs = rmrxs + rmcxs
      gsx = rmrsx + rmcsx
      gss = rmrss + rmcss

c        write (6,*) " MC"
c        write (6,*) " --"
c        write (6,101) rmcxx,rmcxs
c        write (6,101) rmcsx,rmcss
c        write (6,*) " MR"
c        write (6,*) " --"
c        write (6,101) rmrxx,rmrxs
c        write (6,101) rmrsx,rmrss
c        write (6,*) " MT"
c        write (6,*) " --"
c        write (6,101) gxx,gxs
c        write (6,101) gsx,gss

c-----
c Done
c-----

  100 Format (1x,i3,8(1x,f15.10))
  101 Format (9(1x,f15.10))

      Return
      End
