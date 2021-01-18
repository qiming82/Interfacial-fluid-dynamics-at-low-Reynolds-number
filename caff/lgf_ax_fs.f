      subroutine lgf_ax_fs 
     +
     +   (Iopt
     +   ,x,s
     +   ,x0,s0
     +   ,G
     +   ,Gx,Gs
     +   )

c=========================================
c FDLIB, CFDLAB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c=========================================

c-----------------------------------------
c Free-space axisymmetric Green's function.
c of Laplace's equation
c
c  Iopt =  1 compute only the Green's function
c       ne 1 compute Green's function and gradient
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0*pi

c-------- 
c prepare
c-------- 

      Dx  = x-x0
      Dxs = Dx*Dx
      ss0s = (s+s0)**2

      rks = 4.0D0*s*s0/(Dxs+ss0s)

      call ell_int2 (rks,F,E)

c-----------------
c Green's function
c-----------------

      RJ10 = F
      den = Dsqrt(Dxs+ss0s)

      RI10 = 4.0D0*RJ10/den

      G = RI10/pi4

      If(Iopt.eq.1) Go to 99

c---------------------
c compute: I30 and I31
c---------------------

      rksc = 1.0D0-rks
      RJ30 = E/rksc
      RJ31 = (-2.0D0*F+(2.0D0-rks)*E/rksc)/rks
      cf   = 4.0D0/den**3
      RI30 = cf*RJ30
      RI31 = cf*RJ31

c---------
c gradient
c---------

      Gx = - dx * RI30
      Gs = - s*RI30+s0*RI31

      Gx = Gx/pi4
      Gs = Gs/pi4

c-----
c Done
c-----

  99  Continue

      Return
      End
