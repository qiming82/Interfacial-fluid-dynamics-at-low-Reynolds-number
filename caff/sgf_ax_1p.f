      subroutine sgf_ax_1p
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,RL
     +   ,Nsum,Np
     +   ,GXX,GXY
     +   ,GYX,GYY
     +   ,QXXX,QXXY,QXYX,QXYY
     +   ,QYXX,QYXY,QYYX,QYYY
     +   ,PXX,PXY,PYX,PYY
     +   ,Iaxis
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c======================================================
c  Axisymmetric periodic Green's function of Stokes flow
c
c  Symbols:
c  -------
c
c  RL: 		period along the x axis
c  Nsum, Np: 	summation range limits
c
c  y stands for radial distance "sigma"
c
c-------------
c
c  Let (bx, by) be the strength of the point-force ring;
c  the induced velocity is:
c
c  ux(x0) = SXX(x,x0) * bx + SXY(x,x0) * by
c  uy(x0) = SYX(x,x0) * bx + SYY(x,x0) * by
c
c  The kernel of the axisymmetric double-layer potential is:
c
c   Idlpx(x0) = ux * ( Qxxx * vnx + Qxxy * vny)
c             + uy * ( Qxyx * vnx + Qxyy * vny)
c
c   Idlpy(x0) = ux * ( Qyxx * vnx + Qyxy * vny)
c             + uy * ( Qyyx * vnx + Qyyy * vny)
c
c   This is the flow due to a ring distribution of stresslets
c
c   The arguments of Q_abc are (x,x0)
c
c   Summation will be done over x
c
c------------
c
c  Iopt = 1 only the Green's function
c  Iopt = 2          Green's function and the stress tensor
c
c======================================================

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c--------
c prepare
c--------

      subtr = pi4*Y
      RLH   = 0.5D0*RL

c---------------------------
c prepare for fast summation
c using Aitken extrapolation
c---------------------------

      N0 = Nsum
      N1 = N0*Np
      N2 = N1*Np

c-------------
c primary ring
c-------------

      call sgf_ax_fs
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,GXX,GXY
     +   ,GYX,GYY
     +   ,QXXX,QXXY,QXYX,QXYY
     +   ,QYXX,QYXY,QYYX,QYYY
     +   ,PXX,PXY,PYX,PYY
     +   ,Iaxis
     +   )

      GXX = GXX - subtr/Dabs(X+RLH)

c------------------------
c sum over periodic rings
c------------------------

      Xsave = X

      Do 1 Ir=1,N2

      X = Xsave+Ir*RL     ! forward ring

      call sgf_ax_fs
     +
     +   (Iopt
     +   ,X,Y,X0,Y0
     +   ,SXX,SXY
     +   ,SYX,SYY
     +   ,rXXX,rXXY,rXYX,rXYY
     +   ,rYXX,rYXY,rYYX,rYYY
     +   ,eXX,eXY,eYX,eYY
     +   ,Iaxis
     +   )

      SXX = SXX - subtr/Dabs(X+RLH)

      If   (Ir.eq.N0         ! hold for extrapolation
     +  .or.Ir.eq.N1
     +  .or.Ir.eq.N2) then   

       fxx = sxx
       fxy = sxy
       fyx = syx
       fyy = syy

        If(Iopt.eq.2) then  ! stress and pressure
         fxxx = rxxx
         fxxy = rxxy
         fxyx = rxyx
         fxyy = rxyy
         fyxx = ryxx
         fyxy = ryxy
         fyyx = ryyx
         fyyy = ryyy
         uxx  = exx
         uxy  = exy
         uyx  = eyx
         uyy  = eyy
        End If

      End If

      GXX = GXX + SXX
      GXY = GXY + SXY
      GYX = GYX + SYX
      GYY = GYY + SYY

      If(Iopt.eq.2) then

       QXXX = QXXX + RXXX
       QXXY = QXXY + RXXY
       QXYX = QXYX + RXYX
       QXYY = QXYY + RXYY
       QYXX = QYXX + RYXX
       QYXY = QYXY + RYXY
       QYYX = QYYX + RYYX
       QYYY = QYYY + RYYY

       PXX  = PXX+EXX
       PXY  = PXY+EXY
       PYX  = PYX+EYX
       PYY  = PYY+EYY

      End If

      X = Xsave-Ir*RL      ! backward ring

      call sgf_ax_fs
     +
     +   (Iopt
     +   ,X,Y,X0,Y0
     +   ,SXX,SXY
     +   ,SYX,SYY
     +   ,RXXX,RXXY,RXYX,RXYY
     +   ,RYXX,RYXY,RYYX,RYYY
     +   ,EXX,EXY,EYX,EYY
     +   ,IAXIS
     +   )

      SXX = SXX - subtr/Dabs(X+RLH)

      If   (Ir.eq.N0           ! hold for extrapolation
     +  .or.Ir.eq.N1
     +  .or.Ir.eq.N2) then   

       fxx  = fxx  + sxx
       fxy  = fxy  + sxy
       fyx  = fyx  + syx
       fyy  = fyy  + syy

       If(Iopt.eq.2) then
        fxxx = fxxx + rxxx
        fxxy = fxxy + rxxy
        fxyx = fxyx + rxyx
        fxyy = fxyy + rxyy
        fyxx = fyxx + ryxx
        fyxy = fyxy + ryxy
        fyyx = fyyx + ryyx
        fyyy = fyyy + ryyy
        uxx  = uxx  + exx
        uxy  = uxy  + exy
        uyx  = uyx  + eyx
        uyy  = uyy  + eyy
       End If

      End If

      GXX = GXX + SXX
      GXY = GXY + SXY
      GYX = GYX + SYX
      GYY = GYY + SYY

      If(Iopt.eq.2) then
       QXXX = QXXX + RXXX
       QXXY = QXXY + RXXY
       QXYX = QXYX + RXYX
       QXYY = QXYY + RXYY
       QYXX = QYXX + RYXX
       QYXY = QYXY + RYXY
       QYYX = QYYX + RYYX
       QYYY = QYYY + RYYY
       PXX  = PXX+EXX
       PXY  = PXY+EXY
       PYX  = PYX+EYX
       PYY  = PYY+EYY
      End If

      If(Ir.eq.N0) then
        Gxx0  = Gxx - 0.5D0 * fxx
        Gxy0  = Gxy - 0.5D0 * fxy
        Gyx0  = Gyx - 0.5D0 * fyx
        Gyy0  = Gyy - 0.5D0 * fyy
      Else If(Ir.eq.N1) then
        Gxx1  = Gxx - 0.5D0 * fxx
        Gxy1  = Gxy - 0.5D0 * fxy
        Gyx1  = Gyx - 0.5D0 * fyx
        Gyy1  = Gyy - 0.5D0 * fyy
      Else If(Ir.eq.N2) then
        Gxx2  = Gxx - 0.5D0 * fxx
        Gxy2  = Gxy - 0.5D0 * fxy
        Gyx2  = Gyx - 0.5D0 * fyx
        Gyy2  = Gyy - 0.5D0 * fyy
      End If

c---
      If(Iopt.eq.2) then
       If(Ir.eq.N0) then
        Qxxx0 = Qxxx - 0.5D0 * fxxx
        Qxxy0 = Qxxy - 0.5D0 * fxxy
        Qxyx0 = Qxyx - 0.5D0 * fxyx
        Qxyy0 = Qxyy - 0.5D0 * fxyy
        Qyxx0 = Qyxx - 0.5D0 * fyxx
        Qyxy0 = Qyxy - 0.5D0 * fyxy
        Qyyx0 = Qyyx - 0.5D0 * fyyx
        Qyyy0 = Qyyy - 0.5D0 * fyyy
         Pxx0 = Pxx  - 0.5D0 * uxx
         Pxy0 = Pxy  - 0.5D0 * uxy
         Pyx0 = Pyx  - 0.5D0 * uyx
         Pyy0 = Pyy  - 0.5D0 * uyy
       Else If(Ir.eq.N1) then
        Qxxx1 = Qxxx - 0.5D0 * fxxx
        Qxxy1 = Qxxy - 0.5D0 * fxxy
        Qxyx1 = Qxyx - 0.5D0 * fxyx
        Qxyy1 = Qxyy - 0.5D0 * fxyy
        Qyxx1 = Qyxx - 0.5D0 * fyxx
        Qyxy1 = Qyxy - 0.5D0 * fyxy
        Qyyx1 = Qyyx - 0.5D0 * fyyx
        Qyyy1 = Qyyy - 0.5D0 * fyyy
         Pxx1 = Pxx  - 0.5D0 * uxx
         Pxy1 = Pxy  - 0.5D0 * uxy
         Pyx1 = Pyx  - 0.5D0 * uyx
         Pyy1 = Pyy  - 0.5D0 * uyy
       Else If(Ir.eq.N2) then
        Qxxx2 = Qxxx - 0.5D0 *fxxx
        Qxxy2 = Qxxy - 0.5D0 *fxxy
        Qxyx2 = Qxyx - 0.5D0 *fxyx
        Qxyy2 = Qxyy - 0.5D0 *fxyy
        Qyxx2 = Qyxx - 0.5D0 *fyxx
        Qyxy2 = Qyxy - 0.5D0 *fyxy
        Qyyx2 = Qyyx - 0.5D0 *fyyx
        Qyyy2 = Qyyy - 0.5D0 *fyyy
         Pxx2 = Pxx  - 0.5D0 *uxx
         Pxy2 = Pxy  - 0.5D0 *uxy
         Pyx2 = Pyx  - 0.5D0 *uyx
         Pyy2 = Pyy  - 0.5D0 *uyy
       End If
      End If
c---

   1  Continue

      X = Xsave  ! restore

c---
c Aitken extrapolation for the velocity
c---

      GXX = (Gxx0*Gxx2-Gxx1**2)/(Gxx2-2.0D0*Gxx1+Gxx0)
      GXY = (Gxy0*Gxy2-Gxy1**2)/(Gxy2-2.0D0*Gxy1+Gxy0)
      GYX = (Gyx0*Gyx2-Gyx1**2)/(Gyx2-2.0D0*Gyx1+Gyx0)
      GYY = (Gyy0*Gyy2-Gyy1**2)/(Gyy2-2.0D0*Gyy1+Gyy0)

c---
c Aitken extrapolation for the stress and pressure
c---

      If(Iopt.eq.2) then

      Qxxx = (Qxxx0*Qxxx2-Qxxx1**2)/(Qxxx2-2.0D0 *Qxxx1+Qxxx0)
      Qxxy = (Qxxy0*Qxxy2-Qxxy1**2)/(Qxxy2-2.0D0 *Qxxy1+Qxxy0)
      Qxyx = (Qxyx0*Qxyx2-Qxyx1**2)/(Qxyx2-2.0D0 *Qxyx1+Qxyx0)
      Qxyy = (Qxyy0*Qxyy2-Qxyy1**2)/(Qxyy2-2.0D0 *Qxyy1+Qxyy0)

      Qyxx = (Qyxx0*Qyxx2-Qyxx1**2)/(Qyxx2-2.0D0 *Qyxx1+Qyxx0)
      Qyxy = (Qyxy0*Qyxy2-Qyxy1**2)/(Qyxy2-2.0D0 *Qyxy1+Qyxy0)
      Qyyx = (Qyyx0*Qyyx2-Qyyx1**2)/(Qyyx2-2.0D0 *Qyyx1+Qyyx0)
      Qyyy = (Qyyy0*Qyyy2-Qyyy1**2)/(Qyyy2-2.0D0 *Qyyy1+Qyyy0)

      Pxx = (Pxx0*Pxx2-Pxx1**2)/(Pxx2-2.0D0 *Pxx1+Pxx0)
      Pxy = (Pxy0*Pxy2-Pxy1**2)/(Pxy2-2.0D0 *Pxy1+Pxy0)
      Pyx = (Pyx0*Pyx2-Pyx1**2)/(Pyx2-2.0D0 *Pyx1+Pyx0)
      Pyy = (Pyy0*Pyy2-Pyy1**2)/(Pyy2-2.0D0 *Pyy1+Pyy0)

      End If

c---
c optional printing session
c---
c
c     write (6,101) Gxx0,Gxx1,Gxx2,Gxx
c     write (6,101) Gxy0,Gxy1,Gxy2,Gxy
c     write (6,101) Gyx0,Gyx1,Gyx2,Gyx
c     write (6,101) Gyy0,Gyy1,Gyy2,Gyy
c
c     write (6,101) Qxxx0,Qxxx1,Qxxx2,Qxxx
c     write (6,101) Qxxy0,Qxxy1,Qxxy2,Qxxy
c     write (6,101) Qxyx0,Qxyx1,Qxyx2,Qxyx
c     write (6,101) Qxyy0,Qxyy1,Qxyy2,Qxyy
c
c     write (6,101) Qyxx0,Qyxx1,Qyxx2,Qyxx
c     write (6,101) Qyxy0,Qyxy1,Qyxy2,Qyxy
c     write (6,101) Qyyx0,Qyyx1,Qyyx2,Qyyx
c     write (6,101) Qyyy0,Qyyy1,Qyyy2,Qyyy
c
c     write (6,101) Pxx0,Pxx1,Pxx2,Pxx
c     write (6,101) Pxy0,Pxy1,Pxy2,Pxy
c     write (6,101) Pyx0,Pyx1,Pyx2,Pyx
c     write (6,101) Pyy0,Pyy1,Pyy2,Pyy
c---------

c-----
c Done
c-----

 100  format (1x,i3,10(1x,f10.5))
 101  format (10(1x,f12.8))

      Return
      End
