      subroutine sgf_ax_str
     +
     + (Iopt
     + ,X,Y
     + ,X0,Y0
     + ,QXXX,QXXY,QXYX,QXYY
     + ,QYXX,QYXY,QYYX,QYYY
     + ,Iaxis
     + )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c================================================
c  Axisymmetric Green's function of Stokes flow
c  in free space
c
c  Let b be the strength of the point-force ring located at x;
c  then the induced velocity field is:
c
c  ux(x0) = Sxx(x,x0) * bx + Sxs(x,x0) * bs
c  us(x0) = Ssx(x,x0) * bx + Sss(x,x0) * bs
c
c  The kernel of the axisymmetric double-layer potential is:
c
c   Idlpx(x0) = ux * ( Qxxx * vnx + Qxxs * vns)
c             + us * ( Qxsx * vnx + Qxss * vns)
c
c   Idlps(x0) = ux * ( Qsxx * vnx + Qsxs * vns)
c             + us * ( Qssx * vnx + Qsss * vns)
c
c  arguments of Qxxx are (x,x0)
c
c  This is the flow due to a ring distribution of stresslets
c
c  Pij is used for the desingularization of the dlp
c
c------------
c
c  Iopt = 1 produces only the Green's function
c  Iopt = 2 produces the Green's function and the stress tensor
c
c================================================

      Implicit Double Precision (a-h,o-z)

      Parameter (eps=0.00000001)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi

c------------------------
c field point on x axis ?
c------------------------

      Iaxis = 0
      If(y0.lt.eps) Iaxis = 1

c-------
c prepare
c-------

      Y2  = Y*Y
      Y02 = Y0*Y0
      YY2 = Y2+Y02
      YYP = Y*Y0
      DY  = Y-Y0
      DY2 = DY*DY

      if(Iopt.ne.1) then
        Y3  = Y*Y2
        Y66 = 6.0D0*Y
        Y03 = Y0*Y02
      end if
 
      if(Iaxis.eq.0) then   ! off the axis
        Y4  = Y2*Y2
        Y04 = Y02*Y02
        YYR = dsqrt(YYP)
        YY3 = YYR**3
        YY5 = YYR**5
        SY  = Y+Y0
        SY2 = SY*SY
      end if

      DX   = X-X0
      DX2  = DX*DX
      DX4  = DX2*DX2
      DR2  = DX2+DY2
      DR   = dsqrt(DR2)
      DXYY = DX2+YY2

      if(Iopt.ne.1) then
        Y6DX  = Y66*DX
        Y6DX2 = Y66*DX2
        DX3   = DX*DX2
        Y6DX3 = Y66*DX3
      end if

c------------------------
      If(Iaxis.eq.0) then       !  point x0 off the axis
c------------------------

        FC1  = 4.0D0*YYP
        FC2  = DX2+SY2
        RK2  = FC1/FC2
        RK   = dsqrt(RK2)
        RK3  = RK2*RK
        RK4  = RK2**2
        RK5  = RK4*RK
        RK2P = 1.0D0-RK2

        call ell_int (RK2,F,E)

c        RJ10 = 2.0D0*RK*F/YYR
c        RJ11 = RK*(DXYY*F-(DX2+SY2)*E)/YY3
c        RJ30 = 2.0D0*RK*E/(YYR*DR2)
c        RJ31 = RK*(-F+DXYY*E/DR2)/YY3
c        RJ32 = RK*(-DXYY*F+(DX4+2.0*DX2*YY2+Y4+Y04)*E/DR2)/YY5

        if(Iopt.ne.1) then

        RL10 = F
        RL30 = E/RK2P
        RK6  = RK4*RK2
c       RK8  = RK6*RK2
        RL50 =  (2.0D0 * (2.0D0 - RK2)*RL30 - F) / (3.0D0*RK2P)
        RL52 =  (RL50 - RL30) / RK2
        RL54 =  (F-RL50+2.0*RK2*RL52) / RK4
        RL56 = -(E-RL50+3.0*RK2*RL52-3.0*RK4*RL54)/ RK6
c       PREP = ( 2.0*(2.0-RK2)*E - RK2P*RK2P*F ) / 3.0D0
c       RL58 = - (PREP - RL50 + 4.0D0 * RK2*RL52 - 6.0D0 * RK4*RL54
c    +                        + 4.0D0 * RK6*RL56 ) / RK8
        FCTR = RK5/(8.0*YY5)
        RJ50 = FCTR * RL50
        RJ51 = FCTR * (2.0D0*RL52-RL50)
        RJ52 = FCTR * (4.0D0*(RL54-RL52)+RL50)
        RJ53 = FCTR * (8.0D0*RL56 - 12.0D0*RL54
     +                +6.0D0*RL52 - RL50)
c       RJ54 = FCTR * (16.0*RL58 - 32.0*RL56 + 24.0*RL54
c    +                 -8.0*RL52 + RL50)
        End If

c-----------
        Else               !  point x0 on the axis
c-----------

        DR3  = DR2*DR
        DR5  = DR3*DR2
c        RJ10 = PI2/DR
c        RJ11 = 0.0D0
c        RJ30 = PI2/DR3
c        RJ31 = 0.0D0
c        RJ32 = PI/DR3

        IF(Iopt.ne.1) then
         RJ50 = PI2/DR5
         RJ51 = 0.0D0
         RJ52 = PI/DR5
         RJ53 = 0.0D0
c        RJ54 = 3.0*PI/(4.0*DR5)
        End If

c-------------
        End If
c-------------

c---
c Build the Green's function
c---

c     SXX = Y *    (  RJ10 + DX2*RJ30)
c      SXY = Y * DX*(Y*RJ30 - Y0 *RJ31)
c      SYX = Y * DX*(Y*RJ31 - Y0 *RJ30)
c      SYY = Y *    (  RJ11 + YY2*RJ31-YYP*(RJ30+RJ32))
 
c------
c Build the stress tensor
c------
      If(Iopt.ne.1) then

        QXXX = - Y6DX3 * RJ50
        QXXY = - Y6DX2 * (Y*RJ50 - Y0*RJ51)
        QXYX =   QXXY
        QXYY = - Y6DX  * (Y02*RJ52+Y2*RJ50 - 2.0*YYP*RJ51)
        QYXX = - Y6DX2 * (Y  *RJ51-Y0*RJ50)
        QYXY = - Y6DX  * ((Y02+Y2)*RJ51 - YYP*(RJ52+RJ50))
        QYYX = QYXY
        QYYY = - Y66 * ( Y3     *  RJ51
     +                 - Y03    *  RJ52
     +                 - Y*YYP  * (RJ50+2.0*RJ52)
     +                 + Y0*YYP * (RJ53+2.0*RJ51) )
c         PXX = QYXX
c         PXY = QYXY
c         PYX = - Y6DX * (Y2*RJ52 - 2.0D0 * YYP*RJ51 + Y02*RJ50)
c         PYY = - Y66  * (  Y0*YYP * (2.0D0 * RJ52+RJ50)
c     +                   - Y03    * RJ51
c     +                   - Y *YYP * (2.0*RJ51+RJ53)
c     +                   + Y3     * RJ52  )

      End If
c------

c-----
c Done
c-----

      Return
      End
