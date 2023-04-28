!Programa que resuelve las ecuaciones de campo de una estrella de bosones
!        para el caso estático
!
!
module inicon
  implicit none  !DEFINO LAS CONDICIONES INICIALES

  real, parameter ::  m0=0.0d0 , r0=0.0d0  , imp = 0.1       !CAMPO Y MASA
  real, parameter:: n0=0.0d0, a0=0.0d0, b0=0.0d0        !#PARTINCULAS Y METRICA
  real, parameter::  pi=acos(-1.0d0)      !  VALOR DE PI


end module inicon


program campoescalar


          use inicon

          implicit none


          real(8):: rlim, h,  w, w0, sig0, l
          real(kind=8), external :: buscarE, RK4test
          integer:: i






          rlim = 120.0d0
          w0 = 1.2d0
          l= 0.0d0
          sig0 = 0.357d0

         ! do while ( l < 101.0d0 )
                     open(unit=1, file="campoescalar.dat")
                     open(unit=2, file="phi-p.dat")
                     open(unit=3, file="traceT.dat")

                                h = 1.0d-4

                                print*, "entrara a buscar la eigenenergia para:", sig0, "l:", l
                                w = buscarE(rlim, h, w0, sig0, l)

                                print*, "el valor de W que se mandara a la subrutina es :", w




                                call RK4(rlim, h, w, sig0, l)
                                call system("gnuplot -p runge.gnu")
                                w0 = w
           !                     l = l + 10.0d0


          close(1)
          close(2)
          close(3)
           !end do



end program campoescalar


subroutine RK4(rlim, h, w, sig0, l)
           use inicon
           implicit none

           real(8):: phi, phi0, nu ,ld, rho, n, m, u, r, rlim,w ,h, sig0, l
           real(8):: k1, k2, k3, k4, l1, l2, l3, l4
           real(8):: knu1, knu2, knu3, knu4
           real(8):: kld1, kld2, kld3, kld4
           real(8):: km1, km2, km3, km4
           real(8):: kn1, kn2, kn3, kn4
           real(kind=8), external :: f, g, fnu, fld, fm
           integer:: c, uc, phicont
           real(8):: u2, phipi, tem

           r=r0
           phi=sig0
           nu=b0
           ld=a0
           m=m0
           n=n0
           !h = h!*1.0d-1

           u=0.0d0
           c = 0
           uc = 0
           phicont = 0
           u2 = u
           phi0 = phi
           phipi = 0.0d0
           tem = 0.0d0




do while (r<=rlim)


           if ( MOD(c,100)==0 ) then

                      write(1,*) r, phi, u!exp(nu),exp(ld),m,n
                      write(2,*) r, phipi
                      write(3,*) r, tem

           end if
           c=c+1


                 k1 = h*f(r, nu, ld, phi, u,w,l)
                 l1 = h*g(r, nu, ld, phi, u,w,l)
                 knu1 = h*fnu(r, nu, ld, phi, u,w,l)
                 kld1 = h*fld(r, nu, ld, phi, u,w,l)


                 k2 = h*f(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w,l)
                 l2 = h*g(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w,l)
                 knu2 = h*fnu(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,l)
                 kld2 = h*fld(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,l)


                 k3 = h*f(r+0.5d0*h, nu+0.5d0*knu2,ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,l)
                 l3 = h*g(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,l)
                 knu3 = h*fnu(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,u+0.5d0*l2,w,l)
                 kld3 = h*fld(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,l)


                 k4 = h*f(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,l)
                 l4 = h*g(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,l)
                 knu4 = h*fnu(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,l)
                 kld4 = h*fld(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,l)

           km1 = h*fm(r, nu, ld, phi,u,w,l)
           km2 = h*fm(r+0.5*h, nu+0.5*knu1, ld+0.5*kld1, phi+0.5*k1,u+0.5*l1,w,l)
           km3 = h*fm(r+0.5*h, nu+0.5*knu2, ld+0.5*kld2, phi+0.5*k2,  u+0.5*l2,w,l)
           km4 = h*fm(r+0.5*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,l)





           phi = phi + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
           u = u +(1.0d0/6.0d0)*(l1 + 2.0d0*l2 + 2.0*l3 + l4)
           nu = nu + (1.0d0/6.0d0)*(knu1 + 2.0d0*knu2 + 2.0d0*knu3 + knu4)
           ld = ld + (1.0d0/6.0d0)*(kld1 + 2.0d0*kld2 + 2.0d0*kld3 + kld4)
           m = m + (1.0d0/6.0d0)*(km1 + 2.0*km2 + 2.0*km3 + km4)

           phipi = (2.0d0/exp(nu))*(w**2)*phi**2-phi**2 - 0.5d0*l*phi**4
           tem = (1.0d0/exp(nu))*(w**2)*phi**2 - 2.0d0*phi**2 - l*phi**4 - (1.0d0/exp(ld))*u**2
           r=r+h
           if ( u*u2 < 0.0d0 ) then
                      uc = uc +1
           end if
           if ( phi*phi0 < 0.0d0 ) then
                      phicont = phicont +1
           end if
           if (isnan(phi) .OR. isnan(u)) exit
           !if ( 2 < uc) exit
           !if ( 2 <= phicont) exit


           if(phi <0 .AND.abs(phi) > sig0*1.0d-3 ) exit
           if ( r > 10.0 .AND. u > 0.0d0 .AND. abs(phi) > sig0*1.0d-3) exit
           u2 = u
           phi0 = phi


end do
print*, m

end subroutine

real(kind = 8) function RK4test(rlim,h, w, sig0, l)
           use inicon
           implicit none

           real(8):: phi, nu ,ld, rho, u, r, rlim,w ,h, sig0, l
           real(8):: k1, k2, k3, k4, l1, l2, l3, l4
           real(8):: knu1, knu2, knu3, knu4
           real(8):: kld1, kld2, kld3, kld4
           real(8):: phireal, u2
           real(kind=8), external :: f, g, fnu, fld
           integer:: uc, phicont




           r=r0
           phi=sig0
           nu=b0
           ld=a0
           u=0.0d0
           u2 = u
           uc = 0
           phireal = phi
           phicont = 0


do while (r<=rlim)


           k1 = h*f(r, nu, ld, phi, u,w,l)
           l1 = h*g(r, nu, ld, phi, u,w,l)
           knu1 = h*fnu(r, nu, ld, phi, u,w,l)
           kld1 = h*fld(r, nu, ld, phi, u,w,l)


           k2 = h*f(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w,l)
           l2 = h*g(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w,l)
           knu2 = h*fnu(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,l)
           kld2 = h*fld(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,l)


           k3 = h*f(r+0.5d0*h, nu+0.5d0*knu2,ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,l)
           l3 = h*g(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,l)
           knu3 = h*fnu(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,u+0.5d0*l2,w,l)
           kld3 = h*fld(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,l)


           k4 = h*f(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,l)
           l4 = h*g(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,l)
           knu4 = h*fnu(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,l)
           kld4 = h*fld(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,l)




          phi = phi + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
          u = u +(1.0d0/6.0d0)*(l1 + 2.0d0*l2 + 2.0*l3 + l4)
          nu = nu + (1.0d0/6.0d0)*(knu1 + 2.0d0*knu2 + 2.0d0*knu3 + knu4)
          ld = ld + (1.0d0/6.0d0)*(kld1 + 2.0d0*kld2 + 2.0d0*kld3 + kld4)

           r=r+h
           if ( u*u2 < 0.0d0 ) then
                      uc = uc +1
           end if
           if ( phi*phireal < 0.0d0 ) then
                      phicont = phicont +1
           end if

           if (isnan(phi)) exit
           if (isnan(u) .OR. 1.0d4 < abs(u) ) exit
           if ( abs(phi) > sig0*1.0d3 ) exit
           if ( 4 <= uc) exit
           if ( 2 <= phicont) exit


                      phireal=phi
                      u2 = u



end do
           RK4test = phireal
end function


real(kind = 8) function buscarE(rlim, h, w0, sig0 ,l)

           use inicon
           implicit none

           real(8):: rlim,h, phip, phi0,  wp, w0,sig0, l
           integer:: cont, i
           real(8):: p
           real(kind=8), external :: RK4test



           wp = w0
           !print*, wp

                     phip = RK4test(rlim,h, wp, sig0 ,l)
                     print*, wp, phip
                     cont = 0

                     p = 1.0d-3



                     do while ( cont < 20)



                                           phi0 = phip

                                           if ( phip > 0 ) then
                                                      wp = wp + p
                                           elseif ( phip < 0 ) then
                                                      wp = wp - p

                                           end if

                                           phip = RK4test(rlim,h, wp, sig0,l)

                                           !print*,  wp, phip
                                           if ( phip == phi0 .AND. cont >=1) exit


                                           if ( phi0*phip < 0 ) then

                                                      cont = cont +1


                                                                  p = p*1.0d-1

                                                      print*, "el contador va en", cont, wp, p !tv(cont)
                                           end if
                                           !end if



                     end do

                       buscarE = wp
                       print*, "el valor de W dentro de la funcion para encontrar el valor"
                       print*, buscarE


end function

real(kind = 8) function f(r, nu, ld, phi,u,w,l)
  implicit none

  real(8):: r, nu, ld, phi,u,w,l
  f=u
end function

real(kind = 8) function g(r, nu, ld, phi,u,w,l)

  real(8)::  r, nu, ld, phi,u,w,l
  real(kind=8), external :: fnu, fld

  if ( r==0 ) then
             g=(1.0d0/3.0d0)*(phi*(1.-w**2)+l*phi**3)

  else
           g=u*(fld(r, nu, ld, phi,u,w,l)-fnu(r, nu, ld, phi,u,w,l)-(2.0d0/r)) + phi*exp(2.0d0*ld)*(1.0d0-exp(-2.0d0*nu)*w**2)&
             +l*phi**3
           !g=1.1*nu-2.15*phi

  end if

end function

real(kind = 8) function fnu(r, nu, ld, phi,u,w,l)
  implicit none

  real(8):: r, nu, ld, phi,u,w,l

  if ( r==0 ) then
            fnu=0.0d0
  else

  fnu=(r*0.5d0)*(exp(2.0d0*(ld-nu))*(w**2)*(phi**2)+u**2-exp(2.0d0*ld)*(phi**2)-0.5d0*l*exp(2.0d0*ld)*phi**4)&
  +(exp(2.0d0*ld)-1.0d0)/(2.0d0*r)
 end if
 !fnu=nu
end function

real(kind = 8) function fld(r, nu, ld, phi, u, w ,l)
  implicit none

  real(8):: r, nu, ld, phi,u,w,l

  if ( r==0 ) then
             fld=0.
  else
  fld=(0.25d0*r)*(exp(2.0d0*(ld-nu))*(w**2)*(phi**2)+u**2+exp(2.0d0*ld)*(phi**2)+0.5d0*l*exp(2.0d0*ld)*phi**4)&
  -(exp(2.0d0*ld)-1.0d0)/(2.0d0*r)
end if
           !fld=1.1*phi-2.15*nu


end function

real(kind = 8) function fm(r, nu, ld, phi,u,w ,l)
  use inicon
  implicit none

  real(8):: r, nu, ld, phi, u, w, l
  real(kind=8), external :: rho
  fm = 0.5d0*(r**2)*( (exp(-2.0d0*nu)*w**2+1)*phi**2 + 0.5d0*l*phi**4+exp(-2.0d0*ld)*u**2)
end function
