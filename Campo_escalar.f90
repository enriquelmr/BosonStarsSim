
!Programa que resuelve las ecuaciones de campo de una estrella de bosones
!        para el caso estático
!
!
module inicon
  implicit none  !DEFINO LAS CONDICIONES INICIALES

  real, parameter ::  m0=0.0d0 , r0=0.0d0  , imp = 10.0       !CAMPO Y MASA
  real, parameter:: n0=0.0d0, a0=0.0d0, b0=0.0d0        !#PARTINCULAS Y METRICA
  real, parameter::  pi=acos(-1.0)      !  VALOR DE PI


end module inicon


program campoescalar


          use inicon

          implicit none


          real(8):: rlim, h,  w, w0, f0
          real(kind=8), external :: buscarE, RK4test
          real(kind = 8), allocatable, dimension (:) ::  ff!
          integer :: i, nocamp!

         open(unit=2, file="datos.dat")
         open(unit=1, file="campoescalar.dat")

         rlim = 140.0d0
         !f0 = 0.1d0
         w0 =  1.001d0
         nocamp = 23!
         allocate(ff(nocamp))!
         ff(1) = 0.01d0
         ff(2) = 0.1d0
         ff(3) = 0.2d0
         ff(4) = 0.3d0
         ff(5) = 0.4d0
         ff(6) = 0.5d0
         ff(7) = 0.6d0
         ff(8) = 0.8d0
         ff(9) = 1.0d0
         ff(10) = 1.2d0
         ff(11) = 1.4d0
         ff(12) = 1.6d0
         ff(13) = 1.8d0
         ff(14) = 2.0d0
         ff(15) = 2.2d0
         ff(16) = 2.4d0
         ff(17) = 2.6d0
         ff(18) = 2.8d0
         ff(19) = 3.0d0
         ff(20) = 3.2d0
         ff(21) = 3.4d0
         ff(22) = 3.6d0
         ff(23) = 3.8d0
         do i = 1,nocamp , 1
                           f0 = ff(i)
                    if ( f0 > 0.2d0 .AND. f0 <= 2.0d0) then
                               w0 = 1.2d0
                    elseif ( f0 > 1.1d0 ) then
                               w0 = floor(w)
                    end if


                    if ( f0 > 0.0d0 .AND. f0 < 4.0d0) then
                               h = 1.0d-4
                    elseif ( f0 >= 4.0d0 .AND. f0 < 4.5d0)  then
                              h = 1.0d-5
                  elseif ( f0 >= 4.5d0 ) then
                              h = 1.0d-6
                    end if






          print*, "entrara a buscar la eigenenergia para:", f0
          w = buscarE(rlim, h, w0, f0)

         print*, "el valor de W que se mandara a la subrutina es :", w




             call RK4(rlim, h, w, f0)
             if ( 0.0d0 < f0 .AND. f0 < 0.8 ) then
                        f0 = f0 + 0.01d0
            else
                       f0 = f0 + 0.1d0
             end if

          end do
          close(1)
          close(2)

          call system("gnuplot -p runge.gnu")
          deallocate(ff)

end program campoescalar


subroutine RK4(rlim,h,w,f0)
           use inicon
           implicit none

           real(8):: phi, phi0, nu ,ld, rho, n, m, u, r, rlim,w ,h, f0
           real(8):: k1, k2, k3, k4, l1, l2, l3, l4
           real(8):: knu1, knu2, knu3, knu4
           real(8):: kld1, kld2, kld3, kld4
           real(8):: km1, km2, km3, km4
           real(8):: kn1, kn2, kn3, kn4
           integer :: c
           real(kind=8), external :: f, g, fnu, fld, fm, fn


           r=r0
           phi=f0
           nu=b0
           ld=a0
           m=m0
           n=n0
           !h = h!*1.0d-1

           u=0.0d0
           c=0



do while (r<=rlim)

           if ( floor(imp*1d1) == floor(f0*1d1) ) then
           if ( MOD(c,100)==0 ) then

                      write(1,*) r,phi,exp(nu),exp(ld),m,n
           end if
           c=c+1
           end if


                 k1 = h*f(r, nu, ld, phi, u,w)
                 l1 = h*g(r, nu, ld, phi, u,w)
                 knu1 = h*fnu(r, nu, ld, phi, u,w)
                 kld1 = h*fld(r, nu, ld, phi, u,w)


                 k2 = h*f(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w)
                 l2 = h*g(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w)
                 knu2 = h*fnu(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w)
                 kld2 = h*fld(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w)


                 k3 = h*f(r+0.5d0*h, nu+0.5d0*knu2,ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w)
                 l3 = h*g(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w)
                 knu3 = h*fnu(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,u+0.5d0*l2,w)
                 kld3 = h*fld(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w)


                 k4 = h*f(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)
                 l4 = h*g(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)
                 knu4 = h*fnu(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)
                 kld4 = h*fld(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)

           km1 = h*fm(r, nu, ld, phi,u,w)
           km2 = h*fm(r+0.5*h, nu+0.5*knu1, ld+0.5*kld1, phi+0.5*k1,u+0.5*l1,w)
           km3 = h*fm(r+0.5*h, nu+0.5*knu2, ld+0.5*kld2, phi+0.5*k2,  u+0.5*l2,w)
           km4 = h*fm(r+0.5*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)

           kn1 = h*fn(r, nu, ld, phi, u,w)
           kn2 = h*fn(r+0.5*h, nu+0.5*knu1, ld+0.5*kld1, phi+0.5*k1, u,w)
           kn3 = h*fn(r+0.5*h, nu+0.5*knu2, ld+0.5*kld2, phi+0.5*k2, u,w)
           kn4 = h*fn(r+0.5*h, nu+knu3, ld+kld3, phi+k3, u,w)



           phi = phi + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
           u = u +(1.0d0/6.0d0)*(l1 + 2.0d0*l2 + 2.0*l3 + l4)
           nu = nu + (1.0d0/6.0d0)*(knu1 + 2.0d0*knu2 + 2.0d0*knu3 + knu4)
           ld = ld + (1.0d0/6.0d0)*(kld1 + 2.0d0*kld2 + 2.0d0*kld3 + kld4)
           m = m + (1.0d0/6.0d0)*(km1 + 2.0*km2 + 2.0*km3 + km4)
           n = n + (1.0d0/6.0d0)*(kn1 + 2.0*kn2 + 2.0*kn3 + kn4)
           r=r+h


           if (isnan(phi)) exit
           if (u > 0.0d0 .OR. phi < 0.0d0) exit !(phi <0 .AND.abs(phi) > f0*1.0d-3 .OR. u > 0  .AND.abs(phi) > f0*1.0d-3) exit



end do

           write(2,*) f0, w/exp(nu), exp(nu), m,n, (w**2-1.0d0)*0.5*(f0**2),(w**2+1.0d0)*0.5*(f0**2)

end subroutine

real(kind = 8) function RK4test(rlim,h, w, f0)
           use inicon
           implicit none

           real(8):: phi, nu ,ld, rho, u, r, rlim,w ,h, f0
           real(8):: k1, k2, k3, k4, l1, l2, l3, l4
           real(8):: knu1, knu2, knu3, knu4
           real(8):: kld1, kld2, kld3, kld4
           real(8):: phireal
           real(kind=8), external :: f, g, fnu, fld




           r=r0
           phi=f0
           nu=b0
           ld=a0
           u=0.0d0


do while (r<=rlim)


           k1 = h*f(r, nu, ld, phi, u,w)
          l1 = h*g(r, nu, ld, phi, u,w)
          knu1 = h*fnu(r, nu, ld, phi, u,w)
          kld1 = h*fld(r, nu, ld, phi, u,w)


          k2 = h*f(r+0.5*h, nu+0.5*knu1,ld+0.5*kld1,phi+0.5*k1,u+0.5*l1,w)
          l2 = h*g(r+0.5*h, nu+0.5*knu1,ld+0.5*kld1,phi+0.5*k1,u+0.5*l1,w)
          knu2 = h*fnu(r+0.5*h, nu+0.5*knu1, ld+0.5*kld1, phi+0.5*k1,u+0.5*l1,w)
          kld2 = h*fld(r+0.5*h, nu+0.5*knu1, ld+0.5*kld1, phi+0.5*k1,u+0.5*l1,w)


          k3 = h*f(r+0.5*h, nu+0.5*knu2,ld+0.5*kld2, phi+0.5*k2,  u+0.5*l2,w)
          l3 = h*g(r+0.5*h, nu+0.5*knu2, ld+0.5*kld2, phi+0.5*k2,  u+0.5*l2,w)
          knu3 = h*fnu(r+0.5*h, nu+0.5*knu2, ld+0.5*kld2, phi+0.5*k2,u+0.5*l2,w)
          kld3 = h*fld(r+0.5*h, nu+0.5*knu2, ld+0.5*kld2, phi+0.5*k2,  u+0.5*l2,w)


          k4 = h*f(r+0.5*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)
          l4 = h*g(r+0.5*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)
          knu4 = h*fnu(r+0.5*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)
          kld4 = h*fld(r+0.5*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)





          phi = phi + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
          u = u +(1.0d0/6.0d0)*(l1 + 2.0d0*l2 + 2.0*l3 + l4)
          nu = nu + (1.0d0/6.0d0)*(knu1 + 2.0d0*knu2 + 2.0d0*knu3 + knu4)
          ld = ld + (1.0d0/6.0d0)*(kld1 + 2.0d0*kld2 + 2.0d0*kld3 + kld4)

           r=r+h
           if (isnan(phi)) exit
           if ( abs(phi) > f0*1.0d3 ) exit


                      phireal=phi



end do
           RK4test = phireal
end function


real(kind = 8) function buscarE(rlim, h, w0, f0)

           use inicon
           implicit none

           real(8):: rlim,h, phip, phi0,  wp, w0,f0
           integer:: cont, i
           real(8):: p, p1, t
           real(kind=8), external :: RK4test



           wp = w0
           !print*, wp

                     phip = RK4test(rlim,h, wp, f0)
                     print*, wp, phip
                     cont = 0
                     if ( f0 < 0.3 ) then
                                t = 10.0d1
                     elseif ( f0 >= 0.3d0 .AND.  f0 < 2.0d0 ) then
                                t = 10.0d0
                     elseif ( f0 >= 2.0d0 .AND. f0 < 3.0d0) then
                                 t = 1.0d0
                      elseif ( f0 >= 3.0d0 .AND. f0 < 3.6d0) then
                                 t = 1.0d-1
                      elseif ( f0 >= 3.6d0 .AND. f0 < 4.0d0) then
                                 t = 1.0d-2
                      elseif ( f0 >= 4.0d0 .AND. f0 < 4.6d0) then
                                 t = 1.0d-3
                      elseif ( f0 >= 4.6d0 ) then
                                 t = 1.0d-4
                     end if

                     p = 1.0d0/t


                     do while ( cont < 20)



                                           phi0 = phip

                                           if ( phip > 0 ) then
                                                      wp = wp + p
                                           elseif ( phip < 0 ) then
                                                      wp = wp - p

                                           end if

                                           phip = RK4test(rlim,h, wp, f0)

                                           !print*, phip, wp
                                           if ( phip == phi0 .AND. cont >=1) exit

                                           !if ( abs(phip - phi0) ) exit



                                           if ( phi0*phip < 0 ) then

                                                      cont = cont +1
                                                      !call cpu_time(v(cont))
                                                                  t = t*1.0d1
                                                                  p = 1.0d0/t

                                                      print*, "el contador va en", cont, wp, p !tv(cont)

                                           end if



                     end do

                       buscarE = wp
                       print*, "el valor de W dentro de la funcion para encontrar el valor"
                       print*, buscarE


end function

real(kind = 8) function f(r, nu, ld, phi,u,w)
  implicit none

  real(8):: r, nu, ld, phi,u,w
  f=u
end function

real(kind = 8) function g(r, nu, ld, phi,u,w)

  real(8)::  r, nu, ld, phi,u,w
  real(kind=8), external :: fnu, fld

  if ( r==0 ) then
             g=(1/3.0d0)*phi*(1.-w**2)

  else
           g=u*(fld(r, nu, ld, phi,u,w)-fnu(r, nu, ld, phi,u,w)-(2.0d0/r)) + phi*exp(2.0d0*ld)*(1.0d0-exp(-2.0d0*nu)*w**2)
           !g=1.1*nu-2.15*phi

  end if

end function

real(kind = 8) function fnu(r, nu, ld, phi,u,w)
  implicit none

  real(8):: r, nu, ld, phi,u,w

  if ( r==0 ) then
            fnu=1.0d0
  else

  fnu=(r*0.25d0)*(exp(2.0d0*(ld-nu))*(w**2)*(phi**2)+u**2-exp(2.0d0*ld)*(phi**2))+(exp(2.0d0*ld)-1.0d0)/(2.0d0*r)
 end if
 !fnu=nu
end function

real(kind = 8) function fld(r, nu, ld, phi,u,w)
  implicit none

  real(8):: r, nu, ld, phi,u,w

  if ( r==0 ) then
             fld=1.0
  else
  fld=(0.25d0*r)*(exp(2.0d0*(ld-nu))*(w**2)*(phi**2)+u**2+exp(2.0d0*ld)*(phi**2))-(exp(2.0d0*ld)-1.0d0)/(2.0d0*r)
end if
           !fld=1.1*phi-2.15*nu


end function

real(kind = 8) function fm(r, nu, ld, phi,u,w)
  use inicon
  implicit none

  real(8):: r, nu, ld, phi,u,w
  real(kind=8), external :: rho
  fm=4.0d0*pi*(r**2)*0.5d0*rho(r, nu, ld, phi,u, w)
end function

real(kind = 8) function fn(r, nu, ld, phi,u,w)
  use inicon
  implicit none

  real(8):: r, nu, ld, phi,u,w
  fn=4.0d0*pi*w*(r**2)*exp(ld-nu)*(phi**2)
end function

real(kind = 8) function rho(r, nu, ld, phi, u, E)
           implicit none

           real(8):: r, nu, ld, phi,u,E

           rho=(phi**2)*(1.0d0+(E**2)*exp(-2.0d0*nu))+exp(-2.0d0*ld)*u**2
end function rho
