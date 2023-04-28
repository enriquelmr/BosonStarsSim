
!Programa que resuelve las ecuaciones de campo de una estrella de bosones
!        para el caso estático
!
!
module inicon
  implicit none  !DEFINO LAS CONDICIONES INICIALES

  real(8), parameter ::  m0=0.0d0 , r0=0.0d0  , imp = 10.0       !CAMPO Y MASA
  real(8), parameter:: n0=0.0d0, a0=1.0d0, bb0=1.0d0        !#PARTINCULAS Y METRICA
  real(8), parameter::  pi=acos(-1.0)      !  VALOR DE PI
  real(8), parameter:: f0=0.0d0     !derivada del campoescalar
  real(8), parameter:: ml=1d0
  real(8), parameter:: B0=0.05d0, uB0=0.0d0


end module inicon


program campoescalar


          use inicon

          implicit none


          real(8):: rlim, h,  w, w0, u0, l
          real(kind=8), external :: buscarE, RK4test
          integer:: i



         rlim = 100d0

         w0 = 0.1d0!0.8868050452218615d0
         h = 1.0d-3
        l=0 !l es el parametro de autointeraccion
         open(unit=2, file="datos.dat")
         !open(unit=1, file="campoescalar.dat")

         do i = 5, 80
         u0 = (i*1d0)/20d0

          print*, "entrara a buscar la eigenenergia para:", u0
          w =  buscarE(rlim, h, w0, u0,l)

         print*, "el valor de W que se mandara a la subrutina es :", w

             call RK4(rlim, h, w, u0,l)
             write(1,*)  "      "
           end do

          !close(1)
          close(2)

          call system("gnuplot -p runge.gnu")


end program campoescalar


subroutine RK4(rlim,h,w,u0,l)
           use inicon
           implicit none

           real(8):: phi, phi0, nu ,ld, rho, n, m, u, r, rlim,w ,h, u0
           real(8):: k1, k2, k3, k4, l1, l2, l3, l4
           real(8):: knu1, knu2, knu3, knu4
           real(8):: kld1, kld2, kld3, kld4
           real(8):: km1, km2, km3, km4
           real(8):: kn1, kn2, kn3, kn4
           integer :: c
           real(kind=8), external :: f, g, fnu, fld, fm, fn
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----ROTACION-----!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            real(kind=8), external :: X, Y
            real(8):: kb1, kb2, kb3, kB4
            real(8):: lb1, lb2, lb3, lb4
            real(8):: B,uB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!------Autointeraccion--------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            real(8):: l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            B   = B0
            uB  = uB0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           r=r0
           phi=f0
           nu=bb0
           ld=a0
           m=m0
           n=n0
           !h = h!*1.0d-1

           u=u0
           c=0



do while (r<=rlim)


           !if ( MOD(c,100)==0 ) then
            ! write(1,*) r ,phi, nu,ld, B

           !end if
           c=c+1

            !REciordar que l es un paraetro no una variable

                 k1 = h*f(r, nu, ld, phi, u,w)
                 l1 = h*g(r, nu, ld, phi, u,w,B,l)
                 knu1 = h*fnu(r, nu, ld, phi, u,w,B,l)
                 kld1 = h*fld(r, nu, ld, phi, u,w,B,l)
                 !!!!!!!!!!!!!!!!rotacion!!!!!!!!!!!
                 kb1  = h*X(uB)
                 lb1  = h*Y(r, nu, ld, phi, u,w,B,uB,l)


                 k2 = h*f(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w)
                 l2 = h*g(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w,B+0.5*kb1,l)
                 knu2 = h*fnu(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,B+0.5*kb1,l)
                 kld2 = h*fld(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,B+0.5*kb1,l)
                 !!!!!!!!!!!!!!!!rotacion!!!!!!!!!!!
                 kb2  = h*X(uB+0.5*lb1)
                 lb2  = h*Y(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,B+0.5*kb1,uB+0.5*lb1,l)


                 k3 = h*f(r+0.5d0*h, nu+0.5d0*knu2,ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w)
                 l3 = h*g(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,B+0.5*kb2,l)
                 knu3 = h*fnu(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,u+0.5d0*l2,w,B+0.5*kb2,l)
                 kld3 = h*fld(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,B+0.5*kb2,l)

                 !!!!!!!!!!!!!!!!rotacion!!!!!!!!!!!
                 kb3  = h*X(uB+0.5*lb2)
                 lb3  = h*Y(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,B+0.5*kb2,uB+0.5*lb2,l)


                 k4 = h*f(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)
                 l4 = h*g(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,B+kb3,l)
                 knu4 = h*fnu(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,B+kb3,l)
                 kld4 = h*fld(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,B+kb3,l)
                 !!!!!!!!!!!!!!!!rotacion!!!!!!!!!!!
                 kb4  = h*X(uB+lb3)
                 lb4  = h*Y(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,B+kb3,uB+lb3,l)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MASA Y NUMERO DE PARTICULAS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!FUNCIONES INDEPENDIENTES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           km1 = h*fm(r, nu, ld, phi,u,w,B,m,l)
           km2 = h*fm(r+0.5*h, nu+0.5*knu1, ld+0.5*kld1, phi+0.5*k1,u+0.5*l1,w,B,m,l)
           km3 = h*fm(r+0.5*h, nu+0.5*knu2, ld+0.5*kld2, phi+0.5*k2,  u+0.5*l2,w,B,m,l)
           km4 = h*fm(r+0.5*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,B,m,l)

           kn1 = h*fn(r, nu, ld, phi, u,w,B,l)
           kn2 = h*fn(r+0.5*h, nu+0.5*knu1, ld+0.5*kld1, phi+0.5*k1, u,w,B,l)
           kn3 = h*fn(r+0.5*h, nu+0.5*knu2, ld+0.5*kld2, phi+0.5*k2, u,w,B,l)
           kn4 = h*fn(r+0.5*h, nu+knu3, ld+kld3, phi+k3, u,w,B,l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (u0 .lt. 1.4d0 ) then
        if ( r .lt. 14d0 ) then
          phi = phi + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
        else
          phi = 0d0
        end if

      else

        if ( r .lt. 6d0 ) then
          phi = phi + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
        else
          phi = 0d0
        end if

      end if


      if ( r < 1d3 ) then
        ld = ld + (1.0d0/6.0d0)*(kld1 + 2.0d0*kld2 + 2.0d0*kld3 + kld4)
      else
        ld= 1d0

      end if

           u = u +(1.0d0/6.0d0)*(l1 + 2.0d0*l2 + 2.0*l3 + l4)
           nu = nu + (1.0d0/6.0d0)*(knu1 + 2.0d0*knu2 + 2.0d0*knu3 + knu4)

           m = m + (1.0d0/6.0d0)*(km1 + 2.0*km2 + 2.0*km3 + km4)
           n = n + (1.0d0/6.0d0)*(kn1 + 2.0*kn2 + 2.0*kn3 + kn4)
           r=r+h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ROTACION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          B=B+ (1.0d0/6.0d0)*(kB1 + 2.0d0*kB2 + 2.0d0*kB3 + kB4)
          uB = uB +(1.0d0/6.0d0)*(lB1 + 2.0d0*lB2 + 2.0*lB3 + lB4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          ! if (isnan(phi))  exit
           !if (u > 100.0d0 .OR. phi < 0.0d0) exit !(phi <0 .AND.abs(phi) > f0*1.0d-3 .OR. u > 0  .AND.abs(phi) > f0*1.0d-3) exit
           !if (r>40.0 .AND. abs(phi)>10.0) exit


end do

           write(2,*) u0 , w/nu ,  m,n

end subroutine

real(kind = 8) function RK4test(rlim,h, w, u0, ind,l)
           use inicon
           implicit none

           real(8):: phi, nu ,ld, rho, u, r, rlim,w ,h, uCon, u0
           real(8):: k1, k2, k3, k4, l1, l2, l3, l4
           real(8):: knu1, knu2, knu3, knu4
           real(8):: kld1, kld2, kld3, kld4
           real(8):: phireal
           real(kind=8), external :: f, g, fnu, fld
           integer, external :: signo
           integer :: s, ind !indicador de la variable que se solicita 's' es por el signo
           integer :: c
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----ROTACION-----!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     real(kind=8), external :: X, Y
                     real(8):: kb1, kb2, kb3, kB4
                     real(8):: lb1, lb2, lb3, lb4
                     real(8):: B,uB
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!------Autointeraccion--------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          real(8):: l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     B   = B0
                     uB  = uB0
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



           r=r0
           phi=f0
           nu=bb0
           ld=a0
           u=u0

           s = 0
           c  = 0

do while (r<=rlim)


  !REciordar que l es un paraetro no una variable

       k1 = h*f(r, nu, ld, phi, u,w)
       l1 = h*g(r, nu, ld, phi, u,w,B,l)
       knu1 = h*fnu(r, nu, ld, phi, u,w,B,l)
       kld1 = h*fld(r, nu, ld, phi, u,w,B,l)
       !!!!!!!!!!!!!!!!rotacion!!!!!!!!!!!
       kb1  = h*X(uB)
       lb1  = h*Y(r, nu, ld, phi, u,w,B,uB,l)


       k2 = h*f(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w)
       l2 = h*g(r+0.5d0*h, nu+0.5d0*knu1,ld+0.5d0*kld1,phi+0.5d0*k1,u+0.5d0*l1,w,B+0.5*kb1,l)
       knu2 = h*fnu(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,B+0.5*kb1,l)
       kld2 = h*fld(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,B+0.5*kb1,l)
       !!!!!!!!!!!!!!!!rotacion!!!!!!!!!!!
       kb2  = h*X(uB+0.5*lb1)
       lb2  = h*Y(r+0.5d0*h, nu+0.5d0*knu1, ld+0.5d0*kld1, phi+0.5d0*k1,u+0.5d0*l1,w,B+0.5*kb1,uB+0.5*lb1,l)


       k3 = h*f(r+0.5d0*h, nu+0.5d0*knu2,ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w)
       l3 = h*g(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,B+0.5*kb2,l)
       knu3 = h*fnu(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,u+0.5d0*l2,w,B+0.5*kb2,l)
       kld3 = h*fld(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,B+0.5*kb2,l)

       !!!!!!!!!!!!!!!!rotacion!!!!!!!!!!!
       kb3  = h*X(uB+0.5*lb2)
       lb3  = h*Y(r+0.5d0*h, nu+0.5d0*knu2, ld+0.5d0*kld2, phi+0.5d0*k2,  u+0.5d0*l2,w,B+0.5*kb2,uB+0.5*lb2,l)


       k4 = h*f(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w)
       l4 = h*g(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,B+kb3,l)
       knu4 = h*fnu(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,B+kb3,l)
       kld4 = h*fld(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,B+kb3,l)
       !!!!!!!!!!!!!!!!rotacion!!!!!!!!!!!
       kb4  = h*X(uB+lb3)
       lb4  = h*Y(r+0.5d0*h, nu+knu3, ld+kld3, phi+k3, u+l3,w,B+kb3,uB+lb3,l)




          phi = phi + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
          uCon = u
          u = u +(1.0d0/6.0d0)*(l1 + 2.0d0*l2 + 2.0*l3 + l4)
          nu = nu + (1.0d0/6.0d0)*(knu1 + 2.0d0*knu2 + 2.0d0*knu3 + knu4)
          ld = ld + (1.0d0/6.0d0)*(kld1 + 2.0d0*kld2 + 2.0d0*kld3 + kld4)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ROPTACION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          B=B+ (1.0d0/6.0d0)*(kB1 + 2.0d0*kB2 + 2.0d0*kB3 + kB4)
          uB = uB +(1.0d0/6.0d0)*(lB1 + 2.0d0*lB2 + 2.0*lB3 + lB4)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           r=r+h
           if (isnan(phi) .or. abs(phi) .gt. 2.0) exit
           if (r>20.0 .AND. abs(phi)>4.0) exit
           !!!!!################### el cambio en la pendiente##########################

           if ( MOD(c,100)==0 ) then
             phireal = phi


           end if
           c=c+1
           if (  signo(uCon , u) .EQ. 1 .AND. r .gt. 1d-1) then
             s = s + 1
           end if






end do
          if ( ind .EQ. 0) then
           RK4test = phireal
         elseif ( ind .EQ. 1) then
           RK4test = s*1d0
         end if



end function


real(kind = 8) function buscarE(rlim, h, w0, u0,l)

           use inicon
           implicit none

           real(8):: rlim,h, phip, phi0,  wp, w0, u0
           integer:: cont, i, ii, contaux !ii es el contador de los cambios de la derivada
           real(8):: p, p1, t
           real(kind=8), external :: RK4test
           !!!!!!!____Autointeraccion---------!!!!!!!!!!!!!
           real(8):: l



           wp = w0
           !print*, wp

                     phip = RK4test(rlim,h, wp, u0, 0,l)
                     ii = RK4test(rlim,h, wp, u0, 1,l)
                     print*, wp, phip
                     cont = 0
                     contaux = 0


                     p = 1.0d-1

                     print*, p

                     do while ( cont < 30)


!!!!!!!!!!!pemndiente, tomar cada 100 valores de phi como en la subrutina rk4
!!!!!!!!!! no funciona para valores grandes de la pendiente
!!! razon 1.- cruza por el cero, mandar una alerta

                                           if ( ii .LT. 0.98) then
                                             wp = wp + p
                                             !print*, "entro al primer if"
elseif ( 0.99d0 .gt. ii .or. ii .lt. 2.1d0 .OR. ii .gt. 2.1d0 .and. ii .lt. 5.0 .and. contaux .gt. 0 ) then !
                                                                                      !.OR. ii .gt. 2.1d0 .and. ii .lt. 5.0 .and. contaux .gt. 0 condicion para ,0.3
                                                                                    !.or. ii .lt. 2.1d0 ! condiciones para u0 .gt. 0.5
                                             contaux = 1
                                             !print*, "estoy dentro de doble valle"
                                             if (phip .gt. 0.0d0) then
                                               wp = wp + p
                                              ! print*, "estoy arriba"
                                            elseif ( phip .lt. 0d0) then
                                               wp = wp - p
                                               p = p*1d-1
                                               !print*, "estoy abajo"
                                             end if


                                           elseif( ii .gt. 2.1  )then !.AND. contaux .EQ. 0
                                                wp = wp - 10d0*p

                                                p = p*1d-1
                                                !print*, "entre al ultimo if"
                                              end if
                                            ii = RK4test(rlim,h, wp, u0, 1,l)
                                            phi0 = phip
                                            phip = RK4test(rlim,h, wp, u0, 0,l)

                                            if ( phi0*phip < 0 ) then

                                                       cont = cont +1
                                                       !print*, "el contador va en", cont, wp
                                            end if
                                            !if ( wp .gt. 1.0) exit
                                            if ( phip == phi0 .AND. cont >=0) exit

                                            !print*, wp, ii, phip







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

real(kind = 8) function g(r, nu, ld, phi,u,w,B,l)
  use inicon
  implicit none
  real(8)::  r, nu, ld, phi,u,w,B
  real(kind=8), external :: fnu, fld
  real(8):: l

  if ( r==0 ) then
             g=0d0!(1.0/3.0d0)*phi*(1.-w**2)

  else
            g=u*((fld(r, nu, ld, phi,u,w,B,l)/ld)-(fnu(r, nu, ld, phi,u,w,B,l)/nu)-(2.0d0/r)) + phi*(ld**2)*(1.0d0-((w/nu)**2))&
            +phi*(ld**2)*(((ml*(ml+1.0))/(r**2))-(2.0*B*ml*w)/nu)
           !g=1.1*nu-2.15*phi

  end if

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----ROTACION-----!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind = 8) function X(uB)

  implicit none
  real(8):: uB
      X = uB
end function
real(kind = 8) function Y(r, nu, ld, phi,u,w,B,uB,l)
  use inicon
  implicit none
  real(8)::  r, nu, ld, phi,u,w,B,uB , l
  real(kind=8), external :: fnu, fld

  if ( r == 0d0) then
    Y = (ml*phi**2)*(B0*ml-w)
  elseif (r> 0d0 .and. r .lt. 100000d0) then
    Y = uB*(fld(r, nu, ld, phi,u,w,B,l)/ld+fnu(r, nu, ld, phi,u,w,B,l)/nu-4d0/r) &
    -(((2*ld**2)*ml*(phi**2)*(w-B*ml))/(r**2))
  else
    Y = -(4d0*uB)/r
  end if
  end function

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind = 8) function fnu(r, nu, ld, phi,u,w,B,l)

  use inicon
  implicit none
  real(8):: r, nu, ld, phi,u,w,B,l

  if ( r==0d0 ) then
            fnu=0.0d0
  elseif( r < 50)then
  fnu=nu*((r*0.25d0)*(((ld/nu)**2)*(w**2)*(phi**2)+u**2-(ld**2)*(phi**2)+((2.0*B*ml*w*phi**2)/(nu)**2))+(0.5d0)*((ld**2)-1.0d0)/r)
    else
      fnu=0d0
 end if
 !fnu=nu
end function

real(kind = 8) function fld(r, nu, ld, phi,u,w,B,l)
  use inicon
  implicit none
  real(8):: r, nu, ld, phi,u,w,B,l

  if ( r==0 ) then
             fld=0.0
  elseif( r < 40d0)then
  fld=ld*((0.25d0*r)*(((ld/nu)**2)*(w**2)*(phi**2)+u**2+(ld**2)*(phi**2)+((2d0*B*ml*w*phi**2)/(nu**2d0)))-(0.5d0)*((ld**2)-1d0)/r)
else
  fld = -ld*(0.5d0)*((ld**2)-1d0)/r

end if
           !fld=1.1*phi-2.15*nu


end function

real(kind = 8) function fm(r, nu, ld, phi,u,w,B,m,l)
  use inicon
  implicit none

  real(8):: r, nu, ld, phi,u,w, B,m,l
  real(kind=8), external :: rho, fnu, fld
  if ( r==0 ) then
    fm=0
  else
    fm=(r*fld(r, nu, ld, phi,u,w,B,l))/ld**3+m/r
  end if

end function

real(kind = 8) function fn(r, nu, ld, phi,u,w,B,l)
  use inicon
  implicit none

  real(8):: r, nu, ld, phi,u,w, B,l
  fn=0.5*(r**2)*(ld/nu)*(phi**2)*(w+ml*B)
end function

real(kind = 8) function rho(r, nu, ld, phi, u, w,B)
            use inicon
           implicit none

           real(8):: r, nu, ld, phi,u,w, B
           if ( r == 0d0 ) then
             rho =0d0
           else
             rho=(phi**2)*(1.0d0+w**2/nu**2-(6d0*B*ml*w)/nu**2+ml/r)+u**2/ld**2
           end if


end function rho

!!!!#################################################################################
!!!!#################################################################################
integer function signo (a,b)
  real(8):: a,b

  if ( a*b .GT. 0.0 ) then
    signo = 0
  elseif ( a*b .LE. 0.0) then
    signo = 1
  end if
end function
