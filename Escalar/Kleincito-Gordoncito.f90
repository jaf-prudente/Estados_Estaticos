program Escalar

    ! Este programa soluciona el sistema de Dirac con rk4

    implicit none

    !---------------------------------------------------------------
    ! Declaramos variables útiles
    integer Nr, i

    real(8) uno, dos, tres, cuatro, ocho
    real(8) medio, tercio, sexto
    real(8) rmax, dr
    real(8) phi1, w

    ! Declaramos variables del método
    real(8) k1m, k1s, k1phi, k1vphi
    real(8) k2m, k2s, k2phi, k2vphi
    real(8) k3m, k3s, k3phi, k3vphi
    real(8) k4m, k4s, k4phi, k4vphi

    ! Declaramos variables que guarden a las funciones derivadas
    real(8)  dm, ds, dphi, dvphi

    ! Declaramos los arreglos que usaremos
    real(8), allocatable, dimension (:) :: r, m, s, phi, vphi

    ! Iniciamos variables
    uno    = 1.0D0
    dos    = 2.0D0
    tres   = 3.0D0
    cuatro = 4.0D0
    ocho   = 8.0D0
    
    medio  = 1.0D0/2.0D0
    tercio = 1.0D0/3.0D0
    sexto  = 1.0D0/6.0D0

    print *, 'Escribe el valor de phi1'
    read(*,*) phi1
    print *

    print *, 'Escribe el valor de w'
    read(*,*) w
    print *

    print *, 'Escribe el valor de rmax'
    read(*,*) rmax
    print *

    print *, 'Escribe el valor de Nr'
    read(*,*) Nr
    print *

    allocate( r(0:Nr), m(0:Nr), s(0:Nr), phi(0:Nr), vphi(0:Nr))

    dr = rmax/dble(Nr)

    !---------------------------------------------------------------
    ! Llenamos el parámetro
    do i=0, Nr
        r(i) = dble(i)*dr - medio*dr
    end do

    !------------------------------------------------
    ! Definimos las condiciones iniciales (con simetrías).
    m(1)    = (dos*tercio)*phi1*(w**2 + uno)*r(1)**3
    s(1)    = uno + ( w*phi1*r(1) )**2
    phi(1)  = phi1 + (medio*phi1)*(uno - w**2)*r(1)**2
    vphi(1) = phi1*(uno - w**2)*r(1)

    m(0)    = -m(1)
    s(0)    =  s(1)
    phi(0)  =  phi(1)
    vphi(0) = -vphi(1)

    !---------------------------------------------------------------
    ! Hacemos el método de Runge-Kutta en sí
    do i=1, Nr-1

        !------------------------------------------------
        call derivadasE(w, r(i), m(i), s(i), phi(i), vphi(i), &
        dm, ds, dphi, dvphi)

        k1m    = dm
        k1s    = ds
        k1phi  = dphi
        k1vphi = dvphi

        !------------------------------------------------
        call derivadasE(w, r(i) + medio*dr, m(i) + medio*k1m*dr, s(i) + medio*k1s*dr, phi(i) + medio*k1phi*dr, &
        vphi(i) + medio*k1vphi*dr, dm, ds, dphi, dvphi)

        k2m    = dm
        k2s    = ds
        k2phi  = dphi
        k2vphi = dvphi

        !------------------------------------------------
        call derivadasE(w, r(i) + medio*dr, m(i) + medio*k2m*dr, s(i) + medio*k2s*dr, phi(i) + medio*k2phi*dr, &
        vphi(i) + medio*k2vphi*dr, dm, ds, dphi, dvphi)
        
        k3m    = dm
        k3s    = ds
        k3phi  = dphi
        k3vphi = dvphi

        !------------------------------------------------
        call derivadasE(w, r(i) + dr, m(i) + k3m*dr, s(i) + k3s*dr, phi(i) + k3phi*dr, vphi(i) + k3vphi*dr, &
        dm, ds, dphi, dvphi)
        
        k4m    = dm
        k4s    = ds
        k4phi  = dphi
        k4vphi = dvphi

        !------------------------------------------------
        m(i+1)    = m(i) + (k1m + dos*k2m + dos*k3m + k4m)*(dr*sexto)
        s(i+1)    = s(i) + (k1s + dos*k2s + dos*k3s + k4s)*(dr*sexto)
        phi(i+1)  = phi(i) + (k1phi + dos*k2phi + dos*k3phi + k4phi)*(dr*sexto)
        vphi(i+1) = vphi(i) + (k1vphi + dos*k2vphi + dos*k3vphi + k4vphi)*(dr*sexto)

    end do

    !---------------------------------------------------------------
    ! Imprimimos en un archivo
    open(10, file = 'sol.dat', form = 'formatted', status = 'replace')

    do i=1, Nr

        write(10,*) r(i), m(i), s(i), phi(i), vphi(i)

    end do

    close(10)

end program

!---------------------------------------------------------------
! Definimos una subrutina para las derivadas.
subroutine derivadasE(w, r, m, s, phi, vphi, dm, ds, dphi, dvphi)

    implicit none

    real(8) w, r, m, s, phi, vphi, N
    real(8) dN, dm, ds, dphi, dvphi
    real(8) uno, dos

    uno    = 1.0D0
    dos    = 2.0D0
    
    N     = uno - (dos*m)/(r)
    dN    = dos*m/r**2 - (dos/r)*( r**2*(N*vphi**2 + phi**2 + (w**2*phi**2)/(N*s**2)) )

    dm    = r**2*(N*vphi**2 + phi**2 + (w**2*phi**2)/(N*s**2))
    ds    = dos*r*s*( vphi**2 + ((w*phi)/(N*s))**2 )
    dphi  = vphi    
    dvphi = -(dos/r + dN/N + ds/s)*vphi - (w**2/(N*s**2) - uno)*phi/N

end subroutine