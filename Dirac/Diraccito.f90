program Dirac

    ! Este programa soluciona el sistema de Dirac con rk4

    implicit none

    !---------------------------------------------------------------
    ! Declaramos variables útiles
    integer Nr, i

    real(8) uno, dos, tres, cuatro, ocho
    real(8) medio, tercio, sexto
    real(8) rmax, dr
    real(8) F1, w

    ! Declaramos variables del método
    real(8) k1m, k1s, k1F, k1G
    real(8) k2m, k2s, k2F, k2G
    real(8) k3m, k3s, k3F, k3G
    real(8) k4m, k4s, k4F, k4G

    ! Declaramos variables que guarden a las funciones derivadas
    real(8)  dm, ds, dF, dG

    ! Declaramos los arreglos que usaremos
    real(8), allocatable, dimension (:) :: r, m, s, F, G

    ! Iniciamos variables
    uno    = 1.0D0
    dos    = 2.0D0
    tres   = 3.0D0
    cuatro = 4.0D0
    ocho   = 8.0D0
    
    medio  = 1.0D0/2.0D0
    tercio = 1.0D0/3.0D0
    sexto  = 1.0D0/6.0D0

    print *, 'Escribe el valor de F1'
    read(*,*) F1
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

    allocate( r(0:Nr), m(0:Nr), s(0:Nr), F(0:Nr), G(0:Nr))

    dr = rmax/dble(Nr)

    !---------------------------------------------------------------
    ! Llenamos el parámetro
    do i=0, Nr
        r(i) = dble(i)*dr - medio*dr
    end do

    !------------------------------------------------
    ! Definimos las condiciones iniciales (con simetrías).
    m(1) = (dos*tercio)*(F1**2)*(r(1)**3)*w
    s(1) = uno + tercio*(F1**2)*(cuatro*w - uno)*(r(1)**2)
    F(1) = F1*r(1)
    G(1) = tercio*F1*(w - uno)*(r(1)**2)

    m(0) = -m(1)
    s(0) =  s(1)
    F(0) = -F(1)
    G(0) =  G(1)

    !---------------------------------------------------------------
    ! Hacemos el método de Runge-Kutta en sí
    do i=1, Nr-1

        !------------------------------------------------
        call derivadasD(w, r(i), m(i), s(i), F(i), G(i), &
        dm, ds, dF, dG)

        k1m = dm
        k1s = ds
        k1F = dF
        k1G = dG

        !------------------------------------------------
        call derivadasD(w, r(i) + medio*dr, m(i) + medio*k1m*dr, s(i) + medio*k1s*dr, F(i) + medio*k1F*dr, G(i) + medio*k1G*dr, &
        dm, ds, dF, dG)

        k2m = dm
        k2s = ds
        k2F = dF
        k2G = dG

        !------------------------------------------------
        call derivadasD(w, r(i) + medio*dr, m(i) + medio*k2m*dr, s(i) + medio*k2s*dr, F(i) + medio*k2F*dr, G(i) + medio*k2G*dr, &
        dm, ds, dF, dG)
        
        k3m = dm
        k3s = ds
        k3F = dF
        k3G = dG

        !------------------------------------------------
        call derivadasD(w, r(i) + dr, m(i) + k3m*dr, s(i) + k3s*dr, F(i) + k3F*dr, G(i) + k3G*dr, &
        dm, ds, dF, dG)
        
        k4m = dm
        k4s = ds
        k4F = dF
        k4G = dG

        !------------------------------------------------
        m(i+1) = m(i) + (k1m + dos*k2m + dos*k3m + k4m)*(dr*sexto)
        s(i+1) = s(i) + (k1s + dos*k2s + dos*k3s + k4s)*(dr*sexto)
        F(i+1) = F(i) + (k1F + dos*k2F + dos*k3F + k4F)*(dr*sexto)
        G(i+1) = G(i) + (k1G + dos*k2G + dos*k3G + k4G)*(dr*sexto)

    end do

    !---------------------------------------------------------------
    ! Imprimimos en un archivo
    open(10, file = 'sol.dat', form = 'formatted', status = 'replace')

    do i=1, Nr

        write(10,*) r(i), m(i), s(i), F(i), G(i)

    end do

    close(10)

end program

!---------------------------------------------------------------
! Definimos una subrutina para las derivadas.
subroutine derivadasD(w, r, m, s, F, G, dm, ds, dF, dG)

    implicit none

    real(8) w, r, m, s, F, G, N
    real(8) dm, ds, dF, dG
    real(8) uno, dos, cuatro

    uno    = 1.0D0
    dos    = 2.0D0
    cuatro = 4.0D0
    
    N = uno - (dos*m)/(r)

    dm = ( dos*w*( F**2 + G**2 ))/(s)
    ds = ( cuatro*w/(N*r) )*( F**2 + G**2 ) - ( dos*s/(r*sqrt(N)) )*( dos*F*G/r + F**2 - G**2 )
    dF = ( F*G/(r**2*sqrt(N)) + (F**2 - G**2)/(dos*r*sqrt(N)) - m/(r**2*N) + uno/(r*sqrt(N)) )*F - ( uno/(sqrt(N)) + w/(N*s) )*G
    dG = ( F*G/(r**2*sqrt(N)) + (F**2 - G**2)/(dos*r*sqrt(N)) - m/(r**2*N) - uno/(r*sqrt(N)) )*G - ( uno/(sqrt(N)) - w/(N*s) )*F

end subroutine