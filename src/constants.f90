module constants
    integer, parameter :: kind = 8, ckind = 8
    real(kind),parameter :: G_Au = 2.959e-4 !Au^3 /(day*Msun)
    real(kind), parameter :: G=G_Au
    real(kind),parameter :: M_sun = 1.9885e30 ! kg
    real(kind),parameter :: M_earth = 3e-6 ! solar masses
    real(kind),parameter :: M_moon = 3.69e-8 ! solar masses 
    real(kind),parameter :: M_Jup = 9.55e-3 !10 jupiter 
    real(kind), parameter :: AU =1.496e8
    real(kind), parameter :: pi = acos(-1.0_kind)
end module