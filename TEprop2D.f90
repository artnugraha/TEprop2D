!*******************************************************************************
!*******************************************************************************
! Project      : TEprop2D.f90
!===============================================================================
! Purpose      :
! Calculate thermoelectric properties of 2D materials
!-------------------------------------------------------------------------------
! Authors      : ART Nugraha  (nugraha@flex.phys.tohoku.ac.jp)
!                N. T. Hung   (nguyen@flex.phys.tohoku.ac.jp)
! Started      : 2016.12.12
! Latest Vers. : 2018.02.06
!-------------------------------------------------------------------------------
! Overview     :
!  - Read TEprop.inp       --> Main input file
!  - Read band.eig file    --> Energy dispersion and k points
!  - Read linewidth.elself --> Imaginary part of self energy
!  - Calculate group velocity from energy dispersion
!  - Calculate Seebeck coefficient (S)
!  - Calculate electrical conductivity (\sigma)
!  - Calculate power factor (PF)
!  - Data files of S, \sigma, and PF can be as a function of E_F 
!    or as a function of carrier density
!*******************************************************************************
!*******************************************************************************
program TEprop2D
!-------------------------------------------------------------------------------

  implicit none

  character(80)        :: comment, finp0, finp1, finp2
  character(80)        :: fout2, fout3, fout4, fout5
  character(80)        :: fout6, fout7, fout8, fout9
  character(80)        :: fout51, fout61, fout71, fout81, fout91
  integer              :: iks, ikx, iky, nkx, nky, nks, nkpoint, iTempr
  integer              :: iline, nline, ibnd, ibnd0, nbnd, nbndc, nbndvb, nbndcb
  real(8), parameter   :: hpl   = 4.13566766225 ! [1e-15 eV.s] Planck's const.
  real(8), parameter   :: rkB   = 8.617330350D-5 ! Boltzmann constant [eV/K]
  real(8), parameter   :: pi    = 3.14159265358979323846D0
  real(8), parameter   :: zero  = 0.D0
  real(8), allocatable :: Ek(:,:), ImS(:,:), xk(:), yk(:), zk(:)
  real(8), allocatable :: Enk(:,:,:), ImSig(:,:,:), rk(:,:,:)
  real(8), allocatable :: dEdk(:,:,:,:), vg(:,:,:,:), tau(:,:,:)
  real(8)              :: alatt, thick, dummy, vpref, Vunit, hbar, kappa
  real(8)              :: rkBT, En, vx, vy, Efermi0, Efermif, Efermi
  real(8)              :: Efshift, dEfermi, fermiderivs, f0fermi, fermidist
  real(8)              :: sumL0x, sumL0y, sumL1x, sumL1y, sumL2x, sumL2y
  real(8)              :: tauel, Sxx, Syy, sigmax, sigmay, PFx, PFy, ZTx, ZTy
  real(8)              :: kappaex, kappaey, sumfermi, density
  real(8)              :: f0fermi0, sumfermi0, density0

! reading parameter file
  finp0 = "TEprop.inp"
  open(1,file=finp0,status="old")
  read(1,*) finp1
  read(1,*) finp2
  read(1,*) iTempr
  read(1,*) alatt
  read(1,*) thick
  read(1,*) nkx
  read(1,*) nky
  read(1,*) nbnd
  read(1,*) nbndvb
  read(1,*) nbndcb
  read(1,*) nbndc
  read(1,*) Efermi0, Efermif
  read(1,*) kappa
  nks = nkx * nky
  close(1)

  ! open energy band file
  open(10,file=trim(finp1),status="old")
  
  ! initialize index for k points
  iks = 1 

  ! number of lines to read in the energy dispersion file
  nline = 2 * nks

  ! read energy band data sets
  read(10,*) comment
  allocate(xk(nks))
  allocate(yk(nks))
  allocate(zk(nks))
  do iline = 2, nline + 1
     if (mod(iline,2) == 0) then
        ! store k-point coordinates
        read(10,*) xk(iks), yk(iks), zk(iks)
     else
        ! store energy band values
        iks = iks - 1
        read(10,*) dummy !Ek(iks, 1:nbnd)
     end if
     iks = iks + 1
  end do
  close(10)
 
  ! open self energy file
  open(20,file=trim(finp2),status="old")

  ! open output for scattering rate
  fout2 = "scatter-"//trim(finp1)//".dat"  
  open(21,file=fout2)
  
  ! re-initialize index for k points
  iks = 1

  ! number of lines to read in the self energy file
  nline = nbndc * nks

  ! reduced Planck's constant
  hbar  = hpl/(2*pi)

  ! read self energy data sets
  allocate(Ek(nks,nbnd))
  allocate(ImS(nks,nbnd))
  Ek  = zero
  ImS = zero
  read(20,*) comment
  read(20,*) comment
  do iline = 3, nline + 2
     read(20,*) iks, ibnd, Ek(iks,ibnd), ImS(iks,ibnd)
     write(21,"(2F14.6)") Ek(iks,ibnd), 2.D0*ImS(iks,ibnd)/hbar
     if (iline==3) ibnd0 = ibnd
  end do
  close(20); close(21);

  ! convert data from (nks) into (nkx) x (nky) k-point mesh
  allocate(Enk(nbnd,nkx,nky))
  allocate(ImSig(nbnd,nkx,nky))
  allocate(rk(nkx,nky,2))
  Enk = zero
  ImSig = zero
  iks = 1
  do ikx = 1, nkx
     do iky = 1, nky
        rk(ikx,iky,1) = xk(iks)
        rk(ikx,iky,2) = yk(iks)
        do ibnd = 1, nbnd ! the bands to be considered
           Enk(ibnd,ikx,iky)   = Ek(iks,ibnd)
           ImSig(ibnd,ikx,iky) = ImS(iks,ibnd)
        end do
        iks = iks + 1
     end do
  end do  

  ! hbar x group velocity = derivative of energy with respect to k
  ! calculate energy derivative firstly
  allocate(vg(nbnd,nkx,nky,2))
  allocate(dEdk(nbnd,nkx,nky,2))
  vg = zero
  dEdk = zero
  do ibnd = 1, nbnd
     do ikx = 2, nkx - 1
        do iky = 2, nky - 1
           dEdk(ibnd,ikx,iky,1) = (Enk(ibnd,ikx+1,iky)-Enk(ibnd,ikx-1,iky))&
                & / (rk(ikx+1,iky,1)-rk(ikx-1,iky,1))
           dEdk(ibnd,ikx,iky,2) = (Enk(ibnd,ikx,iky+1)-Enk(ibnd,ikx,iky-1))&
                & / (rk(ikx,iky+1,2)-rk(ikx,iky-1,2))
           vg(ibnd,ikx,iky,1)   =  dEdk(ibnd,ikx,iky,1)
           vg(ibnd,ikx,iky,2)   =  dEdk(ibnd,ikx,iky,2)
        end do
     end do
  end do

  ! calculate group velocity with units of [10^5 m/s]
  vpref = alatt / (2*pi*hbar) ! velocity prefactor
  vg    = vpref * vg
  ! write velocity output
  fout3 = "velocity-"//trim(finp1)//".dat"
  open(30,file=fout3)
  do ibnd = 1, nbnd
     do ikx = 2, nkx - 1
        do iky = 2, nky - 1
           write(30,"(I6,4F14.6)") ibnd, rk(ikx,iky,1), rk(ikx,iky,2),&
                & vg(ibnd,ikx,iky,1), vg(ibnd,ikx,iky,2)
        end do
     end do
  end do
  close(30)

  ! calculate relaxation time in units of [10^-12 s] (picosecond)
  ! \tau = \hbar / 2 Im (\Sigma)
  allocate(tau(nbnd,nkx,nky))
  tau = zero ! initialize value for all array components
  do ibnd = 1, nbnd
     do ikx = 1, nkx
        do iky = 1, nky
           tau(ibnd,ikx,iky) = 1.D3*hbar / (2*ImSig(ibnd,ikx,iky))
           ! note that 1.D3 * 1D-15 = 1D-12 s (picosecond)
           ! 1D-15 is from Planck's constant hpl in eV (c.f. hbar = hpl/2pi)
        end do
     end do
  end do

  ! Thermal energy (eV)
  rkBT   = rkB * dble(iTempr)

  ! Prepare output files
  fout4  = "density-"//trim(finp1)//".dat"
  ! additional label "1" means as a function of carrier density
  ! otherwise it's a function of Fermi energy (eV)
  fout5  = "Seebeck-"//trim(finp1)//".dat"
  fout51 = "Seebeck1-"//trim(finp1)//".dat"
  fout6  = "elsigma-"//trim(finp1)//".dat"
  fout61 = "elsigma1-"//trim(finp1)//".dat"
  fout7  = "PF-"//trim(finp1)//".dat"
  fout71 = "PF1-"//trim(finp1)//".dat"
  fout8  = "kappael-"//trim(finp1)//".dat"
  fout81 = "kappael1-"//trim(finp1)//".dat"
  fout9  = "ZT-"//trim(finp1)//".dat"
  fout91 = "ZT1-"//trim(finp1)//".dat"

  open(40,file=fout4)
  open(50,file=fout5)
  open(51,file=fout51)
  open(60,file=fout6)
  open(61,file=fout61)
  open(70,file=fout7)
  open(71,file=fout71)
  open(80,file=fout8)
  open(81,file=fout81)
  open(90,file=fout9)
  open(91,file=fout91)

  write(40,*) "#  Efermi (eV)    Carrier density in 10^21 /cm^3"
  write(50,*) "#  Efermi (eV)    Seebeck coefficient (mV/K) x and y"
  write(51,*) "#  density (10^21 cm-3) S (mV/K) x and y"
  write(60,*) "#  Efermi (eV)    El. conductivity (10^8/ohm.m) x and y"
  write(61,*) "#  density (10^21 cm-3) sigma (10^8/ohm.m) x and y"
  write(70,*) "#  Efermi (eV)    Power factor (W/K^2.m) x and y"
  write(71,*) "#  density (10^21 cm-3) PF (W/K^2.m) x and y"
  write(80,*) "#  Efermi (eV)    El. Therm. Cond. (W/(m.K)) x and y"
  write(81,*) "#  density (10^21 cm-3) kappael (W/(m.K)) x and y"
  write(90,*) "#  Efermi (eV)    ZT at x and y direction"
  write(91,*) "#  density (10^21 cm-3) ZT at x and y "

  Efermi  = Efermi0    ! initialize Fermi energy for data range [-2,2] eV
  Efshift = abs((minval(Ek(:,nbndcb))+maxval(Ek(:,nbndvb)))/2.D0)

  ! for carrier density reference/normalization:
  sumfermi0 = zero
  do ibnd = ibnd0, nbnd ! the bands to be considered (ibnd0 is not always 1)
     do ikx = 2, nkx - 1
        do iky = 2, nky - 1
           f0fermi0  = fermidist((Enk(ibnd,ikx,iky)+Efshift),zero,rkBT)
           sumfermi0 = sumfermi0 + f0fermi0
        end do
     end do
  end do

  nkpoint = (nkx-2)*(nky-2)
  Vunit   = (sqrt(3.D0)/2.D0)*((alatt*1.D-1)**2)*(thick*1.D-1) ! in [nm^3]

  density0 = (2.D0 / (dble(nkpoint)*Vunit) ) * sumfermi0

  ! main part for thermoelectric properties

  do while(Efermi < (Efermif + 0.001))

     ! initialize carrier density
     sumfermi = zero

     ! initialize L0
     sumL0x   = zero
     sumL0y   = zero

     ! initialize L1
     sumL1x   = zero
     sumL1y   = zero

     ! initialize L2
     sumL2x   = zero
     sumL2y   = zero

     ! main loop and summation
     do ibnd = ibnd0, nbnd ! the bands to be considered (ibnd0 is not always 1)
        do ikx = 2, nkx - 1
           do iky = 2, nky - 1

              ! allocate the main variables
              En      = Enk(ibnd,ikx,iky) + Efshift
              vx      = vg(ibnd,ikx,iky,1)
              vy      = vg(ibnd,ikx,iky,2)
              tauel   = tau(ibnd,ikx,iky)
              dEfermi = fermiderivs(En,Efermi,rkBT)
              f0fermi = fermidist(En,Efermi,rkBT)

              ! calculate carrier density
              sumfermi  = sumfermi + f0fermi

              ! x component of L integrals
              sumL0x  = sumL0x + ((vx**2)*tauel*dEfermi)
              sumL1x  = sumL1x + ((En-Efermi)*(vx**2)*tauel*dEfermi)
              sumL2x  = sumL2x + (((En-Efermi)**2)*(vx**2)*tauel*dEfermi)
              
              ! y component of L integrals
              sumL0y  = sumL0y + ((vy**2)*tauel*dEfermi)
              sumL1y  = sumL1y + ((En-Efermi)*(vy**2)*tauel*dEfermi)
              sumL2y  = sumL2y   + (((En-Efermi)**2)*(vy**2)*tauel*dEfermi)

           end do
        end do
     end do

     ! x component
     Sxx     = - (1.D3/dble(iTempr))*(sumL1x/sumL0x) ! S in [mili V/K]
     sigmax  = - 2.D0*(1.60217662D0/(dble(nkpoint)*Vunit)) &
                                   * (sumL0x)        ! sigma in [10^5 / ohm.m]
     PFx     = (Sxx**2) * sigmax                     ! PF in [10^-1 W/mK^2] 
     kappaex = - 2.D0*(1.60217662D0/(dble(nkpoint)*Vunit*dble(iTempr))) &
                 * ((sumL2x)-(sumL1x**2/sumL0x))     ! kappa in [10^6 W/(m.K)]
     ZTx     = (PFx*1D-1) * dble(iTempr) / (kappa + (kappaex*1.D6))

     ! y component
     Syy     = - (1.D3/dble(iTempr))*(sumL1y/sumL0y) ! S in [mili V/K]
     sigmay  = - 2.D0*(1.60217662D0/(dble(nkpoint)*Vunit)) &
                                   * (sumL0y)        ! sigma in [10^5 / ohm.m]
     PFy     = (Syy**2) * sigmay                     ! PF in [10^-1 W/mK^2] 
     kappaey = - 2.D0*(1.60217662D0/(dble(nkpoint)*Vunit*dble(iTempr))) &
                 * ((sumL2y)-(sumL1y**2/sumL0y))     ! kappa in [10^6 W/(m.K)]
     ZTy     = (PFy*1D-1) * dble(iTempr) / (kappa + (kappaey*1.D6))

     ! total carrier density 
     density = ( ((2.D0 / (dble(nkpoint)*Vunit) ) * sumfermi) - density0 )

     ! write density in the output file
     write(40,"(3F14.6)") Efermi, density ! density in units of [10^21 cm^-3]
     
     ! write Seebeck coefficient in the output file
     write(50,"(3F14.6)") Efermi, Sxx, Syy ! in mili V / K

     ! write electrical conductivity in the output file
     write(60,"(3F14.6)") Efermi, sigmax*1D-3, sigmay*1D-3 ! in [10^8 / ohm.m]
     
     ! write power factor in the output file
     write(70,"(3F14.6)") Efermi, PFx*1D-1, PFy*1D-1 ! in [W/K^2 m]

     ! write electronic thermal conductivity in the output file
     write(80,"(3F14.6)") Efermi, kappaex*1D6, kappaey*1D6 ! in [W/(mK)]

     ! write ZT in the output file
     write(90,"(3F14.6)") Efermi, ZTx, ZTy ! dimensionless

     ! write again all data as a function of density
     if (abs(density) > 1D-6) then
        write(51,"(3F14.6)") density, Sxx, Syy 
        write(61,"(3F14.6)") density, sigmax*1D-3, sigmay*1D-3
        write(71,"(3F14.6)") density, PFx*1D-1, PFy*1D-1
        write(81,"(3F14.6)") density, kappaex*1D6, kappaey*1D6
        write(91,"(3F14.6)") density, ZTx, ZTy
     end if

     ! change Fermi energy
     Efermi  = Efermi + 0.01D0
  end do

  close(40); close(50); close(60); close(70); close(80)
  close(61); close(71); close(81); close(81);

  deallocate(xk);  deallocate(yk);  deallocate(zk); deallocate(rk)
  deallocate(vg);  deallocate(tau); deallocate(dEdk)
  deallocate(Ek);  deallocate(Enk)
  deallocate(ImS); deallocate(ImSig)
  
end program TEprop2D
!*******************************************************************************
!*******************************************************************************
real(8) function fermiderivs(E,Ef,rkT)
!-------------------------------------------------------------------------------
  
  implicit none
  real(8), intent(in) :: E, Ef, rkT
  
  ! calculate derivative of Fermi energy
  fermiderivs = -1.D0 / (4*rkT * (cosh((E-Ef)/(2*rkT)))**2 )
  
end function fermiderivs
!*******************************************************************************
!*******************************************************************************
real(8) function fermidist(E,Ef,rkT)
!-------------------------------------------------------------------------------
  
  implicit none
  real(8), intent(in) :: E, Ef, rkT
  
  ! calculate Fermi distribution
  fermidist = 1.D0 / (exp((E-Ef)/rkT) + 1.D0)
  
end function fermidist
!*******************************************************************************
!*******************************************************************************
