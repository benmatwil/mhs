module fftw3

  use iso_c_binding

  implicit none

  include 'fftw3.f03'

end module

module harmonics

  use iso_fortran_env, np => real64
  use fftw3

  implicit none

  real(np), parameter :: pi = acos(-1.0_np), fourpi = 4*pi
  
  integer :: lmax
  real(np) :: rmax

  real(np), parameter :: alpha = 1, d = 0

  ! grid coordinates
  integer :: nrad, ntheta, nphi
  real(np), dimension(:), allocatable :: rads, thetas, phis
  character(2) :: nsplitstr

  ! coefficients required for calculating harmonics
  real(np), dimension(:,:), allocatable :: coeff1, coeff2, dfact
  
  ! all different global arrays for storing harmonics
  real(np), dimension(:,:,:), allocatable :: plm
  real(np), dimension(:,:,:), allocatable :: qlm, dqlm, qlm_sin
  complex(np), dimension(:,:), allocatable :: blm0
  complex(np), dimension(:,:,:), allocatable :: blm, alm

  ! fft variables
  type(c_ptr) :: plan
  complex(np), dimension(:), allocatable :: fft_in, fft_out

  contains

  ! ################################################################################

  subroutine calc_grids()

    ! Sets up grids in r, theta, phi for final field

    implicit none

    integer, parameter :: nsplit = 4
    integer :: ip, ir, it, isplit

#if fft
    if (.true.) then
      ntheta = (lmax+1)+1
      nphi = 2*(lmax+1)+1
      nrad = 0
      do while (exp(pi*nrad/(lmax+1)/2) < rmax)
        nrad = nrad + 1
      enddo
      ! nrad = ceiling(2*(lmax+1)*log(rmax)/pi)

      allocate(rads(nrad), thetas(ntheta), phis(nphi))

      rads = exp((pi/(lmax+1)/2)*real([(ir, ir=0,nrad-1)], np))
      print*, rads(1), rads(nrad)
      rads = (rmax - 1)*(rads - 1)/(rads(nrad) - 1) + 1 ! re-scale to 1 to rmax
      ! rads = 1.5_np*(real([(ir, ir=0,nrad-1)], np))**3/(nrad-1)**3 + 1 ! cubic grid

      thetas = pi*real([(it, it=0,ntheta-1)], np)/(ntheta-1)
      phis = 2*pi*real([(ip, ip=0,nphi-1)], np)/(nphi-1)
    else
      ntheta = (lmax+1)+1
      ntheta = nsplit*ntheta-(nsplit-1)
      nphi = 2*(lmax+1)+1
      nrad = 0
      do while (exp(pi*nrad/(lmax+1)/2) < rmax)
        nrad = nrad + 1
      enddo
      ! nrad = ceiling(2*(lmax+1)*log(rmax)/pi)
      ! nrad = nsplit*nrad-(nsplit-1)

      allocate(rads(nrad), thetas(ntheta), phis(nphi))

      ! do ir = 0, nrad-1, nsplit
      !   rads(ir+1) = exp((pi/(lmax+1)/2)*(ir/nsplit))
      ! enddo
      ! print*, rads(1), rads(nrad)
      ! do ir = 2, nrad, nsplit
      !   do isplit = 0, nsplit-2
      !     rads(ir+isplit) = ((nsplit-isplit-1)*rads(ir-1) + (isplit+1)*rads(ir+nsplit-1))/nsplit
      !   enddo
      ! enddo

      rads = exp((pi/(lmax+1)/2)*real([(ir, ir=0,nrad-1)], np))
      print*, rads(1), rads(nrad)
      rads = (rmax - 1)*(rads - 1)/(rads(nrad) - 1) + 1 ! re-scale to 1 to rmax
      ! rads = 1.5_np*(real([(ir, ir=0,nrad-1)], np))**3/(nrad-1)**3 + 1 ! cubic grid

      thetas = pi*real([(it, it=0,ntheta-1)], np)/(ntheta-1)
      phis = 2*pi*real([(ip, ip=0,nphi-1)], np)/(nphi-1)
    endif
  
#elif analytic

    ! ntheta = (lmax+1)+1
    ! nphi = 2*(lmax+1)+1
    ! nrad = 0
    ! do while (exp(pi*nrad/(lmax+1)/2) < rmax)
    !   nrad = nrad + 1
    ! enddo
    ! ! nrad = ceiling(2*(lmax+1)*log(rmax)/pi)

    ! allocate(rads(nrad), thetas(ntheta), phis(nphi))

    ! rads = exp((pi/(lmax+1)/2)*real([(ir, ir=0,nrad-1)], np))
    ! print*, rads(1), rads(nrad)
    ! rads = (rmax - 1)*(rads - 1)/(rads(nrad) - 1) + 1 ! re-scale to 1 to rmax
    ! ! rads = 1.5_np*(real([(ir, ir=0,nrad-1)], np))**3/(nrad-1)**3 + 1 ! cubic grid

    ! thetas = pi*real([(it, it=0,ntheta-1)], np)/(ntheta-1)
    ! phis = 2*pi*real([(ip, ip=0,nphi-1)], np)/(nphi-1)

    nrad = 0
    do while (exp(pi*nrad/(lmax+1)/2) < rmax)
      nrad = nrad + 1
    enddo
    ntheta = (lmax+1)+1
    nphi = 2*(lmax+1)+1

    nrad = nsplit*nrad-(nsplit-1)
    ntheta = nsplit*ntheta-(nsplit-1)
    nphi = nsplit*nphi-(nsplit-1)

    allocate(rads(nrad), thetas(ntheta), phis(nphi))

    do ir = 0, nrad-1, nsplit
      rads(ir+1) = exp((pi/(lmax+1)/2)*(ir/nsplit))
    enddo
    print*, rads(1), rads(nrad)
    do ir = 2, nrad, nsplit
      do isplit = 0, nsplit-2
        rads(ir+isplit) = ((nsplit-isplit-1)*rads(ir-1) + (isplit+1)*rads(ir+nsplit-1))/nsplit
      enddo
    enddo

    rads = (rmax - 1)*(rads - 1)/(rads(nrad) - 1) + 1

    thetas = pi*real([(it, it=0,ntheta-1)], np)/(ntheta-1)
    phis = 2*pi*real([(ip, ip=0,nphi-1)], np)/(nphi-1)

    nrad = size(rads,1)
    ntheta = size(thetas,1)
    nphi = size(phis,1)

    write(nsplitstr,'(I2.2)') nsplit

#endif

  end

  ! ################################################################################

  subroutine calc_coeffs()

    ! calculates all the coefficients required for both plms and qlms

    implicit none

    integer :: il, im

    allocate(coeff1(0:lmax,0:lmax), coeff2(0:lmax,0:lmax), dfact(0:lmax,0:lmax))
    coeff1 = 0
    coeff2 = 0
    dfact = 0

    !$omp parallel do private(il)
    do im = 0, lmax
      do il = im, lmax
        coeff1(il, im) = sqrt(real(4*il**2 - 1, np))
        coeff2(il, im) = -sqrt(real(2*il + 1, np)*real((il-1)**2 - im**2, np)/(2*il - 3))
        dfact(il, im) = sqrt(real(il**2 - im**2, np))
      enddo
    enddo

  end

  ! ################################################################################

  subroutine calc_plms(xs)

    ! Calculate the Qlms (normalised Plms) for a grid xs
    ! in the range -1 to 1 (cos(theta) = x)

    implicit none

    integer :: il, im, it, ntheta
    real(np) :: xs(:)
    real(np) :: x, y

    if (.not. allocated(coeff1)) call calc_coeffs()

    ntheta = size(xs,1)

    allocate(plm(0:lmax,0:lmax,ntheta))
    plm = 0

    !$omp parallel do private(x, y, im, il)
    do it = 1, ntheta
      x = xs(it)
      y = sqrt(1 - x**2)
      do im = 0, lmax
        if (im == 0) then
          plm(0, im, it) = sqrt(1.0_np/fourpi)
          plm(1, im, it) = sqrt(3.0_np/fourpi)*x
        else
          plm(im, im, it) = -sqrt(real(2*im + 1, np)/(2*im))*y*plm(im-1, im-1, it)
          if (im+1 <= lmax) plm(im+1, im, it) = sqrt(real(2*im + 3, np))*x*plm(im, im, it)
        endif
        do il = im+2, lmax
          if (im < lmax-1) plm(il, im, it) = (x*plm(il-1, im, it)*coeff1(il, im) + &
            plm(il-2, im, it)*coeff2(il, im))/dfact(il, im)
        enddo
      enddo
    enddo

  end

  ! ################################################################################

  subroutine calc_blm_rsun(synmap, nfilter)

    implicit none

    integer :: ia, il, im
    integer, optional :: nfilter
    integer :: ilat, ilon, nlon, nlat
    real(np) :: dt
    real(np), dimension(0:lmax) :: filter
    integer, dimension(0:lmax) :: ls
    complex(np), dimension(:,:), allocatable :: br_blm
    real(np), dimension(:,:), allocatable :: synmap
  
    nlon = size(synmap,1)
    nlat = size(synmap,2)
    
    ia = 10
    do while (exp(-0.25_np*(pi/ia)**2*lmax*(lmax+1)) < 0.1_np)
      ia = ia + 10
    enddo
    ia = ia - 10

    if (present(nfilter)) then
      if (nfilter /= 0) ia = nfilter
    endif
    
    ls = [(il, il=0,lmax)]
    filter = exp(-0.25_np*ls*(ls+1)*(pi/(ia))**2)

    ! not technically calculating the Blms at the correct points... should this be corrected?
    ! may require not using a fft - perhaps quite slow
    allocate(br_blm(0:lmax,nlat))

    allocate(fft_in(nlon), fft_out(nlon))
    plan = fftw_plan_dft_1d(nlon, fft_in,fft_out, fftw_forward,fftw_measure)
    
    !$omp parallel do private(fft_in, fft_out)
    do ilat = 1, nlat
      fft_in = synmap(:,ilat)
      call fftw_execute_dft(plan, fft_in, fft_out)
      br_blm(:,ilat) = fft_out(1:lmax+1)/nlon
    enddo

    call fftw_destroy_plan(plan)
    deallocate(fft_in, fft_out)

    dt = fourpi/nlat

    allocate(blm0(0:lmax, 0:lmax))
    blm0 = 0
    do im = 0, lmax
      do il = im, lmax
        blm0(il, im) = dt*sum(br_blm(im, :)*plm(il, im, :))*filter(il)/il/(il+1)
      enddo
    enddo
    blm0(0,0) = 0 ! otherwise it's infinite (division by zero)

  end subroutine

  ! ################################################################################

  subroutine calc_qlms(theta)

    ! Calculate the Qlms (normalised Plms) for a grid linear in theta
    ! for the range 0 to 2pi

    implicit none

    integer :: il, im, it, ntheta
    real(np) :: theta(:)
    real(np) :: x, y
    real(np) :: c1, c2

    if (.not. allocated(coeff1)) call calc_coeffs()
    
    ! ntheta = 4*(lmax+1)+1
    ntheta = size(theta,1)
    
    allocate(qlm(0:lmax, 0:lmax, ntheta), dqlm(0:lmax, 0:lmax, ntheta))
    qlm = 0
    dqlm = 0
    !$omp parallel do private(x, y, im, il, c1, c2)
    do it = 1, ntheta
      x = cos(theta(it))
      y = sin(theta(it))
      do im = 0, lmax
        if (im == 0) then
          qlm(0, im, it) = sqrt(1.0_np/fourpi)
          dqlm(0, im, it) = 0
          qlm(1, im, it) = sqrt(3.0_np/fourpi)*x
          dqlm(1, im, it) = -sqrt(3.0_np/fourpi)*y
        else
          c1 = -sqrt(real(2*im + 1, np)/(2*im))
          qlm(im, im, it) = c1*y*qlm(im-1, im-1, it)
          dqlm(im, im, it) = c1*(x*qlm(im-1, im-1, it) + y*dqlm(im-1, im-1, it))
          if (im+1 <= lmax) then
            c2 = sqrt(real(2*im + 3, np))
            qlm(im+1, im, it) = c2*x*qlm(im, im, it)
            dqlm(im+1, im, it) = c2*(x*dqlm(im, im, it) - y*qlm(im, im, it))
          endif
        endif
        do il = im+2, lmax
          if (im < lmax-1) then
            qlm(il, im, it) = (x*qlm(il-1, im, it)*coeff1(il, im) + &
              qlm(il-2, im, it)*coeff2(il, im))/dfact(il, im)
            dqlm(il, im, it) = ((x*dqlm(il-1, im, it)-y*qlm(il-1, im, it))*coeff1(il, im) + &
              dqlm(il-2, im, it)*coeff2(il, im))/dfact(il, im)
          endif
        enddo
      enddo
    enddo

    ! add sin dependence for B_phi component
    qlm_sin = qlm
    qlm_sin = 0

    if (ntheta > 1) then
      ! assuming the set of theta end at the poles
      ! !$omp parallel do
      do it = 2, ntheta-1
        qlm_sin(:,1:,it) = qlm(:,1:,it)/sin(theta(it))
      enddo
      qlm_sin(:,1,1) = qlm_sin(:,1,2)
      qlm_sin(:,1,ntheta) = qlm_sin(:,1,ntheta-1)
    else
      qlm_sin(:,1:,1) = qlm(:,1:,1)/sin(theta(1))
    endif

  end

  ! ################################################################################

  subroutine calc_rdep(rads)

    ! extrapolate the blms for each radii

    implicit none

    integer :: ir, il, im, nrad
    complex(np), parameter :: s = cmplx(0,-1,np)
    complex(np) :: fact1
    real(np), dimension(:) :: rads
    real(np), dimension(:), allocatable :: rads1
    real(np) :: rad_fact
    ! integer, dimension(0:lmax) :: ls
    ! real(np), dimension(0:lmax) :: r_rsun, r_rmax, div_fact
    real(np), dimension(:,:), allocatable :: bess1, bess2, dbess1, dbess2
    real(np), dimension(:,:), allocatable :: ibess1
    real(np) :: bess10
    ! real(np) :: rj,ry,rjp,ryp
    complex(np), dimension(0:lmax) :: hfact, bess1fact, bess2fact, div_fact, bess_fact, dbess_fact
    complex(np), dimension(:,:), allocatable :: rdep_blm, rdep_alm

    nrad = size(rads,1)
    allocate(bess1(0:lmax, nrad), bess2(0:lmax, nrad))
    allocate(dbess1(0:lmax, nrad), dbess2(0:lmax, nrad))
    allocate(rdep_blm(0:lmax, nrad), rdep_alm(0:lmax, nrad))

    rads1 = alpha*(rads + d)
    
    ! calculate all required bessel functions
    ! bess1(0,:) = sqrt(2/pi/rads1)*sin(rads1)
    ! bess1(1,:) = sqrt(2/pi/rads1)*(sin(rads1)/rads1 - cos(rads1))

    allocate(ibess1(nrad,3))
    ibess1(:,1) = 0
    ibess1(:,2) = 1
    do il = lmax, 0, -1
      ibess1(:,3) = (2*il + 3)*ibess1(:,2)/rads1 - ibess1(:,1)
      ibess1(:,1:2) = ibess1(:,2:3)
      bess1(il,:) = ibess1(:,3)
    enddo

    do ir = 1, nrad
      bess10 = sqrt(2/pi/rads1(ir))*sin(rads1(ir))
      bess1(:,ir) = bess1(:,ir)*bess10/bess1(0,ir)
    enddo
    
    bess2(0,:) = -sqrt(2/pi/rads1)*cos(rads1)
    bess2(1,:) = -sqrt(2/pi/rads1)*(cos(rads1)/rads1 + sin(rads1))
    do il = 2, lmax
      bess2(il,:) = ((2*il - 1)*bess2(il-1,:) - rads1*bess2(il-2,:))/rads1
    enddo
    ! wronskian if required    
    ! bess1(il,:) = (2/pi/rads1 + bess2(il,:)*bess1(il-1,:))/bess2(il-1,:)

    ! ir = 1
    ! print*, rads1(ir)
    ! do il = 0, lmax
    !   print*, bess1(il,ir), bess2(il,ir)
    ! enddo

    ! and all the derivatives
    dbess1(0,:) = alpha*sqrt(2/pi/rads1)*(cos(rads1) - sin(rads1)/rads1/2)
    dbess2(0,:) = alpha*sqrt(2/pi/rads1)*(sin(rads1) + cos(rads1)/rads1/2)

    do il = 1, lmax
      dbess1(il,:) = alpha*(bess1(il-1,:) - (il + 0.5_np)*bess1(il,:)/rads1)
      dbess2(il,:) = alpha*(bess2(il-1,:) - (il + 0.5_np)*bess2(il,:)/rads1)
    enddo

    ! print*, dbess1(:,nrad)

    fact1 = 1/(rads(nrad) + d)/2 - s*alpha
    ! hfact = (dbess2(:,nrad) + fact1*bess2(:,nrad))/(dbess1(:,nrad) + fact1*bess1(:,nrad))

    ! div_fact = bess2(:, 1) - hfact*bess1(:, 1)

    bess1fact = dbess1(:,nrad) + fact1*bess1(:,nrad)
    bess2fact = dbess2(:,nrad) + fact1*bess2(:,nrad)

    div_fact = bess2(:,1)*bess1fact - bess1(:,1)*bess2fact
    ! print*, bess1(:,1)
    ! print*, '-------------------------------------------'
    ! print*, bess2(:,1)
    ! print*, '-------------------------------------------'
    ! print*, bess1fact
    ! print*, '-------------------------------------------'
    ! print*, bess2fact
    ! print*, '-------------------------------------------'
    ! print*, div_fact

    ! !$omp parallel do private(r_rsun, r_rmax)
    do ir = 1, nrad
      rad_fact = 1/rads(ir)*sqrt((rads(ir) + d)/(1 + d))
      bess_fact = bess2(:, ir)*bess1fact - bess1(:, ir)*bess2fact
      dbess_fact = dbess2(:, ir)*bess1fact - dbess1(:, ir)*bess2fact

      if (ir == nrad) then
        ! print*, bess1(:,ir)
        ! print*, bess2(:,ir)
        ! print*, dbess1(:,ir)
        ! print*, dbess2(:,ir)
        ! print*, bess_fact
        ! print*, '-------------------------------------------'
        ! print*, dbess_fact
        ! print*, '-------------------------------------------'
        ! print*, 1/div_fact
        ! print*, '-------------------------------------------'
        ! print*, rad_fact*(0.5_np*bess_fact/(rads(ir) + d) + dbess_fact)/div_fact
      endif

      ! if (ir == 1) print*, bess_fact/div_fact
      rdep_alm(:, ir) = rad_fact*bess_fact/div_fact
      rdep_blm(:, ir) = rad_fact*(0.5_np*bess_fact/(rads(ir) + d) + dbess_fact)/div_fact
    enddo

    allocate(blm(0:lmax, 0:lmax, nrad))
    ! !$omp parallel do
    do ir = 1, nrad
      blm(:, :, ir) = blm0
    enddo
    !deallocate(blm0)

    alm = blm
    ! !$omp parallel do
    do im = 0, lmax
      blm(:, im, :) = blm(:, im, :)*rdep_blm(:, :)
      alm(:, im, :) = alm(:, im, :)*rdep_alm(:, :)
    enddo

  end

  ! ################################################################################

  subroutine calc_final_field(br, bt, bp)

    implicit none

    integer :: ir, ip, it, il, im, jm
    complex(np) :: mi, expphi
    real(np), dimension(:,:,:), allocatable :: br, bt, bp
    complex(np), dimension(:,:), allocatable :: bbr, bbt, bbp

    allocate(br(nphi, ntheta, nrad), bt(nphi, ntheta, nrad), bp(nphi, ntheta, nrad))

    !$omp parallel workshare
    br = 0
    bt = 0
    bp = 0
    !$omp end parallel workshare

#if analytic

    ! calculating grid analytically
    !$omp parallel do private(ir, it, ip, im, il, mi, jm, expphi)
    do ir = 1, nrad
      print*, ir
      do it = 1, ntheta
        do ip = 1, nphi
          do im = -lmax, lmax
            jm = abs(im)
            mi = cmplx(0.0_np,1.0_np,np)*jm
            expphi = exp(mi*phis(ip))
            if (im >= 0) then
              do il = lmax, abs(im), -1
                br(ip, it, ir) = br(ip, it, ir) + blm(il, jm, ir)*qlm(il, jm, it)*expphi
                bt(ip, it, ir) = bt(ip, it, ir) + alm(il, jm, ir)*dqlm(il, jm, it)*expphi
                bp(ip, it, ir) = bp(ip, it, ir) + mi*alm(il, jm, ir)*qlm_sin(il, jm, it)*expphi
              enddo
            else
              do il = lmax, abs(im), -1
                br(ip, it, ir) = br(ip, it, ir) + conjg(blm(il, jm, ir)*qlm(il, jm, it)*expphi)
                bt(ip, it, ir) = bt(ip, it, ir) + conjg(alm(il, jm, ir)*dqlm(il, jm, it)*expphi)
                bp(ip, it, ir) = bp(ip, it, ir) + conjg(mi*alm(il, jm, ir)*qlm_sin(il, jm, it)*expphi)
              enddo
            endif
          enddo
        enddo
      enddo
    enddo

#elif fft

    allocate(bbr(nphi-1, ntheta), bbt(nphi-1, ntheta), bbp(nphi-1, ntheta))
    
    allocate(fft_in(size(bbr,1)), fft_out(size(bbr,1)))
    plan = fftw_plan_dft_1d(size(fft_in,1), fft_in,fft_out, fftw_backward,fftw_measure)

    !$omp parallel do private(im, mi, il, bbr, bbt, bbp, jm, fft_in, fft_out)
    do ir = 1, nrad
      if (mod(ir, 20) == 0) print*, ir
      bbr = 0
      bbt = 0
      bbp = 0
      do im = lmax, 1, -1
        mi = cmplx(0,im,np)
        jm = im + 1
        do il = lmax, im, -1
          ! sum over l first in preparation for fft
          bbr(jm, :) = bbr(jm, :) + il*(il+1)*alm(il, im, ir)*qlm(il, im, :)/rads(ir)
          bbt(jm, :) = bbt(jm, :) + (mi*alpha*alm(il, im, ir)*qlm_sin(il, im, :) + &
            blm(il, im, ir)*dqlm(il, im, :))
          bbp(jm, :) = bbp(jm, :) + (-alpha*alm(il, im, ir)*dqlm(il, im, :) + &
            mi*blm(il, im, ir)*qlm_sin(il, im, :))
        enddo
        if (im > 0) then
          ! use the conjugation for negative m rather than direct calculation
          bbr(nphi-im, :) = conjg(bbr(jm, :))
          bbt(nphi-im, :) = conjg(bbt(jm, :))
          bbp(nphi-im, :) = conjg(bbp(jm, :))
        endif
      enddo
      ! fft for each constant theta to get field values
      do it = 1, ntheta
        fft_in = bbr(:,it)
        call fftw_execute_dft(plan, fft_in, fft_out)
        br(1:nphi-1,it,ir) = fft_out
        br(nphi,it,ir) = fft_out(1)

        fft_in = bbt(:,it)
        call fftw_execute_dft(plan, fft_in, fft_out)
        bt(1:nphi-1,it,ir) = fft_out
        bt(nphi,it,ir) = fft_out(1)

        fft_in = bbp(:,it)
        call fftw_execute_dft(plan, fft_in, fft_out)
        bp(1:nphi-1,it,ir) = fft_out
        bp(nphi,it,ir) = fft_out(1)
      enddo
    enddo

    call fftw_destroy_plan(plan)
    deallocate(fft_in, fft_out)

#endif

  end

  ! ################################################################################

  function calc_point(r, theta, phi)

    ! need to have calculated blm(rsun) first
    ! r and theta calculated by sending the correct info to other routines

    implicit none

    real(np) :: r, theta, phi, calc_point(3)
    integer :: il, im, jm
    complex(np) :: mi, expphi
    real(np) :: br, bt, bp

    call calc_qlms([theta])
    call calc_rdep([r])

    br = 0
    bt = 0
    bp = 0

    do im = -lmax, lmax
      jm = abs(im)
      mi = cmplx(0.0_np,1.0_np,np)*jm
      expphi = exp(mi*phi)
      if (im >= 0) then
        do il = lmax, abs(im), -1
          br = br + blm(il, jm, 1)*qlm(il, jm, 1)*expphi
          bt = bt + alm(il, jm, 1)*dqlm(il, jm, 1)*expphi
          bp = bp + mi*alm(il, jm, 1)*qlm_sin(il, jm, 1)*expphi
        enddo
      else
        do il = lmax, abs(im), -1
          br = br + conjg(blm(il, jm, 1)*qlm(il, jm, 1)*expphi)
          bt = bt + conjg(alm(il, jm, 1)*dqlm(il, jm, 1)*expphi)
          bp = bp + conjg(mi*alm(il, jm, 1)*qlm_sin(il, jm, 1)*expphi)
        enddo
      endif
    enddo

    deallocate(qlm, dqlm, qlm_sin, alm, blm)

    calc_point = [br, bt, bp]

  end

  ! subroutine bessjy(x,xnu,rj,ry,rjp,ryp)
    
  !   implicit none

  !   real(np), intent(in) :: x, xnu
  !   real(np), intent(out) :: rj, ry, rjp, ryp
  !   integer(int32), parameter :: maxit = 10000
  !   real(np), parameter :: xmin = 2.0_np, eps = 1.0e-16_np, fpmin = 1.0e-300_np
  !   integer(int32) :: i, isign, l, nl
  !   real(np) :: a, b, c, d, del, del1, e, f, fact, fact2, fact3, ff, gam, gam1, gam2,&
  !   gammi, gampl, h, p, pimu, pimu2, q, r, rjl, rjl1, rjmu, rjp1, rjpl, rjtemp, &
  !   ry1, rymu, rymup, rytemp, sum, sum1, w, x2, xi, xi2, xmu, xmu2
  !   complex(np) :: aa, bb, cc, dd, dl, pq

  !   nl=merge(int(xnu+0.5_np), max(0,int(xnu-x+1.5_np)), x < xmin)
  !   xmu=xnu-nl
  !   xmu2=xmu*xmu
  !   xi=1.0_np/x
  !   xi2=2.0_np*xi
  !   w=xi2/pi
  !   isign=1
  !   h=xnu*xi
  !   if (h < fpmin) h=fpmin
  !   b=xi2*xnu
  !   d=0.0
  !   c=h
  !   do i=1,maxit
  !     b=b+xi2
  !     d=b-d
  !     if (abs(d) < fpmin) d=fpmin
  !     c=b-1.0_np/c
  !     if (abs(c) < fpmin) c=fpmin
  !     d=1.0_np/d
  !     del=c*d
  !     h=del*h
  !     if (d < 0.0) isign=-isign
  !     if (abs(del-1.0_np) < eps) exit
  !   end do
  !   if (i > maxit) stop 'x too large in bessjy; try asymptotic expansion'
  !   rjl=isign*fpmin
  !   rjpl=h*rjl
  !   rjl1=rjl
  !   rjp1=rjpl
  !   fact=xnu*xi
  !   do l=nl,1,-1
  !     rjtemp=fact*rjl+rjpl
  !     fact=fact-xi
  !     rjpl=fact*rjtemp-rjl
  !     rjl=rjtemp
  !   end do
  !   if (rjl == 0.0) rjl=eps
  !   f=rjpl/rjl
  !   if (x < xmin) then
  !     x2=0.5_np*x
  !     pimu=pi*xmu
  !     if (abs(pimu) < eps) then
  !       fact=1.0
  !     else
  !       fact=pimu/sin(pimu)
  !     end if
  !     d=-log(x2)
  !     e=xmu*d
  !     if (abs(e) < eps) then
  !       fact2=1.0
  !     else
  !       fact2=sinh(e)/e
  !     end if
  !     call beschb(xmu,gam1,gam2,gampl,gammi)
  !     ff=2.0_np/pi*fact*(gam1*cosh(e)+gam2*fact2*d)
  !     e=exp(e)
  !     p=e/(gampl*pi)
  !     q=1.0_np/(e*pi*gammi)
  !     pimu2=0.5_np*pimu
  !     if (abs(pimu2) < eps) then
  !       fact3=1.0
  !     else
  !       fact3=sin(pimu2)/pimu2
  !     end if
  !     r=pi*pimu2*fact3*fact3
  !     c=1.0
  !     d=-x2*x2
  !     sum=ff+r*q
  !     sum1=p
  !     do i=1,maxit
  !       ff=(i*ff+p+q)/(i*i-xmu2)
  !       c=c*d/i
  !       p=p/(i-xmu)
  !       q=q/(i+xmu)
  !       del=c*(ff+r*q)
  !       sum=sum+del
  !       del1=c*p-i*del
  !       sum1=sum1+del1
  !       if (abs(del) < (1.0_np+abs(sum))*eps) exit
  !     end do
  !     if (i > maxit) stop 'bessy series failed to converge'
  !     rymu=-sum
  !     ry1=-sum1*xi2
  !     rymup=xmu*xi*rymu-ry1
  !     rjmu=w/(rymup-f*rymu)
  !   else
  !     a=0.25_np-xmu2
  !     pq=cmplx(-0.5_np*xi,1.0_np,kind=np)
  !     aa=cmplx(0.0_np,xi*a,kind=np)
  !     bb=cmplx(2.0_np*x,2.0_np,kind=np)
  !     cc=bb+aa/pq
  !     dd=1.0_np/bb
  !     pq=cc*dd*pq
  !     do i=2,maxit
  !       a=a+2*(i-1)
  !       bb=bb+cmplx(0.0_np,2.0_np,kind=np)
  !       dd=a*dd+bb
  !       if (absc(dd) < fpmin) dd=fpmin
  !       cc=bb+a/cc
  !       if (absc(cc) < fpmin) cc=fpmin
  !       dd=1.0_np/dd
  !       dl=cc*dd
  !       pq=pq*dl
  !       if (absc(dl-1.0_np) < eps) exit
  !     end do
  !     if (i > maxit) stop 'cf2 failed in bessjy'
  !     p=real(pq)
  !     q=aimag(pq)
  !     gam=(p-f)/q
  !     rjmu=sqrt(w/((p-f)*gam+q))
  !     rjmu=sign(rjmu,rjl)
  !     rymu=rjmu*gam
  !     rymup=rymu*(p+q/gam)
  !     ry1=xmu*xi*rymu-rymup
  !   end if
  !   fact=rjmu/rjl
  !   rj=rjl1*fact
  !   rjp=rjp1*fact
  !   do i=1,nl
  !     rytemp=(xmu+i)*xi2*ry1-rymu
  !     rymu=ry1
  !     ry1=rytemp
  !   end do
  !   ry=rymu
  !   ryp=xnu*xi*rymu-ry1
  
  ! end subroutine bessjy

  ! function absc(z)

  !   implicit none

  !   complex(np), intent(in) :: z
  !   real(np) :: absc

  !   absc=abs(real(z))+abs(aimag(z))

  ! end function absc

  ! subroutine beschb(x,gam1,gam2,gampl,gammi)
    
  !   implicit none
    
  !   real(np), intent(in) :: x
  !   real(np), intent(out) :: gam1,gam2,gampl,gammi
  !   integer(int32), parameter :: nuse1=7,nuse2=8
  !   !if converting to double precision, set nuse1 = 7, nuse2 = 8.
  !   real(np) :: xx
  !   real(np), dimension(7) :: c1=[-1.142022680371168_np,&
  !   6.5165112670737e-3_np,3.087090173086e-4_np,-3.4706269649e-6_np,&
  !   6.9437664e-9_np,3.67795e-11_np,-1.356e-13_np]
  !   real(np), dimension(8) :: c2=[1.843740587300905_np,&
  !   -7.68528408447867e-2_np,1.2719271366546e-3_np,&
  !   -4.9717367042e-6_np, -3.31261198e-8_np,2.423096e-10_np,&
  !   -1.702e-13_np,-1.49e-15_np]
    
  !   xx=8.0_np*x*x-1.0_np
  !   gam1=chebev(-1.0_np,1.0_np,c1(1:nuse1),xx)
  !   gam2=chebev(-1.0_np,1.0_np,c2(1:nuse2),xx)
  !   gampl=gam2-x*gam1
  !   gammi=gam2+x*gam1

  ! end subroutine beschb

  ! function chebev(a,b,c,x)
    
  !   implicit none
    
  !   real(np), intent(in) :: a,b,x
  !   real(np), dimension(:), intent(in) :: c
  !   real(np) :: chebev_s
  !   integer(int32) :: j,m
  !   real(np) :: d,dd,sv,y,y2
  !   real(np) :: chebev
    
  !   if ((x-a)*(x-b) > 0.0_np) stop 'x not in range in chebev_s'
  !   m=size(c)
  !   d=0.0_np
  !   dd=0.0_np
  !   y=(2.0_np*x-a-b)/(b-a)
  !   y2=2.0_np*y
  !   do j=m,2,-1
  !   sv=d
  !   d=y2*d-dd+c(j)
  !   dd=sv
  !   end do
  !   chebev_s=y*d-dd+0.5_np*c(1)

  ! end function chebev

end module