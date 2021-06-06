    module global

      implicit none

      integer(4),    parameter     ::                    &
          sp = kind(1.0),                                &
          dp = selected_real_kind(2*precision(1.0_sp)),  &
          qp = selected_real_kind(2*precision(1.0_dp))
      integer(4),    parameter     :: kr=qp
      integer(4)                   :: problem

    end module global

    program main

      implicit none

      !$ call omp_set_num_threads(4)

      call drive

    end program main

    subroutine drive

      use global
      !$ use omp_lib, only: omp_set_num_threads

      implicit none

      integer(4)                   :: sol
      integer(4)                   :: jmax
      integer(4),    parameter     :: n=16
      integer(4),    parameter     :: kmax=10000000
      real(kind=kr)                :: c
      real(kind=kr)                :: h
      real(kind=kr)                :: eps
      real(4)                      :: start 
      real(4)                      :: finish

      !$ call omp_set_num_threads(n/2)

      ! Run Problem 1 with coarse mesh
      call cpu_time(start)
      problem=1
      c=1.0_kr
      eps=1.0e-08
      do sol=0,2
        jmax=5
        h=1.0_kr
        call solve_slab(sol,c,n,kmax,jmax,h,eps)
      enddo
      call cpu_time(finish)
      write (0,'(a,f6.3,a)') ' Problem 1 run time ',(finish-start)/60.0,' (sec).'

      ! Run Problem 1 with fine mesh and DD (reference)
      call cpu_time(start)
      eps=1.0e-06
      jmax=10000
      h=0.0005_kr
      call solve_slab(2,c,n,kmax,jmax,h,eps)
      call cpu_time(finish)
      write (0,'(a,f10.3,a)') ' Reference 1 run time ',(finish-start)/60.0,' (sec).'

      ! Run Problem 1 with coarse mesh
      call cpu_time(start)
      problem=2
      c=1.0_kr
      eps=1.0e-08
      do sol=0,2
        jmax=20
        h=1.0_kr
        call solve_slab(sol,c,n,kmax,jmax,h,eps)
      enddo
      call cpu_time(finish)
      write (0,'(a,f6.3,a)') ' Problem 2 run time ',(finish-start)/60.0,' (sec).'

      ! Run Problem 2 with fine mesh and DD (reference)
      call cpu_time(start)
      eps=1.0e-06
      jmax=10000
      h=0.002_kr
      call solve_slab(2,c,n,kmax,jmax,h,eps)
      call cpu_time(finish)
      write (0,'(a,f10.3,a)') ' Reference 1 run time ',(finish-start)/60.0,' (sec).'

    end subroutine drive

    subroutine solve_slab(sol,c,n,kmax,jmax,h,eps)

      use global

      implicit none

      integer(4),    intent(in)    :: sol
      integer(4),    intent(in)    :: n
      integer(4),    intent(in)    :: kmax
      integer(4),    intent(in)    :: jmax
      real(kind=kr), intent(in)    :: c
      real(kind=kr), intent(in)    :: h
      real(kind=kr), intent(in)    :: eps

      integer(4)                   :: j
      integer(4)                   :: bc(2)
      real(kind=kr)                :: mu(n/2)
      real(kind=kr)                :: w (n/2)
      real(kind=kr), allocatable   :: phi (:)
      real(kind=kr), allocatable   :: phi_l(:)
      real(kind=kr), allocatable   :: phi_r(:)
      real(kind=kr)                :: q
      real(kind=kr)                :: sigt
      real(kind=kr)                :: sigs
      character(1)                 :: prob
      character(2)                 :: solopt(0:2)=(/'LD','LC','DD'/)
      character(6)                 :: cells
      character(20)                :: caserun
      character(30)                :: datafile
      character(132)               :: datfmt

    ! dynamic allocation of arrays

      allocate(phi(jmax))
      allocate(phi_l(jmax))
      allocate(phi_r(jmax))
      phi=0.0_kr
      phi_l=0.0_kr
      phi_r=0.0_kr

    ! build source based on options

      if (problem == 1) then
        q=0.01_kr
        bc(1)=0
        bc(2)=1
      elseif (problem == 2) then
        q=10.0_kr
        bc(1)=0
        bc(2)=0
      endif

    ! set cross sections

      sigt=100.0_kr
      sigs=c*sigt

    ! set quadrature

      call quad(n,mu,w)

    ! solve fixed-source problem

      if (sol == 0) then
        call solve_ld(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r)
      elseif (sol == 1) then
        call solve_lc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r)
      elseif (sol == 2) then
        call solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r)
      else
        write(0,'(a)') ' Incorrect solution scheme selected.'
        stop
      endif

    ! write solution into file

      datfmt='(2(es12.5))'
      write(prob,'(i1)') problem
      write(cells,'(i6)') jmax
      write(caserun,'(a)') '-p' // prob // '-' // trim(adjustl(cells)) // '-' // solopt(sol)
      write(datafile,'(a)') 'numave' // trim(adjustl(caserun)) // '.dat'
      open(unit=1,file=datafile,action='write',status='unknown')
      do j=1,jmax
        write(1,datfmt) h*(j-0.5_kr),phi(j)
      enddo
      close(1)
      write(datafile,'(a)') 'numedg' // trim(adjustl(caserun)) // '.dat'
      open(unit=1,file=datafile,action='write',status='unknown')
      do j=1,jmax
        write(1,datfmt) h*(j-1),phi_l(j)
      enddo
      write(1,datfmt) h*(jmax),phi_r(jmax)
      close(1)

    ! clean up arrays

      deallocate(phi)
      deallocate(phi_l)
      deallocate(phi_r)

    end subroutine solve_slab

    subroutine quad(n,mu,w)

      use global

      implicit none

      integer(4),    intent(in)    :: n
      real(kind=kr), intent(out)   :: mu(n/2)
      real(kind=kr), intent(out)   :: w (n/2)

      integer(4)                   :: j
      integer(4),    parameter     :: nmaxp=300
      real(kind=kr)                :: xnew(nmaxp)
      real(kind=kr)                :: wnew(nmaxp)

      xnew=0.0_kr
      wnew=0.0_kr
      call gauleg(-1.0_kr,1.0_kr,xnew,wnew,n)

      do j=1,n/2
        mu(j)=xnew(j)
        w(j) =wnew(j)
      enddo

    end subroutine quad

    subroutine gauleg(x1,x2,x,w,n)

      use global

      implicit none

      integer(4),    intent(in)    :: n
      real(kind=kr), intent(in)    :: x1
      real(kind=kr), intent(in)    :: x2
      real(kind=kr), intent(inout) :: x(n)
      real(kind=kr), intent(inout) :: w(n)

      integer(4)                   :: i
      integer(4)                   :: j
      integer(4)                   :: m
      integer(4)                   :: kount
      integer(4),    parameter     :: nmax=300
      real(kind=kr)                :: xm
      real(kind=kr)                :: xl
      real(kind=kr)                :: p1
      real(kind=kr)                :: p2
      real(kind=kr)                :: p3
      real(kind=kr)                :: pi
      real(kind=kr)                :: pp
      real(kind=kr)                :: z
      real(kind=kr)                :: z1
      real(kind=kr)                :: xtmp(nmax)  ! full set of abscissas
      real(kind=kr)                :: wtmp(nmax)  ! full set of weights
      real(kind=kr), parameter     :: eps=1.0e-30

      pi=4.0_kr*atan(1.0_kr)
      if (n > nmax) then
        write(0,'(a,1i6)') 'Gauss-Leg. integration problem --Increase PARAMETER: NMAX to at least:',n
        stop
      endif

      m=(n+1)/2
      xm=0.5_kr*(x2+x1)
      xl=0.5_kr*(x2-x1)
      do i=1,m
        z=cos(pi*(i-0.25_kr)/(n+0.5_kr))
    1   continue
        p1=1.0_kr
        p2=0.0_kr
        do j=1,n
          p3=p2
          p2=p1
          p1=((2.0_kr*j-1.0_kr)*z*p2-(j-1.0_kr)*p3)/j
        enddo
    !   p1 is now the desired Legendre polynomial. we next compute pp, its derivative,
    !   by a standard relation involving also p2, the polynomial of one lower order.
        pp=n*(z*p1-p2)/(z*z-1.0_kr)
        z1=z
        z=z1-p1/pp
        if (abs(z-z1) > eps) go to 1
        xtmp(i)=    xm-xl*z
        xtmp(n+1-i)=xm+xl*z
    !   the (n+1-i) terms are the symmetric counterparts
        wtmp(i)=2.0_kr*xl/((1.0_kr-z*z)*pp*pp)
        wtmp(n+1-i)=wtmp(i)
      enddo

    ! (half set and assumed symmetric)
      kount=0
      do i=1,n
        if (xtmp(i) >= 0.0_kr) then
          kount=kount+1
          x(kount)=xtmp(i)   ! abscissas
          w(kount)=wtmp(i)   ! weights
        endif
      enddo

    end subroutine gauleg

    subroutine solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r)

      use global

      implicit none

      integer(4),    intent(in)    :: n
      integer(4),    intent(in)    :: jmax
      integer(4),    intent(in)    :: kmax
      integer(4),    intent(in)    :: bc(2)
      real(kind=kr), intent(in)    :: h
      real(kind=kr), intent(in)    :: q
      real(kind=kr), intent(in)    :: eps
      real(kind=kr), intent(in)    :: sigt
      real(kind=kr), intent(in)    :: sigs
      real(kind=kr), intent(in)    :: mu(n/2)
      real(kind=kr), intent(in)    :: w (n/2)
      real(kind=kr), intent(inout) :: phi (jmax)
      real(kind=kr), intent(inout) :: phi_l(jmax)
      real(kind=kr), intent(inout) :: phi_r(jmax)

      integer(4)                   :: j
      integer(4)                   :: k
      integer(4)                   :: kount
      integer(4)                   :: m
      real(kind=kr)                :: tau
      real(kind=kr)                :: norm1
      real(kind=kr)                :: psi
      real(kind=kr)                :: psi_in
      real(kind=kr)                :: psi_bc(n/2)
      real(kind=kr), allocatable   :: c1(:,:)
      real(kind=kr), allocatable   :: c2(:,:)
      real(kind=kr), allocatable   :: phio(:)
      real(kind=kr), allocatable   :: phil(:)
      real(kind=kr), allocatable   :: s(:)

    ! pre-compute coeffs

      allocate(c1(jmax,n/2))
      allocate(c2(jmax,n/2))
      c1=0.0_kr
      c2=0.0_kr

      do m=1,n/2
        do j=1,jmax
          tau=sigt*h/mu(m)
          c1(j,m)=0.5_kr*tau
          c2(j,m)=0.5_kr*h/mu(m)
        enddo
      enddo

    ! solve problem

      allocate(phil(jmax))
      allocate(s(jmax))
      phi =1.0_kr
      phil=0.0_kr

      psi_in=0.0_kr
      psi_bc=0.0_kr
      kount=0

      allocate(phio(jmax))
      do k=1,kmax
        phio=phi
        do j=1,jmax
          if (problem == 1) then
            s(j)=0.5_kr*(sigs*phi(j)+q)
          elseif (problem == 2) then
            if (h*(j-0.5_kr) > 10.0_kr) then
              s(j)=0.5_kr*(sigs*phi(j))
            else
              s(j)=0.5_kr*(0.9_kr*sigs*phi(j)+q)
            endif
          endif
        enddo
        phi  =0.0_kr
        phi_l=0.0_kr
        phi_r=0.0_kr
        !$omp parallel do private(j,psi_in,psi) reduction(+:phi,phi_l,phi_r)
        do m=1,n/2
          psi_in=psi_bc(m) ! left specular bc
          if (bc(1) == 0) psi_in=0.0_kr
          do j=1,jmax
            phi_l(j)=phi_l(j)+psi_in*w(m)
            psi     =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
            phi(j)  =phi(j)+psi*w(m)
            psi_in  =2.0_kr*psi-psi_in
          enddo
          phi_r(jmax)=phi_r(jmax)+2.0_kr*psi_in*w(m)
          if (bc(2) == 0) psi_in=0.0_kr
          do j=jmax,1,-1
            psi     =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
            phi(j)  =phi(j)+psi*w(m)
            psi_in  =2.0_kr*psi-psi_in
            phi_l(j)=phi_l(j)+psi_in*w(m)
          enddo
          psi_bc(m)=psi_in
        enddo
        !$omp end parallel do

        norm1=0.0_kr
        do j=1,jmax
          norm1=norm1+(phi(j)-phio(j))**2
        enddo
        norm1=sqrt(norm1)
        if (norm1 <= eps) exit 
      enddo

      deallocate(c1)
      deallocate(c2)
      deallocate(s)
      deallocate(phio)
      deallocate(phil)

    end subroutine solve_dd

    subroutine solve_ld(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r)

      use global

      implicit none

      integer(4),    intent(in)    :: n
      integer(4),    intent(in)    :: jmax
      integer(4),    intent(in)    :: kmax
      integer(4),    intent(in)    :: bc(2)
      real(kind=kr), intent(in)    :: h
      real(kind=kr), intent(in)    :: q
      real(kind=kr), intent(in)    :: eps
      real(kind=kr), intent(in)    :: sigt
      real(kind=kr), intent(in)    :: sigs
      real(kind=kr), intent(in)    :: mu(n/2)
      real(kind=kr), intent(in)    :: w (n/2)
      real(kind=kr), intent(inout) :: phi (jmax)
      real(kind=kr), intent(inout) :: phi_l(jmax)
      real(kind=kr), intent(inout) :: phi_r(jmax)

      integer(4)                   :: j
      integer(4)                   :: k
      integer(4)                   :: m
      real(kind=kr)                :: tau
      real(kind=kr)                :: theta
      real(kind=kr)                :: norm1
      real(kind=kr)                :: psi
      real(kind=kr)                :: psil
      real(kind=kr)                :: psi_in
      real(kind=kr)                :: psi_out
      real(kind=kr)                :: psi_bc(n/2)
      real(kind=kr), allocatable   :: alpha(:,:)
      real(kind=kr), allocatable   :: c1(:,:)
      real(kind=kr), allocatable   :: c2(:,:)
      real(kind=kr), allocatable   :: phio(:)
      real(kind=kr), allocatable   :: phil(:)
      real(kind=kr), allocatable   :: philo(:)
      real(kind=kr), allocatable   :: s(:)
      real(kind=kr), allocatable   :: sl(:)

    ! lumping parameter (1.0_kr/3.0_kr is LD)
      theta=1.0_kr

    ! pre-compute coeffs

      allocate(alpha(jmax,n/2))
      allocate(c1(jmax,n/2))
      allocate(c2(jmax,n/2))
      alpha=0.0_kr
      c1   =0.0_kr
      c2   =0.0_kr

      do m=1,n/2
        do j=1,jmax
          tau=sigt*h/mu(m)
          alpha(j,m)=1.0_kr/(1.0_kr+6.0_kr/tau)
          alpha(j,m)=1.0_kr/((1.0_kr/theta-3.0_kr)*2.0_kr/tau+1.0_kr/alpha(j,m))
          c1(j,m)   =       (2.0_kr/tau+alpha(j,m)-1.0_kr)
          c2(j,m)   =1.0_kr/(2.0_kr/tau+alpha(j,m)+1.0_kr)
        enddo
      enddo

    ! solve problem

      allocate(phil(jmax))
      allocate(s(jmax))
      allocate(sl(jmax))
      phi =1.0_kr
      phil=0.0_kr

      psi_in=0.0_kr
      psi_bc=1.0_kr

      allocate(phio(jmax))
      allocate(philo(jmax))
      do k=1,kmax
        phio=phi
        philo=phil
        do j=1,jmax
          if (problem == 1) then
            s (j)=0.5_kr*(sigs*phi (j)+q)
            sl(j)=0.5_kr*(sigs*phil(j))
          elseif (problem == 2) then
            if (h*(j-0.5_kr) > 10.0_kr) then
              s (j)=0.5_kr*(sigs*phi (j))
              sl(j)=0.5_kr*(sigs*phil(j))
            else
              s (j)=0.5_kr*(0.9_kr*sigs*phi (j)+q)
              sl(j)=0.5_kr*(0.9_kr*sigs*phil(j))
            endif
          endif
        enddo
        phi=0.0_kr
        phil=0.0_kr
        !$omp parallel do private(j,psi_in,psi_out,psi,psil) reduction(+:phi,phil)
        do m=1,n/2
          psi_in=psi_bc(m) ! left specular bc
          if (bc(1) == 0) psi_in=0.0_kr
          do j=1,jmax
            psi_out=c2(j,m)*(2.0_kr*(s(j)+alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
            psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr-alpha(j,m)*sl(j)/sigt
            psil   =psi_out-psi
            psi_in =psi_out
            phi(j) =phi(j)+psi*w(m)
            phil(j)=phil(j)+psil*w(m)
          enddo
          if (bc(2) == 0) psi_in=0.0_kr
          do j=jmax,1,-1
            psi_out=c2(j,m)*(2.0_kr*(s(j)-alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
            psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr+alpha(j,m)*sl(j)/sigt
            psil   =psi-psi_out
            psi_in =psi_out
            phi(j) =phi(j)+psi*w(m)
            phil(j)=phil(j)+psil*w(m)
          enddo
          psi_bc(m)=psi_in
        enddo
        !$omp end parallel do
        phi_l=phi-phil
        phi_r=phi+phil

        norm1=0.0_kr
        do j=1,jmax
          norm1=norm1+(phi(j)-phio(j))**2 ! +(phil(j)-philo(j))**2 is close to zero
        enddo
        norm1=sqrt(norm1)
        if (norm1 <= eps) exit 
      enddo

      deallocate(alpha)
      deallocate(c1)
      deallocate(c2)
      deallocate(phio)
      deallocate(phil)
      deallocate(philo)
      deallocate(s)
      deallocate(sl)

    end subroutine solve_ld

    subroutine solve_lc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r)

      use global

      implicit none

      integer(4),    intent(in)    :: n
      integer(4),    intent(in)    :: jmax
      integer(4),    intent(in)    :: kmax
      integer(4),    intent(in)    :: bc(2)
      real(kind=kr), intent(in)    :: h
      real(kind=kr), intent(in)    :: q
      real(kind=kr), intent(in)    :: eps
      real(kind=kr), intent(in)    :: sigt
      real(kind=kr), intent(in)    :: sigs
      real(kind=kr), intent(in)    :: mu(n/2)
      real(kind=kr), intent(in)    :: w (n/2)
      real(kind=kr), intent(inout) :: phi (jmax)
      real(kind=kr), intent(inout) :: phi_l(jmax)
      real(kind=kr), intent(inout) :: phi_r(jmax)

      integer(4)                   :: j
      integer(4)                   :: k
      integer(4)                   :: m
      real(kind=kr)                :: tau
      real(kind=kr)                :: theta
      real(kind=kr)                :: tau3
      real(kind=kr)                :: tau5
      real(kind=kr)                :: tau7
      real(kind=kr)                :: norm1
      real(kind=kr)                :: psi
      real(kind=kr)                :: psil
      real(kind=kr)                :: psi_in
      real(kind=kr)                :: psi_out
      real(kind=kr)                :: psi_bc(n/2)
      real(kind=kr), allocatable   :: alpha(:,:)
      real(kind=kr), allocatable   :: rbeta(:,:)
      real(kind=kr), allocatable   :: c1(:,:)
      real(kind=kr), allocatable   :: c2(:,:)
      real(kind=kr), allocatable   :: phio(:)
      real(kind=kr), allocatable   :: phil(:)
      real(kind=kr), allocatable   :: philo(:)
      real(kind=kr), allocatable   :: s(:)
      real(kind=kr), allocatable   :: sl(:)

    ! lumping parameter (1.0_kr/3.0_kr is LD)
      theta=1.0_kr

    ! pre-compute coeffs

      allocate(alpha(jmax,n/2))
      allocate(rbeta(jmax,n/2))
      allocate(c1(jmax,n/2))
      allocate(c2(jmax,n/2))
      alpha=0.0_kr
      rbeta=0.0_kr
      c1   =0.0_kr
      c2   =0.0_kr

      do m=1,n/2
        do j=1,jmax
          tau=sigt*h/mu(m)
          if (tau < 0.01_kr) then
            tau3=tau *tau*tau
            tau5=tau3*tau*tau
            tau7=tau5*tau*tau
            alpha(j,m)=tau/6.0_kr-tau3/360.0_kr+tau5/15120.0_kr-tau7/604800.0_kr
          else
            alpha(j,m)=1.0_kr/tanh(tau/2.0_kr)-2.0_kr/tau
          endif
          rbeta(j,m)=1.0_kr/alpha(j,m)-6.0_kr/tau
          alpha(j,m)=1.0_kr/((1.0_kr/theta-3.0_kr)*2.0_kr/tau+1.0_kr/alpha(j,m))
          c1(j,m)=       (2.0_kr/tau+alpha(j,m)-1.0_kr)
          c2(j,m)=1.0_kr/(2.0_kr/tau+alpha(j,m)+1.0_kr)
        enddo
      enddo

    ! solve problem

      allocate(phil(jmax))
      allocate(s(jmax))
      allocate(sl(jmax))
      phi =1.0_kr
      phil=0.0_kr

      psi_in=0.0_kr
      psi_bc=0.0_kr

      allocate(phio(jmax))
      allocate(philo(jmax))
      do k=1,kmax
        phio=phi
        philo=phil
        do j=1,jmax
          if (problem == 1) then
            s (j)=0.5_kr*(sigs*phi (j)+q)
            sl(j)=0.5_kr*(sigs*phil(j))
          elseif (problem == 2) then
            if (h*(j-0.5_kr) > 10.0_kr) then
              s (j)=0.5_kr*(sigs*phi (j))
              sl(j)=0.5_kr*(sigs*phil(j))
            else
              s (j)=0.5_kr*(0.9_kr*sigs*phi (j)+q)
              sl(j)=0.5_kr*(0.9_kr*sigs*phil(j))
            endif
          endif
        enddo
        phi=0.0_kr
        phil=0.0_kr
        !$omp parallel do private(j,psi_in,psi_out,psi,psil) reduction(+:phi,phil)
        do m=1,n/2
          psi_in=psi_bc(m) ! left specular bc
          if (bc(1) == 0) psi_in=0.0_kr
          do j=1,jmax
            psi_out=c2(j,m)*(2.0_kr*(s(j)+alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
            psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr-alpha(j,m)*sl(j)/sigt
            psil   =((rbeta(j,m)+1.0_kr)*psi_out+(rbeta(j,m)-1.0_kr)*psi_in)/2.0_kr-rbeta(j,m)*psi
            psi_in =psi_out
            phi(j) =phi(j)+psi*w(m)
            phil(j)=phil(j)+psil*w(m)
          enddo
          if (bc(2) == 0) psi_in=0.0_kr
          do j=jmax,1,-1
            psi_out=c2(j,m)*(2.0_kr*(s(j)-alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
            psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr+alpha(j,m)*sl(j)/sigt
            psil   =-((rbeta(j,m)+1.0_kr)*psi_out+(rbeta(j,m)-1.0_kr)*psi_in)/2.0_kr+rbeta(j,m)*psi
            psi_in =psi_out
            phi(j) =phi(j)+psi*w(m)
            phil(j)=phil(j)+psil*w(m)
          enddo
          psi_bc(m)=psi_in
        enddo
        !$omp end parallel do
        phi_l=phi-phil
        phi_r=phi+phil

        norm1=0.0_kr
        do j=1,jmax
          norm1=norm1+(phi(j)-phio(j))**2 ! +(phil(j)-philo(j))**2 is close to zero
        enddo
        norm1=sqrt(norm1)
        if (norm1 <= eps) exit 
      enddo

      deallocate(alpha)
      deallocate(rbeta)
      deallocate(c1)
      deallocate(c2)
      deallocate(phio)
      deallocate(phil)
      deallocate(philo)
      deallocate(s)
      deallocate(sl)

    end subroutine solve_lc
