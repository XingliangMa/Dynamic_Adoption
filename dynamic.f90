!*********************************
!Dynamic Adoption Model
!by Xingliang Ma. ma.xingliang@gmail.com
!*********************************

!---------------------------------
!Version: Oct 30. 2010
!Changes:
! utility function: total variance
! mean profit: back to before
! prior averaged with different sigma_ict (not yet)
! fixed intruments(lagged value)
!Results:
! not good
!---------------------------

!Version: Nov 7. 2010
!Changes:
! utility function
!---------------------------------

!Version: Dec 3. 2010
!Changes:
! initial value: from myopic iterative GMM
!---------------------------------

!Version: June 20, 2011
!Changes: 
! update to the April version: testing the problem of April version.
!---------------------------------

!Version: August 20, 2012
!Changes: 
! update to the April version: new initial values. linear regression
!---------------------------------

!Version: Feburary 16, 2015
!Changes: 
! editorial changes to put it on github
!---------------------------------

!To do:
! parallel compute
!---------------------------------



module Global_Vars
    IMPLICIT NONE

    !data
    integer, parameter :: n = 348*5, ne= n-n/5, MaxIter = 3
    real(8) :: dynamic_data(n,12), Ad_rate(n), Ad_rate_CRD(n), fmsize(n), fmsize_CRD(n), &
               lon(n), lat(n)

    real(8), parameter :: delta=0.96 !discount rate

    !state/control space
    integer,parameter :: na=51, nH=51
    real(8), parameter :: Hmax=1.0d0, Hmin=0.0d0, amax=1.0d0, amin=0.0d0
    real(8) :: a(na), H(nH), deltaH, deltaa

    ! infestation rate z
    integer, parameter :: m = 3, nz=9, Nsim = 20  ! 3*sigma; # of z;  # of simulations

    !GMM
    integer, parameter :: l = 17 ! # of moments
    real(8), dimension(n,l) :: XX  ! data to interact with predict error

    real(8), dimension(l,l) :: W  ! optimal weight matrix
    real(8), dimension(l,1) :: Gn
    real(8), dimension(l,l) :: Ghat_mean

    ! parameter name
    character(len=9) :: names(15)

    integer :: Index_NonZero(n-n/5)

    integer :: MaxNumTheta
    real(8) :: theta_all(MaxIter,15), se_all(MaxIter, 15)
    real(8) :: W_all(MaxIter+1,l,l)
    real(8) :: Ghat_mean_all(MaxIter,l,l)
    real(8) :: Gn_all(MaxIter,l)

    real(8) :: eJAC_tmp(15), eJAC_mean_all(MaxIter, 15)


Contains

FUNCTION eye(m)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: m
    REAL(8) :: eye(m,m)
    INTEGER :: i,j

    DO i=1,m
        DO j=1,m
            IF (i.EQ.j) THEN
                eye(i,j)=1.0d0
            ELSE 
                eye(i,j)=0.0d0
            END IF
        END DO
    END DO

END FUNCTION eye

!initialize the state space
subroutine InitState
    IMPLICIT NONE

    integer :: i,j,t
    !real(8) :: z

    !discretize H
    deltaH = (Hmax-Hmin)/(nH-1)
    H(1) = Hmin
    do i=2, nH
        H(i) = H(i-1) + deltaH
    end do

    !partial-adoption space
    deltaa = (amax-amin)/(na-1)
    a(1) = amin
    do i=2, na
        a(i) = a(i-1) + deltaa
    end do


    open(5, file='data/dynamic_first5_fmsize30_348.txt')
    read(5,*) ( (dynamic_data(i,j), j=1,12), i=1,n)
    close(5)

    Ad_rate = dynamic_data(:,3)
    Ad_rate_CRD = dynamic_data(:,4)

    fmsize = dynamic_data(:,5)/100.0d0
    fmsize_CRD = dynamic_data(:, 6)/100.0d0

    lon = dynamic_data(:,11)/100.0d0
    lat = dynamic_data(:,12)/100.0d0

   
    !for GMM
    XX(:,1) = 1.0d0


    do i = 1,n,5
        XX(i,2) = Ad_rate(i)
        XX(i,3) = Ad_rate_CRD(i)
        do t=1,4
            XX(i+t,2) = Ad_rate(i+t-1) ! last period adoption rate
            XX(i+t,3) = Ad_rate_CRD(i+t-1) ! last period adoption rate
        end do
    end do

    XX(:,4) = fmsize
    XX(:,5) = fmsize_CRD 

    XX(:,6) = XX(:,2)**2 
    XX(:,7) = XX(:,3)**2

    XX(:,8) = fmsize**2
    XX(:,9) = fmsize_CRD ** 2
    
    XX(:,10) = lat
    XX(:,11) = lat**2
    XX(:,12) = lon
    XX(:,13) = lon**2

    ! seed prices
    XX(:, 14:17) = dynamic_data(:,7:10)  

    W = eye(l)
    W_all(1, :, :) = W

    ! parameter name
    names = (/ 'eta_g   ',  'eta_gc  ', 'gamma_0 ', 'sigma_g ', 'sigma_xi', &
        'beta    ',  'b       ', 'lambda  ', 'sigma   ', 'gamma_1 ', &
        'c1_lat  ', 'c2_lat2 ', 'c3_lon  ', 'c4_lon2 ', 'beta_1  '/)

    ! NonZero Index for the prediction error
    j=1
    do i = 1,n,5
        do t=1,4
            Index_NonZero(j) = i+t 
            j=j+1
        end do
    end do


end subroutine InitState

end module Global_Vars

module AR1
    use probability, only :  normcdf => CDF_normal 
    use random, only : sample_uniform, set_seed
    use mytools, only : cumsum
    use Global_Vars, only : m, nz
    implicit none

contains

! compute Tran and invariant prob of a AR1 process using Tauchen(1986)
subroutine Tauchen(lambda, sigma, y, Tran, prob)     
    implicit none

    real(8), intent(in) :: lambda, sigma
    real(8), intent(out) :: y(nz), prob(nz), Tran(nz,nz)

    real(8) :: std_y, ymax, ymin, w
    real(8) :: prob_new(nz), Tol

    integer :: i, j, k


    ! discretize the space
    std_y = sigma/sqrt(1-lambda**2)    ! standard deviation of y_t

    ymax = m * std_y
    ymin = - ymax
    w = (ymax - ymin) / (nz-1)

    y(1) = ymin
    do i = 2, nz
        y(i) = y(i-1) + w
    end do

    ! calculate the transition matrix
    do j=1, nz

        do k =2, nz-1
            Tran(j,k)= normcdf(y(k)-lambda*y(j)+w/2.0d0, 0.0d0, sigma) - normcdf(y(k)-lambda*y(j)-w/2.0d0, 0.0d0, sigma)
        end do
        
        Tran(j,1) = normcdf(y(1)-lambda*y(j)+w/2.0d0, 0.0d0, sigma)
        Tran(j,nz) = 1.0d0 - normcdf(y(nz)-lambda*y(j)-w/2.0d0, 0.0d0, sigma)

    end do

    ! find the invariance probability for each y
    prob = 1.0d0/nz   ! initialize the prob of y
    Tol = 1.0d0
    do while (Tol .ge. 1.0d-6)
        prob_new = matmul(prob, Tran)
        Tol = maxval(abs(prob_new - prob))
        prob = prob_new
    end do

end subroutine Tauchen

!given a population y with density p, draw a random number Z from y
subroutine get_sample(y, p, Z, index_z)
    !use mytools, only : cumsum
    !use random
    implicit none

    real(8), intent(in) :: y(:), p(:)
    real(8), intent(out) :: Z
    integer, intent(out) :: index_z

    real(8) :: p_cum(size(p)), sim
    integer :: np

    integer :: i, j

    np = size(p)

    p_cum = cumsum(p)

    sim = sample_uniform(0.0d0, 1.0d0)

    index_z = 1
    do j = 2, np
        if ( (sim .ge. p_cum(j-1)) .and. (sim .le. p_cum(j)) ) then
            index_z = j
        end if
    end do

    Z = y(index_z)

end subroutine get_sample

! simulate an AR1 process for each farmer at each year
subroutine sim_AR1(y, prob, Tran, z , index_z)
    implicit none

    !y: population of z
    !prob: invariant probability of y
    !Tran: transition matrix
    !z: simulation for each farmer at each year
    real(8), intent(in) :: y(:), prob(:), Tran(:, :)
    real(8), intent(inout) :: z(:, :)

    integer :: index_z(size(z,1), size(z,2))   !index of z

    integer :: i, j, t, Nn, NsimZ

    Nn = size(z,1)
    NsimZ = size(z, 2)

    call set_seed(1234, 4321, 4231, 4123)
    !print*
    do i = 1, Nn, 5

        do j = 1, NsimZ  !simulate each z at steady state
            call get_sample(y, prob, z(i,j), index_z(i,j)) 
        end do

        do t = 1, 4  ! following years

            do j = 1, NsimZ  !simulate each z according to transition probability
                call get_sample(y, Tran(index_z(i+t-1, j), :),  z(i+t, j), index_z(i+t, j) ) 
            end do

        end do

    end do

end subroutine sim_AR1


subroutine sim_AR1_sim(y, prob, Tran, z , index_z)
    implicit none

    !y: population of z
    !prob: invariant probability of y
    !Tran: transition matrix
    !z: simulation for each farmer at each year
    real(8), intent(in) :: y(:), prob(:), Tran(:, :)
    real(8), intent(inout) :: z(:, :)

    integer :: index_z(size(z,1), size(z,2))   !index of z

    integer :: i, j, t, Nn, NsimZ

    Nn = size(z,1)
    NsimZ = size(z, 2)

    call set_seed(1234, 4321, 4231, 4123)
    !print*
    do i = 1, Nn, 5

        do j = 1, NsimZ  !simulate each z at steady state
            call get_sample(y, prob, z(i,j), index_z(i,j)) 
        end do

        do t = 1, 4  ! following years

            do j = 1, NsimZ  !simulate each z according to transition probability
                call get_sample(y, Tran(index_z(i+t-1, j), :),  z(i+t, j), index_z(i+t, j) ) 
            end do

        end do

    end do

end subroutine sim_AR1_sim

end module AR1

module DP
    use Global_Vars
    use AR1

    IMPLICIT NONE

Contains

! Bayesian Updating
function update(sigma_gt_old, G_old, H_old, sigma_g, sigma_xi)
    real(8) :: update
    real(8), intent(in) :: sigma_gt_old, G_old, H_old
    real(8), intent(in) :: sigma_g, sigma_xi

    update = 1.0d0 / ( 1.0d0/sigma_gt_old + G_old/sigma_g + H_old/(sigma_g+sigma_xi))

end function

! payoff
function payoff(a_it, sigma_GM, A_i, eta_gi, eta_gc, sigma_ict, sigma_g, beta_i)
    implicit none
    real(8) :: payoff
    real(8), intent(in) :: a_it, sigma_GM, A_i
    real(8), intent(in) :: eta_gi, eta_gc, sigma_ict, sigma_g, beta_i

    payoff = (eta_gi*a_it -  eta_gc/2.0d0 * a_it**2) - & 
        beta_i/2.0d0 * A_i * ( (a_it**2) * (sigma_GM + sigma_g) + ((1.0d0-a_it)**2) * sigma_ict ) 

end function payoff

! print out results
subroutine myprint(names, theta, se)
    IMPLICIT NONE 
    character(len=*), intent(in) :: names(:)
    real(8), intent(in) :: theta(:)
    real(8), intent(in), optional :: se(:)

    integer :: i

    print*
    do i = 1, size(theta)
        write(*, '(A, A, F12.4)') names(i), '    ', theta(i)
    end do

    ! if printing se
    if (present(se)) then
        do i = 1, size(theta)
            write(*, '(A, A, 2F16.4)') names(i), '    ', theta(i), se(i)
        end do
    end if

end subroutine myprint

subroutine myprint_tex(names, theta, theta_myopic, theta_initial)
    IMPLICIT NONE 
    character(len=*), intent(in) :: names(:)
    real(8), intent(in) :: theta(:)
    real(8), intent(in), optional :: theta_myopic(:)
    real(8), intent(in), optional :: theta_initial(:)

    integer :: i

    ! if printing se
    if (present(theta_myopic) .AND. present(theta_initial)) then
        do i = 1, size(theta)-1
            write(*, '(4A, F12.4, A, F12.4, A, F12.4, A)') & 
            '$\', names(i), '$', ' &', theta(i), ' &', theta_myopic(i), '  &', theta_initial(i), '  \\'
        end do

        write(*, '(4A, F12.4, A, A, A, F12.4, A)') & 
        '$\', names(size(theta)), '$', ' &', theta(size(theta)), ' &', '            ', '  &', theta_initial(size(theta)), '  \\'

    end if

    do i = 1, size(theta)
        write(*, '(4A, F12.4, A)') '$\' , names(i), ' $' ,  '    &', theta(i), ' \\'
    end do


end subroutine myprint_tex


!dynamic obj fuction: put everything in one function
function dynamic_obj_fn(theta)
    implicit none
    real(8), intent(in) :: theta(:)
    real(8) :: dynamic_obj_fn

    !Bayesian beliefs
    real(8), dimension(1:n, na, nH) :: sigma_igt

    !for the dynamic model
    integer, dimension(1:n, na, nH, nz) :: opt_ia
    real(8), dimension(1:n, na, nH, nz) ::  Vmax_next ! Vmax for next year.  (index in current year)
    real(8) :: EVmax_next 
    real(8), dimension(na) :: Vtemp 

    integer :: i,j, t, i_alast, i_Hlast, ia, i_z, itemp(1)

    !for path
    integer :: ia_current(n, Nsim)=0
    real(8) :: a_current(n, Nsim), mean_a_current(n)

    !for gmm
    real(8) :: e(n)
    real(8) :: G(n,l), Ghat(n,l,l)
    real(8) :: obj_tmp(1,1)
    real(8) :: tmp(l,1)

    !theta
    real(8) :: eta_g, eta_gc, gamma_0, sigma_g, sigma_xi, beta, b, beta_1
    real(8) :: lambda, sigma, gamma_1
    real(8) :: beta_i, eta_gi

    ! for infestation rate z
    real(8) :: y(nz), prob(nz), Tran(nz, nz), z(n, Nsim)
    integer :: index_z(n, Nsim)

    real(8) :: sigma_ict
    real(8) :: c1, c2, c3, c4

    ! interpolation
    integer :: ia_1, ia_2, iH_1, iH_2
    real(8) :: a1, H1

    eta_g = theta(1)
    eta_gc = theta(2)

    gamma_0 = theta(3)   !equal to sigma_ict if it is contant

    sigma_g = theta(4)
    sigma_xi = theta(5)

    beta = theta(6)

    b = theta(7)

    lambda = theta(8)
    sigma = theta(9)
    gamma_1 = theta(10)

    ! \eta_gi = \eta_g + c*X_i
    c1 = theta(11)
    c2 = theta(12)
    c3 = theta(13)
    c4 = theta(14)

    beta_1 = theta(15)

    
    ! some parameter restrictions
    if ( sigma_g <= 0.0d0 .or. sigma <= 0.0d0 .or. eta_gc < 0.0d0 )    then
        dynamic_obj_fn =  999999.0d0
        go to 11
    end if

    ! simulate z
    call Tauchen(lambda, sigma, y, Tran, prob)
    call sim_AR1(y, prob, Tran, z, index_z)

    !Bayesian_Updating of sigma_igt
    do i = 1, n, 5

        sigma_ict = gamma_0 + gamma_1 * 0.0d0 ! z( (nz+1)/2 )
        eta_gi = eta_g + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) 
        beta_i = beta + beta_1/fmsize(i)

        !first year prior: based on myopic base, plus a parameter b
        sigma_igt(i, :, :) = (eta_gi - eta_gc*Ad_rate(i) + fmsize(i) * beta_i * sigma_ict) &
             / ( fmsize(i) * beta_i * Ad_rate(i) ) - sigma_ict - sigma_g + b


        do t = 1,4
            do i_alast = 1, na
                do i_Hlast = 1, nH

                    sigma_igt(i+t, i_alast, i_Hlast) = & 
                    update(sigma_igt(i+t-1,i_alast, i_Hlast), &
                    fmsize(i)*a(i_alast), fmsize_CRD(i) * H(i_Hlast), sigma_g, sigma_xi)

                end do
            end do
        end do

    end do


    !dynamic programming
    do i = 5, n, 5
        eta_gi = eta_g + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) 
        beta_i = beta + beta_1/fmsize(i)

        do i_alast = 1, na
            do i_Hlast = 1, nH
                
                do i_z = 1, nz  !infestation rate

                    sigma_ict = gamma_0 + gamma_1 * z(i, i_z)

                    do ia = 1, na
                        Vtemp(ia) = 1.0d0/(1.0d0 - delta) * & 
                         payoff(a(ia), sigma_igt(i,i_alast,i_Hlast), fmsize(i), &
                         eta_gi, eta_gc, sigma_ict, sigma_g, beta_i )
                    end do
                    Vmax_Next(i, i_alast, i_Hlast, i_z) = maxval(Vtemp)

                end do

            end do
        end do
    end do

    !find Vmax and opt control
    do i=1, n, 5
       opt_ia(i, :, :, :) = 0 !first year. won't be used

       eta_gi = eta_g + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) 
       beta_i = beta + beta_1/fmsize(i)


       do t = 4, 1, -1

            do i_alast = 1, na
                do i_Hlast = 1, nH

                    do i_z = 1, nz   ! infestation rate space
                        sigma_ict = gamma_0 + gamma_1 * z(i+t, i_z)
                        EVmax_next = sum( Vmax_next(i+t, ia, i_Hlast, :) * Tran(i_z, :) )   ! EV, be careful of Trans

                        do ia = 1, na
                            Vtemp(ia) = payoff( a(ia), sigma_igt(i+t,i_alast,i_Hlast), fmsize(i), & 
                                 eta_gi, eta_gc, sigma_ict, sigma_g, beta_i)  + & 
                                  delta* EVmax_next           !assume i_H = i_Hlast
                        end do

                        Vmax_next(i+t-1, i_alast, i_Hlast, i_z) = maxval(Vtemp)
                        itemp = maxloc(Vtemp)
                        opt_ia(i+t, i_alast, i_Hlast, i_z) = itemp(1)

                    end do ! i_z


                end do ! i_Hlast
            end do ! i_alast

        end do ! t

    end do

    !simulate a path: get a_current
    do i=1,n, 5
        !at year 2000
        i_alast = floor( (Ad_rate(i) - amin)/deltaa ) + 1
        i_Hlast = floor( (Ad_rate_CRD(i) - Hmin)/deltaH ) + 1 
        
        mean_a_current(i) = Ad_rate(i)
        
        do t = 1,4
            a1 = Ad_rate(i+t-1)
            H1 = Ad_rate_CRD(i+t-1)

            do j = 1, Nsim
                !ia_current(i+t, j) = opt_ia(i+t, i_alast, i_Hlast, index_z(i+t, j))     
                ia_1 = opt_ia(i+t, i_alast, i_Hlast, index_z(i+t, j))     
                ia_2 = opt_ia(i+t, i_alast+1, i_Hlast, index_z(i+t, j))     
                iH_1 = opt_ia(i+t, i_alast, i_Hlast+1, index_z(i+t, j))     
                iH_2 = opt_ia(i+t, i_alast+1, i_Hlast+1, index_z(i+t, j))     

                ! interpolation: check Bill's note Lec 25-28 Part 1
                a_current(i+t, j) = &
                       a(ia_1) * (a(i_alast+1)-a1) * (H(i_Hlast+1)-H1) / (deltaa*deltaH) + & 
                       a(iH_1) * (a(i_alast+1)-a1) * (H1-H(i_Hlast)) / (deltaa*deltaH) + & 
                       a(iH_2) * (a1- a(i_alast)) * (H1-H(i_Hlast)) / (deltaa*deltaH) + & 
                       a(ia_2) * (a1- a(i_alast)) * (H(i_Hlast+1) - H1) / (deltaa*deltaH) 

                !a_current(i+t, j) = a( ia_current(i+t, j) )
            end do
            ! get the mean value
            mean_a_current(i+t) = sum(a_current(i+t, :)) / Nsim

            i_alast = floor( (Ad_rate(i+t) - amin)/deltaa ) + 1
            i_Hlast = floor( (Ad_rate_CRD(i+t) - Hmin)/deltaH ) + 1

        end do

    end do
     
    ! difference 
    e = mean_a_current - Ad_rate

    ! moments
    do j=1,l
        G(:,j) = XX(:,j) * e
    end do

    ! Gn
    Gn(:,1) = sum(G,1)/ne

    ! Ghat
    do i = 1, n
        tmp(:,1) = G(i,:)
        Ghat(i,:, :) = matmul(tmp, transpose(tmp) )
    end do

    Ghat_mean = sum(Ghat,1)/ne

    obj_tmp =  matmul(matmul(transpose(Gn), W), Gn)
    dynamic_obj_fn =  n*obj_tmp(1,1)

11    return

end function dynamic_obj_fn

! minimize
subroutine dynamic_min_simplex
    use minimization, only: bfgs
    use simplex
    use matrix, only : Matrix_Inverse

    implicit none
    real(8) :: obj
    real(8) :: theta0(15), theta1(15)

    integer :: i, j, Iter

    ! with initial values from linear regression model. (fixed effect)
    open(5, file='data/myopic_results_08212012.txt')
    read(5, *) ( theta0(j), j=1, 15)
    close(5)
    call myprint(names, theta0)

    !test
    obj = dynamic_obj_fn(theta0)
    print*
    print*, 'obj value at the initial theta: '
    write(*, '(f16.4)') obj


    ! minimization
    print* 
    print*, 'Minimizing by Nelder-Meade(first stage)......'
    call Nelder_Meade(theta0, 1.0d-6, dynamic_obj_fn, 0 , 100000)
    print*, '   result (first stage): '
    call myprint(names, theta0)
    theta_all(1,:) = theta0

    OPEN(2,FILE='dynamic_theta0_simplex_1.txt')
    write(2, '(15F12.6)') theta0
    close(2) 


    ! Second Stage: with efficient W
    W = Matrix_Inverse( Ghat_mean - matmul(Gn, transpose(Gn)) )
    W_all(2,:,:) = W
    theta1 = theta0
    print* 
    print*, 'Minimizing by Nelder-Meade(second stage)......'
    call Nelder_Meade(theta1, 1.0d-6, dynamic_obj_fn, 0 , 100000)
    print*, ' result (second stage): '
    call myprint(names, theta1)
    theta_all(2,:) = theta1

    OPEN(2,FILE='dynamic_theta0_simplex_2.txt')
    write(2, '(15F12.6)') theta1
    close(2) 

    OPEN(2,FILE='theta_all_backup.txt')
    do j=1, 2
        write(2, '(15F29.19)') theta_all(j,:) 
    end do


   ! iterative GMM
    print* 
    theta0 = theta1
    Iter = 3
    OPEN(2,FILE='theta_all_backup.txt')
 
    do 
        ! update the weight matrix for the next round of minimization
        W = Matrix_Inverse( Ghat_mean - matmul(Gn, transpose(Gn)) )
        W_all(Iter, :, :) = W

        CALL Nelder_Meade(theta1, 1.0d-6, dynamic_obj_fn, 0 , 100000)
        write(*,'(A, I, 2f12.6)') 'GMM stage:', Iter, maxval(abs(theta1-theta0)), theta1(5) - theta0(5)
        theta_all(Iter,:) = theta1

        write(2, '(15F29.19)') theta_all(Iter,:) 

        Iter = Iter + 1

        if ( Iter > MaxIter ) then
            print*, 'Iter exceed Max number of interation '
            exit
        end if

        theta0 = theta1
    end do
    close(2) 

    MaxNumTheta = Iter - 1
    print*, 'Total Number of Interations: ', MaxNumTheta
    OPEN(2,FILE='MaxNumTheta.txt')
    write(2, '(I)') MaxNumTheta
    close(2) 


    print*
    print*,'Output all the theta: '
    OPEN(2,FILE='theta_all.txt')
    do j=1, MaxIter 
        write(2, '(15F29.19)') theta_all(j,:) 
    end do
    close(2) 

end subroutine dynamic_min_simplex

subroutine dynamic_error(num_e, ntheta, theta, e)
    implicit none

    integer :: num_e, ntheta
    real(8) :: theta(:), e(:)

    real(8) :: e_all(n)

    !Bayesian beliefs
    real(8), dimension(1:n, na, nH) :: sigma_igt

    !for the dynamic model
    integer, dimension(1:n, na, nH, nz) :: opt_ia
    real(8), dimension(1:n, na, nH, nz) ::  Vmax_next ! Vmax for next year.  (index in current year)
    real(8) :: EVmax_next 
    real(8), dimension(na) :: Vtemp 

    integer :: i,j, t, i_alast, i_Hlast, ia, i_z, itemp(1)

    !for path
    integer :: ia_current(n, Nsim)=0
    real(8) :: a_current(n, Nsim), mean_a_current(n)

    !for gmm
    !real(8) :: e(n)
    real(8) :: G(n,l), Ghat(n,l,l)
    real(8) :: obj_tmp(1,1)
    real(8) :: tmp(l,1)

    !theta
    real(8) :: eta_g, eta_gc, gamma_0, sigma_g, sigma_xi, beta, b, beta_1
    real(8) :: lambda, sigma, gamma_1
    real(8) :: beta_i, eta_gi

    ! for infestation rate z
    real(8) :: y(nz), prob(nz), Tran(nz, nz), z(n, Nsim)
    integer :: index_z(n, Nsim)

    real(8) :: sigma_ict
    real(8) :: c1, c2, c3, c4

    ! interpolation
    integer :: ia_1, ia_2, iH_1, iH_2
    real(8) :: a1, H1

    eta_g = theta(1)
    eta_gc = theta(2)

    gamma_0 = theta(3)   !equal to sigma_ict if it is contant

    sigma_g = theta(4)
    sigma_xi = theta(5)

    beta = theta(6)

    b = theta(7)

    lambda = theta(8)
    sigma = theta(9)
    gamma_1 = theta(10)

    c1 = theta(11)
    c2 = theta(12)
    c3 = theta(13)
    c4 = theta(14)

    beta_1 = theta(15)

    
    ! simulate z
    call Tauchen(lambda, sigma, y, Tran, prob)
    call sim_AR1(y, prob, Tran, z, index_z)

    !Bayesian_Updating of sigma_igt
    do i = 1, n, 5

        sigma_ict = gamma_0 + gamma_1 * 0.0d0 ! z( (nz+1)/2 )
        eta_gi = eta_g + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) 
        beta_i = beta + beta_1/fmsize(i)

        !first year prior: based on myopic base, plus a parameter b
        sigma_igt(i, :, :) = (eta_gi - eta_gc*Ad_rate(i) + fmsize(i) * beta_i * sigma_ict) &
             / ( fmsize(i) * beta_i * Ad_rate(i) ) - sigma_ict - sigma_g + b


        ! averaged value  (do this later!!!!!)
!        do i_z = 1, Nsim
            !tmp_ict(i_z) = gamma_0 + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) + gamma_1 * z(i, i_z)
            
            !tmp_igt(i_z) = &
            !(eta_c1 * fmsize(i) + 2 * beta_i * sigma_ict + (eta_gc - eta_c1) * fmsize(i) * Ad_rate(i) ) &
                !/ ( 2 * beta_i * Ad_rate(i) )  - sigma_g - tmp_ict(i_z) + b

        !end do

        ! weighted with prob of z
        !sigma_igt(i, :, :) = sum( tmp_igt * prob )

        do t = 1,4
            do i_alast = 1, na
                do i_Hlast = 1, nH

                    sigma_igt(i+t, i_alast, i_Hlast) = & 
                    update(sigma_igt(i+t-1,i_alast, i_Hlast), &
                    fmsize(i)*a(i_alast), fmsize_CRD(i) * H(i_Hlast), sigma_g, sigma_xi)

                end do
            end do
        end do

    end do


    !dynamic programming
    do i = 5, n, 5
        eta_gi = eta_g + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) 
        beta_i = beta + beta_1/fmsize(i)

        do i_alast = 1, na
            do i_Hlast = 1, nH
                
                do i_z = 1, nz  !infestation rate

                    sigma_ict = gamma_0 + gamma_1 * z(i, i_z)

                    do ia = 1, na
                        Vtemp(ia) = 1.0d0/(1.0d0 - delta) * & 
                         payoff(a(ia), sigma_igt(i,i_alast,i_Hlast), fmsize(i), &
                         eta_gi, eta_gc, sigma_ict, sigma_g, beta_i )
                    end do
                    Vmax_Next(i, i_alast, i_Hlast, i_z) = maxval(Vtemp)

                end do

            end do
        end do
    end do

    !find Vmax and opt control
    do i=1, n, 5
       opt_ia(i, :, :, :) = 0 !first year. won't be used

       eta_gi = eta_g + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) 
       beta_i = beta + beta_1/fmsize(i)


       do t = 4, 1, -1

            do i_alast = 1, na
                do i_Hlast = 1, nH

                    do i_z = 1, nz   ! infestation rate space
                        sigma_ict = gamma_0 + gamma_1 * z(i+t, i_z)
                        EVmax_next = sum( Vmax_next(i+t, ia, i_Hlast, :) * Tran(i_z, :) )   ! EV, be careful of Trans

                        do ia = 1, na
                            Vtemp(ia) = payoff( a(ia), sigma_igt(i+t,i_alast,i_Hlast), fmsize(i), & 
                                 eta_gi, eta_gc, sigma_ict, sigma_g, beta_i)  + & 
                                  delta* EVmax_next           !assume i_H = i_Hlast
                        end do

                        Vmax_next(i+t-1, i_alast, i_Hlast, i_z) = maxval(Vtemp)
                        itemp = maxloc(Vtemp)
                        opt_ia(i+t, i_alast, i_Hlast, i_z) = itemp(1)

                    end do ! i_z


                end do ! i_Hlast
            end do ! i_alast

        end do ! t

    end do

    !simulate a path: get a_current
    do i=1,n, 5
        !at year 2000
        i_alast = floor( (Ad_rate(i) - amin)/deltaa ) + 1
        i_Hlast = floor( (Ad_rate_CRD(i) - Hmin)/deltaH ) + 1 
        
        mean_a_current(i) = Ad_rate(i)
        
        do t = 1,4
            a1 = Ad_rate(i+t-1)
            H1 = Ad_rate_CRD(i+t-1)

            do j = 1, Nsim
                !ia_current(i+t, j) = opt_ia(i+t, i_alast, i_Hlast, index_z(i+t, j))     
                ia_1 = opt_ia(i+t, i_alast, i_Hlast, index_z(i+t, j))     
                ia_2 = opt_ia(i+t, i_alast+1, i_Hlast, index_z(i+t, j))     
                iH_1 = opt_ia(i+t, i_alast, i_Hlast+1, index_z(i+t, j))     
                iH_2 = opt_ia(i+t, i_alast+1, i_Hlast+1, index_z(i+t, j))     

                ! interpolation
                a_current(i+t, j) = &
                       a(ia_1) * (a(i_alast+1)-a1) * (H(i_Hlast+1)-H1) / (deltaa*deltaH) + & 
                       a(iH_1) * (a(i_alast+1)-a1) * (H1-H(i_Hlast)) / (deltaa*deltaH) + & 
                       a(iH_2) * (a1- a(i_alast)) * (H1-H(i_Hlast)) / (deltaa*deltaH) + & 
                       a(ia_2) * (a1- a(i_alast)) * (H(i_Hlast+1) - H1) / (deltaa*deltaH) 

                !a_current(i+t, j) = a( ia_current(i+t, j) )
            end do
            ! get the mean value
            mean_a_current(i+t) = sum(a_current(i+t, :)) / Nsim

            i_alast = floor( (Ad_rate(i+t) - amin)/deltaa ) + 1
            i_Hlast = floor( (Ad_rate_CRD(i+t) - Hmin)/deltaH ) + 1

        end do

    end do
     
    ! difference 
    e_all = mean_a_current - Ad_rate
    e = e_all(Index_NonZero)

end subroutine dynamic_error


end module DP

module postEST
    use Global_Vars
    use DP
    IMPLICIT NONE

Contains

subroutine get_Estimation
    IMPLICIT NONE
    integer :: i, j


    open(5, file='theta_all.txt')
    read(5,*) ( (theta_all(i,j), j=1,15), i=1,MaxIter)
    close(5)


end subroutine get_Estimation

SUBROUTINE my_Jacobi_fast(func, nx, nf, FJAC, X, h) ! by centeral diff(omp)
    IMPLICIT NONE 
    integer :: nf, nx 
    real(8) :: X(nx), FJAC(nf, nx) 
    real(8) :: h  
    !real(8) :: h(nx)  

    INTERFACE	                    
        subroutine func(nf, nx, X, F)  
             IMPLICIT NONE		  
             integer :: nf, nx 
             REAL(8) :: X(:), F(:)
         END subroutine func
    END INTERFACE

    integer :: i, j, t
    real(8) :: X1(nx), X2(nx)
    real(8) :: F1(nf,nx), F2(nf,nx)

    X1 = X
    X2 = X

    ! compute the changed function value for each possible variation
    ! And for each possible variation, each function value is changed
    ! independently, so one time for all
    do j = 1, nx ! for each theta

        X1(j) = X(j)*(1 + h)
        call func(nf, nx, X1, F1(:,j))

        X2(j) = X(j)*(1 - h)
        call func(nf, nx, X2, F2(:,j))

        X1 = X
        X2 = X
    end do

    ! Now compute the Jacobian matrix
    do i = 1, nf
        do j = 1, nx
        
            FJAC(i,j) = (F1(i,j) - F2(i,j)) / (2*X(j)*h)
            !FJAC(i,j) = (F1(i) - F2(i)) / (2*h)
            !FJAC(i,j) = (F1(i) - F2(i)) / (2*h(j))

        end do
    end do

END SUBROUTINE my_Jacobi_fast

subroutine get_SE_Petrin(theta, WM, SE)
    use matrix, only : Matrix_Inverse, get_diagonal
    implicit none

    integer, parameter ::  ntheta= 15
    real(8), intent(in) :: theta(:), WM(:,:)
    real(8), intent(out) :: SE(ntheta)


    real(8) :: e_JAC(ne, ntheta), XXnew(ne,l)

    !real(8) :: JAC_EPS= 1.0d-4
    real(8) :: JAC_EPS=  0.1  !percentage change

    integer :: i, j, k
    real(8) :: G_tmp(ne, l, ntheta), G_hat(l, ntheta)
    real(8) :: Sigma_hat(l,l), delta_tmp(ntheta,ntheta)
    real(8) :: Avar(ntheta, ntheta)

    real(8) :: tic, toc


    print*, ' Compute e_JAC:'

    call cpu_time(tic)
    call my_Jacobi_fast(dynamic_error, ntheta, ne, e_JAC, theta, JAC_EPS)
    call cpu_time(toc)
    print *, '  Time to compute e_JAC:', (toc-tic)/60, '  minutes'

    print*, 'e_JAC averaged:'
    write(*, '(15F15.4)') sum(e_JAC, 1)/ne

    eJAC_tmp = sum(e_JAC, 1)/ne

    XXnew(:,:) = XX(Index_NonZero,:)

    !compute G_hat
    do i = 1, ne !loop over each obs
        do j = 1,l
            do k = 1, ntheta
                !G_tmp(i,j,k) = XX(i,j) * e_JAC(i,k)
                G_tmp(i,j,k) = XXnew(i,j) * e_JAC(i,k)
            end do
        end do
    end do

    G_hat = sum(G_tmp, 1)/ne
    !print*, 'G_hat is: '

    Sigma_hat =  matmul(Gn, transpose(Gn))

    delta_tmp = Matrix_Inverse( matmul(transpose(G_hat), G_hat ) )

    AVar = matmul(matmul(matmul(delta_tmp,  transpose(G_hat)), Sigma_hat), matmul(G_hat, delta_tmp) ) 
    se = sqrt(get_diagonal(AVar)/(n))

    call myprint(names, theta, se)

end subroutine get_SE_Petrin

subroutine get_SE(theta, WM, SE)
    use matrix, only : Matrix_Inverse, get_diagonal
    implicit none

    integer, parameter ::  ntheta= 15
    real(8), intent(in) :: theta(:), WM(:,:)
    real(8), intent(out) :: SE(ntheta)


    real(8) :: e_JAC(ne, ntheta), XXnew(ne,l)

    !real(8) :: JAC_EPS= 1.0d-4
    real(8) :: JAC_EPS=  0.05  !percentage change

    integer :: i, j, k
    real(8) :: G_tmp(ne, l, ntheta), G_hat(l, ntheta)
    real(8) :: Sigma_hat(l,l)
    real(8) :: Avar(ntheta, ntheta)

    real(8) :: tic, toc

    print*, ' Compute e_JAC:'

    call cpu_time(tic)
    call my_Jacobi_fast(dynamic_error, ntheta, ne, e_JAC, theta, JAC_EPS)
    call cpu_time(toc)
    print *, '  Time to compute e_JAC:', (toc-tic)/60, '  minutes'

    print*, 'e_JAC averaged:'
    write(*, '(15F15.4)') sum(e_JAC, 1)/ne

    eJAC_tmp = sum(e_JAC, 1)/ne

    XXnew(:,:) = XX(Index_NonZero,:)

    !compute G_hat
    do i = 1, ne !loop over each obs
        do j = 1,l
            do k = 1, ntheta
                !G_tmp(i,j,k) = XX(i,j) * e_JAC(i,k)
                G_tmp(i,j,k) = XXnew(i,j) * e_JAC(i,k)
            end do
        end do
    end do

    G_hat = sum(G_tmp, 1)/ne
    !print*, 'G_hat is: '

    Sigma_hat = WM
    
    

    AVar = Matrix_Inverse( matmul(matmul(transpose(G_hat), Sigma_hat), G_hat ) )
    se = sqrt(get_diagonal(AVar)/n)

    call myprint(names, theta, se)

end subroutine get_SE


! J test
subroutine get_Jtest(theta_hat)
    implicit none 

    real(8) :: theta_hat(:)
    real(8) :: Jtest


    Jtest = dynamic_obj_fn(theta_hat)

    print*
    print*, 'The J test is: '
    write(*,'(f15.4)') Jtest


end subroutine get_Jtest

! predict a path
subroutine adopt_path(theta, e, path)
    implicit none

    real(8), intent(in) :: theta(:)
    real(8), intent(out) :: e(:), path(:)

    !Bayesian beliefs
    real(8), dimension(1:n, na, nH) :: sigma_igt

    !for the dynamic model
    integer, dimension(1:n, na, nH, nz) :: opt_ia
    real(8), dimension(1:n, na, nH, nz) ::  Vmax_next ! Vmax for next year.  (index in current year)
    real(8) :: EVmax_next 
    real(8), dimension(na) :: Vtemp 

    integer :: i,j, t, i_alast, i_Hlast, ia, i_z, itemp(1)

    !for path
    integer :: ia_current(n, Nsim)=0
    real(8) :: a_current(n, Nsim), mean_a_current(n)

    !theta
    real(8) :: eta_g, eta_gc, gamma_0, sigma_g, sigma_xi, beta, b, beta_1
    real(8) :: lambda, sigma, gamma_1
    real(8) :: beta_i, eta_gi

    ! for infestation rate z
    real(8) :: y(nz), prob(nz), Tran(nz, nz), z(n, Nsim)
    integer :: index_z(n, Nsim)

    real(8) :: sigma_ict
    real(8) :: c1, c2, c3, c4

    ! interpolation
    integer :: ia_1, ia_2, iH_1, iH_2
    real(8) :: a1, H1

    eta_g = theta(1)
    eta_gc = theta(2)

    gamma_0 = theta(3)   !equal to sigma_ict if it is contant

    sigma_g = theta(4)
    sigma_xi = theta(5)

    beta = theta(6)

    b = theta(7)

    lambda = theta(8)
    sigma = theta(9)
    gamma_1 = theta(10)

    c1 = theta(11)
    c2 = theta(12)
    c3 = theta(13)
    c4 = theta(14)

    ! simulate z
    call Tauchen(lambda, sigma, y, Tran, prob)
    call sim_AR1(y, prob, Tran, z, index_z)

    !Bayesian_Updating of sigma_igt
    do i = 1, n, 5

        sigma_ict = gamma_0 + gamma_1 * 0.0d0 ! z( (nz+1)/2 )
        eta_gi = eta_g + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) 
        beta_i = beta + beta_1/fmsize(i)

        !first year prior: based on myopic base, plus a parameter b
        sigma_igt(i, :, :) = (eta_gi - eta_gc*Ad_rate(i) + fmsize(i) * beta_i * sigma_ict) &
             / ( fmsize(i) * beta_i * Ad_rate(i) ) - sigma_ict - sigma_g + b


        do t = 1,4
            do i_alast = 1, na
                do i_Hlast = 1, nH

                    sigma_igt(i+t, i_alast, i_Hlast) = & 
                    update(sigma_igt(i+t-1,i_alast, i_Hlast), &
                    fmsize(i)*a(i_alast), fmsize_CRD(i) * H(i_Hlast), sigma_g, sigma_xi)

                end do
            end do
        end do

    end do


    !dynamic programming
    do i = 5, n, 5
        eta_gi = eta_g + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) 
        beta_i = beta + beta_1/fmsize(i)

        do i_alast = 1, na
            do i_Hlast = 1, nH
                
                do i_z = 1, nz  !infestation rate

                    sigma_ict = gamma_0 + gamma_1 * z(i, i_z)

                    do ia = 1, na
                        Vtemp(ia) = 1.0d0/(1.0d0 - delta) * & 
                         payoff(a(ia), sigma_igt(i,i_alast,i_Hlast), fmsize(i), &
                         eta_gi, eta_gc, sigma_ict, sigma_g, beta_i )
                    end do
                    Vmax_Next(i, i_alast, i_Hlast, i_z) = maxval(Vtemp)

                end do

            end do
        end do
    end do

    !find Vmax and opt control
    do i=1, n, 5
       opt_ia(i, :, :, :) = 0 !first year. won't be used

       eta_gi = eta_g + c1*lat(i) + c2*(lat(i)**2) + c3*lon(i) + c4*(lon(i)**2) 
       beta_i = beta + beta_1/fmsize(i)


       do t = 4, 1, -1

            do i_alast = 1, na
                do i_Hlast = 1, nH

                    do i_z = 1, nz   ! infestation rate space
                        sigma_ict = gamma_0 + gamma_1 * z(i+t, i_z)
                        EVmax_next = sum( Vmax_next(i+t, ia, i_Hlast, :) * Tran(i_z, :) )   ! EV, be careful of Trans

                        do ia = 1, na
                            Vtemp(ia) = payoff( a(ia), sigma_igt(i+t,i_alast,i_Hlast), fmsize(i), & 
                                 eta_gi, eta_gc, sigma_ict, sigma_g, beta_i)  + & 
                                  delta* EVmax_next           !assume i_H = i_Hlast
                        end do

                        Vmax_next(i+t-1, i_alast, i_Hlast, i_z) = maxval(Vtemp)
                        itemp = maxloc(Vtemp)
                        opt_ia(i+t, i_alast, i_Hlast, i_z) = itemp(1)

                    end do ! i_z


                end do ! i_Hlast
            end do ! i_alast

        end do ! t

    end do

    !simulate a path: get a_current
    do i=1,n, 5
        !at year 2000
        i_alast = floor( (Ad_rate(i) - amin)/deltaa ) + 1
        i_Hlast = floor( (Ad_rate_CRD(i) - Hmin)/deltaH ) + 1 
        
        mean_a_current(i) = Ad_rate(i)
        
        do t = 1,4
            a1 = Ad_rate(i+t-1)
            H1 = Ad_rate_CRD(i+t-1)

            do j = 1, Nsim
                !ia_current(i+t, j) = opt_ia(i+t, i_alast, i_Hlast, index_z(i+t, j))     
                ia_1 = opt_ia(i+t, i_alast, i_Hlast, index_z(i+t, j))     
                ia_2 = opt_ia(i+t, i_alast+1, i_Hlast, index_z(i+t, j))     
                iH_1 = opt_ia(i+t, i_alast, i_Hlast+1, index_z(i+t, j))     
                iH_2 = opt_ia(i+t, i_alast+1, i_Hlast+1, index_z(i+t, j))     

                ! interpolation: check Bill's note Lec 25-28 Part 1
                a_current(i+t, j) = &
                       a(ia_1) * (a(i_alast+1)-a1) * (H(i_Hlast+1)-H1) / (deltaa*deltaH) + & 
                       a(iH_1) * (a(i_alast+1)-a1) * (H1-H(i_Hlast)) / (deltaa*deltaH) + & 
                       a(iH_2) * (a1- a(i_alast)) * (H1-H(i_Hlast)) / (deltaa*deltaH) + & 
                       a(ia_2) * (a1- a(i_alast)) * (H(i_Hlast+1) - H1) / (deltaa*deltaH) 

                !a_current(i+t, j) = a( ia_current(i+t, j) )
            end do
            ! get the mean value
            mean_a_current(i+t) = sum(a_current(i+t, :)) / Nsim

            i_alast = floor( (Ad_rate(i+t) - amin)/deltaa ) + 1
            i_Hlast = floor( (Ad_rate_CRD(i+t) - Hmin)/deltaH ) + 1

        end do

    end do
     
    ! difference 
    e = mean_a_current - Ad_rate

    path = mean_a_current

    OPEN(2,FILE='adopt_path.txt')
    do i=1, n
        write(2, '(2f10.0, 2F12.6)') dynamic_data(i,1:2), Ad_rate(i), path(i)
    end do
    close(2) 

    call adopt_path_year

end subroutine adopt_path

! to plot the path (year)
subroutine adopt_path_year
    IMPLICIT NONE

    real(8) :: adopt_result(n,4), panel_no(n), year(n), observed_rate(n), predicted_rate(n)
    real(8) :: observed_average_rate(5),  predicted_average_rate(5)

    integer :: Nfarmer = n/5
    integer :: i, j, t

    open(5, file='adopt_path.txt')
    read(5,*) ( (adopt_result(i,j), j=1,4), i=1,n)
    close(5)

    panel_no = adopt_result(:,1)
    year = adopt_result(:,2)
    observed_rate = adopt_result(:,3)
    predicted_rate = adopt_result(:,4)

    observed_average_rate = 0.0d0
    predicted_average_rate = 0.0d0

    do i = 1, n, 5
        do t=1, 5
            observed_average_rate(t) = observed_average_rate(t) + observed_rate(i+t-1)
            predicted_average_rate(t) = predicted_average_rate(t) + predicted_rate(i+t-1)
        end do
    end do

    observed_average_rate =observed_average_rate/Nfarmer
    predicted_average_rate = predicted_average_rate/Nfarmer
        
    print*
    open(2, file = 'result/adopt_path_year.txt')
    do t=1, 5
        write(2, '(I, 2f12.4)') t+2000-1, observed_average_rate(t), predicted_average_rate(t)
        write(*, '(2f12.4)') observed_average_rate(t), predicted_average_rate(t)
    end do

end subroutine adopt_path_year


end module postEST


program DynamicModel
    use Global_Vars
    use DP
    use postEST
    use matrix, only : Matrix_Inverse, get_diagonal

    IMPLICIT NONE

    real(8) :: theta0(15)
    real(8) :: predict_error(n), predict_path(n)

    real(8) :: tic, toc
    integer :: i, j, Iter

    ! initialize states
    call InitState()

    !************************************
    !Estimate Parameters
    !************************************
    call cpu_time(tic)
    call dynamic_min_simplex
    call cpu_time(toc)
    print *, 'Time to compute theta:', (toc-tic)/60, '  minutes'


    !************************************
    !Post Estimation
    !************************************
   call get_Estimation

    ! SE
    W = eye(l)
    W_all(1,:,:) = W
    do Iter=1, MaxIter-1
        call get_Jtest( theta_all(Iter, :) )

        W = Matrix_Inverse( Ghat_mean - matmul(Gn, transpose(Gn)) )
        W_all(Iter+1, :,:) = W
    end do

    !The last iteration
    Iter = MaxIter 
    W = W_all(Iter, :, :)
    call get_SE(theta_all(Iter,:), W, se_all(Iter,:) ) 
    call get_SE_Petrin(theta_all(Iter,:), W, se_all(Iter,:) ) 
    eJAC_mean_all(Iter,:) = eJAC_tmp
    call get_Jtest(theta_all(Iter,:))
    call adopt_path(theta_all(Iter,:), predict_error, predict_path)
    print*, sum(predict_error**2)

end program DynamicModel

