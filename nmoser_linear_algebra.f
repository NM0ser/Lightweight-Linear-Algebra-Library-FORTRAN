c                                                                               
      module mod_user_routines                                                                                                                                 
c
      contains
c
      subroutine buildVoigtCMat(lameShear, lameLambda, cOut, cInvOut)
      implicit none
      double precision, intent(in) :: lameShear, lameLambda
      double precision, intent(inout) :: cOut(6,6), cInvOut(6,6)
c     Builds the elastic stiffness matrix into, cOut, as well as its 
c     inverse, cInvOut. NOTE, the matrix assumes ENGINEERING strains to 
c     stresses! A factor of 2 has been added since engineering 
c     shear strains are not used. In other words, Voigt notation is 
c     adopted.
      integer :: m, n, np
      double precision :: young, poiss
      np = 3
      young = (lameShear*(3.D0*lameLambda+2.D0*lameShear))/
     &        (lameLambda+lameShear)
      poiss = lameLambda/(2.D0*(lameLambda+lameShear))
c      lameShear = young/(2.D0*(1.D0+poiss))
c      lameLambda = young*poiss/((1.D0+poiss)*(1.D0-2.D0*poiss))

      cOut = 0.D0
      do m=1,np
        do n=1,np
          if (m .ne. n) cOut(m,n) = lameLambda
        end do
        cOut(m,m) = 2.D0*lameShear + lameLambda
        cOut(m+3,m+3) = lameShear !NOTE: Lack of 2 here
      end do

      cInvOut = 0.D0
      do m=1,np
        do n=1,np
          if (m .ne. n) cInvOut(m,n) = -poiss/young
        end do
        cInvOut(m,m) = 1.D0/young
        cInvOut(m+3,m+3) = (2.D0/young)*(1.D0 + poiss) !NOTE: The 2
      end do
      end subroutine buildVoigtCMat
c      
      subroutine buildMandelCMat(lameShear, lameLambda, cOut, cInvOut)
      implicit none
      double precision, intent(in) :: lameShear, lameLambda
      double precision, intent(inout) :: cOut(6,6), cInvOut(6,6)
c     Builds the elastic stiffness matrix into, cOut, as well as its 
c     inverse, cInvOut. NOTE, the matrix assumes TENSOR strains to 
c     stresses! A factor of 2 has been added since engineering 
c     shear strains are not used. In other words, Mandel notation is 
c     adopted.
      integer :: m, n, np
      double precision :: young, poiss
      np = 3
      young = (lameShear*(3.D0*lameLambda+2.D0*lameShear))/
     &        (lameLambda+lameShear)
      poiss = lameLambda/(2.D0*(lameLambda+lameShear))
c      lameShear = young/(2.D0*(1.D0+poiss))
c      lameLambda = young*poiss/((1.D0+poiss)*(1.D0-2.D0*poiss))

      cOut = 0.D0
      do m=1,np
        do n=1,np
          if (m .ne. n) cOut(m,n) = lameLambda
        end do
        cOut(m,m) = 2.D0*lameShear + lameLambda
        cOut(m+3,m+3) = 2.D0*lameShear !NOTE: There's a 2 here
      end do

      cInvOut = 0.D0
      do m=1,np
        do n=1,np
          if (m .ne. n) cInvOut(m,n) = -poiss/young
        end do
        cInvOut(m,m) = 1.D0/young
        cInvOut(m+3,m+3) = (1.D0/young)*(1.D0 + poiss) !NOTE: Lack of 2 
      end do
      end subroutine buildMandelCMat

      subroutine buildP4Mat(P4Out)
      implicit none
      double precision, intent(inout) :: P4Out(6,6)
c     Builds the rank-4 deviatoric projection tensor. The output is in
c     P4Out(1:6,1:6). This tensor, when multiplied by a Voigt or Mandel
c     vector representing a rank-2 tensor, say vec(1:6), results in an
c     Voigt or Mandel vector, respectively that corresponds to the
c     deviatoric components of the original rank-2 tensor.
      integer :: m, n, np

      np = 3
      P4Out = 0.D0
      do m=1,np
        do n=1,np
          if (m .ne. n) P4Out(m,n) = -1.D0/3.D0
        end do
        P4Out(m,m) = 2.D0/3.D0
        P4Out(m+3,m+3) = 1.D0
      end do
      end subroutine buildP4Mat
c            
      subroutine mandel2Matrix(manVecIn, matOut)
      implicit none
c     Takes a Mandel vector and converts it back into a 3 X 3 matrix
c     with the SQRT(2) factor taken out of the deviatoric terms. Warp3D
c     follows the same conventions as Abaqus Standard. More 
c     specifically, the order of stress/strain components are:
c     1 - xx
c     2 - yy
c     3 - zz
c     4 - xy
c     5 - xz
c     6 - yz
c     Shear strains are taken to be TENSOR COMPONENTS! The user may
c     need to divide the shear strains by a factor of 2.0 before using
c     this routine if the shear strains were originally given as 
c     engineering components.  
      double precision, parameter :: SQRT2=1.414213562373095048801688724
      double precision, intent(in) :: manVecIn(6)
      double precision, intent(inout) :: matOut(3,3)
c      
      matOut(1,1) = manVecIn(1)
      matOut(2,2) = manVecIn(2)
      matOut(3,3) = manVecIn(3)
      matOut(1,2) = manVecIn(4)/SQRT2
      matOut(2,1) = manVecIn(4)/SQRT2
      matOut(1,3) = manVecIn(5)/SQRT2 ! Abaqus Standard stores
      matOut(3,1) = manVecIn(5)/SQRT2 ! the shear components in a 
      matOut(2,3) = manVecIn(6)/SQRT2 ! different order than Abaqus
      matOut(3,2) = manVecIn(6)/SQRT2 ! Explicit. WATCH OUT
      end subroutine mandel2Matrix
c
      subroutine matrix2Mandel(matIn, manVecOut)
      implicit none
c     Takes a 3 X 3 symmetric matrix and converts into a Mandel vector
c     by adding a SQRT(2) to the shear terms. Warp3D follows the same
c     conventions as Abaqus Standard. More specifically, the order of
c     stress/strain components are:
c     1 - xx
c     2 - yy
c     3 - zz
c     4 - xy
c     5 - xz
c     6 - yz
c     Shear strains are taken to be TENSOR COMPONENTS! The user may
c     need to divide the shear strains by a factor of 2.0 before using
c     this routine if the shear strains were originally given as 
c     engineering components. 
      double precision, parameter :: SQRT2=1.414213562373095048801688724
      double precision, intent(in) :: matIn(3,3)
      double precision, intent(inout) :: manVecOut(6)
c      
      manVecOut(1) = matIn(1,1)
      manVecOut(2) = matIn(2,2)
      manVecOut(3) = matIn(3,3)       ! Abaqus Standard stores
      manVecOut(4) = matIn(1,2)*SQRT2 ! the shear components in a 
      manVecOut(5) = matIn(1,3)*SQRT2 ! different order than Abaqus
      manVecOut(6) = matIn(2,3)*SQRT2 ! Explicit. WATCH OUT
      end subroutine matrix2Mandel
c
      subroutine mandelMat2VoigtMat(matIn, matOut)
      implicit none
c     Takes a matrix in Mandel notation and converts it to Voigt (for 
c     what is usually done with the kinematic terms).    
      double precision, parameter :: SQRT2=1.414213562373095048801688724
      double precision, intent(in) :: matIn(6,6) 
      double precision, intent(inout) :: matOut(6,6)
      integer :: m, n
      double precision :: matTemp(6,6)
c      
      matOut = 0.D0
      matTemp = 0.D0
      do n=1,3
        do m=1,3
          matTemp(m,n) = matIn(m,n)
        end do
      end do

      do n=4,6 
        do m=1,3
          matTemp(m,n) = matIn(m,n)/SQRT2
        end do
      end do

      do n=1,3 
        do m=4,6
          matTemp(m,n) = matIn(m,n)/SQRT2
        end do
      end do

      do n=4,6
        do m=4,6
          matTemp(m,n) = matIn(m,n)/2.D0
        end do
      end do

      matOut(:,:) = matTemp(1:6,1:6)
      end subroutine mandelMat2VoigtMat
c
      subroutine voigtMat2MandelMat(matIn, matOut)
      implicit none
c     Takes a matrix in Voigt notation (for the kinematic terms) and 
c     converts it to Mandel 
      double precision, parameter :: SQRT2=1.414213562373095048801688724
      double precision, intent(in) :: matIn(6,6)
      double precision, intent(inout) :: matOut(6,6)
      double precision :: matTemp(6,6)
      integer :: m, n

      matOut = 0.D0
      matTemp = 0.D0
      do n=1,3
        do m=1,3
          matTemp(m,n) = matIn(m,n)
        end do
      end do

      do n=4,6 
        do m=1,3
          matTemp(m,n) = matIn(m,n)*SQRT2
        end do
      end do

      do n=1,3 
        do m=4,6
          matTemp(m,n) = matIn(m,n)*SQRT2
        end do
      end do

      do n=4,6
        do m=4,6
          matTemp(m,n) = matIn(m,n)*2.D0
        end do
      end do

      matOut(:,:) = matTemp(1:6,1:6)
      end subroutine voigtMat2MandelMat
c
      subroutine mandelVec2VoigtVec(vecIn, vecOut, useTwo)
      implicit none
c     Takes a vector in Mandel notation and converts it to a Voigt 
c     vector. When useTwo < 1, then it is assumed that the Voigt vector
c     is stored in tensor components. Hence, only a factor of sqrt(2) 
c     is needed on the shear terms. If useTwo >= 1, then a factor of two 
c     will be included on the shear terms to accomodate for engineering 
c     (kinematic) components. 
      double precision, parameter :: SQRT2=1.414213562373095048801688724
      integer, intent(in) :: useTwo
      double precision, intent(in) :: vecIn(6)
      double precision, intent(inout) :: vecOut(6)
c
      vecOut = vecIn
      if (useTwo .lt. 1) then
        vecOut(4) = vecOut(4)/SQRT2
        vecOut(5) = vecOut(5)/SQRT2
        vecOut(6) = vecOut(6)/SQRT2
      else
        vecOut(4) = (2.0*vecOut(4))/SQRT2
        vecOut(5) = (2.0*vecOut(5))/SQRT2
        vecOut(6) = (2.0*vecOut(6))/SQRT2
      end if
      end subroutine mandelVec2VoigtVec
c      
      subroutine voigtVec2MandelVec(vecIn, vecOut, useTwo)
      implicit none
c     Takes a vector in Voigt notation and converts it to a Mandel 
c     vector. When useTwo < 1, then it is assumed that the Voigt vector
c     is stored in tensor components. Hence, only a factor of sqrt(2) 
c     is needed on the shear terms. If useTwo >= 1, then a factor of two 
c     will be included on the shear terms to accomodate for engineering 
c     (kinematic) components.
      double precision, parameter :: SQRT2=1.414213562373095048801688724
      integer, intent(in) :: useTwo
      double precision, intent(in) :: vecIn(6)
      double precision, intent(inout) :: vecOut(6)
c
      vecOut = vecIn
      if (useTwo .lt. 1) then
        vecOut(4) = vecOut(4)*SQRT2
        vecOut(5) = vecOut(5)*SQRT2
        vecOut(6) = vecOut(6)*SQRT2
      else
        vecOut(4) = (vecOut(4)/2.0)*SQRT2
        vecOut(5) = (vecOut(5)/2.0)*SQRT2
        vecOut(6) = (vecOut(6)/2.0)*SQRT2
      end if
      end subroutine voigtVec2MandelVec
c      
      subroutine lubksb(a, np, indx, b)
      implicit none
      integer, intent(in) :: np, indx(np)
      double precision, intent(in) :: a(np,np)
      double precision, intent(inout) :: b(np)
c     Solves the set of np linear equations A · X = B. Here a is input, 
c     not as the matrix A but rather as its LU decomposition, determined
c     by the routine ludcmp. indx is input as the permutation vector 
c     returned by ludcmp. b(1:n) is input as the right-hand side vector 
c     B, and returns with the solution vector X. a, n, np, and indx are 
c     not modified by this routine and can be left in place for 
c     successive calls with different right-hand sides b. This routine
c     takes into account the possibility that b will begin with many 
c     zero elements, so it is efficient for use in matrix inversion.
c
      integer :: i, ii, j, ll, n
      double precision :: summ
c      
      n = np
c     When ii is set to a positive value, it will become the index
c     of the first nonvanishing element of b. We now do the forward
c     substitution while also unscrambling the permutation as we go
      ii = 0                      
      do i=1,n                 
        ll = indx(i)              
        summ = b(ll)
        b(ll) = b(i)
        if (ii .ne. 0) then
          do j=ii,i-1
            summ = summ - a(i,j)*b(j)
          end do
        else if (summ .ne. 0.) then
          ii = i   ! A nonzero element is encountered, so from now on                
        endif      ! we will have to do the sums in the loop above.               
        b(i) = summ
      end do
c
      do i=n,1,-1  ! Now do the backsubstitution
        summ = b(i)
        do j=i+1,n
          summ = summ - a(i,j)*b(j)
        end do
        b(i) = summ/a(i,i) ! Store a component of the solution vector X
      end do
      end subroutine lubksb
c      
      subroutine ludcmp(a, np, indx, d, singFlag)
      implicit none
      integer, intent(in) :: np 
      integer, intent(inout) :: indx(np), singFlag
      double precision, intent(inout) :: a(np,np), d 
c  
c     Largest expected n, and a small number.   
      integer, parameter :: NMAX = 500
      double precision, parameter :: TINY = 1.0e-20
c      
c     Given a matrix a(1:np,1:np), with physical dimension np by np,  
c     this routine replaces it by the LU decomposition of a rowwise 
c     permutation of itself. a and n are input. a is output, arranged as 
c     in equation (2.3.14) above; indx(1:n) is an output vector that 
c     records the row permutation effected by the partial pivoting; d is
c     output as ±1 depending on whether the number of row interchanges 
c     was even or odd, respectively. This routine is used in combination
c     with lubksb to solve linear equations or invert a matrix. See
c     Numerical Recipes in FORTRAN 77 for more details on the algorithms

      integer :: i, imax, j, k, n
      double precision :: aamax, dum, summ, vv(NMAX) 
c     vv stores the implicit scaling of each row.
c
      n = np
      d = 1.        ! No row interchanges yet.
      singFlag = 0  ! Set to 1 if singular
      do i=1,n  ! Loop over rows to get the implicit scaling information
        aamax = 0.
        do j=1,n
          if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
        end do
c
        if (aamax .eq. 0.) then ! No nonzero largest element.
          !write (59,*) 'singular matrix in ludcmp' 
          singFlag = 1
          return ! Exit subroutine and delete element
          !read (*,*)
        endif
        vv(i) = 1./aamax  ! Save scaling
      end do
c
      do j=1,n ! This is the loop over columns in Crout's method
        do i=1,j-1
          summ=a(i,j)
          do k=1,i-1
            summ = summ - a(i,k)*a(k,j)
          end do
          a(i,j) = summ
        end do
c
        aamax = 0. ! Initialize for the search for largest pivot element
        do i=j,n
          summ = a(i,j)
          do k=1,j-1
            summ = summ - a(i,k)*a(k,j)
          end do
          a(i,j) = summ
          dum = vv(i)*abs(summ)     ! Figure of merit for the pivot.
          if (dum .ge. aamax) then  ! Is it better than the best so far?
            imax = i
            aamax = dum
          endif
        end do
c  
        if (j .ne. imax) then       ! Do we need to interchange rows?
          do k=1,n               ! Yes, do so...
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          end do
          d = -d                    ! ...and change the parity of d.
          vv(imax) = vv(j)          ! Also interchange the scale factor.
        endif
        indx(j) = imax
        if (a(j,j) .eq. 0.) a(j,j) = TINY
c       If the pivot element is zero, the matrix is singular (at least 
c       to the precision of the algorithm). For some applications on 
c       singular matrices, it is desirable to substitute TINY for zero.
        if (j .ne. n) then ! Now, finally, divide by the pivot element
          dum = 1./a(j,j)
          do i=j+1,n
            a(i,j) = a(i,j)*dum
          end do
        endif
      end do ! Go back for the next column in the reduction.
      end subroutine ludcmp
c
      subroutine calcinv(aMat, np, singFlag, aInv, aLU)
      implicit none
      integer, intent(in) :: np

      integer, intent(inout) :: singFlag
      double precision, intent(in) :: aMat(np,np)
      double precision, intent(inout), optional :: aLU(np,np)
      double precision, intent(inout) :: aInv(np,np)
c     Take the inverse of aMat by solving a linear system equal to the 
c     identity matrix. The inverse is calculated using the LU-
c     decomposition method, which uses ludcmp(...). The resultant LU-
c     decomposed matrix is outputted in aLU. Amat is of size
c     np by np, and the inverse is returned in aInv. If aMat is 
c     singular, and thus no inverse exists, singFlag will be equal to
c     1; implying that if singFlag == 0, the inverse was a success.
      integer :: i, m
      integer :: indx(np)
      double precision :: d, aLU_cp(np,np)

      aInv = 0.
      aLU_cp = aMat
      singFlag = 0.
      m = np

      call ludcmp(aLU_cp, m, indx, d, singFlag)
      if (present(aLU)) then
        aLU = aLU_cp
      end if
      if (singFlag .eq. 1) then
        return
      end if

      call eyeMat(aInv, m)
      do i=1,m
        call lubksb(aLU_cp, m, indx, aInv(1,i))
      end do
      end subroutine calcinv
c
      subroutine luinv(aLU, np, indx, aInv)
      implicit none
      integer, intent(in) :: np, indx(np)   
      double precision, intent(in) :: aLU(np,np)
      double precision, intent(inout) :: aInv(np,np) 
c     Calculates the inverse of aLU and outputs it in aInv. aLU must
c     already be decomposed into an LU-matrix; use ludcmp(...) to do so.
c     indx is an array of integers also outputted by ludcmp(...) that is
c     required here. The size of aLU is expected to np by np, and it is 
c     assumed that aLU is invertible (not singular).  
      integer :: i, m
      double precision :: aLU_cp(np,np)
      m = np

      aLU_cp = aLU
      call eyeMat(aInv, m)
      do i=1,m
        call lubksb(aLU_cp, m, indx, aInv(1,i))
      end do
      end subroutine luinv  
c      
      subroutine ludet(aLU, np, d, aDet)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: aLU(np,np), d
      double precision, intent(inout) :: aDet
c     Works with the output of ludcmp(a, np, indx, d) where matrix 
c     a(1:np) gets destroyed and replaced with a row-permutation of the 
c     LU-decomposition of matrix a(1:np), denoted as input here, 
c     aLU(1:np). Then, d (+1 or -1) is also calculated based on if the
c     number of row permutations are even or odd. This algorithm 
c     attempts to locally  scale the diagonals of aLU while calculating
c     the determinant in order to prevent an overflow. The resultant 
c     determinant is stored into aDet.
c
      integer :: i, m
      double precision :: rho, K, u(np), uMag, nn
c
      m = np
      uMag = 0.
      do i=1,m ! Collect the diagonal components of aLU
        u(i) = aLU(i,i)
      end do
      call vecmag(u, m, uMag)
c
      rho = 1.
      do i=1,m
        rho = rho*(u(i)/uMag)
      end do
      rho = d*rho
c
      nn = np
      K = 0.
      do i=1,m
        K = K + log(uMag)
      end do
      K = exp((1./nn)*K)
c
      aDet = rho*(K**np)
c
c     Original, simpler method, but could cause overflow
c      m = np
c      aDet = d
c      do i=1,m
c        aDet = aDet*aLU(i,i)
c      end do
      end subroutine ludet     
c
      subroutine calcdet(Amat, np, Adet)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: Amat(np,np)
      double precision, intent(inout) :: Adet
c     This is a quicker method of calculating the determinant of
c     either 1) a scaler value (np = 1); 2) a 2D matrix (np = 2);
c     or 3) a 3D matrix (np = 3). Otherwise, the determinant will be
c     calculated from an scaled LU-decomposition method. The 
c     resultant determinant is stored in Adet.
c
      integer :: i, m, failFlag, indxVec(np)
      double precision :: d, Amat_cp(np,np)
      m = np
      Adet = 0.
c
      if (np .eq. 1) then
        Adet = Amat(1,1)
      elseif (np .eq. 2) then
        Adet = Amat(1,1)*Amat(2,2) - Amat(1,2)*Amat(2,1)
      elseif (np .eq. 3) then
        Adet = -Amat(1,3)*Amat(2,2)*Amat(3,1) +
     1         Amat(1,2)*Amat(2,3)*Amat(3,1) +
     2         Amat(1,3)*Amat(2,1)*Amat(3,2) -
     3         Amat(1,1)*Amat(2,3)*Amat(3,2) -
     4         Amat(1,2)*Amat(2,1)*Amat(3,3) +
     5         Amat(1,1)*Amat(2,2)*Amat(3,3)
      else
        Amat_cp = Amat ! Make a copy
        call ludcmp(Amat_cp, np, indxVec, d, failFlag)
        call ludet(Amat_cp, np, d, Adet)
      endif
      end subroutine calcdet
c
      subroutine getInvar(Amat, np, InvarVec)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: Amat(np,np)
      double precision, intent(inout) :: InvarVec(3)
c     Calculates the 3 invariants of matrix Amat(1:np,1:np), and stores
c     them into InvarVec(1:3); InvarVec(1) = I1, InvarVec(2) = I2, and
c     InvarVec(3) = I3.
c
      integer :: i, j, m
      double precision :: I1, I2, summ, I3, Amattr(np,np)
c
      m = np
      I1 = 0.
      I2 = 0.
      I3 = 0.
c
      do i=1,m  ! I1 = trace[Amat]
        I1 = I1 + Amat(i,i)
      end do
c
      call matrixtrnps(Amat, m, m, Amattr)
      call matdbldotmat(Amat, Amattr, m, summ)
      I2 = 0.5*(I1*I1 - summ)
c
      call calcdet(Amat, np, I3)
c
      InvarVec(1) = I1
      InvarVec(2) = I2
      InvarVec(3) = I3
      end subroutine getInvar
c
      subroutine eyeMat(IMat, np)
      implicit none
c     Function eye(np) creates an np by np identity matrix into
c     IMat(1:np, 1:np), which gets destroyed in the process.
      integer, intent(in) :: np
      double precision, intent(inout) :: IMat(np,np)
      integer i, m
c      
      IMat = 0.D0
      m = np
      do i=1,m
        IMat(i,i) = 1.D0
      end do
      end subroutine eyeMat     
c
      subroutine vecdotvec(a, b, np, c)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: a(np), b(np)
      double precision, intent(inout) :: c
c     The inner dot product between vector a(1:np) and b(1:np), and 
c     stores it into c.
c      integer :: i, m
c      m = np
      c = 0.
c
c     do i=1,m
c       c = c + a(i)*b(i)
c     end do
c
      c = dot_product(a,b)
      end subroutine vecdotvec
c
      subroutine vecoutervec(a, b, np, c)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: a(np), b(np)
      double precision, intent(inout) :: c(np,np)
c     Two vectors, a(1:np) and b(1:np), are multiplied using the outer
c     product in order to produce a matrix, c(1:np, 1:np).
      integer :: i, j, m
      c = 0.D0
      m = np
      do j=1,m
        do i=1,m
          c(i,j) = a(i)*b(j)
        end do
      end do
      end subroutine vecoutervec
c
      subroutine vecmag(bvec, np, bmag)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: bvec(np)
      double precision, intent(inout) :: bmag
c     Calculates the vector magnitude of bvec(1:np) and stores it
c     into bmag
c      integer :: i, m
c      m = np
      bmag = 0.

c     do i=1,m
c       bmag = bmag + bvec(i)*bvec(i)
c     end do
c     bmag = sqrt(bmag)
c
      bmag = norm2(bvec)
      end subroutine vecmag      
c
      subroutine matFrobeniusnorm(Amat, np, L2Out)
      implicit none
c     Calculate the Frobenius norm of matrix Amat, and output in L2Out      
      integer, intent(in) :: np
      double precision, intent(in) :: Amat(np,np)
      double precision, intent(inout) :: L2Out
c      
      integer :: i, j, m
      double precision :: summ
c     
      L2Out = 0.D0
      summ = 0.D0
      m = np
      do i=1,m
        do j=1,m
          summ = summ + Amat(j,i)**2
        end do
      end do
c      
      L2Out = sqrt(summ)
      end subroutine matFrobeniusnorm
c
      subroutine matdotvec(Amat, bvec, np, cvec)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: Amat(np,np), bvec(np)
      double precision, intent(inout) :: cvec(np)
c     Calculates A [dot] b = c, or index notation, c_i = A_ij * b_j
c      integer :: i, j, m
c
c      m = np
      cvec = 0.
c     do j=1,m
c       do i=1,m
c         cvec(i) = cvec(i) + Amat(i,j)*bvec(j)
c       end do
c     end do
c
      cvec = MATMUL(Amat, bvec)
      end subroutine matdotvec
c
      subroutine vecdotmat(Amat, bvec, np, cvec)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: Amat(np,np), bvec(np)
      double precision, intent(inout) :: cvec(np)
c     Calculates b [dot] A = c, or index notation, c_j = A_ij * b_i
c      integer :: i, j, m
c      m = np

      cvec = 0.
c     do j=1,m
c       do i=1,m
c         cvec(j) = cvec(j) + Amat(i,j)*bvec(i)
c       end do
c     end do
      cvec = MATMUL(bvec, Amat)
      end subroutine vecdotmat
c      
      subroutine matdotmat(Amat, Bmat, np, Cmat)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) ::  Amat(np,np), Bmat(np,np)
      double precision, intent(inout) :: Cmat(np,np)
c     Calculates the inner product between two square matrices,
c     Amat and Bmat. The resultant matrix is, Cmat(1:np,1:np).
c     In index notation, C_ij = A_ik * B_kj, where it is implied 
c     that there is a sum on k.
c
c      integer :: i, j, k, m
c      m = np
      Cmat = 0.D0
c      
c     do j=1,m
c       do k=1,m
c         do i=1,m
c           Cmat(i,j) = Cmat(i,j) + Amat(i,k)*Bmat(k,j) 
c         end do
c       end do
c     end do
c
      Cmat = matmul(Amat, Bmat)  
      end subroutine matdotmat
c
      subroutine matdbldotmat(Amat, Bmat, np, c)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) ::  Amat(np,np), Bmat(np,np)
      double precision, intent(inout) :: c
c     Calculates the double-dot product between two square matrices,
c     Amat and Bmat. The resultant scalar is, c. In index notation, 
c     c = A_ij * B_ij, where it is implied that there is a sum on both
c     i and j. Note that no transposes are taken.
c
      integer :: i, j, m
      m = np
      c = 0.D0
c      
      do j=1,m
        do i=1,m
          c = c + Amat(i,j)*Bmat(i,j)
        end do
      end do
      end subroutine matdbldotmat      
c
      subroutine matrixtrnps(Amat, rp, cp, Amattr)
      implicit none
      integer, intent(in) :: rp, cp
      double precision, intent(in) :: Amat(rp,cp)
      double precision, intent(inout) :: Amattr(cp,rp)
c     Transposes matrix Amat(1:rp,1:cp) and places it into 
c     Amattr(1:cp,1:rp)... i.e., Amattr_ij = Amat_ji
c      
c      integer :: i, j, m, n
c
      Amattr = 0.D0
c     m = rp
c     n = cp
c     do i=1,m
c       do j=1,n
c         Amattr(j,i) = Amat(i,j)
c       end do 
c     end do
      Amattr = transpose(Amat)
      end subroutine matrixtrnps
c
      subroutine calcPrinVals(Amat, eigsOut)
      implicit none
      double precision, intent(in) :: Amat(3,3)
      double precision, intent(inout) :: eigsOut(3)
c     Calculate the principal values of a 3D symmetric stress or 
c     strain tensor. Amat is the input tensor, assumed to be 3 by 3.
c     The principal values (i.e., eigenvalues) are calculated 
c     manually from the cubic equation based on the invariants. The 
c     three eigenvalues are outputted in eigsOut in descending order.
      integer :: i, j
      double precision, parameter :: PI=3.1415926535897932384626433832
      double precision :: Qc, Rc, I1, I2, I3, RoverQ, theta,
     & e1, e2, e3, tempA 
      double precision :: invarVec(3)
c
      eigsOut = 0. 
      do i=1,3
        do j=1,3
c         Check to ensure symmetry
          if (Amat(i,j) .ne. Amat(j,i)) return
        end do
      end do
c
      call getInvar(Amat, 3, invarVec)
      I1 = invarVec(1)
      I2 = invarVec(2)
      I3 = invarVec(3)
c
      Qc = (3.0*I2 - I1**2)/9.0
      Rc = (2.0*I1**3 - 9.0*I1*I2 + 27.0*I3)/54.0
      RoverQ = Rc/(sqrt(-Qc**3))
      theta = acos(RoverQ)
c
      e1 = 2.0*sqrt(-Qc)*cos(theta/3.0) + I1/3.0
      e2 = 2.0*sqrt(-Qc)*cos((theta + 2.0*PI)/3.0) + I1/3.0
      e3 = 2.0*sqrt(-Qc)*cos((theta + 4.0*PI)/3.0) + I1/3.0
      eigsOut(1) = e1
      eigsOut(2) = e2
      eigsOut(3) = e3
c
      do i=2,3
        tempA = eigsOut(i)
        j = i - 1
        do while (j .ge. 1)
          if (eigsOut(j) .ge. tempA) exit
          eigsOut(j+1) = eigsOut(j)
          j = j - 1
        end do
        eigsOut(j+1) = tempA
      end do
c
      end subroutine calcPrinVals
c
      subroutine tred2(a, np, d, e)
      integer, intent(in) :: np
      double precision, intent(inout) :: a(np,np), d(np), e(np)
c     Householder reduction of a real, symmetric np by np matrix a. On
c     output, a is replaced by the orthogonal matrix Q effecting the
c     transformation. d returns the diagonal elements of the tridiagonal
c     matrix, and e the off-diagonal elements, with e(1) = 0. This 
c     function assumes the user wishes to calculate the eigenvalues and
c     eigenvectors.
      integer :: i, j, k, l, n
      double precision :: f, g, h, hh, scale
c
      n = np
      do i=n,2,-1
        l = i - 1
        h = 0.
        scale=0.
        if (l .gt. 1) then 
          do k=1,l
            scale = scale + abs(a(i,k))
          end do
          if (scale .eq. 0.) then
            e(i) = a(i,l)
          else
            do k=1,l
              a(i,k) = a(i,k)/scale
              h = h + a(i,k)**2
            end do
            f = a(i,l)
            g = -sign(sqrt(h),f)
            e(i) = scale*g
            h = h - f*g
            a(i,l) = f - g
            f = 0.
            do j=1,l
              a(j,i) = a(i,j)/h
              g = 0.
              do k=1,j
                g = g + a(j,k)*a(i,k)
              end do
              do k=j+1,l
                g = g + a(k,j)*a(i,k)
              end do
              e(j) = g/h
              f = f + e(j)*a(i,j)
            end do
            hh=f/(h + h)
            do j=1,l 
              f = a(i,j)
              g = e(j) - hh*f
              e(j) = g
              do k=1,j
                a(j,k) = a(j,k) - f*e(k) - g*a(i,k)
              end do
            end do
          end if
        else
          e(i) = a(i,l)
        end if
        d(i) = h
      end do
      d(1)=0.
      e(1)=0.
      do i=1,n 
        l = i - 1
        if (d(i) .ne. 0.) then
          do j=1,l
            g = 0.
            do k=1,l
              g = g + a(i,k)*a(k,j)
            end do
            do k=1,l
              a(k,j) = a(k,j) - g*a(k,i)
            end do
          end do
        endif
        d(i) = a(i,i)
        a(i,i) = 1. 
        do j = 1,l 
          a(i,j) = 0.
          a(j,i) = 0.
        end do
      end do
      end subroutine tred2
c      
      double precision function pythag(a, b)
      double precision, intent(in) :: a, b
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
      double precision :: p, r, s, t, u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end      
c
      subroutine tqli(d, e, np, z)
      integer, intent(in) :: np
      double precision, intent(inout) :: d(np), e(np), z(np,np)
c     QL algorithm with implicit shifts, to determine the eigenvalues
c     and eigenvectors of a REAL, SYMMETRIC, TRIDIAGONAL matrix. This
c     function works with the output of tred2(...). d is a vector of 
c     length np. On input, its first np elements of the tridiagonal
c     matrix. On output, it returns the eigenvalues. The vector e inputs
c     the sub-diagonal elements of the tridiagonal matrix, with e(1)
c     arbitrary. On output, e is destroyed. This function assumes the 
c     user wishes to calculate the eigenvalues and eigenvectors. The
c     tridiagonal output matrix from tred2(...), called Q, is expected
c     as input for z(np,np). As output, z returns the eigenvectors in 
c     each column. That is, the kth column of z returns the normalized
c     eigenvector corresponding to d(k).
      integer :: i, iter, k, l, m, n
      double precision :: b, c, dd, f, g, p, r, s
      n = np
      do i=2,n
        e(i-1)=e(i)
      end do
      e(n)=0.
      do l=1,n
        iter = 0
    1   do m = l,n-1
          dd = abs(d(m)) + abs(d(m+1))
          if (abs(e(m)) + dd .eq. dd) go to 2
        end do
        m = n
    2   if (m .ne. l) then
          if (iter .eq. 30) then
c           pause 'too many iterations in tqli'            
          end if 
          iter = iter + 1
          g = (d(l+1) - d(l)) / (2.*e(l))
          r = pythag(g, 1.D0)
          g = d(m) - d(l) + e(l)/(g + sign(r,g))
          s = 1.
          c = 1.
          p = 0.
          do i=m-1,l,-1
            f = s*e(i)
            b = c*e(i)
            r = pythag(f, g)
            e(i+1) = r
            if (r .eq. 0.) then
              d(i+1) = d(i+1)-p
              e(m) = 0.
              go to 1
            end if
            s = f/r
            c = g/r
            g = d(i+1) - p
            r = (d(i) - g)*s + 2.*c*b
            p = s*r
            d(i+1) = g + p
            g=c*r - b
            do k=1,n
              f = z(k,i+1)
              z(k,i+1) = s*z(k,i) + c*f
              z(k,i) = c*z(k,i) - s*f
            end do
          end do
          d(l) = d(l) - p
          e(l) = g
          e(m) = 0.
          go to 1
        end if
      end do
      end subroutine tqli
c
      subroutine eigSys(aMat, np, eigVals, eigVecs)
      integer, intent(in) :: np
      double precision, intent(in) :: aMat(np,np)
      double precision, intent(inout) :: eigVals(np), eigVecs(np,np)
c     Calculates the eigenvalues and normalized eigenvectors of matrix
c     aMat(np,np), which must be a real, symmetric matrix. The 
c     eigenvalues are returned in eigVals(np) in descending order. The
c     corresponding eigenvectors are returned in the columns of eigVecs, 
c     such that the kth column of eigVecs returns the normalized
c     eigenvector corresponding to eigVals(k). This function acts as the
c     driver function for calling tred2(...) and tqli(...). Note, that
c     transpose(eigVecs) * aMat * eigVecs  will diagonlize aMat into 
c     its eigenvalues.
      integer i, j, k, n
      double precision :: offDiag(np), p
c      
      n = np
      eigVals = 0.
      eigVecs = aMat
c      
      call tred2(eigVecs, n, eigVals, offDiag)
      call tqli(eigVals, offDiag, n, eigVecs)
c
c     Sort eigenvalues and corresponding eigenvectors
      do i=1,n-1
        k = i
        p = eigVals(i)
        do j=i+1,n
          if (eigVals(j) .ge. p) then
            k = j
            p = eigVals(j)
          end if
        end do
        if (k .ne. i) then
          eigVals(k) = eigVals(i)
          eigVals(i) = p
          do j=1,n
            p = eigVecs(j,i)
            eigVecs(j,i) = eigVecs(j,k)
            eigVecs(j,k) = p
          end do
        end if
      end do

      end subroutine eigSys      
c                                                                               
      end module mod_user_routines                                              