      program main

      use mod_user_routines
      implicit none

c     Local variables
      character(len=80) :: msg
      character(len=1) :: newLine
c
      integer i, j, k, m, n, o
      integer indxVeca(2), indxVecb(3), indxVecc(4)
c
      double precision :: scala, scalb, scalc, scald
      double precision :: vec2a(2), vec3a(3), vec4a(4),
     &                    vec2b(2), vec3b(3), vec4b(4),
     &                    vec2c(2), vec3c(3), vec4c(4),
     &                    vec2d(2), vec3d(3), vec4d(4)
      double precision :: mat2a(2,2), mat3a(3,3), mat4a(4,4),
     &                    mat2b(2,2), mat3b(3,3), mat4b(4,4),
     &                    mat2c(2,2), mat3c(3,3), mat4c(4,4),
     &                    mat2d(2,2), mat3d(3,3), mat4d(4,4)    
c
      msg = ''
      newLine = ' '
      i = 0
      j = 0
      k = 0
      m = 0
      n = 0
      o = 0
      indxVeca = 0
      indxVecb = 0
      indxVecc = 0
      scala = 0.
      scalb = 0.
      scalc = 0.
      scald = 0.
      vec2a = 0.
      vec3a = 0.
      vec4a = 0.
      vec2b = 0.
      vec3b = 0.
      vec4b = 0.
      vec2c = 0.
      vec3c = 0.
      vec4c = 0.
      vec2d = 0.
      vec3d = 0.
      vec4d = 0.
      mat2a = 0.
      mat3a = 0.
      mat4a = 0.
      mat2b = 0.
      mat3b = 0.
      mat4b = 0.
      mat2c = 0.
      mat3c = 0.
      mat4c = 0.
      mat2d = 0.
      mat3d = 0.
      mat4d = 0.
c      
c     File that contains summary of the verification results
      open(69, file = 'verification_summary.txt',
     &  ACTION='READWRITE')
c
c     Testing matrix transpose 
      msg = 'Testing matrix transpose (matrixtrnps)...'
      write(69,999) msg
      scala = 1.0
      do j=1,4
        do i=1,4
          mat4a(i,j) = scala
          scala = scala + 1
        end do
      end do
      call matrixtrnps(mat4a, 4, 4, mat4b)
      scala = -1.0
      do j=1,4
        do i=1,4
          scalb = mat4a(i,j) - mat4b(j,i)
          if (scalb .gt. scala) scala = scalb
        end do
      end do
      write(69,999,advance='no') 'Maximum error = '
      write(69,900) scala
      write(69,999) newLine
c
c     Testing matrix double-dot matrix
      msg = 'Testing matrix double-dot matrix (matdbldotmat)...'
      write(69, 999) msg
      mat3a(1,1) = 1.0
      mat3a(2,1) = 4.0
      mat3a(3,1) = 2.0
      mat3a(1,2) = 2.0
      mat3a(2,2) = 2.0
      mat3a(3,2) = 3.0
      mat3a(1,3) = 3.0
      mat3a(2,3) = 2.0
      mat3a(3,3) = 4.0
      mat3b(1,1) = 1.0
      mat3b(2,1) = 2.0
      mat3b(3,1) = 3.0
      mat3b(1,2) = 4.0
      mat3b(2,2) = 5.0
      mat3b(3,2) = 6.0
      mat3b(1,3) = 7.0
      mat3b(2,3) = 8.0
      mat3b(3,3) = 9.0
      call matdbldotmat(mat3a, mat3b, 3, scala)
      scalb = scala - 124.0
      write(69,999,advance='no') 'Error = '
      write(69,900) scalb
      write(69,999) newLine
c
c     Testing matrix dot matrix
      msg = 'Testing matrix dot matrix (matdotmat)...'
      write(69, 999) msg
      call matdotmat(mat3a, mat3b, 3, mat3c)
      mat3d(1,1) = 14.0
      mat3d(2,1) = 14.0
      mat3d(3,1) = 20.0
      mat3d(1,2) = 32.0
      mat3d(2,2) = 38.0
      mat3d(3,2) = 47.0
      mat3d(1,3) = 50.0
      mat3d(2,3) = 62.0
      mat3d(3,3) = 74.0
      mat3c = mat3c - mat3d
      mat3c = abs(mat3c)
      scala = maxval(mat3c)
      write(69,999,advance='no') 'Maximum Error = '
      write(69,900) scala
      write(69,999) newLine
c
c     Testing vector dot matrix
      msg = 'Testing vector dot matrix (vecdotmat)...'
      write(69, 999) msg
      vec3a(1) = 4.0;
      vec3a(2) = 7.0;
      vec3a(3) = 5.0;
      call vecdotmat(mat3a, vec3a, 3, vec3b)
      vec3c(1) = 42.0
      vec3c(2) = 37.0
      vec3c(3) = 46.0
      vec3d = vec3b - vec3c
      vec3d = abs(vec3d)
      scala = maxval(vec3d)
      write(69,999,advance='no') 'Maximum Error = '
      write(69,900) scala
      write(69,999) newLine
c 
c     Testing matrix dot vector
      msg = 'Testing matrix dot vector (matdotvec)...'
      write(69, 999) msg
      call matdotvec(mat3a, vec3a, 3, vec3b)
      vec3c(1) = 33.0
      vec3c(2) = 40.0
      vec3c(3) = 49.0
      vec3d = vec3b - vec3c
      vec3d = abs(vec3d)
      scala = maxval(vec3d)
      write(69,999,advance='no') 'Maximum Error = '
      write(69,900) scala
      write(69,999) newLine  
c        
c     Testing Frobenius norm of a matrix
      msg = 'Testing Frobenius norm of a matrix (matFrobeniusnorm)...'
      write(69, 999) msg
      call matFrobeniusnorm(mat3a, 3, scala)
      scalb = scala - 8.185352771872450
      write(69,999,advance='no') 'Error = '
      write(69,900) scalb
      write(69,999) newLine
c
c     Testing vector magnitude
      msg = 'Testing vector magnitude (vecmag)...'
      write(69, 999) msg
      call vecmag(vec3a, 3, scala)
      scalb = scala - 9.486832980505138
      write(69,999,advance='no') 'Error = '
      write(69,900) scalb
      write(69,999) newLine 
c
c     Testing vector outer product
      msg = 'Testing vector outer product (vecoutervec)...'
      write(69, 999) msg
      vec3b(1) = 2.0
      vec3b(2) = 3.0
      vec3b(3) = 8.0
      call vecoutervec(vec3a, vec3b, 3, mat3c)
      mat3d(1,1) = 8.0
      mat3d(2,1) = 14.0
      mat3d(3,1) = 10.0
      mat3d(1,2) = 12.0
      mat3d(2,2) = 21.0
      mat3d(3,2) = 15.0
      mat3d(1,3) = 32.0
      mat3d(2,3) = 56.0
      mat3d(3,3) = 40.0
      mat3c = mat3c - mat3d
      mat3c = abs(mat3c)
      scala = maxval(mat3c)
      write(69,999,advance='no') 'Maximum Error = '
      write(69,900) scala
      write(69,999) newLine     
c
c     Testing vector dot vector    
      msg = 'Testing vector dot vector (vecdotvec)...'
      write(69, 999) msg
      call vecdotvec(vec3a, vec3b, 3, scala)
      scalb = scala - 69.0
      write(69,999,advance='no') 'Error = '
      write(69,900) scalb
      write(69,999) newLine 
c
c     Testing the construction of an identity matrix
      msg = 'Testing construction of the identity matrix (eyeMat)...'
      write(69, 999) msg
      call eyeMat(mat3d, 3)
      call printmatrix(mat3d, 'Generated Identity Matrix (3 X 3)', 3)
c
c     Testing the calculation of (symmetric) tensor invariants  
      msg = 'Testing the calculation of (symmetric) tensor invariants'
     &//' (getInvar)...'
      write(69, 999) msg
      mat3c(1,1) = 0.50  
      mat3c(2,1) = 0.30  
      mat3c(3,1) = 0.20  
      mat3c(1,2) = 0.30  
      mat3c(2,2) = -0.20  
      mat3c(3,2) = -0.10  
      mat3c(1,3) = 0.20  
      mat3c(2,3) = -0.10  
      mat3c(3,3) = 0.10
      call getInvar(mat3c, 3, vec3c)
      vec3d(1) = vec3c(1) - 0.40
      vec3d(2) = vec3c(2) - (-0.21)
      vec3d(3) = vec3c(3) - (-0.028)
      vec3d = abs(vec3d)
      scala = maxval(vec3d)
      write(69,999,advance='no') 'Maximum Error = '
      write(69,900) scala
      write(69,999) newLine 
c
c     Testing the calculation of principal values  
      msg = 'Testing the calculation of principal values' 
     &//' (calcPrinVals)...'
      write(69, 999) msg 
      call calcPrinVals(mat3c, vec3a)
      vec3b(1) = 0.655268520795173
      vec3b(2) = 0.115308280957800
      vec3b(3) = -0.370576801752974
      vec3c = vec3a - vec3b
      vec3c = abs(vec3c)
      scala = maxval(vec3c)
      write(69,999,advance='no') 'Maximum Error = '
      write(69,900) scala
      write(69,999) newLine 
c
c     Testing the calculation of eigenvalues and eigenvectors
      msg = 'Testing the calculation of eigenvalues and eigenvectors' 
     &//' (eigSys)...'
      write(69, 999) msg 
      call eigSys(mat3c, 3, vec3a, mat3b)
      vec3c = vec3a - vec3b
      vec3c = abs(vec3c)
      scala = maxval(vec3c)
      write(69,999,advance='no') 'Maximum Eigenvalue Error = '
      write(69,900) scala
      mat3d(1,1) = 0.916135395827559 
      mat3d(2,1) = 0.288850443288549
      mat3d(3,1) = 0.277959273858405
      mat3d(1,2) = 0.140420590746462
      mat3d(2,2) = 0.418221616273271 
      mat3d(3,2) = 0.897425616625795
      mat3d(1,3) = 0.375470363952066 
      mat3d(2,3) = 0.861194577951918 
      mat3d(3,3) = 0.342587076084199
      mat3c = abs(mat3b) - mat3d
      mat3c = abs(mat3c)
      scalb = maxval(mat3c)
      write(69,999,advance='no') 'Maximum Eigenvector Error = '
      write(69,900) scalb
      write(69,999) newLine 
c     
c      call matrixtrnps(mat3b, 3, 3, mat3d)
c      call matdotmat(mat3c, mat3b, 3, mat3a)
c      call matdotmat(mat3d, mat3a, 3, mat3c)
c      call printmatrix(mat3c, 'Diagnolized matrix', 3)
c
c     Testing LU-decomposition linear solver
      msg = 'Testing the LU-decomposition linear solver (ludcmp +'
     &//' lubksb)...'
      write(69, 999) msg
      mat4a(1,1) = 1.0
      mat4a(2,1) = 2.0
      mat4a(3,1) = 3.0
      mat4a(4,1) = -1.0
      mat4a(1,2) = 1.0
      mat4a(2,2) = 1.0
      mat4a(3,2) = -1.0
      mat4a(4,2) = 2.0
      mat4a(1,3) = 0.0
      mat4a(2,3) = -1.0
      mat4a(3,3) = -1.0
      mat4a(4,3) = 3.0
      mat4a(1,4) = 3.0
      mat4a(2,4) = 1.0
      mat4a(3,4) = 2.0
      mat4a(4,4) = -1.0
      mat4b = mat4a
      vec4a(1) = 4.0
      vec4a(2) = 1.0
      vec4a(3) = -3.0
      vec4a(4) = 4.0
      call ludcmp(mat4b, 4, indxVecc, scala, m)
      call lubksb(mat4b, 4, indxVecc, vec4a)
      vec4b(1) = -1.0
      vec4b(2) = 2.0
      vec4b(3) = 0.0
      vec4b(4) = 1.0
      vec4c = vec4b - vec4a
      vec4c = abs(vec4c)
      scalb = maxval(vec4c)
      write(69,999,advance='no') 'Maximum Error = '
      write(69,900) scalb
      write(69,999) newLine
c      
c     Testing the calculation of n-dimensional determinants 
      msg = 'Testing the calculation of n-dimensional determinants'
     &//' (calcdet)...'
      write(69, 999) msg
      mat4a(1,1) = 2.0
      mat4a(2,1) = 4.0
      mat4a(3,1) = 5.0
      mat4a(4,1) = 1.0
      mat4a(1,2) = 5.0
      mat4a(2,2) = 1.0
      mat4a(3,2) = 3.0
      mat4a(4,2) = 0.0
      mat4a(1,3) = 1.0
      mat4a(2,3) = 6.0
      mat4a(3,3) = 7.0
      mat4a(4,3) = 2.0
      mat4a(1,4) = 4.0
      mat4a(2,4) = 3.0
      mat4a(3,4) = 2.0
      mat4a(4,4) = 4.0
      call calcdet(mat4a, 4, scala)
      scalb = scala - (-36.0)
      write(69,999,advance='no') 'Error = '
      write(69,900) scalb
      write(69,999) newLine 
c
c     Testing the calculation of matrix inverse 
      msg = 'Testing the calculation of matrix inverse (luinv)...'
      write(69, 999) msg
      mat4d(1,1) = 0.833333333333332
      mat4d(2,1) = -0.083333333333333
      mat4d(3,1) = -0.583333333333332
      mat4d(4,1) = 0.083333333333333
      mat4d(1,2) = 3.666666666666661
      mat4d(2,2) = -1.166666666666665
      mat4d(3,2) = -2.166666666666663
      mat4d(4,2) = 0.166666666666666
      mat4d(1,3) = -2.611111111111106
      mat4d(2,3) = 0.861111111111110
      mat4d(3,3) = 1.694444444444442
      mat4d(4,3) = -0.194444444444444
      mat4d(1,4) = -2.277777777777775
      mat4d(2,4) = 0.527777777777777
      mat4d(3,4) = 1.361111111111109
      mat4d(4,4) = 0.138888888888889      
c     Method #1
      mat4c = mat4a      
      call ludcmp(mat4c, 4, indxVecc, scala, m)
c      call printmatrix(mat4c, 'LU_mat1', 4)
      if (m .eq. 1) write(69,999) 'Fail! Found singular matrix'
      call luinv(mat4c, 4, indxVecc, mat4b)
      mat4c = mat4b - mat4d
      mat4c = abs(mat4c)
      scalb = maxval(mat4c)
      write(69,999,advance='no') 'Method 1: Maximum Error = '
      write(69,900) scalb
c     Methods #2 and #3
      call calcinv(mat4a, 4, m, mat4b, mat4c)
c      call printmatrix(mat4c, 'LU_mat2', 4)
      if (m .eq. 1) write(69,999) 'Fail! Found singular matrix'
      mat4c = mat4b - mat4d
      mat4c = abs(mat4c)
      scalb = maxval(mat4c)
      write(69,999,advance='no') 'Method 2: Maximum Error = '
      write(69,900) scalb
      call calcinv(mat4a, 4, m, mat4b)
      if (m .eq. 1) write(69,999) 'Fail! Found singular matrix'
      mat4c = mat4b - mat4d
      mat4c = abs(mat4c)
      scalb = maxval(mat4c)
      write(69,999,advance='no') 'Method 3: Maximum Error = '
      write(69,900) scalb
      write(69,999) newLine
c
c
 900  format(E15.5)     
 999  format(A)
c
      close(69)
      end program main
c
c
      subroutine printmatrix(Amat, Aname, np)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: Amat(np,np)
      character(len=*), intent(in) :: Aname
      character(len=1) :: newLine
C     Prints components of Amat(1:np, 1:np) to the screen. Aname is a
C     character string which will be used to denote the matrix name
C     during output.
      integer :: i, j, m
      newLine = ' '

      m = np
      write(69,*) Aname, ' = '
      do i=1,m
        do j=1,m
          write(69,444,advance='no') Amat(i,j)
          write(69,544,advance='no') ' '
        end do
        write(69,544) newLine
      end do
      write(69,544) newLine

 444  format(E15.5)
 544  format(A)
      end subroutine printmatrix

      subroutine printvector(bvec, bname, np)
      implicit none
      integer, intent(in) :: np
      double precision, intent(in) :: bvec(np)
      character(len=*), intent(in) :: bname
      character(len=1) :: newLine
C     Prints components of bvec(1:np) to the screen. bname is a
C     character string which will be used to denote the vector name
C     during output.
      integer :: i, m
      newLine = ' '

      m = np
      write(69,*) bname, ' = '
      do i=1,m
        write(69,445,advance='no') bvec(i)
        write(69,545,advance='no') ' '
      end do
      write(69,545) newLine
      write(69,545) newLine

 445  format(E15.5)
 545  format(A)
      end subroutine printvector