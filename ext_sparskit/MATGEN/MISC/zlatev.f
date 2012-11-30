!!--------------------------------------------------------------------------!
!! ComplexSPARSKIT                                                          ! 
!!                                                                          !
!! ComplexSPARSKIT is derived from the original LGLP library SPARSKIT:      !
!! http://www-users.cs.umn.edu/~saad/software/SPARSKIT/index.html           ! 
!! Copyright (C) 2005, the Regents of the University of Minnesota           !
!!                                                                          !
!! The original library has been modified starting from 2006                !
!! to offer floating complex support. All the additional and                !
!! modified code has been written by:                                       !
!!                                                                          !
!! Alessandro Pecchia, Gabriele Penazzi                                     !  
!! Copyright (C) 2006                                                       !
!!                                                                          !
!! and is released under LGPL 3.0                                           !  
!!                                                                          ! 
!! ComplexSPARSKIT is free software: you can redistribute it and/or modify  !
!! it under the terms of the GNU Lesse General Public License as published  !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with ComplexSPARSKIT.  If not, see                        !
!!  <http://www.gnu.org/licenses/>.                                         ! 
!!--------------------------------------------------------------------------!


      SUBROUTINE MATRF2(M,N,C,INDEX,ALPHA,NN,NZ,A,SNR,RNR,FEJLM)
C--------------------------------------------------------------------
C
C   PURPOSE
C   -------
C   The subroutine generates sparse (rectangular or square) matrices.
C   The dimensions of the matrix and the average number of nonzero
C   elements per row can be specified by the user. Moreover, the user
C   can also change the sparsity pattern and the condition number of the
C   matrix. The non-zero elements of the desired matrix will be
C   accumulated (in an arbitrary order) in the first NZ positions of
C   array A. The column and the row numbers of the non-zero element
C   stored in A(I), I=1,...,NZ, will be found in SNR(I) and RNR(I),
C   respectively. The matrix generated by this subroutine is of the
C   class F(M,N,C,R,ALPHA) (see reference).
C
C   Note: If A is the sparse matrix of type F(M,N,C,R,ALPHA), then
C
C           min|A(i,j)| = 1/ALPHA,
C
C           max|A(i,j)| = max(INDEX*N - N,10*ALPHA).
C
C
C   CONTRIBUTOR: Ernest E. Rothman
C                Cornell Theory Center/Cornell National Supercomputer
C                Facility.
C                e-mail address: BITNET:   eer@cornellf
C                                INTERNET: eer@cornellf.tn.cornell.edu
C
C   minor modifications by Y. Saad. April 26, 1990.
C
C   Note: This subroutine has been copied from the following reference.
C         The allowable array sizes have been changed.
C
C   REFERENCE: Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
C      "A testing Scheme for Subroutines Solving Large Linear Problems",
C      Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
C
C
C   INPUT PARAMETERS
C   ----------------
C   M    - Integer. The number of rows in the desired matrix.
C          N < M+1 < 9000001 must be specified.
C
C   N    - Integer. The number of columns in the desired matrix.
C          21 < N < 9000001 must be specified.
C
C   C    - Integer. The sparsity pattern can be changed by means of this
C          parameter.  10 < C < N-10  must be specified.
C
C   INDEX - Integer.  The average number of non-zero elements per row in
C           the matrix will be equal to INDEX.
C           1 < INDEX < N-C-8 must be specified.
C
C   ALPHA - Real. The condition number of the matrix can be changed
C           BY THIS PARAMETER. ALPHA > 0.0 MUST BE SPECIFIED.
C           If ALPHA is approximately equal to 1.0 then the generated
C           matrix is well-conditioned. Large values of ALPHA will
C           usually produce ill-conditioned matrices. Note that no
C           round-off errors during the computations in this subroutine
C           are made if ALPHA = 2**I (where I is an arbitrary integer
C           which produces numbers in the machine range).
C
C   NN    - Integer. The length of arrays A, RNR, and SNR (see below).
C           INDEX*M+109 < NN < 9000001 must be specified.
C
C
C   OUTPUT PARAMETERS
C   -----------------
C   NZ    - Integer. The number of non-zero elements in the matrix.
C
C   A(NN) - Real array. The non-zero elements of the matrix generated
C           are accumulated in the first NZ locations of array A.
C
C   SNR(NN) - INTEGER array. The column number of the non-zero element
C           kept in A(I), I=1,...NZ, is stored in SNR(I).
C
C   RNR(NN) - Integer array. The row number of the non-zero element
C           kept in A(I), I=1,...NZ, is stored in RNR(I).
C
C   FEJLM - Integer. FEJLM=0 indicates that the call is successful.
C           Error diagnostics are given by means of positive values of
C           this parameter as follows:
C             FEJLM = 1    -  N       is out of range.
C             FEJLM = 2    -  M       is out of range.
C             FEJLM = 3    -  C       is out of range.
C             FEJLM = 4    -  INDEX   is out of range.
C             FEJLM = 5    -  NN      is out of range.
C             FEJLM = 7    -  ALPHA   is out of range.
C
C
C
C
      REAL*8 A, ALPHA, ALPHA1
      INTEGER M, N, NZ, C, NN, FEJLM, M1, NZ1, RR1, RR2, RR3, K
      INTEGER M2, N2
      INTEGER SNR, RNR
      DIMENSION A(NN), SNR(NN), RNR(NN)
      M1 = M
      FEJLM = 0
      NZ1 = INDEX*M + 110
      K = 1
      ALPHA1 = ALPHA
      INDEX1 = INDEX - 1
C
C  Check the parameters.
C
      IF(N.GE.22) GO TO 1
2     FEJLM = 1
      RETURN
1     IF(N.GT.9000000) GO TO 2
      IF(M.GE.N) GO TO 3
4     FEJLM = 2
      RETURN
3     IF(M.GT.9000000) GO TO 4
      IF(C.LT.11)GO TO 6
      IF(N-C.GE.11)GO TO 5
6     FEJLM = 3
      RETURN
5     IF(INDEX.LT.1) GO TO 12
      IF(N-C-INDEX.GE.9)GO TO 13
12    FEJLM = 4
13    IF(NN.GE.NZ1)GO TO 7
8     FEJLM = 5
      RETURN
7     IF(NN.GT.9000000)GO TO 8
      IF(ALPHA.GT.0.0)GO TO 9
      FEJLM = 6
      RETURN
9     CONTINUE
C
C  End of the error check. Begin to generate the non-zero elements of
C  the required matrix.
C
      DO 20 I=1,N
      A(I) = 1.0d0
      SNR(I) = I
20    RNR(I) = I
      NZ = N
      J1 = 1
      IF(INDEX1.EQ.0) GO TO 81
      DO 21 J = 1,INDEX1
      J1 = -J1
      DO 22 I=1,N
      A(NZ+I) = dfloat(J1*J*I)
      IF(I+C+J-1.LE.N)SNR(NZ+I) = I + C + J - 1
      IF(I+C+J-1.GT.N)SNR(NZ+I) = C + I + J - 1 - N
22    RNR(NZ + I) = I
21    NZ = NZ + N
81    RR1 = 10
      RR2 = NZ
      RR3 = 1
25    CONTINUE
      DO 26 I=1,RR1
      A(RR2 + I) = ALPHA*dfloat(I)
      SNR(RR2+I) = N - RR1 + I
      RNR(RR2+I) = RR3
26    CONTINUE
      IF(RR1.EQ.1) GO TO 27
      RR2 = RR2 + RR1
      RR1 = RR1 - 1
      RR3 = RR3 + 1
      GO TO 25
27    NZ = NZ + 55
29    M1 = M1 - N
      ALPHA = 1.0d0/ALPHA
      IF(M1.LE.0) GO TO 28
      N2 = K*N
      IF(M1.GE.N)M2 = N
      IF(M1.LT.N)M2 = M1
      DO 30 I=1,M2
      A(NZ+I) = ALPHA*dfloat(K+1)
      SNR(NZ + I) = I
30    RNR(NZ + I) = N2 + I
      NZ = NZ + M2
      IF(INDEX1.EQ.0) GO TO 82
      J1 = 1
      DO 41 J = 1,INDEX1
      J1 = -J1
      DO 42 I = 1,M2
      A(NZ+I) = ALPHA*dFLOAT(J*J1)*(dfloat((K+1)*I)+1.0d0)
      IF(I+C+J-1.LE.N)SNR(NZ+I) = I + C + J - 1
      IF(I+C+J-1.GT.N)SNR(NZ+I) = C + I + J - 1 - N
42    RNR(NZ + I) = N2 + I
41    NZ = NZ +M2
82    K = K + 1
      GO TO 29
28    CONTINUE
      ALPHA = 1.0d0/ALPHA1
      RR1 = 1
      RR2 = NZ
35    CONTINUE
      DO 36 I = 1,RR1
      A(RR2+I) = ALPHA*dfloat(RR1+1-I)
      SNR(RR2+I) = I
      RNR(RR2+I) = N - 10 + RR1
36    CONTINUE
      IF(RR1.EQ.10) GO TO 34
      RR2 = RR2 + RR1
      RR1 = RR1 + 1
      GO TO 35
34    NZ = NZ + 55
      ALPHA = ALPHA1
      RETURN
      END
      SUBROUTINE DCN(AR,IA,JA,N,NE,IC,NN,IERR)
C-----------------------------------------------------------------------
C
C   PURPOSE
C   -------
C   The subroutine generates sparse (square) matrices of the type
C   D(N,C).  This type of matrix has the following characteristics:
C   1's in the diagonal, three bands at the distance C above the
C   diagonal (and reappearing cyclicly under it), and a 10 x 10
C   triangle of elements in the upper right-hand corner.
C   Different software libraries require different storage schemes.
C   This subroutine generates the matrix in the storage  by
C   indices mode.
C
C
C   Note: If A is the sparse matrix of type D(N,C), then
C
C       min|A(i,j)| = 1,     max|A(i,j)| = max(1000,N + 1)
C
C
C
C   CONTRIBUTOR: Ernest E. Rothman
C                Cornell Theory Center/Cornell National Supercomputer
C                Facility.
C                e-mail address: BITNET:   eer@cornellf
C                                INTERNET: eer@cornellf.tn.cornell.edu
C
C
C   REFERENCE
C   ---------
C   1) Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
C      "A Testing Scheme for Subroutines Solving Large Linear Problems",
C       Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
C   2) Osterby, Ole and Zletev, Zahari;
C      "Direct Methods for Sparse Matrices";
C       Springer-Verlag 1983.
C
C
C
C   INPUT PARAMETERS
C   ----------------
C   N    - Integer. The size of the square matrix.
C          N > 13 must be specified.
C
C   NN   - Integer. The dimension of integer arrays IA and JA and 
C          real array AR. Must be at least NE.
C
C   IC   - Integer. The sparsity pattern can be changed by means of this
C          parameter.  0 < IC < N-12  must be specified.
C
C
C   OUTPUT PARAMETERS
C   -----------------
C   NE   - Integer. The number of nonzero elements in the sparse matrix
C          of the type D(N,C). NE = 4*N + 55.
C
C   AR(NN) - Real array. (Double precision)
C            Stored entries of a sparse matrix to be generated by this
C            subroutine.
C            NN is greater then or equal to, NE, the number of
C            nonzeros including a mandatory diagonal entry for
C            each row. Entries are stored by indices.
C
C   IA(NN) - Integer array.
C            Pointers to specify rows for the stored nonzero entries
C            in AR.
C
C   JA(NN) - Integer array.
C            Pointers to specify columns for the stored nonzero entries
C            in AR.
C
C   IERR   - Error parameter is returned as zero on successful
C             execution of the subroutine.
C             Error diagnostics are given by means of positive values
C             of this parameter as follows:
C             IERR = 1    -  N       is out of range.
C             IERR = 2    -  IC      is out of range.
C             IERR = 3    -  NN      is out of range.
C
C----------------------------------------------------------------------
C
      real*8 ar(nn)
      integer ia(nn), ja(nn), ierr
      ierr = 0
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Check the input parameters:
c
      if(n.le.13)then
         ierr = 1
         return
      endif
      if(ic .le. 0 .or. ic .ge. n-12)then
         ierr = 2
         return
      endif
      ne = 4*n+55 
      if(nn.lt.ne)then
         ierr = 3
         return
      endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c Begin to generate the nonzero elements as well as the row and column
c pointers:
c
      do 20 i=1,n
        ar(i) = 1.0d0
        ia(i) = i
        ja(i) = i
20    continue
      ilast = n
      do 30 i=1,n-ic
        it = ilast + i
        ar(it) = 1.0 + dfloat(i)
        ia(it) = i
        ja(it) = i+ic
30    continue
      ilast = ilast + n-ic
      do 40 i=1,n-ic-1
        it = ilast + i
        ar(it) = -dfloat(i)
        ia(it) = i
        ja(it) = i+ic+1
40    continue
      ilast = ilast + n-ic-1
      do 50 i=1,n-ic-2
        it = ilast + i
        ar(it) = 16.0d0
        ia(it) = i
        ja(it) = i+ic+2
50    continue
      ilast = ilast + n-ic-2
      icount = 0
      do 70 j=1,10
        do 60 i=1,11-j
         icount = icount + 1
         it = ilast + icount
         ar(it) = 100.0d0 * dfloat(j)
         ia(it) = i
         ja(it) = n-11+i+j
60    continue
70    continue
      icount = 0
      ilast = 55 + ilast
      do 80 i=n-ic+1,n
        icount = icount + 1
        it = ilast + icount
        ar(it) = 1.0d0 + dfloat(i)
        ia(it) = i
        ja(it) = i-n+ic
80    continue
      ilast = ilast + ic
      icount = 0
      do 90 i=n-ic,n
        icount = icount + 1
        it = ilast + icount
        ar(it) = -dfloat(i)
        ia(it) = i
        ja(it) = i-n+ic+1
90    continue
      ilast = ilast + ic + 1
      icount = 0
      do 100 i=n-ic-1,n
        icount = icount + 1
        it = ilast + icount
        ar(it) = 16.0d0
        ia(it) = i
        ja(it) = i-n+ic+2
100   continue
c     ilast = ilast + ic + 2
c     if(ilast.ne.4*n+55) then
c     write(*,*)' ilast equal to ', ilast
c     write(*,*)' ILAST, the number of nonzeros, should = ', 4*n + 55
c     stop
c     end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
      SUBROUTINE ECN(N,IC,NE,IA,JA,AR,NN,IERR)
C----------------------------------------------------------------------
C
C   PURPOSE
C   -------
C   The subroutine generates sparse (square) matrices of the type
C   E(N,C).  This type of matrix has the following characteristics:
C   Symmetric, positive-definite, N x N matrices with 4 in the diagonal
C   and -1 in the two sidediagonal and in the two bands at the distance
C   C from the diagonal. These matrices are similar to matrices obtained
C   from using the five-point formula in the discretization of the
C   elliptic PDE.
C
C
C   Note: If A is the sparse matrix of type E(N,C), then
C
C       min|A(i,j)| = 1,     max|A(i,j)| = 4
C
C
C
C   CONTRIBUTOR: Ernest E. Rothman
C                Cornell Theory Center/Cornell National Supercomputer
C                Facility.
C                e-mail address: BITNET:   eer@cornellf
C                                INTERNET: eer@cornellf.tn.cornell.edu
C
C
C   REFERENCE
C   ---------
C   1) Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
C      "A Testing Scheme for Subroutines Solving Large Linear Problems",
C       Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
C   2) Osterby, Ole and Zletev, Zahari;
C      "Direct Methods for Sparse Matrices";
C       Springer-Verlag 1983.
C
C
C
C   INPUT PARAMETERS
C   ----------------
C   N    - Integer. The size of the square matrix.
C          N > 2 must be specified.
C
C   NN   - Integer. The dimension of integer arrays IA and JA and 
C          real array AR. Must be at least NE.
C
C   NN  - Integer. The dimension of integer array JA. Must be at least
C          NE.
C
C   IC   - Integer. The sparsity pattern can be changed by means of this
C          parameter.  1 < IC < N   must be specified.
C
C
C
C   OUTPUT PARAMETERS
C   -----------------
C   NE   - Integer. The number of nonzero elements in the sparse matrix
C          of the type E(N,C). NE = 5*N - 2*IC - 2 . 
C
C   AR(NN)  - Real array.
C             Stored entries of the sparse matrix A.
C             NE is the number of nonzeros including a mandatory
C             diagonal entry for each row.
C
C   IA(NN)  - Integer array.(Double precision)
C             Pointers to specify rows for the stored nonzero entries
C             in AR.
C
C   JA(NN) - Integer array.
C             Pointers to specify columns for the stored nonzero entries
C             in AR.
C
C   IERR    - Error parameter is returned as zero on successful
C             execution of the subroutine.
C             Error diagnostics are given by means of positive values
C             of this parameter as follows:
C             IERR = 1    -  N       is out of range.
C             IERR = 2    -  IC      is out of range.
C             IERR = 3    -  NN      is out of range.
C
C---------------------------------------------------------------------
C
C
      real*8 ar(nn)
      integer ia(nn), ja(nn), n, ne, ierr
      ierr = 0
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
Check the input parameters:
c
      if(n.le.2)then
         ierr = 1
         return
      endif
      if(ic.le.1.or.ic.ge.n)then
         ierr = 2
         return
      endif
c
      ne = 5*n-2*ic-2 
      if(nn.lt.ne)then
         ierr = 3
         return
      endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c Begin to generate the nonzero elements as well as the row and column
c pointers:
c
      do 20 i=1,n
      ar(i) = 4.0d0
      ia(i) = i
      ja(i) = i
20    continue
      ilast = n
      do 30 i=1,n-1
      it = ilast + i
      ar(it) = -1.0d0
      ia(it) = i+1
      ja(it) = i
30    continue
      ilast = ilast + n - 1
      do 40 i=1,n-1
      it = ilast + i
      ar(it) = -1.0d0
      ia(it) = i
      ja(it) = i+1
40    continue
      ilast = ilast + n-1
      do 50 i=1,n-ic
      it = ilast + i
      ar(it) = -1.0d0
      ia(it) = i+ic
      ja(it) = i
50    continue
      ilast = ilast + n-ic
      do 60 I=1,n-ic
      it = ilast + i
      ar(it) = -1.0d0
      ia(it) = i
      ja(it) = i+ic
60    continue
c      ilast = ilast + n-ic
c      if(ilast.ne.5*n-2*ic-2) then
c      write(*,*)' ilast equal to ', ilast
c      write(*,*)' ILAST, the no. of nonzeros, should = ', 5*n-2*ic-2
c      stop
c      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end

