*****************************************************************************
* Subroutines required for the earthquake inversion project
*****************************************************************************
* In unix:   R CMD SHLIB xxx.f    creates xxx.so
*            nm -g xxx.so         to see what objects are exported
* All routines must return values via their arguments: 
* i.e. use SUBROUTINEs only
*
* Use only:  INTEGER          as.integer()
*            REAL*8           as.double()
*            CHARACTER*255    as.character()
*
* NB: Fortran and R both store arrays in column major order: i.e with the
*     leftmost index varying fastest
*
*****************************************************************************
* Routine functions
*****************************************************************************
*-----------------------------------------------------------------------------
      SUBROUTINE SFLOOR(X,IVAL)

* Largest integer less than X

      IMPLICIT NONE
* Input variables
      REAL*8  X           ! value
* Ouput variables
      INTEGER IVAL        ! FLOOR(X)

      
      IVAL = DINT(X)
      IF(X.LT.0.AND.X.NE.IVAL) IVAL = IVAL-1

      RETURN
      END

*****************************************************************************
* Vector and matrix operations
*****************************************************************************
*-----------------------------------------------------------------------------
      SUBROUTINE VECTOR_COPY(AVEC1,N,AVEC2)

* Copy vector AVEC1 into AVEC2

      IMPLICIT NONE
* Input variables
      INTEGER N                      ! vector dimension
      REAL*8  AVEC1(N)               ! vector to copy
* Output variables
      REAL*8  AVEC2(N)               ! output vector
* Local variables
      INTEGER I                      ! loop control

      DO I = 1,N
         AVEC2(I) = AVEC1(I)
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE MATRIX_COPY(AMAT1,N1,N2,AMAT2)

* Copy matrix AMAT1 into AMAT2

      IMPLICIT NONE
* Input variables
      INTEGER N1, N2                 ! matrix dimensions
      REAL*8  AMAT1(N1,N2)           ! matrix to copy
* Output variables
      REAL*8  AMAT2(N1,N2)           ! output matrix
* Local variables
      INTEGER I, J                   ! loop control

      DO I = 1,N1
         DO J = 1,N2
            AMAT2(I,J) = AMAT1(I,J)
         END DO
      END DO

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE SCALAR_PRODUCT(A,B,N,ADOTB)

* Scalar product of N-vectors A and B

      IMPLICIT NONE
* Input variables
      INTEGER N                 ! vector dimension
      REAL*8  A(N), B(N)        ! input vectors
* Output variables
      REAL*8  ADOTB             ! output value
* Local variables
      INTEGER I                 ! loop control

      ADOTB = 0.0D0
      DO I = 1,N
         ADOTB = ADOTB + A(I)*B(I)
      END DO

      RETURN
      END


*-----------------------------------------------------------------------------
      SUBROUTINE MAT_TRANSPOSE(MAT,N,M,TMAT)

* Transpose a NxM matrix

      IMPLICIT NONE
* Input variables
      INTEGER N, M        ! matrix dimensions
      REAL*8  MAT(N,M)    ! matrix
* Output variables
      REAL*8  TMAT(M,N)   ! transposed matrix
* Local variables
      INTEGER I, J        ! loop control

      DO I = 1,M
         DO J = 1,N
            TMAT(I,J) = MAT(J,I)
         END DO
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE SQMAT_TRANSPOSE_INPLACE(MAT,N)

* Transpose a square NxN matrix - replace the input with the output

      IMPLICIT NONE
* Input variables
      INTEGER N           ! matrix dimension
* Modified input variables
      REAL*8  MAT(N,N)    ! matrix
* Local variables
      INTEGER I, J        ! loop control
      REAL*8  TEMP        ! temporary storage

      DO I = 2,N
        DO J = 1,I-1
           TEMP = MAT(I,J)
           MAT(I,J) = MAT(J,I)
           MAT(J,I) = TEMP
        END DO
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE MATRIX_MULTIPLY(MAT1,MAT2,N1,N12,N2,OUTMAT)

* Multiply: OUTMAT = MAT1 %*% MAT2

      IMPLICIT NONE
* Input variables
      INTEGER N1, N12, N2       ! array dimensions
      REAL*8  MAT1(N1,N12),     ! input arrays
     +        MAT2(N12,N2)      !
* Output variables
      REAL*8  OUTMAT(N1,N2)     ! output array
* Local variables
      INTEGER I, J, K           ! loop control
      REAL*8  TEMP              ! sum storage
      
      DO I = 1,N1
         DO J = 1,N2
            TEMP = 0.0
            DO K = 1,N12
               TEMP = TEMP + MAT1(I,K)*MAT2(K,J)
            END DO
            OUTMAT(I,J) = TEMP
         END DO
      END DO

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE ROWSWAP(RMAT,J0)

* Reverse row J0, and swap the other two rows
* Do nothing if J0=0

* Input variables
      INTEGER J0         ! row to reverse
      REAL*8  RMAT(3,3)  ! matrix to work on
* Local variables
      INTEGER I          ! loop control
      INTEGER J1, J2     ! other columns
      REAL*8  TEMP       ! temporary storage

      IF(J0.EQ.0) RETURN ! do nothing if J0=0

      J1 = J0-1
      J2 = J0+1
      IF(J1.LT.1) J1 = J1 + 3
      IF(J2.GT.3) J2 = J2 - 3

      DO I = 1,3
         TEMP = RMAT(J1,I)
         RMAT(J1,I) = RMAT(J2,I)
         RMAT(J2,I) = TEMP
         RMAT(J0,I) =-RMAT(J0,I)
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE COLSWAP(RMAT,J0)

* Reverse column J0, and swap the other two columns
* Do nothing if J0=0

* Input variables
      INTEGER J0         ! column to reverse
      REAL*8  RMAT(3,3)  ! matrix to work on
* Local variables
      INTEGER I          ! loop control
      INTEGER J1, J2     ! other columns
      REAL*8  TEMP       ! temporary storage

      IF(J0.EQ.0) RETURN ! do nothing if J0=0

      J1 = J0-1
      J2 = J0+1
      IF(J1.LT.1) J1 = J1 + 3
      IF(J2.GT.3) J2 = J2 - 3

      DO I = 1,3
         TEMP = RMAT(I,J1)
         RMAT(I,J1) = RMAT(I,J2)
         RMAT(I,J2) = TEMP
         RMAT(I,J0) =-RMAT(I,J0)
      END DO

      RETURN
      END

*****************************************************************************
* Operations with Euler angles
*****************************************************************************
*-----------------------------------------------------------------------------
      SUBROUTINE NVEC_PHIV2(PHIV2,NVEC)
     
* Unit vector NVEC corresponding to the Euler angles PHIV2=(phi,theta)

      IMPLICIT NONE
* Input variables
      REAL*8  PHIV2(2)          ! Euler angles
* Output variables
      REAL*8  NVEC(3)           ! Unit vector
* Local variables
      REAL*8  ST                ! sin(theta)

      ST = DSIN(PHIV2(2))
      NVEC(1) = ST*DCOS(PHIV2(1))
      NVEC(2) = ST*DSIN(PHIV2(1))
      NVEC(3) = DCOS(PHIV2(2))

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE EXTRACT_PHIV2_NVEC(NVEC,PHIV2)
     
* Get the azimuth (phi) and colatitude=polar distance (theta) of a unit vector

      IMPLICIT NONE
* Input variables
      REAL*8  PHIV2(2)          ! Euler angles
* Output variables
      REAL*8  NVEC(3)           ! Unit vector

      PHIV2(1) = DATAN2( NVEC(2), NVEC(1) )  ! Phi
      PHIV2(2) = DACOS(NVEC(3))              ! Theta

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE EXCHANGE_AXIS(AMAT,J)

* Exchange about axis J: negate column J, interchange the other two columns

      IMPLICIT NONE
* Input variables
      REAL*8  AMAT(3,3)              ! matrix to operate on
      INTEGER J                      ! axis
* Output variables
* Local variables
      INTEGER J1, J2                 ! other axes
      INTEGER I                      ! loop control
      REAL*8  TEMP                   ! temporary storage

      J1 = J-1
      J2 = J+1
      IF(J1.LT.1) J1 = 3
      IF(J2.GT.3) J2 = 1
      DO I = 1,3
         AMAT(I,J) = -AMAT(I,J)
         TEMP = AMAT(I,J1)
         AMAT(I,J1) = AMAT(I,J2)
         AMAT(I,J2) = TEMP
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE ROTMAT_PHIV(PHIV,RMAT)

* Calculate the rotation matrix RMAT for Euler angles PHIV

      IMPLICIT NONE
* Input variables
      REAL*8  PHIV(3)           ! Rotation angles
* Output variables
      REAL*8  RMAT(3,3)         ! Rotation matrix
* Local variables
      REAL*8  CP, SP            ! cos and sin phi
      REAL*8  CT, ST            ! cos and sin theta
      REAL*8  CS, SS            ! cos and sin psi

      CP = DCOS(PHIV(1))
      SP = DSIN(PHIV(1))
      CT = DCOS(PHIV(2))
      ST = DSIN(PHIV(2))
      CS = DCOS(PHIV(3))
      SS = DSIN(PHIV(3))

      RMAT(1,1) =  CP*CT*CS-SP*SS
      RMAT(1,2) = -CP*CT*SS-SP*CS
      RMAT(1,3) =  CP*ST

      RMAT(2,1) =  SP*CT*CS+CP*SS
      RMAT(2,2) = -SP*CT*SS+CP*CS
      RMAT(2,3) =  SP*ST

      RMAT(3,1) = -ST*CS
      RMAT(3,2) =  ST*SS
      RMAT(3,3) =  CT

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE EXTRACT_PHIV_ROTMAT(RMAT,PHIV)

* Calculate the rotation matrix RMAT for Euler angles PHIV

      IMPLICIT NONE
* Input variables
      REAL*8  RMAT(3,3)         ! Rotation matrix
* Output variables
      REAL*8  PHIV(3)           ! Rotation angles
* Local variables
* Parameters
      REAL*8  PI
      PARAMETER(PI=3.141592653589793D0)

      PHIV(2) = DACOS( RMAT(3,3) ) ! Theta
      IF(ABS(RMAT(3,3)).EQ.1)THEN
         PHIV(1) = DATAN2(-RMAT(1,2), RMAT(2,2) ) ! Phi
         PHIV(3) = 0.D0                           ! Psi
      ELSE
         PHIV(1) = DATAN2( RMAT(2,3), RMAT(1,3) ) ! Phi
         PHIV(3) = DATAN2( RMAT(3,2),-RMAT(3,1) ) ! Psi
      END IF

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE FLIP_PHIV(NPHIV,PHIV,PHIV_ALT)

* Convert PHIV into its alternative representation
* swap the x and z axes, and reverse the y axis

      IMPLICIT NONE
* Input variables
      INTEGER NPHIV             ! Number of observations
      REAL*8  PHIV(3,NPHIV)     ! Original data
* Output variables
      REAL*8  PHIV_ALT(3,NPHIV) ! Alternative data
* Local variables
      INTEGER I, J              ! Loop control
      REAL*8  RMAT1(3,3)        ! Original rotation matrix
      REAL*8  RMAT2(3,3)        ! Flipped rotation matrix

      DO I = 1,NPHIV
          CALL ROTMAT_PHIV(PHIV(1,I),RMAT1)
          DO J = 1,3
             RMAT2(J,1) =  RMAT1(J,3)
             RMAT2(J,2) = -RMAT1(J,2)
             RMAT2(J,3) =  RMAT1(J,1)
          END DO
          CALL EXTRACT_PHIV_ROTMAT(RMAT2,PHIV_ALT(1,I))
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE ROTANGLE_RMAT(RMAT,RANGLE)

* Get the rotation angle of rotation matrix RMAT

      IMPLICIT NONE
* Input variables
      REAL*8  RMAT(3,3)         ! Rotation matrix
* Output variables
      REAL*8  RANGLE            ! Rotation angle
* Local variables
      INTEGER I                 ! Loop control
      REAL*8  TT, T1            ! Temporary storage

* Trace
      TT = 0.0D0
      DO I = 1,3
         TT = TT + RMAT(I,I)
      END DO
      T1 = DMAX1(-2.0D0, DMIN1(2.0D0, TT-1.0D0))
* Rotation angle
      RANGLE = DATAN2( SQRT(4.0D0-T1*T1), T1 )

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE FLIP_SDR(NSDR,SDR,SDR_ALT)

* Convert SDR into its alternative representation

      IMPLICIT NONE
* Input variables
      INTEGER NSDR              ! Number of observations
      REAL*8  SDR(3,NSDR)       ! Original data
* Output variables
      REAL*8  SDR_ALT(3,NSDR)   ! Alternative data
* Local variables
      REAL*8  PHIV(3)           ! Original data as PHIV
      REAL*8  PHIV_ALT(3)       ! Converted data
      LOGICAL RESTRICT_DIP      ! Should the DIP be restricted?
      INTEGER I                 ! Loop control
* Data statements
      DATA    RESTRICT_DIP  /.TRUE./

      DO I = 1,NSDR
         CALL CONVERT_SDR_PHIV(1,SDR(1,I),PHIV)
         CALL FLIP_PHIV(1,PHIV,PHIV_ALT)
         CALL CONVERT_PHIV_SDR(1,PHIV_ALT,SDR_ALT(1,I),RESTRICT_DIP)
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE CONVERT_SDR_PHIV(NSDR,SDR,PHIV)

* Convert SDR into PHIV

      IMPLICIT NONE
* Input variables
      INTEGER NSDR              ! Number of observations
      REAL*8  SDR(3,NSDR)       ! Original data as SDR
* Output variables
      REAL*8  PHIV(3,NSDR)      ! Data converted to PHIV
* Local variables
      INTEGER I, J              ! Loop control
      INTEGER ITEMP             ! Temporary storage
* Parameters
      REAL*8  PI
      PARAMETER(PI=3.141592653589793D0)

      DO I = 1,NSDR
         PHIV(1,I) = SDR(1,I) + PI/2.
         PHIV(2,I) = PI - SDR(2,I)
         PHIV(3,I) = SDR(3,I) - PI/2.
* Put PHIV in the range ([0,2*pi], [0,pi], [0,2*pi])
         DO J = 1,3
            CALL SFLOOR(PHIV(J,I)/(2.*PI),ITEMP)
            PHIV(J,I) = PHIV(J,I) - 2*PI*ITEMP
         END DO
         IF(PHIV(2,I).GT.PI)THEN 
            PHIV(2,I) = 2*PI - PHIV(2,I)
            PHIV(1,I) = PHIV(1,I) + PI
            IF(PHIV(1,I).GE.2.*PI) PHIV(1,I) = PHIV(1,I) - 2.*PI
            PHIV(3,I) = PHIV(3,I) + PI
            IF(PHIV(3,I).GE.2.*PI) PHIV(3,I) = PHIV(3,I) - 2.*PI
         END IF
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE CONVERT_PHIV_SDR(NPHIV,PHIV,SDR,RESTRICT_DIP)

* Convert PHIV into SDR

      IMPLICIT NONE
* Input variables
      INTEGER NPHIV             ! Number of observations
      REAL*8  PHIV(3,NPHIV)     ! Original data as PHIV
      LOGICAL RESTRICT_DIP      ! Should the range of dip be restricted?
* Output variables
      REAL*8  SDR(3,NPHIV)      ! Data converted to SDR
* Local variables
      INTEGER I, J              ! Loop control
      INTEGER ITEMP             ! Temporary storage
* Parameters
      REAL*8  PI
      PARAMETER(PI=3.141592653589793D0)

      DO I = 1,NPHIV
         SDR(1,I) = PHIV(1,I) - PI/2.
         SDR(2,I) = PI - PHIV(2,I)
         SDR(3,I) = PHIV(3,I) + PI/2.
* Put SDR in the ranges ([0,2*pi], [0,pi], [0,2*pi])
         DO J = 1,3
            CALL SFLOOR(SDR(J,I)/(2.*PI),ITEMP)
            SDR(J,I) = SDR(J,I) - 2*PI*ITEMP
         END DO
         IF(SDR(2,I).GT.PI)THEN
            SDR(2,I) = 2*PI - SDR(2,I)
            SDR(1,I) = SDR(1,I) + PI
            IF(SDR(1,I).GE.2.*PI) SDR(1,I) = SDR(1,I) - 2.*PI
            SDR(3,I) = SDR(3,I) + PI
            IF(SDR(3,I).GE.2.*PI) SDR(3,I) = SDR(3,I) - 2.*PI
         END IF
         IF(RESTRICT_DIP.AND.SDR(2,I).GT.PI/2.)THEN
* restrict the (SDR) vector so that the dip is in the range [0,pi/2]
            SDR(1,I) = SDR(1,I) + PI
            IF(SDR(1,I).GE.2.*PI) SDR(1,I) = SDR(1,I) - 2.*PI
            SDR(2,I) = PI - SDR(2,I)
            SDR(3,I) = 2*PI - SDR(3,I)
         END IF
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE ROTMAT_SDR(SDR,RMAT)

* Calculate the rotation matrix RMAT for SDR

      IMPLICIT NONE
* Input variables
      REAL*8  SDR(3)            ! SDR angles
* Local variables
      REAL*8  PHIV(3)           ! Euler angles
* Output variables
      REAL*8  RMAT(3,3)         ! Rotation matrix

      CALL CONVERT_SDR_PHIV(1,SDR,PHIV)
      CALL ROTMAT_PHIV(PHIV,RMAT)

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE RESTRICT_PHIV_THETA(PHIV)

* Restrict a set of angles (phiv)=(phi,theta,psi)
* so that theta is in the range [0,pi/2]

      IMPLICIT NONE
* Input/Output variables
      REAL*8  PHIV(3)    ! Euler angles
* Local variables
      INTEGER I          ! Loop control
      INTEGER ITEMP      ! Output from SFLOOR()
* Parameters
      REAL*8  PI
      PARAMETER(PI=3.141592653589793D0)

* First put all angles in the range [0,2*pi]
      DO I = 1,3
         CALL SFLOOR(PHIV(I)/(2*PI),ITEMP)
         PHIV(I) = PHIV(I) - 2*PI*ITEMP
      END DO
* Check for theta value in the range [pi,2*pi]
      IF(PHIV(2).GT.PI)THEN
         ! add pi to phi
         PHIV(1) = PI + PHIV(1)
         CALL SFLOOR(PHIV(I)/(2*PI),ITEMP)
         PHIV(1) = PHIV(I) - 2*PI*ITEMP
         PHIV(2) = 2*PI - PHIV(2)
      END IF
* Now check for theta in the range [pi/2,pi]
      IF(PHIV(2).GT.PI/2.)THEN
         PHIV(1) = PI + PHIV(1)
         PHIV(2) = PI - PHIV(2)
         PHIV(3) = 2*PI - PHIV(3)
      END IF
* Finally check for psi greater than pi
*      IF(PHIV(3).GT.PI)THEN
*         PHIV(3) = PHIV(3) - PI
*      END IF
* And once more check that all angles are in [0,2*pi]
      DO I = 1,3
         CALL SFLOOR(PHIV(I)/(2*PI),ITEMP)
         PHIV(I) = PHIV(I) - 2*PI*ITEMP
      END DO

      RETURN
      END


*-----------------------------------------------------------------------------
      SUBROUTINE INTERP1D(X,F,NX,XMIN,DX,FVAL)

* Interpolation on a 1D regular grid XMIN + (0:(NX-1))*DX

      IMPLICIT NONE
* Input variables
      REAL*8  X        ! location at which interpolation is to take place
      INTEGER NX       ! number of grid points
      REAL*8  XMIN, DX ! min X value and stepsize
      REAL*8  F(NX)    ! array of function values
* Output variables
      REAL*8  FVAL     ! output function value
* Local variables
      INTEGER I1       ! interval number
      REAL*8  U        ! proportion of the interval

      I1 = MAX(1,MIN(NX-1, INT((X-XMIN)/DX)+1))
      U = (X-XMIN-(I1-1)*DX)/DX
      FVAL = F(I1)*(1-U) + F(I1+1)*U

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE INTERP2D(X,Y,F,NX,XMIN,DX,NY,YMIN,DY,FVAL)

* Interpolation on a 2D regular grid XMIN + (0:(NX-1))*DX
*                                    YMIN + (0:(NY-1))*DY

      IMPLICIT NONE
* Input variables
      REAL*8  X, Y     ! location at which interpolation is to take place
      INTEGER NX       ! number of grid points
      REAL*8  XMIN, DX ! min X value and stepsize
      INTEGER NY       ! number of grid points
      REAL*8  YMIN, DY ! min Y value and stepsize
      REAL*8  F(NX,NY) ! array of function values
* Output variables
      REAL*8  FVAL     ! output function value
* Local variables
      INTEGER I1, J1   ! interval number
      REAL*8  U, V     ! proportion of the interval

      I1 = MAX(1,MIN(NX-1, INT((X-XMIN)/DX)+1))
      J1 = MAX(1,MIN(NY-1, INT((Y-YMIN)/DY)+1))
      U = (X-XMIN-(I1-1)*DX)/DX
      V = (Y-YMIN-(J1-1)*DY)/DY
      FVAL =  F(I1,J1)*(1-U)*(1-V) + F(I1+1,J1)*U*(1-V) 
     +      + F(I1,J1+1)*(1-U)*V   + F(I1+1,J1+1)*U*V

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE INTERP3D(X,Y,Z,F,NX,XMIN,DX,NY,YMIN,DY,NZ,ZMIN,DZ,FVAL)

* Interpolation on a 3D regular grid XMIN + (0:(NX-1))*DX
*                                    YMIN + (0:(NY-1))*DY
*                                    ZMIN + (0:(NZ-1))*DZ

      IMPLICIT NONE
* Input variables
      REAL*8  X, Y, Z     ! location at which interpolation is to take place
      INTEGER NX          ! number of grid points
      REAL*8  XMIN, DX    ! min X value and stepsize
      INTEGER NY          ! number of grid points
      REAL*8  YMIN, DY    ! min Y value and stepsize
      INTEGER NZ          ! number of grid points
      REAL*8  ZMIN, DZ    ! min Z value and stepsize
      REAL*8  F(NX,NY,NZ) ! array of function values
* Output variables
      REAL*8  FVAL        ! output function value
* Local variables
      INTEGER I1, J1, K1  ! interval number
      REAL*8  U, V, W     ! proportion of the interval

      I1 = MAX(1,MIN(NX-1, INT((X-XMIN)/DX)+1))
      J1 = MAX(1,MIN(NY-1, INT((Y-YMIN)/DY)+1))
      K1 = MAX(1,MIN(NZ-1, INT((Z-ZMIN)/DZ)+1))
      U = (X-XMIN-(I1-1)*DX)/DX
      V = (Y-YMIN-(J1-1)*DY)/DY
      W = (Z-ZMIN-(K1-1)*DZ)/DZ
      FVAL =  F(I1,J1,K1)*(1-U)*(1-V)*(1-W) 
     +      + F(I1+1,J1,K1)*U*(1-V)*(1-W)
     +      + F(I1,J1+1,K1)*(1-U)*V*(1-W)   
     +      + F(I1+1,J1+1,K1)*U*V*(1-W)
     +      + F(I1,J1,K1+1)*(1-U)*(1-V)*W
     +      + F(I1+1,J1,K1+1)*U*(1-V)*W
     +      + F(I1,J1+1,K1+1)*(1-U)*V*W
     +      + F(I1+1,J1+1,K1+1)*U*V*W

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE INTERP3D_1(X,Y,Z,F,N0,I0,
     +                      NX,XMIN,DX,NY,YMIN,DY,NZ,ZMIN,DZ,FVAL)

* Interpolation on a 3D regular grid XMIN + (0:(NX-1))*DX
*                                    YMIN + (0:(NY-1))*DY
*                                    ZMIN + (0:(NZ-1))*DZ
* This version has an initial index on the table which is held constant

      IMPLICIT NONE
* Input variables
      REAL*8  X, Y, Z     ! location at which interpolation is to take place
      INTEGER N0, I0      ! grid size and location of first index
      INTEGER NX          ! number of grid points
      REAL*8  XMIN, DX    ! min X value and stepsize
      INTEGER NY          ! number of grid points
      REAL*8  YMIN, DY    ! min Y value and stepsize
      INTEGER NZ          ! number of grid points
      REAL*8  ZMIN, DZ    ! min Z value and stepsize
      REAL*8  F(N0,NX,NY,NZ) ! array of function values
* Output variables
      REAL*8  FVAL        ! output function value
* Local variables
      INTEGER I1, J1, K1  ! interval number
      REAL*8  U, V, W     ! proportion of the interval

      I1 = MAX(1,MIN(NX-1, INT((X-XMIN)/DX)+1))
      J1 = MAX(1,MIN(NY-1, INT((Y-YMIN)/DY)+1))
      K1 = MAX(1,MIN(NZ-1, INT((Z-ZMIN)/DZ)+1))
      U = (X-XMIN-(I1-1)*DX)/DX
      V = (Y-YMIN-(J1-1)*DY)/DY
      W = (Z-ZMIN-(K1-1)*DZ)/DZ
      FVAL =  F(I0,I1,J1,K1)*(1-U)*(1-V)*(1-W) 
     +      + F(I0,I1+1,J1,K1)*U*(1-V)*(1-W)
     +      + F(I0,I1,J1+1,K1)*(1-U)*V*(1-W)   
     +      + F(I0,I1+1,J1+1,K1)*U*V*(1-W)
     +      + F(I0,I1,J1,K1+1)*(1-U)*(1-V)*W
     +      + F(I0,I1+1,J1,K1+1)*U*(1-V)*W
     +      + F(I0,I1,J1+1,K1+1)*(1-U)*V*W
     +      + F(I0,I1+1,J1+1,K1+1)*U*V*W

      RETURN
      END

*****************************************************************************
* Functions for earthquake statistical model
*****************************************************************************
*-----------------------------------------------------------------------------
      SUBROUTINE WALLACE_BOTT_MAT_NVEC(NVEC,NU,PHIV_GS,RMAT_GS,WBMAT)

* Wallace-Bott Rotation Matrix R_SF given the fault normal NVEC
* and the earth model (NU,PHIV_GS,RMAT_GS)

      IMPLICIT NONE
* Input variables
      REAL*8  NVEC(3)           ! Fault normal
      REAL*8  NU                ! Triaxiality of stress tensor
      REAL*8  PHIV_GS(3)        ! Euler angles of stress tensor
      REAL*8  RMAT_GS(3,3)      ! Rotation matrix for stress tensor
* Output variables
      REAL*8  WBMAT(3,3)        ! Wallace Bott Matrix RSF for this fault normal
* Local variables
      INTEGER I                 ! Lopp control
      REAL*8  KAPPA             ! kappa of the unit vector + earth model
      REAL*8  LAMBDA            ! normalisation

      
      KAPPA = NVEC(1)*NVEC(1) + NU*NVEC(2)*NVEC(2)
      LAMBDA = 1/DSQRT( KAPPA*(1-KAPPA) - NU*(1-NU)*NVEC(2)*NVEC(2)  )
* Slip vector (1)
      WBMAT(1,1) =  LAMBDA* (1-KAPPA)*NVEC(1)
      WBMAT(2,1) = -LAMBDA*(KAPPA-NU)*NVEC(2)
      WBMAT(3,1) = -LAMBDA*     KAPPA*NVEC(3)
* Auxiliary vector (2)
      WBMAT(1,2) = -LAMBDA*    NU*NVEC(2)*NVEC(3)
      WBMAT(2,2) =  LAMBDA       *NVEC(1)*NVEC(3)
      WBMAT(3,2) = -LAMBDA*(1-NU)*NVEC(1)*NVEC(2)
* Normal vector (3)
      DO I = 1,3
         WBMAT(I,3) = NVEC(I)
      END DO

      RETURN
      END

*****************************************************************************
* Routines for integration and evaluation of the Q function
*****************************************************************************
*-----------------------------------------------------------------------------
      SUBROUTINE FTEST(NDIM,X,NUMFUN,FVAL)

* It must have parameters (NDIM,X,NUMFUN,FUNVLS)
* Ready for integration by ADBAYS()

* Test routine

      IMPLICIT NONE
* Input variables
      INTEGER NDIM               ! dimension of argument X
      INTEGER NUMFUN             ! dimension of output value FTEST()
      REAL*8  X(NDIM)            ! argument X
* Output variables
      REAL*8  FVAL(NUMFUN)       ! value

      FVAL(1) = DSIN(X(1))*X(2) - X(1)*X(1)

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE INTEGRATE_FTEST(NDIM, A, B, MINPTS, MAXPTS, 
     +                        EPSABS, EPSREL, 
     +                        RESULT, ABSERR, NEVAL, IFAIL, 
     +                        NWMAX, WORK)

* Integrate the test routine FTEST using the adaptive integrator ADBAYS

      IMPLICIT NONE
* Input variables
      INTEGER NDIM               ! dimension of the function FTEST
      REAL*8  A(NDIM), B(NDIM)   ! lower and upper bounds
      INTEGER MINPTS, MAXPTS     ! min and max number of points
      REAL*8  EPSABS, EPSREL     ! absolute and relative error requested
      INTEGER NWMAX              ! size of workspace provided
      REAL*8  WORK(NWMAX)        ! workspace
* Output variables
      REAL*8  RESULT(1)          ! result of integral
      REAL*8  ABSERR             ! achieved absolute error
      INTEGER NEVAL              ! number of function evaluations
      INTEGER IFAIL              ! error condition
* Local variables
      INTEGER  NUMFUN            ! number of functions to evaluate=1
      INTEGER  KEY               ! type of integration (set to 1)
      INTEGER  RESTAR            ! restart integral (set to 0)
      PARAMETER(NUMFUN=1, KEY=1, RESTAR=0)!
      INTEGER  NUM, MAXSUB       ! intermediate results
      INTEGER  NW                ! size of workspace required
* External declarations
      EXTERNAL FTEST             ! the function to be integrated

      IF(MINPTS.LT.0) MINPTS = 100
* Calculate size of workspace required
      IF(NDIM.EQ.1)THEN
         NUM = 15
      ELSE IF(KEY.EQ.0)THEN
         IF(NDIM.LT.12)THEN
            NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
         ELSE
            NUM = 1 + 2*NDIM*(NDIM+4)
         END IF
      ELSE IF(KEY.EQ.1)THEN 
         NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
      ELSE IF(KEY.EQ.2)THEN
         NUM = 1 + 4*NDIM + 6*NDIM*NDIM 
     +           + 4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
      ELSE 
         NUM = 1 + 2*NDIM*(NDIM+4)
      END IF
      MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1

      NW = MAXSUB*( 2*NDIM + 2*NUMFUN + 2 ) + 7*NUMFUN + NDIM
      IF(MAXPTS.LT.NUM) MAXPTS = NUM
      IF(MAXPTS.LT.MINPTS) MAXPTS = MINPTS

      CALL ADBAYS( NDIM, NUMFUN, A, B, MINPTS, MAXPTS, FTEST,
     +     EPSABS, EPSREL, KEY, NW, RESTAR, RESULT, ABSERR, NEVAL,
     +     IFAIL, WORK)

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE QINTEGRAND(NDIM,X,NUMFUN,FVAL,
     +                      FLIP,TAU,EPAR,HAS_RMAT_GS,RMAT_GS)

* It must have parameters (NDIM,X,NUMFUN,FUNVLS)

* Test routine

      IMPLICIT NONE
* Input variables
      INTEGER NDIM               ! dimension of argument X
      INTEGER NUMFUN             ! dimension of output value FTEST()
      REAL*8  X(NDIM)            ! argument X = PHICT
      INTEGER FLIP               ! flip the observation?
      REAL*8  TAU                ! precision of the observation
      REAL*8  EPAR(4)            ! (NU,PHIV_S(3))
      INTEGER HAS_RMAT_GS        ! is RMAT_GS calculated already?
      REAL*8  RMAT_GS(3,3)       ! rotation matrix E(PHIV_S)
* Output variables
      REAL*8  FVAL(NUMFUN)       ! value
* Local variables
      REAL*8  NU                 ! stress ratio
      REAL*8  PHIV_S(3)          ! orientation of the stress tensor
      REAL*8  PHIV_F(3)          ! focal mechanism in S coordinates
      REAL*8  SINPHI, COSPHI     ! temporary storage
      REAL*8  NUM, DEN           ! temporary storage
      REAL*8  RMAT_SF(3,3)       ! rotation matrix E(PHIV_F)
      REAL*8  TMAT(3,3)          ! temporary rotation matrix
      INTEGER I                  ! loop control

      NU = EPAR(1)
      DO I = 1,3
         PHIV_S(I) = EPAR(I+1)
      END DO
      SINPHI = DSIN(X(1))
      COSPHI = DCOS(X(1))
      IF(HAS_RMAT_GS.EQ.0) CALL ROTMAT_PHIV(PHIV_S,RMAT_GS)

      PHIV_F(1) = X(1)
      PHIV_F(2) = DACOS(X(2))
      NUM = -(1-NU)*SINPHI*COSPHI
      DEN = (COSPHI*COSPHI + NU*SINPHI*SINPHI)*X(2)
      PHIV_F(3) = ATAN2(NUM,DEN)
      CALL ROTMAT_PHIV(PHIV_F,RMAT_SF)
      CALL MATRIX_MULTIPLY(RMAT_GS,RMAT_SF,3,3,3,TMAT)

      IF(FLIP.EQ.1)THEN ! flip the observation
         FVAL(1) = TAU*(TMAT(3,1)-TMAT(2,2)+TMAT(1,3)-3) 
      ELSE 
         FVAL(1) = TAU*(TMAT(1,1)+TMAT(2,2)+TMAT(3,3)-3) 
      END IF
      FVAL(1) = DEXP(FVAL(1))

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE QINTEGRAND_ADBAYS(NDIM,X,NUMFUN,FVAL)

* Version of QINTEGRAND to be used with ADBAYS(): 
* the parameters TAU,EPAR,HAS_RMAT_GS,RMAT_GS 
* are passed through a common block

* Ready for integration by ADBAYS()

* Test routine

      IMPLICIT NONE
* Input variables
      INTEGER NDIM               ! dimension of argument X
      INTEGER NUMFUN             ! dimension of output value FTEST()
      REAL*8  X(NDIM)            ! argument X = PHICT
* Output variables
      REAL*8  FVAL(NUMFUN)       ! value
* Local variables
* Common block variables
      INTEGER C_FLIP             ! flip the observation?
      INTEGER C_HAS_RMAT_GS      ! is RMAT_GS calculated already?
      REAL*8  C_RMAT_GS(3,3)     ! rotation matrix E(PHIV_S)
      REAL*8  C_TAU              ! precision of the observation
      REAL*8  C_EPAR(4)          ! (NU,PHIV_S(3))
      COMMON/QINTEGRAND_ADBAYS_CMN/C_FLIP, C_HAS_RMAT_GS, C_RMAT_GS,
     +                             C_TAU, C_EPAR
      
      CALL QINTEGRAND(NDIM,X,NUMFUN,FVAL,
     +                C_FLIP,C_TAU,C_EPAR,C_HAS_RMAT_GS,C_RMAT_GS)

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE INTEGRATE_QINTEGRAND(NDIM, A, B, MINPTS, MAXPTS, 
     +                        EPSABS, EPSREL, 
     +                        RESULT, ABSERR, NEVAL, IFAIL, 
     +                        NWMAX, WORK,
     +                        FLIP,TAU,EPAR,HAS_RMAT_GS,RMAT_GS)

* Integrate the function QINTEGRAND using the adaptive integrator ADBAYS

      IMPLICIT NONE
* Input variables
      INTEGER NDIM               ! dimension of the function QINTEGRAND
      REAL*8  A(NDIM), B(NDIM)   ! lower and upper bounds
      INTEGER MINPTS, MAXPTS     ! min and max number of points
      REAL*8  EPSABS, EPSREL     ! absolute and relative error requested
      INTEGER NWMAX              ! size of workspace provided
      REAL*8  WORK(NWMAX)        ! workspace
      INTEGER FLIP               ! flip the observation?
      REAL*8  TAU                ! precision of the observation
      REAL*8  EPAR(4)            ! (NU,PHIV_S(3))
      INTEGER HAS_RMAT_GS        ! is RMAT_GS calculated already?
      REAL*8  RMAT_GS(3,3)       ! rotation matrix E(PHIV_S)
* Output variables
      REAL*8  RESULT(1)          ! result of integral
      REAL*8  ABSERR             ! achieved absolute error
      INTEGER NEVAL              ! number of function evaluations
      INTEGER IFAIL              ! error condition
* Local variables
      INTEGER  NUMFUN            ! number of functions to evaluate=1
      INTEGER  KEY               ! type of integration (set to 1)
      INTEGER  RESTAR            ! restart integral (set to 0)
      PARAMETER(NUMFUN=1, KEY=1, RESTAR=0)!
      INTEGER  NUM, MAXSUB       ! intermediate results
      INTEGER  NW                ! size of workspace required
      INTEGER  I, J              ! loop control
* External declarations
      EXTERNAL QINTEGRAND_ADBAYS ! the function to be integrated
* Common block variables
      INTEGER C_FLIP             ! flip the observation?
      INTEGER C_HAS_RMAT_GS      ! is RMAT_GS calculated already?
      REAL*8  C_RMAT_GS(3,3)     ! rotation matrix E(PHIV_S)
      REAL*8  C_TAU              ! precision of the observation
      REAL*8  C_EPAR(4)          ! (NU,PHIV_S(3))
      COMMON/QINTEGRAND_ADBAYS_CMN/C_FLIP, C_HAS_RMAT_GS, C_RMAT_GS,
     +                             C_TAU, C_EPAR

* Store the common block variables
      C_FLIP = FLIP
      C_TAU = TAU
      DO I = 1,3
         C_EPAR(I) = EPAR(I)
      END DO
      C_HAS_RMAT_GS = HAS_RMAT_GS
      DO I = 1,3
         DO J = 1,3
            C_RMAT_GS(I,J) = RMAT_GS(I,J)
         END DO
      END DO

      IF(MINPTS.LT.0) MINPTS = 100
* Calculate size of workspace required
      IF(NDIM.EQ.1)THEN
         NUM = 15
      ELSE IF(KEY.EQ.0)THEN
         IF(NDIM.LT.12)THEN
            NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
         ELSE
            NUM = 1 + 2*NDIM*(NDIM+4)
         END IF
      ELSE IF(KEY.EQ.1)THEN 
         NUM = 1 + 2*NDIM*(NDIM+3) + 2**NDIM
      ELSE IF(KEY.EQ.2)THEN
         NUM = 1 + 4*NDIM + 6*NDIM*NDIM 
     +           + 4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
      ELSE 
         NUM = 1 + 2*NDIM*(NDIM+4)
      END IF
      MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1

      NW = MAXSUB*( 2*NDIM + 2*NUMFUN + 2 ) + 7*NUMFUN + NDIM
      IF(MAXPTS.LT.NUM) MAXPTS = NUM
      IF(MAXPTS.LT.MINPTS) MAXPTS = MINPTS

      CALL ADBAYS( NDIM, NUMFUN, A, B, MINPTS, MAXPTS, 
     +     QINTEGRAND_ADBAYS,
     +     EPSABS, EPSREL, KEY, NW, RESTAR, RESULT, ABSERR, NEVAL,
     +     IFAIL, WORK)

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE QFUNC(N,EPAR,TAU,FLIP,FVAL,ABSERR,NEVAL,IFAIL,
     +                 EPSABS,EPSREL,NWMAX,WORK)

* Evaluate the function QFUNC using the adaptive integrator ADBAYS

      IMPLICIT NONE
* Input variables
      INTEGER N                  ! number of points to evaluate
      REAL*8  EPAR(4,N)          ! (NU,PHIV_S(3))
      REAL*8  TAU                ! precision of the observation      
      INTEGER FLIP               ! flip the observation?
      REAL*8  EPSABS, EPSREL     ! absolute and relative error requested
      INTEGER NWMAX              ! size of workspace provided
      REAL*8  WORK(NWMAX)        ! workspace
* Output variables
      REAL*8  FVAL(N)            ! result of integral
      REAL*8  ABSERR(N)          ! achieved absolute error
      INTEGER NEVAL(N)           ! number of function evaluations
      INTEGER IFAIL(N)           ! error condition
* Local variables
      INTEGER NDIM               ! dimension of the function QINTEGRAND
      PARAMETER(NDIM=2)
      REAL*8  A(NDIM), B(NDIM)   ! lower and upper bounds
      INTEGER MINPTS, MAXPTS     ! min and max number of points
      INTEGER HAS_RMAT_GS        ! is RMAT_GS calculated already?
      REAL*8  RMAT_GS(3,3)       ! rotation matrix E(PHIV_S)
      INTEGER  I                 ! loop control
* Parameters
      REAL*8  PI
      PARAMETER(PI=3.141592653589793D0)

      HAS_RMAT_GS = 1
      A(1) = 0D0   ! range of phi integration
      B(1) = 2*PI
      A(2) = -1D0  ! range of cos(theta) integration
      B(2) = +1D0
      MINPTS = 100
      MAXPTS = 10000
      DO I = 1,N
         CALL ROTMAT_PHIV(EPAR(2,I),RMAT_GS)
         CALL INTEGRATE_QINTEGRAND(NDIM, A, B, MINPTS, MAXPTS, 
     +                             EPSABS, EPSREL, 
     +                             FVAL(I), ABSERR(I), 
     +                             NEVAL(I), IFAIL(I), 
     +                             NWMAX, WORK,
     +                             FLIP,TAU,EPAR(1,I),
     +                             HAS_RMAT_GS,RMAT_GS)
      END DO

      RETURN
      END
*-----------------------------------------------------------------------------
      SUBROUTINE RETABULATE_QFUNC(PHIV_G1,FLIP,EXCHANGE,
     +                            QVALTAB1,QVALTAB0,
     +                            NNU, NUMIN, DNU,
     +                            NPHI, PHIMIN, DPHI,
     +                            NTHETA, THETAMIN, DTHETA,
     +                            NPSI, PSIMIN, DPSI)

* Evaluate the posterior distribution QVALTAB for an arbitrary observation 
* PHIV_G1 using the tabulation QKTAB0 for a single observation with
* PHIV_G=(0,0,0) with the same errors (TAU)

* To flip the observation round input PHIV_G1=(0,pi/2,pi) or
* specify FLIP=1

      IMPLICIT NONE
* Input variables
      REAL*8  PHIV_G1(3)                      ! observation of interest
      INTEGER FLIP                            ! flip the observation?
      INTEGER EXCHANGE                        ! =0(no change),1(swap 2-3),2(swap 1-3)
      INTEGER NNU, NPHI, NTHETA, NPSI         ! array sizes
      REAL*8  NUMIN, PHIMIN, THETAMIN, PSIMIN ! vector min values
      REAL*8  DNU, DPHI, DTHETA, DPSI         ! vector stepsizes
      REAL*8  QVALTAB0(NNU,NPHI,NTHETA,NPSI)  ! input table
* Output variables
      REAL*8  QVALTAB1(NNU,NPHI,NTHETA,NPSI)  ! new table
* Local variables
      INTEGER INU, IPHI, ITHETA, IPSI         ! loop control
      REAL*8  TRMAT_GF1(3,3)                  ! (R_GF1)^T
      REAL*8  PHIV0(3), PHIV1(3)              ! vector storage of coordinates
      REAL*8  EMAT0(3,3), EMAT1(3,3)          ! rotation matrices

* Rotation matrix of the observation
      CALL ROTMAT_PHIV(PHIV_G1,TRMAT_GF1)
      CALL SQMAT_TRANSPOSE_INPLACE(TRMAT_GF1,3)
      IF(FLIP.EQ.1)THEN
         ! swap the first and third rows, negate the second row
         CALL ROWSWAP(TRMAT_GF1,2)
      END IF

      DO INU = 1,NNU
         DO IPSI = 1,NPSI
            PHIV1(3) = PSIMIN + (IPSI-1)*DPSI
            DO ITHETA = 1,NTHETA
               PHIV1(2) = THETAMIN + (ITHETA-1)*DTHETA
               DO IPHI = 1,NPHI
                  PHIV1(1) = PHIMIN + (IPHI-1)*DPHI

* Rotation matrix of the current location EMAT1(PHIV1)
                  CALL ROTMAT_PHIV(PHIV1,EMAT1)
* Rotation matrix in existing tabulation EMAT0(PHIV2)
                  CALL MATRIX_MULTIPLY(TRMAT_GF1,EMAT1,3,3,3,EMAT0)
* Swap axes EXCHANGE=0 (no swap) =1 (swap y-z) =2 (swap x-z)
                  CALL COLSWAP(EMAT0,EXCHANGE)
* Extract the angles PHIV2 and ensure they all lie in the range (0,pi)
                  CALL EXTRACT_PHIV_ROTMAT(EMAT0,PHIV0)
                  CALL RESTRICT_PHIV_THETA(PHIV0)
* Interpolate
                  CALL INTERP3D_1(PHIV0(1),PHIV0(2),PHIV0(3),
     +                            QVALTAB0, NNU, INU,
     +                            NPHI, PHIMIN, DPHI,
     +                            NTHETA, THETAMIN, DTHETA,
     +                            NPSI, PSIMIN, DPSI,
     +                            QVALTAB1(INU,IPHI,ITHETA,IPSI))

               END DO
            END DO
         END DO
      END DO

      RETURN
      END

*-----------------------------------------------------------------------------
      SUBROUTINE RETABULATE_QPHIV_FUNC(PHIV_G1,FLIP,EXCHANGE,
     +                                 QPHIV_VALTAB1,QPHIV_VALTAB0,
     +                                 NPHI, PHIMIN, DPHI,
     +                                 NTHETA, THETAMIN, DTHETA,
     +                                 NPSI, PSIMIN, DPSI)

* Evaluate the posterior distribution QPHIV_VALTAB for an arbitrary observation 
* PHIV_G1 using the tabulation QPHIV_VALTAB0 for a single observation with
* PHIV_G=(0,0,0) with the same errors (TAU)

* This does the same job as RETABULATE_QFUNC - but for an 
* array without the NU dimension

* To flip the observation round input PHIV_G1=(0,pi/2,pi) or
* specify FLIP=1

      IMPLICIT NONE
* Input variables
      REAL*8  PHIV_G1(3)                      ! observation of interest
      INTEGER FLIP                            ! flip the observation?
      INTEGER EXCHANGE                        ! =0(no change),1(swap 2-3),2(swap 1-3)
      INTEGER NPHI, NTHETA, NPSI              ! array sizes
      REAL*8  PHIMIN, THETAMIN, PSIMIN        ! vector min values
      REAL*8  DPHI, DTHETA, DPSI              ! vector stepsizes
      REAL*8  QPHIV_VALTAB0(NPHI,NTHETA,NPSI) ! input table
* Output variables
      REAL*8  QPHIV_VALTAB1(NPHI,NTHETA,NPSI) ! new table
* Local variables
      INTEGER IPHI, ITHETA, IPSI              ! loop control
      REAL*8  TRMAT_GF1(3,3)                  ! (R_GF1)^T
      REAL*8  PHIV0(3), PHIV1(3)              ! vector storage of coordinates
      REAL*8  EMAT0(3,3), EMAT1(3,3)          ! rotation matrices

* Rotation matrix of the observation
      CALL ROTMAT_PHIV(PHIV_G1,TRMAT_GF1)
      CALL SQMAT_TRANSPOSE_INPLACE(TRMAT_GF1,3)
      IF(FLIP.EQ.1)THEN
         ! swap the first and third rows, negate the second row
         CALL ROWSWAP(TRMAT_GF1,2)
      END IF

      DO IPSI = 1,NPSI
         PHIV1(3) = PSIMIN + (IPSI-1)*DPSI
         DO ITHETA = 1,NTHETA
            PHIV1(2) = THETAMIN + (ITHETA-1)*DTHETA
            DO IPHI = 1,NPHI
               PHIV1(1) = PHIMIN + (IPHI-1)*DPHI

* Rotation matrix of the current location EMAT1(PHIV1)
               CALL ROTMAT_PHIV(PHIV1,EMAT1)
* Rotation matrix in existing tabulation EMAT0(PHIV2)
               CALL MATRIX_MULTIPLY(TRMAT_GF1,EMAT1,3,3,3,EMAT0)
* Swap axes EXCHANGE=0 (no swap) =1 (swap y-z) =2 (swap x-z)
               CALL COLSWAP(EMAT0,EXCHANGE)
* Extract the angles PHIV2 and ensure they all lie in the range (0,pi)
               CALL EXTRACT_PHIV_ROTMAT(EMAT0,PHIV0)
               CALL RESTRICT_PHIV_THETA(PHIV0)
* Interpolate
               CALL INTERP3D(PHIV0(1),PHIV0(2),PHIV0(3),
     +                       QPHIV_VALTAB0, 
     +                       NPHI, PHIMIN, DPHI,
     +                       NTHETA, THETAMIN, DTHETA,
     +                       NPSI, PSIMIN, DPSI,
     +                       QPHIV_VALTAB1(IPHI,ITHETA,IPSI))

            END DO
         END DO
      END DO

      RETURN
      END

*****************************************************************************
      SUBROUTINE ARRAY_COLUMN_BOXCAR(N1,N2,A,MMAX,S,WRAP,IDEND,W,XTMP)

* Apply a triangular boxcar filter to each column of matrix A(N1,N2)
* Construct weights using the vector S

      IMPLICIT NONE
* Input variables
      INTEGER N1, N2
      REAL*8  A(N1,N2)
      INTEGER MMAX
      INTEGER S(N2)
      INTEGER WRAP, IDEND
      REAL*8  W(MMAX)
      REAL*8  XTMP(N1)
* Local variables
      INTEGER I1, I2, IW, NW

      DO I2 = 1,N2
         IF(S(I2).GT.0)THEN
            NW = 2*S(I2)+1
            DO IW = 1,S(I2)
               W(IW) = DBLE(IW)
               W(NW+1-IW) = DBLE(IW)
            END DO
            W(S(I2)+1) = S(I2)+1.0D0
            CALL BOXCAR(N1,A(1,I2),NW,W,WRAP,IDEND,XTMP)
            DO I1 = 1,N1
               A(I1,I2) = XTMP(I1)
            END DO
         END IF
      END DO

      RETURN
      END

*****************************************************************************
      SUBROUTINE BOXCAR(N,X,NW,W,WRAP,IDEND,XOUT)

* Apply the boxcar filter with weights W() to X() - return as XOUT()
* If WRAP=1 then wrap the function - noting in addition that
* if IDEND=1 then the first and last elements of X() should be 
* identical
      
      IMPLICIT NONE
* Input variables
      INTEGER N
      REAL*8  X(N)
      INTEGER NW
      REAL*8  W(NW)
      INTEGER WRAP, IDEND
* Output variables
      REAL*8  XOUT(N)
* Local variables
      INTEGER I, IW, S, I1, N1
      REAL*8  TEMP

      S = INT(NW/2)
      IF(2*S.EQ.NW) WRITE(6,*) 'Error in length of W'
      IF(S.GT.INT(N/2)) WRITE(6,*) 'W too wide'
      IF(NW.EQ.1)THEN ! no smoothing required
         DO I = 1,N
            XOUT(I) = X(I)
         END DO
         RETURN
      END IF

      TEMP = 0.0D0
      DO IW = 1,NW
         TEMP = TEMP + W(IW)
         IF(W(IW).LT.0) WRITE(6,*) 'negative W() element found'
      END DO
      IF(TEMP.EQ.0) WRITE(6,*) 'weights sum to zero'
      DO IW = 1,NW
         W(IW) = W(IW)/TEMP
      END DO

      IF(WRAP.EQ.1.AND.IDEND.EQ.1)THEN
         N1 = N
      ELSE
         N1 = N-1
      END IF

      DO I = S+1,N1-S
         XOUT(I) = 0.0D0
         DO IW = 1,NW
            I1 = I-S-1+IW
            XOUT(I) = XOUT(I) + W(IW)*X(I1)
         END DO
      END DO

      IF(WRAP.EQ.1)THEN
         DO I = 1,S
            XOUT(I) = 0.0D0
            DO IW = 1,S+1-I
               I1 = N1-S+I-1+IW
               XOUT(I) = XOUT(I) + W(IW)*X(I1)
            END DO
            DO IW = S+2-I,NW
               I1 = I-S-1+IW
               XOUT(I) = XOUT(I) + W(IW)*X(I1)
            END DO
         END DO
         DO I = N1-S+1,N1
            XOUT(I) = 0.0D0
            DO IW = 1,N1-I+S+1
               I1 = I-S-1+IW
               XOUT(I) = XOUT(I) + W(IW)*X(I1)
            END DO
            DO IW = N1-I+S+2,NW
               I1 = -N1+I-S-1+IW
               XOUT(I) = XOUT(I) + W(IW)*X(I1)
            END DO
         END DO
         IF(IDEND.EQ.1) XOUT(N) = XOUT(1)
      ELSE
         DO I = 1,S
            XOUT(I) = 0.0D0
            DO IW = 1,S+1-I
               I1 = 1
               XOUT(I) = XOUT(I) + W(IW)*X(I1)
            END DO
            DO IW = S+2-I,NW
               I1 = I-S-1+IW
               XOUT(I) = XOUT(I) + W(IW)*X(I1)
            END DO
         END DO
         DO I = N1-S+1,N1
            XOUT(I) = 0.0D0
            DO IW = 1,N1-I+S+1
               I1 = I-S-1+IW
               XOUT(I) = XOUT(I) + W(IW)*X(I1)
            END DO
            DO IW = N1-I+S+2,NW
               I1 = N1
               XOUT(I) = XOUT(I) + W(IW)*X(I1)
            END DO
         END DO
      END IF

      RETURN
      END
*****************************************************************************
