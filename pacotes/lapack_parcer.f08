module lapack_parcer
    use lapack
    
    implicit none
    
    !#####################################
    !            Module Usage:
    !#####################################
    private
    
    !The following subroutines/variables are the only
    !ones accessible outside this module
    public:: lapack_inverse, lapack_lin_sys, lapack_eig
    !-------------------------------------
    
contains
    
     
!##############################################################################################################################   
!Calculate the Inverse of a square Matrix using Lapack routines
subroutine lapack_inverse(n,A,InvA)
    implicit none
    !Entrada:
    integer:: n
    double precision, allocatable, dimension (:,:):: A, InvA
    !Local:
    integer:: ierr
    integer, allocatable, dimension (:):: ww1i
    double precision, allocatable, dimension (:):: ww1d
    
    InvA = A
    !Allocate variables needed for Lapack's DGETRF and DGETRI subroutines
    allocate(ww1i(n),ww1d(n))
    ww1d=0.d0; ww1i=0
    
    ierr=0
    call DGETRF(n,n,InvA,n,ww1i,ierr)
    if (ierr/=0) then
        write(*,'(1x,a,i3)')'ERROR IN LU FACTORIZATION. ERROR TRACKER:',ierr
        stop
    end if
    ierr=0
    call DGETRI(n,InvA,n,ww1i,ww1d,n,ierr)
    if (ierr/=0) then
        write(*,'(1x,a,i3)')'ERROR IN INVERSE MATRIX. ERROR TRACKER:',ierr
        stop
    end if
    deallocate(ww1i,ww1d)
    
end subroutine lapack_inverse
!##############################################################################################################################   


!##############################################################################################################################   
!Solve the Linear System [A].{X} = {B} using Lapack routines
subroutine lapack_lin_sys(n,A,B,X)
    implicit none
    !Entrada:
    integer:: n
    double precision, allocatable, dimension (:,:):: A
    double precision, allocatable, dimension (:):: B, X
    !Local:
    integer:: ierr
    double precision, allocatable, dimension (:,:)::A2
    integer, allocatable, dimension (:):: ww1i
    
    !Allocate variables needed for Lapack's DGESV subroutine
    allocate(ww1i(n))
    ww1i=0; ierr=0
    !Allocate Matrices to keep A and B unchanged
    allocate(A2(n,n))
    A2=A
    X=B
    
    !Solve linear system
    call DGESV(n, 1, A2, n, ww1i, X, n, ierr )
    if (ierr/=0) then
        write(*,'(1x,a,i3)')'ERROR IN LINEAR SOLVER. ERROR TRACKER:',ierr
        stop
    end if
    
    deallocate(ww1i,A2)

end subroutine lapack_lin_sys
!##############################################################################################################################   


!##############################################################################################################################   
!Solve the Generalized Eigenvalue Problem [A].{v} = lambda.[B].{v} using Lapack routines
subroutine lapack_eig(n,A,B,EigVal,EigVet)
    implicit none
    !Entrada:
    integer:: n
    double precision, allocatable, dimension (:,:):: A, B, EigVet
    double precision, allocatable, dimension (:):: EigVal
    !Local:
    integer:: ierr
    double precision, allocatable, dimension (:,:):: A2, B2, EigVet_L
    double precision, allocatable, dimension (:):: ALPHAR, ALPHAI, BETAf, ww1d
    
    
    !Allocate variables needed for Lapack's DGGEV subroutine
    allocate(ALPHAR(n),ALPHAI(n),BETAf(n),ww1d(8*n),EigVet_L(n,n))
    ALPHAR=0.d0; ALPHAI=0.d0; BETAf=0.d0; ww1d=0.d0; EigVet_L = 0.d0
    
    !Allocate Matrices to keep A and B unchanged
    allocate(A2(n,n),B2(n,n))
    A2=A
    B2=B
    ierr=0
    
    !Solve EigenValue problem
    CALL DGGEV('N','V',n,A2,n,B2,n,ALPHAR,ALPHAI,BETAf,EigVet_L,n,EigVet,n,ww1d,8*n,ierr)
    
    !Verify Solver result to inform user if needed
    if (ierr/=0) then
        write(*,*)'ERROR IN EIGEN SOLVER ERROR TRACKER:',ierr
        stop
    end if
    
    !Calculate EigenValues
    EigVal = ALPHAR/BETAf
    
    deallocate(ALPHAR,ALPHAI,BETAf,ww1d,EigVet_L,A2,B2)
    
end subroutine lapack_eig
!##############################################################################################################################   

    
end module lapack_parcer