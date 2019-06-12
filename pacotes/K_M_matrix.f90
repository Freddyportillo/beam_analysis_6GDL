module K_M_matrix

    implicit none
    save
    public :: KM_matrix

    contains
  
    subroutine KM_matrix (prop,mat_conect,lengELE,nb_ele,K_global, M_global)
    
    integer :: ngdl_global, ngdl_ele, ii, status
    integer, intent(in), dimension (:,:) :: mat_conect
    integer, intent(in) :: nb_ele
    double precision, intent(in) :: lengELE
    integer, allocatable, dimension (:,:) :: ident, mat_transf, transmat
    double precision :: pi = 3.14159
    double precision :: E, I, S, a, aux1, aux2, aux3, kcc, L,Rho, ni, G, J, r
    double precision, dimension (12,12) :: Kc, Kf, M1, Ms
    double precision,allocatable, dimension (:,:) :: Kc_global, Kf_global, M1_global, Ms_global, Kc1, Kf1, M1_1, Ms1
    double precision, intent(out), allocatable, dimension (:,:) :: K_global, M_global
    double precision, intent(in), dimension (1,5) :: prop
    
    ngdl_global = (nb_ele+1)*6
    ngdl_ele = 12
    L = lengELE
    !prop = [E I Rho ni S]
    r = sqrt(S/pi)
    E = prop(1,1)
    I = prop(1,2)
    Rho = prop(1,3)
    ni = prop(1,4) 
    S = prop(1,5)   
    kcc = 6*(1+ni)/(7+6*ni)
    G = E/(3*(1+ni))
    J = pi/2*r**4
    a = 12*E*I/(G*kcc*S*L**2)

    aux1=E*I/((1+a)*L**3);
    aux2=G*J/L;
    aux3=E*S/L;
    !-------------------------------------------------------------------------------------------------------------------------
    !MATRIZ DE RIGIDEZ CLÁSSICA
    Kc(1,:) = (/12.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-6.0d0*L,-12.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-6.0d0*L/)       
    Kc(2,:) = (/0.0d0,aux3/aux1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-aux3/aux1,0.0d0,0.0d0,0.0d0,0.0d0/)        
    Kc(3,:) = (/0.0d0,0.0d0,12.0d0,6.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,-12.0d0,6.0d0*L,0.0d0,0.0d0/)         
    Kc(4,:) = (/0.0d0,0.0d0,6.0d0*L,(4.0d0+a)*L**2,0.0d0,0.0d0,0.0d0,0.0d0,-6.0d0*L,(2.0d0-a)*L**2,0.0d0,0.0d0/)          
    Kc(5,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,aux2/aux1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-aux2/aux1,0.0d0/)        
    Kc(6,:) = (/-6.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,(4.0d0+a)*L**2,6.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,(2.0d0-a)*L**2/)   
    Kc(7,:) = (/-12.0d0,0.0d0,0.0d0,0.0d0,0.0d0,6.0d0*L,12.0d0,0.0d0,0.0d0,0.0d0,0.0d0,6.0d0*L/)        
    Kc(8,:) = (/0.0d0,-aux3/aux1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,aux3/aux1,0.0d0,0.0d0,0.0d0,0.0d0/)          
    Kc(9,:) = (/0.0d0,0.0d0,-12.0d0,-6.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,12.0d0,-6.0d0*L,0.0d0,0.0d0/)           
    Kc(10,:) = (/0.0d0,0.0d0,6.0d0*L,(2.0d0-a)*L**2,0.0d0,0.0d0,0.0d0,0.0d0,-6.0d0*L,(4.0d0+a)*L**2,0.0d0,0.0d0/)          
    Kc(11,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,-aux2/aux1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,aux2/aux1,0.0d0/)        
    Kc(12,:) = (/-6.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,(2.0d0-a)*L**2,6.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,(4.0d0+a)*L**2/)

    Kc = aux1*Kc
    !-------------------------------------------------------------------------------------------------------------------------
    !MATRIZ DE RIGIDEZ GEOMÉTRICA
    Kf(1,:) = (/36.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0*L,-36.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0*L/)
    Kf(2,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)   
    Kf(3,:) = (/0.0d0,0.0d0,36.0d0,3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,-36.0d0,3.0d0*L,0.0d0,0.0d0/)   
    Kf(4,:) = (/0.0d0,0.0d0,3.0d0*L,4.0d0*L**2,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0*L,-L**2,0.0d0,0.0d0/)   
    Kf(5,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)   
    Kf(6,:) = (/-3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0*L**2,3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,-L**2/)
    Kf(7,:) = (/-36.0d0,0.0d0,0.0d0,0.0d0,0.0d0,3.0d0*L,36.0d0,0.0d0,0.0d0,0.0d0,0.0d0,3.0d0*L/)     
    Kf(8,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)      
    Kf(9,:) = (/0.0d0,0.0d0,-36.0d0,-3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,36.0d0,-3.0d0*L,0.0d0,0.0d0/)         
    Kf(10,:) = (/0.0d0,0.0d0,3.0d0*L,-L**2,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0*L,4.0d0*L**2,0.0d0,0.0d0/)       
    Kf(11,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)       
    Kf(12,:) = (/-3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,-L**2,3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0*L**2/)
  
    Kf = (1/(30*L))*Kf
    !-------------------------------------------------------------------------------------------------------------------------
    !MATRIZ DE MASSA
    aux1=Rho*S*L/420
    aux2=Rho*J*L/6
    M1(1,:) = (/156.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-22.0d0*L,54.0d0,0.0d0,0.0d0,0.0d0,0.0d0,13.0d0*L/)     
    M1(2,:) = (/0.0d0,140.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,70.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)         
    M1(3,:) = (/0.0d0,0.0d0,156.0d0,22.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,54.0d0,-13.0d0*L,0.0d0,0.0d0/)         
    M1(4,:) = (/0.0d0,0.0d0,22.0d0*L,4.0d0*L**2,0.0d0,0.0d0,0.0d0,0.0d0,13.0d0*L,-3.0d0*L**2,0.0d0,0.0d0/)          
    M1(5,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,2.0d0*aux2/aux1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,aux2/aux1/)         
    M1(6,:) = (/-22.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0*L**2,-13.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,-3*L**2/)    
    M1(7,:) = (/54.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-13.0d0*L,156.0d0,0.0d0,0.0d0,0.0d0,0.0d0,22.0d0*L/)        
    M1(8,:) = (/0.0d0,70.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,140.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)           
    M1(9,:) = (/0.0d0,0.0d0,54.0d0,13.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,156.0d0,-22.0d0*L,0.0d0,0.0d0/)          
    M1(10,:) =(/0.0d0,0.0d0,-13.0d0*L,-3.0d0*L**2,0.0d0,0.0d0,0.0d0,0.0d0,-22.0d0*L,4.0d0*L**2,0.0d0,0.0d0/)      
    M1(11,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,aux2/aux1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,2.0d0*aux2/aux1,0.0d0/)    
    M1(12,:) = (/13.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0*L**2,22.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0*L**2/)
    
    M1 = aux1*M1
    !-------------------------------------------------------------------------------------------------------------------------

    Ms(1,:) = (/36.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0*L,-36.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0*L/)     
    Ms(2,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)         
    Ms(3,:) = (/0.0d0,0.0d0,36.0d0,3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,-36.0d0,3.0d0*L,0.0d0,0.0d0/)        
    Ms(4,:) = (/0.0d0,0.0d0,3.0d0*L,4.0d0*L**2,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0*L,-L**2,0.0d0,0.0d0/)          
    Ms(5,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)          
    Ms(6,:) = (/-3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0*L**2,3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,-L**2/)      
    Ms(7,:) = (/-36.0d0,0.0d0,0.0d0,0.0d0,0.0d0,3.0d0*L,36.0d0,0.0d0,0.0d0,0.0d0,0.0d0,3.0d0*L/)     
    Ms(8,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)     
    Ms(9,:) = (/0.0d0,0.0d0,-36.0d0,-3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,36.0d0,-3.0d0*L,0.0d0,0.0d0/)         
    Ms(10,:) = (/0.0d0,0.0d0,3.0d0*L,-L**2,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0*L,4.0d0*L**2,0.0d0,0.0d0/)           
    Ms(11,:) = (/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/)         
    Ms(12,:) = (/-3.0d0*L,0.0d0,0.0d0,0.0d0,0.0d0,-L**2,3*L,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0*L**2/)
    
    Ms = (Rho*I/(30*L))*Ms
    !-------------------------------------------------------------------------------------------------------------------------
    
    
    !Matriz identidade
    allocate (ident(ngdl_global,ngdl_global), STAT = status)
    ident = 0.0d0
    do ii=1,ngdl_global
      ident(ii,ii) = 1
    end do
  
    allocate(mat_transf(ngdl_ele,ngdl_global), transmat(ngdl_global, ngdl_ele), Kc_global(ngdl_global,ngdl_global), STAT = status)
    allocate(Kf_global(ngdl_global,ngdl_global), K_global(ngdl_global,ngdl_global), STAT = status)
    allocate(M1_global(ngdl_global,ngdl_global), Ms_global(ngdl_global,ngdl_global), STAT = status)
    allocate(M_global(ngdl_global,ngdl_global), STAT = status)
    allocate(Kf1(ngdl_global,ngdl_global), Kc1(ngdl_global,ngdl_global), STAT = status)
    allocate(M1_1(ngdl_global,ngdl_global), Ms1(ngdl_global,ngdl_global), STAT = status)
    Kc_global = 0.0d0
    Kc_global = 0.0d0
    M1_global = 0.0d0
    Ms_global = 0.0d0
    K_global = 0.0d0
    M_global = 0.0d0
    
    do ii = 1, nb_ele
        mat_transf (1,:) = ident(mat_conect(ii,1),:)
        mat_transf (2,:) = ident(mat_conect(ii,2),:)
        mat_transf (3,:) = ident(mat_conect(ii,3),:)
        mat_transf (4,:) = ident(mat_conect(ii,4),:)
        mat_transf (5,:) = ident(mat_conect(ii,5),:)
        mat_transf (6,:) = ident(mat_conect(ii,6),:)
        mat_transf (7,:) = ident(mat_conect(ii,7),:)
        mat_transf (8,:) = ident(mat_conect(ii,8),:)
        mat_transf (9,:) = ident(mat_conect(ii,9),:)
        mat_transf (10,:) = ident(mat_conect(ii,10),:)
        mat_transf (11,:) = ident(mat_conect(ii,11),:)
        mat_transf (12,:) = ident(mat_conect(ii,12),:)
        transmat = transpose(mat_transf)
        call matmu(transmat,Kc,mat_transf,Kc1)
        call matmu(transmat,Kf,mat_transf,Kf1)
        call matmu(transmat,M1,mat_transf,M1_1)
        call matmu(transmat,Ms,mat_transf,Ms1)

        Kc_global = Kc_global + Kc1
        Kf_global = Kf_global + Kf1
        M1_global = M1_global + M1_1
        Ms_global = Ms_global + Ms1

    enddo   
      K_global = Kc_global + Kf_global
      M_global = M1_global + Ms_global

    deallocate(mat_transf, transmat, Kc_global, Kf_global, ident, STAT = status)
    

end subroutine KM_matrix

subroutine matmu (A,B,C,D)

integer, dimension (36,12) :: A
double precision, dimension (36,12) :: mul
double precision, dimension (12,12) :: B
integer, dimension (12,36) :: C
double precision, dimension (36,36) :: D

mul = matmul(A,B)
D = matmul (mul,C)

end subroutine matmu


end module K_M_matrix
