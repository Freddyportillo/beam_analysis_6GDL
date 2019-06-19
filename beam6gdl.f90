program beam6gdl
use K_M_matrix
use lapack_parcer
    implicit none

integer :: ii, jj, n, k, nb_ele, nb_nodes, ngdl_ele, ngdl_global, ngdl_node, status, n_modos
integer, allocatable, dimension (:,:) :: mat_conect, cc, fa, ident
double precision :: lengBEAM, lengELE, I, r
double precision, dimension(1,5) :: prop
double precision, allocatable, dimension (:,:) ::K_global, M_global, eigvet, K_ll, M_ll, K_ll_mat, modos
double precision, allocatable, dimension (:) :: eigval, freqs, FN, K_ll_diag, freq
integer, allocatable, dimension (:) :: gdl_livres
double precision :: pi=3.14159265359


!prop = [E I Rho ni A(m^2)]


prop(1,:) = (/2.10*10**11, 4.166/10**9, 7.8*10**3, 0.29, 0.0005/)

lengBEAM = 2.0
!DEFINA O NO. DE ELEMENTOS
nb_ele = 5
nb_nodes = nb_ele+1
ngdl_node = 6
ngdl_ele = ngdl_node*2
ngdl_global = nb_nodes*ngdl_node
lengELE = lengBEAM/(nb_ele)

!condições de contornos para uma viga monoengastada
allocate (cc(ngdl_node,2))
do ii = 1,ngdl_node
    cc(ii,1) = ii
    cc(ii,2) = 0.0d0
enddo

!aplicação das forças
allocate (fa(ngdl_global-6,2))
do ii = ngdl_node+1,ngdl_global
    fa(ii-6,1) = ii
    fa(ii-6,2) = 0.0d0
enddo

!Gerando a matriz conectividade
allocate (mat_conect(nb_ele, ngdl_ele), STAT = status)
k = 1
do ii = 1,nb_ele
    do jj = 1,ngdl_ele
        mat_conect(ii,jj) = k
        k = k+1
    enddo
    k=k-6
enddo

!Gerando as Matrizes elementares

allocate (K_global(ngdl_global,ngdl_global), STAT = status)
allocate (M_global(ngdl_global,ngdl_global), STAT = status)

K_global = 0.0d0

call KM_matrix (prop,mat_conect,lengELE,nb_ele,K_global,M_global)

print*, '-------------------------------------'

n = ngdl_global-6
allocate(gdl_livres(n), STAT = status)
gdl_livres = fa(:,1)

allocate (eigval(n), eigvet(n,n), freqs(10), freq(n), STAT = status)
allocate (K_ll(n,n), M_ll(n,n), STAT = status)

K_ll = K_global(gdl_livres(1):gdl_livres(n),gdl_livres(1):gdl_livres(n))
M_ll = M_global(gdl_livres(1):gdl_livres(n),gdl_livres(1):gdl_livres(n))

call lapack_eig (n,K_ll,M_ll,EigVal,EigVet)

OPEN (UNIT=11, FILE='K_ll.dat', STATUS='replace')
OPEN (UNIT=10, FILE='M_ll.dat', STATUS='replace')
do ii=1,n
    write(11,*) K_ll(ii,:)
    write(10,*) M_ll(ii,:)
end do

print*, 'As frequencias naturais sao (Hz): '
freq = sqrt(abs(eigval))/(2.0d0*pi)
n_modos = 5
print*,'------------'


do ii = 1, 10
    freqs(ii) = minval(freq)
    freq(minloc(freq)) = 999999.0d0
end do

print*, freqs

allocate(modos(ngdl_global-6,n_modos))
modos = eigvet(:,1:n_modos)

OPEN (UNIT=12, FILE='modos.dat', STATUS='replace')
OPEN (UNIT=13, FILE='modos_exp.dat', STATUS='replace')

do ii=1,n
    write(12,*) modos(ii,:)
    write(13,*) eigvet(ii,:)

end do

end program