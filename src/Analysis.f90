MODULE Analysis
     IMPLICIT NONE
     PUBLIC Read_Input_File, Set_DOF_Number, Calculate_Index, Assemble_Kt, Assemble_Load
     PUBLIC Solve_Equations, Displacement_Stress, Print_Matrix
     PUBLIC COLSOL	
	
     CONTAINS	
! ---------------------------------------------------------------------------------------	
! read input file
      subroutine Read_Input_File()
         real(8) :: thickness
         character :: bufs
         integer :: bufi, i

         ! read nodal information
         read(1,*) bufs; read(1,*) node_n; read(1,*) bufs
         allocate(Node(node_n)) ! allocate a node array  
         do i=1, node_n
            read(1,*) bufi, Node(i).x(1:2), Node(i).bc(3:5), Node(i).pm(3:5)
            Node(i).bc(1:2) = 1 ! fixed DOFs not used in plate problems.
            Node(i).bc(6) = 1   ! fixed DOFs not used in plate problems.
      
      
         end do

         ! read elemental information
         read(1,*) bufs; read(1,*) element_n; read(1,*) bufs
         allocate(Element(element_n)) ! allocate a element array
         do i=1, element_n
            read(1,*) bufi, Element(i).cn(:), Element(i).q(:)
         end do

         ! read properties
         read(1,*) bufs; read(1,*) thickness
         Element(:).thickness = thickness
         read(1,*) bufs; read(1,*) Material.Young
         read(1,*) bufs; read(1,*) Material.Poisson
         read(1,*) bufs; read(1,*) Material.Shear
      end subroutine Read_Input_File

! ---------------------------------------------------------------------------------------

! calculate number of total DOF, number of free DOF, number of fixed DOF
! and assign equation numbers to DOFs
      subroutine Set_DOF_Number(tn, fn, cn)
         ! number of total DOF, number of free DOF, number of fixed DOF
         integer, intent(out) :: tn, fn, cn  
         integer :: i, j
   
         write(3, *) "EQUATION NUMBER"
         write(3, *) "---------------"
         write(3, *) "    node   dof    eqn"
   
         tn = node_n * 6  ! number of total DOF
         fn = 0; cn = 0
         do i=1, node_n
            do j=1, 6
               if (Node(i).bc(j) == 0) then 
                  fn = fn + 1
                  Node(i).eq_n(j) = fn
                  write(3, "(3i7)") i, j, fn
               else
                  cn = cn + 1
                  Node(i).eq_n(j) = tn - cn + 1
               end if
            end do
         end do
         write(3, *)      
      end subroutine Set_DOF_Number

! ---------------------------------------------------------------------------------------

! calculate column index for skyline solver
      subroutine Calculate_Index()   
         integer :: column_h(total_dof_n) ! column height
         integer :: a_index(6*4)          ! index for assemblage
         integer :: en                    ! number of element nodes
         integer :: i, j, k               ! index for iteration
         integer :: buf_sum

         ! allocate c_index array
         allocate ( c_index(free_dof_n+1) ) 
   
         ! column height
         column_h(:) = 0
         do i=1, element_n
            ! assemblage index
            en = 4
            do j=1, 6
               do k=1, en
                 a_index( en*j+k-en ) = node( element(i).cn(k) ).eq_n(j)
               end do
            end do      
            ! column height   
            do k=1, 6*en
               do j=1, 6*en
                  if(a_index(j)<=free_dof_n .and. a_index(k)<=free_dof_n) then
                     if( a_index(j) < a_index(k) ) then
                        if( a_index(k)-a_index(j) > column_h(a_index(k)) ) then
                           column_h(a_index(k)) = a_index(k) - a_index(j)
                        end if
                     end if
                  end if
               end do
            end do  
         end do

         ! c_index array
         buf_sum = 1
         do i=1, free_dof_n
            c_index(i) = buf_sum
            buf_sum = buf_sum + Column_H(i) + 1
         end do
         c_index(free_dof_n+1) = buf_sum
         write(3,"(a18, i4, a3)") "REQUIRED MEMORY =", buf_sum*8/1000000, " MB"
         write(3,*)
      end subroutine Calculate_Index

! ---------------------------------------------------------------------------------------

! assemble total stiffness matrix by using equation numbers
      subroutine Assemble_Kt()
         double precision :: eNode(4,2) ! nodal position of 4-node element (x,y)
         double precision :: Ke(8,8)    ! stifness matrix of element
         integer :: a_index(8)          ! assemblage index
         integer :: i, j, k, address

         Kt(:) = 0.0d0
  
         do i=1, element_n
            ! nodal position of element
            do j=1, 4
              do k=1, 2
                 eNode(j, k) = Node( Element(i).cn(j) ).x(k)
              end do
            end do  
            ! calculate stiffness matrix of element
            call Plane_Stiffness(Material.Young, Material.Poisson, Material.Shear, Element(i).thickness, eNode, Ke)
            write(2,"(a24,i4)") " STIFFNESS of ELEMENT : ", i
            call Print_Matrix(Ke)      
            ! assemblage index
            do j=1, 2
               do k=1, 4
                  a_index( 4*j+k-4 ) = node( element(i).cn(k) ).eq_n(j)
               end do
            end do
            ! assemble total stiffness matrix
            do j=1, 8
               do k=1, 8
                  if(a_index(j)<=free_dof_n .and. a_index(k)<=free_dof_n) then
                     if( a_index(j) <= a_index(k) ) then
                        address = c_index(a_index(k)) + a_index(k) - a_index(j)
                        Kt(address) = Kt(address) + Ke(j, k)
                     end if
                  end if
               end do
            end do
         end do  
      end subroutine Assemble_Kt

! ---------------------------------------------------------------------------------------

! assemble load vector			
      subroutine Assemble_Load()						
         real(8) :: eNode(4,2)
         real(8) :: NodalLoad(8) ! equivalent nodal load
         integer :: i, j, k

         R(:) = 0

         ! assemble load vector for nodal load
         do i=1, node_n
            do j=1, 2
               R(Node(i).eq_n(j)) = Node(i).pm(j)
            end do
         end do
  
         ! assemble load vector for body force
         do i=1, element_n
            ! nodal position of element
            do j=1, 4
              do k=1, 2
                 eNode(j, k) = Node( Element(i).cn(j) ).x(k)
              end do
            end do
            ! calculate equivalent nodal load from body force
            call Plane_Load(eNode, Element(i).q, NodalLoad)
            ! assemble load vector
            do j=1, 2
               do k=1, 4
                  R( Node( Element(i).cn(k) ).eq_n(j) ) = R( Node( Element(i).cn(k) ).eq_n(j) ) &
                                                  + NodalLoad(4*j+k-4)
               end do   	 
            end do
         end do
      end subroutine Assemble_Load

! ---------------------------------------------------------------------------------------

! solve linear equations
      subroutine Solve_Equations()
         integer :: i

         U(:) = 0.0d0; i=3
         U(1:free_dof_n) = R(1:free_dof_n)
         call COLSOL(Kt(1:c_index(free_dof_n+1)-1), U(1:free_dof_n), c_index(1:free_dof_n+1), &
                     free_dof_n,  c_index(free_dof_n+1)-1, free_dof_n+1, 1, i)
         call COLSOL(Kt(1:c_index(free_dof_n+1)-1), U(1:free_dof_n), c_index(1:free_dof_n+1), &
                     free_dof_n,  c_index(free_dof_n+1)-1, free_dof_n+1, 2, i) 
      end subroutine Solve_Equations

! ---------------------------------------------------------------------------------------

! calculate stress and print solutions
      subroutine Displacement_Stress()
         double precision :: eNode(4,2)  ! nodal position of element
         double precision :: displace(8) ! nodal displacement vector of element
   
         double precision :: Stress(4,3) ! Sxx, Syy, Sxy in Gauss points or nodal positions (2*2)
         integer :: i, j, k
   
         ! print strain energy
         write(3,"(a17, E14.6)") "STRAIN ENERGY = ", 0.5d0*dot_product(R(1:free_dof_n), U(1:free_dof_n))
         write(4,*) element_n

         ! print nodal displacement
         write(3,*)
         write(3,*) "DISPLACEMENT "
         write(3,*) "------------"
         write(3,*) "  Node      Dx         Dy     "
         do i=1, node_n
            write(3,"(1x,i4,2x,2(1P,E11.3))") i, U(Node(i).eq_n(1)), U(Node(i).eq_n(2))
         end do
         write(3,*)

         do i=1, element_n
            ! nodal position of element
            do j=1, 4
              do k=1, 2
                 eNode(j, k) = Node( Element(i).cn(j) ).x(k)
              end do
            end do
            ! displacement vector of element 
            do j=1, 2
               do k=1, 4
                  displace(4*j+k-4) = U( Node( Element(i).cn(k) ).eq_n(j) )
               end do
            end do
            ! calculate stress of element
            call Plane_Stress(Material.Young, Material.Poisson, eNode, displace, Stress)
            ! print stress
            write(3,"(a21,i4)") " STRESS of ELEMENT : ", i
            write(3,*) "------------------------" 
            write(3,*) " Position       Sxx        Syy        Sxy     "
            do j=1, 4
               write(3,"(1x,i4,5x, 3x,3(1P,E11.3))") j, Stress(j,:)
            end do
            write(3,*)
           ! print deformed shape and stress for post processing by using MATLAB 
            write(4,"(1x,28(1P,E13.5))") eNode(1,:), displace(1), displace(5), Stress(1,:),&
                                         eNode(2,:), displace(2), displace(6), Stress(2,:),&
                                         eNode(3,:), displace(3), displace(7), Stress(3,:),&
                                         eNode(4,:), displace(4), displace(8), Stress(4,:)
         end do
      end subroutine Displacement_Stress

! ---------------------------------------------------------------------------------------

! print 8x8 matrix
      subroutine Print_Matrix(M)
         double precision, intent(in) :: M(:,:)
         integer :: i

         write(2,*) "---------------------------"
         do i=1, 8
            write(2,"(8E12.4)" ) M(i,:)
         end do
         write(2,*)
      end subroutine Print_Matrix

! ---------------------------------------------------------------------------------------

   subroutine COLSOL(A,V,MAXA,NN,NWK,NNM,KKK,IOUT)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO SOLVE FINITE ELEMENT STATIC EQUILIBRIUM EQUATIONS IN    . 
! .        CORE, USING COMPACTED STORAGE AND COLUMN REDUCTION SCHEME  . 
! .                                                                   . 
! .  - - INPUT VARIABLES - -                                          . 
! .        A(NWK)    = STIFFNESS MATRIX STORED IN COMPACTED FORM      . 
! .        V(NN)     = RIGHT-HAND-SIDE LOAD VECTOR                    . 
! .        MAXA(NNM) = VECTOR CONTAINING ADDRESSES OF DIAGONAL        . 
! .                    ELEMENTS OF STIFFNESS MATRIX IN A              . 
! .        NN        = NUMBER OF EQUATIONS                            . 
! .        NWK       = NUMBER OF ELEMENTS BELOW SKYLINE OF MATRIX     . 
! .        NNM       = NN + 1                                         . 
! .        KKK       = INPUT FLAG                                     . 
! .            EQ. 1   TRIANGULARIZATION OF STIFFNESS MATRIX          . 
! .            EQ. 2   REDUCTION AND BACK-SUBSTITUTION OF LOAD VECTOR . 
! .        IOUT      = UNIT NUMBER USED FOR OUTPUT                    . 
! .                                                                   . 
! .  - - OUTPUT - -                                                   . 
! .        A(NWK)    = D AND L - FACTORS OF STIFFNESS MATRIX          . 
! .        V(NN)     = DISPLACEMENT VECTOR                            . 
! .                                                                   . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     . 
! .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      . 
! .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     . 
! .   SINGLE PRECISION ARITHMETIC.                                    . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      double precision, intent(inout) :: A(NWK), V(NN)
      integer,          intent(in)    :: MAXA(NNM)                                  
!                                                                       
!     PERFORM L*D*L(T) FACTORIZATION OF STIFFNESS MATRIX                
!                                                                       
      IF (KKK-2) 40,150,150                                             
   40 DO 140 N=1,NN                                                     
      KN=MAXA(N)                                                        
      KL=KN + 1                                                         
      KU=MAXA(N+1) - 1                                                  
      KH=KU - KL                                                        
      IF (KH) 110,90,50                                                 
   50 K=N - KH                                                          
      IC=0                                                              
      KLT=KU                                                            
      DO 80 J=1,KH                                                      
      IC=IC + 1                                                         
      KLT=KLT - 1                                                       
      KI=MAXA(K)                                                        
      ND=MAXA(K+1) - KI - 1                                             
      IF (ND) 80,80,60                                                  
   60 KK=MIN0(IC,ND)                                                    
      C=0.                                                              
      DO 70 L=1,KK                                                      
   70 C=C + A(KI+L)*A(KLT+L)                                            
      A(KLT)=A(KLT) - C                                                 
   80 K=K + 1                                                           
   90 K=N                                                               
      B=0.                                                              
      DO 100 KK=KL,KU                                                   
      K=K - 1                                                           
      KI=MAXA(K)                                                        
      C=A(KK)/A(KI)                                                     
      B=B + C*A(KK)                                                     
  100 A(KK)=C                                                           
      A(KN)=A(KN) - B                                                   
  110 IF (A(KN)) 120,120,140                                            
  120 WRITE (iout,2000) N,A(KN)
      WRITE (*,2000) N,A(KN)
      pause                              
      GO TO 800                                                         
  140 CONTINUE                                                          
      GO TO 900                                                         
!                                                                       
!     REDUCE RIGHT-HAND-SIDE LOAD VECTOR                                
!                                                                       
  150 DO 180 N=1,NN                                                     
      KL=MAXA(N) + 1                                                    
      KU=MAXA(N+1) - 1                                                  
      IF (KU-KL) 180,160,160                                            
  160 K=N                                                               
      C=0.                                                              
      DO 170 KK=KL,KU                                                   
      K=K - 1                                                           
  170 C=C + A(KK)*V(K)                                                  
      V(N)=V(N) - C                                                     
  180 CONTINUE                                                          
!                                                                       
!     BACK-SUBSTITUTE                                                   
!                                                                       
      DO 200 N=1,NN                                                     
      K=MAXA(N)                                                         
  200 V(N)=V(N)/A(K)                                                    
      IF (NN.EQ.1) GO TO 900                                            
      N=NN                                                              
      DO 230 L=2,NN                                                     
      KL=MAXA(N) + 1                                                    
      KU=MAXA(N+1) - 1                                                  
      IF (KU-KL) 230,210,210                                            
  210 K=N                                                               
      DO 220 KK=KL,KU                                                   
      K=K - 1                                                           
  220 V(K)=V(K) - A(KK)*V(N)                                            
  230 N=N - 1                                                           
      GO TO 900                                                         
!                                                                       
  800 STOP                                                              
  900 RETURN                                                            
!                                                                       
 2000 FORMAT (//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//,  &  
                ' NONPOSITIVE PIVOT FOR EQUATION ',I8,//,             &  
                ' PIVOT = ',E20.12 )                                    
 
   end subroutine COLSOL
! ---------------------------------------------------------------------------------------   
END MODULE Analysis
