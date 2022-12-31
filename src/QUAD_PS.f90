module QUAD_PS  ! -------------- Begin of Module -----------------------------------------
   implicit none
   
   public Plane_Stiffness, Plane_Stress, Plane_Load
   
   private Gauss22_Point, Shape_Function ,dHrs_Matrix
   private Strain_displacement, dHxy_Matrix, Material_law

contains  ! ------------------------------ Internal Procedures --------------------------

! ---------------------------------------------------------------------------------------

! calculate stiffness matrix of plane stress element
! DOF vector (u1,u2,u3,u4, v1,v2,v3,v4)
      subroutine Plane_Stiffness(Young, Poisson, thick, node, Ke)
         real(8), intent(in)  :: Young, Poisson  ! Young's modulus, Poison's ratio
         real(8), intent(in)  :: thick           ! thickness
         real(8), intent(in)  :: node(4,2)       ! nodal position of element
         real(8), intent(out) :: Ke(8,8)         ! stiffness matrix
  
         real(8) :: r, s, weight                 ! integration point, weight factor
         real(8) :: det_j                        ! determinant of Jacobian
         real(8) :: B(3,8)                       ! B-matrix (strain-displacement matrix)
         real(8) :: C(3,3)                       ! material matrix (material law)
   
         integer :: i, j

         ! calculate a material matrix
         call Material_Law(Young, Poisson, C)

         ! numerical integration
         Ke(:,:) = 0.0d0  
         do i=1, 2
            do j=1, 2
               ! Gauss interation point and weight factor
               call Gauss22_Point(i, j, r, s, weight)
               ! calculate determinant of Jacobian and B-matrix
               call Strain_displacement(r, s, node, det_j, B)
               ! integration
               Ke = Ke + weight * thick * matmul(matmul(transpose(B),C),B) * det_j  
            end do  
         end do
      end subroutine Plane_Stiffness

! ---------------------------------------------------------------------------------------

! calculate stress (Sx, Sy, Sxy)
      subroutine Plane_Stress(Young, Poisson, Node, Displace, Stress)
         real(8), intent(in)  :: Young, Poisson  ! Young's modulus, Possion's ratio
         real(8), intent(in)  :: node(4,2)       ! nodal position of element
         real(8), intent(in)  :: displace(8)     ! displacement vector
         real(8), intent(out) :: Stress(4,3)     ! stress (Sx, Sy, Sxy)
         real(8) :: C(3,3), B(3,8)
         real(8) :: buf_det
         real(8) :: r, s, weight
         integer :: i, j, r4(4), s4(4)
         data r4 / 1.0d0, -1.0d0, -1.0d0,  1.0d0 /
         data s4 / 1.0d0,  1.0d0, -1.0d0, -1.0d0 /
 
         ! calculate a material matrix
         call Material_law(Young, Poisson, C)
   
         do i=1, 2
            do j=1, 2
               ! set postions where stresses are out
               !call Gauss22_Point(i, j, r, s, weight) ! at Gauss points 
               r = r4((i-1)*2+j); s=s4((i-1)*2+j)      ! at nodal positions of element
               ! calculate B-matrix
               call strain_displacement(r, s, node, buf_det, B)
               ! calculate stresses
               Stress(i*2+j-2,:) = matmul( matmul(C,B), displace )
            end do
         end do
      end subroutine Plane_Stress

! ---------------------------------------------------------------------------------------

! equivalent nodal loads
      subroutine Plane_Load(node, q, nodal_load)
         real(8), intent(in)  :: node(4,2)      ! node position
         real(8), intent(in)  :: q(2)           ! body forces in x- and y-directions
         real(8), intent(out) :: nodal_load(8)  ! element load vector
 
         real(8) :: H(4)                        ! shape functions
         real(8) :: r, s, weight                ! Gauss point, weight factor
         real(8) :: det_j                       ! determinant of Jacobian
         integer          :: i, j

         ! numerical integration
         nodal_load(:) = 0.0d0
         do i=1, 2
            do j=1, 2
               ! Gauss points
               call Gauss22_Point(i, j, r, s, weight)      
               ! determinant of Jacobian and shape function (H)					   
               call H_Matrix(r, s, node, H, det_j)      
               ! equivalent nodal load vector
               nodal_load(1:4) = nodal_load(1:4) + weight*dabs(det_j)*H(:)*q(1)
               nodal_load(5:8) = nodal_load(5:8) + weight*dabs(det_j)*H(:)*q(2)
            end do
         end do
      end subroutine Plane_Load

! ---------------------------------------------------------------------------------------

! Gauss integration point (2*2) and weight factor
      subroutine Gauss22_Point(i, j, r, s, weight)
         integer, intent(in)           :: i, j
         real(8), intent(out) :: r, s, weight

         real(8) :: GaussPoint(2), w(2)
         data GaussPoint / -0.577350269d0, 0.577350269d0 / 
         data w          /  1.000000000d0, 1.000000000d0 / 
   
         weight = w(i)*w(j) 
         r = GaussPoint(i)
         s = GaussPoint(j)   
      end subroutine Gauss22_Point

! ---------------------------------------------------------------------------------------

! calculate shape functions
      subroutine Shape_Function(r, s, H)
         real(8), intent(in)  :: r, s  ! natural coordinate
         real(8), intent(out) :: H(4)  ! shape functions

         H(1) = 0.25d0 * (1.0d0 - r) * (1.0d0 - s)
         H(2) = 0.25d0 * (1.0d0 + r) * (1.0d0 - s)
         H(3) = 0.25d0 * (1.0d0 + r) * (1.0d0 + s)
         H(4) = 0.25d0 * (1.0d0 - r) * (1.0d0 + s)
      end subroutine Shape_Function

! ---------------------------------------------------------------------------------------

! shape function (H matrix), determinant of Jacobian
      subroutine H_Matrix(r, s, node, H, det_j)
         real(8), intent(in)  :: r, s
         real(8), intent(in)  :: node(4,2)
         real(8), intent(out) :: H(4), det_j

         real(8) :: dHrs(2,4), Jacob(2,2)

         ! shape function (H)
         call Shape_Function(r, s, H)
         ! dHrs = derivative of shape functions (dH/dr, dH/ds)
         call dHrs_Matrix(r, s, dHrs)
         ! Jacob = Jacobian matrix
         Jacob = matmul(dHrs, node)         
         ! determinant of jacobian matrix
         det_j = Jacob(1,1)*Jacob(2,2) - Jacob(1,2)*Jacob(2,1)
      end subroutine H_Matrix

! ---------------------------------------------------------------------------------------

! derivatives of shape functions. dHrs(1,:)=dH/dr, dHrs(2,:)=dH/ds
      subroutine dHrs_Matrix(r, s, dHrs)
         real(8), intent(in)  :: r, s
         real(8), intent(out) :: dHrs(2,4)

         dHrs(1,1) = -0.25d0 * (1.0d0 - s)
         dHrs(1,2) =  0.25d0 * (1.0d0 - s)
         dHrs(1,3) =  0.25d0 * (1.0d0 + s)
         dHrs(1,4) = -0.25d0 * (1.0d0 + s)

         dHrs(2,1) = -0.25d0 * (1.0d0 - r)
         dHrs(2,2) = -0.25d0 * (1.0d0 + r)
         dHrs(2,3) =  0.25d0 * (1.0d0 + r)
         dHrs(2,4) =  0.25d0 * (1.0d0 - r)
      end subroutine dHrs_Matrix

! ---------------------------------------------------------------------------------------

! dHxy matirx. dHxy(1,:)=dH/dx, dHxy(2,:)=dH/dy
      subroutine dHxy_Matrix(r, s, node, det_j, dHxy)
         real(8), intent(in)  :: r, s, Node(4,2)
         real(8), intent(out) :: det_j, dHxy(2,4) 

         real(8) :: buf, dHrs(2,4), Jacob(2,2)

         ! dHrs = derivative of shape function (dH/dr, dH/ds)
         call dHrs_Matrix(r, s, dHrs)

         ! Jacob = Jacobian Matrix
         Jacob = matmul(dHrs, node) 
        
         ! Jacob => inverse of jacobian Matrix
         det_j = Jacob(1,1)*Jacob(2,2) - Jacob(1,2)*Jacob(2,1)
         Jacob(1,2) = -Jacob(1,2)
         Jacob(2,1) = -Jacob(2,1)
         buf = Jacob(1,1)
         Jacob(1,1) = Jacob(2,2)
         Jacob(2,2) = buf
         Jacob = Jacob / det_j
 
         ! dHxy(1,:)=dH/dx, dHxy(2,:)=dH/dy
         dHxy = matmul(Jacob, dHrs)
      end subroutine dHxy_Matrix

! ---------------------------------------------------------------------------------------

! B-matrix (strain-displacement matrix)
      subroutine Strain_displacement(r, s, node, det_j, B)
         real(8), intent(in)  :: r, s, node(4,2)
         real(8), intent(out) :: det_j, B(3,8)

         real(8) :: dHxy(2,4)

         ! calculate dHxy. dHxy(1,:)=dH/dx, dH(2,:)=dHxy/dy
         call dHxy_Matrix(r, s, node, det_j, dHxy)

         ! B-matrix
         B(:,:)   = 0.0d0
         B(1,1:4) = dHxy(1,:)
         B(2,5:8) = dHxy(2,:)
         B(3,1:4) = dHxy(2,:)
         B(3,5:8) = dHxy(1,:)
      end subroutine Strain_displacement

! ---------------------------------------------------------------------------------------

! material law for plane stress condition
      subroutine Material_law(Young, Poisson, C)
         real(8), intent(in)  :: Young, Poisson  ! material constants
         real(8), intent(out) :: C(3,3)          ! material matrix

         C(:,:) = 0.0d0
         C(1,1) = Young / (1.0d0 - Poisson**2)
         C(2,2) = C(1,1)
         C(1,2) = Poisson * C(1,1)
         C(2,1) = C(1,2)
         C(3,3) = 0.5d0* Young / (1.0d0 + Poisson)
      end subroutine Material_law

! ---------------------------------------------------------------------------------------

end Module QUAD_PS  ! --------------- End of Module --------------------------------------
