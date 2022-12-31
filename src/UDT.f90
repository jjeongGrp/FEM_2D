! ---------------------------------------------------------------------------------------
!            Module for defining user-defined types (UDT)
! ---------------------------------------------------------------------------------------
 
module UDT  ! -------------- Begin of Module ---------------------------------------

   ! define a type for node
   type :: NodeType
      real(8) :: x(3)       ! nodal position (x, y, z)
      real(8) :: pm(6)      ! nodal force (Px, Py, Pz, Mx, My, Mz)
      integer :: bc(6)      ! displacement BC (u, v, w, Rx, Ry, Rz) (1=fixed, 0=free)
      integer :: eq_n(6)    ! equation number (u, v, w, Rx, Ry, Rz)
   end type NodeType

   ! define a type for element 
   type :: ElementType
      integer :: cn(4)      ! connectivity
      real(8) :: thickness  ! thickness
      real(8) :: q(2)       ! distributed load in x- and y-directions
   end type ElementType

   ! define a type for material
   type :: MaterialType
      real(8) :: Young      ! Young's modulus
      real(8) :: Poisson    ! Poison's ratio
   end type MaterialType 

end module UDT  ! --------------- End of Module ----------------------------------
