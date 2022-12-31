      program FEM_2D  ! ------------------------- Begin of Program ---------------------------
   
         use UDT
         use Analysis
         use QUAD_PS

         implicit none
  
         ! define internal variables
         integer :: node_n, element_n                             ! number of nodes, number of elements
         type (NodeType),    allocatable, dimension(:) :: Node    ! node
         type (ElementType), allocatable, dimension(:) :: Element ! element
         type (MaterialType) :: Material                          ! material

         integer :: total_dof_n, free_dof_n, fixed_dof_n          ! number of DOFs (total, free, fixed)
         real(8), allocatable, dimension(:) :: Kt        ! stiffness vector
         real(8), allocatable, dimension(:) :: U         ! diplacement vector
         real(8), allocatable, dimension(:) :: R         ! load vector
         integer,          allocatable, dimension(:) :: c_index   ! column index

         call Main ! main subroutine

      contains  
! ------------------------------ Internal Procedures --------------------------
! ---------------------------------------------------------------------------------------

! main subroutine
      subroutine Main()
         character(len=20) :: filename
   
         ! input file name
         write(*, "(' 1/ Input File Name (.txt): ',$)")
         read(*,*) filename
   
         ! open files
         open(unit=1, file=trim(filename)//".txt",     form="formatted")
         open(unit=2, file=trim(filename)//"_out.txt", form="formatted")
         open(unit=3, file=trim(filename)//"_res.txt", form="formatted")
         open(unit=4, file=trim(filename)//"_pos.txt", form="formatted")

         ! read input file    
         print *, "2/ Reading Input File";
         call Read_Input_File
   
         ! calculate # of total DOF, # of free DOF, # of fixed DOF
         ! and assign equation numbers
         call Set_DOF_Number(total_dof_n, free_dof_n, fixed_dof_n)
   
         ! calculate column index
         call Calculate_Index
   
         ! allocate memory
         allocate( Kt( c_index(free_dof_n+1)-1 ) ) ! total stiffness vector (Kt)
         allocate( U(total_dof_n) )                ! displacement vector
         allocate( R(total_dof_n) )                ! load vector
   
         ! assemble total stiffness matrix
         print *, "3/ Assembling Stiffness and Load"
         call Assemble_Kt

         ! assemble load vector
         call Assemble_Load

         ! slove linear system (find U in K*U=R)
         print *, "4/ Solving Linear System"
         call Solve_Equations 

         ! calculate stress and print solutions
         print *, "5/ Printing Output Files"
         call Displacement_Stress

         ! deallocate memory
         deallocate(Node, Element, Kt, R, U, c_index)

         ! close files
         close(unit=1); close(unit=2); close(unit=3); close(unit=4)
         print *, "6/ Completed !"; print *, ""; pause
      end subroutine Main
! ---------------------------------------------------------------------------------------


end program FEM_2D  ! ------------------------- End of Program -------------------------
