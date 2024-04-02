module c_interface
   use, intrinsic :: iso_c_binding, only: C_INT,C_CHAR,C_NULL_CHAR
   implicit none

   !> Interface to C mkdir
   interface 
      ! function mkdirf(path) bind(C,name="c_mkdir755")
      !    ! use iso_c_binding, only: C_INT, C_CHAR
      !    import
      !    implicit none
      !    integer(C_INT) :: mkdir
      !    character(C_CHAR) :: path
      ! end function mkdirf
      subroutine F_mkdir(path) bind(C,name="c_mkdir755noreturn")
         ! use iso_c_binding, only: C_INT, C_CHAR
         import
         implicit none
         character(C_CHAR) :: path
      end subroutine F_mkdir
   end interface

   contains

      subroutine mkdir(path)
         implicit none
         character(len=*) :: path
         character(len=len_trim(path)+1) :: c_path

         c_path=trim(path)//C_NULL_CHAR
         call F_mkdir(c_path)
      end subroutine mkdir


end module c_interface
