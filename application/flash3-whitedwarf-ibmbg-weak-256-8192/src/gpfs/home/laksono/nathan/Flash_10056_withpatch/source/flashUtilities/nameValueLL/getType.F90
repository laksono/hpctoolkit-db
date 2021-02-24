!!****if* source/flashUtilities/nameValueLL/getType
!!
!! NAME
!!   getType
!!
!! SYNOPSIS
!!   getType(context_type(INOUT) :: context,
!!           character(IN)       :: name,
!!           integer(OUT)        :: name_type)
!!
!! DESCRIPTION: 
!!    Given the name in a linked list and a context to search, find
!!    the type (real, integer, etc.) of the parameter.  If the
!!    parameter is not found, return a value which indicates that
!!    the parameter is invalid.
!!
!! ARGUMENTS
!!    context:         name of context
!!    name:          name (in)
!!    name_type:     type of parameter (out)
!!
!!***

subroutine getType (context, name, name_type)
  
  use nameValueLL_data !, ONLY: context_type, nameValueLL_find, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type

  implicit none

  type (context_type), intent(inout)       :: context
  character(len=*), intent(in)          :: name
  integer, intent(out)                  :: name_type
  type (real_list_type), pointer:: real_test
  type (int_list_type), pointer :: int_test
  type (str_list_type), pointer :: str_test
  type (log_list_type), pointer :: log_test

  
  call nameValueLL_find (context, name, real_test)
  call nameValueLL_find (context, name, int_test)
  call nameValueLL_find (context, name, str_test)
  call nameValueLL_find (context, name, log_test)
  
  name_type = name_invalid
  !               Only one of the following should test true if the filters
  !               in add_*_
  if (associated(real_test)) name_type = name_real
  if (associated(int_test))  name_type = name_int
  if (associated(str_test))  name_type = name_str
  if (associated(log_test))  name_type = name_log

  return
end subroutine getType
