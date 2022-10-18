program test
   use mod_shell_correction
   implicit none
   integer :: z, n, a
   real, dimension(4) :: bt
   real :: xshell
   bt = 0.
   z = 82; n = 126
   xshell = get_e_shell(z,n,bt(2),bt(3),bt(4))
   a = z+n
   write(*,*)xshell, ' <- old value'
   write(*,*)z,n,a
end program test
