!# 1 "fox_init_module.f90"
MODULE fox_init_module
!# 8 "fox_init_module.f90"
IMPLICIT NONE 
PRIVATE
PUBLIC     :: fox_init 
!# 12 "fox_init_module.f90"
CONTAINS 
   SUBROUTINE fox_init() 
      INTEGER   :: errcodes(3)
!# 26 "fox_init_module.f90"
         errcodes(:) = 0
!# 28 "fox_init_module.f90"
   END SUBROUTINE fox_init
END MODULE fox_init_module
