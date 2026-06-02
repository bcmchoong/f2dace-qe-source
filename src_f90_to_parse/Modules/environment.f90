!# 1 "environment.f90"
!
! Copyright (C) 2002-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!# 11 "environment.f90"
!
!==-----------------------------------------------------------------------==!
MODULE environment
  !==-----------------------------------------------------------------------==!
  !! Environment management.
  USE kinds, ONLY: DP
  USE io_files, ONLY: crash_file, nd_nmbr
  USE io_global, ONLY: stdout, meta_ionode
  USE mp_world,  ONLY: nproc, nnode
  USE mp_images, ONLY: me_image, my_image_id, root_image, nimage, &
      nproc_image
  USE mp_pools,  ONLY: npool
  USE mp_bands,  ONLY: ntask_groups, nproc_bgrp, nbgrp, nyfft
  USE global_version, ONLY: version_number
  USE fox_init_module, ONLY: fox_init
  USE command_line_options, ONLY : nmany_
  USE clib_wrappers, ONLY : get_mem_avail
!# 32 "environment.f90"
  IMPLICIT NONE
!# 34 "environment.f90"
  ! ... code name
  CHARACTER(LEN=20) :: code = 'notset'
!# 37 "environment.f90"
  SAVE
!# 39 "environment.f90"
  PRIVATE
!# 41 "environment.f90"
  PUBLIC :: environment_start
  PUBLIC :: environment_end
  PUBLIC :: opening_message
  PUBLIC :: parallel_info
  PUBLIC :: print_cuda_info
!# 47 "environment.f90"
  !==-----------------------------------------------------------------------==!
CONTAINS
  !==-----------------------------------------------------------------------==!
!# 51 "environment.f90"
  SUBROUTINE environment_start( code_in )
!# 53 "environment.f90"
    CHARACTER(LEN=*), INTENT(IN) :: code_in
!# 55 "environment.f90"
    LOGICAL           :: exst, debug = .false.
    CHARACTER(LEN=80) :: code_version, uname
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    CHARACTER(LEN=3)           :: env_maxdepth
    INTEGER :: ios, crashunit, max_depth 
!# 68 "environment.f90"
    ! ... use ".FALSE." to disable all clocks except the total cpu time clock
    ! ... use ".TRUE."  to enable clocks
    CALL init_clocks(.TRUE.) 
    code = code_in
    CALL start_clock( TRIM(code) )
!# 74 "environment.f90"
    code_version = TRIM (code) // " v." // TRIM (version_number)
!# 76 "environment.f90"
    ! ... for compatibility with PWSCF
!# 80 "environment.f90"
    nd_nmbr = ' '
!# 83 "environment.f90"
    IF( meta_ionode ) THEN
!# 85 "environment.f90"
       ! ...  search for file CRASH and delete it
!# 87 "environment.f90"
       INQUIRE( FILE=TRIM(crash_file), EXIST=exst )
       IF( exst ) THEN
          OPEN( NEWUNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD',IOSTAT=ios )
          IF (ios==0) THEN
             CLOSE( UNIT=crashunit, STATUS='DELETE', IOSTAT=ios )
          ELSE
             WRITE(stdout,'(5x,"Remark: CRASH file could not be deleted")')
          END IF
       END IF
!# 97 "environment.f90"
    ELSE
       ! ... one processor per image (other than meta_ionode)
       ! ... or, for debugging purposes, all processors,
       ! ... open their own standard output file
!#define DEBUG
!# 105 "environment.f90"
       IF (me_image == root_image .OR. debug ) THEN
          uname = 'out.' // trim(int_to_char( my_image_id )) // '_' // &
               trim(int_to_char( me_image))
          OPEN ( unit = stdout, file = TRIM(uname),status='unknown')
       ELSE
!# 113 "environment.f90"
          OPEN ( unit = stdout, file='/dev/null', status='unknown' )
!# 115 "environment.f90"
       END IF
!# 117 "environment.f90"
    END IF
    !
    CALL opening_message( code_version )
!# 123 "environment.f90"
    CALL serial_info()
!# 125 "environment.f90"
    CALL fox_init()
!# 129 "environment.f90"
    !
    WRITE(stdout,'(5x, I0, A, A)') get_mem_avail()/1024, &
                &" MiB available memory on the printing compute node ", &
                &"when the environment starts"
    WRITE(stdout, *)
  END SUBROUTINE environment_start
!# 136 "environment.f90"
  !==-----------------------------------------------------------------------==!
!# 138 "environment.f90"
  SUBROUTINE environment_end( code_in )
    ! next line for back-compatibility
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: code_in
!# 144 "environment.f90"
    IF ( meta_ionode ) WRITE( stdout, * )
!# 146 "environment.f90"
    IF ( PRESENT(code_in) ) THEN
       code = code_in
    END IF
    IF ( code == 'notset' ) THEN
       WRITE( stdout,'(5X,"WARNING: environment_end needs a call to environment_start at the beginning of the run")')
    ELSE
       CALL stop_clock(  TRIM(code) )
       CALL print_clock( TRIM(code) )
    END IF
    CALL closing_message( )
!# 157 "environment.f90"
    IF( meta_ionode ) THEN
       WRITE( stdout,'(A)')      '   JOB DONE.'
       WRITE( stdout,3335)
    END IF
3335 FORMAT('=',78('-'),'=')
    FLUSH(stdout)
!# 164 "environment.f90"
    RETURN
  END SUBROUTINE environment_end
!# 167 "environment.f90"
  !==-----------------------------------------------------------------------==!
!# 169 "environment.f90"
  SUBROUTINE opening_message( code_version )
!# 171 "environment.f90"
    CHARACTER(LEN=*), INTENT(IN) :: code_version
    CHARACTER(LEN=9)  :: cdate, ctime
!# 174 "environment.f90"
    CALL date_and_tim( cdate, ctime )
    !
    WRITE( stdout, '(/5X,"Program ",A," starts on ",A9," at ",A9)' ) &
         TRIM(code_version), cdate, ctime
!# 188 "environment.f90"
    !
    WRITE( stdout, '(/5X,"This program is part of the open-source Quantum ",&
         &    "ESPRESSO suite", &
         &/5X,"for quantum simulation of materials; please cite",   &
         &/9X,"""P. Giannozzi et al., J. Phys.:Condens. Matter 21 ",&
         &    "395502 (2009);", &
         &/9X,"""P. Giannozzi et al., J. Phys.:Condens. Matter 29 ",&
         &    "465901 (2017);", &
         &/9X,"""P. Giannozzi et al., J. Chem. Phys. 152 ",&
         &    "154105 (2020);", &
         &/9X," URL http://www.quantum-espresso.org"", ", &
         &/5X,"in publications or presentations arising from this work. More details at",&
         &/5x,"http://www.quantum-espresso.org/quote")' )
    RETURN
  END SUBROUTINE opening_message
!# 204 "environment.f90"
  !==-----------------------------------------------------------------------==!
!# 206 "environment.f90"
  SUBROUTINE closing_message( )
!# 208 "environment.f90"
    CHARACTER(LEN=9)  :: cdate, ctime
    CHARACTER(LEN=80) :: time_str
!# 211 "environment.f90"
    CALL date_and_tim( cdate, ctime )
!# 213 "environment.f90"
    time_str = 'This run was terminated on:  ' // ctime // ' ' // cdate
!# 215 "environment.f90"
    IF( meta_ionode ) THEN
       WRITE( stdout,*)
       WRITE( stdout,3334) time_str
       WRITE( stdout,3335)
    END IF
!# 221 "environment.f90"
3334 FORMAT(3X,A60,/)
3335 FORMAT('=',78('-'),'=')
!# 224 "environment.f90"
    RETURN
  END SUBROUTINE closing_message
!# 227 "environment.f90"
  !==-----------------------------------------------------------------------==!
  SUBROUTINE parallel_info ( )
    !
!# 241 "environment.f90"
    WRITE( stdout, '(/5X,"Parallel version (MPI), running on ",&
         &I5," processors")' ) nproc
!# 244 "environment.f90"
    !
    WRITE( stdout, '(/5X,"MPI processes distributed on ",&
         &I5," nodes")' ) nnode
    IF ( nimage > 1 ) WRITE( stdout, &
         '(5X,"path-images division:  nimage    = ",I7)' ) nimage
    IF ( npool > 1 ) WRITE( stdout, &
         '(5X,"K-points division:     npool     = ",I7)' ) npool
    IF ( nbgrp > 1 ) WRITE( stdout, &
         '(5X,"band groups division:  nbgrp     = ",I7)' ) nbgrp
    IF ( nproc_bgrp > 1 ) WRITE( stdout, &
         '(5X,"R & G space division:  proc/nbgrp/npool/nimage = ",I7)' ) nproc_bgrp
    IF ( nyfft > 1 ) WRITE( stdout, &
         '(5X,"wavefunctions fft division:  Y-proc x Z-proc = ",2I7)' ) &
         nyfft, nproc_bgrp / nyfft
    IF ( ntask_groups > 1 ) WRITE( stdout, &
         '(5X,"wavefunctions fft division:  task group distribution",/,34X,"#TG    x Z-proc = ",2I7)' ) &
         ntask_groups, nproc_bgrp / ntask_groups
    IF ( nmany_ > 1) WRITE( stdout, '(5X,"FFT bands division:     nmany     = ",I7)' ) nmany_
    !
  END SUBROUTINE parallel_info
!# 265 "environment.f90"
  !==-----------------------------------------------------------------------==!
  SUBROUTINE serial_info ( )
    !
!# 271 "environment.f90"
    !
!# 277 "environment.f90"
    WRITE( stdout, '(/5X,"Serial version")' )
!# 279 "environment.f90"
    !
  END SUBROUTINE serial_info
!
!-----------------------------------------------------------------------
SUBROUTINE print_cuda_info(check_use_gpu) 
  !-----------------------------------------------------------------------
  !
  USE io_global,       ONLY : stdout
  USE control_flags,   ONLY : use_gpu_=> use_gpu, iverbosity
  USE mp_world,        ONLY : nnode, nproc
  USE mp,              ONLY : mp_sum, mp_max
!# 293 "environment.f90"
  !
  IMPLICIT NONE
  !
  LOGICAL, OPTIONAL,INTENT(IN)  :: check_use_gpu 
  !! if present and trues the internal variable use_gpu is checked
!# 353 "environment.f90"
  !
END SUBROUTINE print_cuda_info
!# 356 "environment.f90"
  !==-----------------------------------------------------------------------==!
END MODULE environment
!==-----------------------------------------------------------------------==!
