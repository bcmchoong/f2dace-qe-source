!# 1 "casino2upf.f90"
!
! Copyright (C) 2008 Simon Binnie
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM casino2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in CASINO tabulated
  !     format to unified pseudopotential format
!# 16 "casino2upf.f90"
  USE casino_pp
  USE write_upf_new, ONLY : write_upf
  USE pseudo_types,  ONLY : deallocate_pseudo_upf, pseudo_upf
  ! 
  IMPLICIT NONE
  !
  CHARACTER(len=256) :: pp_data
  CHARACTER(len=256) :: upf_file
  CHARACTER(len=256), ALLOCATABLE:: wavefile(:)
  INTEGER, ALLOCATABLE :: waveunit(:)
  INTEGER nofiles, i, ios, pp_unit
  TYPE(pseudo_upf)      :: upf_out
!# 29 "casino2upf.f90"
  NAMELIST / inputpp / &
       pp_data,        &         !CASINO pp filename
       upf_file,        &         !output file
       tn_grid,        &         !.true. if Trail and Needs grid is used
       tn_prefac,      &
       xmin,           &         !xmin for standard QE grid
       dx                        !dx for Trail and Needs and standard QE
                                 !grid
   pp_data= 'pp.data'
   upf_file= 'out.UPF'
!# 40 "casino2upf.f90"
   WRITE(0,*) 'CASINO2UPF Converter'
!# 42 "casino2upf.f90"
   READ(*,inputpp,iostat=ios)
   READ(*,*,iostat=ios) nofiles
!# 45 "casino2upf.f90"
   ALLOCATE(wavefile(nofiles), waveunit(nofiles))
!# 47 "casino2upf.f90"
   !Now read in the awfn file names and open the files
!# 49 "casino2upf.f90"
   DO i=1,nofiles
      READ(*,*,iostat=ios) wavefile(:)
      OPEN(newunit=waveunit(i),file=trim(wavefile(i)),&
           status='old',form='formatted', iostat=ios)
      CALL upf_error ('casino2upf', 'cannot read file'//trim(wavefile(i)), ios)
   ENDDO
!# 56 "casino2upf.f90"
   OPEN(newunit=pp_unit,file=trim(pp_data),status='old',form='formatted', iostat=ios)
   CALL upf_error ('casino2upf', 'cannot read file'//trim(pp_data), ios)
!# 59 "casino2upf.f90"
   CALL read_casino(pp_unit,nofiles, waveunit)
!# 61 "casino2upf.f90"
   CLOSE (unit=pp_unit)
   DO i=1,nofiles
      CLOSE (waveunit(i))
   ENDDO
!# 66 "casino2upf.f90"
   DEALLOCATE( wavefile, waveunit )
!# 68 "casino2upf.f90"
   ! convert variables read from CASINO format into those needed
   ! by the upf format - add missing quantities
!# 71 "casino2upf.f90"
   CALL convert_casino(upf_out)
!# 73 "casino2upf.f90"
   PRINT '(''Output PP file in UPF format :  '',a)', upf_file
   CALL write_upf(filename = TRIM(upf_file), UPF = upf_out, SCHEMA = 'V2') 
   CALL  deallocate_pseudo_upf( upf_out )
!# 77 "casino2upf.f90"
   STOP
!# 79 "casino2upf.f90"
END PROGRAM casino2upf
