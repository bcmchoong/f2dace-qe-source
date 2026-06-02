!# 1 "tg_gather.f90"
!----------------------------------------------------------------------------------------------------------------
!-real version 
SUBROUTINE tg_gather( dffts, v, tg_v )
  !
  USE fft_param
  USE fft_types,      ONLY : fft_type_descriptor
!# 8 "tg_gather.f90"
  IMPLICIT NONE
!# 10 "tg_gather.f90"
  TYPE(fft_type_descriptor), INTENT(in) :: dffts
  REAL(DP), INTENT(IN)  :: v(dffts%nnr)
  REAL(DP), INTENT(OUT) :: tg_v(dffts%nnr_tg)
!# 14 "tg_gather.f90"
  INTEGER :: nxyp, ir3, off, tg_off
  INTEGER :: i, nsiz, ierr
!# 17 "tg_gather.f90"
  nxyp   = dffts%nr1x*dffts%my_nr2p
  !
  !  The potential in v is distributed so that each Z-plane is shared among nproc2 processors.
  !  We collect the data of whole planes in tg_v to be used with task group distributed wfcs.
  !
  tg_v(:) = (0.d0,0.d0)
  do ir3 =1, dffts%my_nr3p
     off    = dffts%nr1x*dffts%my_nr2p*(ir3-1)
     tg_off = dffts%nr1x*dffts%nr2x   *(ir3-1) + dffts%nr1x*dffts%my_i0r2p
     tg_v(tg_off+1:tg_off+nxyp) = v(off+1:off+nxyp)
  end do
  !write (6,*) ' tg_v ', dffts%my_i0r2p, dffts%my_nr2p
  !write (6,'(20f12.7)') (v(dffts%my_i0r2p+i+dffts%nr1x*(i-1)), i=1,dffts%my_nr2p)
  !write (6,'(20f12.7)') (tg_v(i+dffts%nr1x*(i-1)), i=1,dffts%nr2x)
!# 41 "tg_gather.f90"
  !write (6,'(20f12.7)') (tg_v(i+dffts%nr1x*(i-1)), i=1,dffts%nr1x)
  RETURN
END SUBROUTINE tg_gather
!# 45 "tg_gather.f90"
!-complex version of previous routine
SUBROUTINE tg_cgather( dffts, v, tg_v )
  !
  USE fft_param
  USE fft_types,      ONLY : fft_type_descriptor
!# 51 "tg_gather.f90"
  IMPLICIT NONE
!# 53 "tg_gather.f90"
  TYPE(fft_type_descriptor), INTENT(in) :: dffts
  COMPLEX(DP), INTENT(IN)  :: v(dffts%nnr)
  COMPLEX(DP), INTENT(OUT) :: tg_v(dffts%nnr_tg)
!# 57 "tg_gather.f90"
  INTEGER :: nxyp, ir3, off, tg_off
  INTEGER :: i, nsiz, ierr
!# 60 "tg_gather.f90"
  nxyp   = dffts%nr1x*dffts%my_nr2p
  !
  !  The potential in v is distributed so that each Z-plane is shared among nproc2 processors.
  !  We collect the data of whole planes in tg_v to be used with task group distributed wfcs.
  !
  tg_v(:) = (0.d0,0.d0)
  do ir3 =1, dffts%my_nr3p
     off    = dffts%nr1x*dffts%my_nr2p*(ir3-1)
     tg_off = dffts%nr1x*dffts%nr2x   *(ir3-1) + dffts%nr1x*dffts%my_i0r2p
     tg_v(tg_off+1:tg_off+nxyp) = v(off+1:off+nxyp)
  end do
  !write (6,*) ' tg_v ', dffts%my_i0r2p, dffts%my_nr2p
  !write (6,'(20f12.7)') (v(dffts%my_i0r2p+i+dffts%nr1x*(i-1)), i=1,dffts%my_nr2p)
  !write (6,'(20f12.7)') (tg_v(i+dffts%nr1x*(i-1)), i=1,dffts%nr2x)
!# 84 "tg_gather.f90"
  !write (6,'(20f12.7)') (tg_v(i+dffts%nr1x*(i-1)), i=1,dffts%nr1x)
  RETURN
END SUBROUTINE tg_cgather
