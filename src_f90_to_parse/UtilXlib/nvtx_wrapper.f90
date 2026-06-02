!# 1 "nvtx_wrapper.f90"
!MIT License
!# 3 "nvtx_wrapper.f90"
!Copyright (c) 2019 maxcuda
!# 5 "nvtx_wrapper.f90"
!This module has been downloaded and adapted from
!   https://github.com/maxcuda/NVTX_example     
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!# 15 "nvtx_wrapper.f90"
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!# 18 "nvtx_wrapper.f90"
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!# 27 "nvtx_wrapper.f90"
! ----
! nvtx
! ----
!# 31 "nvtx_wrapper.f90"
module nvtx
  use iso_c_binding
!# 36 "nvtx_wrapper.f90"
  implicit none
!# 76 "nvtx_wrapper.f90"
contains
!# 78 "nvtx_wrapper.f90"
  subroutine nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
!# 98 "nvtx_wrapper.f90"
  end subroutine nvtxStartRange
!# 100 "nvtx_wrapper.f90"
  subroutine nvtxStartRangeAsync(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
!# 116 "nvtx_wrapper.f90"
  end subroutine nvtxStartRangeAsync
!# 119 "nvtx_wrapper.f90"
  subroutine nvtxEndRange
!# 127 "nvtx_wrapper.f90"
  end subroutine nvtxEndRange
!# 129 "nvtx_wrapper.f90"
  subroutine nvtxEndRangeAsync
!# 133 "nvtx_wrapper.f90"
  end subroutine nvtxEndRangeAsync
!# 135 "nvtx_wrapper.f90"
end module nvtx
