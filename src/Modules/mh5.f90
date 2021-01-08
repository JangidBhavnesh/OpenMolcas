!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! Generig interfaces for HDF5

module mh5

implicit none
private
public :: mh5_create_file, mh5_open_file_rw, mh5_open_file_r, mh5_close_file, mh5_is_hdf5, mh5_open_group, mh5_close_group, &
          mh5_exists_dset, mh5_exists_attr, mh5_open_dset, mh5_close_dset, mh5_open_attr, mh5_close_attr, mh5_create_attr_int, &
          mh5_create_attr_real, mh5_create_attr_str, mh5_put_attr, mh5_get_attr, mh5_init_attr, mh5_fetch_attr, &
          mh5_create_dset_int, mh5_create_dset_real, mh5_create_dset_str, mh5_put_dset, mh5_get_dset, mh5_init_dset, &
          mh5_fetch_dset, mh5_resize_dset, mh5_get_dset_dims
! these are needed because assumed-size arguments match only one-dimensional arrays
public :: mh5_put_dset_array_int, mh5_put_dset_array_real, mh5_fetch_dset_array_real, mh5_fetch_dset_array_str


! generic interface for dataset creation

interface

  logical function mh5_is_hdf5 (name)
  implicit none
  character(len=*) :: name
  end function

  function mh5_create_file (filename) result (lu)
  implicit none
  character(len=*) :: filename
  integer :: lu
  end function

  function mh5_open_file_rw (filename) result (lu)
  implicit none
  character(len=*) :: filename
  integer :: lu
  end function

  function mh5_open_file_r (filename) result (lu)
  implicit none
  character(len=*) :: filename
  integer :: lu
  end function

  subroutine mh5_close_file (lu)
  implicit none
  integer :: lu
  end subroutine

! open/close group

  function mh5_open_group (lu, name) result (id)
  implicit none
  character(len=*) :: name
  integer :: lu, id
  end function

  subroutine mh5_close_group (id)
  implicit none
  integer :: id
  end subroutine

! open/close attribute

  logical function mh5_exists_attr (id, name)
  implicit none
  integer :: id
  character(len=*) :: name
  end function

  function mh5_open_attr (lu, attrname) result (attrid)
  implicit none
  integer :: lu
  character(len=*) :: attrname
  integer :: attrid
  end function

  subroutine mh5_close_attr (attrid)
  implicit none
  integer :: attrid
  end subroutine

! open/close dataset

  logical function mh5_exists_dset (id, name)
  implicit none
  integer :: id
  character(len=*) :: name
  end function

  function mh5_open_dset (lu, dsetname) result (dsetid)
  implicit none
  integer :: lu
  character(len=*) :: dsetname
  integer :: dsetid
  end function

  subroutine mh5_close_dset (dsetid)
  implicit none
  integer :: dsetid
  end subroutine

end interface

!=================
! ATTRIBUTES
!=================

interface mh5_create_attr_int

  function mh5_create_attr_scalar_int (lu, name) result (dset_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: dset_id
  end function

  function mh5_create_attr_array_int (lu, name, rank, dims) result (attr_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  integer :: attr_id
  end function

end interface mh5_create_attr_int

interface mh5_create_attr_real

  function mh5_create_attr_scalar_real (lu, name) result (dset_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: dset_id
  end function

  function mh5_create_attr_array_real (lu, name, rank, dims) result (attr_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  integer :: attr_id
  end function

end interface mh5_create_attr_real

interface mh5_create_attr_str

  function mh5_create_attr_scalar_str (lu, name) result (dset_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: dset_id
  end function

  function mh5_create_attr_array_str (lu, name, rank, dims, size) result (attr_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  integer :: size
  integer :: attr_id
  end function

end interface mh5_create_attr_str

interface mh5_put_attr

  subroutine mh5_put_attr_scalar_int (attr_id, value)
  implicit none
  integer :: attr_id
  integer :: value
  end subroutine

  subroutine mh5_put_attr_scalar_real (attr_id, value)
  implicit none
  integer :: attr_id
  real*8 :: value
  end subroutine

  subroutine mh5_put_attr_scalar_str (attr_id, value)
  implicit none
  integer :: attr_id
  character(len=*) :: value
  end subroutine

  subroutine mh5_put_attr_array_int (attr_id, buffer)
  implicit none
  integer :: attr_id
  integer :: buffer(*)
  end subroutine

    subroutine mh5_put_attr_array_real (attr_id, buffer)
  implicit none
  integer :: attr_id
  real*8 :: buffer(*)
  end subroutine

  subroutine mh5_put_attr_array_str (attr_id, buffer)
  implicit none
  integer :: attr_id
  character :: buffer(*)
  end subroutine

end interface mh5_put_attr

interface mh5_get_attr

  subroutine mh5_get_attr_scalar_int (attr_id, value)
  implicit none
  integer :: attr_id
  integer :: value
  end subroutine

  subroutine mh5_get_attr_scalar_real (attr_id, value)
  implicit none
  integer :: attr_id
  real*8 :: value
  end subroutine

  subroutine mh5_get_attr_scalar_str (attr_id, value)
  implicit none
  integer :: attr_id
  character(len=*) :: value
  end subroutine

  subroutine mh5_get_attr_array_int (attr_id, buffer)
  implicit none
  integer :: attr_id
  integer :: buffer(*)
  end subroutine

  subroutine mh5_get_attr_array_real (attr_id, buffer)
  implicit none
  integer :: attr_id
  real*8 :: buffer(*)
  end subroutine

  subroutine mh5_get_attr_array_str (attr_id, buffer)
  implicit none
  integer :: attr_id
  character :: buffer(*)
  end subroutine

end interface mh5_get_attr

interface mh5_init_attr

  subroutine mh5_init_attr_scalar_int (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: value
  end subroutine

  subroutine mh5_init_attr_scalar_real (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  real*8 :: value
  end subroutine

  subroutine mh5_init_attr_scalar_str (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  character(len=*) :: value
  end subroutine

  subroutine mh5_init_attr_array_int (lu, name, rank, dims, buffer)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  integer :: buffer(*)
  end subroutine

  subroutine mh5_init_attr_array_real (lu, name, rank, dims, buffer)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  real*8 :: buffer(*)
  end subroutine

  subroutine mh5_init_attr_array_str (lu, name, rank, dims, buffer, size)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  character :: buffer(*)
  integer :: size
  end subroutine

end interface mh5_init_attr

interface mh5_fetch_attr

  subroutine mh5_fetch_attr_scalar_int (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: value
  end subroutine

  subroutine mh5_fetch_attr_scalar_real (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  real*8 :: value
  end subroutine

  subroutine mh5_fetch_attr_scalar_str (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  character(len=*) :: value
  end subroutine

  subroutine mh5_fetch_attr_array_int (lu, name, buffer)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: buffer(*)
  end subroutine

  subroutine mh5_fetch_attr_array_real (lu, name, buffer)
  implicit none
  integer :: lu
  character(len=*) :: name
  real*8 :: buffer(*)
  end subroutine

  subroutine mh5_fetch_attr_array_str (lu, name, buffer)
  implicit none
  integer :: lu
  character(len=*) :: name
  character :: buffer(*)
  end subroutine

end interface mh5_fetch_attr

!=================
! DATASETS
!=================

interface mh5_create_dset_int

  function mh5_create_dset_scalar_int (lu, name) result (dset_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: dset_id
  end function

  function mh5_create_dset_array_int (lu, name, rank, dims, dyn) result (dset_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  logical, optional :: dyn
  integer :: dset_id
  end function

end interface mh5_create_dset_int

interface mh5_create_dset_real

  function mh5_create_dset_scalar_real (lu, name) result (dset_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: dset_id
  end function

  function mh5_create_dset_array_real (lu, name, rank, dims, dyn) result (dset_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  logical, optional :: dyn
  integer :: dset_id
  end function

end interface mh5_create_dset_real

interface mh5_create_dset_str

  function mh5_create_dset_scalar_str (lu, name) result (dset_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: dset_id
  end function

  function mh5_create_dset_array_str (lu, name, rank, dims, size, dyn) result (dset_id)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  integer :: size
  logical, optional :: dyn
  integer :: dset_id
  end function

end interface mh5_create_dset_str

interface mh5_put_dset

  subroutine mh5_put_dset_scalar_int (dset_id, value)
  implicit none
  integer :: dset_id
  integer :: value
  end subroutine

  subroutine mh5_put_dset_scalar_real (dset_id, value)
  implicit none
  integer :: dset_id
  real*8 :: value
  end subroutine

  subroutine mh5_put_dset_scalar_str (dset_id, value)
  implicit none
  integer :: dset_id
  character(len=*) :: value
  end subroutine

  subroutine mh5_put_dset_array_int (dset_id, buffer, exts, offs)
  implicit none
  integer :: dset_id
  integer :: buffer(*)
  integer, optional :: exts(*), offs(*)
  end subroutine

  subroutine mh5_put_dset_array_real (dset_id, buffer, exts, offs)
  implicit none
  integer :: dset_id
  real*8 :: buffer(*)
  integer, optional :: exts(*), offs(*)
  end subroutine

  subroutine mh5_put_dset_array_str (dset_id, buffer, exts, offs)
  implicit none
  integer :: dset_id
  character :: buffer(*)
  integer, optional :: exts(*), offs(*)
  end subroutine

end interface mh5_put_dset

interface mh5_get_dset

  subroutine mh5_get_dset_scalar_int (dset_id, value)
  implicit none
  integer :: dset_id
  integer :: value
  end subroutine

  subroutine mh5_get_dset_scalar_real (dset_id, value)
  implicit none
  integer :: dset_id
  real*8 :: value
  end subroutine

  subroutine mh5_get_dset_scalar_str (dset_id, value)
  implicit none
  integer :: dset_id
  character(len=*) :: value
  end subroutine

  subroutine mh5_get_dset_array_int (dset_id, buffer, exts, offs)
  implicit none
  integer :: dset_id
  integer :: buffer(*)
  integer, optional :: exts(*), offs(*)
  end subroutine

  subroutine mh5_get_dset_array_real (dset_id, buffer, exts, offs)
  implicit none
  integer :: dset_id
  real*8 :: buffer(*)
  integer, optional :: exts(*), offs(*)
  end subroutine

  subroutine mh5_get_dset_array_str (dset_id, buffer, exts, offs)
  implicit none
  integer :: dset_id
  character :: buffer(*)
  integer, optional :: exts(*), offs(*)
  end subroutine

end interface mh5_get_dset

interface mh5_init_dset

  subroutine mh5_init_dset_scalar_int (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: value
  end subroutine

  subroutine mh5_init_dset_scalar_real (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  real*8 :: value
  end subroutine

  subroutine mh5_init_dset_scalar_str (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  character(len=*) :: value
  end subroutine

  subroutine mh5_init_dset_array_int (lu, name, rank, dims, buffer, dyn)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  integer :: buffer(*)
  logical, optional :: dyn
  end subroutine

  subroutine mh5_init_dset_array_real (lu, name, rank, dims, buffer, dyn)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  real*8 :: buffer(*)
  logical, optional :: dyn
  end subroutine

  subroutine mh5_init_dset_array_str (lu, name, rank, dims, buffer, size, dyn)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: rank
  integer :: dims(*)
  character :: buffer(*)
  integer :: size
  logical, optional :: dyn
  end subroutine

end interface mh5_init_dset

interface mh5_fetch_dset

  subroutine mh5_fetch_dset_scalar_int (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: value
  end subroutine

  subroutine mh5_fetch_dset_scalar_real (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  real*8 :: value
  end subroutine

  subroutine mh5_fetch_dset_scalar_str (lu, name, value)
  implicit none
  integer :: lu
  character(len=*) :: name
  character(len=*) :: value
  end subroutine

  subroutine mh5_fetch_dset_array_int (lu, name, buffer, exts, offs)
  implicit none
  integer :: lu
  character(len=*) :: name
  integer :: buffer(*)
  integer, optional :: exts(*), offs(*)
  end subroutine

  subroutine mh5_fetch_dset_array_real (lu, name, buffer, exts, offs)
  implicit none
  integer :: lu
  character(len=*) :: name
  real*8 :: buffer(*)
  integer, optional :: exts(*), offs(*)
  end subroutine

  subroutine mh5_fetch_dset_array_str (lu, name, buffer, exts, offs)
  implicit none
  integer :: lu
  character(len=*) :: name
  character :: buffer(*)
  integer, optional :: exts(*), offs(*)
  end subroutine

end interface mh5_fetch_dset

interface mh5_resize_dset

  subroutine mh5_extend_dset_array (dset_id, dims)
  implicit none
  integer :: dset_id
  integer :: dims(*)
  end subroutine

end interface mh5_resize_dset

interface mh5_get_dset_dims

  subroutine mh5_get_dset_array_dims(dset_id, dims)
  implicit none
  integer :: dset_id
  integer :: dims(*)
  end subroutine

end interface mh5_get_dset_dims

end module mh5
