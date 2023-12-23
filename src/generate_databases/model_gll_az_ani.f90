!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
!
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!
!--------------------------------------------------------------------------------------------------

  subroutine model_gll_az_ani(myrank,nspec,LOCAL_PATH)

  use constants, only: NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IMAIN,MAX_STRING_LEN,CUSTOM_REAL,IIN,IANISOTROPY_MODEL1

  use generate_databases_par, only: ATTENUATION

  use create_regions_mesh_ext_par, only: rhostore,rho_vp,rho_vs,&
                                        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                        c36store,c44store,c45store,c46store,c55store,c56store,c66store

  use shared_parameters, only: ADIOS_FOR_MESH,HDF5_ENABLED

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH
  integer :: ispec,i,j,k

  ! local parameters
  real, dimension(:,:,:,:),allocatable :: gs_read, gc_read
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname_lp,filename
  real(kind=CUSTOM_REAL) :: vp,vs,rho,c11,c12,c13,c14,c15,c16,c22,c23,c24,&
                            c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,&
                            d44, d45, d55

  ! select routine for file i/o format
  if (ADIOS_FOR_MESH) then
    ! ADIOS
    call model_gll_adios(myrank,nspec,LOCAL_PATH)
    ! all done
    return
  else if (HDF5_ENABLED) then
    ! not implemented yet
    stop 'HDF5_ENABLED not supported yet for model_gll() routine, please return without flag...'
  else
    ! default binary
    ! implemented here below, continue
    continue
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using GLL model from: ',trim(LOCAL_PATH)
  endif

  ! processors name
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)// '/' //'proc',myrank,'_'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! if only vp structure is available (as is often the case in exploration seismology),
  !!! use lines for vp only

  ! gs
  allocate(gs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 648')
  if (ier /= 0) stop 'error allocating array gs_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: gc.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'gs_read.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading gc.bin file'
  endif

  read(IIN) gs_read
  close(IIN)

  ! gc
  allocate(gc_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 649')
  if (ier /= 0) stop 'error allocating array gc_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: gc_read.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'gc_read.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading gc_read.bin file'
  endif

  read(IIN) gc_read
  close(IIN)

  call synchronize_all()

  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! calculate vp, vs, rho
          rho = rhostore(i,j,k,ispec)
          vp = rho_vp(i,j,k,ispec)/rho
          vs = rho_vs(i,j,k,ispec)/rho

          ! calculate iso terms
          call model_aniso(IANISOTROPY_MODEL1,rho,vp,vs, &
                          c11,c12,c13,c14,c15,c16, &
                          c22,c23,c24,c25,c26,c33, &
                          c34,c35,c36,c44,c45,c46,c55,c56,c66)
          
          ! calculate ani terms
          d44 = rho*vs*vs - gc_read(i,j,k,ispec)
          d55 = rho*vs*vs + gc_read(i,j,k,ispec)
          d45 = -gs_read(i,j,k,ispec)
          c44 = d55
          c45 = - d45
          c55 = d44

          ! overwrite elastic tensor
          c11store(i,j,k,ispec) = c11
          c12store(i,j,k,ispec) = c12
          c13store(i,j,k,ispec) = c13
          c14store(i,j,k,ispec) = c14
          c15store(i,j,k,ispec) = c15
          c16store(i,j,k,ispec) = c16
          c22store(i,j,k,ispec) = c22
          c23store(i,j,k,ispec) = c23
          c24store(i,j,k,ispec) = c24
          c25store(i,j,k,ispec) = c25
          c26store(i,j,k,ispec) = c26
          c33store(i,j,k,ispec) = c33
          c34store(i,j,k,ispec) = c34
          c35store(i,j,k,ispec) = c35
          c36store(i,j,k,ispec) = c36
          c44store(i,j,k,ispec) = c44
          c45store(i,j,k,ispec) = c45
          c46store(i,j,k,ispec) = c46
          c55store(i,j,k,ispec) = c55
          c56store(i,j,k,ispec) = c56
          c66store(i,j,k,ispec) = c66
        enddo
      enddo
    enddo
  enddo



  end subroutine model_gll_az_ani

