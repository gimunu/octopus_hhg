!! Copyright (C) 2006-2011 M. Marques, U. De Giovannini
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: this.F90 9210 2012-07-19 17:10:36Z umberto $

#include "global.h"

module harmonic_spect_m
  use batch_m
  use comm_m
  use cube_m
  use cube_function_m
  use datasets_m
  use density_m
  use fft_m
  use fourier_space_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use io_function_m
  use io_m
  use io_binary_m
  use math_m
  use mesh_m
  use mesh_cube_parallel_map_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use output_m
  use parser_m
  use profiling_m
  use qshepmod_m
  use restart_m
  use simul_box_m
  use states_io_m
  use states_m
  use string_m
  use system_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                        &
    harmonic_spect_t,              &
    harmonic_spect_init,           &
    harmonic_spect_end,            &
    harmonic_spect_calc,           &
    harmonic_spect_checkpoint,     &
    harmonic_spect_restart_read

  type harmonic_spect_t !< Angular-resolved Harmonic Spectrum structure
    
    logical        :: calc = .false.         !< enable HS calculation 
    
    type(mesh_t), pointer :: mesh => NULL()  
    type(grid_t), pointer :: gr => NULL()  
    
    type(cube_t)   :: cube
    type(mesh_cube_parallel_map_t) :: mesh_cube_map  !< The parallel map
    
    type(cube_function_t) :: Jkint(3) !< Time integrated FS density current
        
    FLOAT          :: dt
    integer        :: what                  !< what do we want to calculate
    
  end type harmonic_spect_t

  integer, parameter ::     &
    HS_NONE = 0, &
    HS_AR   = 2, &
    HS_TOT  = 4


contains

  ! ---------------------------------------------------------
  subroutine harmonic_spect_nullify(this)
    type(harmonic_spect_t), intent(out) :: this
    PUSH_SUB(harmonic_spect_nullify)

    this%mesh  => null()
    this%gr  => null()

    POP_SUB(harmonic_spect_nullify)
  end subroutine harmonic_spect_nullify

  ! ---------------------------------------------------------
  subroutine harmonic_spect_init(this, gr, dt)
    type(harmonic_spect_t), intent(out) :: this
    type(grid_t), target,   intent(in)    :: gr
    FLOAT,                  intent(in)    :: dt


    integer :: hs_flags
    character(len=256) :: lab, str
    integer :: n_block, ib, ii
    type(block_t) :: blk
    FLOAT         :: memory 
    integer :: i, box(MAX_DIM)

    PUSH_SUB(harmonic_spect_init)
    
    call harmonic_spect_nullify(this)


    this%calc = .false.

    !%Variable HarmonicSpectrum 
    !%Type flag
    !%Default no
    !%Section Time-Dependent::HarmonicSpectrum 
    !%Description
    !%This variable controls the method used for the calculation of
    !%the harmonic spectrum from the current density. 
    !% Another possiblity is to calculate the from the multipoles or acceleration 
    !% with the <tt>oct-harmonic-spectrum</tt> utility.
    !% <tt>HarmonicSpectrum = angular_resolved + total</tt>
    !%Option none 0
    !% The harmonic spectrum is not calculated. This is the default.
    !%Option angular_resolved 2
    !% Angular resolved harmonic spectrum calculated from the current density.
    !%Option total 4
    !% The angle integrated quantity.
    !%End

    call parse_integer(datasets_check('HarmonicSpectrum'), HS_NONE, hs_flags)
    if(.not.varinfo_valid_option('HarmonicSpectrum', hs_flags, is_flag = .true.)) then
      call input_error('HarmonicSpectrum')
    end if
    
    this%what = hs_flags
    
    !Header HarmonicSpectrum info
    if( hs_flags > 0) then
      !init values
      this%calc = .true.       
      this%gr => gr
      this%mesh => gr%mesh
      this%dt = dt
      box(:) = this%mesh%idx%ll(:)

      write(str, '(a,i5)') 'Harmonic Spectrum'
      call messages_print_stress(stdout, trim(str))
     
    

!       !%Variable HarmonicSpectrumGrid
!       !%Type block
!       !%Section Time-Dependent::HarmonicSpectrum
!       !%Description
!       !% Specify the grid that we want to use. 
!       !% 
!       !% <tt>%HarmonicSpectrumGrid 
!       !% <br>&nbsp;&nbsp; 'omega' | min | max | spacing 
!       !% <br>&nbsp;&nbsp; 'theta' | min | max | spacing  (radiants)
!       !% <br>&nbsp;&nbsp; 'phi'   | min | max | spacing  (radiants)
!       !% <br>%</tt>
!       !%
!       !%End
!       if (parse_block(datasets_check('HarmonicSpectrumGrid'), blk) < 0)then
!         call input_error('HarmonicSpectrumGrid')
! 
!       else 
!         n_block = parse_block_n(blk)
!       
!         do ib = 1, n_block
!           call parse_block_string(blk, ib-1, 0, lab)
!           select case (lab)
!             case ('omega')
!               call parse_block_float(blk, ib-1, 1, this%omega(1), units_inp%energy)
!               call parse_block_float(blk, ib-1, 2, this%omega(2), units_inp%energy)
!               call parse_block_float(blk, ib-1, 3, this%omega(3), units_inp%energy)
! 
!             case ('theta')
!               call parse_block_float(blk, ib-1, 1, this%theta(1), unit_one)
!               call parse_block_float(blk, ib-1, 2, this%theta(2), unit_one)
!               call parse_block_float(blk, ib-1, 3, this%theta(3), unit_one)
! 
!             case ('phi')
!               call parse_block_float(blk, ib-1, 1, this%phi(1), unit_one)
!               call parse_block_float(blk, ib-1, 2, this%phi(2), unit_one)
!               call parse_block_float(blk, ib-1, 3, this%phi(3), unit_one)
! 
!             case default  
!               call input_error('HarmonicSpectrumGrid')
! 
!           end select
!         end do   
!       end if
!     
!       this%nw = nint(abs((this%omega(2)-this%omega(1))/this%omega(3)))
!       this%nt = nint(abs((this%theta(2)-this%theta(1))/this%theta(3)))
!       this%np = nint(abs((this%phi(2)-this%phi(1))/this%phi(3)))
!     
!     
!       memory =  M_TWO * REAL_PRECISION * dble(this%nw * this%nt *this%np)
!       memory =  memory/CNST(1024.0)**2
!     
!       write(message(1), '(a)') "Grid:"
!       call messages_info(1)
!     
!       write(message(1),'(3a)') '  Omega [', trim(units_abbrev(units_out%energy)), ']   = ('
!       do ii = 1, 3
!         if(ii == 2) write(message(1), '(2a)') trim(message(1)), ','
!         if(ii == 3) write(message(1), '(2a)') trim(message(1)), ';'
!         write(message(1), '(a,f6.3)') trim(message(1)), units_from_atomic(units_out%energy, this%omega(ii))
!       end do
!       write(message(1), '(2a,i10)') trim(message(1)), ')   # points  = ',  this%nw
! 
!       write(message(2),'(3a)') '  Theta [', 'Rad', '] = ('
!       do ii = 1, 3
!         if(ii == 2) write(message(2), '(2a)') trim(message(2)), ','
!         if(ii == 3) write(message(2), '(2a)') trim(message(2)), ';'
!         write(message(2), '(a,f6.3)') trim(message(2)), units_from_atomic(units_out%energy, this%theta(ii))
!       end do
!       write(message(2), '(2a,i10)') trim(message(2)), ')   # points  = ',  this%nt
! 
!       write(message(3),'(3a)') '  Phi   [', 'Rad', '] = ('
!       do ii = 1, 3
!         if(ii == 2) write(message(3), '(2a)') trim(message(3)), ','
!         if(ii == 3) write(message(3), '(2a)') trim(message(3)), ';'
!         write(message(3), '(a,f6.3)') trim(message(3)), units_from_atomic(units_out%energy, this%phi(ii))
!       end do
!       write(message(3), '(2a,i10)') trim(message(3)), ')   # points  = ',  this%np
!     
!       write(message(4),'(a,i10)') '  # total mesh (nw * nt * np) = ', this%np * this%nw * this%nt * 3
!       write(message(5),'(a,f10.1,a)')'  memory =', memory, ' [Mb]'
! 
!       call messages_info(5)
    

      !ALLOCATION 

      call cube_init(this%cube, box, this%mesh%sb, fft_type = FFT_REAL, verbose = .true.)
      call cube_init_fs_coords(this%cube, this%mesh%spacing(1:3), this%mesh%sb%dim)
      if ( this%mesh%parallel_in_domains .and. this%cube%parallel_in_domains) then
        call mesh_cube_parallel_map_init(this%mesh_cube_map, this%mesh, this%cube)
      end if      
      
      
      
      
      do i = 1, 3
      	call cube_function_null(this%Jkint(i))    
      	call cube_function_alloc_FS(this%cube, this%Jkint(i), force_alloc = .true.)   
        this%Jkint(i)%Fs = M_z0     
      end do 


     !Footer HarmonicSpectrum info
      call messages_print_stress(stdout)
    end if 

    POP_SUB(harmonic_spect_init)
  end subroutine harmonic_spect_init


  ! ---------------------------------------------------------
  subroutine harmonic_spect_end(this)
    type(harmonic_spect_t), intent(inout) :: this
    
    integer :: i

    PUSH_SUB(harmonic_spect_end)
    if(this%calc) then
      do i = 1, 3
        call cube_function_free_fs(this%cube, this%Jkint(i))        
      end do

      call cube_end(this%cube)   
      
      call harmonic_spect_nullify(this)
    end if
    POP_SUB(harmonic_spect_end)
  end subroutine harmonic_spect_end


  ! ---------------------------------------------------------
  subroutine harmonic_spect_calc(this, st, ii)
    type(harmonic_spect_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st
    integer,             intent(in)    :: ii

    FLOAT, allocatable :: Js(:,:,:), J(:,:)
    integer :: iw, it, ip, ir, dir
    FLOAT   :: NN(1:3), omega, theta, phi, t
    CMPLX   :: nnj(1:3) 
    type(cube_function_t) :: cfJ
  


    PUSH_SUB(harmonic_spect_calc)
    
    ASSERT(this%calc)
    ASSERT(states_are_complex(st))

  	call cube_function_null(cfJ)    
  	call dcube_function_alloc_RS(this%cube, cfJ)   
  	call  cube_function_alloc_FS(this%cube, cfJ)   

    SAFE_ALLOCATE(J(1:this%mesh%np, 1:3))
    
    SAFE_ALLOCATE(Js(1:this%mesh%np, 1:3, 1:st%d%nspin))      
    call states_calc_quantities(this%gr%der, st, paramagnetic_current = Js )
    J(:,:) = sum(Js(:,:,1:st%d%nspin)) ! sum up all the spin channels
    SAFE_DEALLOCATE_A(Js)

    t = ii * this%dt
    
    do dir= 1, 3
      
      if (this%cube%parallel_in_domains) then
        call dmesh_to_cube_parallel(this%mesh, J(:,dir), this%cube, cfJ, this%mesh_cube_map)
      else
        if(this%mesh%parallel_in_domains) then
          call dmesh_to_cube(this%mesh, J(:,dir), this%cube, cfJ, local = .true.)
        else 
          call dmesh_to_cube(this%mesh, J(:,dir), this%cube, cfJ)
        end if
      end if
      
      call dcube_function_rs2fs(this%cube, cfJ)
      
      call apply_phase(cfJ)
      
      this%Jkint(dir)%FS(:,:,:) = this%Jkint(dir)%FS(:,:,:) + cfJ%FS(:,:,:)
    end do

   


    call dcube_function_free_RS(this%cube, cfJ)
    call  cube_function_free_FS(this%cube, cfJ)
    SAFE_DEALLOCATE_A(J)
    
    POP_SUB(harmonic_spect_calc)
    
    contains 
    
    subroutine apply_phase(cfJ)
      type(cube_function_t), intent(inout) :: cfJ
      
      integer :: ix,iy,iz
      FLOAT :: KK(3), K
      
      do ix = 1, this%cube%fs_n(1)
        KK(1) = this%cube%k(ix,1)
        do iy = 1, this%cube%fs_n(2)
          KK(2) = this%cube%k(iy,2)
          do iz = 1, this%cube%fs_n(3)
            KK(3) = this%cube%k(iz,3)
            
            K = sqrt(KK(1)**2 + KK(2)**2 + KK(3)**2)
            cfJ%FS(ix, iy, iz) = cfJ%FS(ix, iy, iz) *  this%dt * exp(M_zI * this%dt * K * P_C)
            
          end do
        end do
      end do
      
      
    end subroutine apply_phase
    
  end subroutine harmonic_spect_calc

  ! ---------------------------------------------------------
  subroutine harmonic_spect_checkpoint(this)
    type(harmonic_spect_t), intent(in) :: this
    
    PUSH_SUB(harmonic_spect_checkpoint)
    
    if(this%calc) then
      call harmonic_spect_restart_write(this)
      call harmonic_spect_out(this)
    end if
    
    POP_SUB(harmonic_spect_checkpoint)
  end subroutine harmonic_spect_checkpoint

  ! ---------------------------------------------------------
  ! g(k) = |n x (n x J(k))|^2
  ! ---------------------------------------------------------
  subroutine harmonic_spect_gk(this, gk)
    type(harmonic_spect_t), intent(in) :: this
    FLOAT,                  intent(out):: gk(:,:,:)

    integer :: ix,iy,iz
    FLOAT :: KK(3), K, N(3)
    
    PUSH_SUB(harmonic_spect_gk)

    gk(:,:,:) = M_ZERO
          
    do ix = 1, this%cube%fs_n(1)
      KK(1) = this%cube%k(ix,1)
      do iy = 1, this%cube%fs_n(2)
        KK(2) = this%cube%k(iy,2)
        do iz = 1, this%cube%fs_n(3)
          KK(3) = this%cube%k(iz,3)
            
          K = sqrt(KK(1)**2 + KK(2)**2 + KK(3)**2)
          N(:) = KK(:)/K
          
          gk(ix, iy, iz) = abs(this%Jkint(1)%FS(ix, iy, iz))**2 &
                         + abs(this%Jkint(2)%FS(ix, iy, iz))**2 &
                         + abs(this%Jkint(3)%FS(ix, iy, iz))**2
            
          gk(ix, iy, iz) = gk(ix, iy, iz) &
                         - abs(N(1) * this%Jkint(1)%FS(ix, iy, iz))**2 &
                         - abs(N(2) * this%Jkint(2)%FS(ix, iy, iz))**2 &
                         - abs(N(3) * this%Jkint(3)%FS(ix, iy, iz))**2 
            
          gk(ix, iy, iz) = gk(ix, iy, iz) * K**2 
            
        end do
      end do
    end do
    
    gk = gk/ (CNST(4.0) * M_PI**2 * P_C)
    

    POP_SUB(harmonic_spect_gk)  
  end subroutine harmonic_spect_gk  

  ! ---------------------------------------------------------
  subroutine harmonic_spect_out(this)
    type(harmonic_spect_t),      intent(in) :: this

    
    FLOAT, allocatable :: gk(:,:,:)
    integer :: ix,iy,iz
    FLOAT :: KK(3), K, N(3)
    
    
    PUSH_SUB(harmonic_spect_out)

    ASSERT(this%calc)
    
    SAFE_ALLOCATE(gk(1:this%cube%fs_n(1),1:this%cube%fs_n(2),1:this%cube%fs_n(3)))
    
    call harmonic_spect_gk(this, gk)
    
    call harmonic_spect_total(this,"td.general/harmonic_spect.tot", gk)
    
    SAFE_DEALLOCATE_A(gk)

    POP_SUB(harmonic_spect_out)
  end subroutine harmonic_spect_out

  ! ---------------------------------------------------------
  subroutine harmonic_spect_restart_write(this)
    type(harmonic_spect_t),    intent(in) :: this

    character(len=256) :: filename, path, number
    integer :: dir, ierr, npoints
    
    PUSH_SUB(harmonic_spect_restart_write)

    ASSERT(this%calc)
    
    path = trim(restart_dir)//'td/'
    
    npoints = this%cube%fs_n(1)*this%cube%fs_n(2)*this%cube%fs_n(3)
    
    do dir = 1, 3
      write(number,'(i1.1)') dir
      filename = trim(path)//'hs_spect.'//number
      call io_binary_write(filename, npoints, this%Jkint(dir)%FS(:,:,:), ierr)    

      if(ierr > 0) then
        message(1) = "Failed to write file "//trim(filename)
        call messages_fatal(1)
      end if
      
    end do


    POP_SUB(harmonic_spect_restart_write)
  end subroutine harmonic_spect_restart_write

  ! ---------------------------------------------------------
  subroutine harmonic_spect_restart_read(this)
    type(harmonic_spect_t),    intent(inout) :: this

    character(len=256) :: filename, path, number
    integer :: dir, ierr, npoints

    PUSH_SUB(harmonic_spect_restart_read)

    ASSERT(this%calc)
 
    path = trim(restart_dir)//'td/'
    
    npoints = this%cube%fs_n(1)*this%cube%fs_n(2)*this%cube%fs_n(3)
    
    do dir = 1, 3
      write(number,'(i1.1)') dir
      filename = trim(path)//'hs_spect.'//number      
      call io_binary_read(filename, npoints, this%Jkint(dir)%FS(:,:,:), ierr)    

      if(ierr > 0) then
        message(1) = "Failed to write file "//trim(filename)
        call messages_fatal(1)
      end if
      
    end do


    POP_SUB(harmonic_spect_restart_read)
  end subroutine harmonic_spect_restart_read


  ! ---------------------------------------------------------
  subroutine harmonic_spect_init_write(this, mesh, st)
    type(harmonic_spect_t),    intent(in)  :: this
    type(mesh_t),   intent(in)  :: mesh
    type(states_t), intent(in)  :: st


    PUSH_SUB(harmonic_spect_init_write)

 

    POP_SUB(harmonic_spect_init_write)
  end subroutine harmonic_spect_init_write

  ! ---------------------------------------------------------
  subroutine harmonic_spect_total(this, file, gk, interpolate)
    type(harmonic_spect_t), intent(in) :: this
    character(len=*),       intent(in) :: file
    FLOAT,                  intent(in) :: gk(:,:,:)
    logical, optional,      intent(in) :: interpolate

    integer :: ist, ik, ii, ix, iy, iz, iunit,idim
    FLOAT ::  KK(3),vec

    integer :: nn
    FLOAT  :: step, DE
    FLOAT, allocatable :: npoints(:), out(:)

    ! needed for interpolation in 2D and 3D 
    FLOAT, pointer :: cube_f(:)
    type(qshep_t) :: interp

    FLOAT :: Dtheta, Dphi, theta, phi,EE , Emax, Edir(3)
    integer :: np, Ntheta, Nphi, ith, iph
    integer :: ll(3), dim




    PUSH_SUB(harmonic_spect_total)


    ll = 1
    step = M_ZERO
    dim = this%mesh%sb%dim
    do ii = 1, dim
      ll(ii) = size(gk,ii) 
      step= step + abs(this%cube%k(2, ii) - this%cube%k(1, ii))**2 
      Edir(ii) = maxval(abs(this%cube%k(1:this%cube%fs_n(ii), ii)))**2/M_TWO
    end do
    step = step /M_TWO
    
    
    Emax = maxval(Edir(1:dim))

!     print *,"Emax", Emax, "DE", step

    nn  = int(Emax/step)


    Ntheta = 360
    Dtheta = M_TWO*M_PI/Ntheta

    Nphi = 180
    Dphi = M_PI/Nphi

    SAFE_ALLOCATE(out(1:nn))
    out = M_ZERO

    SAFE_ALLOCATE(npoints(1:nn))
    npoints = M_ZERO

    !in 1D we do not interpolate 
    if ( (.not. optional_default(interpolate, .false.))  ) then 

      do ix = 1,ll(1)
        KK(1) = this%cube%k(ix,1)
        do iy = 1, ll(2)
          KK(2) = this%cube%k(iy,2)
          do iz = 1, ll(3)
            KK(3) = this%cube%k(iz,3)

            if(KK(1).ne.0 .or. KK(2).ne.0 .or. KK(3).ne.0) then
              ! the power spectrum
              vec = sum(KK(1:dim)**2) / M_TWO
              ii = int(vec / step) + 1

              if(ii <= nn) then

                out(ii) = out(ii) + gk(ix,iy,iz)
                npoints(ii) = npoints(ii) + M_ONE

              end if
            end if

          end do
        end do
      end do
      
      

    end if
      ! Interpolate the output
 !    else
 
!       call PES_mask_interpolator_init(PESK, Lk, dim, cube_f, interp)
! 
!       select case(dim)
!       case(2)
! 
!         do ii = 1, nn
!           EE = (ii-1)*step
!           do ith = 0, Ntheta
!             theta = ith*Dtheta
!             KK(1) = sqrt(2*EE)*cos(theta) 
!             KK(2) = sqrt(2*EE)*sin(theta)
!             out(ii) = out(ii) + qshep_interpolate(interp, cube_f, KK(1:2))
!           end do
!         end do
! 
!         out = out * Dtheta
! 
!       case(3)
! 
!         do ii = 1, nn
!           EE = (ii-1)*step
!           do ith = 0, Ntheta
!             theta = ith * Dtheta
!             do iph = 0, Nphi
!               phi = iph * Dphi
! 
!               KK(2) = sqrt(M_TWO*EE)*sin(phi)*cos(theta) 
!               KK(3) = sqrt(M_TWO*EE)*sin(phi)*sin(theta)
!               KK(1) = sqrt(M_TWO*EE)*cos(phi)
! 
!               out(ii) = out(ii) + qshep_interpolate(interp, cube_f, KK(1:3))
!             end do
!           end do
!           out(ii) = out(ii) * sqrt(M_TWO*EE)    
!         end do
! 
!         out = out * Dtheta * Dphi
! 
!       end select
! 
!       call PES_mask_interpolator_end(cube_f, interp)
! 
!     end if


!     if (interpolate) then 
!       call PES_mask_write_power_total(file, step, out)
!     else 
!       call PES_mask_write_power_total(file, step, out, npoints)
!     end if


    iunit = io_open(file, action='write')

    !!Header
    write(iunit, '(a)') '##################################################'
    write(iunit, '(a1,a18,a18)') '#', str_center("E", 18), str_center("S(E)", 18)
    write(iunit, '(a1,a18,a18)') &
      '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 18), &
      str_center('[1/' //trim(units_abbrev(units_out%energy))//']', 18)
    write(iunit, '(a)') '##################################################'
 

    do ii = 1, nn
!       if(present(npoints)) then

        if(npoints(ii) > 0) then
          write(iunit, '(es19.12,2x,es19.12,2x,es19.12)')  units_from_atomic(units_out%energy, (ii - 1) * step), out(ii), npoints(ii)
        end if

!       else
!         write(iunit, '(es19.12,2x,es19.12)')  units_from_atomic(units_out%energy, (ii - 1) * step), pes(ii)
!       end if
    end do
  
    call io_close(iunit)



    SAFE_DEALLOCATE_A(out)
    SAFE_DEALLOCATE_A(npoints)

    POP_SUB(harmonic_spect_total)


  end subroutine harmonic_spect_total



end module harmonic_spect_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
