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
  use datasets_m
  use density_m
  use fourier_space_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use io_function_m
  use io_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use output_m
  use parser_m
  use profiling_m
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
    harmonic_spect_out

  type harmonic_spect_t !< Angular-resolved Harmonic Spectrum structure
    
    logical        :: calc = .false.         !< enable HS calculation 
    
    type(mesh_t), pointer :: mesh => NULL()  
    type(grid_t), pointer :: gr => NULL()  
    
    CMPLX, pointer :: spect(:,:,:,:) => NULL() !< The spectrum as a function of omega, theta and phi 
    integer        :: nw
    integer        :: nt 
    integer        :: np


    FLOAT          :: omega(1:3)            !< omega: 1 - min, 2 - max, 3 - delta 
    FLOAT          :: theta(1:3)            !< theta: 1 - min, 2 - max, 3 - delta 
    FLOAT          :: phi(1:3)              !< phi: 1 - min, 2 - max, 3 - delta 
    
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

    this%spect => null()
    this%mesh  => null()
    this%gr  => null()

    POP_SUB(harmonic_spect_nullify)
  end subroutine harmonic_spect_nullify

  ! ---------------------------------------------------------
  subroutine harmonic_spect_init(this, gr, dt )
    type(harmonic_spect_t), intent(out) :: this
    type(grid_t), target,   intent(in)    :: gr
    FLOAT,                  intent(in)    :: dt


    integer :: hs_flags
    character(len=256) :: lab, str
    integer :: n_block, ib, ii
    type(block_t) :: blk
    FLOAT         :: memory 

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

      write(str, '(a,i5)') 'Harmonic Spectrum'
      call messages_print_stress(stdout, trim(str))
     
    

      !%Variable HarmonicSpectrumGrid
      !%Type block
      !%Section Time-Dependent::HarmonicSpectrum
      !%Description
      !% Specify the grid that we want to use. 
      !% 
      !% <tt>%HarmonicSpectrumGrid 
      !% <br>&nbsp;&nbsp; 'omega' | min | max | spacing 
      !% <br>&nbsp;&nbsp; 'theta' | min | max | spacing  (radiants)
      !% <br>&nbsp;&nbsp; 'phi'   | min | max | spacing  (radiants)
      !% <br>%</tt>
      !%
      !%End
      if (parse_block(datasets_check('HarmonicSpectrumGrid'), blk) < 0)then
        call input_error('HarmonicSpectrumGrid')

      else 
        n_block = parse_block_n(blk)
      
        do ib = 1, n_block
          call parse_block_string(blk, ib-1, 0, lab)
          select case (lab)
            case ('omega')
              call parse_block_float(blk, ib-1, 1, this%omega(1), units_inp%energy)
              call parse_block_float(blk, ib-1, 2, this%omega(2), units_inp%energy)
              call parse_block_float(blk, ib-1, 3, this%omega(3), units_inp%energy)

            case ('theta')
              call parse_block_float(blk, ib-1, 1, this%theta(1), unit_one)
              call parse_block_float(blk, ib-1, 2, this%theta(2), unit_one)
              call parse_block_float(blk, ib-1, 3, this%theta(3), unit_one)

            case ('phi')
              call parse_block_float(blk, ib-1, 1, this%phi(1), unit_one)
              call parse_block_float(blk, ib-1, 2, this%phi(2), unit_one)
              call parse_block_float(blk, ib-1, 3, this%phi(3), unit_one)

            case default  
              call input_error('HarmonicSpectrumGrid')

          end select
        end do   
      end if
    
      this%nw = nint(abs((this%omega(2)-this%omega(1))/this%omega(3)))
      this%nt = nint(abs((this%theta(2)-this%theta(1))/this%theta(3)))
      this%np = nint(abs((this%phi(2)-this%phi(1))/this%phi(3)))
    
    
      memory =  M_TWO * REAL_PRECISION * dble(this%nw * this%nt *this%np)
      memory =  memory/CNST(1024.0)**2
    
      write(message(1), '(a)') "Grid:"
      call messages_info(1)
    
      write(message(1),'(3a)') '  Omega [', trim(units_abbrev(units_out%energy)), '] = ('
      do ii = 1, 3
        if(ii == 2) write(message(1), '(2a)') trim(message(1)), ','
        if(ii == 3) write(message(1), '(2a)') trim(message(1)), ';'
        write(message(1), '(a,f6.3)') trim(message(1)), units_from_atomic(units_out%energy, this%omega(ii))
      end do
      write(message(1), '(2a,i10)') trim(message(1)), ')   # points  = ',  this%nw

      write(message(2),'(3a)') '  Theta [', trim(units_abbrev(units_out%energy)), '] = ('
      do ii = 1, 3
        if(ii == 2) write(message(2), '(2a)') trim(message(2)), ','
        if(ii == 3) write(message(2), '(2a)') trim(message(2)), ';'
        write(message(2), '(a,f6.3)') trim(message(2)), units_from_atomic(units_out%energy, this%theta(ii))
      end do
      write(message(2), '(2a,i10)') trim(message(2)), ')   # points  = ',  this%nt

      write(message(3),'(3a)') '  Phi   [', trim(units_abbrev(units_out%energy)), '] = ('
      do ii = 1, 3
        if(ii == 2) write(message(3), '(2a)') trim(message(3)), ','
        if(ii == 3) write(message(3), '(2a)') trim(message(3)), ';'
        write(message(3), '(a,f6.3)') trim(message(3)), units_from_atomic(units_out%energy, this%phi(ii))
      end do
      write(message(3), '(2a,i10)') trim(message(3)), ')   # points  = ',  this%np
    
      write(message(4),'(a,i10)') '  # total mesh (nw * nt * np) = ', this%np * this%nw * this%nt * 3
      write(message(5),'(a,f10.1,a)')'  memory =', memory, ' [Mb]'

      call messages_info(5)
    

      !ALLOCATION 
      SAFE_ALLOCATE(this%spect(1:this%nw, 1:this%nt, 1:this%np, 1:3))
      this%spect = M_z0



     !Footer HarmonicSpectrum info
      call messages_print_stress(stdout)
    end if 

    POP_SUB(harmonic_spect_init)
  end subroutine harmonic_spect_init


  ! ---------------------------------------------------------
  subroutine harmonic_spect_end(this)
    type(harmonic_spect_t), intent(inout) :: this

    PUSH_SUB(harmonic_spect_end)
    if(this%calc) then
      SAFE_DEALLOCATE_P(this%spect)  
      call harmonic_spect_nullify(this)
    end if
    POP_SUB(harmonic_spect_end)
  end subroutine harmonic_spect_end


  ! ---------------------------------------------------------
  subroutine harmonic_spect_calc(this, st, ii)
    type(harmonic_spect_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st
    integer,             intent(in)    :: ii

    FLOAT, allocatable :: J(:,:,:), nnJx(:,:), nx(:)
    integer :: iw, it, ip, ir, dir
    FLOAT   :: NN(1:3), omega, theta, phi, t
    CMPLX   :: nnj(1:3) 
  


    PUSH_SUB(harmonic_spect_calc)
    
    ASSERT(this%calc)
    ASSERT(states_are_complex(st))

    SAFE_ALLOCATE(J(1:this%mesh%np, 1:3, 1:st%d%nspin))
    SAFE_ALLOCATE(nnJx(1:this%mesh%np, 1:3))
    SAFE_ALLOCATE(nx(1:this%mesh%np))
    
    call states_calc_quantities(this%gr%der, st, paramagnetic_current = J)

    t = ii * this%dt

    do iw = 1, this%nw
      omega = this%omega(1) + iw * this%omega(3)
      do it = 1, this%nt
        theta = this%theta(1) + it * this%theta(3)
        do ip = 1, this%np 
          phi = this%phi(1) + ip * this%phi(3)
        
          
          NN(1) = sin(theta) * cos(phi) 
          NN(2) = sin(theta) * sin(phi)
          NN(3) = cos(theta)
 
          
          do ir = 1, this%mesh%np
            nx(ir) = ddot_product(NN, mesh_x_global(this%mesh, ir))
            nnJx(:,:) = M_ZERO 
            do dir = 1, st%d%nspin         
              nnJx(ir,:) = nnJx(ir,:) + dcross_product(NN, &
                                              dcross_product(NN, J(ir, :, dir) ) )
            end do
          end do
          
          do dir= 1, 3
            nnJ(dir) = zmf_integrate(this%mesh, exp( M_zI * omega * nx(:)/P_C) * nnJx(:,dir))
          end do
                     
          this%spect(iw, it, ip, :) = this%spect(iw, it, ip, :) + exp( - M_zI * omega * t ) * nnJ(:)
            
        end do
      end do
    end do


    SAFE_DEALLOCATE_A(J)
    SAFE_DEALLOCATE_A(nnJx)
    SAFE_DEALLOCATE_A(nx)
    
    POP_SUB(harmonic_spect_calc)
  end subroutine harmonic_spect_calc


  ! ---------------------------------------------------------
  subroutine harmonic_spect_out(this, outp)
    type(harmonic_spect_t),      intent(inout) :: this
    type(output_t),   intent(in)    :: outp

    PUSH_SUB(harmonic_spect_out)

    ASSERT(this%calc)
    



    POP_SUB(harmonic_spect_out)
  end subroutine harmonic_spect_out

  ! ---------------------------------------------------------
  subroutine harmonic_spect_restart_write(this)
    type(harmonic_spect_t),    intent(in) :: this

    PUSH_SUB(harmonic_spect_restart_write)

    ASSERT(this%calc)


    POP_SUB(harmonic_spect_restart_write)
  end subroutine harmonic_spect_restart_write

  ! ---------------------------------------------------------
  subroutine harmonic_spect_restart_read(this)
    type(harmonic_spect_t),    intent(inout) :: this

    PUSH_SUB(harmonic_spect_restart_read)

    ASSERT(this%calc)
 

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


end module harmonic_spect_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
