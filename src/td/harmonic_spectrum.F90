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
  use nfft_m
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
    CMPLX, pointer        :: J(:,:,:) => NULL()!< the curent density at different timed used to calcualte dJ/dt     
    
    type(cube_function_t) :: cftmp(3)    
        
    FLOAT          :: dt
    integer        :: what                  !< what do we want to calculate
    integer        :: how                   !< how to do the calculation
    
    type(fft_t)    :: fft
    
  end type harmonic_spect_t

  integer, parameter ::     &
    HS_NONE = 0, &
    HS_AR   = 2, &
    HS_TOT  = 4

  integer, parameter ::     &
    HS_FROM_J      = 1, &
    HS_FROM_DJDT   = 2


contains

  ! ---------------------------------------------------------
  subroutine harmonic_spect_nullify(this)
    type(harmonic_spect_t), intent(out) :: this
    PUSH_SUB(harmonic_spect_nullify)

    this%mesh  => null()
    this%gr    => null()
    this%J     => null()

    POP_SUB(harmonic_spect_nullify)
  end subroutine harmonic_spect_nullify

  ! ---------------------------------------------------------
  subroutine harmonic_spect_init(this, gr, dt)
    type(harmonic_spect_t), intent(out) :: this
    type(grid_t), target,   intent(in)    :: gr
    FLOAT,                  intent(in)    :: dt


    integer :: hs_flags
    character(len=256) :: lab, str
    integer :: ii, dim, idim
    type(block_t) :: blk
    FLOAT         :: memory, omegaMax, omegaMin, temp(3), WMax, Sfactor
    integer :: i, box(MAX_DIM)
    logical :: use_nfft

    !FOR NFFT
    FLOAT, allocatable :: XX(:,:)
    logical :: optimize(3)
    integer :: optimize_parity(3)

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
    !%Option yes 1
    !% Calcualte all.
    !%Option integrated_current 2
    !% Angular resolved harmonic spectrum calculated from the current density.
    !%Option total 4
    !% The angle integrated quantity.
    !%Option all_hs 8
    !% Calcualte all.
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
      dim = this%mesh%sb%dim

      write(str, '(a,i5)') 'Harmonic Spectrum'
      call messages_print_stress(stdout, trim(str))
     
      !%Variable HarmonicSpectrumHow 
      !%Type integer
      !%Default use_J
      !%Section Time-Dependent::HarmonicSpectrumHow 
      !%Description
      !%Option use_J 1
      !% n x (n x J).
      !%Option use_dJdt 2
      !% n x (n x dJ/dt).
      !%End
    
      call parse_integer(datasets_check('HarmonicSpectrumHow'), HS_FROM_J, this%how)
      if(.not.varinfo_valid_option('HarmonicSpectrumHow', this%how)) then
        call input_error('HarmonicSpectrumHow')
      end if
      call messages_print_var_option(stdout, "HarmonicSpectrumHow", this%how)
    
    
      !%Variable HarmonicSpectrumOmegaMax 
      !%Type float
      !%Section Time-Dependent::HarmonicSpectrum 
      !%Description
      !% Set a target maximum harmonic frequency.      
      !% Requires NFFT libraries. 
      !%End
      call parse_float(datasets_check('HarmonicSpectrumOmegaMax'), CNST(-1.0), Wmax)
      Wmax = units_to_atomic(units_inp%energy, Wmax)
      Sfactor = - M_ONE
      use_nfft = .false.
      if(Wmax .gt. M_ZERO) use_nfft = .true. 
#ifndef HAVE_NFFT
      if(use_nfft) then
        write(message(1),'(a)')'HarmonicSpectrumOmegaMax needs NFFT library.'
        call messages_fatal(1)
      end if 
#endif
      
         

      !ALLOCATION 
      
      if(this%how .eq. HS_FROM_DJDT) then
        SAFE_ALLOCATE(this%J(1:this%mesh%np, 1:3, 2))
        this%J = M_z0
      end if
      
      if (.not. use_nfft) then
        call cube_init(this%cube, box, this%mesh%sb, fft_type = FFT_COMPLEX, verbose = .true.)
        
      else
        
        
        ! we just add 2 points for the enlarged region
        box(1:dim) = box(1:dim) + 2 

  #ifdef HAVE_NFFT    
        !Set NFFT defaults to values that are optimal for PES (at least for the cases I have tested)
        !These values are overridden by the NFFT options in the input file 
        this%fft%nfft%set_defaults = .true.
        this%fft%nfft%guru = .true.
        this%fft%nfft%mm = 2 
        this%fft%nfft%sigma = CNST(1.1)
        this%fft%nfft%precompute = NFFT_PRE_PSI
  #endif
      
        ! These options should not affect NFFT scheme  
        optimize(1:3) = .false.
        optimize(this%mesh%sb%periodic_dim+1:this%mesh%sb%dim) = .true.
        optimize_parity(1:this%mesh%sb%periodic_dim) = 0
        optimize_parity(this%mesh%sb%periodic_dim+1:this%mesh%sb%dim) = 1

        call fft_init(this%fft, box, dim, FFT_COMPLEX, FFTLIB_NFFT, optimize, optimize_parity )
          
        Sfactor =  (maxval(M_PI / (this%mesh%spacing(1:dim)))*P_C) / Wmax        
        
        SAFE_ALLOCATE(XX(1:maxval(box(:)),3))
    
        !Generate the NFFT-enlarged node grid
        if (Sfactor > 0) then
          do idim = 1, dim
            do ii=2, box(idim)-1 
              XX(ii,idim)= (ii - int(box(idim)/2) -1)*this%mesh%spacing(idim)
            end do
            XX(1,idim)= (-int(box(idim)/2)) * this%mesh%spacing(idim) * Sfactor 
            XX(box(idim),idim)= (int(box(idim)/2)) * this%mesh%spacing(idim) * Sfactor 
          end do
      
        else
          do idim = 1, dim
            do ii=1, box(dim) 
              XX(ii,dim)= (ii - int(box(dim)/2) -1) * this%mesh%spacing(1)
            end do
          end do
        end if

        do idim = dim + 1, 3
          XX(:,idim) = XX(:,1)
        end do
    
        !Set the node points and precompute the NFFT plan
        call fft_init_stage1(this%fft, XX)

        SAFE_DEALLOCATE_A(XX)

        call cube_init(this%cube, box, this%mesh%sb)  
        SAFE_ALLOCATE(this%cube%fft)
        this%cube%fft = this%fft
        call fft_get_dims(this%cube%fft, this%cube%rs_n_global, this%cube%fs_n_global, this%cube%rs_n, this%cube%fs_n, &
             this%cube%rs_istart, this%cube%fs_istart)

        write(message(1),'(a, f10.3)') '  NFFT Scaling factor: ', Sfactor
        call messages_info(1)
      end if 
      
      
      call cube_init_fs_coords(this%cube, this%mesh%spacing(1:3), this%mesh%sb%dim, Sfactor)
      if ( this%mesh%parallel_in_domains .and. this%cube%parallel_in_domains) then
        call mesh_cube_parallel_map_init(this%mesh_cube_map, this%mesh, this%cube)
      end if      
      
      
      

      do ii = 1, dim
        temp(ii) = maxval(abs(this%cube%k(1:this%cube%fs_n(ii), ii)))*P_C
      end do
      omegaMax = maxval(temp(1:dim))
      omegaMin = sqrt(sum(abs(this%cube%k(2, 1:dim) - this%cube%k(1, 1:dim))**2))*P_C
      write(message(1),'(3a)') '  (Max Omega, Delta Omega)  [', trim(units_abbrev(units_out%energy)), ']   =  ('
      write(message(1),'(a, f10.3, a)') trim(message(1)), units_from_atomic(units_out%energy, omegaMax),', '
      write(message(1),'(a, f10.6, a)') trim(message(1)), units_from_atomic(units_out%energy, omegaMin) ,') '
      call messages_info(1)
      
      
      
      
      do i = 1, 3
      	call cube_function_null(this%Jkint(i))    
      	call cube_function_alloc_FS(this%cube, this%Jkint(i), force_alloc = .true.)   
        this%Jkint(i)%Fs = M_z0     
        
      	call cube_function_null(this%cftmp(i))    
      	call cube_function_alloc_FS(this%cube, this%cftmp(i), force_alloc = .true.)   
        this%cftmp(i)%Fs = M_z0     
        
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
        call cube_function_free_fs(this%cube, this%cftmp(i))        
      end do

      call cube_end(this%cube)   
      
      call harmonic_spect_nullify(this)
      
      if(this%how .eq. HS_FROM_DJDT) then
        SAFE_DEALLOCATE_P(this%J)
      end if
      
    end if
    POP_SUB(harmonic_spect_end)
  end subroutine harmonic_spect_end


  ! ---------------------------------------------------------
  subroutine harmonic_spect_calc(this, st, iter)
    type(harmonic_spect_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st
    integer,             intent(in)    :: iter

    FLOAT, allocatable :: Js(:,:,:)
    CMPLX, allocatable :: J(:,:)
    integer :: iw, it, ip, ir, dir
    FLOAT   :: NN(1:3), omega, theta, phi
    CMPLX   :: nnj(1:3) 
    type(cube_function_t) :: cfJ
  


    PUSH_SUB(harmonic_spect_calc)
    
    ASSERT(this%calc)
    ASSERT(states_are_complex(st))

  	call cube_function_null(cfJ)    
  	call zcube_function_alloc_RS(this%cube, cfJ)   
  	call  cube_function_alloc_FS(this%cube, cfJ)   

    SAFE_ALLOCATE(J(1:this%mesh%np, 1:3))

    SAFE_ALLOCATE(Js(1:this%mesh%np, 1:MAX_DIM, 1:st%d%nspin))    
    Js = M_ZERO
    J = M_z0  
    call states_calc_quantities(this%gr%der, st, paramagnetic_current = Js )

    J(:,:) = Js(:,:,1) ! sum up all the spin channels
    if(this%how == HS_FROM_DJDT) then
      this%J(:,:,2) = this%J(:,:,1)      
      this%J(:,:,1) = J(:, :)
      J(:,:) = (this%J(:,:,1) - this%J(:,:,2)) / this%dt 
    end if

    SAFE_DEALLOCATE_A(Js)

 
    
    do dir= 1, 3

      call zmesh_to_cube(this%mesh, J(:,dir), this%cube, cfJ, local = .true.)
      
      call zcube_function_rs2fs(this%cube, cfJ)
      
      this%cftmp(dir)%FS = cfJ%FS
      
      call apply_phase(cfJ)
      
      this%Jkint(dir)%FS(:,:,:) = this%Jkint(dir)%FS(:,:,:) + this%dt * cfJ%FS(:,:,:)
!       this%Jkint(dir)%FS(:,:,:) = cfJ%zRS(:,:,:)

    end do

   


    call zcube_function_free_RS(this%cube, cfJ)
    call  cube_function_free_FS(this%cube, cfJ)
    SAFE_DEALLOCATE_A(J)
    
    POP_SUB(harmonic_spect_calc)
    
    contains 
    
    subroutine apply_phase(cfJ)
      type(cube_function_t), intent(inout) :: cfJ
      
      integer :: ix,iy,iz
      FLOAT :: KK(3), K, t

      t = (iter - 1) * this%dt
      
      do ix = 1, this%cube%fs_n(1)
        KK(1) = this%cube%k(ix,1)
        do iy = 1, this%cube%fs_n(2)
          KK(2) = this%cube%k(iy,2)
          do iz = 1, this%cube%fs_n(3)
            KK(3) = this%cube%k(iz,3)
            
            K = sqrt(sum(KK(1:3)**2))
            cfJ%FS(ix, iy, iz) = cfJ%FS(ix, iy, iz) * exp( M_zI * t * K * P_C)
            
          end do
        end do
      end do
      
      
    end subroutine apply_phase
    
  end subroutine harmonic_spect_calc


  !***********************************************************
  ! OUTPUT ROUTINES
  !***********************************************************

  ! ---------------------------------------------------------
  subroutine harmonic_spect_checkpoint(this, iter)
    type(harmonic_spect_t), intent(in) :: this
    integer,                intent(in) :: iter
    
    PUSH_SUB(harmonic_spect_checkpoint)
    
    if(this%calc .and. mpi_grp_is_root(mpi_world)) then
      call harmonic_spect_restart_write(this)
      call harmonic_spect_out(this, iter)
    end if
    
    POP_SUB(harmonic_spect_checkpoint)
  end subroutine harmonic_spect_checkpoint

  ! ---------------------------------------------------------
  ! g(k) = |n x (n x J(k))|^2 = |J(k)|^2 - |n.J(k)|^2 
  ! with 
  ! n = k/|k|
  ! ---------------------------------------------------------
  subroutine harmonic_spect_gk(this, gk)
    type(harmonic_spect_t), intent(in) :: this
    FLOAT,                  intent(out):: gk(:,:,:)

    integer :: ix,iy,iz, idim
    FLOAT :: KK(3), K, scale
    CMPLX :: N(3), P(3), J(3)
    
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

          J(1) = this%Jkint(1)%FS(ix, iy, iz)
          J(2) = this%Jkint(2)%FS(ix, iy, iz)
          J(3) = this%Jkint(3)%FS(ix, iy, iz)
          P = zcross_product(N, zcross_product(N, J))
          gk(ix,iy,iz) = sum(abs(P(1:3))**2)
           
!           gk(ix, iy, iz) = abs(this%Jkint(1)%FS(ix, iy, iz))**2 &
!                          + abs(this%Jkint(2)%FS(ix, iy, iz))**2 &
!                          + abs(this%Jkint(3)%FS(ix, iy, iz))**2
!             
!           gk(ix, iy, iz) = gk(ix, iy, iz) &
!                          - abs(sum(N(1:3) * this%Jkint(1:3)%FS(ix, iy, iz)))**2 !&
! !                          - abs(N(2) * this%Jkint(2)%FS(ix, iy, iz))**2 &
! !                          - abs(N(3) * this%Jkint(3)%FS(ix, iy, iz))**2 
            
          if(this%how .eq. HS_FROM_J) gk(ix, iy, iz) = gk(ix, iy, iz) * K**2 
            
        end do
      end do
    end do
        
   
    
    gk(:,:,:) = gk(:,:,:)/ (CNST(4.0) * M_PI**2 * P_C)
    
    if(this%how .eq. HS_FROM_J) gk(:,:,:) = gk(:,:,:)/P_C**2
    
    ! This is needed in order to normalize the Fourier integral 
    scale = M_ONE
    do idim=1, this%mesh%sb%dim
      scale = scale *( this%mesh%spacing(idim)/sqrt(M_TWO*M_PI))**2
    end do
    
    gk = gk * scale
    

    POP_SUB(harmonic_spect_gk)  
  end subroutine harmonic_spect_gk  

  ! ---------------------------------------------------------
  subroutine harmonic_spect_write_cf(cf, Lk, dim, file, cut, cutdim)
    FLOAT,            intent(in) :: cf(:,:,:)
    FLOAT,            intent(in) :: Lk(:,:)
    integer,          intent(in) :: dim
    character(len=*), intent(in) :: file
    integer,          intent(in) :: cut
    integer,          intent(in) :: cutdim
  
       
    character(len=10) :: dirlab(3), cutlab(3)   
    character(len=256) :: filename, path
    integer :: iunit, ll(3), ii, ix, iy, iz
    integer, allocatable :: idx(:,:)
    FLOAT, allocatable   :: Lk_(:,:)
    FLOAT                :: KK(3), temp, Retemp, Imtemp
    CMPLX                :: tmp
    
    PUSH_SUB(harmonic_spect_write_cf)

    ASSERT(cut <= 3 .and. cut >=1)
    ASSERT(cutdim <= 3 .and. cutdim >=1)
    ASSERT(size(cf, 1) == size(Lk,1))
    
    select case (cutdim)
      case (1)
        cutlab = (/'ky=0.kz=0','kx=0.kz=0','kx=0.ky=0'/)
        
      case (2)
        cutlab = (/'kx=0','ky=0','kz=0'/)
    
      case (3)
    end select

    
    filename = trim(file)//'.'//trim(cutlab(cut))

    iunit = io_open(filename, action='write')
  
    ll = 1
    do ii = 1, dim
      ll(ii) = size(cf, ii) 
    end do
    

    SAFE_ALLOCATE(idx(maxval(ll(:)), 3))
    SAFE_ALLOCATE(Lk_(maxval(ll(:)), 3))
  

    do ii = 1, 3
      Lk_(:,ii) = Lk(:,ii)
      call sort(Lk_(1:ll(ii), ii), idx(1:ll(ii), ii)) !We need to sort the k-vectors in order to dump in gnuplot format
    end do  
!     print *,"idx(1)", idx(:,1)
!     print *,"idx(2)", idx(1:3,2), "ll(dir)/2 + 1=", int(ll(dir)/2) + 1,"idx", idx(ll(dir)/2 + 1, 3)
!     print *,"idx(3)", idx(1:3,3), "ll(dir)/2 + 1=", int(ll(dir)/2) + 1,"idx", idx(ll(dir)/2 + 1, 3)

    select case (cutdim)
      case (1)
        do ix = 1, ll(1)
          KK(1) = Lk_(ix, 1)
      
            select case (cut)
            
              case (1)
!             print *,ix,"cf",idx(ix, 1), idx(ll(2)/2 + 1, 2), idx(ll(3)/2 + 1, 3)

                tmp = cf(idx(ix, 1), idx(ll(2)/2 + 1, 2), idx(ll(3)/2 + 1, 3))
    
              case (2)
                tmp = cf(idx(ll(1)/2 + 1, 1), idx(ix, 2), idx(ll(3)/2 + 1, 3))    
      
              case (3)
                tmp = cf(idx(ll(1)/2 + 1, 1), idx(ll(2)/2 + 1, 2), idx(ix, 3))
    
            end select
            
            write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &
                    units_from_atomic(sqrt(units_out%energy), KK(1)),&
                    real(tmp), aimag(tmp)
       
        end do
    

      case (2)
        do ix = 1, ll(1)
          KK(1) = Lk_(ix, 1)
          do iy = 1, ll(2)
            KK(2) = Lk_(iy, 2)
      

            select case (cut)
              case (1)
                temp = cf(idx(ll(1)/2 + 1, 1), idx(ix, 2), idx(iy, 3))    
    
              case (2)
                temp = cf(idx(ix, 1), idx(ll(2)/2 + 1, 2), idx(iy, 3))    
      
              case (3)
                temp = cf(idx(ix, 1), idx(iy, 2), idx(ll(3)/2 + 1, 3))
    
            end select
        
            write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &
                    units_from_atomic(sqrt(units_out%energy), KK(1)),&
                    units_from_atomic(sqrt(units_out%energy), KK(2)),&
                    temp
 
       
          end do
          write(iunit, *)  
    
        end do
        
    end select
  
    call io_close(iunit)

    

    
    POP_SUB(harmonic_spect_write_cf)
  end subroutine harmonic_spect_write_cf


  ! ---------------------------------------------------------
  subroutine harmonic_spect_write_gk(this, gk)
    type(harmonic_spect_t), intent(in) :: this
    FLOAT,                  intent(in) :: gk(:,:,:)
       
    character(len=256) :: filename, path
    integer :: ierr, npoints
    
    PUSH_SUB(harmonic_spect_write_gk)

    path ='td.general/'
    
    npoints = this%cube%fs_n(1)*this%cube%fs_n(2)*this%cube%fs_n(3)
    
    filename = trim(path)//'harmonic_spect_gk.obf'
    call io_binary_write(filename, npoints, gk(:,:,:), ierr)    

    if(ierr > 0) then
      message(1) = "Failed to write file "//trim(filename)
      call messages_fatal(1)
    end if
    
    
    POP_SUB(harmonic_spect_write_gk)
  end subroutine harmonic_spect_write_gk

  ! ---------------------------------------------------------
  subroutine harmonic_spect_read_gk(gk, npoints)
    FLOAT,   intent(out) :: gk(:,:,:)
    integer, intent(in) :: npoints
       
    character(len=256) :: filename, path
    integer :: ierr
    
    PUSH_SUB(harmonic_spect_read_gk)

    path ='td.general/'
    
    filename = trim(path)//'harmonic_spect_gk.obf'
    call io_binary_read(filename, npoints, gk(:,:,:), ierr)    

    if(ierr > 0) then
      message(1) = "Failed to write file "//trim(filename)
      call messages_fatal(1)
    end if
    
    
    POP_SUB(harmonic_spect_read_gk)
  end subroutine harmonic_spect_read_gk


  ! ---------------------------------------------------------
  subroutine harmonic_spect_out(this, iter)
    type(harmonic_spect_t),      intent(in) :: this
    integer,                     intent(in) :: iter

    
    FLOAT, allocatable :: gk(:,:,:)
    character(len=256) :: file
    
    
    PUSH_SUB(harmonic_spect_out)

    ASSERT(this%calc)
    
    SAFE_ALLOCATE(gk(1:this%cube%fs_n(1),1:this%cube%fs_n(2),1:this%cube%fs_n(3)))
    
    if(this%mesh%sb%dim == 1) then
      gk = abs(this%Jkint(1)%FS(:,:,:))**2
    else
      call harmonic_spect_gk(this, gk)
    end if
    
    
    file ='td.general/harmonic_spect_j-x'
    call harmonic_spect_write_cf(abs(this%Jkint(1)%FS), this%cube%k, this%mesh%sb%dim, file, cut = 1, cutdim = 1)

    file ='td.general/harmonic_spect_gk'
    call harmonic_spect_write_cf(gk, this%cube%k, this%mesh%sb%dim, file, cut = 3, cutdim = 2)

    
!     write(file, '(a,i7.7)') "td.", iter
!     file=trim(file)//'/hs_jk-'    
!     call harmonic_spect_write_cf(abs(this%cftmp(1)%FS), this%cube%k, this%mesh%sb%dim, dir = 1, cut = 1, cutdim = 1, file = file)
    
    call harmonic_spect_write_gk(this, gk)
    
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
      filename = trim(path)//'hs_j-'//trim(number)//'.obf'
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
      filename = trim(path)//'hs_j-'//trim(number)//'.obf'      
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
  subroutine harmonic_spect_total(this, file, gk, bounds, interpolate)
    type(harmonic_spect_t), intent(in) :: this
    character(len=*),       intent(in) :: file
    FLOAT,                  intent(in) :: gk(:,:,:)
    FLOAT, optional,        intent(in) :: bounds(3)
    logical, optional,      intent(in) :: interpolate

    integer :: ist, ik, ii, ix, iy, iz, iunit,idim
    FLOAT ::  KK(3),vec

    integer :: nn
    FLOAT  :: step, DE
    FLOAT, allocatable :: npoints(:), out(:)

    ! needed for interpolation in 2D and 3D 
    FLOAT, pointer :: cube_f(:)
    type(qshep_t) :: interp

    FLOAT :: Dtheta, Dphi, theta, phi,EE , Wmin, Wmax, Wdir(3)
    integer :: np, Ntheta, Nphi, ith, iph
    integer :: ll(3), dim


    PUSH_SUB(harmonic_spect_total)


    dim = this%mesh%sb%dim
    ll = 1
    do ii = 1, dim
      ll(ii) = size(gk,ii) 
    end do
    
    
    if(.not. present(bounds)) then 
      step = M_ZERO
      do ii = 1, dim
        step= step + abs(this%cube%k(2, ii) - this%cube%k(1, ii))**2 
        Wdir(ii) = maxval(abs(this%cube%k(1:this%cube%fs_n(ii), ii)))*P_C
      end do
      step = sqrt(step)*P_C
      Wmax = maxval(Wdir(1:dim))
    else
      Wmin = bounds(1)
      Wmax = bounds(2)
      step = bounds(3)
    end if

!      print *,"-- Wmax", Wmax, "DW", step, "dim", dim, "ll", ll

    nn  = int(Wmax/step)


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
              vec = sqrt(sum(KK(1:dim)**2))*P_C
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
    write(iunit, '(a1,a18,a18)') '#', str_center("w", 20), str_center("H(w)", 20)
    write(iunit, '(a1,a18,a18)') &
      '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 20), &
      str_center('[' //trim(units_abbrev(units_out%length))//'/' &
          //trim(units_abbrev(units_out%time**2))//']', 20)
    write(iunit, '(a)') '##################################################'
 

    do ii = 1, nn
!       if(present(npoints)) then

        if(npoints(ii) > M_ZERO) then
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
