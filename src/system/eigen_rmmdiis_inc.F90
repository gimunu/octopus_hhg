!! Copyright (C) 2009 X. Andrade
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
!! $Id: eigen_rmmdiis_inc.F90 6342 2010-03-04 01:25:36Z dstrubbe $

! ---------------------------------------------------------
!> See http://prola.aps.org/abstract/PRB/v54/i16/p11169_1
subroutine X(eigensolver_rmmdiis) (gr, st, hm, pre, tol, niter, converged, ik, diff, blocksize)
  type(grid_t),           intent(in)    :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(preconditioner_t), intent(in)    :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  integer,                intent(in)    :: ik
  FLOAT,                  intent(out)   :: diff(1:st%nst)
  integer,                intent(in)    :: blocksize

  R_TYPE, allocatable :: res(:, :, :, :), tmp(:, :)
  R_TYPE, allocatable :: psi(:, :, :, :)
  R_TYPE, allocatable :: mm(:, :, :, :), evec(:, :, :)
  R_TYPE, allocatable :: eigen(:)
  FLOAT,  allocatable :: eval(:, :)
  FLOAT, allocatable :: lambda(:), nrm(:)
  integer :: ist, minst, idim, ip, ii, iter, nops, maxst, jj, bsize, ib
  R_TYPE :: ca, cb, cc
  R_TYPE, allocatable :: fr(:, :)
  type(profile_t), save :: prof
  type(batch_t), allocatable :: psib(:), resb(:)
  type(batch_t) :: psibit, resbit, stpsib
  integer, allocatable :: done(:), last(:)
  logical, allocatable :: failed(:)

  PUSH_SUB(X(eigensolver_rmmdiis))

  SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim, 1:niter, 1:blocksize))
  SAFE_ALLOCATE(res(1:gr%mesh%np_part, 1:st%d%dim, 1:niter, 1:blocksize))
  SAFE_ALLOCATE(tmp(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(lambda(1:st%nst))
  SAFE_ALLOCATE(psib(1:niter))
  SAFE_ALLOCATE(resb(1:niter))
  SAFE_ALLOCATE(done(1:blocksize))
  SAFE_ALLOCATE(last(1:blocksize))
  SAFE_ALLOCATE(failed(1:blocksize))
  SAFE_ALLOCATE(fr(1:4, 1:blocksize))
  SAFE_ALLOCATE(nrm(1:blocksize))

  nops = 0

  call profiling_in(prof, "RMMDIIS")

  failed = .false.

  do ib = st%block_start, st%block_end
    minst = states_block_min(st, ib)
    maxst = states_block_max(st, ib)
    bsize = maxst - minst + 1

    call batch_init(psib(1), st%d%dim, minst, maxst, psi(:, :, 1, :))
    call batch_init(resb(1), st%d%dim, minst, maxst, res(:, :, 1, :))

    if(batch_is_packed(st%psib(ib, ik))) then 
      call batch_pack(psib(1), copy = .false.)
      call batch_pack(resb(1), copy = .false.)
    end if

    call batch_copy_data(gr%mesh%np, st%psib(ib, ik), psib(1))

    call X(hamiltonian_apply_batch)(hm, gr%der, psib(1), resb(1), ik)
    nops = nops + bsize

    SAFE_ALLOCATE(eigen(1:bsize))

    call X(mesh_batch_dotp_vector)(gr%mesh, psib(1), resb(1), eigen)

    st%eigenval(minst:maxst, ik) = eigen(1:bsize)

    SAFE_DEALLOCATE_A(eigen)

    call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), psib(1), resb(1))

    done = 0

    call mesh_batch_nrm2(gr%mesh, resb(1), nrm)

    if(batch_is_packed(st%psib(ib, ik))) then 
      call batch_unpack(psib(1))
      call batch_unpack(resb(1))
    end if

    do ii = 1, bsize
      if(nrm(ii) < tol) done(ii) = 1
    end do

    if(all(done(1:bsize) /= 0)) cycle

    ! initialize the remaining batch objects
    do iter = 2, niter
      call batch_init(psib(iter), st%d%dim, maxst - minst + 1)
      call batch_init(resb(iter), st%d%dim, maxst - minst + 1)
      
      do ist = minst, maxst
        ii = ist - minst + 1
        call batch_add_state(psib(iter), ist, psi(:, :, iter, ii))
        call batch_add_state(resb(iter), ist, res(:, :, iter, ii))
      end do
    end do
    
    ! get lambda 
    call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb(1), psib(2))
    call X(hamiltonian_apply_batch)(hm, gr%der, psib(2), resb(2), ik)
    nops = nops + bsize

    call X(mesh_batch_dotp_vector)(gr%mesh, psib(2), psib(2), fr(1, :), reduce = .false.)
    call X(mesh_batch_dotp_vector)(gr%mesh, psib(1), psib(2), fr(2, :), reduce = .false.)
    call X(mesh_batch_dotp_vector)(gr%mesh, psib(2), resb(2), fr(3, :), reduce = .false.)
    call X(mesh_batch_dotp_vector)(gr%mesh, psib(1), resb(2), fr(4, :), reduce = .false.)

    if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, fr)

    do ist = minst, maxst
      ii = ist - minst + 1

      ca = fr(1, ii)*fr(4, ii) - fr(3, ii)*fr(2, ii)
      cb = fr(3, ii) - st%eigenval(ist, ik)*fr(1, ii)
      cc = st%eigenval(ist, ik)*fr(2, ii) - fr(4, ii)

      lambda(ist) = CNST(2.0)*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))

      ! restrict the value of lambda to be between 0.1 and 1.0
      if(abs(lambda(ist)) > CNST(1.0)) lambda(ist) = lambda(ist)/abs(lambda(ist))
      if(abs(lambda(ist)) < CNST(0.1)) lambda(ist) = CNST(0.1)*lambda(ist)/abs(lambda(ist))
    end do

    do iter = 2, niter

      ! for iter == 2 the preconditioning was done already
      if(iter > 2) call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb(iter - 1), psib(iter))

      do ist = minst, maxst
        ii = ist - minst + 1

        ! predict by jacobi
        forall(idim = 1:st%d%dim, ip = 1:gr%mesh%np)
          psi(ip, idim, iter, ii) = lambda(ist)*psi(ip, idim, iter, ii) + psi(ip, idim, iter - 1, ii)
        end forall
      end do

      ! calculate the residual
      call X(hamiltonian_apply_batch)(hm, gr%der, psib(iter), resb(iter), ik)
      nops = nops + bsize

      SAFE_ALLOCATE(mm(1:iter, 1:iter, 1:2, 1:bsize))
      SAFE_ALLOCATE(evec(1:iter, 1:1, 1:blocksize))
      SAFE_ALLOCATE(eval(1:iter, 1:blocksize))

      call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), psib(iter), resb(iter))

      do ist = minst, maxst
        ii = ist - minst + 1
        
        ! perform the diis correction

        !  mm_1(i, j) = <res_i|res_j>
        call batch_init(resbit, st%d%dim, 1, iter, res(:, :, :, ii))
        call X(mesh_batch_dotp_matrix)(gr%mesh, resbit, resbit, mm(:, :, 1, ii), symm = .true., reduce = .false.)
        call batch_end(resbit)

        !  mm_2(i, j) = <psi_i|psi_j>
        call batch_init(psibit, st%d%dim, 1, iter, psi(:, :, :, ii))
        call X(mesh_batch_dotp_self)(gr%mesh, psibit, mm(:, :, 2, ii), reduce = .false.)
        call batch_end(psibit)

      end do

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, mm)

      do ist = minst, maxst
        ii = ist - minst + 1

        failed(ii) = .false.
        call lalg_lowest_geneigensolve(1, iter, mm(:, :, 1, ii), mm(:, :, 2, ii), eval(:, ii), evec(:, :, ii), bof = failed(ii))
        if(failed(ii)) then
          last(ii) = iter - 1
          cycle
        end if

        ! generate the new vector and the new residual (the residual
        ! might be recalculated instead but that seems to be a bit
        ! slower).
        do idim = 1, st%d%dim
          call lalg_scal(gr%mesh%np, evec(iter, 1, ii), psi(:, idim, iter, ii))
          call lalg_scal(gr%mesh%np, evec(iter, 1, ii), res(:, idim, iter, ii))
        end do

        do jj = 1, iter - 1
          do idim = 1, st%d%dim
            call lalg_axpy(gr%mesh%np, evec(jj, 1, ii), psi(:, idim, jj, ii), psi(:, idim, iter, ii))
            call lalg_axpy(gr%mesh%np, evec(jj, 1, ii), res(:, idim, jj, ii), res(:, idim, iter, ii))
          end do
        end do

      end do

      SAFE_DEALLOCATE_A(mm)
      SAFE_DEALLOCATE_A(eval)      
      SAFE_DEALLOCATE_A(evec)

    end do

    ! end with a trial move
    call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb(niter), resb(niter - 1))

    do ist = minst, maxst
      ii = ist - minst + 1

      if(.not. failed(ii)) then

        forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np)
          res(ip, idim, niter - 1, ii) = psi(ip, idim, niter, ii) + lambda(ist)*res(ip, idim, niter - 1, ii)
        end forall
        
        call states_set_state(st, gr%mesh, ist, ik, res(:, :, niter - 1, ii))
      else
        call states_set_state(st, gr%mesh, ist, ik, psi(:, :, last(ii), ii))
      end if

      if(mpi_grp_is_root(mpi_world)) then
        call loct_progress_bar(st%nst * (ik - 1) +  ist, st%nst*st%d%nik)
      end if

    end do

    do iter = 1, niter
      call batch_end(psib(iter))
      call batch_end(resb(iter))
    end do

  end do

  call profiling_out(prof)

  call X(states_orthogonalization_full)(st, gr%mesh, ik)

  ! recalculate the eigenvalues and residuals

  SAFE_ALLOCATE(eigen(st%st_start:st%st_end))

  do ib = st%block_start, st%block_end
    minst = states_block_min(st, ib)
    maxst = states_block_max(st, ib)

    call batch_copy(st%psib(ib, ik), resb(1), reference = .false.)

    call X(hamiltonian_apply_batch)(hm, gr%der, st%psib(ib, ik), resb(1), ik)
    call X(mesh_batch_dotp_vector)(gr%der%mesh, st%psib(ib, ik), resb(1), eigen(minst:maxst))

    st%eigenval(minst:maxst, ik) = eigen(minst:maxst)

    call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), st%psib(ib, ik), resb(1))

    call X(mesh_batch_dotp_vector)(gr%der%mesh, resb(1), resb(1), eigen(minst:maxst))

    diff(minst:maxst) = sqrt(abs(eigen(minst:maxst)))

    call batch_end(resb(1), copy = .false.)

    nops = nops + maxst - minst + 1
  end do

  SAFE_DEALLOCATE_A(eigen)

  niter = nops

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(psib)
  SAFE_DEALLOCATE_A(resb)
  SAFE_DEALLOCATE_A(done)
  SAFE_DEALLOCATE_A(last)
  SAFE_DEALLOCATE_A(failed)
  SAFE_DEALLOCATE_A(fr)
  SAFE_DEALLOCATE_A(nrm)

  POP_SUB(X(eigensolver_rmmdiis))

end subroutine X(eigensolver_rmmdiis)

! ---------------------------------------------------------
subroutine X(eigensolver_rmmdiis_min) (gr, st, hm, pre, tol, niter, converged, ik, blocksize)
  type(grid_t),           intent(in)    :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(preconditioner_t), intent(in)    :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  integer,                intent(in)    :: ik
  integer,                intent(in)    :: blocksize

  integer, parameter :: sweeps = 5
  integer, parameter :: sd_steps = 2

  integer :: isd, ist, sst, est, ib, ii
  R_TYPE  :: ca, cb, cc
  R_TYPE, allocatable :: res(:, :, :), lambda(:)
  R_TYPE, allocatable :: kres(:, :, :)
  R_TYPE, allocatable :: me1(:, :), me2(:, :)

  type(batch_t) :: resb, kresb

  PUSH_SUB(X(eigensolver_rmmdiis_min))

  SAFE_ALLOCATE(me1(1:2, 1:st%d%block_size))
  SAFE_ALLOCATE(me2(1:4, 1:st%d%block_size))
  SAFE_ALLOCATE(res(1:gr%mesh%np_part, 1:st%d%dim, 1:st%d%block_size))
  SAFE_ALLOCATE(kres(1:gr%mesh%np_part, 1:st%d%dim, 1:st%d%block_size))
  SAFE_ALLOCATE(lambda(1:st%nst))

  niter = 0

  do ib = st%block_start, st%block_end
    sst = states_block_min(st, ib)
    est = states_block_max(st, ib)

    call batch_copy(st%psib(ib, ik), resb, reference= .false.)
    call batch_copy(st%psib(ib, ik), kresb, reference= .false.)

    do isd = 1, sd_steps

      call X(hamiltonian_apply_batch)(hm, gr%der, st%psib(ib, ik), resb, ik)

      call X(mesh_batch_dotp_vector)(gr%mesh, st%psib(ib, ik), resb, me1(1, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, st%psib(ib, ik), st%psib(ib, ik), me1(2, :), reduce = .false.)

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, me1, (/2, st%d%block_size/))

      forall(ist = sst:est) st%eigenval(ist, ik) = me1(1, ist - sst + 1)/me1(2, ist - sst + 1)
 
      call batch_axpy(gr%mesh%np, -st%eigenval(:, ik), st%psib(ib, ik), resb)

      call X(preconditioner_apply_batch)(pre, gr, hm, ik, resb, kresb)

      call X(hamiltonian_apply_batch)(hm, gr%der, kresb, resb, ik)

      niter = niter + 2*(est - sst + 1)

      call X(mesh_batch_dotp_vector)(gr%mesh, kresb, kresb, me2(1, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, st%psib(ib, ik),  kresb, me2(2, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, kresb, resb,  me2(3, :), reduce = .false.)
      call X(mesh_batch_dotp_vector)(gr%mesh, st%psib(ib, ik),  resb,  me2(4, :), reduce = .false.)

      if(gr%mesh%parallel_in_domains) call comm_allreduce(gr%mesh%mpi_grp%comm, me2, (/4, st%d%block_size/))

      do ist = sst, est
        ii = ist - sst + 1

        ca = me2(1, ii) * me2(4, ii) - me2(3, ii) * me2(2, ii)
        cb = me1(2, ii) * me2(3, ii) - me1(1, ii) * me2(1, ii)
        cc = me1(1, ii) * me2(2, ii) - me1(2, ii) * me2(4, ii)

        lambda(ist) = CNST(2.0)*cc/(cb + sqrt(cb**2 - CNST(4.0)*ca*cc))

      end do

      call batch_axpy(gr%mesh%np, lambda, kresb, st%psib(ib, ik))

    end do

    call batch_end(resb, copy = .false.)
    call batch_end(kresb, copy = .false.)

    if(mpi_grp_is_root(mpi_world)) then
      call loct_progress_bar(st%nst*(ik - 1) +  est, st%nst*st%d%nik)
    end if

  end do

  call X(states_orthogonalization_full)(st, gr%mesh, ik)

  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(me1)
  SAFE_DEALLOCATE_A(me2)
  SAFE_DEALLOCATE_A(res)
  SAFE_DEALLOCATE_A(kres)

  POP_SUB(X(eigensolver_rmmdiis_min))

end subroutine X(eigensolver_rmmdiis_min)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
