! Joint collaboration between Kuan-Hsiang Hsu and Chunjing Jia
! Revised to calculate Matrix with imaginary elements
!	Nov 1, 2024

module NumOfOrbitalAndElectrons
   integer*8 :: N, N3d, N2p, Ndim
   integer*8 :: nup, ndn, Hsize
   integer*8 :: Nsite, Nkpt
endmodule NumOfOrbitalAndElectrons

module ModelParas
   double precision :: Utratio, ttprimeratio, tttprimeratio
   double precision :: U, E_site(1:30)
   double precision :: Uc, t, tt, ttt, miu
   double precision :: Ud, Up, UQ, E_d, E_pxy, E_pz, tpd, tpp, tpdz, tpzdz, tppz
   double precision :: U_eg, U_t2g, Upp, E_eg, E_t2g, cfs
   double precision :: rA, rB, rC, F_0, F_2, G_1, G_3, Dq, Ds, Dt
   double precision :: U_rest(5,5,5,5), U_pddp(3,5,5,3), U_dpdp(5,3,5,3)
   double precision :: U_ext(5,5), J_ext(5,5)
   double precision :: phase(30,30)
   double complex :: U_SO_2p(3,2,3,2)
   double precision :: xi_2p
   double complex :: disp(0:19,3), den(0:19), delta(5,2), lambda ! Delta are site locations
endmodule ModelParas

module MPIParas
   integer  ::   comm, myid, nprocs
   integer  ::   source, dest, tagmpi, ierr, tag
   integer*8 ::  nloc, localstart, localend
   integer*8 ::  nloc_array(1024)
endmodule MPIParas

module ConstantParas
   integer*8 :: HsizeEst=25000000, nzeprEst= 400
   integer*8 :: NNN=1000, niter_CFE = 150, niter_CG = 150
   double precision :: tol_CG=0.0001, tol_ED = 0.000001
   double precision :: pi=3.1415926535, sqrthalf = 0.707106781
   double precision :: epsilone_CFE=0.3d0,  epsilone_CG = 1.0d0
endmodule ConstantParas

module BettsCluster
   integer*8 :: Nmore
   integer*8 :: rpoint(0:19,2), runit
   integer*8 :: kpoint(0:19,2), kunit
   integer*8 :: site(0:19,4)
   integer*8 :: nQpoints, Kmap(1:19)
   double complex :: QPhase(0:19,0:8)
   double precision :: Qpoint(0:19,2)
   integer*8 :: Kmapsize
   integer :: twistx, twisty, twist1, twist2
endmodule BettsCluster

module PolarParas
   double precision :: k_theta_in(0:15), k_theta_out(0:15)
   double precision :: k_phi_in(0:15), k_phi_out(0:15)
   double precision :: p_theta_in(0:15), p_phi_in(0:15)
   double precision :: p_theta_out(0:15), p_phi_out(0:15)
   double precision :: Gn_in(3,0:15), Gn_out(3,0:15)
   character*1 :: PorS_in, PorS_out
endmodule PolarParas

module ScanRegion
   integer*8 :: divX, divY
   integer*8 :: eneX, eneY, eneY_start, eneY_read
   double precision :: startX, endX, startY, endY
endmodule ScanRegion

module Measurements
   integer :: charge_gap, opt_cond, occ
   integer :: static_corr, dynamic_corr
   integer :: nres_raman, res_raman, eff_raman
   integer :: a1g, a2g, b1g, b2g
endmodule Measurements


Program Main
use NumOfOrbitalAndElectrons; use ModelParas; use BettsCluster
use ConstantParas; use PolarParas; use ScanRegion; use MPIParas
use mpi; use Measurements
implicit none
! include 'mpif.h'
integer*8 :: SprSize, SprSize_f, ksize_f
integer*8 :: ii, jj, kk, mm, nn, i, j, sz, site1, site2
integer*8 :: Ediv, Eloss_div
integer*8 :: thread_num, ctag
integer*8 :: qx, qy, gsmom, gssz, szmax
integer*8 :: nupinp, ndninp, Hsizeinp
integer*8, external :: factorial
integer*8, external :: omp_get_num_procs
integer*8, external :: ksubtract, ksum
integer*8, allocatable :: IndexI(:), IndexJ(:)
integer*8, allocatable :: IndexI_f(:), IndexJ_f(:)
integer*8, allocatable :: Hsp(:), Hsp_f(:)
double precision :: maxE, minE, maxE_loss, minE_loss
double precision :: gsen
double precision :: time1, time2, rtemp_1, tsum, E_arr(50)
double precision, external :: Lorentzian
double precision :: E, E_0, Em1_0, Ep1_0, temp_read
double precision :: specX(6001), specY(6001), specY_read(6001), read_ev(2000)
double complex :: gamma, z, trace, sisj, ninj, docc
double complex, external :: lsf, lsf_chiral, lsf_test
double complex, allocatable :: sparseH(:), sparseH_f(:) ! Complex value in Hamiltonian
double complex, allocatable :: specraman(:), specopt(:), ni(:)
double complex, allocatable :: tempv(:), tempv1(:), tempvcmplx(:), tempv1cmplx(:)
double complex, allocatable :: H_0(:), H_f(:) ! Complex eigenvectors
integer*8 :: qq, qqk, nqq, ferr, order, symm, skip
integer :: diag, nevin, evsize, ee, iostat
integer*8, external :: nchoosek
character(1000) :: filename, symmname, iomsg
double precision, allocatable :: H_0_copy(:)
logical, allocatable :: H_0_max(:)

      call MPI_INIT( ierr )
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
      if(myid.eq.0) write(*,*) '# of processors = ', nprocs

open(unit=10101,file='input',Status='old');
read(10101,'(3I8)') N3d, Ndim, Nmore
read(10101,'(2I8)') nup, ndn
read(10101,'(3F8.2)') Utratio, ttprimeratio, tttprimeratio
read(10101,'(2F8.2)') maxE_loss,minE_loss
read(10101,'(2F8.2)') maxE,minE
read(10101,'(2I8)') Ediv,Eloss_div
read(10101,'(2I8)') niter_CFE, niter_CG
read(10101,'(1F8.6)') tol_CG
read(10101,'(2F8.2)') epsilone_CG, epsilone_CFE
read(10101,'(2I8)') diag, nevin
read(10101,'(5I4)') occ, charge_gap, opt_cond, static_corr, dynamic_corr
read (10101, '(3I4)') nres_raman, res_raman, eff_raman
read (10101, '(4I4)') a1g, a2g, b1g, b2g
read (10101, '(2I4)') twistx, twisty
close(10101)
tol_CG = 0.0001;

if(myid.eq.0) then
write(*,*) ''
write(*,*) '     N3d     Ndim    Nmore'
write(*,'(3I8)') N3d, Ndim, Nmore
write(*,*) ''
write(*,*) '     nup     ndn'
write(*,'(2I8)') nup, ndn
write(*,*) ''
if (Ndim.eq.5) then
   write(*,*) '     U/t     delta/t_p, ???'
else
   write(*,*) '     U/t     t/t_p, t/t_pp'
endif
write(*,'(2F8.4)') Utratio, ttprimeratio, tttprimeratio
write(*,*) ''
write(*,*) 'maxE_loss minE_loss'
write(*,'(2F8.2)') maxE_loss,minE_loss
write(*,*) ''
write(*,*) '    maxE    minE'
write(*,'(2F8.2)') maxE,minE
write(*,*) ''
write(*,*) '    Ediv Eloss_div'
write(*,'(2I8)') Ediv,Eloss_div
write(*,*) ''
write(*,*) '  CFE      CG'
write(*,'(2I8)') niter_CFE, niter_CG
write(*,*) ''
write(*,*) '  tol_CG   tol_ED'
write(*,'(2F8.2)') tol_CG, tol_ED
write(*,*) ''
write(*,*) 'epsilone_CG, CFE'
write(*,'(2F8.2)') epsilone_CG, epsilone_CFE
write(*,*) ''
write(*,*) 'diag, nevin'
write(*,'(2I4)') diag, nevin
write(*,*) ''
write(*,*) 'occ, charge_gap, opt_cond, static_corr, dynamic_corr'
write(*,'(4I4)') occ, charge_gap, opt_cond, static_corr, dynamic_corr
write(*,*) ''
write(*,*) 'nres_raman, res_raman, eff_raman'
write(*,'(3I8)') nres_raman, res_raman, eff_raman
write(*,*) ''
write(*,*) 'a1g, a2g, b1g, b2g'
write(*,'(4I4)') a1g, a2g, b1g, b2g
write(*,*) ''
write(*,*) 'twistx, twisty'
write(*,'(2I4)') twistx, twisty
write(*,*) ''
write(*,*) '******************************'
write(*,*) ''

endif

skip = 0
if (epsilone_CFE.eq.0) then
   epsilone_CFE = 0.1*4/Utratio
   if (myid.eq.0) write(*,*) "NEW epsilone_CFE: ", epsilone_CFE
endif
if (Ediv.eq.1.and.minE.ne.maxE) then
   maxE = minE
   if (myid.eq.0) write(*,*) "Doing only 1 incident energy: ", maxE
endif
if (twistx .lt. 1) then
   if (myid.eq.0) write(*,*) "WARNING: twistx < 1, setting to 1"
   twistx = 1
endif
if (twisty .lt. 1) then
   if (myid.eq.0) write(*,*) "WARNING: twisty < 1, setting to 1"
   twisty = 1
endif
! Negative maxE loss means this is in J unit
!============== The Step 0: preparation ==============
call SetModelParameters
if(myid.eq.0) write(*,*) 'SetModelPara and Setnloc are done'
if (maxE_loss.le.0) then
   maxE_loss = -maxE_loss * 4 / U
   if (myid.eq.0) write(*,*) "new maxE_loss in t: ", maxE_loss
endif
if (maxE.le.0) then
   maxE = -maxE * 4 / U
   if (myid.eq.0) write(*,*) "new maxE in t: ", maxE
endif
N2p=0
N=N3d+N2p
startX = minE_loss; endX = maxE_loss; divX = Eloss_div-1
startY = minE; endY = maxE; divY = Ediv-1
eneY_start = 0
if (minE .eq. maxE) then 
   divY = 1
   eneY_start = 1
endif
if (N3d.ge.18) then 
   nzeprEst= 200
   HsizeEst = nchoosek(N3d,nup) * nchoosek(N3d,ndn) / N3d * 1.1 * Nsite
   if(myid.eq.0) then 
      write(*,*) 'Large Matrix, re-estimating HsizeEst: ',HsizeEst
      write(*,*) 'HsizeEst/nprocs*nzeprEst = ', HsizeEst/nprocs*nzeprEst/Nsite
   end if 
else
   if(myid.eq.0) then 
      write(*,*) 'HsizeEst, estimate: ',HsizeEst, nchoosek(N3d,nup) * nchoosek(N3d,ndn) / N3d * 1.1 * Nsite
      write(*,*) 'HsizeEst/nprocs*nzeprEst = ', HsizeEst/nprocs*nzeprEst
   end if 
endif
allocate(Hsp(HsizeEst))
allocate(H_0(HsizeEst))
allocate(tempv(HsizeEst))
allocate(tempv1(HsizeEst))
allocate(H_0_copy(HsizeEst))
allocate(H_0_max(HsizeEst))
! allocate(Hsp_f(HsizeEst))
! allocate(IndexI_f(HsizeEst/nprocs*nzeprEst*2))
! allocate(IndexJ_f(HsizeEst/nprocs*nzeprEst*2))
! allocate(sparseH_f(HsizeEst/nprocs*nzeprEst))
if (N.lt.18) then
   allocate(IndexI(HsizeEst/nprocs*nzeprEst/Nsite))
   allocate(IndexJ(HsizeEst/nprocs*nzeprEst/Nsite))
   allocate(sparseH(HsizeEst/nprocs*nzeprEst/Nsite))
   if (charge_gap.ne.0 .or. static_corr.ne.0 .or. occ.ne.0 .or. dynamic_corr.ne.0) then
      if(myid.eq.0) write(*,*) 'Charge Gap/Static correlation/Occupation, allocating non 0 q point matrices'
      allocate(Hsp_f(HsizeEst))
      allocate(IndexI_f(HsizeEst/nprocs*nzeprEst))
      allocate(IndexJ_f(HsizeEst/nprocs*nzeprEst))
      allocate(sparseH_f(HsizeEst/nprocs*nzeprEst))
      Hsp_f=0
      SprSize_f = HsizeEst/nprocs*nzeprEst
   endif
else
   if (static_corr.ne.0.or.dynamic_corr.ne.0) then
      if(myid.eq.0) write(*,*) 'Large matrix, doing only non-q point calculation'
      allocate(Hsp_f(HsizeEst))
      Hsp_f = 0
      ! Hamiltonian not needed for static correlation
      ! allocate(IndexI_f(HsizeEst/nprocs*nzeprEst/Nsite))
      ! allocate(IndexJ_f(HsizeEst/nprocs*nzeprEst/Nsite))
      ! allocate(sparseH_f(HsizeEst/nprocs*nzeprEst/Nsite))
      nres_raman = 0
      res_raman = 0 
      eff_raman = 0
      diag = 0 ! Diag set to 0 to save memory
      skip = 1
   else
      if(myid.eq.0) write(*,*) 'Large matrix, doing only q point calculation'
      allocate(IndexI(HsizeEst/nprocs*nzeprEst/Nsite))
      allocate(IndexJ(HsizeEst/nprocs*nzeprEst/Nsite))
      allocate(sparseH(HsizeEst/nprocs*nzeprEst/Nsite))
   endif 
endif
SprSize = HsizeEst/nprocs*nzeprEst
if(myid.eq.1) then 
   write(*,*) 'Allocations are done'
   write(*,*) 'Array size (Gb): ', Real(HsizeEst)/nprocs*nzeprEst *8/1e9
   write(*,*) 'Upperbound of IndexI', Ubound(IndexI,1,16)
   write(*,*) 'Size of IndexI', Size(IndexI,1,16)
endif
E_0=0.0d0;
specX=0.0d0;
specY=0.0d0;
specY_read=0.0d0;
!!!
if (skip.ne.1) then
   IndexI = 0
   if(myid.eq.0) write(*,*) "Allocating indexI"
   IndexJ = 0
   if(myid.eq.0) write(*,*) "Allocating indexJ"
   sparseH = 0
   if(myid.eq.0) write(*,*) "Allocating sparseH"
endif
Hsp = 0
if(myid.eq.0) write(*,*) "Allocating Hsp"
H_0 = 0
H_0_copy = 0
H_0_max = .True.
if(myid.eq.0) write(*,*) "Allocating H0"
tempv = 0
if(myid.eq.0) write(*,*) "Allocating tempv"
tempv1 = 0
if(myid.eq.0) write(*,*) "Allocating tempv1"

!!!
! stop
! if (myid.eq.0) then
!    gamma = lsf_test(int(1,8),30)
! endif 

! write(*,*) "tt", tt
! if (myid.eq.0) then
!    do kk = 0,N-1
!       do qq = 0,N-1
!          qqk = ksubtract(int(0,8),kk)
!          qqk = ksubtract(qqk,qq)
!          write(*,*) kpoint(kk,1),kpoint(kk,2),kpoint(qq,1),kpoint(qq,2),kpoint(qqk,1),kpoint(qqk,2)
!          ! gamma = lsf(qq,qqk,int(1,8))
!          gamma = lsf_chiral(kk,qq,qqk,1)
!       enddo
!    enddo
! endif
! stop
! Find lowest eigenvalue in each spin sector
gsen = 10000
if (diag.eq.2) then
   if (myid.eq.0) then 
      call system("rm eigval_sectors.dat")
      open(111,file='eigval_sectors.dat',access='append')
      write(111,*) "    Sz,    qx,    qy,     E_0"
      close(111)
   endif
   nupinp = nup
   ndninp = ndn
   do twist1 = 0, twistx-1
      do twist2 = 0, twisty-1
         if (myid.eq.0) write(*,*) 'Finding ground state energy in different sectors'
         ! This is for mom = 0, sz = 0
         if (myid.eq.0) then
            write(*,*) "Doing twist, x, y:", twist1, twist2
            open(111,file='eigval_sectors.dat',access='append')
            write(111,*) "twist, x, y:", twist1, twist2
            close(111)
         endif
         ! Loop through different sector
         if (N3d.ge.18) then
            szmax = 0
         else 
            szmax = 0
         endif
         do sz = 0,szmax
            nup = nupinp + sz
            ndn = ndninp - sz
            if (myid.eq.0) write(*,*) "Doing sz sector:", sz
            call MPI_BCAST(nup, 1, MPI_INTEGER8, 0, comm, ierr)
            call MPI_BCAST(ndn, 1, MPI_INTEGER8, 0, comm, ierr)
            call SetModelParameters
            do qqk = 1, Kmapsize
               ! Loop over interesting Kpoints
               qq = Kmap(qqk)
               qx = kpoint(qq,1)
               qy = kpoint(qq,2)
               if (myid.eq.0) write(*,*) "Doing mom:", qx, qy
               Hsp = 0.0
               IndexI = 0.0
               IndexJ = 0.0
               sparseH = 0.0
               SprSize = HsizeEst/nprocs*nzeprEst
               ! if (myid.eq.0) write(*,*) 'Upperbound of array', Ubound(IndexJ)
               ! if (myid.eq.0) write(*,*) 'Size of array', Size(IndexJ)
               call GenHsp_kspace(qq, Hsp, Hsize)
               call MPI_BCAST(Hsize, 1, MPI_INTEGER8, 0, comm, ierr)
               call MPI_BCAST(Hsp, Hsize, MPI_INTEGER8, 0, comm, ierr)
               call Setnloc(Hsize)
               if(myid.eq.0) write(*,*) 'Setnloc is done'
               call GenMatrix_kspace(Hsp, Hsize, IndexI, IndexJ, sparseH, SprSize)
               if(myid.eq.0) write(*,*) 'GenMatrix is done and SprSize = ', SprSize
               call ED_CMPLX_PARPACK(Hsize, SprSize, IndexI, IndexJ, sparseH, E_0, H_0, nevin)
               if(myid.eq.0) write(*,*) 'ED is done and E_0 = ', E_0
               if (E_0.le.gsen) then
                  if (myid.eq.0) call system("mv eigenvector_1.h5 eigenvector_1-lowest.h5")
                  gsen = E_0
                  gsmom = qq
                  gssz = sz
               endif
               if (myid.eq.0) then
                  open(unit=112,file='eigenvalue.dat',Status='old');
                  do ii = 1,nevin
                     read(112,'(F20.5)') E_arr(ii)
                  enddo
                  close(112)
                  open(111,file='eigval_sectors.dat',access='append')
                  write(111,*) sz, qx, qy, E_arr(1), E_arr(2), E_arr(3), E_arr(4)
                  close(111)
               endif 
            enddo
         enddo
      enddo
   enddo
   if (myid.eq.0) then
      call system("rm targetev.dat targetsector.dat")
      write(*,*) "gs energy, mom, sz: ", gsen, gsmom, gssz
      open(unit=113,file='targetev.dat',access='append');
         write(113,*) gsen
      close(113)
      open(unit=114,file='targetsector.dat',access='append');
         write(114,'(2I4)') gsmom, gssz
      close(114)
      call system("mv eigenvector_1-lowest.h5 eigenvector_1.h5")
   endif
   diag = 0 ! Eigenvector and eigenvalues does not need to be recalculated
   stop
endif


do twist1 = 0, twistx-1
   do twist2 = 0, twisty-1
      call SetModelParameters
      ! if diag = 0/1, read in momentum and sz sector 
      if (diag.le.0.or.diag.eq.1) then
         if (myid.eq.0) then
            open(unit=10103, file='targetsector.dat', status='old', action='read',iostat=ferr)
            if (ferr == 0) then
               read(10103,'(2I4)') gsmom, gssz
               write(*,*) 'Doing momentum sector: ', kpoint(gsmom,1), kpoint(gsmom,2)
               write(*,*) 'Doing gssz sector: ', gssz
               close(10103)
            else
               write(*,*) 'targetsector.dat not found, doing mom = (0,0), sz = 0'
               gsmom = 0
               gssz = 0
            endif
         end if
      endif
      
      if (myid.eq.0) write(*,*) "Doing twist, x, y:", twist1, twist2
      !============== The Step 1: Diagonalize ==============
      ! Charge Gap, E(N+1), E(N-1)
      nupinp = nup
      ndninp = ndn
      if (charge_gap.ne.0) then
         if(myid.eq.0) then
            write(*,*) '=============================='
            write(*,*) '      Measure Charge Gap      '
            write(*,*) '=============================='
            ! Calculate N-1
            Hsize = 0
            nup = nup
            ndn = ndn - 1
            call GenHsp_kspace(gsmom, Hsp_f, Hsize)
            write(*,*) 'GenHsp is done, N-1 Hsize = ', Hsize
         endif
         call MPI_BCAST(nup, 1, MPI_INTEGER8, 0, comm, ierr)
         call MPI_BCAST(ndn, 1, MPI_INTEGER8, 0, comm, ierr)
         call MPI_BCAST(Hsize, 1, MPI_INTEGER8, 0, comm, ierr)
         call MPI_BCAST(Hsp_f, Hsize, MPI_INTEGER8, 0, comm, ierr)

         call Setnloc(Hsize)
         if(myid.eq.0) write(*,*) 'Setnloc is done'
         call GenMatrix_kspace(Hsp_f, Hsize, IndexI_f, IndexJ_f, sparseH_f, SprSize_f)
         if(myid.eq.0) write(*,*) 'GenMatrix is done and SprSize = ', SprSize_f
         call ED_CMPLX_PARPACK(Hsize, SprSize_f, IndexI_f, IndexJ_f, sparseH_f, Em1_0, tempv, 2)
         if(myid.eq.0) write(*,*) 'ED is done and E_0 = ', Em1_0
         
         if(myid.eq.0) then
            ! Calculate N+1
            Hsize = 0
            nup = nupinp
            ndn = ndninp + 1
            call GenHsp_kspace(gsmom, Hsp_f, Hsize)
            write(*,*) 'GenHsp is done, N+1 Hsize = ', Hsize
         endif
         call MPI_BCAST(nup, 1, MPI_INTEGER8, 0, comm, ierr)
         call MPI_BCAST(ndn, 1, MPI_INTEGER8, 0, comm, ierr)
         call MPI_BCAST(Hsize, 1, MPI_INTEGER8, 0, comm, ierr)
         call MPI_BCAST(Hsp_f, Hsize, MPI_INTEGER8, 0, comm, ierr)
         SprSize_f = HsizeEst/nprocs*nzeprEst
         ! I need to check if odd number is correct
         call Setnloc(Hsize)
         if(myid.eq.0) write(*,*) 'Setnloc is done'
         call GenMatrix_kspace(Hsp_f, Hsize, IndexI_f, IndexJ_f, sparseH_f, SprSize_f)
         if(myid.eq.0) write(*,*) 'GenMatrix is done and SprSize = ', SprSize_f
         call ED_CMPLX_PARPACK(Hsize, SprSize_f, IndexI_f, IndexJ_f, sparseH_f, Ep1_0, tempv, 2)
         if(myid.eq.0) write(*,*) 'ED is done and N+1 E_0 = ', Ep1_0
      endif

      ! Calculate Ground state
      ! 0 stands for (0,0) momentum
      Hsize = 0
      nup = nupinp + gssz
      ndn = ndninp - gssz
      if (diag.eq.-1) gsmom = -1
      if(myid.eq.0) then
         write(*,*) 'gsmom, nup, ndn', gsmom, nup, ndn
         call GenHsp_kspace(gsmom, Hsp, Hsize)
         write(*,*) 'GenHsp is done, Hsize = ', Hsize
      endif
      if (diag.eq.-1) gsmom = 0 !This is a bit irresponsible manuver.
      call MPI_BCAST(Hsize, 1, MPI_INTEGER8, 0, comm, ierr)
      call MPI_BCAST(Hsp, Hsize, MPI_INTEGER8, 0, comm, ierr)
      call Setnloc(Hsize)
      if(myid.eq.0) write(*,*) 'Setnloc is done'
      if (skip.ne.1) then
         call GenMatrix_kspace(Hsp, Hsize, IndexI, IndexJ, sparseH, SprSize)
         if(myid.eq.0) write(*,*) 'GenMatrix is done and SprSize = ', SprSize
      else
         if(myid.eq.0) write(*,*) 'Skipping genmatrix'
      end if

      if (diag.eq.0) then
         if (myid.eq.0) then
            ! Eigenvector naming format: eigenvector_1xy, x=twistx, y=twisty
            if (twist1.eq.0 .and. twist2.eq.0) then
               call readhdf5cmplx(1,Hsize,H_0)
            else 
               write(*,*) "Cannot reread twist eigenvectors for now, stopping..."
               stop
            endif
            open(unit=10102, file='targetev.dat', status='old', action='read',iostat=ferr)
            if (ferr == 0) then
               read(10102,'(F20.16)') E_0
               write(*,*) 'ED is skipped and E_0 = ', E_0
               close(10102)
            else
               write(*,*) 'diag = 0 but no targetev.dat found, stopping...'
               stop
            endif
         end if
         call MPI_BCAST(E_0, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
         call MPI_BCAST(H_0, Hsize, MPI_DOUBLE_COMPLEX, 0, comm, ierr)
      else if (diag.eq.1) then
         if (myid.eq.0) write(*,*) 'Diagonalizing using PARPACK'
         call ED_CMPLX_PARPACK(Hsize, SprSize, IndexI, IndexJ, sparseH, E_0, H_0, nevin)
         if(myid.eq.0) write(*,*) 'ED is done and E_0 = ', E_0
         if(myid.eq.0) then 
            do ii = 1,Hsize
               H_0_copy(ii) = abs(H_0(ii))
            enddo
            do ii = 1,15
               nn = MAXLOC(H_0_copy,Hsize,H_0_max)
               write(*,*) "MAXloc: ", nn, Hsp(nn), H_0(nn)
               H_0_max(nn) = .False.
            enddo
         endif
      else if (diag.eq.-1) then
         if (myid.eq.0) write(*,*) 'Diagonalizing using PARPACK, with no spin symmetry'
         ! call sparse_trace(SprSize, IndexI, IndexJ, sparseH, trace)
         ! call sparse_symmetric(SprSize, IndexI, IndexJ, sparseH, Hsize)
         call ED_CMPLX_PARPACK(Hsize, SprSize, IndexI, IndexJ, sparseH, E_0, H_0, nevin)
         if (myid.eq.0) write(*,*) 'ED is done and E_0 = ', E_0
         if (myid.eq.0) write(*,*) 'trace = ', trace
         if (myid.eq.0) write(*,*) 'Stopping......'
         stop
      else 
         if (myid.eq.0) write(*,*) 'Diagonalize the matrix again at correct sector'
         call ED_CMPLX_PARPACK(Hsize, SprSize, IndexI, IndexJ, sparseH, E_0, H_0, nevin)
         if(myid.eq.0) write(*,*) 'ED is done and E_0 = ', E_0
         ! Loop through different sector
      endif

      !============== The Step 1.A: Measure some properties ==============
      ! Charge Gap
      if (charge_gap.ne.0) then
         open(114, file='charge_gap.dat')
         if(myid.eq.0) write(*,*) 'E_0, E_(N-1), E_(N+1): ', E_0, Em1_0, Ep1_0
         if(myid.eq.0) write(*,*) 'Charge Gap: ', (Em1_0+Ep1_0)-2*E_0
         if(myid.eq.0) write(114,*) 'Charge Gap: ', (Em1_0+Ep1_0)-2*E_0
         close(114)
      endif
      if (myid.eq.0) then
         if (diag.eq.0) then
            evsize = 0
            open(unit=10102, file='targetev.dat', status='old', action='read',iostat=ferr)
            ! Check if the file opened successfully
            if (ferr == 0) then
               do
                  read(10102,'(1F20.16)',iostat=ferr) E
                  if (ferr /= 0) exit ! Exit loop on end of file or error
                  evsize = evsize + 1
                  read_ev(evsize) = E
               end do
               if (evsize.eq.0) then
                  evsize = 1
                  read_ev(1) = E_0
               end if
               close(10102)
            else
               write(*,*) 'Error opening targetev.dat'
            end if
         else
            evsize = 1
            read_ev(1) = E_0
         end if
      end if
      call MPI_BCAST(evsize, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(read_ev, size(read_ev), MPI_DOUBLE_PRECISION, 0, comm, ierr)

      do ee = 1, evsize
         if (myid.eq.0) then
            call readhdf5cmplx(ee,Hsize,H_0)
            write(*,*) 'Calculating properties for energy: ', read_ev(ee)
         end if
         call MPI_BCAST(H_0, Hsize, MPI_DOUBLE_COMPLEX, 0, comm, ierr)
         ! Measure Occupation of each state
         if (occ.ne.0) then
            if(myid.eq.0) then
               allocate(ni(N))
               write(*,*) "==========================" 
               write(*,*) "Calculate Site Occupations" 
               open(17,file='occupation.dat')
               do kk = 0, Nkpt-1
                  ! There are some problems here with multiple site
                  ! Measure up spin occupation
                  qq = ksubtract(kk,gsmom) !!! gsmom + qq = kk
                  qx = kpoint(qq,1)
                  qy = kpoint(qq,2)
                  nup = nup - 1
                  call GenHsp_kspace(qq, Hsp_f, ksize_f) 
                  write(*,*) 'GenHsp up spin is done, Hsize = ', ksize_f, qq
                  tempv = 0.0d0
                  ctag = 1
                  call GenDensityState(qq, H_0, tempv, Hsp, Hsp_f, Hsize, ksize_f, ctag)
                  ni(qq+1) = ni(qq+1) + DOT_PRODUCT(tempv, tempv)
                  write(*,*) "Spin up:", DOT_PRODUCT(tempv, tempv)
                  write(17,*) qx, qy, DOT_PRODUCT(tempv, tempv)
                  nup = nup + 1
                  ! Measure down spin occupation
                  ndn = ndn - 1
                  Hsp_f = 0
                  ksize_f = 0
                  call GenHsp_kspace(qq, Hsp_f, ksize_f) 
                  write(*,*) 'GenHsp down spin is done, Hsize = ', ksize_f, qq
                  tempv = 0.0d0
                  ctag = 2
                  write(*,*) N, nup, ndn
                  call GenDensityState(qq, H_0, tempv, Hsp, Hsp_f, Hsize, ksize_f, ctag)
                  ni(qq+1) = ni(qq+1) + DOT_PRODUCT(tempv, tempv)
                  write(*,*) "Spin down:", qq, DOT_PRODUCT(tempv, tempv)
                  write(17,*) qx, qy, DOT_PRODUCT(tempv, tempv)
                  ndn = ndn + 1
               enddo
               close(17)
               write(*,*) "==========================" 
               deallocate(ni)
            endif
         endif

        ! Static Spin/Charge Correlation
         if ((static_corr.ne.0)) then
            if(myid.eq.0) write(*,*) "Calculate Static Correlations"
            Hsp_f = 0.0d0
            tempv = 0.0d0
            ! call system("rm ev1_ninj_ks.dat")
            ! call system("rm ev1_sisj_ks.dat")
            if (ee.lt.10) then 
               write(filename, '(A,I1,A)') 'ev',ee,'_sisj_ks.dat'
            else 
               write(filename, '(A,I2,A)') 'ev',ee,'_sisj_ks.dat'
            endif
            open(115, file=filename,access='append')
            if (ee.lt.10) then 
               write(filename, '(A,I1,A)') 'ev',ee,'_ninj_ks.dat'
            else 
               write(filename, '(A,I2,A)') 'ev',ee,'_ninj_ks.dat'
            endif
            open(116, file=filename,access='append')
            if (ee.lt.10) then 
               write(filename, '(A,I1,A)') 'ev',ee,'_docc.dat'
            else 
               write(filename, '(A,I2,A)') 'ev',ee,'_docc.dat'
            endif
            open(117, file=filename,access='append')
            
            ! Double Occupancy <nini>
            ! do  site1 = 1,Nsite
            tempv = CMPLX(0.0d0,0.0d0)
            docc = CMPLX(0.0d0,0.0d0)
            call GenDoubleOccupiedState(H_0, tempv, Hsp, Hsize, int(1,8), Nsite, docc)
            if (myid.eq.0) then
               write(*,*) "Double Occupancy: ", REAL(docc), docc
               write(117,*) REAL(docc)
            endif
            close(117)
            ! enddo
            
            do qq = 0, Nkpt-1
               kk = ksum(gsmom,qq*Nsite) !!! gsmom + qq = kk ==> kk is transferred momentum
               if (myid.eq.0) then
                  call GenHsp_kspace(kk, Hsp_f, ksize_f) 
                  write(*,*) 'GenHsp up spin is done, Hsize = ', ksize_f, qq
                  write(*,*) 'gsmom   : ', gsmom, kpoint(gsmom,1), kpoint(gsmom,2)
                  write(*,*) 'qq      : ', qq*Nsite, kpoint(qq*Nsite,1), kpoint(qq*Nsite,2)
                  write(*,*) 'gsmom+qq: ', kk, kpoint(kk,1), kpoint(kk,2)
               endif
               call MPI_BCAST(ksize_f, 1, MPI_INTEGER8, 0, comm, ierr)
               call MPI_BCAST(Hsp_f, ksize_f, MPI_INTEGER8, 0, comm, ierr)

               ! Spin-Spin <SiSj>
               do site1 = 1,Nsite
                  tempv = 0.0d0 
                  call GenSpinChargeCorrelatedState(qq*Nsite, H_0, tempv, Hsp, Hsp_f, Hsize, ksize_f, int(1,8), site1, site1)
                  do site2 = 1,Nsite
                     tempv1 = 0.0d0
                     sisj = 0
                     call GenSpinChargeCorrelatedState(qq*Nsite, H_0, tempv1, Hsp, Hsp_f, Hsize, ksize_f, int(1,8), site2, site2)
                     do ii = 1,ksize_f
                        sisj = sisj + tempv(ii)*CONJG(tempv1(ii))/4/N
                     enddo
                     ! sisj = ZDOTC(tempv, tempv1)/4/N
                     if(myid.eq.0) then
                        ! write(*,*) 'Static Spin correlation', kpoint(qq*Nsite,1), kpoint(qq*Nsite,2), site1, site2, sisj
                        write(115,*)  kpoint(qq*Nsite,1), kpoint(qq*Nsite,2), site1, site2, REAL(sisj), AIMAG(sisj)
                     endif 
                  enddo
               enddo
               ! Charge-Charge <ninj>-!!!!<ni><nj>!!!! This part not yet calculated
               ! do site1 = 1,Nsite
               !    tempv = 0.0d0 
               !    call GenSpinChargeCorrelatedState(qq*Nsite, H_0, tempv, Hsp, Hsp_f, Hsize, ksize_f, int(2,8), site1, site1)
               !    do site2 = 1,Nsite
               !       tempv1 = 0.0d0
               !       call GenSpinChargeCorrelatedState(qq*Nsite, H_0, tempv1, Hsp, Hsp_f, Hsize, ksize_f, int(2,8), site2, site2)
               !       ninj = DOT_PRODUCT(tempv, tempv1)/N
               !       if(myid.eq.0) then
               !          ! write(*,*) 'Static Charge correlation', kpoint(qq,1), kpoint(qq,2), site1, site2, ninj
               !          write(116,*)  kpoint(qq*Nsite,1), kpoint(qq*Nsite,2), site1, site2, REAL(ninj), AIMAG(ninj)
               !       endif 
               !    enddo
               ! enddo
            enddo
            close(115)
            close(116)
         endif
         
         ! Dynamic Spin/Charge Correlation
         if (dynamic_corr.ne.0) then
            if(myid.eq.0) write(*,*) "Calculating Dynamic Correlations"
            do qqk = 1, Kmapsize
               ! Loop over interesting Kpoints
               Hsp_f = 0.0
               kk = Kmap(qqk) ! Final k point
               qq = ksubtract(kk,gsmom) ! gsmom + qq = kk
               qx = kpoint(qq,1)
               qy = kpoint(qq,2)
               if(myid.eq.0) then 
                  write(*,*) '========================================'
                  write(*,*) 'Calculating Dynamic Correlation for kpoint: ', qx, qy
                  call GenHsp_kspace(qq, Hsp_f, ksize_f) 
                  write(*,*) 'GenHsp is done, Hsize = ', ksize_f, qq
               endif
               call MPI_BCAST(ksize_f, 1, MPI_INTEGER8, 0, comm, ierr)
               call MPI_BCAST(Hsp_f, ksize_f, MPI_INTEGER8, 0, comm, ierr)
               call Setnloc(ksize_f)
               if(myid.eq.0) write(*,*) 'Setnloc is done'
               SprSize_f = HsizeEst/nprocs*nzeprEst
               call GenMatrix_kspace(Hsp_f, ksize_f, IndexI_f, IndexJ_f, sparseH_f, SprSize_f)
               call MPI_Barrier(comm, ierr)
               if(myid.eq.0) write(*,*) 'GenMatrix is done and SprSize = ', SprSize
               call MPI_Barrier(comm, ierr)
               ! Calculate S(q,w)
               if(myid.eq.0) write(*,*) 'Calculating S_q(w)'
               tempv = 0.0d0
               call GenSpinChargeCorrelatedState(qq, H_0, tempv, Hsp, Hsp_f, Hsize, ksize_f, int(1,8), int(1,8), Nsite)
               if (ee.lt.10) then 
                  write(filename, '(A,I1,A,I1,A,I1,A)') 'ev',ee,'_sqw_',qx,'_',qy,'.dat'
               else 
                  write(filename, '(A,I2,A,I1,A,I1,A)') 'ev',ee,'_sqw_',qx,'_',qy,'.dat'
               endif
               call PContFracExpanNRES(ksize_f, SprSize_f, tempv, E_0, IndexI_f, IndexJ_f, sparseH_f, specX, specY)
               open(131, file=filename)
               if(myid.eq.0) then
                  do mm=1, divY+1
                     write(131,*) specX(mm), specY(mm)
                  enddo
               endif
               close(131)
               ! Calculate N(q,w)
               if(myid.eq.0) write(*,*) 'Calculating N_q(w)'
               tempv = 0.0d0
               call GenSpinChargeCorrelatedState(qq, H_0, tempv, Hsp, Hsp_f, Hsize, ksize_f, int(2,8), int(1,8), Nsite)
               if (ee.lt.10) then 
                  write(filename, '(A,I1,A,I1,A,I1,A)') 'ev',ee,'_nqw_',qx,'_',qy,'.dat'
               else 
                  write(filename, '(A,I2,A,I1,A,I1,A)') 'ev',ee,'_nqw_',qx,'_',qy,'.dat'
               endif
               call PContFracExpanNRES(ksize_f, SprSize_f, tempv1, E_0, IndexI_f, IndexJ_f, sparseH_f, specX, specY)
               open(132, file=filename)
               if(myid.eq.0) then
                  do mm=1, divY+1
                     write(132,*) specX(mm), specY(mm)
                  enddo
               endif
               close(132)
            enddo
         endif
         if (eff_raman.ne.0) then
            if(myid.eq.0) write(*,*) "Calculating Effective Raman Operator"
            ! Calculate O(A1g), O(B1g), O(B2g)
            qq = gsmom ! gsmom = qq
            do symm = 0,3
               if ((symm.eq.0 .and. a1g.ne.0) .or. (symm.eq.2 .and. b1g.ne.0)&
               .or. (symm.eq.3 .and. b2g.ne.0)) then
                  if (myid.eq.0) then
                     if (symm.eq.0) then 
                        write(*,*) "O(A1g): ctag (a1g=0,2): ", ctag
                        write(symmname,'(A)') "_OA1g"
                     else if (symm.eq.1) then 
                        ! This is omitted for now
                        write(*,*) "O(A2g): ctag (a2g=12): ", ctag
                        write(symmname,'(A)') "_OA2g"
                     else if (symm.eq.2) then 
                        write(*,*) "O(B1g): ctag (b1g=20,22): ", ctag
                        write(symmname,'(A)') "_OB1g"
                     else if (symm.eq.3) then 
                        write(*,*) "O(B2g): ctag (b2g=30,32): ", ctag
                        write(symmname,'(A)') "_OB2g"
                     endif
                  endif
                  do order = 0, 1
                     tempv = 0.0d0
                     tempv1 = 0.0d0
                     ctag = order*2 + symm*10
                     if (ctag.ne.10) then
                        if (myid.eq.0) write(*,*) "ctag: ", ctag
                        call MPI_Barrier(comm, ierr)
                        call Gen2SpinCorrelatedState(qq, H_0, tempv1, Hsp, Hsp, Hsize, Hsize, ctag)
                        call MPI_Barrier(comm, ierr)
                        call MPI_Allreduce(tempv1, tempv, Hsize, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
                        if (ee.lt.10) then 
                           write(filename, '(A,I1,A5,I1,A,I1,A,I1,A)') 'ev',ee,symmname,order*2,'_',qx,'_',qy,'.dat'
                        else 
                           write(filename, '(A,I1,A5,I1,A,I1,A,I1,A)') 'ev',ee,symmname,order*2,'_',qx,'_',qy,'.dat'
                        endif
                        specX = 0
                        specY = 0
                        open(133, file=filename)
                        call PContFracExpanNRESNgs(Hsize, SprSize, tempv, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
                        if(myid.eq.0) then
                           do mm=1, divY+1
                              write(133,*) specX(mm), specY(mm)
                           enddo
                        endif
                        close(133)   
                     endif         
                  enddo
               else if(myid.eq.0) then
                  write(*,*) "2 spin Correlation unavailable for selected symmetry"
               endif 
            enddo
            
            ! Calculate S3(q,w)
            allocate(tempv1cmplx(HsizeEst))
            allocate(tempvcmplx(HsizeEst))
            tempvcmplx = 0.0d0
            tempv1cmplx = 0.0d0
            qq = gsmom ! gsmom = qq
            if (a2g.ne.0) then
               do order = 2,2
                  ctag = order * 2
                  if (myid.eq.0) write(*,*) "S3qw: ctag (a2g=2,4): ", ctag
                  call Gen3SpinCorrelatedState(qq, H_0, tempv1cmplx, Hsp, Hsp, Hsize, Hsize, ctag)
                  call MPI_Barrier(comm, ierr)
                  !!!!! THERE IS A PROBLEM IN ALLREDUCE
                  call MPI_Allreduce(tempv1cmplx, tempvcmplx, Hsize, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
                  if (ee.lt.10) then 
                     write(filename, '(A,I1,A,I1,A,I1,A,I1,A)') 'ev',ee,'_OA2g',ctag,'c_',qx,'_',qy,'.dat'
                  else 
                     write(filename, '(A,I2,A,I1,A,I1,A,I1,A)') 'ev',ee,'_OA2g',ctag,'c_',qx,'_',qy,'.dat'
                  endif
                  if (myid.eq.0) then
                     do i = 1,Hsize
                        z = z + conjg(tempvcmplx(i))*tempvcmplx(i)
                     enddo
                     write(*,*) "ID, Hsize, sum: ", myid,Hsize,z
                  endif
                  specX = 0
                  specY = 0
                  open(134, file=filename)
                  call PContFracExpanNRES(Hsize, SprSize, tempvcmplx, E_0, IndexI, IndexJ, sparseH, specX, specY)
                  if(myid.eq.0) then
                     do mm=1, divY+1
                        write(134,*) specX(mm), specY(mm)
                     enddo
                  endif
                  close(134)
               end do
            else
               if(myid.eq.0) write(*,*) "chiral Correlation unavailable for selected symmetry"
            endif
         endif
         !!!!! Do we need to find optimal momentum point for lowest GS energy??????
      enddo

      !============== The Step 2.A: Calculate Optical Conductivity ==============
      !! Reload Ground State Eigenvector
      if (diag.eq.0) then
         if (myid.eq.0) then
            call readhdf5cmplx(1,Hsize,H_0)
            open(unit=10102, file='targetev.dat', status='old', action='read',iostat=ferr)
            if (ferr == 0) then
               read(10102,'(1F20.16)') E_0
               write(*,*) 'ED is skipped and E_0 = ', E_0
               close(10102)
            endif
         end if
         call MPI_BCAST(E_0, 1, MPI_DOUBLE_COMPLEX, 0, comm, ierr)
         call MPI_BCAST(H_0, Hsize, MPI_DOUBLE_COMPLEX, 0, comm, ierr)
      end if

      if (opt_cond.ne.0) then
         if(myid.eq.0) then 
            write(*,*) "Measuring optical conductivity" 
         endif
         ! Calculate jx
         allocate(specopt(HsizeEst))
         open(127, file='opt_jx.dat',access='append')
         tempv=0.0d0
         tempv=H_0;
         ctag = 1;
         call GenCurrentState(tempv, Hsp, Hsize, ctag);
         specopt = tempv
         call PContFracExpanNRESNgs(Hsize, SprSize, specopt, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
         if(myid.eq.0) then
            do mm=1, divY+1
               ! write(*,*) specX(mm), specY(mm)
               ! specY(mm) = specY(mm)/specX(mm)
               specY(mm) = specY(mm)
               write(127,*) specX(mm), specY(mm)
            enddo
         endif
         ! Calculate jy
         open(128, file='opt_jy.dat',access='append')
         tempv=0.0d0
         tempv=H_0;
         ctag = 2;
         call GenCurrentState(tempv, Hsp, Hsize, ctag);
         if(myid.eq.0) then
            tsum = 0
            do mm = 1,Hsize
               tsum = tsum + tempv(mm)*H_0(mm)
            enddo
            ! write(*,*) "tsum 1: ", tsum, tsum*tsum
            tsum = 0
            do mm = 1,Hsize
               tsum = tsum + tempv(mm)*tempv(mm)
            enddo
            ! write(*,*) "tsum 2: ", tsum, tsum*tsum
            write(*,*) "zero frequency optical conductivity", DOT_PRODUCT(tempv, H_0)*DOT_PRODUCT(tempv, H_0)
         endif 
         specopt = tempv
         call PContFracExpanNRESNgs(Hsize, SprSize, specopt, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
         if(myid.eq.0) then
            do mm=1, divY+1
               ! write(*,*) specX(mm), specY(mm)
               specY(mm) = specY(mm)/specX(mm)
               write(128,*) specX(mm), specY(mm)
            enddo
         endif

         close(127)
         close(128)
      endif
      !============== The Step 2.B: Calculate NRES ==============
      if (nres_raman.ne.0) then
         if(myid.eq.0) then 
            write(*,*) "Calculate non-resonant Raman" 
         endif
         allocate(specraman(HsizeEst))
         tempv=0.0d0
         ! Forgive me using GenCurrentState, it's easier
         ! Calculate A1g cos(kx)+cos(ky) 
         if (a1g.ne.0) then
            open(124, file='NRESRamanA1g.dat',access='append')
            tempv=H_0;
            ctag = 3;
            call GenCurrentState(tempv, Hsp, Hsize, ctag);
            if(myid.eq.0) then
               tsum = 0
               do mm = 1,Hsize
                  tsum = tsum + tempv(mm)*H_0(mm)
               enddo
               write(*,*) "tsum A1g NRES: ", tsum, tsum*tsum
            endif
            specraman = tempv
            call PContFracExpanNRESNgs(Hsize, SprSize, specraman, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
            ! call PContFracExpanNRES(Hsize, SprSize, specraman, E_0, IndexI, IndexJ, sparseH, specX, specY)
            if(myid.eq.0) then
               do mm=1, divY+1
                  write(124,*) specX(mm), specY(mm)
               enddo
            endif
            close(124)
         endif

         ! Calculate B1g cos(kx)-cos(ky) 
         if (b1g.ne.0) then
            tempv=H_0;
            ctag = 4;
            call GenCurrentState(tempv, Hsp, Hsize, ctag);
            if(myid.eq.0) then
               if (twist1.ne.0.or.twist2.ne.0) then
                  open(unit=125, file='NRESRamanB1g.dat', status='old', action='read',iostat=ferr)
                  do mm=1, divY+1
                     read(125,*) temp_read, specY_read(mm)
                  enddo
                  close(125)
               else
                  specY_read = 0.0d0
               endif
            endif
            specraman = tempv
            call PContFracExpanNRESNgs(Hsize, SprSize, specraman, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
            if(myid.eq.0) then
               ! Delete NRES file
               filename = 'NRESRamanE21.dat'
               if (Ndim.eq.2) filename = 'NRESRamanB1g.dat'
               write(*,*) "B1g NRES finished, writing file..."
               open(125, file=filename, status='old', iostat=ferr)
               if (ferr == 0) close(125, status='delete')
               open(125, file=filename, access='append', action='write')
               do mm=1, divY+1
                  specY(mm) = specY(mm) + specY_read(mm)
                  if (twist1.eq.(twistx-1).and.twist2.eq.(twisty-1)) specY(mm) = specY(mm)/dfloat(twistx*twisty)
                  write(125,*) specX(mm), specY(mm)
               enddo
               close(125)
            endif
         endif
         ! Calculate B2g cos(kx)-cos(ky) 
         if (b2g.ne.0) then
            filename = 'NRESRamanE22.dat'
            if (Ndim.eq.2) filename = 'NRESRamanB2g.dat'
            open(126, file=filename,access='append')
            tempv=H_0;
            ctag = 5;
            call GenCurrentState(tempv, Hsp, Hsize, ctag);
            if(myid.eq.0) then
               tsum = 0
               do mm = 1,Hsize
                  tsum = tsum + tempv(mm)*H_0(mm)
               enddo
               write(*,*) "tsum B2g NRES: ", tsum, tsum*tsum
            endif
            specraman = tempv
            call PContFracExpanNRESNgs(Hsize, SprSize, specraman, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
            if(myid.eq.0) then
               do mm=1, divY+1
                  write(126,*) specX(mm), specY(mm)
               enddo
            endif
            close(126)
         endif
         deallocate(specraman)
      endif
      !============== The Step 2.C: Calculate RES ==============
      if (res_raman.ne.0) then
         if(myid.eq.0) then 
            write(*,*) "Calculate resonant Raman" 
         endif
         allocate(specraman(Hsize))
         tempv=0.0d0
         tempv1=0.0d0
         do eneY=eneY_start, divY ! loop over intermediate state energy
            z = CMPLX(dble(eneY)/dble(divY)*(endY-startY)+startY+E_0, epsilone_CG) 
            if(myid.eq.0) then
               write(*,*) '*************************************'
               write(*,*) '  incident energy index =  ', dble(eneY)/dble(divY)*(endY-startY)+startY
               write(*,*) '*************************************'
            endif 

            if (a2g.ne.0.or.b2g.ne.0) then
               ! j_y j_x ! Generate <jx|0) and <j_y|0)
               tempv=H_0;
               ctag = 1;
               call GenCurrentState(tempv, Hsp, Hsize, ctag); ! tag = 1 for x; tag = 2 for y
               call PBiCGS(z, SprSize, tempv, IndexI, IndexJ, sparseH);
               ctag = 2;
               call GenCurrentState(tempv, Hsp, Hsize, ctag); ! tag = 1 for x; tag = 2 for y
               ! j_x j_y
               tempv1=H_0;
               ctag = 2;
               call GenCurrentState(tempv1, Hsp, Hsize, ctag); ! tag = 1 for x; tag = 2 for y
               call PBiCGS(z, SprSize, tempv1, IndexI, IndexJ, sparseH);
               ctag = 1;
               call GenCurrentState(tempv1, Hsp, Hsize, ctag); ! tag = 1 for x; tag = 2 for y
            endif
            
            ! Raman A2g: j_x j_y - j_y j_x
            if (a2g.ne.0) then
               ! Read in file
               if(myid.eq.0) then
                  if (twist1.ne.0.or.twist2.ne.0) then
                     open(unit=130, file='RamanA2g.dat', status='old', action='read',iostat=ferr)
                     ! SKIP irrelevant energy
                     if (eneY.ne.eneY_start) then
                        do mm=0,divX*(eneY-eneY_start)
                           read(130,*)
                        enddo
                     endif
                     do mm=1, divX+1
                        read(130,*) temp_read, temp_read, specY_read(mm)
                     enddo
                     close(130)
                  else
                     specY_read = 0.0d0
                  endif
               endif
               specraman = tempv1 - tempv;
               call PContFracExpanNgs(Hsize, SprSize, specraman, H_0, E_0, IndexI, &
                                       IndexJ, sparseH, specX, specY)
               if(myid.eq.0) then
                  open(120, file='RamanA2g-tmp.dat',access='append')
                  write(*,*) 'ready to write A2g'
                  do mm=1, divX+1
                     specY(mm) = specY(mm) + specY_read(mm)
                     if (twist1.eq.(twistx-1).and.twist2.eq.(twisty-1)) specY(mm) = specY(mm)/dfloat(twistx*twisty)
                     write(120,*) dble(eneY)/dble(divY)*(endY-startY)+startY, specX(mm), specY(mm)
                  enddo
                  close(120)
                  if (eneY.eq.divY) call system("mv RamanA2g-tmp.dat RamanA2g.dat")
               endif
            endif 

            ! Raman B2g: j_x j_y + j_y j_x
            if (b2g.ne.0) then
               filename = 'RamanE22.dat'
               if (Ndim.eq.2) filename = 'RamanB2g.dat'
               open(121, file=filename,access='append')
               specraman = tempv + tempv1;
               call PContFracExpanNgs(Hsize, SprSize, specraman, H_0, E_0, IndexI, &
                  IndexJ, sparseH, specX, specY)
               ! call PContFracExpan(Hsize, SprSize, specraman, E_0, IndexI, &
               !    IndexJ, sparseH, specX, specY)
               if(myid.eq.0) then
                  write(*,*) 'ready to write B2g/Eg2'
                  do mm=1, divX+1
                     write(121,*) dble(eneY)/dble(divY)*(endY-startY)+startY, specX(mm), specY(mm)
                  enddo
               endif
               close(121)
            endif

            if (a1g.ne.0.or.b1g.ne.0) then
               ! j_x j_x
               tempv=H_0;
               ctag = 1;
               call GenCurrentState(tempv, Hsp, Hsize, ctag); ! tag = 1 for x; tag = 2 for y
               call PBiCGS(z, SprSize, tempv, IndexI, IndexJ, sparseH);
               ctag = 1;
               call GenCurrentState(tempv, Hsp, Hsize, ctag); ! tag = 1 for x; tag = 2 for y
               
               ! j_y j_y
               tempv1=H_0;
               ctag = 2;
               call GenCurrentState(tempv1, Hsp, Hsize, ctag); ! tag = 1 for x; tag = 2 for y
               call PBiCGS(z, SprSize, tempv1, IndexI, IndexJ, sparseH);
               ctag = 2;
               call GenCurrentState(tempv1, Hsp, Hsize, ctag); ! tag = 1 for x; tag = 2 for y
            endif
            
            ! Raman B1g: j_x j_x - j_y j_y
            if (b1g.ne.0) then
               ! Read in file
               if(myid.eq.0) then
                  if (twist1.ne.0.or.twist2.ne.0) then
                     open(unit=133, file='RamanB1g.dat', status='old', action='read',iostat=ferr)
                     ! SKIP irrelevant energy
                     if (eneY.ne.eneY_start) then
                        do mm=0,divX*(eneY-eneY_start)
                           read(133,*)
                        enddo
                     endif
                     do mm=1, divX+1
                        read(133,*) temp_read, temp_read, specY_read(mm)
                     enddo
                     close(133)
                  else
                     specY_read = 0.0d0
                  endif
               endif
               specraman = tempv1 - tempv;
               call PContFracExpanNgs(Hsize, SprSize, specraman, H_0, E_0, IndexI, &
                  IndexJ, sparseH, specX, specY)
               if(myid.eq.0) then
                  filename = 'RamanE21-tmp.dat'
                  write(*,*) 'ready to write B1g/Eg1'
                  if (Ndim.eq.2) filename = 'RamanB1g-tmp.dat'
                  open(123, file=filename,access='append')
                  do mm=1, divX+1
                     specY(mm) = specY(mm) + specY_read(mm)
                     if (twist1.eq.(twistx-1).and.twist2.eq.(twisty-1)) specY(mm) = specY(mm)/dfloat(twistx*twisty)
                     write(123,*) dble(eneY)/dble(divY)*(endY-startY)+startY, specX(mm), specY(mm)
                  enddo
                  close(123)
                  if (eneY.eq.divY) then
                     if (Ndim.eq.2) then
                        call system("mv RamanB1g-tmp.dat RamanB1g.dat")
                     else 
                        call system("mv RamanE21-tmp.dat RamanE21.dat")
                     endif
                  endif
               endif
            endif

            ! Raman A1g: j_x j_x + j_y j_y
            if (a1g.ne.0) then
               open(122, file='RamanA1g.dat',access='append')
               specraman = tempv + tempv1;
               call PContFracExpanNgs(Hsize, SprSize, specraman, H_0, E_0, IndexI, &
                  IndexJ, sparseH, specX, specY)
               if(myid.eq.0) then
                  write(*,*) 'ready to write A1g'
                  do mm=1, divX+1
                     ! write(*,*) dble(eneY)/dble(divY)*(endY-startY)+startY, specX(mm), specY(mm)
                     write(122,*) dble(eneY)/dble(divY)*(endY-startY)+startY, specX(mm), specY(mm)
                  enddo
               endif
               close(122)
            endif
         enddo
         deallocate(specraman)
      endif
   enddo
enddo
end program

subroutine BiCGS(z, Hsizet, SprSizet, groundH, IndexI, IndexJ, sparseH)
use NumOfOrbitalAndElectrons; use ConstantParas; use ScanRegion
implicit none

double complex, INTENT(IN) :: z
integer*8, INTENT(IN) :: SprSizet, Hsizet
double complex, DIMENSION(Hsizet), INTENT(INOUT) :: groundH
integer*8, DIMENSION(SprSizet), INTENT(IN) :: IndexI, IndexJ
double complex, DIMENSION(SprSizet), INTENT(IN) :: sparseH

double complex, allocatable:: x0(:), xn(:)
double complex, allocatable:: r0(:), rn(:), r00(:)
double complex, allocatable:: p0(:), pn(:)
double complex, allocatable:: v0(:), vn(:)
double complex, allocatable:: s(:), at(:)
double complex, allocatable:: dvtemp(:)
double complex :: alpha, rho0, rhon
double complex :: myrhon, mytsum, tsum, mytsum1, tsum1
double complex :: betha, omega0, omegan
integer*8 :: spr_size, nloc, ii, jj
double precision :: time1, time2
     !z = (0.0d0, 0.0d0)
   !****************        STEP 3.A        *****************
   !********              Initialization           **********

nloc = Hsizet
spr_size = SprSizet

allocate(x0(nloc),xn(nloc))
allocate(p0(nloc),pn(nloc))
allocate(v0(nloc),vn(nloc))
allocate(r0(nloc),rn(nloc),r00(nloc))
allocate(s(nloc),at(nloc))
allocate(dvtemp(nloc))

  x0(:)=0.0d0;
 ! r0(:)=groundH(:)-A(:,:)*x0(:);
  r0(:)=groundH(:);
  r00(:)=r0(:);
  rho0=1.0d0;
  alpha=1.0d0;
  omega0=1.0d0;
  v0(:)=0.0d0;
  p0(:)=0.0d0;

!if( myid.eq.0) then
!do ii=1,nloc
!   write(*,*) groundH(ii)
!enddo
!endif

   !****************        STEP 3.B        *****************
   !*********         Conjugate Gradient     *********************

call CPU_time(time1)

  do ii=1,niter_CG  !CG loop
   rhon = DOT_PRODUCT(r00,r0);
   !rhon = (0.0d0,0.0d0)
   !call MPI_Barrier(comm, rc)
   !call MPI_Allreduce(myrhon, rhon, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   betha = (rhon/rho0)*(alpha/omega0);

!   write(*,*) 'rhon', rhon
!   write(*,*) 'rho0', rho0
!   write(*,*) 'betha', betha

   pn(1:nloc) = r0(1:nloc) + betha*(p0(1:nloc)-omega0*v0(1:nloc))
   !call MPI_Barrier(comm, rc)
   !call MPI_Gather(pn,nloc,MPI_DOUBLE_COMPLEX,vtemp,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Gather(pn,nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(vtemp,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(phi,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)

   !call MPI_Allgather(pn,nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,comm,rc)

   vn(:)=(0.0d0,0.0d0)
   do jj=1, spr_size
      vn(IndexI(jj))=vn(IndexI(jj))-pn(IndexJ(jj))*sparseH(jj)
   enddo
!   do jj=1, spr_size_CH(CHsite)
!      vn(rIndex(jj,CHsite)-localstart+1)=vn(rIndex(jj,CHsite)-localstart+1)-pn(rIndex(jj,CHsite)-localstart+1)*rsparseH(jj,CHsite)
!   enddo
   do jj=1,nloc
      vn(jj)=vn(jj)+pn(jj)*z
   enddo
   tsum = DOT_PRODUCT(r00,vn)
   !call MPI_Barrier(comm, rc)
   !call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   alpha = rhon/tsum

!   write(*,*) 'a rhon', rhon
!   write(*,*) 'a tsum', tsum

   s(1:nloc) = r0(1:nloc) - alpha*vn(1:nloc)
   at(:)=(0.0d0,0.0d0)
   !call MPI_Allgather(s,nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,comm,rc)

   do jj=1, spr_size
      at(IndexI(jj))=at(IndexI(jj))-s(IndexJ(jj))*sparseH(jj)
   enddo
!   do jj=1, spr_size_CH(CHsite)
!      at(rIndex(jj,CHsite)-localstart+1)=at(rIndex(jj,CHsite)-localstart+1)-s(rIndex(jj,CHsite)-localstart+1)*rsparseH(jj,CHsite)
!   enddo
   do jj=1,nloc
      at(jj)=at(jj)+s(jj)*z
   enddo
   tsum = DOT_PRODUCT(at,s)
   tsum1 = DOT_PRODUCT(at,at)
   !call MPI_Barrier(comm, rc)
   !call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   !call MPI_Allreduce(mytsum1, tsum1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   omegan = tsum/tsum1
   xn(:) = x0(:)+alpha*pn(:)+omegan*s(:)
   rn(:) = s(:) - omegan*at(:)

!   write(*,*) 'alpha', alpha
!   write(*,*) 'betha', betha
!   write(*,*) 'omegan', omegan

 ! xn=x0+alpha*pn+omegan*s;
  ! dvtemp(1:nloc)= x0(1:nloc)+alpha*vn(1:nloc)+omegan*at(1:nloc)-groundH(1:nloc);
   dvtemp(:)=0.0d0
   !call MPI_Barrier(comm, rc)
   !call MPI_Gather(xn(1),nloc,MPI_DOUBLE_COMPLEX,vtemp,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Gather(xn(1),nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(vtemp,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Bcast(phi,ReducedHsize,MPI_DOUBLE_COMPLEX,0,comm,rc)
   !call MPI_Allgather(xn(1),nloc,MPI_DOUBLE_COMPLEX,phi,nloc,MPI_DOUBLE_COMPLEX,comm,rc)

   do jj=1, spr_size
      dvtemp(IndexI(jj))=dvtemp(IndexI(jj))-xn(IndexJ(jj))*sparseH(jj)
   enddo
!   do jj=1, spr_size_CH(CHsite)
!      dvtemp(rIndex(jj,CHsite)-localstart+1)=dvtemp(rIndex(jj,CHsite)-localstart+1)-xn(rIndex(jj,CHsite)-localstart+1)*rsparseH(jj,CHsite)
!   enddo
   do jj=1,nloc
      dvtemp(jj)=dvtemp(jj)+xn(jj)*z
   enddo
   do jj=1,nloc
      dvtemp(jj)=dvtemp(jj)-groundH(jj)
   enddo

   tsum = (0.0d0,0.0d0)
   do jj=1,nloc
      tsum = tsum + CONJG(dvtemp(jj))*dvtemp(jj)
   enddo
   !call MPI_Barrier(comm, rc)
   !call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, rc)
   write(*,*) ii, tsum
   if(real(tsum).lt.tol_CG) then
      exit
   endif

   x0(:)=xn(:)
   rho0=rhon;
   p0(:)=pn(:)
   r0(:)=rn(:)
   omega0=omegan;
   v0(:)=vn(:)

  enddo  !CG loop

call CPU_time(time2)
 write(*,*) 'CG finishes at:', time2, 'secs for ', ii,' iters'
 if(ii.gt.niter_CG) write(*,*) 'BLOWUP TOO BAD', real(tsum)

groundH(:)=xn(:)

deallocate(x0,xn,p0,pn,v0,vn,r0,rn,r00,s,at,dvtemp)

end subroutine

subroutine sparse_symmetric(SprSizet,IndexIt,IndexJt,sparseHt,ksize)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas;
use MPIParas; use mpi;
implicit none
! Use with caution
integer*8, INTENT(IN) :: SprSizet
integer*8, DIMENSION(SprSizet), INTENT(IN) :: IndexIt, IndexJt
double complex, DIMENSION(SprSizet), INTENT(IN) :: sparseHt
integer*8, INTENT(IN) :: ksize
double complex :: mat_loc(ksize,ksize), mat(ksize,ksize)
integer :: kk,ii,jj

mat = 0.0d0
write(*,*) "Checking Hermitianess: ", ksize, SprSizet
if (ksize.le.1000) then
   do kk = 1,SprSizet
      mat_loc(indexIt(kk),indexJt(kk)) = mat_loc(indexIt(kk),indexJt(kk)) +sparseHt(kk)
   enddo 
   call MPI_Barrier(comm, ierr)
   call MPI_Allreduce(mat_loc, mat, ksize*ksize, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
   if (myid.eq.0) then
      do ii = 1,ksize
         do jj = ii+1,ksize
            if (abs(mat(ii,jj) - conjg(mat(jj,ii))).ge.1E-6) then
               write(*,*) "Matrix not Hermitian: "
               write(*,*) "Entry=",ii,jj,mat(ii,jj)
               write(*,*) "Entry=",jj,ii,mat(jj,ii)
               stop
            endif
         enddo
      enddo
      call system("rm mat.dat")
      open(unit=133, file='mat.dat',Access = 'append')
      do ii = 1,ksize
         do jj = 1,ksize
            write(133,*) REAL(mat(ii,jj)),IMAG(mat(ii,jj))
         enddo
      enddo
      close(133)
   endif
else
   if (myid.eq.0) write(*,*) "Matrix too large, not checking if it's symmetric"
endif
end subroutine

subroutine sparse_trace(SprSizet,IndexIt,IndexJt,sparseHt,trace)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas;
use MPIParas; use mpi;
implicit none
integer*8, INTENT(IN) :: SprSizet
integer*8, DIMENSION(SprSizet), INTENT(IN) :: IndexIt, IndexJt
double complex, DIMENSION(SprSizet), INTENT(IN) :: sparseHt
double complex, INTENT(OUT) :: trace
double complex :: trace_loc
integer :: kk

trace_loc = 0
trace = 0
do kk = 1,SprSizet
   if (IndexIt(kk).eq.IndexJt(kk).and.ABS(sparseHt(kk)).ge.1E-5) then
      ! write(*,*) IndexIt(kk), IndexJt(kk), sparseHt(kk)
      trace_loc = trace_loc + sparseHt(kk)
   endif
enddo 
call MPI_Barrier(comm, ierr)
call MPI_Allreduce(trace_loc, trace, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
end subroutine

subroutine GenDensityState(qq, groundH, midH, listqpt0, listqpt, &
        ksize0, ksize, ctag)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas;
use ModelParas
implicit none
! 
!   Calculate |n_k|phi> for each site
!   Spin up: ctag = 1, Spin down: ctag = 2
INTEGER*8, INTENT(IN) :: qq, ksize, ksize0, ctag
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double complex, DIMENSION(ksize0), INTENT(INOUT) :: groundH
double complex, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit
integer*8 :: jj, l
integer*8 :: gup, gdn, mup, mdn, midsign
double precision :: ptx, pty
write(*,*) "ksize f", ksize
do jj=1,ksize0           !!!!!change from groundstate basis to first excited state basis
   gup=listqpt0(jj)/2**N
   gdn=mod(listqpt0(jj),2**N)
   if (ctag.eq.1) then
      if(BTEST(gup,qq)) then 
         mup=IBCLR(gup,qq)
         mdn=gdn
         midsign=(-1)**(sumbeforebit(mup,qq)+sumeverybit(mdn))
         call BinarySearch(listqpt,ksize,mup*(2**N)+mdn,l)
         if(l.eq.-1) then
            write(*,*) "state: ", jj, mup*(2**N)+mdn
            write(*,*) "bits: ", POPCNT(mup*(2**N)+mdn), POPCNT(mup*(2**N)), POPCNT(mdn)
            write(*,*) 'No overlap between groundstate and N-1 state'
            write(*,*) 'kk=',qq
            stop
         endif
         midH(l)=midH(l)+midsign*groundH(jj) 
      endif
   endif 
   if (ctag.eq.2) then
      if(BTEST(gdn,qq)) then 
         mup=gup
         mdn=IBCLR(gdn,qq)
         midsign=(-1)**(sumbeforebit(mdn,qq))
         call BinarySearch(listqpt,ksize,mup*(2**N)+mdn,l)
         if(l.eq.-1) then
            write(*,*) "state: ", mup*(2**N)+mdn
            write(*,*) "bits: ", POPCNT(mup*(2**N)), mdn
            write(*,*) 'No overlap between groundstate and N-1 state'
            write(*,*) 'kk=',qq
            stop
         endif
         midH(l)=midH(l)+midsign*groundH(jj) 
      endif
   endif
enddo

end subroutine

subroutine GenCurrentState(v0, listqpt, ksize, ctag) ! tag = 1 for x; tag = 2 for y
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas;
use ModelParas; 
use MPIParas; use mpi;
implicit none
! 
!   j | n > = \sum_{k\sigma} 2t sink_x c^{\dagger}_{k\sigma}c_{k\sigma} 
!
double complex, DIMENSION(ksize), INTENT(INOUT) :: v0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
INTEGER*8, INTENT(IN) :: ksize
integer*8, INTENT(IN) :: ctag

integer*8 :: jj, kk, kkp
integer*8 :: kks, qqs, l
integer*8 :: gup, gdn, site1, site2
integer*8 :: initstatup,initstatdn
integer*8 :: afterstatup,afterstatdn,initsign,aftersign
double complex, allocatable :: midH(:), finalH(:)
double complex :: prefac(0:25,3)
double precision :: ptx, pty
integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum

prefac = 0.0d0
do kk=0,Nkpt-1
   if (Ndim.eq.2) then
      ptx = kpoint(kk*Nsite,1)*2.0d0*pi/kunit + dfloat(twist1)*2.0d0*pi/twistx;
      pty = kpoint(kk*Nsite,2)*2.0d0*pi/kunit + dfloat(twist2)*2.0d0*pi/twisty;
   else
      ptx = kpoint(kk*Nsite,1)*2.0d0*pi/kunit + twist1*2.0d0*pi/twistx;
      pty = (kpoint(kk*Nsite,1)*0.5d0+kpoint(kk*Nsite,2))*4.0d0/sqrt(3.0)*pi/kunit + twist2*2.0d0*pi/twisty;
   endif
   selectcase(Ndim)
   !!!!! In band dispersion, there is always a factor of -2t in front of the matrix
   case(1)
      selectcase (ctag)
      case(1) !jx
         prefac(kk,1) = t*(sin(ptx) + sin(ptx/2)*cos(sqrt(3.0)*pty/2)) &
                        + 3*tt*(sin(3*ptx/2)*cos(sqrt(3.0)*pty/2))
      case(2) !jy
         prefac(kk,1) = t*(sqrt(3.0)*cos(ptx/2)*sin(sqrt(3.0)*pty/2)) &
                        + sqrt(3.0)*tt*(sin(sqrt(3.0)*pty)+cos(ptx*3/2)*sin(pty*sqrt(3.0)/2))
      case(3) !A1g xx+yy
         prefac(kk,1) = t*(2*cos(ptx)+4*cos(ptx/2)*cos(sqrt(3.0)*pty/2)) &
                        + 6*tt*(cos(sqrt(3.0)*pty)+2*cos(3*ptx/2)*cos(sqrt(3.0)*pty/2))
      case(4) !B1g xx-yy
         prefac(kk,1) = t*(2*cos(ptx)-2*cos(ptx/2)*cos(sqrt(3.0)*pty/2)) &
                        ! + 3*tt*(cos(3*ptx/2)*cos(sqrt(3.0)*pty/2)-cos(sqrt(3.0)*pty))
                        + 6*tt*(cos(3*ptx/2)*cos(sqrt(3.0)*pty/2)-cos(sqrt(3.0)*pty))
      case(5) !B2g xy+xy
         prefac(kk,1) = -2*sqrt(3.0)*t*sin(ptx/2)*sin(sqrt(3.0)*pty/2) &
                        -6*sqrt(3.0)*tt*sin(3*ptx/2)*sin(sqrt(3.0)*pty/2)
      case DEFAULT
         ! if (myid.eq.0) write(*,*) "invalid tag in Gen Current"
         stop
      endselect
   case(2)
      selectcase(ctag)
      case(1) !jx Current operators are off by a factor of 2
         prefac(kk,1) = sin(ptx) + 2*ttprimeratio*sin(ptx)*cos(pty) &
                        + 2*tttprimeratio*sin(2*ptx)
      case(2) !jy
         prefac(kk,1) = sin(pty) + 2*ttprimeratio*cos(ptx)*sin(pty)&
                        + 2*tttprimeratio*sin(2*pty)
      case(3) !A1g
         prefac(kk,1) = 0.5*(cos(ptx)+cos(pty))
      case(4) !B1g
         prefac(kk,1) = 0.5*(cos(ptx)-cos(pty))
      case(5) !B2g
         prefac(kk,1) = sin(ptx)*sin(pty)
      case DEFAULT
         ! if (myid.eq.0) write(*,*) "invalid tag in Gen Current"
         stop
      endselect
   case(3)
      selectcase(ctag)
      case(1) ! jx I guess it's 0...?
         prefac(kk*Nsite,2)   = exp(CMPLX(0,-pty/sqrt(3.0)/2)) * sin(ptx/2)
         prefac(kk*Nsite+1,1) = exp(CMPLX(0,pty/sqrt(3.0)/2)) * sin(ptx/2)
      case(2) !jy
         prefac(kk*Nsite,2)   =  CMPLX(0,1/sqrt(3.0)) * exp(CMPLX(0,pty/sqrt(3.0))) &
                                    * (-1+cos(ptx/2)*exp(CMPLX(0,-pty*sqrt(3.0))/2))
         prefac(kk*Nsite+1,1) = -CMPLX(0,1/sqrt(3.0)) * exp(CMPLX(0,-pty/sqrt(3.0))) * &
                                 (-1+cos(ptx/2)*exp(CMPLX(0,pty*sqrt(3.0))/2))
      case(3) !A1g
         prefac(kk*Nsite,2)   = 2.0/3.0 * exp(CMPLX(0,pty/sqrt(3.0))) * (1+cos(ptx/2)*exp(CMPLX(0,-pty*sqrt(3.0))/2))
         prefac(kk*Nsite+1,1) = 2.0/3.0 * exp(CMPLX(0,-pty/sqrt(3.0))) * (1+cos(ptx/2)*exp(CMPLX(0,pty*sqrt(3.0))/2))
      case(4) !B1g
         if (myid.eq.0) then
            write(*,*) ptx, pty
            write(*,*) exp(CMPLX(0,pty/sqrt(3.0))) * (-2+cos(ptx/2)*exp(CMPLX(0,-pty*sqrt(3.0))/2))
            write(*,*) exp(CMPLX(0,-pty/sqrt(3.0))) * (-2+cos(ptx/2)*exp(CMPLX(0,pty*sqrt(3.0))/2))
         endif
         prefac(kk*Nsite,2)   = 1.0/3.0 * exp(CMPLX(0,pty/sqrt(3.0))) * (-2+cos(ptx/2)*exp(CMPLX(0,-pty*sqrt(3.0))/2))
         prefac(kk*Nsite+1,1) = 1.0/3.0 * exp(CMPLX(0,-pty/sqrt(3.0))) * (-2+cos(ptx/2)*exp(CMPLX(0,pty*sqrt(3.0))/2))
      case(5) !B2g
         prefac(kk*Nsite,2)   = -CMPLX(0,1/sqrt(3.0)) * exp(CMPLX(0,-pty/sqrt(3.0)/2)) * sin(ptx/2)
         prefac(kk*Nsite+1,1) =  CMPLX(0,1/sqrt(3.0)) * exp(CMPLX(0,pty/sqrt(3.0)/2)) * sin(ptx/2)
      case DEFAULT
         ! if (myid.eq.0) write(*,*) "invalid tag in Gen Current"
         stop
      endselect
      if (myid.eq.0) then
         if (ctag.eq.1) write(*,*) "jx"
         if (ctag.eq.2) write(*,*) "jy"
         write(*,*) "Ksite: ", kk*Nsite, kk*Nsite+1
         write(*,*) "kpoint1, disp: ", ptx,pty,prefac(kk*Nsite,1),prefac(kk*Nsite,2)
         write(*,*) "kpoint2, disp: ", ptx,pty,prefac(kk*Nsite+1,1),prefac(kk*Nsite+1,2)
         write(*,*) "---------------------------------------------------------------------------------"
      endif
   case(4)
      selectcase(ctag)
      case(1) 
         prefac(kk*Nsite,2)   = -0.5*sin(-ptx/4+sqrt(3.0)/4*pty)
         prefac(kk*Nsite,3)   =  sin(ptx/2)
         prefac(kk*Nsite+1,1) = -0.5*sin(-ptx/4+sqrt(3.0)/4*pty)
         prefac(kk*Nsite+1,3) =  0.5*sin(ptx/4+sqrt(3.0)/4*pty)
         prefac(kk*Nsite+2,1) =  sin(ptx/2)
         prefac(kk*Nsite+2,2) =  0.5*sin(ptx/4+sqrt(3.0)/4*pty)
      case(2) !jy
         prefac(kk*Nsite,2)   =  sqrt(3.0)/2*sin(-ptx/4+sqrt(3.0)/4*pty)
         prefac(kk*Nsite+1,1) =  sqrt(3.0)/2*sin(-ptx/4+sqrt(3.0)/4*pty)
         prefac(kk*Nsite+1,3) =  sqrt(3.0)/2*sin( ptx/4+sqrt(3.0)/4*pty)
         prefac(kk*Nsite+2,2) =  sqrt(3.0)/2*sin( ptx/4+sqrt(3.0)/4*pty)
      case(3) !A1g
         prefac(kk*Nsite,2)   = 0.5*cos(ptx/4-sqrt(3.0)/4*pty)
         prefac(kk*Nsite,3)   = 0.5*cos(ptx/2)
         prefac(kk*Nsite+1,1) = 0.5*cos(ptx/4-sqrt(3.0)/4*pty)
         prefac(kk*Nsite+1,3) = 0.5*cos(ptx/4+sqrt(3.0)/4*pty)
         prefac(kk*Nsite+2,1) = 0.5*cos(ptx/2)
         prefac(kk*Nsite+2,2) = 0.5*cos(ptx/4+sqrt(3.0)/4*pty)
      case(4) !B1g
         prefac(kk*Nsite,2)   = -0.25*cos(ptx/4-sqrt(3.0)/4*pty)
         prefac(kk*Nsite,3)   =  0.5*cos(ptx/2)
         prefac(kk*Nsite+1,1) = -0.25*cos(ptx/4-sqrt(3.0)/4*pty)
         prefac(kk*Nsite+1,3) = -0.25*cos(ptx/4+sqrt(3.0)/4*pty)
         prefac(kk*Nsite+2,1) =  0.5*cos(ptx/2)
         prefac(kk*Nsite+2,2) = -0.25*cos(ptx/4+sqrt(3.0)/4*pty)
      case(5) !B2g
         prefac(kk*Nsite,2)   = -sqrt(3.0)/4*cos(ptx/4-sqrt(3.0)/4*pty)
         prefac(kk*Nsite+1,1) = -sqrt(3.0)/4*cos(ptx/4-sqrt(3.0)/4*pty)
         prefac(kk*Nsite+1,3) =  sqrt(3.0)/4*cos(ptx/4+sqrt(3.0)/4*pty)
         prefac(kk*Nsite+2,2) =  sqrt(3.0)/4*cos(ptx/4+sqrt(3.0)/4*pty)
      case DEFAULT
         ! if (myid.eq.0) write(*,*) "invalid tag in Gen Current"
         stop
      endselect
      if (myid.eq.0) then
         if (ctag.eq.1) write(*,*) "jx"
         if (ctag.eq.2) write(*,*) "jy"
         write(*,*) "Ksite: ", kk*Nsite, kk*Nsite+1, kk*Nsite+2
         write(*,*) "kpoint1, disp: ", ptx,pty,prefac(kk*Nsite,1),prefac(kk*Nsite,2),prefac(kk*Nsite,3)
         write(*,*) "kpoint2, disp: ", ptx,pty,prefac(kk*Nsite+1,1),prefac(kk*Nsite+1,2),prefac(kk*Nsite+1,3)
         write(*,*) "kpoint3, disp: ", ptx,pty,prefac(kk*Nsite+2,1),prefac(kk*Nsite+2,2),prefac(kk*Nsite+2,3)
         write(*,*) "---------------------------------------------------------------------------------"
      endif
   endselect
enddo

allocate(midH(ksize))
midH = 0.0
! do jj= 1,ksize  !loop B
do jj = localstart,localend
   initstatup=listqpt(jj)/2**N;
   initstatdn=mod(listqpt(jj),2**N);
   ! if (myid.eq.1) WRITE(*,'(B64)') initstatup, initstatdn
   do kk=0,Nkpt-1
      ! Insert dispersion relation here
      do site1=1,Nsite
         do site2=1,Nsite
            kks=kk*Nsite+site1-1
            qqs=kk*Nsite+site2-1
            if (kks.ge.N3d.or.qqs.ge.N3d) then
               if (myid.eq.1) write(*,*) "kpoint out of bounds: ", kks, qqs
            endif
            ! if (myid.eq.1) write(*,*) kks,qqs,site1,site2
            ! SPIN UP
            if(BTEST(initstatup,kks)) then       
               afterstatdn=initstatdn
               afterstatup=IBCLR(initstatup,kks)
               aftersign=((-1)**(sumbeforebit(afterstatup,kks)+sumeverybit(initstatdn)))
               if(.not.(BTEST(afterstatup,qqs))) then
                  afterstatup=IBSET(afterstatup,qqs)
                  aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,qqs)+sumeverybit(afterstatdn)))
                  call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                  if(l.eq.-1) stop
                  midH(l) = midH(l)+aftersign*prefac(kks,site2)*v0(jj)
               endif
            endif
            ! SPIN DOWN, why is the fermion sign different?????   
            if(BTEST(initstatdn,kks)) then
               afterstatdn=IBCLR(initstatdn,kks)
               aftersign=((-1)**(sumbeforebit(initstatdn,kks)))
               if(.not.(BTEST(afterstatdn,qqs))) then
                  afterstatdn=IBSET(afterstatdn,qqs)
                  aftersign=aftersign*((-1)**(sumbeforebit(afterstatdn,qqs)))
                  afterstatup=initstatup
                  call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                  if(l.eq.-1) stop
                  midH(l) = midH(l)+aftersign*prefac(kks,site2)*v0(jj)
                  ! if (myid.eq.1) write(*,*) "Result down: ",jj,aftersign*disp(kks,site2),H_value(H(l)),aftersign
               endif
            endif
         enddo
      enddo   
   enddo
enddo
call MPI_Barrier(comm, ierr)
! allocate(finalH(ksize))
! finalH = 0.0d0
v0 = 0.0d0
do jj=1,ksize
   v0(jj) = midH(jj);
enddo
midH = 0.0d0
call MPI_Allreduce(v0, midH, ksize, MPI_DOUBLE_COMPLEX, &
                     MPI_SUM, comm, ierr)
do jj=1,ksize
   v0(jj) = midH(jj);
enddo
call MPI_Barrier(comm, ierr)
! deallocate(midH,finalH)
deallocate(midH)
end subroutine

subroutine GenDoubleOccupiedState(groundH, midH, listqpt, ksize, start_site, end_site, docc)
use NumOfOrbitalAndElectrons; use BettsCluster
use ModelParas; use ConstantParas;
use MPIParas; use mpi;
implicit none    
!******************************************************************************
!                \---
!                 \   +   +
! Double Occ:     /  c   c    c  c
!                /--- k+q k'-q k' k
!                q,k,k'
!******************************************************************************
INTEGER*8, INTENT(IN) :: ksize, start_site, end_site
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double complex, INTENT(INOUT) :: docc
double complex, DIMENSION(ksize), INTENT(IN) :: groundH
double complex, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8 :: jj,kk,l,kkp,qq,ii,pp,kks,kkps,qqs
integer*8 :: temp,tempi,site1,site2
integer*8 :: initstatup,initstatdn
integer*8 :: afterstatup,afterstatdn,initsign,aftersign
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp1,iktemp2,iktemp
integer*8 :: ktempx,ktempy,ktempx1,ktempx2,ktempy1,ktempy2,ktempx3,ktempy3
integer*8, external :: sumeverybit, sumbeforebit, find_kpt
double precision :: time1, time2, ptx, pty
double complex :: fphase
double complex, allocatable :: finalH(:)


if (start_site.gt.Nsite.or.end_site.gt.Nsite) then
   write(*,*) "Correlation: start/end site too large", start_site, end_site, Nsite
   stop
endif

ktempx1=0
ktempy1=0
ktempx2=0
ktempy2=0
midH = 0.0d0
do jj= localstart, localend
   initstatup=listqpt(jj)/2**N;
   initstatdn=mod(listqpt(jj),2**N);
   do kk=0,Nkpt-1 
      do kkp=0,Nkpt-1
         do qq=0,Nkpt-1
            do site1=start_site,end_site
               kks = kk*Nsite+site1-1
               kkps = kkp*Nsite+site1-1
               qqs = qq*Nsite+site1-1
               ! if (jj.eq.1) then
               !    write(*,*) kks,kkps,qqs
               ! endif
               if(BTEST(initstatup,kks).and.BTEST(initstatdn,kkps)) then
                  ktempx1 = kpoint(kks,1) + kpoint(qqs,1)
                  ktempy1 = kpoint(kks,2) + kpoint(qqs,2)
                  ktempx2 = kpoint(kkps,1) - kpoint(qqs,1) + kunit
                  ktempy2 = kpoint(kkps,2) - kpoint(qqs,2) + kunit
                  iktemp1 = find_kpt(ktempx1,ktempy1,site1-1,Nsite) ! k +q
                  iktemp2 = find_kpt(ktempx2,ktempy2,site1-1,Nsite) ! k'-q
                  ! A phase is introduced when Fourier transform/Should be multiples of 2 pi
                  ktempx3 = kpoint(kks,1) + kpoint(kkps,1) - kpoint(iktemp1,1) - kpoint(iktemp2,1)
                  ktempy3 = kpoint(kks,2) + kpoint(kkps,2) - kpoint(iktemp1,2) - kpoint(iktemp2,2)
                  if (Ndim.ne.2) then
                     pty = (ktempx3*0.5d0+ktempy3)*4.0d0/sqrt(3.0)*pi/kunit*2
                     ptx = ktempx3*2.0d0*pi/kunit*2
                  else
                     ptx = ktempx3*2.0d0*pi/kunit 
                     pty = ktempy3*2.0d0*pi/kunit
                  endif
                  fphase = exp(CMPLX(0,1)*(ptx * delta(site1,1) + pty * delta(site1,2)))
                  ! fphase = 1.0d0
                  if (MOD(iktemp1,Nsite).ne.(site1-1)) stop
                  if (MOD(iktemp2,Nsite).ne.(site1-1)) stop
                  afterstatup=IBCLR(initstatup,kks)
                  aftersign=((-1)**(sumbeforebit(afterstatup,kks)+sumeverybit(initstatdn)))
                  afterstatdn=IBCLR(initstatdn,kkps)
                  aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,kkps))
                  if(.not.BTEST(afterstatup,iktemp1).and..not.BTEST(afterstatdn,iktemp2)) then
                     afterstatdn=IBSET(afterstatdn,iktemp2)
                     aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,iktemp2))
                     afterstatup=IBSET(afterstatup,iktemp1)
                     aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,iktemp1)+sumeverybit(afterstatdn)))
                     call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                     if(l.eq.-1) stop
                     midH(l) = midH(l) + aftersign/dfloat(Nkpt)*groundH(jj)*fphase
                     !!
                     ! if (ABS(aftersign/dfloat(Nkpt)*fphase).ge.0.1.and.jj.eq.20) then
                     !    write(*,*) ktempx1,ktempy1,ktempx2,ktempy2,ktempx3,ktempy3
                     !    write(*,*) "FPHASE: ", site1, fphase, aftersign/dfloat(Nkpt)*fphase, midH(l), groundH(jj)
                     ! endif
                     !!
                  endif
                  iktemp1 = 0
                  iktemp2 = 0
               endif
            enddo
         enddo
      enddo
   enddo
enddo

! midH = conjg(midH)
! docc = docc + DOT_PRODUCT(groundH,midH)
call MPI_Barrier(comm, ierr)
allocate(finalH(ksize))
finalH = 0.0
call MPI_Allreduce(midH, finalH, ksize, MPI_DOUBLE_COMPLEX, &
                     MPI_SUM, comm, ierr)

do jj=1,ksize
   midH(jj) = finalH(jj);
enddo
call MPI_Barrier(comm, ierr)
deallocate(finalH)

do ii = 1,ksize
   docc = docc + CONJG(groundH(ii))*midH(ii)
enddo
if (myid.eq.0) write(*,*) "DOCC: ", docc
docc = docc/N

end subroutine

! Need to add number of sites
subroutine GenSpinChargeCorrelatedState(qq, groundH, midH, listqpt0, listqpt, &
                  ksize0, ksize, ctag, start_site, end_site)
use NumOfOrbitalAndElectrons; use BettsCluster
use ModelParas; use ConstantParas;
use MPIParas; use mpi;
implicit none

!                    MINUS: Spin-Spin, PLUS: Charge-Charge
!                    -----
!                     \       +              +   
!                     /    ( C     C      +- C    C      ) | phi 0 >
!                    -----    k+q,up  k,up   k+q,dn k,dn
!                       k
!                    PRL 112, 156402 (2014)

INTEGER*8, INTENT(IN) :: qq, ksize0, ksize, ctag, start_site, end_site
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double complex, DIMENSION(ksize0), INTENT(IN) :: groundH
double complex, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj, kk, kkp, kkps, kks, site1
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp, kstat, ktempx, ktempy
double precision :: ptx, pty
double complex :: fphase
double complex, allocatable :: finalH(:)

midH(:)=0.0d0

! A prefactor might be needed to account for multsite correlation
do kkp= 0, Nkpt-1
   do site1 = start_site, end_site
      ! write(*,*) 'ktemp: ', ktempx, ktempy
      kkps = kkp*Nsite+site1-1
      kks = ksum(kkps,qq)
      ! code to get kkp
      ! write(*,*) 'kk    = ', kkps, kpoint(kkps,1), kpoint(kkps,2)
      ! write(*,*) 'qq    = ', qq, kpoint(qq,1), kpoint(qq,2)
      ! write(*,*) 'kk+qq = ', kks, kpoint(kks,1), kpoint(kks,2)
      ktempx = kpoint(kkps,1) + kpoint(qq,1) - kpoint(kks,1)
      ktempy = kpoint(kkps,2) + kpoint(qq,2) - kpoint(kks,2)
      if (Ndim.ne.2) then
         pty = (ktempx*0.5d0+ktempy)*4.0d0/sqrt(3.0)*pi/kunit*2
         ptx = ktempx*2.0d0*pi/kunit*2
      else
         ptx = ktempx*2.0d0*pi/kunit 
         pty = ktempy*2.0d0*pi/kunit
      endif
      fphase = exp(CMPLX(0,1)*(ptx * delta(site1,1) + pty * delta(site1,2)))
      ! if (myid.eq.0) write(*,*) 'site, fphase: ', site1, fphase
      ! do jj=1,ksize0
      do jj = localstart,localend
         groundstatup=listqpt0(jj)/2**N
         groundstatdn=mod(listqpt0(jj),2**N)
         !! up spin
         if(BTEST(groundstatup,kkps).eqv..true.) then
            midstatup=IBCLR(groundstatup,kkps)
            midstatdn=groundstatdn
            midsign=(-1)**(sumbeforebit(midstatup,kkps)+sumeverybit(midstatdn))
            if(BTEST(midstatup,kks).eqv..false.) then
               afterstatup=IBSET(midstatup,kks)
               afterstatdn=midstatdn
               aftersign=midsign*(-1)**(sumbeforebit(midstatup,kks)+sumeverybit(midstatdn))
               call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
               if(l.eq.-1) then
                  write(*,*) 'No overlap between groundstate and N-1 state, SUP'
                  write(*,*) 'kk=',kkps,'qq=',qq,'kk+qq=',kks
                  stop
                  continue
               endif
               midH(l) = midH(l) + fphase * aftersign * groundH(jj)
            endif
         endif
         !! down spin
         if(BTEST(groundstatdn,kkps).eqv..true.) then
            midstatdn=IBCLR(groundstatdn,kkps)
            midstatup=groundstatup
            midsign=(-1)**(sumbeforebit(midstatdn,kkps))
            if(BTEST(midstatdn,kks).eqv..false.) then
               afterstatdn=IBSET(midstatdn,kks)
               afterstatup=midstatup
               aftersign=midsign*(-1)**(sumbeforebit(midstatdn,kks))
               call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
               if(l.eq.-1) then
                  write(*,*) 'No overlap between groundstate and N-1 state, SDN'
                  write(*,*) 'kk=',kkps,'qq=',qq,'kk+qq=',kks
                  stop
                  continue
               endif
               if (ctag.eq.1) midH(l) = midH(l) - fphase * aftersign * groundH(jj)
               if (ctag.eq.2) midH(l) = midH(l) + fphase * aftersign * groundH(jj)
            endif
         endif
      enddo
   enddo
enddo

call MPI_Barrier(comm, ierr)
allocate(finalH(ksize))
finalH = 0.0d0
call MPI_Allreduce(midH, finalH, ksize, MPI_DOUBLE_COMPLEX, &
                     MPI_SUM, comm, ierr)
do jj=1,ksize
   midH(jj) = finalH(jj);
enddo
call MPI_Barrier(comm, ierr)
deallocate(finalH)

end subroutine

subroutine Gen2SpinCorrelatedState(qq, groundH, midH, listqpt0, listqpt, &
        ksize0, ksize, ctag)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas
use MPIParas
implicit none

INTEGER*8, INTENT(IN) :: qq, ksize0, ksize, ctag
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double complex, DIMENSION(ksize0), INTENT(IN) :: groundH
double complex, DIMENSION(ksize), INTENT(OUT) :: midH

!  This function parallelize over k points naively
!  ctag determines the Raman vertex that will be calculated
! next***               
                        !            -----          -----
                        !             \              \                      
                        !             /    < phi 0 | /    S S    | phi n >
                        !            -----          -----  k q-k 
!********               !                             k
integer*8, external :: ksubtract, ksum
double complex, external :: lsf
integer*8 :: kk, q2qk, q1k
integer*8 :: ksize1
integer*8, allocatable :: listqpt1(:)
double complex:: prefac
double complex, allocatable :: tempv1(:), tempv2(:), tempvc(:)

allocate(listqpt1(HsizeEst))
allocate(tempv2(ksize))


! Sum over k, S_k,S_q-k
do kk = 0, N-1
   if ((myid-1).eq.MOD(kk,nprocs)) then
      q2qk = ksubtract(qq, kk)
      prefac = lsf(kk,q2qk,ctag)
      ! write(*,*) "id, kk, q2qk, lsf, ctag: ", myid, kk, q2qk, prefac, ctag
      ! Sz component 2, (Sum_q2: c+_q2+q-k_up*c_q2_up - c+_q2+q-k_dn*c_q2_dn)
      listqpt1 = 0
      call GenHsp_kspace(q2qk, listqpt1, ksize1)
      allocate(tempv1(ksize1))
      call GenSpinChargeCorrelatedState(q2qk , groundH, tempv1, listqpt0, listqpt1, ksize0, ksize1, int(1,8), int(1,8), Nsite)
      if (myid.eq.1) write(*,*) 'Sz component 1 done'    
      ! Sz component 1, (Sum_q1: c+_q1+k_up*c_q1_up - c+_q1+k_dn*c_q1_dn)
      call GenSpinChargeCorrelatedState(kk, tempv1, tempv2, listqpt1, listqpt, ksize1, ksize, int(1,8), int(1,8), Nsite)
      midH = midH + tempv2 * real(prefac) / 4
      if (myid.eq.1) write(*,*) 'Sz component 2 done'    
      
      ! Sum_spin: c+_q1_up*c_q1+k_dn*c+_q2_dn*c_q2+q-k_up
      tempv2 = 0
      call GenFourParticle(kk, q2qk, groundH, tempv2, listqpt0, listqpt, ksize0, ksize, int(1,8))
      midH = midH + tempv2 * real(prefac) / 2
      if (myid.eq.1) write(*,*) 'Sz*Sz - Sy*Sy done'    
      deallocate(tempv1)
   endif
enddo

deallocate(listqpt1)
deallocate(tempv2)

end subroutine

subroutine Gen3SpinCorrelatedState(qq, groundH, midH, listqpt0, listqpt, &
        ksize0, ksize, ctag)
use NumOfOrbitalAndElectrons; use BettsCluster; use ConstantParas
use MPIParas
implicit none

INTEGER*8, INTENT(IN) :: qq, ksize0, ksize, ctag
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double complex, DIMENSION(ksize0), INTENT(IN) :: groundH
double complex, DIMENSION(ksize), INTENT(OUT) :: midH

!  This function parallelize over k points naively
!  ctag determines the Raman vertex that will be calculated
! next***               
                        !            -----          -----
                        !             \              \                      
                        !             /    < phi 0 | /    S (S X S      )| phi n >
                        !            -----          -----  k  k'  q-k-k' 
!********               !                             k

! Sum_spin: c+_q1_up*c_q1+k_dn*c+_q2_dn*c_q2+q-k_up
! d^{+}_q2,s d_{q2+qqt2},sbar d^{+}_q1,sbar d_{q1+qqt1},
integer*8, external :: ksubtract, ksum
double complex, external :: lsf_chiral
double complex :: prefac, sum
integer*8 :: kk, kkp, q3qkkp, q2qk, q1kp, ntag
integer*8 :: ksize1, ksize2, i
integer*8, allocatable :: listqpt1(:)
double complex, allocatable :: tempv1(:), tempv2(:)

allocate(listqpt1(HsizeEst))
allocate(tempv2(ksize))

do kk = 0, N-1
   do kkp = 0, N-1
      if ((myid-1).eq.MOD(kk*N+kkp,nprocs)) then
         listqpt1 = 0
         tempv2 = 0.0   
         q2qk = ksubtract(qq,kk) ! q-k
         q3qkkp = ksubtract(q2qk,kkp) ! q-k-k'
         prefac = lsf_chiral(kk,kkp,q3qkkp,ctag)
         write(*,*) "id, k, k', prefac: ", myid, kk, kkp, prefac
         ! k,k',q-k-k'            
         ! First Permutation: S_k,S_k',S_q-k-k'
         call GenHsp_kspace(q3qkkp, listqpt1, ksize1)
         allocate(tempv1(ksize1))            
         call GenSpinChargeCorrelatedState(q3qkkp, groundH, tempv1, listqpt0, listqpt1, ksize0, ksize1, int(1,8), int(1,8), Nsite)
         call GenFourParticle(kk, kkp, tempv1, tempv2, listqpt1, &
                              listqpt, ksize1, ksize, int(2,8))
         midH = midH + tempv2 * prefac
         deallocate(tempv1)
         if (myid.eq.1) write(*,*) "processor 1 finishes perm 1"
         ! Second Permutation: S_q-k-k',S_k,S_k'
         listqpt1 = 0
         tempv2 = 0  
         call GenHsp_kspace(kkp, listqpt1, ksize1)
         allocate(tempv1(ksize1))            
         call GenSpinChargeCorrelatedState(kkp, groundH, tempv1, listqpt0, listqpt1, ksize0, ksize1, int(1,8), int(1,8), Nsite)
         call GenFourParticle(q3qkkp, kk, tempv1, tempv2, listqpt1, &
                              listqpt, ksize1, ksize, int(2,8))
         midH = midH + tempv2 * prefac
         deallocate(tempv1)
         if (myid.eq.1) write(*,*) "processor 1 finishes perm 2"
         ! Third Permutation: S_k',S_q-k-k'S_k
         listqpt1 = 0
         tempv2 = 0
         call GenHsp_kspace(kk, listqpt1, ksize1)
         allocate(tempv1(ksize1))            
         call GenSpinChargeCorrelatedState(kk, groundH, tempv1, listqpt0, listqpt1, ksize0, ksize1, int(1,8), int(1,8), Nsite)
         call GenFourParticle(kkp, q3qkkp, tempv1, tempv2, listqpt1, &
                              listqpt, ksize1, ksize, int(2,8))
         midH = midH + tempv2 * prefac
         deallocate(tempv1)
         if (myid.eq.1) write(*,*) "processor 1 finishes perm 3"
      endif
   enddo
enddo
deallocate(listqpt1, tempv2)
do i = 1,ksize
   sum = sum + conjg(midH(i))*midH(i)
enddo
write(*,*) "id, sum: ", myid,sum
end subroutine


subroutine GenFourParticle(q2qk, kk, groundH, midH, listqpt0, listqpt, ksize0, ksize, ntag)
use NumOfOrbitalAndElectrons; use BettsCluster
implicit none

INTEGER*8, INTENT(IN) :: kk, q2qk, ksize0, ksize, ntag
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double complex, DIMENSION(ksize0), INTENT(IN) :: groundH
double complex, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit, ksubtract, ksum
integer*8 :: jj
integer*8 :: qq1s, qq2s, qq1, qq2
integer*8 :: initstatup, initstatdn
integer*8 :: afterstatup,afterstatdn,aftersign,l


midH(:)=0.0d0
do jj= 1, ksize0  !loop B
   initstatup=listqpt0(jj)/2**N;
   initstatdn=mod(listqpt0(jj),2**N);
   do qq1=0,N-1 ! creation operator 1 c_q1
   do qq2=0,N-1 ! creation operator 2 c_q2
      qq1s = ksubtract(qq1, kk) ! Annihilation operator 1 c_q1-k
      qq2s = ksubtract(qq2, q2qk) ! Annihilation operator 2 c_q2-q+k ==> Transfer Momentum = q
      if(BTEST(initstatdn,qq1s)) then   !dn
         afterstatdn=IBCLR(initstatdn,qq1s)
         aftersign=((-1)**sumbeforebit(initstatdn,qq1s)) 
      if(.not.BTEST(initstatup,qq1)) then  !up
         afterstatup=IBSET(initstatup,qq1)
         aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,qq1)+sumeverybit(afterstatdn)))
      if(BTEST(afterstatup,qq2s)) then   !up
         afterstatup=IBCLR(afterstatup,qq2s)
         aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,qq2s)+sumeverybit(afterstatdn)))
      if(.not.BTEST(afterstatdn,qq2)) then !dn
         afterstatdn=IBSET(afterstatdn,qq2)
         aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,qq2))

         call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
         if(l.eq.-1) then
            write(*,*) '4 particle state 1', POPCNT(afterstatup*(2**N)+afterstatdn)
            write(*,*) 'No overlap between groundstate and N-1 state'
            write(*,*) 'kk & qq-kk =',kpoint(kk,1), kpoint(kk,2),",", kpoint(q2qk,1), kpoint(q2qk,2)
            write(*,*) 'q1 =', kpoint(qq1,1), kpoint(qq1,2)
            write(*,*) 'q2 =', kpoint(qq2,1), kpoint(qq2,2)
            write(*,*) 'q1+k =', kpoint(qq1s,1), kpoint(qq1s,2)
            write(*,*) 'q2+q-k =', kpoint(qq2s,1), kpoint(qq2s,2)
            stop
         endif
         midH(l)=midH(l)+aftersign*groundH(jj)
      endif                 
      endif                 
      endif                 
      endif                 
      if(BTEST(initstatup,qq1s)) then  !up
         afterstatup=IBCLR(initstatup,qq1s)
         aftersign=((-1)**(sumbeforebit(afterstatup,qq1s)+sumeverybit(initstatdn)))
      if(.not.BTEST(initstatdn,qq1)) then   !dn
         afterstatdn=IBSET(initstatdn,qq1)
         aftersign=aftersign*((-1)**sumbeforebit(initstatdn,qq1))
      if(BTEST(afterstatdn,qq2s)) then !dn
         afterstatdn=IBCLR(afterstatdn,qq2s)
         aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,qq2s))
      if(.not.BTEST(afterstatup,qq2)) then   !up
         afterstatup=IBSET(afterstatup,qq2)
         aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,qq2)+sumeverybit(afterstatdn)))

         call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
         if(l.eq.-1) then
            write(*,*) '4 particle state 2', POPCNT(afterstatup*(2**N)+afterstatdn)
            write(*,*) 'No overlap between groundstate and N-1 state'
            write(*,*) 'kk & qq-kk =',kpoint(kk,1), kpoint(kk,2),",", kpoint(q2qk,1), kpoint(q2qk,2)
            write(*,*) 'q1 =', kpoint(qq1,1), kpoint(qq1,2)
            write(*,*) 'q2 =', kpoint(qq2,1), kpoint(qq2,2)
            write(*,*) 'q1+k =', kpoint(qq1s,1), kpoint(qq1s,2)
            write(*,*) 'q1+q-k =', kpoint(qq2s,1), kpoint(qq2s,2)
            stop
         endif
         if (ntag.eq.1) midH(l)=midH(l)+aftersign*groundH(jj)
         if (ntag.eq.2) midH(l)=midH(l)-aftersign*groundH(jj)
      endif                 
      endif                 
      endif                 
      endif        
   enddo
   enddo
enddo !loop B

end subroutine


subroutine GenChargeCorrelatedState(qq, groundH, midH, listqpt0, listqpt, &
        ksize0, ksize)
use NumOfOrbitalAndElectrons; use BettsCluster
implicit none

INTEGER*8, INTENT(IN) :: qq, ksize0, ksize
integer*8, DIMENSION(ksize0), INTENT(IN) :: listqpt0
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
double complex, DIMENSION(ksize0), INTENT(IN) :: groundH
double complex, DIMENSION(ksize), INTENT(OUT) :: midH

integer*8, external :: sumbeforebit, sumeverybit
integer*8 :: jj, kk, kkp
integer*8 :: groundstatup, groundstatdn
integer*8 :: px, py
integer*8 :: afterstatup,afterstatdn,initsign,aftersign,l
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp, kstat, ktempx, ktempy


px = kpoint(Kmap(qq),1)
py = kpoint(Kmap(qq),2)


midH=dcmplx(0.0d0,0.0d0)
do jj=1,ksize0           !!!!!change from groundstate basis to first excited state basis

   groundstatup=listqpt0(jj)/2**N
   groundstatdn=mod(listqpt0(jj),2**N)

! next***               !ii  sum over n
                        !            -----         -----
                        !             \             \     +                 2
                        !             /    < phi0 | /    C  C     | phi n >
                        !            -----         -----  k  k-q
!********               !              n             k
   do kk=0,N-1
      ! code to get kkp
      iktemp=0
      kstat=1
      ktempx = mod(kpoint(kk,1)-px+kunit,kunit)
      ktempy = mod(kpoint(kk,2)-py+kunit,kunit)
      do while (kstat.ne.0)
         if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy)) then
            kstat=0
         else
            iktemp=iktemp+1
         endif
         if(iktemp.eq.N) then
            write(*,*) 'iktemp out of bounds, spin up'
            stop
         endif
      enddo
      kkp=iktemp
      if(BTEST(groundstatup,kkp).eqv..true.) then
         midstatup=IBCLR(groundstatup,kkp)
         midstatdn=groundstatdn
         midsign=(-1)**(sumbeforebit(midstatup,kkp)+sumeverybit(midstatdn))
         if(BTEST(midstatup,kk).eqv..false.) then
            afterstatup=IBSET(midstatup,kk)
            afterstatdn=midstatdn
            aftersign=midsign*(-1)**(sumbeforebit(midstatup,kk)+sumeverybit(midstatdn))
            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) "state: ", afterstatup*(2**N)+afterstatdn
               write(*,*) "bits: ", POPCNT(afterstatup*(2**N)+afterstatdn)
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               write(*,*) '============================================'
               return
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)
         endif
      endif

      if(BTEST(groundstatdn,kkp).eqv..true.) then
         midstatdn=IBCLR(groundstatdn,kkp)
         midstatup=groundstatup
         midsign=(-1)**(sumbeforebit(midstatdn,kkp))

         if(BTEST(midstatdn,kk).eqv..false.) then
            afterstatdn=IBSET(midstatdn,kk)
            afterstatup=midstatup
            aftersign=midsign*(-1)**(sumbeforebit(midstatdn,kk))

            call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
            if(l.eq.-1) then
               write(*,*) 'No overlap between groundstate and N-1 state'
               write(*,*) 'kk=',kk,'kkp=',kkp
               write(*,*) '============================================'
               return
            endif
            midH(l)=midH(l)+aftersign*groundH(jj)
         endif
      endif
   enddo
enddo
end subroutine

integer*8 function find_kpt(kx,ky,kptin,inc)
use ConstantParas; use NumOfOrbitalAndElectrons
use BettsCluster; use MPIParas
use mpi
implicit none

INTEGER*8, INTENT(IN) :: kx, ky, kptin, inc
integer*8 :: kstat, iktemp, ktempx, ktempy

! write(*,*) "kx,ky", kx,ky
ktempx = mod(kx,kunit)
ktempy = mod(ky,kunit)
kstat = 1
iktemp = kptin
! write(*,*) "mod: kx,ky", ktempx,ktempy
do while (kstat.ne.0)
   ! Find momentum based on the unit vector * slanted 
   if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy) .or. &
      (kpoint(iktemp,1).eq.(ktempx+kunit)).and.(kpoint(iktemp,2).eq.ktempy)) then
      kstat=0
   else
      iktemp=iktemp+inc !Increment per site
   endif
   if(iktemp.ge.N) then
      write(*,*) 'find kpt, cannot find kpoint'
      write(*,*) kx, ky, ktempx, ktempy
      write(*,*) kptin, inc, iktemp
      write(*,*) kpoint(iktemp,1), kpoint(iktemp,2)
      call MPI_ABORT(MPI_COMM_WORLD,comm,ierr)
      stop
   endif
enddo
find_kpt = iktemp
! write(*,*) "iktemp", iktemp
end function


subroutine GenHsp_kspace(kfixed, listqpt, ksize)
use ConstantParas; use NumOfOrbitalAndElectrons
use BettsCluster
implicit none

INTEGER*8, INTENT(IN):: kfixed
INTEGER*8, DIMENSION(HsizeEst), INTENT(OUT):: listqpt
INTEGER*8, INTENT(OUT):: ksize

integer*8 :: iiup,iidn, Hsizet
integer*8, allocatable :: listup(:,:),listdn(:,:),listversup(:),listversdn(:)
! integer*8 :: listup(65536,2),listdn(65536,2),listversup(65536),listversdn(65536)
integer*8 :: ktempx,ktempy, px, py, kfixed_final
integer*8 :: jj, ii, i, kstat, iq, snum
integer*8 :: iktemp1,iktemp2,iktemp
integer*8 :: kstat1,kstat2
integer*8 :: temp,tempi
integer*8 :: iup,idn, jup, jdn
integer*8 :: initstatup, initstatdn
integer*8, external :: sumeverybit, factorial, find_kpt
integer*8 :: ktempx1,ktempx2,ktempy1,ktempy2

allocate(listup(2**N,2),listdn(2**N,2),listversup(2**N),listversdn(2**N))

!! Modifying kfixed for multiple sites
if (kfixed.ne.-1) then
   snum = MOD(kfixed,Nsite)
   kfixed_final = kfixed - snum
else 
   kfixed_final = 0
endif 
! write(*,*) "New kfixed: ", kfixed_final

iiup=0;iidn=0
do ii=0,2**N-1
   listup(ii+1,1)=0;listup(ii+1,2)=0;listdn(ii+1,1)=0;listdn(ii+1,2)=0;
   listversup(ii+1)=0;listversdn(ii+1)=0
end do
ktempx=0
ktempy=0
iktemp=0
kstat=1
do iup=0,2**N-1
   if(sumeverybit(iup).eq.nup) then                !find one up state
      do i=0,N-1
         if(BTEST(iup,i)) then
            ktempx = ktempx + kpoint(i,1)
            ktempy = ktempy + kpoint(i,2)
         endif
      enddo
      ! write(*,*) "iup number: ",iup
      ! write(*,'(B64)') iup
      iktemp = find_kpt(ktempx,ktempy,int(0,8),Nsite)
      listup(iiup+1,1)=iup;
      listup(iiup+1,2)=iktemp;
      listversup(iup)=iiup;
      iiup=iiup+1;
      ktempx=0;
      ktempy=0;
      iktemp=0;
      kstat=1;
   endif
enddo
!********************************************************************************
! This loop is superfluous at Sz=0
!********************************************************************************

ktempx=0
ktempy=0
iktemp=0
kstat=1

do idn=0,2**N-1
   if(sumeverybit(idn).eq.ndn) then                !find one down state
      do i=0,N-1
         if(BTEST(idn,i)) then
            ktempx = ktempx + kpoint(i,1)
            ktempy = ktempy + kpoint(i,2)
         endif
      enddo
      iktemp = find_kpt(ktempx,ktempy,int(0,8),Nsite)
      listdn(iidn+1,1)=idn;
      listdn(iidn+1,2)=iktemp;
      listversdn(idn)=iidn;
      iidn=iidn+1;
      ktempx=0;
      ktempy=0;
      iktemp=0;
      kstat=1;
   endif
enddo

Hsizet = (factorial(N)/factorial(nup)/factorial(N-nup))*&
(factorial(N)/factorial(ndn)/factorial(N-ndn))

ktempx=0
ktempy=0
iktemp=0
kstat=1
iq=0
listqpt=0
do jj=0,Hsizet-1         !loop B
   jdn=mod(jj,iidn);
   jup=jj/iidn;
   initstatup=listup(jup+1,1);
   initstatdn=listdn(jdn+1,1);
   ktempx=kpoint(listup(jup+1,2),1)+kpoint(listdn(jdn+1,2),1)
   ktempy=kpoint(listup(jup+1,2),2)+kpoint(listdn(jdn+1,2),2)
   iktemp = find_kpt(ktempx,ktempy,int(0,8),Nsite)
!   if(iktemp.eq.Kmap(kfixed)) then
   if (kfixed.ne.-1) then 
      if(iktemp.eq.kfixed_final) then !! Calculate for non (0,0)
         iq=iq+1
         listqpt(iq) = initstatup*(2**N)+initstatdn
      endif
   else
      iq=iq+1
      listqpt(iq) = initstatup*(2**N)+initstatdn
   endif
   iktemp=0
   kstat=1
enddo !loop B

!write(*,*) 'listup&dn set up OK', myid
ksize = iq
! write(*,*) 'Inside GenHsp kfixed = ', kfixed
! write(*,*) 'ksize = ', ksize
! write(*,*) 'px and py = ', kpoint(kfixed, 1), kpoint(kfixed, 2)
! write(*,*) 'px = ', kpoint(Kmap(kfixed), 1)
! write(*,*) 'py = ', kpoint(Kmap(kfixed), 2)
! write(*,*) 'Kmap = ', Kmap(kfixed)

! do i = 1,iq
!    ! write(*,'(B64)') listqpt(i)
!    write(*,*) listqpt(i)
! enddo
deallocate(listup,listdn,listversup,listversdn)
end subroutine




subroutine GenMatrix_kspace(listqpt, ksize, IndexIt, IndexJt, sparseHt, SprSizet)
use NumOfOrbitalAndElectrons; use ModelParas; use MPIParas; use ConstantParas
use BettsCluster
implicit none

integer*8, INTENT(IN) :: ksize
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
integer*8, INTENT(INOUT) :: SprSizet
integer*8, DIMENSION(SprSizet), INTENT(OUT) :: IndexIt, IndexJt
double complex, DIMENSION(SprSizet), INTENT(OUT) :: sparseHt

! integer*8 :: H_index(1:50000)
integer*8 :: jj,kk,l,kkp,qq,ii,pp,kks,kkps,qqs,Nsite_local
integer*8 :: temp,tempi,site1,site2
integer*8 :: initstatup,initstatdn
integer*8 :: afterstatup,afterstatdn,initsign,aftersign
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp1,iktemp2,iktemp
integer*8 :: kstat1,kstat2
integer*8 :: ktempx,ktempy,ktempx1,ktempx2,ktempy1,ktempy2,ktempx3,ktempy3
integer*8 :: printfreq
integer*8, external :: sumeverybit, sumbeforebit, find_kpt
integer*8, allocatable :: H(:), H_index(:)
! double precision :: H_value(1:50000)
double complex, allocatable :: H_value(:)
double complex :: fphase, dsum
double precision :: time1, time2, ptx, pty
!!!! Check MOD
!!!******************Change the parameters here****************

allocate(H(ksize))
if (N3d.ge.18) then  
   if (myid.eq.0) write(*,*) 'Large Matrix, allocate more array...'
   allocate(H_value(HsizeEst/100),H_index(HsizeEst/100))
else
   allocate(H_value(50000))
   allocate(H_index(50000))
endif
! write(*,*) 'H allocation done, myid = ', myid
! if (myid.eq.0) write(*,*) 'Upperbound of array', Ubound(IndexJt)
! if (myid.eq.0) write(*,*) 'Size of array', Size(IndexJt)

H=0
H_value=0
H_index=0
tempi=0
printfreq = SprSizet/50
if (Ndim.eq.5) then 
   Nsite_local = 1
else 
   Nsite_local = Nsite
endif
call CPU_time(time1)
do jj= localstart, localend  !loop B
   if(mod(jj,printfreq).eq.0) write(*,*) 'myid index_percent',myid,1.0*(jj-localstart)/(nloc)
   temp=0
   initstatup=listqpt(jj)/2**N;
   initstatdn=mod(listqpt(jj),2**N);
   ! if (myid.eq.0) WRITE(*,'(B64)') initstatup, initstatdn
   !******************************************************************************
   ! Site energy Delta * sum c+k ck (Use delta = t' as stand in)
   ! Count how many fermions reside on "oxygen"
   !******************************************************************************
   ! if (Ndim.eq.5) then
   !    dsum = 0
   !    do kk=0,Nkpt-1
   !       kks=kk*2+1
   !       if(BTEST(initstatup,kks)) dsum = dsum + ttprimeratio
   !       if(BTEST(initstatdn,kks)) dsum = dsum + ttprimeratio
   !    enddo
   !    if(H(jj).eq.0) then
   !       temp=temp+1
   !       H(jj)=temp
   !       H_value(temp)=dsum
   !       H_index(temp)=jj
   !    else
   !       H_value(H(jj))=H_value(H(jj))+dsum
   !    endif
   ! endif 
   do kk=0,Nkpt-1
      !**************************************************************************
      ! KE in kspace, e(k)n_k, Momentum doesn't mix, plugin dispersion relation
      ! Twisted average boundary condition for hopping, averaging dispersion
      ! ^ This is performed in Bett's cluster
      !**************************************************************************
      do site1=1,Nsite
         do site2=1,Nsite
            kks=kk*Nsite+site1-1
            qqs=kk*Nsite+site2-1
            ! SPIN UP !
            if(BTEST(initstatup,kks)) then       
               afterstatdn=initstatdn
               afterstatup=IBCLR(initstatup,kks)
               aftersign=((-1)**(sumbeforebit(afterstatup,kks)+sumeverybit(initstatdn)))
               if(.not.(BTEST(afterstatup,qqs))) then
                  afterstatup=IBSET(afterstatup,qqs)
                  aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,qqs)+sumeverybit(afterstatdn)))
                  call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                  if(l.eq.-1) stop
                  if(H(l).eq.0) then
                     temp=temp+1
                     H(l)=temp
                     H_value(temp)=aftersign*disp(kks,site2)
                     H_index(temp)=l
                  else
                     H_value(H(l))=H_value(H(l))+aftersign*disp(kks,site2)
                  endif
                  ! if (myid.eq.0) write(*,*) "Result down: ",jj,l,kks,qqs
                  ! if (myid.eq.0) write(*,*) "Result down: ",aftersign*disp(kks,site2),aftersign
               endif
            endif
            ! SPIN DOWN
            if(BTEST(initstatdn,kks)) then
               afterstatdn=IBCLR(initstatdn,kks)
               aftersign=((-1)**(sumbeforebit(initstatdn,kks)))
               if(.not.(BTEST(afterstatdn,qqs))) then
                  afterstatdn=IBSET(afterstatdn,qqs)
                  aftersign=aftersign*((-1)**(sumbeforebit(afterstatdn,qqs)))
                  afterstatup=initstatup
                  call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                  if(l.eq.-1) stop
                  if(H(l).eq.0) then
                     temp=temp+1
                     H(l)=temp
                     H_value(temp)=aftersign*disp(kks,site2)
                     H_index(temp)=l
                  else
                     H_value(H(l))=H_value(H(l))+aftersign*disp(kks,site2)
                  endif
               endif
            endif
         enddo
      enddo   
      !******************************************************************************
      !               \---
      !                \   +   +
      ! Hubbard U:  U  /  c   c    c  c
      !               /--- k+q k'-q k' k
      !               q,k,k'
      !******************************************************************************
      ktempx1=0
      ktempy1=0
      ktempx2=0
      ktempy2=0
      do kkp=0,Nkpt-1
         do qq=0,Nkpt-1
            do site1=1,Nsite_local
               kks = kk*Nsite+site1-1
               kkps = kkp*Nsite+site1-1
               qqs = qq*Nsite+site1-1
               if(BTEST(initstatup,kks).and.BTEST(initstatdn,kkps)) then
                  ktempx1 = kpoint(kks,1) + kpoint(qqs,1)
                  ktempy1 = kpoint(kks,2) + kpoint(qqs,2)
                  ktempx2 = kpoint(kkps,1) - kpoint(qqs,1) + kunit
                  ktempy2 = kpoint(kkps,2) - kpoint(qqs,2) + kunit
                  iktemp1 = find_kpt(ktempx1,ktempy1,site1-1,Nsite)
                  iktemp2 = find_kpt(ktempx2,ktempy2,site1-1,Nsite)
                  ! A phase is introduced when Fourier transform/Should be multiples of 2 pi
                  ktempx3 = kpoint(kks,1) + kpoint(kkps,1) - kpoint(iktemp1,1) - kpoint(iktemp2,1)
                  ktempy3 = kpoint(kks,2) + kpoint(kkps,2) - kpoint(iktemp1,2) - kpoint(iktemp2,2)
                  if (Ndim.eq.2.or.Ndim.eq.5) then
                     ptx = ktempx3*2.0d0*pi/kunit 
                     pty = ktempy3*2.0d0*pi/kunit
                  else
                     pty = (ktempx3*0.5d0+ktempy3)*4.0d0/sqrt(3.0)*pi/kunit
                     ptx = ktempx3*2.0d0*pi/kunit
                  endif
                  fphase = exp(CMPLX(0,1)*(ptx * delta(site1,1) + pty * delta(site1,2)))
                  ! if (Ndim.eq.0) fphase = 1
                  if (MOD(iktemp1,Nsite).ne.(site1-1)) stop
                  if (MOD(iktemp2,Nsite).ne.(site1-1)) stop
                  afterstatup=IBCLR(initstatup,kks)
                  aftersign=((-1)**(sumbeforebit(afterstatup,kks)+sumeverybit(initstatdn)))
                  afterstatdn=IBCLR(initstatdn,kkps)
                  aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,kkps))
                  if(.not.BTEST(afterstatup,iktemp1).and..not.BTEST(afterstatdn,iktemp2)) then
                     afterstatdn=IBSET(afterstatdn,iktemp2)
                     aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,iktemp2))
                     afterstatup=IBSET(afterstatup,iktemp1)
                     aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,iktemp1)+sumeverybit(afterstatdn)))
                     call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                     if(l.eq.-1) stop
                     if(H(l).eq.0) then
                        temp=temp+1
                        H(l)=temp
                        H_value(temp)=aftersign*fphase*U/dfloat(Nkpt)
                        H_index(temp)=l
                     else
                        H_value(H(l))=H_value(H(l))+aftersign*fphase*U/dfloat(Nkpt)
                     endif
                  endif
                  iktemp1 = 0
                  iktemp2 = 0
               endif
            enddo
         enddo
      enddo
   enddo
   ! stop
   do ii=1,temp
      if(H_value(ii).ne.0) then
         tempi = tempi + 1
         IndexJt(tempi)  = H_index(ii) !still on the old numbering
         IndexIt(tempi)  = jj
         sparseHt(tempi) = H_value(ii)
      endif
      H(H_index(ii))=0
   enddo
   !H_index = 0
   !H_value = 0
   !H = 0
enddo   !loop B

SprSizet = tempi

call CPU_time(time2)
call MPI_Barrier(comm, rc)
if ( myid .eq. 0 .or. myid.eq.nprocs-1) then
   print *,'The time to Generate this Matrix is about ', time2-time1
   write(*,*) 'spr_size in GenMatrix subroutine =',SprSizet
   write(*,*) 'localstart&end',localstart,localend
endif

deallocate(H,H_value,H_index)

end subroutine





subroutine SetModelParameters
use ModelParas; use NumOfOrbitalAndElectrons; use BettsCluster
use ConstantParas; use MPIParas
implicit none

integer*8 :: ii, kk, NdimCluster
double precision :: ptx, pty, kxx, kyy
!The Hamiltonain:
!include 'Coulombint.f90'
!define U_ext(:,:), J_ext(:,:)

t=1.0d0; tt=t*ttprimeratio; ttt=t*tttprimeratio
U=t*Utratio;
! t = 0.0d0
miu = 0.0

if (myid.eq.0) then
   write(*,*) 't, tt, ttt, U'
   write(*,'(2F8.2)') t,tt,ttt,U
endif

rA=0.0d0; rB=0.1365d0; rC=0.5093;
Uc=-2.0d0;
!include 'multiplet.f90'
!define E_site(:), U(:), phase(:)
Ud=rA+4.0d0*rB+3.0d0*rC
cfs=1.8d0;
Ds= 0.25d0;
Dt= 0.10d0;
E_d=0.0d0;
lambda=0.0d0
NdimCluster = Ndim
if (NdimCluster .eq. 1) then 
   NdimCluster = 2
   if (myid.eq.0) write(*,*) "Selecting triangular lattice"
endif

if (NdimCluster .eq. 5) then 
   NdimCluster = 3
   if (myid.eq.0) write(*,*) "Selecting effective 2 band model"
endif

selectcase (NdimCluster)
case(2)
Nsite = 1
if (myid.eq.0) write(*,*) "Selecting square/triangular lattice"
selectcase(N3d)
case(4) ! Simple test case
   kpoint(0,1)=0;kpoint(0,2)=0;
   kpoint(1,1)=2;kpoint(1,2)=1;
   kpoint(2,1)=0;kpoint(2,2)=2;
   kpoint(3,1)=2;kpoint(3,2)=3;
   kunit=4
   Kmap(1) = 0
   Kmap(2) = 2
   Kmapsize = 2
case(8)
!************************************************************
! K-points for the 8-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 |     x     |     x
!                 |      (1,3)|      (3,3)
!                 |           |
!                 x-----------x------
!                 |(0,2)      |(2,2)
!                 |           |
!                 |     x     |     x
!                 |      (1,1)|      (3,1)
!                 |           |
!                 x-----------x------
!                 |(0,0)      |(2,0)
!************************************************************
kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=2;kpoint(1,2)=0;
kpoint(2,1)=0;kpoint(2,2)=2;
kpoint(3,1)=2;kpoint(3,2)=2;
kpoint(4,1)=1;kpoint(4,2)=1;
kpoint(5,1)=3;kpoint(5,2)=1;
kpoint(6,1)=1;kpoint(6,2)=3;
kpoint(7,1)=3;kpoint(7,2)=3;

kunit=4

nQpoints=4
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 3
Kmap(4) = 4
Kmapsize = 4


case(10)
!************************************************************
! K-points for the 10B-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 p-----x-----x-----x-----x-----p
!                 |                    b  |
!                 |                       |
!                 x     x     p     x     x     x
!                 |  b                    |
!                 |                       |
!                 x     x     x     x     p     x
!                 |              b        |
!                 |                       |
!                 x     p     x     x     x     x
!                 |                       |  b
!                 |                       | 
!                 x     x     x     p     x     x
!                 |        b              |
!                 |                       |
!                 p-----x-----x-----x-----x-----x

!************************************************************

kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=3;kpoint(1,2)=1;
kpoint(2,1)=2;kpoint(2,2)=4;
kpoint(3,1)=5;kpoint(3,2)=5;
kpoint(4,1)=6;kpoint(4,2)=2;
kpoint(5,1)=9;kpoint(5,2)=3;
kpoint(6,1)=8;kpoint(6,2)=6;
kpoint(7,1)=1;kpoint(7,2)=7;
kpoint(8,1)=4;kpoint(8,2)=8;
kpoint(9,1)=7;kpoint(9,2)=9;

kunit=10

nQpoints=4
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
Kmap(4) = 3
Kmapsize = 4




case(12)
   selectcase(Nmore)
   case(3)
      if (myid.eq.0) write(*,*) "selected 12C cluster"
      kpoint(0,1)=0;   kpoint(0,2)=0
      kpoint(1,1)=3;   kpoint(1,2)=0
      kpoint(2,1)=1;   kpoint(2,2)=1
      kpoint(3,1)=4;   kpoint(3,2)=1
      kpoint(4,1)=2;   kpoint(4,2)=2
      kpoint(5,1)=5;   kpoint(5,2)=2
      kpoint(6,1)=0;   kpoint(6,2)=3
      kpoint(7,1)=3;   kpoint(7,2)=3
      kpoint(8,1)=1;   kpoint(8,2)=4
      kpoint(9,1)=4;   kpoint(9,2)=4
      kpoint(10,1)=2;   kpoint(10,2)=5
      kpoint(11,1)=5;   kpoint(11,2)=5
      kunit=6
      ! if (Ndim.eq.1) then 
      !!! This is high symmetry path for triangular lattice
      nQpoints = 4
      Kmap(1) = 0
      Kmap(2) = 2
      Kmap(3) = 4
      Kmap(4) = 1
      Kmapsize = 4

   case(5)
      ! 4x3 real space cluster
      kpoint(0,1)=0;   kpoint(0,2)=0
      kpoint(1,1)=3;   kpoint(1,2)=0
      kpoint(2,1)=6;   kpoint(2,2)=0
      kpoint(3,1)=9;   kpoint(3,2)=0
      kpoint(4,1)=0;   kpoint(4,2)=4
      kpoint(5,1)=3;   kpoint(5,2)=4
      kpoint(6,1)=6;   kpoint(6,2)=4
      kpoint(7,1)=9;   kpoint(7,2)=4
      kpoint(8,1)=0;   kpoint(8,2)=8
      kpoint(9,1)=3;   kpoint(9,2)=8
      kpoint(10,1)=6;   kpoint(10,2)=8
      kpoint(11,1)=9;   kpoint(11,2)=8

      kunit=12

      nQpoints=6
      Kmap(1) = 0
      Kmap(2) = 1
      Kmap(3) = 2
      Kmap(4) = 3
      Kmap(5) = 6
      Kmap(6) = 7
      Kmapsize = 6

   case(4)
      !************************************************************
      ! K-points for the 12D-site cluster with (+) integer values
      ! with a factor of (pi/2) removed
      !
      !                 x-----x-----x-----x-----x-----x
      !                 |                 |      
      !                 |                 |      
      !                 x     x     x     x     x     x
      !                 |                 |      
      !                 |                 |      
      !                 B     x     x     x     x     x
      !                 |        _        |      
      !                 |                 |      
      !                 x-----x-----x-----B-----x-----x
      !                 |                 |      
      !                 B                 |      
      !                 x     x     x     x     x     x
      !                 |                 |        
      !                 |        P        |       
      !                 x     x     x     x     x     x
      !                 |                 B      
      !                 |                 |      
      !                 B-----x-----x-----x-----x-----x

      !************************************************************

      kpoint(0,1)=0;   kpoint(0,2)=0
      kpoint(1,1)=9;   kpoint(1,2)=1
      kpoint(2,1)=6;   kpoint(2,2)=2
      kpoint(3,1)=3;   kpoint(3,2)=3
      kpoint(4,1)=0;   kpoint(4,2)=4
      kpoint(5,1)=9;   kpoint(5,2)=5
      kpoint(6,1)=6;   kpoint(6,2)=6
      kpoint(7,1)=3;   kpoint(7,2)=7
      kpoint(8,1)=0;   kpoint(8,2)=8
      kpoint(9,1)=9;   kpoint(9,2)=9
      kpoint(10,1)=6;   kpoint(10,2)=10
      kpoint(11,1)=3;   kpoint(11,2)=11
      kunit=12

      nQpoints=5
      Kmap(1) = 0
      Kmap(2) = 2
      Kmap(3) = 3
      Kmap(4) = 4
      Kmap(5) = 6
      Kmapsize = 5
   endselect
case(16)
   selectcase(Nmore)
   case(2)
!************************************************************
! K-points for the 16B-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 x     x     x     x     x
!                 |                       |
!                 |                       |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!                 |(0,3)(1,3) (2,3) (3,3) |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(0,2)(1,2) (2,2) (3,2) |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(0,1)(1,1) (2,1) (3,1) |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(0,0)(1,0) (2,0) (3,0) |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!************************************************************
if (myid.eq.0) write(*,*) "selected 16B cluster"
kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=1;kpoint(1,2)=0;
kpoint(2,1)=2;kpoint(2,2)=0;
kpoint(3,1)=3;kpoint(3,2)=0;
kpoint(4,1)=0;kpoint(4,2)=1;
kpoint(5,1)=1;kpoint(5,2)=1;
kpoint(6,1)=2;kpoint(6,2)=1;
kpoint(7,1)=3;kpoint(7,2)=1;
kpoint(8,1)=0;kpoint(8,2)=2;
kpoint(9,1)=1;kpoint(9,2)=2;
kpoint(10,1)=2;kpoint(10,2)=2;
kpoint(11,1)=3;kpoint(11,2)=2;
kpoint(12,1)=0;kpoint(12,2)=3;
kpoint(13,1)=1;kpoint(13,2)=3;
kpoint(14,1)=2;kpoint(14,2)=3;
kpoint(15,1)=3;kpoint(15,2)=3;

kunit=4
if (Ndim.eq.2) then
   nQpoints=6
   Kmap(1) = 0
   Kmap(2) = 1
   Kmap(3) = 2
   Kmap(4) = 5
   Kmap(5) = 6
   Kmap(6) = 10
   Kmapsize = 6
else
   nQpoints=4
   Kmap(1) = 0
   Kmap(2) = 1
   Kmap(3) = 2
   Kmap(4) = 5
   Kmapsize = 4
endif 
   case(1)
!************************************************************
! K-points for the 16A-site cluster with (+) integer values
! with a factor of (pi/4) removed
!
!                 |  x     x     x     x     x
!                 |                       |
!                 |                       |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!                 |  (1,6)(3,6)(5,6)(7,6) |
!                 |  x     x     x     x  |  x
!                 |                       |
!                 |(0,4)(2,4)(4,4)(6,4)   |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(1,2)(3,2) (5,2) (7,2) |
!                 |  x     x     x     x  |  x
!                 |                       |
!                 |(0,0)(2,0) (4,0) (6,0) |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!************************************************************
if (myid.eq.0) write(*,*) "selected 16A cluster"
kpoint(0,1)=0;   kpoint(0,2)=0
kpoint(1,1)=4;   kpoint(1,2)=0
kpoint(2,1)=2;   kpoint(2,2)=1
kpoint(3,1)=6;   kpoint(3,2)=1
kpoint(4,1)=0;   kpoint(4,2)=2
kpoint(5,1)=4;   kpoint(5,2)=2
kpoint(6,1)=2;   kpoint(6,2)=3
kpoint(7,1)=6;   kpoint(7,2)=3
kpoint(8,1)=0;   kpoint(8,2)=4
kpoint(9,1)=4;   kpoint(9,2)=4
kpoint(10,1)=2;   kpoint(10,2)=5
kpoint(11,1)=6;   kpoint(11,2)=5
kpoint(12,1)=0;   kpoint(12,2)=6
kpoint(13,1)=4;   kpoint(13,2)=6
kpoint(14,1)=2;   kpoint(14,2)=7
kpoint(15,1)=6;   kpoint(15,2)=7
kunit=8

if (Ndim.eq.2) then
   Kmap(1) = 0
   Kmap(2) = 1
   Kmap(3) = 2
   Kmap(4) = 4
   Kmap(5) = 5
   Kmap(6) = 9
   Kmapsize = 6
else
   Kmap(1) = 0
   Kmap(2) = 4
   Kmap(3) = 8
   Kmap(4) = 6
   Kmap(5) = 2
   Kmapsize = 5
endif

! What is the initeresting Kpoint in 16A cluster
case(3)
!************************************************************
! K-points for the 16C-site cluster with (+) integer values
! with a factor of (pi/8) removed
!
!                 |                               |
!                 x---------------x---------------x
!                 |     x(3,14)         x(11,14)  |
!                 |           x(6,12)         x(14,12)
!                 | x(1,10)         x(9,10)       |
!                 |       x(4,8)          x(12,8) |
!                 |             x(7,6)          x(15,6)
!                 |   x(2,4)          x(10,4)     |
!                 |         x(5,2)          x(13,2)
!                 x(0,0)----------x(8,0)----------x
!                 |                               |
!************************************************************
if (myid.eq.0) write(*,*) "selected 16C cluster"

kpoint(0,1)=0;   kpoint(0,2)=0
kpoint(1,1)=10;   kpoint(1,2)=1
kpoint(2,1)=4;   kpoint(2,2)=2
kpoint(3,1)=14;   kpoint(3,2)=3
kpoint(4,1)=8;   kpoint(4,2)=4
kpoint(5,1)=2;   kpoint(5,2)=5
kpoint(6,1)=12;   kpoint(6,2)=6
kpoint(7,1)=6;   kpoint(7,2)=7
kpoint(8,1)=0;   kpoint(8,2)=8
kpoint(9,1)=10;   kpoint(9,2)=9
kpoint(10,1)=4;   kpoint(10,2)=10
kpoint(11,1)=14;   kpoint(11,2)=11
kpoint(12,1)=8;   kpoint(12,2)=12
kpoint(13,1)=2;   kpoint(13,2)=13
kpoint(14,1)=12;   kpoint(14,2)=14
kpoint(15,1)=6;   kpoint(15,2)=15
kunit=16

nQpoints=7
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
Kmap(4) = 4
Kmap(5) = 5
Kmap(6) = 6
Kmap(7) = 8
Kmapsize = 7

case(4)
!************************************************************
! K-points for the 16D-site cluster with (+) integer values
! with a factor of (pi/8) removed
!
!                 |                               |
!                 x-------x-------x-------x-------x
!                 |      (3,12)  (7,12)  (11,12) (15,12)
!                 |     x       x       x       x |
!                 |    (2,8)   (6,8)   (10,8)  (14,8)
!                 |   x       x       x       x   |
!                 |  (1,4)   (5,4)   (9,4)   (13,4)  
!                 | x       x       x       x     |
!                 |(0,0)   (4,0)   (8,0)   (12,0) |
!                 x-------x-------x-------x-------x
!                 |                               |
!************************************************************
if (myid.eq.0) write(*,*) "selected 16D cluster"
kpoint(0,1)=0;   kpoint(0,2)=0
kpoint(1,1)=12;   kpoint(1,2)=1
kpoint(2,1)=8;   kpoint(2,2)=2
kpoint(3,1)=4;   kpoint(3,2)=3
kpoint(4,1)=0;   kpoint(4,2)=4
kpoint(5,1)=12;   kpoint(5,2)=5
kpoint(6,1)=8;   kpoint(6,2)=6
kpoint(7,1)=4;   kpoint(7,2)=7
kpoint(8,1)=0;   kpoint(8,2)=8
kpoint(9,1)=12;   kpoint(9,2)=9
kpoint(10,1)=8;   kpoint(10,2)=10
kpoint(11,1)=4;   kpoint(11,2)=11
kpoint(12,1)=0;   kpoint(12,2)=12
kpoint(13,1)=12;   kpoint(13,2)=13
kpoint(14,1)=8;   kpoint(14,2)=14
kpoint(15,1)=4;   kpoint(15,2)=15
kunit=16

Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
Kmap(4) = 4
Kmap(5) = 5
Kmap(6) = 9
Kmapsize = 6
   endselect
case(18)
   selectcase(Nmore)
   case(1)
!************************************************************
! K-points for the 18A-site cluster with (+) integer values
! with a factor of (pi/3) removed
!
!                 |    x         x         x
!                 |                       
!                 |                       
!                 x---------x---------x---------x
!                 |    (1,5)     (3,5)     (5,5)|
!                 |    x         x         x    |
!                 |(0,4)     (2,4)     (4,4)    |
!                 x---------x---------x---------x
!                 |    (1,3)     (3,3)     (5,3)|
!                 |    x         x         x    |
!                 |(0,2)     (2,2)     (4,2)    |
!                 x---------x---------x---------x
!                 |    (1,1)     (3,1)     (5,1)|
!                 |    x         x         x    |
!                 |(0,0)     (2,0)     (4,0)    |
!                 x---------x---------x---------x
!                 |                             |
!************************************************************
kpoint(0,1)=0;   kpoint(0,2)=0;
kpoint(1,1)=2;   kpoint(1,2)=0;
kpoint(2,1)=4;   kpoint(2,2)=0;
kpoint(3,1)=1;   kpoint(3,2)=1;
kpoint(4,1)=3;   kpoint(4,2)=1;
kpoint(5,1)=5;   kpoint(5,2)=1;
kpoint(6,1)=0;   kpoint(6,2)=2;
kpoint(7,1)=2;   kpoint(7,2)=2;
kpoint(8,1)=4;   kpoint(8,2)=2;
kpoint(9,1)=1;   kpoint(9,2)=3;
kpoint(10,1)=3;  kpoint(10,2)=3;
kpoint(11,1)=5;  kpoint(11,2)=3;
kpoint(12,1)=0;  kpoint(12,2)=4;
kpoint(13,1)=2;  kpoint(13,2)=4;
kpoint(14,1)=4;  kpoint(14,2)=4;
kpoint(15,1)=1;  kpoint(15,2)=5;
kpoint(16,1)=3;  kpoint(16,2)=5;
kpoint(17,1)=5;  kpoint(17,2)=5;
kunit = 6

! This is a bit off center
if (Ndim.eq.2) then
   nQpoints = 7
   Kmap(1) = 0
   Kmap(2) = 1
   Kmap(3) = 2
   Kmap(4) = 3
   Kmap(5) = 4
   Kmap(6) = 7
   Kmap(7) = 10
   Kmapsize = 7
else
   nQpoints = 5
   Kmap(1) = 0
   Kmap(2) = 15
   Kmap(3) = 13
   Kmap(4) = 10
   Kmap(5) = 14
   Kmapsize = 5
endif

case(2)
!************************************************************
! K-points for the 18B-site cluster with (+) integer values
! with a factor of (pi/9) removed
!                 |                       |
!                 x-----------x-----------x
!                 |   x (3,16)     x (12,16)
!                 |       x (6,14)     x (15,14)
!                 x-(0,12)----x-(9,12)----x
!                 |   x (3,10)     x (12,10)
!                 |       x (6,8)     x (15,8)
!                 x-(0,6)-----x-(9,6)-----x
!                 |   x (3,4)     x (12,4)|
!                 |       x (6,2)     x (15,2)
!                 x-(0,0)-----x-(9,0)-----x
!                 |                       |
!************************************************************

kpoint(0,1)=0;   kpoint(0,2)=0
kpoint(1,1)=11;   kpoint(1,2)=1
kpoint(2,1)=4;   kpoint(2,2)=2
kpoint(3,1)=15;   kpoint(3,2)=3
kpoint(4,1)=8;   kpoint(4,2)=4
kpoint(5,1)=1;   kpoint(5,2)=5
kpoint(6,1)=12;   kpoint(6,2)=6
kpoint(7,1)=5;   kpoint(7,2)=7
kpoint(8,1)=16;   kpoint(8,2)=8
kpoint(9,1)=9;   kpoint(9,2)=9
kpoint(10,1)=2;   kpoint(10,2)=10
kpoint(11,1)=13;   kpoint(11,2)=11
kpoint(12,1)=6;   kpoint(12,2)=12
kpoint(13,1)=17;   kpoint(13,2)=13
kpoint(14,1)=10;   kpoint(14,2)=14
kpoint(15,1)=3;   kpoint(15,2)=15
kpoint(16,1)=14;   kpoint(16,2)=16
kpoint(17,1)=7;   kpoint(17,2)=17
kunit=18

if (Ndim.eq.2) then
   nQpoints = 2
   Kmap(1) = 0
   Kmap(2) = 9
   Kmapsize = 2   
else
   nQpoints = 4
   Kmap(1) = 0
   Kmap(2) = 15
   Kmap(3) = 12
   Kmap(4) = 9
   Kmapsize = 4   
endif

case(3)
!************************************************************
! K-points for the 18C-site cluster with (+) integer values
! with a factor of (pi/9) removed
!                 |                                   |
!                 x-----------------x-----------------x
!                 |       x (4,16)          x (13,16) |
!                 |               x (8,14)          x (17,14)
!                 |     x (3,12)          x (12,12)   |         
!                 |             x (7,10)          x (16,10)
!                 |   x (2,8)           x (11,8)      |
!                 |           x (6,6)           x (15,6)
!                 | x (1,4)           x (10,4)        |
!                 |         x (5,2)           x (14,2)|
!                 x-(0,0)----------x-(9,0)------------x
!                 |                                   |
!************************************************************
kpoint(0,1)=0;   kpoint(0,2)=0
kpoint(1,1)=4;   kpoint(1,2)=1
kpoint(2,1)=8;   kpoint(2,2)=2
kpoint(3,1)=12;   kpoint(3,2)=3
kpoint(4,1)=16;   kpoint(4,2)=4
kpoint(5,1)=2;   kpoint(5,2)=5
kpoint(6,1)=6;   kpoint(6,2)=6
kpoint(7,1)=10;   kpoint(7,2)=7
kpoint(8,1)=14;   kpoint(8,2)=8
kpoint(9,1)=0;   kpoint(9,2)=9
kpoint(10,1)=4;   kpoint(10,2)=10
kpoint(11,1)=8;   kpoint(11,2)=11
kpoint(12,1)=12;   kpoint(12,2)=12
kpoint(13,1)=16;   kpoint(13,2)=13
kpoint(14,1)=2;   kpoint(14,2)=14
kpoint(15,1)=6;   kpoint(15,2)=15
kpoint(16,1)=10;   kpoint(16,2)=16
kpoint(17,1)=14;   kpoint(17,2)=17
kunit=18
if (Ndim.eq.2) then
   nQpoints = 4
   Kmap(1) = 0
   Kmap(2) = 1
   Kmap(3) = 6
   Kmap(4) = 9
   Kmapsize = 4
else
   nQpoints = 4
   Kmap(1) = 0
   Kmap(2) = 15
   Kmap(3) = 6
   Kmap(4) = 9
   Kmapsize = 4
endif

case(4)
! K-points for the 18D-site cluster with (+) integer values
! with a factor of (pi/9) removed
kpoint(0,1)=0;   kpoint(0,2)=0
kpoint(1,1)=3;   kpoint(1,2)=1
kpoint(2,1)=6;   kpoint(2,2)=2
kpoint(3,1)=9;   kpoint(3,2)=3
kpoint(4,1)=12;   kpoint(4,2)=4
kpoint(5,1)=15;   kpoint(5,2)=5
kpoint(6,1)=0;   kpoint(6,2)=6
kpoint(7,1)=3;   kpoint(7,2)=7
kpoint(8,1)=6;   kpoint(8,2)=8
kpoint(9,1)=9;   kpoint(9,2)=9
kpoint(10,1)=12;   kpoint(10,2)=10
kpoint(11,1)=15;   kpoint(11,2)=11
kpoint(12,1)=0;   kpoint(12,2)=12
kpoint(13,1)=3;   kpoint(13,2)=13
kpoint(14,1)=6;   kpoint(14,2)=14
kpoint(15,1)=9;   kpoint(15,2)=15
kpoint(16,1)=12;   kpoint(16,2)=16
kpoint(17,1)=15;   kpoint(17,2)=17
kunit=18
if (Ndim.eq.2) then
   nQpoints = 3
   Kmap(1) = 0
   Kmap(2) = 6
   Kmap(3) = 9
   Kmapsize = 3
else
   nQpoints = 3
   Kmap(1) = 0
   Kmap(2) = 6
   Kmap(3) = 9
   Kmapsize = 3
endif

endselect

case(20)
! KPOINTS for 20 site lattice are likely incorrect
   selectcase(Nmore)
   case(1)
!************************************************************
! K-points for the 20A-site cluster with (+) integer values
! with a factor of (pi/5) removed
!                 |                             |
!                 x--------------x--------------x
!                 |     x(2,9)   |     x(7,9)   |
!                 |           x(4,8)         x (9,8)
!                 |  x(1,7)      |  x(6,7)      |
!                 |        x(3,6)|        x(8,6)| 
!                 x(0,5)---------x(5,5)---------x
!                 |     x(2,4)   |     x(7,4)   |
!                 |           x(4,3)         x (9,3)
!                 |  x(1,2)      |  x(6,2)      |
!                 |        x(3,1)|        x(8,1)| 
!                 x(0,0)---------x(5,0)---------x
!                 |                             |
!************************************************************
kpoint(0,1)=0;   kpoint(0,2)=0;
kpoint(1,1)=5;   kpoint(1,2)=0;
kpoint(2,1)=3;   kpoint(2,2)=1;
kpoint(3,1)=8;   kpoint(3,2)=1;
kpoint(4,1)=1;   kpoint(4,2)=2;
kpoint(5,1)=6;   kpoint(5,2)=2;
kpoint(6,1)=4;   kpoint(6,2)=3;
kpoint(7,1)=9;   kpoint(7,2)=3;
kpoint(8,1)=2;   kpoint(8,2)=4;
kpoint(9,1)=7;   kpoint(9,2)=4;
kpoint(10,1)=0;  kpoint(10,2)=5;
kpoint(11,1)=5;  kpoint(11,2)=5;
kpoint(12,1)=3;  kpoint(12,2)=6;
kpoint(13,1)=8;  kpoint(13,2)=6;
kpoint(14,1)=1;  kpoint(14,2)=7;
kpoint(15,1)=6;  kpoint(15,2)=7;
kpoint(16,1)=4;  kpoint(16,2)=8;
kpoint(17,1)=9;  kpoint(17,2)=8;
kpoint(18,1)=2;  kpoint(18,2)=9;
kpoint(19,1)=7;  kpoint(19,2)=9;
kunit = 10
   endselect
case DEFAULT
   write(*,*) "Specified Square lattice size unavailable, stopping..."
   stop
endselect
case(3)
   Nsite = 2
   if (myid.eq.0) write(*,*) "Selecting Honeycomb Lattice"
   selectcase(N3d)
   case(6)
      if (myid.eq.0) write(*,*) "3x1 Honeycomb Lattice"
      kpoint(0,1)=0;kpoint(0,2)=0;
      kpoint(1,1)=0;kpoint(1,2)=0;
      kpoint(2,1)=1;kpoint(2,2)=0;
      kpoint(3,1)=1;kpoint(3,2)=0;
      kpoint(4,1)=2;kpoint(4,2)=0;
      kpoint(5,1)=2;kpoint(5,2)=0;
      kunit=3

      nQpoints=2
      Kmap(1) = 0
      Kmap(2) = 2
      Kmapsize = 2
   case(8)
      if (myid.eq.0) write(*,*) "2x2 Honeycomb Lattice"
      kpoint(0,1)=0;kpoint(0,2)=0;
      kpoint(1,1)=0;kpoint(1,2)=0;
      kpoint(2,1)=1;kpoint(2,2)=0;
      kpoint(3,1)=1;kpoint(3,2)=0;
      kpoint(4,1)=0;kpoint(4,2)=1;
      kpoint(5,1)=0;kpoint(5,2)=1;
      kpoint(6,1)=1;kpoint(6,2)=1;
      kpoint(7,1)=1;kpoint(7,2)=1;
      kunit=2

      nQpoints=2
      Kmap(1) = 0
      Kmap(2) = 2
      Kmapsize = 2
   case(12) 
      selectcase (Nmore)
      case(1)
         if (myid.eq.0) write(*,*) "12A Honeycomb Lattice"
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=2;kpoint(2,2)=0;
         kpoint(3,1)=2;kpoint(3,2)=0;
         kpoint(4,1)=4;kpoint(4,2)=0;
         kpoint(5,1)=4;kpoint(5,2)=0;
         kpoint(6,1)=1;kpoint(6,2)=3;
         kpoint(7,1)=1;kpoint(7,2)=3;
         kpoint(8,1)=3;kpoint(8,2)=3;
         kpoint(9,1)=3;kpoint(9,2)=3;
         kpoint(10,1)=5;kpoint(10,2)=3;
         kpoint(11,1)=5;kpoint(11,2)=3;
         kunit=6

         nQpoints=2
         Kmap(1) = 0
         Kmap(2) = 8
         Kmapsize = 2
      case(2)
         if (myid.eq.0) write(*,*) "3 x 2 Honeycomb Lattice"
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=2;kpoint(2,2)=0;
         kpoint(3,1)=2;kpoint(3,2)=0;
         kpoint(4,1)=4;kpoint(4,2)=0;
         kpoint(5,1)=4;kpoint(5,2)=0;
         kpoint(6,1)=0;kpoint(6,2)=3;
         kpoint(7,1)=0;kpoint(7,2)=3;
         kpoint(8,1)=2;kpoint(8,2)=3;
         kpoint(9,1)=2;kpoint(9,2)=3;
         kpoint(10,1)=4;kpoint(10,2)=3;
         kpoint(11,1)=4;kpoint(11,2)=3;
         kunit=6
         nQpoints=2
         Kmap(1) = 0
         Kmap(2) = 6
         Kmapsize = 2
      endselect 
   case(16)
      selectcase (Nmore)
      case(1)
         ! 8A Honeycomb Lattice
         if (myid.eq.0) write(*,*) "8A Betts Honeycomb Lattice"
         ! kpoint(0,1)=0;kpoint(0,2)=0;
         ! kpoint(1,1)=0;kpoint(1,2)=0;
         ! kpoint(2,1)=2;kpoint(2,2)=0;
         ! kpoint(3,1)=2;kpoint(3,2)=0;
         ! kpoint(4,1)=0;kpoint(4,2)=2;
         ! kpoint(5,1)=0;kpoint(5,2)=2;
         ! kpoint(6,1)=2;kpoint(6,2)=2;
         ! kpoint(7,1)=2;kpoint(7,2)=2;
         ! kpoint(8,1)=1;kpoint(8,2)=1;
         ! kpoint(9,1)=1;kpoint(9,2)=1;
         ! kpoint(10,1)=3;kpoint(10,2)=1;
         ! kpoint(11,1)=3;kpoint(11,2)=1;
         ! kpoint(12,1)=1;kpoint(12,2)=3;
         ! kpoint(13,1)=1;kpoint(13,2)=3;
         ! kpoint(14,1)=3;kpoint(14,2)=3;
         ! kpoint(15,1)=3;kpoint(15,2)=3;


         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=2;kpoint(2,2)=0;
         kpoint(3,1)=2;kpoint(3,2)=0;
         kpoint(4,1)=1;kpoint(4,2)=1;
         kpoint(5,1)=1;kpoint(5,2)=1;
         kpoint(6,1)=3;kpoint(6,2)=1;
         kpoint(7,1)=3;kpoint(7,2)=1;
         kpoint(8,1)=0;kpoint(8,2)=2;
         kpoint(9,1)=0;kpoint(9,2)=2;
         kpoint(10,1)=2;kpoint(10,2)=2;
         kpoint(11,1)=2;kpoint(11,2)=2;
         kpoint(12,1)=1;kpoint(12,2)=3;
         kpoint(13,1)=1;kpoint(13,2)=3;
         kpoint(14,1)=3;kpoint(14,2)=3;
         kpoint(15,1)=3;kpoint(15,2)=3;

         kunit=4

         nQpoints=4
         Kmap(1) = 0
         Kmap(2) = 2
         Kmap(3) = 4
         Kmap(4) = 10
         Kmapsize = 4
      endselect
   case(18)
      selectcase (Nmore)
      case(1)
         !3x3 Honeycomb Lattice
         if (myid.eq.0) write(*,*) "3 x 3 Honeycomb Lattice"
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=1;kpoint(2,2)=0;
         kpoint(3,1)=1;kpoint(3,2)=0;
         kpoint(4,1)=2;kpoint(4,2)=0;
         kpoint(5,1)=2;kpoint(5,2)=0;
         kpoint(6,1)=0;kpoint(6,2)=1;
         kpoint(7,1)=0;kpoint(7,2)=1;
         kpoint(8,1)=1;kpoint(8,2)=1;
         kpoint(9,1)=1;kpoint(9,2)=1;
         kpoint(10,1)=2;kpoint(10,2)=1;
         kpoint(11,1)=2;kpoint(11,2)=1;
         kpoint(12,1)=0;kpoint(12,2)=2;
         kpoint(13,1)=0;kpoint(13,2)=2;
         kpoint(14,1)=1;kpoint(14,2)=2;
         kpoint(15,1)=1;kpoint(15,2)=2;
         kpoint(16,1)=2;kpoint(16,2)=2;
         kpoint(17,1)=2;kpoint(17,2)=2;

         kunit=3

         Kmap(1) = 0
         Kmap(2) = 8
         Kmapsize = 2
      endselect
   case DEFAULT
      write(*,*) "Specified Triangular lattice size unavailable, stopping..."
      stop
   endselect
case (4)
   Nsite = 3
   if (myid.eq.0) write(*,*) "Selecting Kagome Lattice"
   selectcase(N3d)
   case(3)
      kpoint(0,1)=0;kpoint(0,2)=0;
      kpoint(1,1)=0;kpoint(1,2)=0;
      kpoint(2,1)=0;kpoint(2,2)=0;
      kunit=1
      Kmap(1) = 0
      Kmapsize = 1
   case(6) ! 1x2x3 Kagome Lattice
      selectcase (Nmore)
      case(1) !6x
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=0;kpoint(2,2)=0;
         kpoint(3,1)=1;kpoint(3,2)=0;
         kpoint(4,1)=1;kpoint(4,2)=0;
         kpoint(5,1)=1;kpoint(5,2)=0;
         kunit=2

         nQpoints=2
         Kmap(1) = 0
         Kmap(2) = 3
         Kmapsize = 2
      case(2) !6y
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=0;kpoint(2,2)=0;
         kpoint(3,1)=0;kpoint(3,2)=1;
         kpoint(4,1)=0;kpoint(4,2)=1;
         kpoint(5,1)=0;kpoint(5,2)=1;
         kunit=2

         nQpoints=2
         Kmap(1) = 0
         Kmap(2) = 3
         Kmapsize = 2
      case(3) !6c [4,0],[2,2]
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=0;kpoint(2,2)=0;
         kpoint(3,1)=1;kpoint(3,2)=1;
         kpoint(4,1)=1;kpoint(4,2)=1;
         kpoint(5,1)=1;kpoint(5,2)=1;
         kunit=2

         nQpoints=2
         Kmap(1) = 0
         Kmap(2) = 3
         Kmapsize = 2
      endselect
   case(12) !2x2x3 Kagome Lattice
      selectcase (Nmore)
      case(1)
         if (myid.eq.0) write(*,*) "12A Kagome Lattice"
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=0;kpoint(2,2)=0;
         kpoint(3,1)=2;kpoint(3,2)=1;
         kpoint(4,1)=2;kpoint(4,2)=1;
         kpoint(5,1)=2;kpoint(5,2)=1;
         kpoint(6,1)=0;kpoint(6,2)=2;
         kpoint(7,1)=0;kpoint(7,2)=2;
         kpoint(8,1)=0;kpoint(8,2)=2;
         kpoint(9,1)=2;kpoint(9,2)=3;
         kpoint(10,1)=2;kpoint(10,2)=3;
         kpoint(11,1)=2;kpoint(11,2)=3;

         kunit = 4
         nQpoints= 4
         Kmap(1) = 0
         Kmap(2) = 3
         Kmap(3) = 6
         Kmap(4) = 9
         Kmapsize = 4
      case(2)
         if (myid.eq.0) write(*,*) "2 x 2 Kagome Lattice"
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=0;kpoint(2,2)=0;
         kpoint(3,1)=1;kpoint(3,2)=0;
         kpoint(4,1)=1;kpoint(4,2)=0;
         kpoint(5,1)=1;kpoint(5,2)=0;
         kpoint(6,1)=0;kpoint(6,2)=1;
         kpoint(7,1)=0;kpoint(7,2)=1;
         kpoint(8,1)=0;kpoint(8,2)=1;
         kpoint(9,1)=1;kpoint(9,2)=1;
         kpoint(10,1)=1;kpoint(10,2)=1;
         kpoint(11,1)=1;kpoint(11,2)=1;
         kunit=2
         nQpoints=4
         Kmap(1) = 0
         Kmap(2) = 3
         Kmap(3) = 6
         Kmap(4) = 9
         Kmapsize = 4
      endselect
   case(18) 
      selectcase (Nmore)
      case(1)
         !18A Kagome Lattice
         if (myid.eq.0) write(*,*) "18A Kagome Lattice"
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=0;kpoint(2,2)=0;
         kpoint(3,1)=1;kpoint(3,2)=1;
         kpoint(4,1)=1;kpoint(4,2)=1;
         kpoint(5,1)=1;kpoint(5,2)=1;
         kpoint(6,1)=2;kpoint(6,2)=2;
         kpoint(7,1)=2;kpoint(7,2)=2;
         kpoint(8,1)=2;kpoint(8,2)=2; 
         kpoint(9,1)=3;kpoint(9,2)=3;
         kpoint(10,1)=3;kpoint(10,2)=3;
         kpoint(11,1)=3;kpoint(11,2)=3;
         kpoint(12,1)=4;kpoint(12,2)=4;
         kpoint(13,1)=4;kpoint(13,2)=4;
         kpoint(14,1)=4;kpoint(14,2)=4;
         kpoint(15,1)=5;kpoint(15,2)=5;
         kpoint(16,1)=5;kpoint(16,2)=5;
         kpoint(17,1)=5;kpoint(17,2)=5;

         kunit=6

         Kmap(1) = 0
         Kmap(2) = 9
         Kmap(3) = 12
         Kmapsize = 3
      case(2)
         !2x3x3 Kagome Lattice
         if (myid.eq.0) write(*,*) "18B Kagome Lattice"
         kpoint(0,1)=0;kpoint(0,2)=0;
         kpoint(1,1)=0;kpoint(1,2)=0;
         kpoint(2,1)=0;kpoint(2,2)=0;
         kpoint(3,1)=2;kpoint(3,2)=0;
         kpoint(4,1)=2;kpoint(4,2)=0;
         kpoint(5,1)=2;kpoint(5,2)=0;
         kpoint(6,1)=4;kpoint(6,2)=0;
         kpoint(7,1)=4;kpoint(7,2)=0;
         kpoint(8,1)=4;kpoint(8,2)=0; 
         kpoint(9,1)=1;kpoint(9,2)=3;
         kpoint(10,1)=1;kpoint(10,2)=3;
         kpoint(11,1)=1;kpoint(11,2)=3;
         kpoint(12,1)=3;kpoint(12,2)=3;
         kpoint(13,1)=3;kpoint(13,2)=3;
         kpoint(14,1)=3;kpoint(14,2)=3;
         kpoint(15,1)=5;kpoint(15,2)=3;
         kpoint(16,1)=5;kpoint(16,2)=3;
         kpoint(17,1)=5;kpoint(17,2)=3;

         kunit=6

         Kmap(1) = 0
         Kmap(2) = 12
         Kmapsize = 2
      endselect
   case DEFAULT
      write(*,*) "Specified Kagome lattice size unavailable, stopping..."
      stop
   endselect
case DEFAULT
   write(*,*) "Specified geomtery unavailable, stopping..."
   stop
endselect
if (myid.eq.0) write(*,*) 't and tt and ttt and miu', t, tt, ttt, miu
if (myid.eq.0) write(*,*) 'total averaged twist', twistx*twisty
Nkpt = INT(N3d/Nsite)
if (myid.eq.0) write(*,*) 'NKPT, Nsite: ', Nkpt, Nsite
selectcase(Ndim)
case(1) ! Triangular Lattice
   delta = 0.0d0
   if (myid.eq.0) write(*,*) 'Triangular Lattice Dispersion'
   do kk=0,Nkpt-1
      ! Real Space Unit cell is taken to be (1,0), (-1/2,sqrt(3)/2)
      ! (kx,ky) ==> (sqrt(3)/2*kx,0.5kx+ky)
      ptx = kpoint(kk,1)*2.0d0*pi/kunit + twist1*2.0d0*pi/twistx;
      pty = (kpoint(kk,1)*0.5d0+kpoint(kk,2))*4.0d0/sqrt(3.0)*pi/kunit + twist2*2.0d0*pi/twisty;
      ! Twisted boundary condition, disp shift by 2*pi/n
      disp(kk,1) = -2*t*(cos(ptx)+2*cos(ptx/2)*cos(pty*sqrt(3.0)/2)) &
                  -2*tt*(cos(sqrt(3.0)*pty)+2*cos(3*ptx/2)*cos(pty*sqrt(3.0)/2))
      if (myid.eq.0) write(*,*) "kpoint in: ", kpoint(kk,1),kpoint(kk,2), kunit
      if (myid.eq.0) write(*,*) "kpoint, disp: ", ptx,pty,disp(kk,1)
   enddo
case(2) ! Square Lattice
   delta = 0.0d0
   if (myid.eq.0) write(*,*) 'Square Lattice Dispersion'
   do kk=0,Nkpt-1
      ptx = kpoint(kk,1)*2.0d0*pi/kunit + twist1*2.0d0*pi/twistx;
      pty = kpoint(kk,2)*2.0d0*pi/kunit + twist2*2.0d0*pi/twisty;
      ! Twisted boundary condition, disp shift by 2*pi/n
      disp(kk,1) = -2*t*(cos(ptx)+cos(pty))&
                  -4*tt*(cos(ptx)*cos(pty))&
                  -2*ttt*(cos(2*ptx)+cos(2*pty))&
                  -miu
      den(kk) = (cos(ptx*2.0d0*pi/kunit)-cos(pty*2.0d0*pi/kunit))
      if (myid.eq.0) write(*,*) "kpoint, disp: ", ptx,pty,disp(kk,1)
   enddo
case(3) ! Honeycomb Lattice
   delta = 0.0d0
   delta(2,1) = 0.0
   delta(2,2) = sqrt(3.0)/3.0
   do kk=0,Nkpt-1
      ptx = kpoint(kk*Nsite,1)*2.0d0*pi/kunit + twist1*2.0d0*pi/twistx;
      pty = (kpoint(kk*Nsite,1)*0.5d0+kpoint(kk*Nsite,2))*4.0d0/sqrt(3.0)*pi/kunit + twist2*2.0d0*pi/twisty;
      ! bond length is 1/sqrt(3) of unit cell
      disp(kk*Nsite,2) = -t * (exp(CMPLX(0,-pty/sqrt(3.0)))+exp(CMPLX(0,ptx/2+pty/2/sqrt(3.0))) &
                               + exp(CMPLX(0,-ptx/2+pty/2/sqrt(3.0))))
      disp(kk*Nsite+1,1) = -t * (exp(CMPLX(0,pty/sqrt(3.0)))+exp(CMPLX(0,ptx/2-pty/2/sqrt(3.0))) &
                               + exp(CMPLX(0,-ptx/2-pty/2/sqrt(3.0))))
      if (myid.eq.0) write(*,*) "kpoint, disp: ", ptx,pty,disp(kk*Nsite,1),disp(kk*Nsite,2)
      if (myid.eq.0) write(*,*) "kpoint, disp: ", ptx,pty,disp(kk*Nsite+1,1),disp(kk*Nsite+1,2)
   enddo
case(4) ! Kagome Lattice, 3 points on a triangle
   delta = 0.0d0
   delta(2,1) = -0.25
   delta(2,2) = sqrt(3.0)/4
   delta(3,1) = 0.5
   delta(3,2) = 0
   do kk=0,Nkpt-1
      ptx = kpoint(kk*Nsite,1)*2.0d0*pi/kunit + twist1*2.0d0*pi/twistx;
      pty = (kpoint(kk*Nsite,1)*0.5d0+kpoint(kk*Nsite,2))*4.0d0/sqrt(3.0)*pi/kunit + twist2*2.0d0*pi/twisty;
      ! Twisted boundary condition, disp shift by 2*pi/n
      disp(kk*Nsite,2) = -2*t * (cos(-ptx/4+sqrt(3.0)/4*pty))
      disp(kk*Nsite,3) = -2*t * (cos(ptx/2))
      disp(kk*Nsite+1,1) = -2*t * (cos(-ptx/4+sqrt(3.0)/4*pty))
      disp(kk*Nsite+1,3) = -2*t * (cos(ptx/4+sqrt(3.0)/4*pty))
      disp(kk*Nsite+2,1) = -2*t * (cos(ptx/2))
      disp(kk*Nsite+2,2) = -2*t * (cos(ptx/4+sqrt(3.0)/4*pty))
      if (myid.eq.0) write(*,*) "Ksite: ", kk*Nsite, kk*Nsite+1, kk*Nsite+2
      if (myid.eq.0) write(*,*) "kpoint1, disp: ", ptx,pty,disp(kk*Nsite,1),disp(kk*Nsite,2),disp(kk*Nsite,3)
      if (myid.eq.0) write(*,*) "kpoint2, disp: ", ptx,pty,disp(kk*Nsite+1,1),disp(kk*Nsite+1,2),disp(kk*Nsite+1,3)
      if (myid.eq.0) write(*,*) "kpoint3, disp: ", ptx,pty,disp(kk*Nsite+2,1),disp(kk*Nsite+2,2),disp(kk*Nsite+2,3)
      if (myid.eq.0) write(*,*) "---------------------------------------------------------------------------------"
   enddo
case(5) ! 2 band Hubbard model, PRB 82 064513 (2010)
   delta = 0.0d0
   delta(2,1) = 0.5
   delta(2,2) = 0.5
   do kk=0,Nkpt-1
      ptx = kpoint(kk*Nsite,1)*2.0d0*pi/kunit + twist1*2.0d0*pi/twistx;
      pty = kpoint(kk*Nsite,2)*2.0d0*pi/kunit + twist2*2.0d0*pi/twisty;
      ! bond length is 1/sqrt(3) of unit cell
      disp(kk*Nsite,2) = 2*t*sqrt(sin(ptx/2)*sin(ptx/2)+sin(pty/2)*sin(pty/2))
      disp(kk*Nsite+1,1) = 2*t*sqrt(sin(ptx/2)*sin(ptx/2)+sin(pty/2)*sin(pty/2))
      disp(kk*Nsite,1) = 0
      disp(kk*Nsite+1,2) = ttprimeratio
      if (myid.eq.0) write(*,*) "kpoint, disp: ", ptx,pty,disp(kk*Nsite,1),disp(kk*Nsite,2)
      if (myid.eq.0) write(*,*) "kpoint, disp: ", ptx,pty,disp(kk*Nsite+1,1),disp(kk*Nsite+1,2)
   enddo
endselect
end subroutine

!*****************************************************


subroutine ContFracExpan(Hsize, SprSize, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
use ScanRegion; use ConstantParas
implicit none
INTEGER*8, INTENT(IN) :: Hsize, SprSize
INTEGER*8, DIMENSION(SprSize), INTENT(IN) :: IndexI, IndexJ
DOUBLE PRECISION, INTENT(IN) :: E_0
DOUBLE PRECISION, DIMENSION(Hsize),INTENT(IN) :: H_0
DOUBLE PRECISION, DIMENSION(SprSize), INTENT(IN) :: sparseH
DOUBLE PRECISION, DIMENSION(divX+1), INTENT(OUT) :: specX, specY

integer*8 :: ii,jj,kk
double precision :: mysum, factor
double precision, allocatable :: alpha(:), betha(:)
double precision, allocatable:: phi(:), phil(:),phip(:), phipp(:)
double complex :: z
double complex :: Intensity(divX+1)


   allocate(phi(Hsize), phipp(Hsize))
   allocate(phil(Hsize), phip(Hsize))
   allocate(alpha(niter_CFE),betha(niter_CFE))

   phi(1:Hsize) = H_0(1:Hsize)
   ! Normalize the matrix?
   mysum = 0.0;
   do jj=1, Hsize
        mysum = mysum + phi(jj)*phi(jj)
   enddo
   factor = mysum
   do jj=1, Hsize
        phi(jj) = phi(jj)/sqrt(mysum)
   enddo
   write(*,*) 'sum=',mysum
   !--------------------------------------------------

   do ii = 1, niter_CFE

!       According to different paper, the index of betha 
!       is different 
!       here, betha runs from 2 to niter
!             alpha runs from 1 to niter

      phip=0.0d0
      do kk=1,SprSize
         phip(IndexI(kk))=phip(IndexI(kk))+sparseH(kk)*phi(IndexJ(kk))
      enddo

      if(ii.ne.1) then
         do kk=1, Hsize
           phip(kk) = phip(kk) - betha(ii)*phil(kk)
         enddo
      endif

      alpha(ii)=0.0d0
      do kk=1, Hsize
        alpha(ii) = alpha(ii) + phip(kk)*phi(kk)
      enddo
      write(*,*) 'alpha', alpha(ii)

      do kk=1, Hsize
        phipp(kk) = phip(kk) - alpha(ii)*phi(kk)
      enddo
   if(ii.ne.niter_CFE) then
      betha(ii+1) = 0.0d0
      do kk=1, Hsize
        betha(ii+1) = betha(ii+1) + phipp(kk)**2
      enddo
      betha(ii+1) = sqrt(betha(ii+1))
      write(*,*) 'betha', betha(ii+1)
      do kk=1, Hsize
        phil(kk) = phi(kk)
      enddo
      do kk=1, Hsize
        phi(kk) = phipp(kk) / betha(ii+1)
      enddo
   endif

   enddo
   deallocate(phi, phil, phip, phipp)

do ii=0, divX
   z = CMPLX(dble(ii)/divX*(endX-startX)+startX+E_0,epsilone_CFE)
   Intensity(ii+1)=z-alpha(niter_CFE)
   do jj=1,niter_CFE-1
      Intensity(ii+1)=z-alpha(niter_CFE-jj)-betha(niter_CFE-jj+1)**2/Intensity(ii+1)
   enddo
   specX(ii+1)=dble(ii)/divX*(endX-startX)+startX
   specY(ii+1)=-1*1/pi*AIMAG(factor/Intensity(ii+1))
enddo

deallocate(alpha,betha)

end subroutine


subroutine ChenCmplxContFracExpan(Hsize_3d, SprSize_3d, H_e, groundE, IndexI, IndexJ, sparseH, specX, specY)
use NumOfOrbitalAndElectrons; use ConstantParas; use ScanRegion
implicit none
! Repurpose this to pass in real Hamiltonian 
double precision, INTENT(IN) :: groundE
integer*8, INTENT(IN) :: Hsize_3d, SprSize_3d
double complex, DIMENSION(Hsize_3d), INTENT(IN) :: H_e
integer*8, DIMENSION(SprSize_3d), INTENT(IN) :: IndexI, IndexJ
double precision, DIMENSION(SprSize_3d), INTENT(IN) :: sparseH
double precision, DIMENSION(divX+1), INTENT(OUT) :: specX
double precision, DIMENSION(divX+1), INTENT(OUT) :: specY

integer*8 :: ii,jj,kk, i, j, k
integer*8 :: Itemp1, Itemp2, Itemp3
integer*8 :: spr_size, iter, iter_num
double precision :: rtemp1, rtemp2, rtemp3
double precision :: b_n
double precision :: time1, time2, time3
!!double precision, allocatable :: sparseH(:)
!!for complex Hamiltonian matrix:
double precision :: a_n
double complex :: ctemp1, ctemp2, ctemp3, ctemp4
double complex :: con_ct1, con_ct2, con_ct3
double complex, allocatable :: phi_n1(:), phi_nm(:), phi_n0(:)
double complex, allocatable :: carray_1(:), carray_2(:), carray_3(:)
double precision, allocatable :: work_an(:)
double precision, allocatable ::work_bn(:)
double precision :: Emin, Emax
integer*8 :: energy_pt, Hsizet

double complex :: Gamma

Gamma=dcmplx(0.0d0,epsilone_CFE)
write(*,*) 'starts to in ChenCmplx'
call CPU_time(time1)

spr_size=SprSize_3d;
Hsizet = Hsize_3d
write(*,*) 'Hsizet = ', Hsizet
!!========================================================================================
!!Part I: Construct the tri-diagonal matrix and obtain the continued faction coefficients:
!!========================================================================================

allocate(carray_1(Hsizet)); carray_1=0.0d0;
allocate(phi_n0(Hsizet)); phi_n0=0.0d0;
allocate(phi_n1(Hsizet)); phi_n1=0.0d0;
allocate(phi_nm(Hsizet)); phi_nm=0.0d0;
write(*,*) 'allocation is done'

!!For CFE: using the O|GS> as the initial starting vector.
phi_n0(1:Hsizet) = H_e(1:Hsizet)
iter_num=niter_CFE;
allocate(work_an(iter_num)); work_an=0.0d0;
allocate(work_bn(iter_num)); work_bn=0.0d0;

!!Specify the order of continued fraction expansions.
!!Note: Usually iter_num=50 already gives truncation error < 10-6;
!!      For iter_num > 100 orthogonality may leak and result in "NaN".
!!************************************************

!phi = initial
!do
!  c(:)=c(:)*Matrix
!  ct = c(:)*phi

iter=0; j=0;
do while(j.eq.0)

   call CPU_time(time2)

   iter=iter+1

   !!Constructing H|phi_n0>: later want to do in one step so only 3 vectors are needed.
   carray_1=0.0d0;
   do ii=1, spr_size
      carray_1(IndexI(ii))=  carray_1(IndexI(ii))+sparseH(ii)*phi_n0(IndexJ(ii));
   enddo

   !!Calculate <phi_n0|H|phi_n0>:
   ctemp1=0.0d0;
   rtemp1=0.0d0;
   do ii=1, Hsizet
      ctemp1=ctemp1+conjg(carray_1(ii))*phi_n0(ii)
      rtemp1=rtemp1+conjg(phi_n0(ii))*phi_n0(ii)
   enddo

   a_n=ctemp1/rtemp1;


   if(iter.eq.1) b_n=0.0d0;
   phi_n1=carray_1-a_n*phi_n0-b_n*phi_nm;

   !write(*,*) iter, a_n, b_n, rtemp1
   work_an(iter)=a_n;
   work_bn(iter)=b_n;

   rtemp2=0.0d0;
   do ii=1, Hsizet
      rtemp2=rtemp2+phi_n1(ii)*conjg(phi_n1(ii))
   enddo
   b_n=rtemp2/rtemp1;

   phi_nm=phi_n0
   phi_n0=phi_n1;
   phi_n1=0.0d0;

   if(iter.eq.iter_num) j=3

   call CPU_time(time3)
   do ii=1,Hsizet
      phi_n0(ii)=phi_n0(ii)/sqrt(rtemp2)
      phi_nm(ii)=phi_nm(ii)/sqrt(rtemp2)
   enddo


enddo

!!=======================================================
!!Part II: Calculating the continued fraction expansion:
!!=======================================================

!!Note and Suggestion:
!!Specify Emin to be the ground energy;
!!Emax-Emin would be the energy range of interest (for spectra).

!Emin=   -9.2832;
Emin = startX
Emax = endX
!!How many points in the spectra; 1001 is goodl
!!this in principle does not affect the time of the code.
energy_pt=divX+1;

!!The Lorenztian broadening of the spectra:
!Gamma=(0.0d0, 0.3d0);

carray_1=0.0d0;
rtemp1=0.0d0;

!open(unit=123, file='OGS.txt',status='old')
!do ii=1, Hsizet
!   read(123,*) carray_1(ii)
!   rtemp1=rtemp1+carray_1(ii)*conjg(carray_1(ii))
!enddo
carray_1(1:Hsizet) = H_e(1:Hsizet)
do ii=1, Hsizet
   rtemp1 = rtemp1 + carray_1(ii)*conjg(carray_1(ii))
enddo

!open(unit=1111,file='check_CFE.txt')


!do ii=1, energy_pt
do eneX=0,divX
   !rtemp2= Emin+ (Emax-Emin)*(ii-1)/(energy_pt-1)
   !ctemp1= rtemp2 + Gamma
   rtemp2 = dble(eneX)/dble(divX)*(endX-startX)+startX
   ctemp1 = rtemp2 + groundE + Gamma
   j=0
   do jj=iter_num, 2, -1

      j=j+1
      if(j.eq.1) ctemp4=0.0d0;

      ctemp3=0.0d0;
      ctemp3= work_bn(jj)/(ctemp1-work_an(jj) -ctemp4);
      ctemp4= ctemp3;

   enddo
   con_ct1=rtemp1/(ctemp1-work_an(1)-ctemp4)


   j=0
   do jj=iter_num-5, 2, -1

      j=j+1
      if(j.eq.1) ctemp4=0.0d0;

      ctemp3=0.0d0;
      ctemp3= work_bn(jj)/(ctemp1-work_an(jj) -ctemp4);
      ctemp4= ctemp3;

   enddo
   con_ct2=rtemp1/(ctemp1-work_an(1)-ctemp4)

   write(*,*) rtemp2, imag(con_ct1)/(-pi), imag(con_ct2)/(-pi)
   specX(eneX+1) = rtemp2
   specY(eneX+1) = imag(con_ct1)/(-pi)

enddo
!close(1111)


call CPU_time(time2)
write(*,*) 'Complex CFE (secs):', time2-time1
deallocate(carray_1, phi_n0, phi_n1, phi_nm, work_an, work_bn)

end subroutine

!********************************************************

subroutine BinarySearch(listqpt,ksize,statupdn,l)
implicit none
INTEGER*8, INTENT(IN) :: ksize
INTEGER*8, DIMENSION(1:ksize), INTENT(IN) :: listqpt
INTEGER*8, INTENT(IN) :: statupdn
INTEGER*8, INTENT(OUT):: l

integer*8:: head, tail, middle

head=1; tail=ksize; middle=(head+tail)/2

do while((listqpt(middle).ne.statupdn).and.(head.le.tail)) !! Less than or Less than equal
   if(statupdn.gt.listqpt(middle)) then
      head = middle+1
   else
      tail = middle-1
   endif
   middle = (head+tail)/2
enddo

if(listqpt(middle).eq.statupdn) then
   l=middle;
else
   l=-1;
endif

end

!*****************************************************

integer*8 function factorial(tempn)
implicit none
integer*8:: tempn,tempn1

factorial=1
tempn1=tempn
do while (tempn1.ne.0)
        factorial=factorial*tempn1;
        tempn1=tempn1-1;
enddo
return
end function

!********************************************************

INTEGER*8 function sumeverybit(tempa)
implicit none
integer*8::tempa,tempa1,tempsum

tempsum=0
tempa1=tempa
do while (tempa1.ne.0)
        tempsum=tempsum+mod(tempa1,2)
        tempa1=ISHFT(tempa1,-1)
enddo
sumeverybit=tempsum
return
end function

!********************************************************

integer*8 function sumbeforebit(temps,ks)
implicit none
integer*8::temps,temps1,is,ks

sumbeforebit=0;

if(ks.eq.0) return

temps1=temps
do is=0,ks-1
        if(BTEST(temps1,is).eqv..true.) then
                sumbeforebit=sumbeforebit+1
        endif
enddo
return
end function

subroutine Setnloc(Hsize)
use MPIParas
use mpi
implicit none
! include 'mpif.h'
INTEGER*8, INTENT(IN) :: Hsize
integer  ::   status(MPI_STATUS_SIZE)

if(mod(Hsize,nprocs).eq.0) then
   nloc = Hsize/nprocs
else
   nloc = Hsize/nprocs + 1
endif
localstart = 1+nloc*myid
localend = nloc*(myid+1)
if(myid.eq.nprocs-1) then
   localend=Hsize
   nloc=localend-localstart+1
endif

nloc_array=0;
nloc_array(myid+1)=nloc;
tag=1
if (myid.eq.0) then
   do source=1, nprocs-1
      call MPI_RECV(nloc_array(source+1), 1, MPI_INTEGER, source, tag, comm, status, ierr)
   enddo
else
   call MPI_SEND(nloc_array(myid+1), 1, MPI_INTEGER, 0, tag, comm, ierr)
endif
call MPI_BCAST(nloc_array(1), nprocs, MPI_INTEGER, 0, comm, ierr)

!if(myid.eq.0) then
!   do source = 1, nprocs
!      write(*,*) 'nloc_array', source-1, nloc_array(source)
!   enddo
!endif

end subroutine

subroutine PBiCGS(z, SprSizet, groundH, IndexI, IndexJ, sparseH)
use NumOfOrbitalAndElectrons; use ConstantParas; use ScanRegion
use MPIParas
use mpi
implicit none
! include 'mpif.h'

integer  ::   status(MPI_STATUS_SIZE)
double complex, INTENT(IN) :: z
integer*8, INTENT(IN) :: SprSizet
double complex, DIMENSION(Hsize), INTENT(INOUT) :: groundH
integer*8, DIMENSION(SprSizet), INTENT(IN) :: IndexI, IndexJ
double complex, DIMENSION(SprSizet), INTENT(IN) :: sparseH

double complex, allocatable:: x0(:), xn(:)
double complex, allocatable:: r0(:), rn(:), r00(:)
double complex, allocatable:: p0(:), pn(:)
double complex, allocatable:: v0(:), vn(:)
double complex, allocatable:: s(:), at(:)
double complex, allocatable:: dvtemp(:), phi(:)
double complex, allocatable :: rvtemp(:)
double complex :: alpha, rho0, rhon
double complex :: myrhon, mytsum, tsum, mytsum1, tsum1
double complex :: betha, omega0, omegan
integer*8 :: spr_size, ii, jj
integer*8 :: nloct, lstart, lend, Hsizet
double precision :: time1, time2
     !z = (0.0d0, 0.0d0)
   !****************        STEP 3.A        *****************
   !********              Initialization           **********

spr_size = SprSizet
nloct = nloc
lstart = localstart
lend = localend
Hsizet = Hsize

allocate(x0(nloct),xn(nloct))
allocate(p0(nloct),pn(nloct))
allocate(v0(nloct),vn(nloct))
allocate(r0(nloct),rn(nloct),r00(nloct))
allocate(s(nloct),at(nloct))
allocate(dvtemp(nloct))
allocate(rvtemp(nloct))

allocate(phi(Hsizet))

  x0(:)=0.0d0;
  r0(1:nloct)=groundH(lstart:lend);
  r00(:)=r0(:);
  rho0=1.0d0;
  alpha=1.0d0;
  omega0=1.0d0;
  v0(:)=0.0d0;
  p0(:)=0.0d0;


   !****************        STEP 3.B        *****************
   !*********         Conjugate Gradient     *********************

call CPU_time(time1)

  do ii=1,niter_CG  !CG loop

   ! myrhon = DOT_PRODUCT(CONJG(r00),r0);
   myrhon = (0.0d0,0.0d0)
   do jj=1,nloct
      myrhon = myrhon + CONJG(r00(jj))*r0(jj)
   enddo
   rhon = (0.0d0,0.0d0)
   call MPI_Barrier(comm, ierr)
   call MPI_Allreduce(myrhon, rhon, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
   betha = (rhon/rho0)*(alpha/omega0);
   pn(1:nloc) = r0(1:nloc) + betha*(p0(1:nloc)-omega0*v0(1:nloc))
   call MPI_Barrier(comm, ierr)

   !call MPI_Allgather(pn,nloct,MPI_DOUBLE_COMPLEX,phi,nloct,MPI_DOUBLE_COMPLEX,comm,ierr)
   if(myid.eq.0) then
      phi(1:nloc) = pn(1:nloc);
      do source =1,nprocs-1
         call MPI_RECV(phi( sum(nloc_array(1:source))+1), &
              nloc_array(source+1), MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
      enddo
   else
      call MPI_SEND(pn, nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
   endif
   call MPI_BCAST(phi, Hsize, MPI_DOUBLE_COMPLEX, 0, comm, ierr)

   vn(:)=(0.0d0,0.0d0)
   do jj=1, spr_size
      vn(IndexI(jj)-lstart+1)=vn(IndexI(jj)-lstart+1)-phi(IndexJ(jj))*sparseH(jj)
   enddo
   do jj=1,nloct
      vn(jj)=vn(jj)+pn(jj)*z
   enddo
   ! mytsum = DOT_PRODUCT(conjg(r00),vn)
   mytsum = (0.0d0,0.0d0)
   do jj=1,nloct
      mytsum = mytsum + CONJG(r00(jj))*vn(jj)
   enddo
   call MPI_Barrier(comm, ierr)
   call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
   alpha = rhon/tsum
   s(1:nloct) = r0(1:nloct) - alpha*vn(1:nloct)
   at(:)=(0.0d0,0.0d0)
   call MPI_Barrier(comm, ierr)
   !call MPI_Allgather(s,nloct,MPI_DOUBLE_COMPLEX,phi,nloct,MPI_DOUBLE_COMPLEX,comm,ierr)
   if(myid.eq.0) then
      phi(1:nloc) = s(1:nloc);
      do source =1,nprocs-1
         call MPI_RECV(phi( sum(nloc_array(1:source))+1), &
              nloc_array(source+1), MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
      enddo
   else
      call MPI_SEND(s, nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
   endif
   call MPI_BCAST(phi, Hsize, MPI_DOUBLE_COMPLEX, 0, comm, ierr)

   do jj=1, spr_size
      at(IndexI(jj)-lstart+1)=at(IndexI(jj)-lstart+1)-phi(IndexJ(jj))*sparseH(jj)
   enddo
   do jj=1,nloct
      at(jj)=at(jj)+s(jj)*z
   enddo


   ! mytsum = DOT_PRODUCT(conjg(at),s)
   ! mytsum1 = DOT_PRODUCT(conjg(at),at)
   mytsum = (0.0d0,0.0d0)
   mytsum1 = (0.0d0,0.0d0)
   do jj=1,nloct
      mytsum = mytsum + CONJG(at(jj))*s(jj)
   enddo
   do jj=1,nloct
      mytsum1 = mytsum1 + CONJG(at(jj))*at(jj)
   enddo
   call MPI_Barrier(comm, ierr)
   call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
   call MPI_Allreduce(mytsum1, tsum1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
   omegan = tsum/tsum1
   xn(:) = x0(:)+alpha*pn(:)+omegan*s(:)
   rn(:) = s(:) - omegan*at(:)


  ! xn=x0+alpha*pn+omegan*s;
  ! dvtemp(1:nloc)= x0(1:nloc)+alpha*vn(1:nloc)+omegan*at(1:nloc)-groundH(1:nloc);
   dvtemp(:)=0.0d0
   call MPI_Barrier(comm, ierr)
   !call MPI_Allgather(xn(1),nloct,MPI_DOUBLE_COMPLEX,phi,nloct,MPI_DOUBLE_COMPLEX,comm,ierr)
   if(myid.eq.0) then
      phi(1:nloc) = xn(1:nloc);
      do source =1,nprocs-1
         call MPI_RECV(phi( sum(nloc_array(1:source))+1), &
              nloc_array(source+1), MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
      enddo
   else
      call MPI_SEND(xn, nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
   endif
   call MPI_BCAST(phi, Hsize, MPI_DOUBLE_COMPLEX, 0, comm, ierr)

   do jj=1, spr_size
      dvtemp(IndexI(jj)-lstart+1)=dvtemp(IndexI(jj)-lstart+1)-phi(IndexJ(jj))*sparseH(jj)
   enddo
   do jj=1,nloct
      dvtemp(jj)=dvtemp(jj)+xn(jj)*z
   enddo
   do jj=1,nloct
      dvtemp(jj)=dvtemp(jj)-groundH(jj+lstart-1)
   enddo

   mytsum = (0.0d0,0.0d0)
   do jj=1,nloct
      mytsum = mytsum + CONJG(dvtemp(jj))*dvtemp(jj)
   enddo
   call MPI_Barrier(comm, ierr)
   call MPI_Allreduce(mytsum, tsum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
   if(myid.eq.0) write(*,*) ii, tsum
   if(real(tsum).lt.tol_CG) then
      exit
   endif

   x0(:)=xn(:)
   rho0=rhon;
   p0(:)=pn(:)
   r0(:)=rn(:)
   omega0=omegan;
   v0(:)=vn(:)

  enddo  !CG loop

call CPU_time(time2)
 if(myid.eq.0) write(*,*) 'CG finishes at:', time2, 'secs for ', ii,' iters'
 if(ii.gt.niter_CG) write(*,*) 'BLOWUP TOO BAD', real(tsum)

rvtemp(:)=xn(:)
!call MPI_Allgather(rvtemp,nloct,MPI_DOUBLE_PRECISION,groundH,nloct,MPI_DOUBLE_PRECISION,comm,ierr)
   if(myid.eq.0) then
      groundH(1:nloc) = rvtemp(1:nloc);
      do source =1,nprocs-1
         call MPI_RECV(groundH( sum(nloc_array(1:source))+1), &
              nloc_array(source+1), MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
      enddo
   else
      call MPI_SEND(rvtemp, nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
   endif
   call MPI_BCAST(groundH, Hsize, MPI_DOUBLE_COMPLEX, 0, comm, ierr)

   do jj=1, spr_size
      dvtemp(IndexI(jj)-lstart+1)=dvtemp(IndexI(jj)-lstart+1)-phi(IndexJ(jj))*sparseH(jj)
   enddo

deallocate(x0,xn,p0,pn,v0,vn,r0,rn,r00,s,at,dvtemp,phi)

end subroutine

! znaupd/zneupd needed for Honeycomb Lattice
subroutine ED_CMPLX_PARPACK(Hsize, SprSize, IndexI, IndexJ, sparseH, E_0, H_0, nevin)
use MPIParas; use ConstantParas
use mpi;
implicit none
! include 'mpif.h'

integer  ::   status(MPI_STATUS_SIZE), nevin
INTEGER*8, INTENT(IN) :: Hsize, SprSize
INTEGER*8, DIMENSION(SprSize), INTENT(IN) :: IndexI, IndexJ
double complex, DIMENSION(SprSize), INTENT(IN) :: sparseH
double precision, INTENT(OUT) :: E_0
double complex, DIMENSION(Hsize), INTENT(OUT) :: H_0 !Ground state Eigenvectors

character        bmat*1, which*2
integer       :: ndigit, logfil, msaupd!For MPI use.
integer         :: ldv, ldz, ido, ncv, nev, lworkl, idxsize, iev
integer         :: info, j, nconv, maxitr, mode, ishfts
integer*8       :: flag, tempuse1
double precision :: time1, time2
double precision :: sigma, zero, factor
logical         :: rvec
double complex, allocatable :: temparray_in(:), temparray_gather(:)
double complex, allocatable ::  v(:,:), workl(:), workd(:), d(:), resid(:), z(:,:), workev(:)
double precision, allocatable :: rwork(:)
logical, allocatable :: select(:)
integer, allocatable :: idx(:)
integer         :: iparam(11), ipntr(14)
integer*8 :: ksize, spr_size, ii, jj, acount

comm = MPI_COMM_WORLD
ksize = Hsize
spr_size = SprSize
nev = nevin
ncv = 2*nev+10
allocate(idx(nev))
if(myid.eq.0) write(*,*) 'number of eigenvectors', nev

!**************************PARPACK internal parameters***********************
!     %--------------------------------------------------%
!     | The work array WORKL is used in PSSAUPD as       |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in PSSAUPD to start the Arnoldi        |
!     | iteration.                                       |
!     %--------------------------------------------------%

        zero = 0.0d0;
        bmat = 'I'
        which = 'SR'!Keep the smallest Real Part
        info = 0
        ido = 0

      ndigit = -3
      logfil = 6
      msaupd = 1

!**************************************************************************************        

!Set up the size of the partition of the matrix through nodes:

ido = 0
info = 0
lworkl = 3*ncv**2+5*ncv+8
ldv=nloc
allocate(v(ldv,ncv))
allocate(workd(3*ldv))
allocate(resid(ldv))
ldz=nloc
allocate(z(ldz, nev))
allocate(workev(3*ncv))
allocate(d(nev))
allocate(select(ncv))
allocate(rwork(ncv))
allocate(workl(lworkl))
workl=0;
allocate(temparray_in(ksize));
allocate(temparray_gather(nloc));

!********************************Message Passing Interface*****************************   

!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of PSSAUPD is used    |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in PSSAUPD.                         |
!     %---------------------------------------------------%

      ishfts = 1
      maxitr = 500!Maximum iteration times; change it to 500 for large size problems.
      mode   = 1

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!********************************************************************************************



!Read in the matrix:
acount=0
call MPI_Barrier(comm, ierr)
!********************************************************************************************

!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%

!!!!!!! TODO: kick out this goto syntax
 10   continue

!        %---------------------------------------------%
!        | Repeatedly call the routine PSSAUPD and take| 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
call pznaupd ( comm, ido, bmat, nloc, which, nev, tol_ED, resid, ncv, v, ldv,&
                iparam, ipntr, workd, workl, lworkl, rwork, info )
acount=acount+1
if(myid.eq.0) then
write(*,*) 'the ',acount,'th iteration by pznaupd'
endif

if (ido .eq. -1 .or. ido .eq. 1) then

!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%                    
                        
!*********************************************************************************
tag=1
   
        temparray_in(1:ksize)=0;

        temparray_in(localstart: localend) =workd(ipntr(1):ipntr(1)+nloc-1)
        temparray_gather(1:nloc)=workd(ipntr(1):ipntr(1)+nloc-1)

        tag=1
        if (myid.eq.0) then
           do source=1, nprocs-1
              call MPI_RECV( temparray_in(sum(nloc_array(1:source))+1), nloc_array(source+1), &
                              MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
           enddo
        else
           call MPI_SEND( temparray_gather(1), nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
        endif

        call MPI_BCAST(temparray_in(1), ksize, MPI_DOUBLE_COMPLEX, 0, comm, ierr)

        workd(ipntr(2):ipntr(2)+nloc-1)=0.0d0;
        do ii=1, spr_size
                workd(ipntr(2)+IndexI(ii)-localstart)=workd(ipntr(2)+IndexI(ii)-localstart)+temparray_in(IndexJ(ii))*sparseH(ii)
        enddo


call MPI_Barrier(comm,ierr)
!*********************************************************************************

go to 10 ! L O O P   B A C K to call PZNAUPD again.
end if

!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%

      if ( info .lt. 0 ) then

!        %--------------------------%
!        | Error message. Check the |
!        | documentation in PSSAUPD.|
!        %--------------------------%

         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _saupd, info = ', info
            print *, ' Check documentation in _saupd '
            print *, ' '
         endif

      else

!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using PSSEUPD.               |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%

         rvec = .true.

         if(myid.eq.0) then
            write(*,*) 'ready for pzneupd'
            write(*,*) 'ncv=',ncv
            write(*,*) 'converged Ritz value=',iparam(5)
         endif
         call pzneupd ( comm, rvec, 'A', select, d, z, ldz, sigma, workev, bmat, nloc, which, &
            nev, tol_ED, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info) !d,z eigenval,vec

!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%

         if ( info .ne. 0) then !Error condition
            if ( myid .eq. 0 ) then
                print *, ' '
                print *, ' Error with _neupd, info = ', info
                print *, ' Check the documentation of _neupd. '
                print *, ' nev,ncv,nloc', nev,ncv,nloc
                print *, ' '
            endif                       
         end if

         !*****Print out the eigenvalues and eigenvectors***** 
         ! Lattices used in this program should have only real values!
         if (myid.eq.0) then
            call system("rm eigenvalue.dat")
            open(16,file='eigenvalue.dat',Access = 'append')
            do ii=1, nev
               write(16,*) REAL(d(ii))
               write(*,*) d(ii)
            end do
            close(16)
            call getEigvalindex(nev,REAL(d),idxsize,idx)
         endif
         E_0 = REAL(d(4));
         call MPI_BCAST(idxsize, 1, MPI_INTEGER8, 0, comm, ierr)
         if (myid.eq.0) write(*,*) "index size: ", idxsize
         ! Printing Eigenvector
         open(17,file='write_ev.dat')
         if (idxsize.ne.0) then
            do iev=1,idxsize
               tag=1
               ii = idx(iev)
               call MPI_BCAST(ii, 1, MPI_INTEGER8, 0, comm, ierr)
               if (myid.eq.0) then
                  temparray_in=0;
                  temparray_in(1:nloc)=z(1:nloc,ii)
                  do source=1, nprocs-1
                     call MPI_RECV(temparray_in( sum(nloc_array(1:source))+1), &
                        nloc_array(source+1), MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
                  enddo
               else
                  call MPI_SEND(z(1,ii), nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
               endif
               ! Printing Eigenvector
               if (myid.eq.0) then
                  factor = 0
                  do jj=1, ksize
                     factor = factor + CONJG(temparray_in(jj))*temparray_in(jj)
                  enddo
                  factor = sqrt(factor)
                  do jj=1, ksize
                     temparray_in(jj) = temparray_in(jj)/factor
                     if (idx(iev).eq.1) H_0(jj) = temparray_in(jj)
                  enddo
                  write(*,*) "writing eigenvector at energy, factor: ", d(idx(iev)), factor
                  call write2hdf5cmplx(iev,ksize,temparray_in)
                  write(17,*) iev,idx(iev),REAL(d(idx(iev)))
               endif
               call MPI_BCAST(H_0, Hsize, MPI_DOUBLE_COMPLEX, 0, comm, ierr)
            enddo
         end if
         close(17)
!        ***********************************************************************************************

!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%

         if (myid .eq. 0)then
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
               print *, ' '
            end if

            print *, ' '
            print *, '_NDRV1 '
            print *, '====== '
            print *, ' '
            print *, ' Size of the matrix is ', Hsize
            print *, ' The number of processors is ', nprocs
            print *, ' The number of Ritz values requested is ', nev
            print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
            print *, ' What portion of the spectrum: ', which
!            print *, ' The number of converged Ritz values is ', nconv 
            print *, ' The number of Implicit Arnoldi update',' iterations taken is ', iparam(3)
            print *, ' The number of OP*x is ', iparam(9)
            print *, ' The convergence criterion is ', tol_ED
            print *, ' '
         endif

      end if!end if for info.eq.0

call CPU_time(time2)

if (myid.eq.0) then
write(*,*) 'Total CPU time for ED:', time2-time1, 'seconds.'
end if

deallocate(v,workd,resid,z,workev,d,select,rwork,workl)
end subroutine

subroutine ED_PARPACK(Hsize, SprSize, IndexI, IndexJ, sparseH, E_0, H_0, nevin)
use MPIParas; use ConstantParas
use mpi;
implicit none
! include 'mpif.h'

integer  ::   status(MPI_STATUS_SIZE), nevin
INTEGER*8, INTENT(IN) :: Hsize, SprSize
INTEGER*8, DIMENSION(SprSize), INTENT(IN) :: IndexI, IndexJ
double precision, DIMENSION(SprSize), INTENT(IN) :: sparseH
double precision, INTENT(OUT) :: E_0
double precision, DIMENSION(Hsize), INTENT(OUT) :: H_0 !Ground state Eigenvectors

character        bmat*1, which*2
integer       :: ndigit, logfil, msaupd!For MPI use.
integer         :: ldv, ldz, ido, ncv, nev, lworkl, idxsize, iev
integer         :: info, j, nconv, maxitr, mode, ishfts
integer*8       :: flag, tempuse1
double precision :: time1, time2
double precision :: sigma, zero, factor
double precision, allocatable :: temparray_in(:), temparray_gather(:)
logical         :: rvec
double precision, allocatable ::  v(:,:), workl(:), workd(:), d(:), resid(:), z(:,:)
logical, allocatable :: select(:)
integer, allocatable :: idx(:)
integer         :: iparam(11), ipntr(11)
integer*8 :: ksize, spr_size, ii, jj, acount

comm = MPI_COMM_WORLD
ksize = Hsize
spr_size = SprSize
nev = nevin
ncv = 2*nev+10
allocate(idx(nev))
if(myid.eq.0) write(*,*) 'number of eigenvectors', nev

!**************************PARPACK internal parameters***********************
!     %--------------------------------------------------%
!     | The work array WORKL is used in PSSAUPD as       |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in PSSAUPD to start the Arnoldi        |
!     | iteration.                                       |
!     %--------------------------------------------------%

        zero = 0.0d0;
        bmat = 'I'
        which = 'SA'!Keep the smallest part of the eignvalue spectrum.
        info = 0
        ido = 0

      ndigit = -3
      logfil = 6
      msaupd = 1

!**************************************************************************************        

!Set up the size of the partition of the matrix through nodes:

lworkl = ncv*(ncv+8)
ldv=nloc
allocate(v(ldv,ncv))
allocate(workd(3*ldv))
allocate(resid(ldv))
ldz=nloc
allocate(z(ldz, nev))
allocate(d(nev))
allocate(select(ncv))
allocate(workl(lworkl))
workl=0;
allocate(temparray_in(ksize));
allocate(temparray_gather(nloc));

!********************************Message Passing Interface*****************************   

!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of PSSAUPD is used    |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in PSSAUPD.                         |
!     %---------------------------------------------------%

      ishfts = 1
      maxitr = 600!Maximum iteration times; change it to 500 for large size problems.
      mode   = 1

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!********************************************************************************************



!Read in the matrix:
acount=0
call MPI_Barrier(comm, ierr)
!********************************************************************************************

!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%

!!!!!!! TODO: kick out this goto syntax
 10   continue

!        %---------------------------------------------%
!        | Repeatedly call the routine PSSAUPD and take| 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
call pdsaupd ( comm, ido, bmat, nloc, which, nev, tol_ED, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
acount=acount+1
if(myid.eq.0) then
write(*,*) 'the ',acount,'th iteration by pdsaupd'
endif

if (ido .eq. -1 .or. ido .eq. 1) then

!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%                    
                        
!*********************************************************************************
tag=1

        temparray_in(1:ksize)=0;

        temparray_in(localstart: localend) =workd(ipntr(1):ipntr(1)+nloc-1)
        temparray_gather(1:nloc)=workd(ipntr(1):ipntr(1)+nloc-1)

        tag=1
        if (myid.eq.0) then
           do source=1, nprocs -1
              call MPI_RECV( temparray_in(sum(nloc_array(1:source))+1), nloc_array(source+1), &
                              MPI_DOUBLE_PRECISION, source, tag, comm, status, ierr)
           enddo
        else
           call MPI_SEND( temparray_gather(1), nloc, MPI_DOUBLE_PRECISION, 0, tag, comm, ierr)
        endif

        call MPI_BCAST(temparray_in(1), ksize, MPI_DOUBLE_PRECISION, 0, comm, ierr)

        workd(ipntr(2):ipntr(2)+nloc-1)=0.0d0;
        do ii=1, spr_size
                workd(ipntr(2)+IndexI(ii)-localstart)=workd(ipntr(2)+IndexI(ii)-localstart)+temparray_in(IndexJ(ii))*sparseH(ii)
        enddo


call MPI_Barrier(comm,ierr)
!*********************************************************************************

go to 10 ! L O O P   B A C K to call PSSAUPD again.
end if

!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%

      if ( info .lt. 0 ) then

!        %--------------------------%
!        | Error message. Check the |
!        | documentation in PSSAUPD.|
!        %--------------------------%

         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _saupd, info = ', info
            print *, ' Check documentation in _saupd '
            print *, ' '
         endif

      else

!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using PSSEUPD.               |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%

         rvec = .true.

         if(myid.eq.0) then
         write(*,*) 'ready for pdseupd'
         write(*,*) 'ncv=',ncv
         endif
         call pdseupd ( comm, rvec, 'A', select, d, z, ldz, sigma, bmat, nloc, which, &
            nev, tol_ED, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info) !d,z eigenval,vec

!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | Corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%

         if ( info .ne. 0) then!Error condition: Check the documentation of PSSEUPD.
            if ( myid .eq. 0 ) then
                print *, ' '
                print *, ' Error with _seupd, info = ', info
                print *, ' Check the documentation of _seupd. '
                print *, ' nev,ncv,nloc', nev,ncv,nloc
                print *, ' '
            endif                       
         end if

         !*****Print out the eigenvalues and eigenvectors***** 
         if (myid.eq.0) then
            call system("rm eigenvalue.dat")
            open(16,file='eigenvalue.dat',Access = 'append')
            do ii=1, nev
               write(16,*) d(ii)
               write(*,*) d(ii)
            end do
            close(16)
            call getEigvalindex(nev,d,idxsize,idx)
         endif
         E_0 = d(1);
         call MPI_BCAST(idxsize, 1, MPI_INTEGER8, 0, comm, ierr)
         if (myid.eq.0) write(*,*) "index size: ", idxsize
         ! Printing Eigenvector
         open(17,file='write_ev.dat')
         if (idxsize.ne.0) then
            do iev=1, idxsize
               tag=1
               ii = idx(iev)
               call MPI_BCAST(ii, 1, MPI_INTEGER8, 0, comm, ierr)
               if (myid.eq.0) then
                  temparray_in=0;
                  temparray_in(1:nloc)=z(1:nloc,ii)
                  do source=1, nprocs-1
                     call MPI_RECV(temparray_in( sum(nloc_array(1:source))+1), &
                        nloc_array(source+1), MPI_DOUBLE_PRECISION, source, tag, comm, status, ierr)
                  enddo
               else
                  call MPI_SEND(z(1,ii), nloc, MPI_DOUBLE_PRECISION, 0, tag, comm, ierr)
               endif
               ! Printing Eigenvector
               if (myid.eq.0) then
                  factor = 0
                  do jj=1, ksize
                     factor = factor + temparray_in(jj)*temparray_in(jj)
                  enddo
                  factor = sqrt(factor)
                  do jj=1, ksize
                     temparray_in(jj) = temparray_in(jj)/factor
                     if (idx(iev).eq.1) H_0(jj) = temparray_in(jj)
                  enddo
                  write(*,*) "writing eigenvector at energy, factor: ", d(idx(iev)), factor
                  call write2hdf5(iev,ksize,temparray_in)
                  write(17,*) iev,idx(iev),d(idx(iev))
               endif
               call MPI_BCAST(H_0, Hsize, MPI_DOUBLE_PRECISION, 0, comm, ierr)
            enddo
         end if
         close(17)
!        ***********************************************************************************************

!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%

         if (myid .eq. 0)then
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
               print *, ' '
            end if

            print *, ' '
            print *, '_SDRV1 '
            print *, '====== '
            print *, ' '
            print *, ' Size of the matrix is ', Hsize
            print *, ' The number of processors is ', nprocs
            print *, ' The number of Ritz values requested is ', nev
            print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
            print *, ' What portion of the spectrum: ', which
!            print *, ' The number of converged Ritz values is ', nconv 
            print *, ' The number of Implicit Arnoldi update',' iterations taken is ', iparam(3)
            print *, ' The number of OP*x is ', iparam(9)
            print *, ' The convergence criterion is ', tol_ED
            print *, ' '
         endif

      end if!end if for info.eq.0

call CPU_time(time2)

if (myid.eq.0) then
write(*,*) 'Total CPU time for ED:', time2-time1, 'seconds.'
end if

end subroutine


subroutine readhdf5(idx,hsize,v)
use hdf5
implicit none
integer, INTENT(IN) :: idx
integer*8, INTENT(IN):: hsize
DOUBLE PRECISION, DIMENSION(hsize), INTENT(OUT) :: v

integer*8 :: file_id, dset_id, dspace_id
integer(HSIZE_T), dimension(1) :: dims, maxdims
character(80) :: filename
integer :: i, rank, error

rank = 1
! dims(1) = hsize
if (idx.lt.10) then
   write(filename, '(A,I1,A)') 'eigenvector_',idx,'.h5'
else 
   write(filename, '(A,I2,A)') 'eigenvector_',idx,'.h5'
endif
call h5open_f(error)
call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error)
call h5dopen_f(file_id,"eigvec",dset_id,error)
call h5dget_space_f(dset_id,dspace_id,error)
call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)
if (dims(1).ne.hsize) then
   write(*,*) "Incorrect wavefunction size: ", dims(1), hsize
endif
call h5dread_f(dset_id,h5t_native_double,v,dims,error)
call h5sclose_f(dspace_id,error)
call h5dclose_f(dset_id,error)
call h5fclose_f(file_id,error)
end subroutine

subroutine write2hdf5(idx,hsize,v)
use hdf5
implicit none
integer, INTENT(IN) :: idx
integer*8, INTENT(IN):: hsize
DOUBLE PRECISION, DIMENSION(hsize), INTENT(IN) :: v

integer*8 :: file_id, dset_id, dspace_id
integer(HSIZE_T), dimension(1) :: dims
character(80) :: filename
integer :: i, rank, error

rank = 1
dims(1) = hsize
if (idx.lt.10) then
   write(filename, '(A,I1,A)') 'eigenvector_',idx,'.h5'
else 
   write(filename, '(A,I2,A)') 'eigenvector_',idx,'.h5'
endif
call h5open_f(error)
call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,error)
call h5screate_simple_f(rank,dims,dspace_id,error)
call h5dcreate_f(file_id,"eigvec",h5t_native_double,dspace_id,dset_id,error)
call h5dwrite_f(dset_id,h5t_native_double,v,dims,error)
call h5dclose_f(dset_id,error)
call h5sclose_f(dspace_id,error)
call h5fclose_f(file_id,error)
end subroutine


subroutine readhdf5cmplx(idx,hsize,v)
use hdf5
implicit none
integer, INTENT(IN) :: idx
integer*8, INTENT(IN):: hsize
DOUBLE COMPLEX, DIMENSION(hsize), INTENT(OUT) :: v

integer*8 :: file_id, dset_id, dspace_id
integer(HSIZE_T), dimension(1) :: dims, maxdims
character(80) :: filename
integer :: i, rank, error
double precision, dimension(hsize*2) :: v_read

rank = 1
! dims(1) = hsize
if (idx.lt.10) then
   write(filename, '(A,I1,A)') 'eigenvector_',idx,'.h5'
else 
   write(filename, '(A,I2,A)') 'eigenvector_',idx,'.h5'
endif
call h5open_f(error)
call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error)
call h5dopen_f(file_id,"eigvec",dset_id,error)
call h5dget_space_f(dset_id,dspace_id,error)
call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)
if (dims(1).ne.hsize*2) then
   write(*,*) "Incorrect wavefunction size: ", dims(1), hsize*2
endif
call h5dread_f(dset_id,h5t_native_double,v_read,dims,error)
call h5sclose_f(dspace_id,error)
call h5dclose_f(dset_id,error)
call h5fclose_f(file_id,error)

do i = 1,hsize
   v(i) = CMPLX(v_read(2*i-1),v_read(2*i))
enddo
end subroutine

subroutine write2hdf5cmplx(idx,hsize,v)
use hdf5
implicit none
integer, INTENT(IN) :: idx
integer*8, INTENT(IN):: hsize
DOUBLE COMPLEX, DIMENSION(hsize), INTENT(IN) :: v

integer*8 :: file_id, dset_id, dspace_id
integer(HSIZE_T), dimension(1) :: dims
character(80) :: filename
integer :: i, rank, error
double precision, dimension(hsize*2) :: v_write

!! Break this down to 2 
rank = 1
dims(1) = hsize*2
do i = 1,hsize
   v_write(2*i-1) = REAL(v(i))
   v_write(2*i) = IMAG(v(i))
enddo

if (idx.lt.10) then
   write(filename, '(A,I1,A)') 'eigenvector_',idx,'.h5'
else 
   write(filename, '(A,I2,A)') 'eigenvector_',idx,'.h5'
endif
call h5open_f(error)
call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,error)
call h5screate_simple_f(rank,dims,dspace_id,error)
call h5dcreate_f(file_id,"eigvec",h5t_native_double,dspace_id,dset_id,error)
call h5dwrite_f(dset_id,h5t_native_double,v_write,dims,error)
call h5dclose_f(dset_id,error)
call h5sclose_f(dspace_id,error)
call h5fclose_f(file_id,error)
end subroutine

subroutine getEigvalindex(nev,eigval,idxsize,idx)
use ConstantParas;
integer, INTENT(IN) :: nev
DOUBLE PRECISION, DIMENSION(nev), INTENT(IN) :: eigval
integer, INTENT(OUT) :: idxsize
integer, DIMENSION(nev), INTENT(OUT) :: idx
! Get index of target eigenvalues for writing eigenvectors
integer*8 :: i, j, f
double precision :: e


idxsize = 0
open(unit=10102,file='targetev.dat',Access = 'append',Status='old',iostat=i);
! Check if the file opened successfully
if (i == 0) then
   do
      read(10102,'(1F20.16)',iostat=i) e
      if (i /= 0) exit ! Exit loop on end of file or error
      idxsize = idxsize + 1
      f = 0
      do j = 1, nev
         if (abs(eigval(j)-e) < 1e-5) then
            idx(idxsize) = j
            f = 1
         end if
      enddo
      if (f.eq.0) then
         write(*,*) "target eigenvalue not found", e
         idxsize = idxsize - 1
      else
         write(*,*) idxsize, e, idx(idxsize)
      end if
   end do 
   if (idxsize.eq.0) then
      idxsize = 1
      idx(1) = 1 ! At least return ground state eigenvalue
      write(*,*) idxsize, idx(1), eigval(1)
   end if
   close(10102)
else
   ! Error opening file
   write(*,*) 'targetev.dat dos not exist'
   idxsize = 1
   idx(1) = 1 ! At least return ground state eigenvalue
end if
end subroutine

subroutine PContFracExpan(Hsize, SprSize, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
use ScanRegion; use ConstantParas; use MPIParas
use mpi
implicit none
! include 'mpif.h'

INTEGER*8, INTENT(IN) :: Hsize, SprSize !SprSize is input final states?
INTEGER*8, DIMENSION(SprSize), INTENT(IN) :: IndexI, IndexJ
DOUBLE PRECISION, INTENT(IN) :: E_0
DOUBLE COMPLEX, DIMENSION(Hsize),INTENT(IN) :: H_0
DOUBLE COMPLEX, DIMENSION(SprSize), INTENT(IN) :: sparseH
DOUBLE PRECISION, DIMENSION(divX+1), INTENT(OUT) :: specX, specY

integer  ::   status(MPI_STATUS_SIZE)
integer*8 :: ii,jj,kk
double complex :: mysum, tsum, factor
double complex, allocatable :: alpha(:), betha(:)
double complex, allocatable :: myalpha(:), mybetha(:)
double complex, allocatable:: phi(:), phil(:),phip(:), phipp(:)
double complex :: z
double complex :: Intensity(divX+1)

   comm = MPI_COMM_WORLD
   allocate(phi(Hsize), phipp(nloc))
   allocate(phil(nloc), phip(nloc))
   allocate(alpha(niter_CFE),betha(niter_CFE))
   allocate(myalpha(niter_CFE),mybetha(niter_CFE))

   phi(localstart:localend) = H_0(localstart:localend)

   mysum=0.0d0
   do jj=localstart, localend
        mysum = mysum + conjg(phi(jj))*phi(jj)
   enddo
   call MPI_Allreduce(mysum, tsum, 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
   factor = tsum
   tsum = sqrt(tsum)
   do jj=localstart, localend
        phi(jj) = phi(jj)/tsum
   enddo
   if(myid.eq.0) write(*,*) 'sum=',tsum
   !--------------------------------------------------

   do ii = 1, niter_CFE

!       According to different paper, the index of betha 
!       is different 
!       here, betha runs from 2 to niter
!             alpha runs from 1 to niter

      phip=0.0d0
      !use phip as temprary
      phip(1:nloc)=phi(localstart:localend)
      !call MPI_Gather(phip,nloc,MPI_DOUBLE_PRECISION,phi,nloc,MPI_DOUBLE_PRECISION,0,comm,ierr)  
      tag = 1;
      if(myid.eq.0) then
         phi(1:nloc)=phip(1:nloc)
         do source=1, nprocs-1
            call MPI_RECV(phi( sum(nloc_array(1:source))+1), &
                 nloc_array(source+1), MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
         enddo
      else
         call MPI_SEND(phip(1), nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
      endif
      call MPI_Bcast(phi,Hsize,MPI_DOUBLE_COMPLEX,0,comm,ierr)

      phip=0.0d0
      do kk=1,SprSize
         phip(IndexI(kk)-localstart+1)=phip(IndexI(kk)-localstart+1)+sparseH(kk)*phi(IndexJ(kk))
      enddo

      if(ii.ne.1) then
         do kk=1, nloc
           phip(kk) = phip(kk) - betha(ii)*phil(kk)
         enddo
      endif

      myalpha(ii)=0.0d0
      do kk=1, nloc
        myalpha(ii) = myalpha(ii) + conjg(phi(kk+localstart-1))*phip(kk)
      enddo
      call MPI_Allreduce(myalpha(ii), alpha(ii), 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
      !if(myid.eq.0) write(*,*) 'alpha', alpha(ii)
      do kk=1, nloc
        phipp(kk) = phip(kk) - alpha(ii)*phi(kk+localstart-1)
      enddo

   if(ii.ne.niter_CFE) then
      mybetha(ii+1) = 0.0d0
      do kk=1, nloc
        mybetha(ii+1) = mybetha(ii+1) + conjg(phipp(kk))*phipp(kk)
      enddo
      call MPI_Allreduce(mybetha(ii+1), betha(ii+1), 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
      betha(ii+1) = sqrt(betha(ii+1))
      !if(myid.eq.0) write(*,*) 'betha', betha(ii+1)
      do kk=1, nloc
        phil(kk) = phi(kk+localstart-1)
      enddo
      do kk=1, nloc
        phi(kk+localstart-1) = phipp(kk) / betha(ii+1)
      enddo
   endif

enddo

!***************alpha and betha ready******************
!                  compute I(w) 
 

do ii=0, divX
   z = CMPLX(dble(ii)/divX*(endX-startX)+startX+E_0,epsilone_CFE)
   Intensity(ii+1)=z-alpha(niter_CFE)
   do jj=1,niter_CFE-1
      Intensity(ii+1)=z-alpha(niter_CFE-jj)-betha(niter_CFE-jj+1)**2/Intensity(ii+1)
   enddo
   specX(ii+1)=dble(ii)/divX*(endX-startX)+startX
   specY(ii+1)=-1*1/pi*AIMAG(factor/Intensity(ii+1))
enddo

call MPI_Barrier(comm, ierr)
!write(*,*) 'inside CFE done myid',myid
deallocate(phi, phil, phip, phipp)
deallocate(alpha,betha,myalpha,mybetha)

end subroutine 

subroutine PContFracExpanNgs(Hsize, SprSize, H_0, groundH, E_0, IndexI, IndexJ, sparseH, specX, specY)
use ScanRegion; use ConstantParas; use MPIParas
use mpi
implicit none
! include 'mpif.h'

INTEGER*8, INTENT(IN) :: Hsize, SprSize !SprSize is input final states?
INTEGER*8, DIMENSION(SprSize), INTENT(IN) :: IndexI, IndexJ
DOUBLE PRECISION, INTENT(IN) :: E_0
DOUBLE COMPLEX, DIMENSION(Hsize),INTENT(IN) :: H_0, groundH
DOUBLE COMPLEX, DIMENSION(SprSize), INTENT(IN) :: sparseH
DOUBLE PRECISION, DIMENSION(divX+1), INTENT(OUT) :: specX, specY

integer  ::   status(MPI_STATUS_SIZE)
integer*8 :: ii,jj,kk
double complex :: mysum, tsum, factor, gssum, tgssum
double complex, allocatable :: alpha(:), betha(:)
double complex, allocatable :: myalpha(:), mybetha(:)
double complex, allocatable:: phi(:), phil(:),phip(:), phipp(:)
double complex :: z
double complex :: Intensity(divX+1)

   comm = MPI_COMM_WORLD
   allocate(phi(Hsize), phipp(nloc))
   allocate(phil(nloc), phip(nloc))
   allocate(alpha(niter_CFE),betha(niter_CFE))
   allocate(myalpha(niter_CFE),mybetha(niter_CFE))

   phi(localstart:localend) = H_0(localstart:localend)

   mysum=0.0d0
   do jj=localstart, localend
        mysum = mysum + conjg(phi(jj))*phi(jj)
   enddo
   call MPI_Allreduce(mysum, tsum, 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
   factor = tsum
   tsum = sqrt(tsum)
   do jj=localstart, localend
        phi(jj) = phi(jj)/tsum
   enddo
   if(myid.eq.0) write(*,*) 'sum=',tsum
   ! Calculate subtraction of ground state <H_0|phi_0><phi_0|H_0>
   gssum = 0.0d0
   tgssum = 0.0d0
   do jj = 1, Hsize
      gssum = gssum + conjg(H_0(jj)) * groundH(jj)
   enddo
   ! call MPI_Allreduce(gssum, tgssum, 1, MPI_DOUBLE_PRECISION,MPI_sum,comm,ierr)
   if(myid.eq.0) write(*,*) 'gssum=', gssum, conjg(gssum)*gssum
   tgssum = conjg(gssum) * gssum
   !--------------------------------------------------

   do ii = 1, niter_CFE

!       According to different paper, the index of betha 
!       is different 
!       here, betha runs from 2 to niter
!             alpha runs from 1 to niter

      phip=0.0d0
      !use phip as temprary
      phip(1:nloc)=phi(localstart:localend)
      !call MPI_Gather(phip,nloc,MPI_DOUBLE_PRECISION,phi,nloc,MPI_DOUBLE_PRECISION,0,comm,ierr)  
      tag = 1;
      if(myid.eq.0) then
         phi(1:nloc)=phip(1:nloc)
         do source=1, nprocs-1
            call MPI_RECV(phi( sum(nloc_array(1:source))+1), &
                 nloc_array(source+1), MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
         enddo
      else
         call MPI_SEND(phip(1), nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
      endif
      call MPI_Bcast(phi,Hsize,MPI_DOUBLE_COMPLEX,0,comm,ierr)

      phip=0.0d0
      do kk=1,SprSize
         phip(IndexI(kk)-localstart+1)=phip(IndexI(kk)-localstart+1)+sparseH(kk)*phi(IndexJ(kk))
      enddo

      if(ii.ne.1) then
         do kk=1, nloc
           phip(kk) = phip(kk) - betha(ii)*phil(kk)
         enddo
      endif

      myalpha(ii)=0.0d0
      do kk=1, nloc
        myalpha(ii) = myalpha(ii) + conjg(phi(kk+localstart-1))*phip(kk)
      enddo
      call MPI_Allreduce(myalpha(ii), alpha(ii), 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
      !if(myid.eq.0) write(*,*) 'alpha', alpha(ii)
      do kk=1, nloc
        phipp(kk) = phip(kk) - alpha(ii)*phi(kk+localstart-1)
      enddo

   if(ii.ne.niter_CFE) then
      mybetha(ii+1) = 0.0d0
      do kk=1, nloc
        mybetha(ii+1) = mybetha(ii+1) + conjg(phipp(kk))*phipp(kk)
      enddo
      call MPI_Allreduce(mybetha(ii+1), betha(ii+1), 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
      betha(ii+1) = sqrt(betha(ii+1))
      !if(myid.eq.0) write(*,*) 'betha', betha(ii+1)
      do kk=1, nloc
        phil(kk) = phi(kk+localstart-1)
      enddo
      do kk=1, nloc
        phi(kk+localstart-1) = phipp(kk) / betha(ii+1)
      enddo
   endif

enddo

!***************alpha and betha ready******************
!                  compute I(w) 
 
if (myid.eq.0) then
   ! write(*,*) "factor, gssum: ", factor, tgssum
   do ii=0, divX
      z = CMPLX(dble(ii)/divX*(endX-startX)+startX+E_0,epsilone_CFE)
      Intensity(ii+1)=z-alpha(niter_CFE)
      do jj=1,niter_CFE-1
         Intensity(ii+1)=z-alpha(niter_CFE-jj)-betha(niter_CFE-jj+1)**2/Intensity(ii+1)
      enddo
      ! if (ii.eq.0) write(*,*) "z, int: ", z, Intensity(ii+1)
      specX(ii+1)=dble(ii)/divX*(endX-startX)+startX
      specY(ii+1)=-1*1/pi*AIMAG(factor/Intensity(ii+1) - &
            tgssum/(cmplx(dble(ii)/divX*(endX-startX)+startX,epsilone_CFE)))
      ! write(*,*) "spec, gs, diff: ", factor/Intensity(ii+1), tgssum/(cmplx(dble(ii)/divX*(endX-startX)+startX,epsilone_CFE)) &
      !             , factor/Intensity(ii+1) - tgssum/(cmplx(dble(ii)/divX*(endX-startX)+startX,epsilone_CFE))
   enddo
endif

call MPI_Barrier(comm, ierr)
!write(*,*) 'inside CFE done myid',myid
deallocate(phi, phil, phip, phipp)
deallocate(alpha,betha,myalpha,mybetha)

end subroutine 


subroutine PContFracExpanNRES(Hsize, SprSize, H_0, E_0, IndexI, IndexJ, sparseH, specX, specY)
use ScanRegion; use ConstantParas; use MPIParas
use mpi
implicit none
! include 'mpif.h'

INTEGER*8, INTENT(IN) :: Hsize, SprSize !SprSize is input final states?
INTEGER*8, DIMENSION(SprSize), INTENT(IN) :: IndexI, IndexJ
DOUBLE PRECISION, INTENT(IN) :: E_0
DOUBLE COMPLEX, DIMENSION(Hsize),INTENT(IN) :: H_0
DOUBLE COMPLEX, DIMENSION(SprSize), INTENT(IN) :: sparseH
DOUBLE PRECISION, DIMENSION(divY+1), INTENT(OUT) :: specX, specY

integer  ::   status(MPI_STATUS_SIZE)
integer*8 :: ii,jj,kk
double complex :: mysum, tsum, factor
double complex, allocatable :: alpha(:), betha(:)
double complex, allocatable :: myalpha(:), mybetha(:)
double complex, allocatable:: phi(:), phil(:),phip(:), phipp(:)
double complex :: z
double complex :: Intensity(divY+1)

   allocate(phi(Hsize), phipp(nloc))
   allocate(phil(nloc), phip(nloc))
   allocate(alpha(niter_CFE),betha(niter_CFE))
   allocate(myalpha(niter_CFE),mybetha(niter_CFE))

   phi(localstart:localend) = H_0(localstart:localend)

   mysum=0.0d0
   do jj=localstart, localend
        mysum = mysum + conjg(phi(jj))*phi(jj)
   enddo
   call MPI_Allreduce(mysum, tsum, 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
   factor = tsum
   tsum = sqrt(tsum)
   do jj=localstart, localend
        phi(jj) = phi(jj)/tsum
   enddo
   if(myid.eq.0) write(*,*) 'sum=',tsum
   !--------------------------------------------------

   do ii = 1, niter_CFE

!       According to different paper, the index of betha 
!       is different 
!       here, betha runs from 2 to niter
!             alpha runs from 1 to niter

      phip=0.0d0
      !use phip as temprary
      phip(1:nloc)=phi(localstart:localend)
      !call MPI_Gather(phip,nloc,MPI_DOUBLE_PRECISION,phi,nloc,MPI_DOUBLE_PRECISION,0,comm,ierr)  
      tag = 1;
      if(myid.eq.0) then
         phi(1:nloc)=phip(1:nloc)
         do source=1, nprocs-1
            call MPI_RECV(phi( sum(nloc_array(1:source))+1), &
                 nloc_array(source+1), MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
         enddo
      else
         call MPI_SEND(phip(1), nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
      endif
      call MPI_Bcast(phi,Hsize,MPI_DOUBLE_COMPLEX,0,comm,ierr)

      phip=0.0d0
      do kk=1,SprSize
         phip(IndexI(kk)-localstart+1)=phip(IndexI(kk)-localstart+1)+sparseH(kk)*phi(IndexJ(kk))
      enddo

      if(ii.ne.1) then
         do kk=1, nloc
           phip(kk) = phip(kk) - betha(ii)*phil(kk)
         enddo
      endif

      myalpha(ii)=0.0d0
      do kk=1, nloc
        myalpha(ii) = myalpha(ii) + conjg(phi(kk+localstart-1))*phip(kk)
      enddo
      call MPI_Allreduce(myalpha(ii), alpha(ii), 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
      !if(myid.eq.0) write(*,*) 'alpha', alpha(ii)
      do kk=1, nloc
        phipp(kk) = phip(kk) - alpha(ii)*phi(kk+localstart-1)
      enddo

   if(ii.ne.niter_CFE) then
      mybetha(ii+1) = 0.0d0
      do kk=1, nloc
        mybetha(ii+1) = mybetha(ii+1) + conjg(phipp(kk))*phipp(kk)
      enddo
      call MPI_Allreduce(mybetha(ii+1), betha(ii+1), 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
      betha(ii+1) = sqrt(betha(ii+1))
      !if(myid.eq.0) write(*,*) 'betha', betha(ii+1)
      do kk=1, nloc
        phil(kk) = phi(kk+localstart-1)
      enddo
      do kk=1, nloc
        phi(kk+localstart-1) = phipp(kk) / betha(ii+1)
      enddo
   endif

enddo

!***************alpha and betha ready******************
!                  compute I(w) 

call MPI_Barrier(comm, ierr)
if(myid.eq.0) then 
   do ii=0, divY
      z = CMPLX(dble(ii)/divY*(endY-startY)+startY+E_0,epsilone_CFE)
      Intensity(ii+1)=z-alpha(niter_CFE)
      do jj=1,niter_CFE-1
         Intensity(ii+1)=z-alpha(niter_CFE-jj)-betha(niter_CFE-jj+1)**2/Intensity(ii+1)
      enddo
      specX(ii+1)=dble(ii)/divY*(endY-startY)+startY
      specY(ii+1)=-1*1/pi*AIMAG(factor/Intensity(ii+1))
   enddo
   ! write(*,*) 'inside CFE done myid',myid
endif
call MPI_Barrier(comm, ierr)
deallocate(phi, phil, phip, phipp)
deallocate(alpha,betha,myalpha,mybetha)

end subroutine 

subroutine PContFracExpanNRESNgs(Hsize, SprSize, H_0, groundH, E_0, IndexI, IndexJ, sparseH, specX, specY)
use ScanRegion; use ConstantParas; use MPIParas
use mpi
implicit none
! include 'mpif.h'

INTEGER*8, INTENT(IN) :: Hsize, SprSize !SprSize is input final states?
INTEGER*8, DIMENSION(SprSize), INTENT(IN) :: IndexI, IndexJ
DOUBLE PRECISION, INTENT(IN) :: E_0
DOUBLE COMPLEX, DIMENSION(Hsize),INTENT(IN) :: H_0, groundH
DOUBLE COMPLEX, DIMENSION(SprSize), INTENT(IN) :: sparseH
DOUBLE PRECISION, DIMENSION(divY+1), INTENT(OUT) :: specX, specY

integer  ::   status(MPI_STATUS_SIZE)
integer*8 :: ii,jj,kk
double complex :: mysum, tsum, factor, gssum, tgssum
double complex, allocatable :: alpha(:), betha(:)
double complex, allocatable :: myalpha(:), mybetha(:)
double complex, allocatable:: phi(:), phil(:),phip(:), phipp(:)
double complex :: z
double complex :: Intensity(divY+1)


   allocate(phi(Hsize), phipp(nloc))
   allocate(phil(nloc), phip(nloc))
   allocate(alpha(niter_CFE),betha(niter_CFE))
   allocate(myalpha(niter_CFE),mybetha(niter_CFE))

   phi(localstart:localend) = H_0(localstart:localend)

   mysum=0.0d0
   do jj=localstart, localend
        mysum = mysum + conjg(phi(jj))*phi(jj)
   enddo
   call MPI_Allreduce(mysum, tsum, 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
   factor = tsum
   tsum = sqrt(tsum)
   do jj=localstart, localend
        phi(jj) = phi(jj)/tsum
   enddo
   if(myid.eq.0) write(*,*) 'sum=',tsum
   ! Calculate subtraction of ground state <H_0|phi_0><phi_0|H_0>
   gssum = 0.0d0
   tgssum = 0.0d0
   do jj = 1, Hsize
      gssum = gssum + conjg(H_0(jj))*groundH(jj)
   enddo
   ! call MPI_Allreduce(gssum, tgssum, 1, MPI_DOUBLE_PRECISION,MPI_sum,comm,ierr)
   if(myid.eq.0) write(*,*) 'gssum=', gssum, conjg(gssum)*gssum
   tgssum = conjg(gssum) * gssum
   !--------------------------------------------------

   do ii = 1, niter_CFE

!       According to different paper, the index of betha 
!       is different 
!       here, betha runs from 2 to niter
!             alpha runs from 1 to niter

      phip=0.0d0
      !use phip as temprary
      phip(1:nloc)=phi(localstart:localend)
      !call MPI_Gather(phip,nloc,MPI_DOUBLE_PRECISION,phi,nloc,MPI_DOUBLE_PRECISION,0,comm,ierr)  
      tag = 1;
      if(myid.eq.0) then
         phi(1:nloc)=phip(1:nloc)
         do source=1, nprocs-1
            call MPI_RECV(phi( sum(nloc_array(1:source))+1), &
                 nloc_array(source+1), MPI_DOUBLE_COMPLEX, source, tag, comm, status, ierr)
         enddo
      else
         call MPI_SEND(phip(1), nloc, MPI_DOUBLE_COMPLEX, 0, tag, comm, ierr)
      endif
      call MPI_Bcast(phi,Hsize,MPI_DOUBLE_COMPLEX,0,comm,ierr)

      phip=0.0d0
      do kk=1,SprSize
         phip(IndexI(kk)-localstart+1)=phip(IndexI(kk)-localstart+1)+sparseH(kk)*phi(IndexJ(kk))
      enddo

      if(ii.ne.1) then
         do kk=1, nloc
           phip(kk) = phip(kk) - betha(ii)*phil(kk)
         enddo
      endif

      myalpha(ii)=0.0d0
      do kk=1, nloc
        myalpha(ii) = myalpha(ii) + conjg(phi(kk+localstart-1))*phip(kk)
      enddo
      call MPI_Allreduce(myalpha(ii), alpha(ii), 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
      !if(myid.eq.0) write(*,*) 'alpha', alpha(ii)
      do kk=1, nloc
        phipp(kk) = phip(kk) - alpha(ii)*phi(kk+localstart-1)
      enddo

   if(ii.ne.niter_CFE) then
      mybetha(ii+1) = 0.0d0
      do kk=1, nloc
        mybetha(ii+1) = mybetha(ii+1) + conjg(phipp(kk))*phipp(kk)
      enddo
      call MPI_Allreduce(mybetha(ii+1), betha(ii+1), 1, MPI_DOUBLE_COMPLEX,MPI_sum,comm,ierr)
      betha(ii+1) = sqrt(betha(ii+1))
      !if(myid.eq.0) write(*,*) 'betha', betha(ii+1)
      do kk=1, nloc
        phil(kk) = phi(kk+localstart-1)
      enddo
      do kk=1, nloc
        phi(kk+localstart-1) = phipp(kk) / betha(ii+1)
      enddo
   endif

enddo


! if (myid .eq. 0) then
!    do ii = 1, niter_CFE
!       write(*,*) "i, alpha, beta: ", ii, alpha(ii), betha(ii)
!    enddo
! endif

!***************alpha and betha ready******************
!                  compute I(w) 
 
call MPI_Barrier(comm, ierr)
if (myid .eq. 0) then
   ! write(*,*) "factor, gssum: ", factor, tgssum
   do ii=0, divY
      z = CMPLX(dble(ii)/divY*(endY-startY)+startY+E_0,epsilone_CFE)
      Intensity(ii+1)=z-alpha(niter_CFE)
      do jj=1,niter_CFE-1
         Intensity(ii+1)=z-alpha(niter_CFE-jj)-betha(niter_CFE-jj+1)**2/Intensity(ii+1)
         ! write(*,*) ii, Intensity(ii+1)
      enddo
      ! if (ii.eq.0) write(*,*) "z, int: ", CMPLX(dble(ii)/divY*(endY-startY)+startY+E_0,epsilone_CFE), Intensity(ii+1)
      ! write(*,*) "spec, gs: ", factor/Intensity(ii+1), tgssum/(cmplx(dble(ii)/divY*(endY-startY)+startY,epsilone_CFE))
      specX(ii+1)=dble(ii)/divY*(endY-startY)+startY
      specY(ii+1)=-1*1/pi*AIMAG(factor/Intensity(ii+1)- &
         tgssum/(cmplx(dble(ii)/divY*(endY-startY)+startY,epsilone_CFE)))
   enddo
endif

if (myid .eq. 0) write(*,*) 'inside CFE-Ngs done myid',myid
deallocate(phi, phil, phip, phipp)
deallocate(alpha,betha,myalpha,mybetha)

end subroutine 



!**********************************************************
integer*8 function ksubtract(kk, kkp)
use BettsCluster; use NumOfOrbitalAndElectrons
implicit none
integer*8 :: kk, kkp
integer*8 :: iktemp, kstat, ktempx, ktempy
integer*8 :: px, py

px = kpoint(kkp,1)
py = kpoint(kkp,2)
      ! code to get kkp
      iktemp=0
      kstat=1
      ktempx = mod(kpoint(kk,1)-px+2*kunit,kunit)
      ktempy = mod(kpoint(kk,2)-py+2*kunit,kunit)
      do while (kstat.ne.0)
         if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy)) then
            kstat=0
         else
            iktemp=iktemp+1
         endif
         if(iktemp.eq.N) then
            write(*,*) 'iktemp out of bounds, spin up'
            stop
         endif
      enddo
      ksubtract=iktemp
end function


!**********************************************************
integer*8 function ksum(kk, kkp)
use BettsCluster; use NumOfOrbitalAndElectrons
implicit none
integer*8 :: kk, kkp
integer*8 :: iktemp, kstat, ktempx, ktempy
integer*8 :: px, py

   px = kpoint(kkp,1)
   py = kpoint(kkp,2)
   ! code to get kkp
   iktemp=mod(kk,Nsite)
   ! write(*,*) "iktemp: ",iktemp
   kstat=1
   ktempx = mod(kpoint(kk,1)+px+kunit,kunit)
   ktempy = mod(kpoint(kk,2)+py+kunit,kunit)
   do while (kstat.ne.0)
      if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy)) then
         kstat=0
      else
         iktemp=iktemp+Nsite
      endif
      if(iktemp.eq.N) then
         write(*,*) 'iktemp out of bounds, ksum'
         stop
      endif
   enddo
   ksum=iktemp
   ! write(*,*) "ksum: ",ksum
end function

!**********************************************************
double complex function r2dotk2(ntag,qq1,qq2,r1x,r1y,r2x,r2y)
use BettsCluster; use NumOfOrbitalAndElectrons;
use ConstantParas; use MPIParas; use ModelParas;
implicit none
! e^(ik1*r1+ik2*r2), include reflections
integer*8, external :: ksubtract, ksum
integer :: i, j, x, y, ntag
integer*8 :: qq1, qq2
integer :: r1x, r1y, r2x, r2y
double precision :: k1x, k1y, k2x, k2y
double complex :: dotsum
integer, dimension(2) :: indx = [1,-1]

k1x = real(kpoint(qq1,1))*2.0d0*pi/kunit
k1y = real(kpoint(qq1,2))*2.0d0*pi/kunit
k2x = real(kpoint(qq2,1))*2.0d0*pi/kunit
k2y = real(kpoint(qq2,2))*2.0d0*pi/kunit
dotsum = cmplx(0.0d0,0.0d0)
! Only perform a mirror flip
! (x,y) -> (-x,-y)
dotsum = dotsum + cmplx(cos(k1x*real(r1x) + k1y*real(r1y) + k2x*real(r2x) &
                           + k2y*real(r2y)) &
                  , sin(k1x*real(r1x) + k1y*real(r1y) + k2x*real(r2x) &
                           + k2y*real(r2y)))
   
dotsum = dotsum + cmplx(cos(k1x*real(-r1x) + k1y*real(-r1y) + k2x*real(-r2x) &
                           + k2y*real(-r2y)) &
                  , sin(k1x*real(-r1x) + k1y*real(-r1y) + k2x*real(-r2x) &
                           + k2y*real(-r2y)))

r2dotk2 = dotsum
end function

!**********************************************************
double complex function lsf(q1,q2,ctag)
use BettsCluster; use NumOfOrbitalAndElectrons;
use ConstantParas; use MPIParas; use ModelParas;
implicit none
! Lattice Structure Factor
integer*8, external :: ksubtract, ksum
integer*8 :: q1, q2, ctag
double complex, external :: r2dotk2
double complex :: prefac

prefac = cmplx(0.0d0,0.0d0)
if (ctag.eq.0) then !A1g, n = 0
   prefac = prefac + r2dotk2(0,q1,q2,0,0,0,1) + r2dotk2(0,q1,q2,0,0,1,0) &
               + 2*tt*tt*r2dotk2(0,q1,q2,0,0,1,1) + 2*tt*tt*r2dotk2(0,q1,q2,0,0,1,-1) &
               + 2*ttt*ttt*r2dotk2(0,q1,q2,0,0,0,2) + 2*ttt*ttt*r2dotk2(0,q1,q2,0,0,2,0)
else if (ctag.eq.20) then !B1g, n = 0
   prefac = prefac - r2dotk2(0,q1,q2,0,0,0,1) + r2dotk2(0,q1,q2,0,0,1,0) &
               - 2*ttt*ttt*r2dotk2(0,q1,q2,0,0,0,2) + 2*ttt*ttt*r2dotk2(0,q1,q2,0,0,2,0)
else if (ctag.eq.30) then !B2g, n = 0
   prefac = prefac + 2*tt*tt*r2dotk2(0,q1,q2,0,0,1,-1) - 2*tt*tt*r2dotk2(0,q1,q2,0,0,1,1) 
else if (ctag.eq.2) then !A1g, n = 2
   prefac = prefac +4*r2dotk2(0,q1,q2,0,0,2,0)+4*r2dotk2(0,q1,q2,0,0,0,2) & 
                  +12*r2dotk2(0,q1,q2,0,0,1,1)+12*r2dotk2(0,q1,q2,0,0,-1,1)
else if (ctag.eq.22) then !B1g, n = 2
   prefac = prefac + 32*r2dotk2(0,q1,q2,0,0,1,0)- 32*r2dotk2(0,q1,q2,0,0,0,1) & 
                     + 4*r2dotk2(0,q1,q2,0,0,2,0) - 4*r2dotk2(0,q1,q2,0,0,0,2) !&
else if (ctag.eq.32) then !B2g, n = 2
   prefac = prefac + 12 *r2dotk2(0,q1,q2,0,0,1,1)-12*r2dotk2(0,q1,q2,0,0,-1,1)
else
   write(*,*) "invalid prefactor: ", ctag
endif 
if (myid.eq.0) write(*,*) 'q1, q2, tag, lsf', q1, q2, ctag, prefac

lsf = prefac
end function

!**********************************************************
double complex function rsdotks(qq1,qq2,qq3,r1x,r1y,r2x,r2y,r3x,r3y)
use BettsCluster; use NumOfOrbitalAndElectrons;
use ConstantParas; use MPIParas; use ModelParas;
implicit none
! e^(ik1*r1+ik2*r2), include reflections
integer*8, external :: ksubtract, ksum
integer :: i, j, x, y
integer*8 :: qq1, qq2, qq3
integer :: r1x, r1y, r2x, r2y, r3x, r3y
double precision :: k1x, k1y, k2x, k2y, k3x, k3y
double complex :: dotsum
integer, dimension(2) :: indx = [1,-1]

k1x = real(kpoint(qq1,1))*2.0d0*pi/kunit
k1y = real(kpoint(qq1,2))*2.0d0*pi/kunit
k2x = real(kpoint(qq2,1))*2.0d0*pi/kunit
k2y = real(kpoint(qq2,2))*2.0d0*pi/kunit
k3x = real(kpoint(qq3,1))*2.0d0*pi/kunit
k3y = real(kpoint(qq3,2))*2.0d0*pi/kunit
dotsum = cmplx(0.0d0,0.0d0)
! For chiral operators, mirror operator flips sign
! (x,y) -> (-y,x) -> (-x,-y) -> (y,-x)
dotsum = dotsum + cmplx(cos(k1x*real(r1x) + k1y*real(r1y) + k2x*real(r2x) &
                           + k2y*real(r2y) + k3x*real(r3x) + k3y*real(r3y)) &
                  , sin(k1x*real(r1x) + k1y*real(r1y) + k2x*real(r2x) &
                           + k2y*real(r2y) + k3x*real(r3x) + k3y*real(r3y)))
   
dotsum = dotsum + cmplx(cos(k1x*real(-r1x) + k1y*real(-r1y) + k2x*real(-r2x) &
                          + k2y*real(-r2y) + k3x*real(-r3x) + k3y*real(-r3y)) &
                 , sin(k1x*real(-r1x) + k1y*real(-r1y) + k2x*real(-r2x) &
                          + k2y*real(-r2y) + k3x*real(-r3x) + k3y*real(-r3y)))

dotsum = dotsum + cmplx(cos(k1x*real(r1y) + k1y*real(-r1x) + k2x*real(r2y) &
                          + k2y*real(-r2x) + k3x*real(r3y) + k3y*real(-r3x)) &
                 , sin(k1x*real(r1y) + k1y*real(-r1x) + k2x*real(r2y) &
                          + k2y*real(-r2x) + k3x*real(r3y) + k3y*real(-r3x)))

dotsum = dotsum + cmplx(cos(k1x*real(-r1y) + k1y*real(r1x) + k2x*real(-r2y) &
                          + k2y*real(r2x) + k3x*real(-r3y) + k3y*real(r3x)) &
                 , sin(k1x*real(-r1y) + k1y*real(r1x) + k2x*real(-r2y) &
                          + k2y*real(r2x) + k3x*real(-r3y) + k3y*real(r3x)))
rsdotks = dotsum
end function


!**********************************************************
double complex function lsf_chiral(q1,q2,q3,ctag)
use BettsCluster; use NumOfOrbitalAndElectrons;
use ConstantParas; use MPIParas; use ModelParas;
implicit none
! Chiral structure factor, (q-k-k')r1+k'r2
integer*8, external :: ksubtract, ksum
integer*8 :: q1, q2, q3, ctag
double complex, external :: rsdotks
double complex :: lsf

lsf = cmplx(0.0d0,0.0d0)
if (ctag.eq.2) then
   ! Hopping in 1 Unit cell
   lsf = lsf + (+8*rsdotks(q1,q2,q3,0,0,1,0,0,1)+24*rsdotks(q1,q2,q3,0,0,1,0,1,1) & 
               -24*rsdotks(q1,q2,q3,0,0,0,1,1,1)-8*rsdotks(q1,q2,q3,1,0,0,1,1,1)) !&
   !          * tt * tt
   ! ! Hopping in 2 cell, middle, 2t'
   lsf = lsf + (-4*rsdotks(q1,q2,q3,0,0,1,0,0,1)+12*rsdotks(q1,q2,q3,0,0,0,1,1,1) & 
               -4*rsdotks(q1,q2,q3,0,0,0,1,-1,0)+4*rsdotks(q1,q2,q3,1,0,0,1,-1,0) & 
               -12*rsdotks(q1,q2,q3,0,0,0,1,-1,1)-4*rsdotks(q1,q2,q3,0,0,1,1,-1,1)) !&
   !             * tt * tt
   ! ! Hopping in 2 cell, corner 2t'
   lsf = lsf + (-16*rsdotks(q1,q2,q3,0,0,1,0,1,1)+16*rsdotks(q1,q2,q3,0,0,0,1,1,1)) !&
   !             * tt * tt *2
   ! lsf = & !4*rsdotks(q1,q2,q3,0,0,1,0,0,1)-4*rsdotks(q1,q2,q3,0,0,0,1,-1,0) !& 
   !       ! +4*rsdotks(q1,q2,q3,1,0,0,1,-1,0)-4*rsdotks(q1,q2,q3,0,0,1,1,-1,1)! & 
   !       +4*rsdotks(q1,q2,q3,0,0,0,1,1,1)-8*rsdotks(q1,q2,q3,1,0,0,1,1,1) & 
   !       -12*rsdotks(q1,q2,q3,0,0,0,1,-1,1)+8*rsdotks(q1,q2,q3,0,0,1,0,1,1)
   ! Hopping in 2 cell, middle, t' + t''
   ! lsf = lsf - (16*rsdotks(q1,q2,q3,0,0,1,0,0,1)-8*rsdotks(q1,q2,q3,0,0,0,1,1,1) & 
   !             +16*rsdotks(q1,q2,q3,0,0,0,1,-1,0)+12*rsdotks(q1,q2,q3,1,0,0,1,-1,0) & 
   !             -4*rsdotks(q1,q2,q3,0,0,1,1,-1,0)+6*rsdotks(q1,q2,q3,1,0,1,1,-1,0) & 
   !             -4*rsdotks(q1,q2,q3,0,0,1,0,-1,1)+8*rsdotks(q1,q2,q3,0,0,0,1,-1,1) & 
   !             +2*rsdotks(q1,q2,q3,1,0,1,1,-1,1)-6*rsdotks(q1,q2,q3,1,0,-1,0,-1,1) & 
   !             -2*rsdotks(q1,q2,q3,1,1,-1,0,-1,1)) * tt * ttt
   ! lsf = lsf - (-40*rsdotks(q1,q2,q3,0,0,1,0,0,1)-8*rsdotks(q1,q2,q3,0,0,2,0,0,1) & 
   !             +6*rsdotks(q1,q2,q3,1,0,2,0,0,1)-12*rsdotks(q1,q2,q3,0,0,1,0,1,1) & 
   !             -8*rsdotks(q1,q2,q3,0,0,2,0,1,1)+8*rsdotks(q1,q2,q3,1,0,2,0,1,1) & 
   !             +12*rsdotks(q1,q2,q3,0,0,0,1,1,1)+2*rsdotks(q1,q2,q3,2,0,0,1,1,1) & 
   !             +2*rsdotks(q1,q2,q3,1,0,2,0,2,1)-2*rsdotks(q1,q2,q3,1,0,0,1,2,1) & 
   !             +2*rsdotks(q1,q2,q3,2,0,1,1,2,1)-8*rsdotks(q1,q2,q3,0,0,1,0,0,2) & 
   !             -6*rsdotks(q1,q2,q3,1,0,0,1,0,2)-8*rsdotks(q1,q2,q3,0,0,1,1,0,2) & 
   !             -2*rsdotks(q1,q2,q3,1,0,1,1,0,2)+8*rsdotks(q1,q2,q3,0,1,1,1,0,2) & 
   !             -2*rsdotks(q1,q2,q3,1,0,0,1,1,2)-2*rsdotks(q1,q2,q3,0,1,0,2,1,2) & 
   !             +2*rsdotks(q1,q2,q3,1,1,0,2,1,2)) * tt * ttt
   ! n = 4
else if (ctag.eq.4) then
   ! lsf = lsf +34*rsdotks(q1,q2,q3,0,0,1,0,0,1)+9*rsdotks(q1,q2,q3,0,0,1,0,2,1) & 
   !           +9*rsdotks(q1,q2,q3,0,0,1,2,0,1)+30*rsdotks(q1,q2,q3,0,1,2,1,0,1)&
   !           +9*rsdotks(q1,q2,q3,0,0,1,0,0,2)+9*rsdotks(q1,q2,q3,0,0,2,0,0,1)

   lsf = lsf +20*rsdotks(q1,q2,q3,0,0,1,0,0,1)+2*rsdotks(q1,q2,q3,0,0,1,0,1,1) & 
               -2*rsdotks(q1,q2,q3,0,0,0,1,1,1)+8*rsdotks(q1,q2,q3,1,0,0,1,1,1)
   lsf = lsf -4*rsdotks(q1,q2,q3,0,0,1,0,0,1)-12*rsdotks(q1,q2,q3,0,0,1,0,1,1) & 
            +8*rsdotks(q1,q2,q3,0,0,0,1,1,1)+2*rsdotks(q1,q2,q3,1,0,0,1,1,1) & 
            -4*rsdotks(q1,q2,q3,0,0,0,1,-1,0)-6*rsdotks(q1,q2,q3,1,0,0,1,-1,0) & 
            +10*rsdotks(q1,q2,q3,0,0,1,1,-1,0)-8*rsdotks(q1,q2,q3,0,1,1,1,-1,0) & 
            +10*rsdotks(q1,q2,q3,0,0,1,0,-1,1)-8*rsdotks(q1,q2,q3,0,0,0,1,-1,1) & 
            +8*rsdotks(q1,q2,q3,1,0,0,1,-1,1)+4*rsdotks(q1,q2,q3,0,0,1,1,-1,1) & 
            -2*rsdotks(q1,q2,q3,1,0,1,1,-1,1)+12*rsdotks(q1,q2,q3,0,0,-1,0,-1,1) & 
            +2*rsdotks(q1,q2,q3,0,1,-1,0,-1,1)+2*rsdotks(q1,q2,q3,1,1,-1,0,-1,1)
   lsf = lsf +16*rsdotks(q1,q2,q3,0,0,1,0,0,1)+6*rsdotks(q1,q2,q3,0,0,2,0,0,1) & 
               -3*rsdotks(q1,q2,q3,1,0,2,0,0,1)+10*rsdotks(q1,q2,q3,0,0,1,0,1,1) & 
               +4*rsdotks(q1,q2,q3,0,0,2,0,1,1)+2*rsdotks(q1,q2,q3,1,0,2,0,1,1) & 
               -10*rsdotks(q1,q2,q3,0,0,0,1,1,1)+10*rsdotks(q1,q2,q3,1,0,0,1,1,1) & 
               +6*rsdotks(q1,q2,q3,2,0,0,1,1,1)-2*rsdotks(q1,q2,q3,0,0,1,0,2,1) & 
               -rsdotks(q1,q2,q3,1,0,2,0,2,1)-12*rsdotks(q1,q2,q3,1,0,0,1,2,1) & 
               -5*rsdotks(q1,q2,q3,2,0,0,1,2,1)-2*rsdotks(q1,q2,q3,0,0,1,1,2,1) & 
               -13*rsdotks(q1,q2,q3,1,0,1,1,2,1)-8*rsdotks(q1,q2,q3,2,0,1,1,2,1) & 
               +6*rsdotks(q1,q2,q3,0,0,1,0,0,2)+3*rsdotks(q1,q2,q3,1,0,0,1,0,2) & 
               +4*rsdotks(q1,q2,q3,0,0,1,1,0,2)-6*rsdotks(q1,q2,q3,1,0,1,1,0,2) & 
               +2*rsdotks(q1,q2,q3,0,1,1,1,0,2)+2*rsdotks(q1,q2,q3,0,0,0,1,1,2) & 
               -12*rsdotks(q1,q2,q3,1,0,0,1,1,2)+2*rsdotks(q1,q2,q3,0,0,1,1,1,2) & 
               +13*rsdotks(q1,q2,q3,0,1,1,1,1,2)-5*rsdotks(q1,q2,q3,1,0,0,2,1,2) & 
               +rsdotks(q1,q2,q3,0,1,0,2,1,2)-8*rsdotks(q1,q2,q3,1,1,0,2,1,2)
else
   write(*,*) "invalid prefactor (Has to be: 2,4)", ctag
endif 
if (myid.eq.0) write(*,*) 'q1, q2, q3, tag, lsf', q1, q2, q3, ctag, lsf
lsf_chiral = lsf
end function

integer*8 function nchoosek(n,k)
integer*8, intent(in) :: n,k
integer*8 :: i, numerator, denominator
numerator = 1
denominator = 1
! Calculate C(n,k)
do i = 1, min(k, n - k)
   numerator = numerator * (n - i + 1)
   denominator = denominator * i
end do
nchoosek = numerator / denominator
end function
