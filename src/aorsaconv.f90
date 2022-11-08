! Convert results into a single netcdf file named fpm.nc
! Files read:
! - fpm
! - transform.out
! - movie_wdot
! - out_edotj

program transformconv
use netcdf

implicit none
integer :: iargc
integer, parameter :: CHARLEN = 128
character(len=CHARLEN) :: file, fileOut="fpm.nc"

integer :: NULL(0)
integer, parameter :: LUNERR = 0
integer, parameter :: LUNOUT = 6
integer :: i, j, k, n
integer :: istat, istat200

integer :: id_nc, id_xdim, id_ydim, id_phidim, id_rhodim, id_spdim, &
      id_mcapdim, id_x, id_capr, id_y, id_phiprime, id_psilim, id_rt, &
      id_b0, id_psi, id_bmod, id_rhon, id_xnavg, id_mcap, &
      id_reealpha, id_imealpha, id_reebeta, id_imebeta, &
      id_reeb, id_imeb, id_rene, id_imne, id_rekparae, id_imkparae, &
      id_ptot, id_pcrto2, id_pcito2, &
      id_redotjavg, id_wdotavg, id_wdot_ql, &
      id_wdote, id_wdot1, id_wdot2, id_wdot3, &
      id_edotje, id_edotj1, id_edotj2, id_edotj3, &
      id_rho, id_thetap, id_bxn, id_byn, id_bzn

integer, dimension(4) :: sz

integer :: nnodex, nnodey, nnodephi, nnoderho, nmcap, mcap, imcap
integer :: nsp = 4, nfp
real :: psilim, rt, b0, ptot, pcrto2, pcito2
real, dimension(:), allocatable :: x, y, capr, phiprime
real, dimension(:), allocatable :: rhon, xnavg
real, dimension(:, :), allocatable:: work2
real, dimension(:, :, :), allocatable:: work3
real, dimension(:, :, :), allocatable :: bxn, byn, bzn, bmod, psi, rho, thetap
complex, dimension(:, :, :), allocatable :: cwork

! n = iargc()
! if (n >= 1) then
!     call getarg(1, fileOut)
! endif

309 format(10i10)
310 format(1p6e12.4)

! create netcdf file
istat = nf90_create(fileOut, NF90_CLOBBER, id_nc)
if (istat /= NF90_NOERR) then
    write(LUNERR, *) trim(nf90_strerror(istat))
    stop
endif

! --------------------------------------------
! variable definition
! --------------------------------------------

! -- fpm --
file = "fpm"

open(200, file=file, status="old", iostat=istat200)
if (istat200 == 0) then

    read(200, 309) nmcap
    read(200, 309) nnodex, nnodey, nnodephi
    nnoderho = nnodex / 2.

    ! dimensions
    istat = nf90_def_dim(id_nc, "x", nnodex, id_xdim)
    istat = nf90_def_dim(id_nc, "y", nnodey, id_ydim)
    istat = nf90_def_dim(id_nc, "phi", nnodephi, id_phidim)
    istat = nf90_def_dim(id_nc, "rho", nnoderho, id_rhodim)
    istat = nf90_def_dim(id_nc, "sp", nsp, id_spdim)
    istat = nf90_def_dim(id_nc, "mcap", nmcap, id_mcapdim)

    ! variables
    istat = nf90_def_var(id_nc, "x", NF90_DOUBLE, id_xdim, id_x)
    istat = nf90_def_var(id_nc, "capr", NF90_DOUBLE, id_xdim, id_capr)
    istat = nf90_def_var(id_nc, "y", NF90_DOUBLE, id_ydim, id_y)
    istat = nf90_def_var(id_nc, "phiprime", NF90_DOUBLE, id_phidim, id_phiprime)

    istat = nf90_def_var(id_nc, "psilim", NF90_DOUBLE, NULL, id_psilim)
    istat = nf90_def_var(id_nc, "rt", NF90_DOUBLE, NULL, id_rt)
    istat = nf90_def_var(id_nc, "b0", NF90_DOUBLE, NULL, id_b0)

    sz = (/id_xdim, id_ydim, id_phidim, id_mcapdim/)
    istat = nf90_def_var(id_nc, "psi", NF90_DOUBLE, sz(1:3), id_psi)
    istat = nf90_def_var(id_nc, "bmod", NF90_DOUBLE, sz(1:3), id_bmod)

    istat = nf90_def_var(id_nc, "rhon", NF90_DOUBLE, id_rhodim, id_rhon)
    istat = nf90_def_var(id_nc, "xnavg", NF90_DOUBLE, id_rhodim, id_xnavg)

    istat = nf90_def_var(id_nc, "mcap", NF90_INT, id_mcapdim, id_mcap)

    istat = nf90_def_var(id_nc, "Reealpha", NF90_DOUBLE, sz, id_reealpha)
    istat = nf90_def_var(id_nc, "Imealpha", NF90_DOUBLE, sz, id_imealpha)
    istat = nf90_def_var(id_nc, "Reebeta", NF90_DOUBLE, sz, id_reebeta)
    istat = nf90_def_var(id_nc, "Imebeta", NF90_DOUBLE, sz, id_imebeta)
    istat = nf90_def_var(id_nc, "Reeb", NF90_DOUBLE, sz, id_reeb)
    istat = nf90_def_var(id_nc, "Imeb", NF90_DOUBLE, sz, id_imeb)
    istat = nf90_def_var(id_nc, "Rene", NF90_DOUBLE, sz, id_rene)
    istat = nf90_def_var(id_nc, "Imne", NF90_DOUBLE, sz, id_imne)
    istat = nf90_def_var(id_nc, "Rekparae", NF90_DOUBLE, sz, id_rekparae)
    istat = nf90_def_var(id_nc, "Imkparae", NF90_DOUBLE, sz, id_imkparae)

    istat = nf90_def_var(id_nc, "ptot", NF90_DOUBLE, id_mcapdim, id_ptot)
    istat = nf90_def_var(id_nc, "pcrto2", NF90_DOUBLE, id_mcapdim, id_pcrto2)
    istat = nf90_def_var(id_nc, "pcito2", NF90_DOUBLE, id_mcapdim, id_pcito2)

    istat = nf90_def_var(id_nc, "redotjavg", NF90_DOUBLE, &
        & (/id_rhodim, id_spdim, id_mcapdim/), id_redotjavg)

    istat = nf90_def_var(id_nc, "wdotavg", NF90_DOUBLE, &
        & (/id_rhodim, id_spdim, id_mcapdim/), id_wdotavg)
    istat = nf90_def_var(id_nc, "wdot_ql", NF90_DOUBLE, &
        & (/id_rhodim, id_spdim, id_mcapdim/), id_wdot_ql)
else
    write(LUNERR, *) "Error opening "//file
endif

! ------------------------------------------------------------------

! -- transform.out --
file = "transform.out"

open(100, file=file, status="old", iostat=istat)
if (istat /= 0) then
    write(LUNERR, *) "Error opening "//file
    stop
endif

read(100, 309) nnodex, nnodey, nnodephi

! dimensions
if (istat200 /= 0) then
    istat = nf90_def_dim(id_nc, "x", nnodex, id_xdim)
    istat = nf90_def_dim(id_nc, "y", nnodey, id_ydim)
    istat = nf90_def_dim(id_nc, "phi", nnodephi, id_phidim)
endif

! variables
istat = nf90_def_var(id_nc, "rt", NF90_DOUBLE, NULL, id_rt)
istat = nf90_def_var(id_nc, "b0", NF90_DOUBLE, NULL, id_b0)

istat = nf90_def_var(id_nc, "bxn", NF90_DOUBLE, sz(1:3), id_bxn)
istat = nf90_def_var(id_nc, "byn", NF90_DOUBLE, sz(1:3), id_byn)
istat = nf90_def_var(id_nc, "bzn", NF90_DOUBLE, sz(1:3), id_bzn)
istat = nf90_def_var(id_nc, "bmod", NF90_DOUBLE, sz(1:3), id_bmod)
istat = nf90_def_var(id_nc, "psi", NF90_DOUBLE, sz(1:3), id_psi)
istat = nf90_def_var(id_nc, "rho", NF90_DOUBLE, sz(1:3), id_rho)
istat = nf90_def_var(id_nc, "thetap", NF90_DOUBLE, sz(1:3), id_thetap)
! ------------------------------------------------------------------


! -- movie_wdot --
file = "movie_wdot"

open(101, file=file, status="old", iostat=istat)
if (istat /= 0) then
    write(LUNERR, *) "Error opening "//file
    stop
endif

read(101, 309) nmcap
read(101, 309) nnodex, nnodey, nnodephi

! dimensions
if (istat200 /= 0) then
    istat = nf90_def_dim(id_nc, "mcap", nmcap, id_mcapdim)
endif

! variables
istat = nf90_def_var(id_nc, "wdote", NF90_DOUBLE, sz, id_wdote)
istat = nf90_def_var(id_nc, "wdot1", NF90_DOUBLE, sz, id_wdot1)
istat = nf90_def_var(id_nc, "wdot2", NF90_DOUBLE, sz, id_wdot2)
istat = nf90_def_var(id_nc, "wdot3", NF90_DOUBLE, sz, id_wdot3)
! ------------------------------------------------------------------


! -- out_edotj --
file = "out_edotj"

open(102, file=file, status="old", iostat=istat)
if (istat /= 0) then
    write(LUNERR, *) "Error opening "//file
    stop
endif

read(102, 309) nmcap
read(102, 309) nnodex, nnodey, nnodephi


! variables
istat = nf90_def_var(id_nc, "edotje", NF90_DOUBLE, sz, id_edotje)
istat = nf90_def_var(id_nc, "edotj1", NF90_DOUBLE, sz, id_edotj1)
istat = nf90_def_var(id_nc, "edotj2", NF90_DOUBLE, sz, id_edotj2)
istat = nf90_def_var(id_nc, "edotj3", NF90_DOUBLE, sz, id_edotj3)
! ------------------------------------------------------------------


! -- end definitions --
istat = nf90_enddef(id_nc)
if (istat /= NF90_NOERR) then
    write(LUNERR, *) trim(nf90_strerror(istat))
    istat = nf90_close(id_nc)
    stop
endif

! --------------------------------------------
! fpm
! --------------------------------------------
if (istat200 == 0) then

    allocate(x(nnodex))
    allocate(y(nnodey))
    allocate(capr(nnodex))
    allocate(phiprime(nnodephi))
    allocate(psi(nnodex, nnodey, nnodephi))
    allocate(bmod(nnodex, nnodey, nnodephi))
    read(200, 310) (x(i), i = 1, nnodex)
    read(200, 310) (y(j), j = 1, nnodey)
    read(200, 310) (capr(i), i = 1, nnodex)
    read(200, 310) (phiprime(i), i = 1, nnodephi)
    read(200, 310) psilim, rt, b0
    read(200, 310) (((psi(i, j, k),  i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
    read(200, 310) (((bmod(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)

    read(200, 309) nnoderho

    allocate(rhon(nnoderho))
    allocate(xnavg(nnoderho))

    read(200, 310) (rhon(i), i = 1, nnoderho)
    read(200, 310) (xnavg(i), i = 1, nnoderho)

    istat = nf90_put_var(id_nc, id_x, x)
    istat = nf90_put_var(id_nc, id_capr, capr)
    istat = nf90_put_var(id_nc, id_y, y)
    istat = nf90_put_var(id_nc, id_phiprime, phiprime)
    istat = nf90_put_var(id_nc, id_psilim, psilim)
    istat = nf90_put_var(id_nc, id_rt, rt)
    istat = nf90_put_var(id_nc, id_b0, b0)
    istat = nf90_put_var(id_nc, id_psi, psi)
    istat = nf90_put_var(id_nc, id_bmod, bmod)
    istat = nf90_put_var(id_nc, id_rhon, rhon)
    istat = nf90_put_var(id_nc, id_xnavg, xnavg)

    allocate(cwork(nnodex, nnodey, nnodephi))
    allocate(work2(nnoderho, nsp))

    do imcap = 1, nmcap
        read(200, 309) mcap
        istat = nf90_put_var(id_nc, id_mcap, (/mcap/), start=(/imcap/), count=(/1/))
        read(200, 310) (((cwork(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
        istat = nf90_put_var(id_nc, id_reealpha, real(cwork),  & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
        istat = nf90_put_var(id_nc, id_imealpha, aimag(cwork), & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
        read(200, 310) (((cwork(i, j, k), i = 1, nnodex), j = 1, nnodey), k = 1, nnodephi)
        istat = nf90_put_var(id_nc, id_reebeta, real(cwork),  & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
        istat = nf90_put_var(id_nc, id_imebeta, aimag(cwork), & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
        read(200, 310) (((cwork(i, j, k), i = 1, nnodex), j = 1, nnodey), k = 1, nnodephi)
        istat = nf90_put_var(id_nc, id_reeb, real(cwork),  & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
        istat = nf90_put_var(id_nc, id_imeb, aimag(cwork), & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
        read(200, 310) (((cwork(i, j, k), i = 1, nnodex), j = 1, nnodey), k = 1, nnodephi)
        istat = nf90_put_var(id_nc, id_rene, real(cwork),  & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
        istat = nf90_put_var(id_nc, id_imne, aimag(cwork), & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))

        ! AORSA isn't saving kparae to fpm
        ! read(200, 310) (((cwork(i, j, k), i = 1, nnodex), j = 1, nnodey), k = 1, nnodephi)
        ! istat = nf90_put_var(id_nc, id_rekparae, real(cwork),  & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
        ! istat = nf90_put_var(id_nc, id_imkparae, aimag(cwork), & start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))

        read(200, 310) ptot, pcrto2, pcito2
        istat = nf90_put_var(id_nc, id_ptot, (/ptot/), & start=(/imcap/), count=(/1/))
        istat = nf90_put_var(id_nc, id_pcrto2, (/pcrto2/), & start=(/imcap/), count=(/1/))
        istat = nf90_put_var(id_nc, id_pcito2, (/pcito2/), & start=(/imcap/), count=(/1/))

        do j = 1, nsp
            read(200, 310) (work2(i, j), i = 1, nnoderho)
        enddo
        istat = nf90_put_var(id_nc, id_redotjavg, work2, & start=(/1, 1, imcap/), count=(/nnoderho, nsp, 1/))

        do j = 1, nsp
            read(200, *) (work2(i, j), i = 1, nnoderho)
        enddo
        istat = nf90_put_var(id_nc, id_wdotavg, work2, & start=(/1, 1, imcap/), count=(/nnoderho, nsp, 1/))

        do j = 1, nsp
            read(200, 310) (work2(i, j), i = 1, nnoderho)
        enddo
        istat = nf90_put_var(id_nc, id_wdot_ql, work2, & start=(/1, 1, imcap/), count=(/nnoderho, nsp, 1/))
    enddo
    close(200)
    endif

! --------------------------------------------
! transform.out
! --------------------------------------------

allocate(bxn(nnodex, nnodey, nnodephi))
allocate(byn(nnodex, nnodey, nnodephi))
allocate(bzn(nnodex, nnodey, nnodephi))
allocate(rho(nnodex, nnodey, nnodephi))
allocate(thetap(nnodex, nnodey, nnodephi))

if (istat200 /= 0) then
    allocate(bmod(nnodex, nnodey, nnodephi))
    allocate(psi(nnodex, nnodey, nnodephi))
endif

read(100, 310) rt, b0
read(100, 310) (((bxn(i,j,k),    i=1,nnodex), j=1,nnodey), k=1,nnodephi)
read(100, 310) (((byn(i,j,k),    i=1,nnodex), j=1,nnodey), k=1,nnodephi)
read(100, 310) (((bzn(i,j,k),    i=1,nnodex), j=1,nnodey), k=1,nnodephi)
read(100, 310) (((bmod(i,j,k),   i=1,nnodex), j=1,nnodey), k=1,nnodephi)
read(100, 310) (((psi(i,j,k),    i=1,nnodex), j=1,nnodey), k=1,nnodephi)
read(100, 310) (((rho(i,j,k),    i=1,nnodex), j=1,nnodey), k=1,nnodephi)
read(100, 310) (((thetap(i,j,k), i=1,nnodex), j=1,nnodey), k=1,nnodephi)
read(100, *) nfp
close(100)

istat = nf90_put_var(id_nc, id_rt, rt)
istat = nf90_put_var(id_nc, id_b0, b0)
istat = nf90_put_var(id_nc, id_bxn, bxn)
istat = nf90_put_var(id_nc, id_byn, byn)
istat = nf90_put_var(id_nc, id_bzn, bzn)
istat = nf90_put_var(id_nc, id_bmod, bmod)
istat = nf90_put_var(id_nc, id_psi, psi)
istat = nf90_put_var(id_nc, id_rho, rho)
istat = nf90_put_var(id_nc, id_thetap, thetap)


! --------------------------------------------
! movie_wdot
! --------------------------------------------

allocate(work3(nnodex, nnodey, nnodephi))

do imcap = 1, nmcap
    read(101, 310) (((work3(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
    istat = nf90_put_var(id_nc, id_wdote, work3,  start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
    read(101, 310) (((work3(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
    istat = nf90_put_var(id_nc, id_wdot1, work3,  start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
    read(101, 310) (((work3(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
    istat = nf90_put_var(id_nc, id_wdot2, work3,  start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
    read(101, 310) (((work3(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
    istat = nf90_put_var(id_nc, id_wdot3, work3,  start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
enddo
close(101)

! --------------------------------------------
! out_edotj
! --------------------------------------------

do imcap = 1, nmcap
    read(102, 310) (((work3(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
    istat = nf90_put_var(id_nc, id_edotje, work3,  start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
    read(102, 310) (((work3(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
    istat = nf90_put_var(id_nc, id_edotj1, work3,  start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
    read(102, 310) (((work3(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
    istat = nf90_put_var(id_nc, id_edotj2, work3,  start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
    read(102, 310) (((work3(i, j, k), i = 1, nnodex),j = 1,nnodey), k = 1, nnodephi)
    istat = nf90_put_var(id_nc, id_edotj3, work3,  start=(/1, 1, 1, imcap/), count=(/nnodex, nnodey, nnodephi, 1/))
enddo
close(102)


! --------------------------------------------
! close netcdf file
! --------------------------------------------

istat = nf90_close(id_nc)

end program
