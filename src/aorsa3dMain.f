
*------------------------------------------------------------------------------
*     This version (4/01/03: newlab3) does not solve equations for the points
*     outside psi = 1.
*------------------------------------------------------------------------------


*------------------------------------------------------------------------
*     Stix2 frame using SCALAPACK with new faster loading by E. D`Azevedo
*------------------------------------------------------------------------

      program aorsa3dMain

      use ql_myra_mod
      use profile_mod

      use aorsa3din_mod

      implicit none

      integer n_upper, n_lower, m_upper, m_lower, l_upper, l_lower

      integer i0, j0, k0
      real rhomin

      integer iflag
      integer imcap
      integer nxmx, nymx, nphimx
      integer ier, l, irnc
      integer nrow, ncol, norder

      integer kplot1, kplot2, kplot3, kplot4, ndkplot

      integer nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2

      integer nkx1, nkx2, nldim, nldim3, nky1, nky2, iant, nphi1, nphi2

      real :: kperp_max
      real  kperp_max_actual

      real ttotal, time0, time00, time, dummy, second1
      real time1, time2, tmin1, tmin2, dum
      real tmin, gflops, gflopsp, ops, theta_antr

      real signb, sqx

      real capa, capx, capy, z, hz, coshz, sinhz, r2,
     .    bcapx, bcapy, bcapz, xkh

      integer nmodesmax, mmodesmax, lmodesmax, nrhomax
      real xkh0, xlp
      real epsl, r, angle, sinang, cosang

      real xkx_cutoff, xky_cutoff, xkz_cutoff

      real UminPara, UmaxPara
      real vce_mks, vc1_mks, vc2_mks, vc3_mks

      integer, parameter :: lmaxdim = 25

      integer nbessj
      integer, parameter :: n_psi_dim = 200
      integer :: n_psi
      real :: rho_a(n_psi_dim)

*     ----------------------
*     96 x 96 x 64 modes max
*     ----------------------
      parameter (nmodesmax = 96)
      parameter (mmodesmax = 96)
      parameter (lmodesmax = 160)


      parameter (nxmx = nmodesmax)
      parameter (nymx = mmodesmax)
      parameter (nphimx = lmodesmax)

      parameter (nrhomax = nmodesmax)

      parameter (nkdim1 = - nmodesmax / 2)
      parameter (nkdim2 =   nmodesmax / 2)

      parameter (mkdim1 = - mmodesmax / 2)
      parameter (mkdim2 =   mmodesmax / 2)

      parameter (lkdim1 = - lmodesmax / 2)
      parameter (lkdim2 =   lmodesmax / 2)

      parameter (nldim  = nxmx * nymx * nphimx)
      parameter (nldim3 = 3 * nldim)


      real, dimension(:),   allocatable :: UPERP, UPARA
      real, dimension(:),   allocatable :: VPERP, VPARA
      real, dimension(:),   allocatable :: UPERP_work, UPARA_work

      real, dimension(:,:), allocatable :: DFDUPERe, DFDUPARe
      real, dimension(:,:), allocatable :: DFDUPER1, DFDUPAR1
      real, dimension(:,:), allocatable :: DFDUPER2, DFDUPAR2
      real, dimension(:,:), allocatable :: DFDUPER3, DFDUPAR3

      real, dimension(:,:), allocatable :: f

      real, dimension(:,:,:), allocatable :: fe_cql_cart
      real, dimension(:,:,:), allocatable :: dfe_cql_uprp
      real, dimension(:,:,:), allocatable :: dfe_cql_uprl

      real, dimension(:,:,:), allocatable :: f1_cql_cart
      real, dimension(:,:,:), allocatable :: df1_cql_uprp
      real, dimension(:,:,:), allocatable :: df1_cql_uprl

      real, dimension(:,:,:), allocatable :: f2_cql_cart
      real, dimension(:,:,:), allocatable :: df2_cql_uprp
      real, dimension(:,:,:), allocatable :: df2_cql_uprl

      real, dimension(:,:,:), allocatable :: f3_cql_cart
      real, dimension(:,:,:), allocatable :: df3_cql_uprp
      real, dimension(:,:,:), allocatable :: df3_cql_uprl

      real, dimension(:,:,:), allocatable :: bqlavg_e
      real, dimension(:,:,:), allocatable :: cqlavg_e
      real, dimension(:,:,:), allocatable :: eqlavg_e
      real, dimension(:,:,:), allocatable :: fqlavg_e

      real, dimension(:,:,:), allocatable :: bqlavg_i1
      real, dimension(:,:,:), allocatable :: cqlavg_i1
      real, dimension(:,:,:), allocatable :: eqlavg_i1
      real, dimension(:,:,:), allocatable :: fqlavg_i1

      real, dimension(:,:,:), allocatable :: bqlavg_i2
      real, dimension(:,:,:), allocatable :: cqlavg_i2
      real, dimension(:,:,:), allocatable :: eqlavg_i2
      real, dimension(:,:,:), allocatable :: fqlavg_i2

      real, dimension(:,:,:), allocatable :: bqlavg_i3
      real, dimension(:,:,:), allocatable :: cqlavg_i3
      real, dimension(:,:,:), allocatable :: eqlavg_i3
      real, dimension(:,:,:), allocatable :: fqlavg_i3

c***  3-D arrays:
      complex, dimension(:,:,:), allocatable :: xjx, xjy, xjz

      real, dimension(:,:,:), allocatable :: psi, thetap, rho, xkte,
     .     xkti, xkti2, xkti3, xnea, xn1a, xn2a, xn3a,
     .     omgce, omgci1, omgci2, omgci3,
     .     omgpe2, omgp12, omgp22, omgp32, bmod, bmod_mid

      complex, dimension(:,:,:), allocatable :: xjpxe, xjpye, xjpze,
     .     xjpx1, xjpy1, xjpz1, xjpx2, xjpy2, xjpz2,
     .     xjpx3, xjpy3, xjpz3, xjpx, xjpy, xjpz

      real, dimension(:,:,:), allocatable :: pcre, pcim, redotj,
     .     redotje, redotj1, redotj2,redotj3,redotjt

      real, dimension(:,:,:), allocatable :: bxn, byn, bzn

      real, dimension(:,:,:), allocatable :: uxx, uxy, uxz,
     .     uyx, uyy, uyz, uzx, uzy, uzz
      real, dimension(:,:,:), allocatable :: dxuxx, dxuxy, dxuxz,
     .     dxuyx, dxuyy, dxuyz,
     .     dxuzx, dxuzy, dxuzz
      real, dimension(:,:,:), allocatable :: dyuxx, dyuxy, dyuxz,
     .     dyuyx, dyuyy, dyuyz,
     .     dyuzx, dyuzy, dyuzz
      real, dimension(:,:,:), allocatable :: dyyuxx, dyyuxy, dyyuxz,
     .     dyyuyx, dyyuyy, dyyuyz,
     .     dyyuzx, dyyuzy, dyyuzz
      real, dimension(:,:,:), allocatable :: dxyuxx, dxyuxy, dxyuxz,
     .     dxyuyx, dxyuyy, dxyuyz,
     .     dxyuzx, dxyuzy, dxyuzz
      real, dimension(:,:,:), allocatable :: dxxuxx, dxxuxy, dxxuxz,
     .     dxxuyx, dxxuyy, dxxuyz,
     .     dxxuzx, dxxuzy, dxxuzz
      real, dimension(:,:,:), allocatable :: dphiuxx, dphiuxy, dphiuxz,
     .     dphiuyx, dphiuyy, dphiuyz,
     .     dphiuzx, dphiuzy, dphiuzz
      real, dimension(:,:,:), allocatable :: d2phiuxx,d2phiuxy,d2phiuxz,
     .     d2phiuyx, d2phiuyy, d2phiuyz,
     .     d2phiuzx, d2phiuzy, d2phiuzz
      real, dimension(:,:,:), allocatable :: dxphiuxx,dxphiuxy,dxphiuxz,
     .     dxphiuyx, dxphiuyy, dxphiuyz,
     .     dxphiuzx, dxphiuzy, dxphiuzz
      real, dimension(:,:,:), allocatable :: dyphiuxx,dyphiuxy,dyphiuxz,
     .     dyphiuyx, dyphiuyy, dyphiuyz,
     .     dyphiuzx, dyphiuzy, dyphiuzz

      real, dimension(:,:,:), allocatable :: gradprlb


      complex, dimension(:,:,:), allocatable :: ealphak, ebetak, ebk

      complex, dimension(:,:,:), allocatable :: ex, ey, ez
      complex, dimension(:,:,:), allocatable :: ealpha, ebeta, eb
      complex, dimension(:,:,:), allocatable :: ntilda_e

      complex, dimension(:,:), allocatable :: xx, yy, zz

      complex one

      integer info
      complex zi, cexpkxkykz
      complex norm2, alpha
      integer incX, global_row, global_col

      complex fdk, fek, ffk, fgk, fak, fpk, frk, fqk, fsk

      complex arg, xil(0:100), xilp(0:100)
      complex bl(100)


      real xkxsav(nkdim1 : nkdim2),
     .     xkysav(mkdim1 : mkdim2),
     .     xkzsav(lkdim1 : lkdim2),
     .     ptot, ptot_wdot, pcito2, pcrto2, powtot, pscale,
     .     t1, gaussian, frho, q07qa, psimax,
     .     rhomax, rhoant, gausspsi, dpsiant,
     .     deltay, shapey, yant_max, zstrap, zn, rant1
      real phiant, xjymax, modxjy
      complex shapephi(nphimx)

      real ealphakmod(nkdim1 : nkdim2, mkdim1 : mkdim2),
     .     ebetakmod(nkdim1 : nkdim2, mkdim1 : mkdim2),
     .     ebkmod(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real ealphakmod1(lkdim1 : lkdim2),
     .     ebetakmod1(lkdim1 : lkdim2),
     .     ebkmod1(lkdim1 : lkdim2)

      integer nnodex, nnodey, nnodephi, i, j, k,
     .    jequat, icenter, np

      integer nnoderho

      integer n, m, nphi

      real xant, btor, rhoplasm,
     .    reomg1, reomg2, reomg3,
     .    phimax,
     .    rhonorm, xiota0, rholim, psi1,
     .    xwleft, xwright, psiright,
     .    psi_lim
      real telimj, tilimj, ti2limj, ti3limj
      real xkphi, xkphi0, xnphi0
      real rplasm, rlim

      real xkthrho, wphase, xkthdx,
     .   v0i, vthe, vphase, rnz, xn2, eta2,
     .   xk0, eta1, xn1, xmi1, xmh, xme, qi1, xmi2,
     .   xmi3, t0i, t0i2, t0i3, t0e, q, clight, xmu0, xlnlam,
     .   omgrf, xmax, qe, qi2, pi, two_pi, eps0, xn3, qi3,
     .   costh, sinth, rnx,  rny, rnphi

      real xjantx, xjanty, xjantz, xjant

      complex xb, xc, xd

      real xprimec(nxmx), caprc(nxmx), xcourse(nxmx), capr(nxmx),
     .   xprime(nxmx), x(nxmx), dx, dxc

      real rhon(nrhomax), drho, wdoti1avg(nrhomax), wdoti2avg(nrhomax),
     .     wdoti3avg(nrhomax), wdoteavg(nrhomax), wdotavg(nrhomax),
     .     dvol(nrhomax), bmod_min(nrhomax)
      real wdote_ql(nrhomax),  wdoti1_ql(nrhomax), wdoti2_ql(nrhomax),
     .     wdoti3_ql(nrhomax)
      real xnavg(nrhomax), xn1avg(nrhomax), xn2avg(nrhomax),
     .     xn3avg(nrhomax)
      real xkteavg(nrhomax), xktiavg(nrhomax), xkti2avg(nrhomax),
     .     xkti3avg(nrhomax)
      real qhat
      real dldbavg(nrhomax)
      real redotjeavg(nrhomax), redotj1avg(nrhomax),
     .     redotj2avg(nrhomax), redotj3avg(nrhomax)

      real vol(nrhomax), fvol(nrhomax)

      real theta(nxmx, nymx), theta0(nxmx, nymx),
     .     bx, by, bz

      real bphi, bth, br


      real dbdx, dbdy, dbdphi

      complex adisp(nxmx, nymx), acold(nxmx, nymx)

      complex pc

      real, dimension(:,:,:), allocatable :: reson1, reson2

      real, dimension(:,:,:), allocatable :: wdote, wdoti1, wdoti2,
     .     wdoti3, wdot

      real, dimension(:,:,:), allocatable :: fx0e, fx0i1, fx0i2, fx0i3,
     .     fy0e, fy0i1, fy0i2, fy0i3,
     .     fz0e, fz0i1, fz0i2, fz0i3, fx0, fy0, fz0

      real denom, fz0tot


      real pcedotj1, pcedotje, pcedotj2, pcedotjt, pcedotj3
      real pedotj1, pedotje, pedotj2, pedotjt, pedotj3

      real pi1, pi2, pi3, pe, pt
      real pcti1, pcte, pcti2, pctt, pcti3

      complex sigexx, sigexy, sigexz,
     .        sigeyx, sigeyy, sigeyz,
     .        sigezx, sigezy, sigezz

      complex sig1xx, sig1xy, sig1xz,
     .        sig1yx, sig1yy, sig1yz,
     .        sig1zx, sig1zy, sig1zz

      complex sig2xx, sig2xy, sig2xz,
     .        sig2yx, sig2yy, sig2yz,
     .        sig2zx, sig2zy, sig2zz

      complex sig3xx, sig3xy, sig3xz,
     .        sig3yx, sig3yy, sig3yz,
     .        sig3zx, sig3zy, sig3zz


      complex
     .     sigxx, sigxy, sigxz,
     .     sigyx, sigyy, sigyz,
     .     sigzx, sigzy, sigzz

      complex
     .     xkxx, xkxy, xkxz,
     .     xkyx, xkyy, xkyz,
     .     xkzx, xkzy, xkzz
      complex scold

      complex
     .     dxx, dxy, dxz,
     .     dyx, dyy, dyz,
     .     dzx, dzy, dzz

      real yprimec(nymx), ycourse(nymx),
     .     yprime(nymx), y(nymx), dy, dyc

      real phiprime(nphimx), phi(nphimx), phi0(nphimx),
     .     dphi



      common/sigcom/zi, eps0, v0i, omgrf

* ------------------------------
* storage for parallel scalapack
* ------------------------------
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, M_, N_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_, LLD_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

      integer            p_amat_dim, p_brhs_dim, p_ipiv_dim
      integer            p_amat_size

      complex, allocatable :: p_amat(:)
      complex, allocatable :: p_brhs(:)
      integer, allocatable :: p_ipiv(:)

      integer icnc
      logical ismine

      integer numroc
      external numroc
      integer lld,nrow_local,ncol_local

      character*4 suffix

      integer mb,nb,myid,nproc,myrow,mycol
      integer icontxt
      integer lrindx,lcindx,rsrc,csrc,ipos
      integer desc_amat(dlen_), desc_brhs(dlen_)

      integer rsrc1, csrc1, irnc1, icnc1

!efd-begin

*       ---------------------------------------------------
*       local variables for transformation to real space
*       ---------------------------------------------------
        logical, parameter :: use_Btmp_all = .true.
        integer :: ib,jb
        integer :: Locp, Locq, irow, icol, nia,nja
        integer :: niastart, niaend
        integer :: num_org, num_new, niasize
        integer :: org_nrow, org_ncol, mmb, nnb
        integer :: new_nrow, new_ncol
        integer :: num_inside, mm, ni, nn, niu, miu

        integer, dimension(:), allocatable :: new_to_org
        integer, dimension(:,:), allocatable :: niabegin_all
        integer, dimension(:,:), allocatable :: isize_all
        integer, dimension(:,:,:), allocatable :: descBtmp_all
        integer, dimension(:,:,:), allocatable :: descbrhs_all
        integer, dimension(DLEN_) :: descBtmp, descbrhs

        complex, dimension(:), allocatable :: brhs_tmp,brhs2
        complex, dimension(:,:), allocatable :: brhs
        complex, dimension(:,:), allocatable :: Btmp
        complex, dimension(:), allocatable :: row,rowk

        complex, dimension(:,:), allocatable :: inv_xx,inv_xx_t
        complex, dimension(:,:), allocatable :: inv_yy,inv_yy_t
        complex, dimension(:,:), allocatable :: inv_zz,inv_zz_t


        integer, dimension(:,:,:), allocatable :: mask
        real :: psimask



      integer idebug
      parameter(idebug=0)
      integer ia, ja
      integer, allocatable, dimension(:) ::
     .         itable, jtable, ktable,
     .         ntable, mtable, nphitable
      logical isok
      integer undefined
      parameter(undefined=-987654321)

      integer indxg2p, indxl2g
      external indxg2p, indxl2g

*     ---------------
*     allocate arrays
*     ---------------

      allocate( xjx(nxmx, nymx, nphimx),
     .        xjy(nxmx, nymx, nphimx),
     .        xjz(nxmx, nymx, nphimx) )

      allocate( psi(nxmx, nymx, nphimx),
     .     thetap(nxmx, nymx, nphimx),
     .     rho(nxmx, nymx, nphimx),
     .     xkte(nxmx, nymx, nphimx),
     .     xkti(nxmx, nymx, nphimx),
     .     xkti2(nxmx, nymx, nphimx),
     .     xkti3(nxmx, nymx, nphimx),
     .     xnea(nxmx, nymx, nphimx),
     .     xn1a(nxmx, nymx, nphimx),
     .     xn2a(nxmx, nymx, nphimx),
     .     xn3a(nxmx, nymx, nphimx),
     .     omgce(nxmx, nymx, nphimx),
     .     omgci1(nxmx, nymx, nphimx),
     .     omgci2(nxmx, nymx, nphimx),
     .     omgci3(nxmx, nymx, nphimx),
     .     omgpe2(nxmx, nymx, nphimx),
     .     omgp12(nxmx, nymx, nphimx),
     .     omgp22(nxmx, nymx, nphimx),
     .     omgp32(nxmx, nymx, nphimx),
     .     bmod(nxmx, nymx, nphimx),
     .     bmod_mid(nxmx, nymx, nphimx) )


      allocate( xjpxe(nxmx, nymx, nphimx),
     .        xjpye(nxmx, nymx, nphimx),
     .        xjpze(nxmx, nymx, nphimx),
     .        xjpx1(nxmx, nymx, nphimx),
     .        xjpy1(nxmx, nymx, nphimx),
     .        xjpz1(nxmx, nymx, nphimx),
     .        xjpx2(nxmx, nymx, nphimx),
     .        xjpy2(nxmx, nymx, nphimx),
     .        xjpz2(nxmx, nymx, nphimx),
     .        xjpx3(nxmx, nymx, nphimx),
     .        xjpy3(nxmx, nymx, nphimx),
     .        xjpz3(nxmx, nymx, nphimx),
     .        xjpx(nxmx, nymx, nphimx),
     .        xjpy(nxmx, nymx, nphimx),
     .        xjpz(nxmx, nymx, nphimx) )

      allocate( pcre(nxmx, nymx, nphimx),
     .     pcim(nxmx, nymx, nphimx),
     .     redotj(nxmx, nymx, nphimx),
     .     redotje(nxmx, nymx, nphimx),
     .     redotj1(nxmx, nymx, nphimx),
     .     redotj2(nxmx, nymx, nphimx),
     .     redotj3(nxmx, nymx, nphimx),
     .     redotjt(nxmx, nymx, nphimx) )


      allocate( bxn(nxmx, nymx, nphimx),
     .     byn(nxmx, nymx, nphimx),
     .     bzn(nxmx, nymx, nphimx) )

      allocate( uxx(nxmx, nymx, nphimx), uxy(nxmx, nymx, nphimx),
     .     uxz(nxmx, nymx, nphimx),
     .     uyx(nxmx, nymx, nphimx), uyy(nxmx, nymx, nphimx),
     .     uyz(nxmx, nymx, nphimx),
     .     uzx(nxmx, nymx, nphimx), uzy(nxmx, nymx, nphimx),
     .     uzz(nxmx, nymx, nphimx) )

      allocate( dxuxx(nxmx, nymx, nphimx), dxuxy(nxmx, nymx, nphimx),
     .     dxuxz(nxmx, nymx, nphimx),
     .     dxuyx(nxmx, nymx, nphimx), dxuyy(nxmx, nymx, nphimx),
     .     dxuyz(nxmx, nymx, nphimx),
     .     dxuzx(nxmx, nymx, nphimx), dxuzy(nxmx, nymx, nphimx),
     .     dxuzz(nxmx, nymx, nphimx) )

      allocate( dyuxx(nxmx, nymx, nphimx), dyuxy(nxmx, nymx, nphimx),
     .     dyuxz(nxmx, nymx, nphimx),
     .     dyuyx(nxmx, nymx, nphimx), dyuyy(nxmx, nymx, nphimx),
     .     dyuyz(nxmx, nymx, nphimx),
     .     dyuzx(nxmx, nymx, nphimx), dyuzy(nxmx, nymx, nphimx),
     .     dyuzz(nxmx, nymx, nphimx) )


      allocate( dyyuxx(nxmx, nymx, nphimx), dyyuxy(nxmx, nymx, nphimx),
     .     dyyuxz(nxmx, nymx, nphimx),
     .     dyyuyx(nxmx, nymx, nphimx), dyyuyy(nxmx, nymx, nphimx),
     .     dyyuyz(nxmx, nymx, nphimx),
     .     dyyuzx(nxmx, nymx, nphimx), dyyuzy(nxmx, nymx, nphimx),
     .     dyyuzz(nxmx, nymx, nphimx) )

      allocate( dxyuxx(nxmx, nymx, nphimx), dxyuxy(nxmx, nymx, nphimx),
     .     dxyuxz(nxmx, nymx, nphimx),
     .     dxyuyx(nxmx, nymx, nphimx), dxyuyy(nxmx, nymx, nphimx),
     .     dxyuyz(nxmx, nymx, nphimx),
     .     dxyuzx(nxmx, nymx, nphimx), dxyuzy(nxmx, nymx, nphimx),
     .     dxyuzz(nxmx, nymx, nphimx) )

      allocate( dxxuxx(nxmx, nymx, nphimx), dxxuxy(nxmx, nymx, nphimx),
     .     dxxuxz(nxmx, nymx, nphimx),
     .     dxxuyx(nxmx, nymx, nphimx), dxxuyy(nxmx, nymx, nphimx),
     .     dxxuyz(nxmx, nymx, nphimx),
     .     dxxuzx(nxmx, nymx, nphimx), dxxuzy(nxmx, nymx, nphimx),
     .     dxxuzz(nxmx, nymx, nphimx) )

      allocate( dphiuxx(nxmx, nymx, nphimx), dphiuxy(nxmx, nymx,nphimx),
     .     dphiuxz(nxmx, nymx, nphimx),
     .     dphiuyx(nxmx, nymx, nphimx), dphiuyy(nxmx, nymx, nphimx),
     .     dphiuyz(nxmx, nymx, nphimx),
     .     dphiuzx(nxmx, nymx, nphimx), dphiuzy(nxmx, nymx, nphimx),
     .     dphiuzz(nxmx, nymx, nphimx) )


      allocate( d2phiuxx(nxmx, nymx, nphimx),d2phiuxy(nxmx,nymx,nphimx),
     .     d2phiuxz(nxmx, nymx, nphimx),
     .     d2phiuyx(nxmx, nymx, nphimx), d2phiuyy(nxmx, nymx, nphimx),
     .     d2phiuyz(nxmx, nymx, nphimx),
     .     d2phiuzx(nxmx, nymx, nphimx), d2phiuzy(nxmx, nymx, nphimx),
     .     d2phiuzz(nxmx, nymx, nphimx) )

      allocate( dxphiuxx(nxmx, nymx, nphimx),dxphiuxy(nxmx,nymx,nphimx),
     .     dxphiuxz(nxmx, nymx, nphimx),
     .     dxphiuyx(nxmx, nymx, nphimx), dxphiuyy(nxmx, nymx, nphimx),
     .     dxphiuyz(nxmx, nymx, nphimx),
     .     dxphiuzx(nxmx, nymx, nphimx), dxphiuzy(nxmx, nymx, nphimx),
     .     dxphiuzz(nxmx, nymx, nphimx) )

      allocate( dyphiuxx(nxmx, nymx, nphimx),dyphiuxy(nxmx,nymx,nphimx),
     .     dyphiuxz(nxmx, nymx, nphimx),
     .     dyphiuyx(nxmx, nymx, nphimx), dyphiuyy(nxmx, nymx, nphimx),
     .     dyphiuyz(nxmx, nymx, nphimx),
     .     dyphiuzx(nxmx, nymx, nphimx), dyphiuzy(nxmx, nymx, nphimx),
     .     dyphiuzz(nxmx, nymx, nphimx) )

      allocate( gradprlb(nxmx, nymx, nphimx) )


      allocate( ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2,lkdim1:lkdim2),
     .        ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2),
     .        ebk(nkdim1 : nkdim2, mkdim1 : mkdim2, lkdim1 : lkdim2) )

      allocate( ex(nxmx, nymx, nphimx),
     .        ey(nxmx, nymx, nphimx),
     .        ez(nxmx, nymx, nphimx) )

      allocate( ealpha(nxmx, nymx, nphimx),
     .        ebeta(nxmx, nymx, nphimx),
     .        eb(nxmx, nymx, nphimx) )

      allocate( ntilda_e(nxmx, nymx, nphimx) )


      allocate( xx(nkdim1 : nkdim2, 1 : nxmx),
     .        yy(mkdim1 : mkdim2, 1 : nymx),
     .        zz(lkdim1 : lkdim2, 1 : nphimx) )

      allocate( reson1(nxmx, nymx, nphimx),
     .     reson2(nxmx, nymx, nphimx) )

      allocate( wdote(nxmx, nymx, nphimx), wdoti1(nxmx, nymx, nphimx),
     .     wdoti2(nxmx, nymx, nphimx), wdoti3(nxmx, nymx, nphimx),
     .     wdot(nxmx, nymx, nphimx) )


      allocate( fx0e(nxmx, nymx, nphimx), fx0i1(nxmx, nymx, nphimx),
     .     fx0i2(nxmx, nymx, nphimx), fx0i3(nxmx, nymx, nphimx),
     .     fy0e(nxmx, nymx, nphimx), fy0i1(nxmx, nymx, nphimx),
     .     fy0i2(nxmx, nymx, nphimx), fy0i3(nxmx, nymx, nphimx),
     .     fz0e(nxmx, nymx, nphimx), fz0i1(nxmx, nymx, nphimx),
     .     fz0i2(nxmx, nymx, nphimx), fz0i3(nxmx, nymx, nphimx),
     .     fx0(nxmx, nymx, nphimx), fy0(nxmx, nymx, nphimx),
     .     fz0(nxmx, nymx, nphimx) )

* --------------------------
* setup parallel environment
* --------------------------

* -----------------------------------------
* some environment may not require mpi_init
* -----------------------------------------
      call mpi_init( info )

*     ------------------------
*     Read namelist input data
*     ------------------------
      open(unit=63, file='aorsa3d.in', status='old', form='formatted')
      rewind(63)
      read (63, aorsa3din)


*     ---------------------------------------------------------
*     if (mcap_number .gt. 1), set prfin=0.0 so that pscale=1.0
*     ---------------------------------------------------------
      ! if(mcap_number .gt. 1) prfin = 0.0

      if (myid.eq.0) then
         open(unit=15,file='out15',status='unknown',form='formatted')
         open(unit=28,file='out28',status='unknown',form='formatted')
         open(unit=38,file='out38',status='unknown',form='formatted')
         open(unit=39, file='out_edotj', status='unknown',
     .      form='formatted')
         open(unit=53, file='fields_realSpace', status='unknown',
     .        form='formatted')
         open(unit=54, file='fields_fourier', status='unknown',
     .        form='formatted')
         open(unit=69,file='movie_wdot',status='unknown',
     .        form='formatted')

         write(39, 309) mcap_number
         write(39, 309) nmodesx, nmodesy, nmodesphi

         write(69, 309) mcap_number
         write(69, 309) nmodesx, nmodesy, nmodesphi
      endif

      idiag = nmodesx / 2
      jdiag = nmodesy / 2
      kdiag = nmodesphi / 2

      nbessj = lmax + 2

      ! TODO: nnodephi used before defined ??
      if (kplot .eq. 0) kplot = nnodephi / 2 + 1


!       --------------------------
!       determine points in plasma
!       --------------------------
        psimask = psilim


      call blacs_pinfo( myid, nproc )
      if (nproc .lt. 1) then
         write(6,*) '** blacs_pinfo returns:myid,nproc',myid,nproc
         stop '** blacs not setup '
      endif
      call blacs_get(-1,0,icontxt)

*     --------------------
*     setup processor grid
*     --------------------
*     -------------------------------------------------------
*     1 x n processor grid would simplify pivoting
*     n x 1 processor grid would simplify eqnx,eqny,eqnz
*     sqrt(n) x sqrt(n) normally would give best performance
*     -------------------------------------------------------


      call blacs_gridinit( icontxt, 'column-major', nprow,npcol)
      call blacs_gridinfo( icontxt, nprow, npcol, myrow, mycol)

      time00 = second1(dummy)

*     ----------------------------------------------------------
*     simplifying assumption that  all processors can open
*     same files and perform read
*     no need to re-broadcast input data
*     ----------------------------------------------------------
      if (myid .eq. 0) then
          write(6,*) 'blacs started: nprow,npcol,nproc ',
     &            nprow,npcol,nproc
      endif


*     --------------------------------------------------
*     may need to open different files like 'out38.001'
*     for processor '001'
*     --------------------------------------------------
      suffix = '.000'
      suffix(4:4) = char(ichar('0') + mod( myid,    10))
      suffix(3:3) = char(ichar('0') + mod( myid/10, 10))
      suffix(2:2) = char(ichar('0') + mod( myid/100,10))


      nkx2 = nmodesx / 2
      nkx1 = - nmodesx / 2 + 1
      nnodex = nmodesx
      nnoderho = nnodex / 2

      nky2 = nmodesy / 2
      nky1 = - nmodesy / 2 + 1
      nnodey = nmodesy

      nphi2 = nmodesphi / 2
      nphi1 = - nmodesphi / 2 + 1
      nnodephi = nmodesphi

      ndkplot = nnodephi / 4

      jequat  = nnodey / 2
      icenter = nnodex / 2

      if (qavg0 .ne. 0.0) xiota0 = 1./qavg0

      q = 1.6e-19
      if(te0 .eq. 0.0)te0 = ti0

      t0e = te0 * q
      t0i = ti0 * q
      t0i2 = ti02 * q
      t0i3 = ti03 * q

      qhat = qavg0

      xme = 9.11e-31
      xmh = 1.67e-27
      xmi1 = amu1 * xmh
      xmi2 = amu2 * xmh
      xmi3 = amu3 * xmh
      qi1 = z1 * q
      qi2 = z2 * q
      qi3 = z3 * q
      qe = -q
      zi  = cmplx(0.0,1.0)
      one = cmplx(1.0,0.0)
      eps0 = 8.85e-12
      pi = 3.141592654
      two_pi = 2.0 * pi
      xlnlam = 20.0
      xmu0 = 1.26e-06
      clight = 1.0 / sqrt(eps0 * xmu0)

      theta_antr = theta_ant / 180. * pi


      omgrf = two_pi * xnurf
      xk0 = omgrf / clight
      v0i = sqrt(2.0 * t0i / xmi1)

      xkphi0 = nphiant / rt
      xnphi0 = xkphi0 / xk0
      vphase = omgrf / xkphi0
      vthe = sqrt(t0e / xme)
      wphase = 0.0
      if(vthe .ne. 0.0)wphase = vphase / vthe

      xkthrho = 0.2
      xkthdx = 1.0


      rnz = xkphi0 / xk0
      xn1 = xn0 / z1 * (1.0 - z2 * eta - z3 * eta3)
      eta1 = xn1 / xn0
      xn2 = xn0 * eta
      eta2 = xn2 / xn0
      xn3 = xn0 * eta3

      if(myid .eq. 0)then

         write(15, *) "AORSA3D:"
         write(15, *) "mcap_number = ", mcap_number
         write(15, *) "mcap = ", mcap

         write (15, *)
         write (15, *) "xkperp_cutoff = ", xkperp_cutoff
         write (15, *) "damping       = ", damping

         write (15, *)
         write (15, *) "nnodex  = ", nnodex
         write (15, *) "nnodey  = ", nnodey
         write (15, *) "nnodephi = ", nnodephi

         write (15, *)
         write (15, *) "nprow   = ", nprow
         write (15, *) "npcol   = ", npcol
         write (15, *)
         write (15, *) "prfin     = ", prfin
         write (15, *)
         write (15, *) "delta0  = ", delta0
         write (15, *) "xnuomg  = ", xnuomg
         write (15, *)
         write (15, *) "lmax    = ", lmax
         write (15, *) "zeffcd  = ", zeffcd
         write (15, *)

         write(15, 162)
         write(15, 7013)nmodesx, nmodesy, nmodesphi
         write(15, 7113)nwdot
         write(15, 7213)nnodecx, nnodecy

         write(15, 1815)t0e
         write(15, 1821)t0i

         write(15, 7217)qavg0
         write(15, 7014)lmax
         write(15, 7015)ibessel
         write(15, 7115)nzfun
         write(15, 7116)xnuomg
         write(15, 7217)qavg0
         write(15, 1321)wphase
         write(15, 1323)vthe
         write(15, 1021)rnz
         write(15, 1812)rt
         write(15, 1822)aplasm
         write(15, 1823)rant
         write(15, 1809)b0
         write(15, 6813)xn0
         write(15, 1813)xn1
         write(15, 1814)xn2
         write(15, 1834)xn3

         write(15, 1012)omgrf
         write(15, 1009)xnurf
         write(15, 1321)wphase
         write(15, 1322)vphase
         write(15, 1323)vthe
         write(15, 1714)xk0
         write(15, 1016)xnuead
         write(15, 1017)xnu1ad
         write(15, 1018)xnu2ad
         write(15, 1020)nnodex, nnodey

         write(15, *) "iprofile = ", iprofile
         write(15, *) "psilim = ", psilim
      end if

      telimj   = telim  * q
      tilimj   = tilim  * q
      ti2limj  = ti2lim * q
      ti3limj  = ti3lim * q


      t1 = second1(dummy)

*     ------------------------------------
*     Stellarator flux surfaces from VMEC
*     ------------------------------------
      if (igeom .eq. 5) then
         call vmec_setup(myid, nproc, iwout, wout,
     .       nxmx, nymx, nphimx, nnodex, nnodey, nnodephi,
     .       capr, y, phi, rwleft, rwright, ytop, ybottom, rt, b0,
     .       bxn, byn, bzn, bmod, psi, rho, thetap, nfp)

         if (wout_bscale .lt. 0.0) then
            do k = 1, nnodephi
               do i = 1, nnodex
                  do j = 1, nnodey
                     bxn(i, j, k) = -bxn(i, j, k)
                     byn(i, j, k) = -byn(i, j, k)
                     bzn(i, j, k) = -bzn(i, j, k)
                  enddo
               enddo
            enddo
         endif
         if (abs(wout_bscale) .ne. 1.0) then
            do k = 1, nnodephi
               do i = 1, nnodex
                  do j = 1, nnodey
                     bmod(i, j, k) = bmod(i, j, k) * abs(wout_bscale)
                  enddo
               enddo
            enddo

            b0 = b0 * abs(wout_bscale)
         endif

         xkphi0 = nphiant / rt
         xnphi0 = xkphi0 / xk0
         vphase = omgrf / xkphi0

         wphase = 0.0
         if(vthe .ne. 0.0)wphase = vphase / vthe

         rnz = xkphi0 / xk0

         if(myid .eq. 0)then
            write(15, 1812)rt
            write(15, 1809)b0
         endif

      end if

      time = second1(dummy) - t1
      tmin = time / 60.
      if(myid .eq. 0)then
         write(15, 386) tmin
      end if
  386 format('time to call vmec_setup = ',f9.3, 4h min)



*     --------------------------------
*     Distribution function
*     --------------------------------

!     -----------------
!     allocate arrays
!     -----------------
      allocate( UPERP(nuper) )
      allocate( UPARA(nupar) )
      allocate( VPERP(nuper) )
      allocate( VPARA(nupar) )
      allocate( UPERP_work(nuper) )
      allocate( UPARA_work(nupar) )


      allocate( f(nuper, nupar) )

      allocate( dfdupere(nuper, nupar) )
      allocate( dfdupare(nuper, nupar) )

      allocate( dfduper1(nuper, nupar) )
      allocate( dfdupar1(nuper, nupar) )

      allocate( dfduper2(nuper, nupar) )
      allocate( dfdupar2(nuper, nupar) )

      allocate( dfduper3(nuper, nupar) )
      allocate( dfdupar3(nuper, nupar) )


      UPERP = 0.0
      UPARA = 0.0
      UPERP_work = 0.0
      UPARA_work = 0.0

      f = 0.0

      dfdupere = 0.0
      dfdupare = 0.0

      dfduper1 = 0.0
      dfdupar1 = 0.0

      dfduper2 = 0.0
      dfdupar2 = 0.0

      dfduper3 = 0.0
      dfdupar3 = 0.0


      if(ndisti1   .ne. 0 .or.
     .   ndisti2   .ne. 0 .or.
     .   ndisti3   .ne. 0 .or.
     .   ndiste    .ne. 0)   then  ! ndist!=0


         allocate( fe_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( dfe_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( dfe_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f1_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df1_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df1_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f2_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df2_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df2_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f3_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df3_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df3_cql_uprl(nuper, nupar, n_psi_dim) )



         fe_cql_cart = 0.0
         dfe_cql_uprp = 0.0
         dfe_cql_uprl = 0.0

         f1_cql_cart = 0.0
         df1_cql_uprp = 0.0
         df1_cql_uprl = 0.0

         f2_cql_cart = 0.0
         df2_cql_uprp = 0.0
         df2_cql_uprl = 0.0

         f3_cql_cart = 0.0
         df3_cql_uprp = 0.0
         df3_cql_uprl = 0.0

         call blacs_barrier(icontxt, 'All')

         ! Common for three cases
         n_psi = nnoderho
         UminPara = -1.0
         UmaxPara = 1.0
         do n = 1, n_psi
            rho_a(n) = real(n-1) / real(n_psi-1)
         end do
         do n = 1, NUPER
            UPERP(n) = real(n-1) / real(NUPER-1)
         end do
         do m = 1, NUPAR
            UPARA(m) = UminPara
     .          + (UmaxPara - UminPara) * (real(m-1)/real(NUPAR-1))
         end do

         if (ndiste .eq. 1) then
            vce_mks = 3.0 * sqrt(2.0 * te0 * q / xme)
            do n = 1, NUPER
               do m = 1, NUPAR
                  do l = 1, n_psi
                     fe_cql_cart(n, m, l)
     .                   = exp(-(uperp(n)**2 + upara(m)**2) * (3.0**2))
     .                   * ((3.0 / sqrt(pi))**3)
                     dfe_cql_uprp(n, m, l)
     .                   = fe_cql_cart(n, m, l)
     .                   * (3.0**2) * (-2.0 * uperp(n))
                     dfe_cql_uprl(n, m, l)
     .                   = fe_cql_cart(n, m, l)
     .                   * (3.0**2) * (-2.0 * upara(m))
                  end do
               end do
            end do

         end if

         call blacs_barrier(icontxt, 'All')

         if (ndisti1 .eq. 1) then
            vc1_mks = 3.0 * sqrt(2.0 * ti0 * q / xmi1)
            do n = 1, NUPER
               do m = 1, NUPAR
                  do l = 1, n_psi
                     f1_cql_cart(n, m, l)
     .                   = exp(-(uperp(n)**2 + upara(m)**2) * (3.0**2))
     .                   * ((3.0 / sqrt(pi))**3)
                     df1_cql_uprp(n, m, l)
     .                   = f1_cql_cart(n, m, l)
     .                   * (3.0**2) * (-2.0 * uperp(n))
                     df1_cql_uprl(n, m, l)
     .                   = f1_cql_cart(n, m, l)
     .                   * (3.0**2) * (-2.0 * upara(m))
                  end do
               end do
            end do
         end if

         call blacs_barrier(icontxt, 'All')

         if (ndisti2 .eq. 1) then
            vc2_mks = 3.0 * sqrt(2.0 * ti02 * q / xmi2)
            do n = 1, NUPER
               do m = 1, NUPAR
                  do l = 1, n_psi
                     f2_cql_cart(n, m, l)
     .                   = exp(-(uperp(n)**2 + upara(m)**2) * (3.0**2))
     .                   * ((3.0 / sqrt(pi))**3)
                     df2_cql_uprp(n, m, l)
     .                   = f2_cql_cart(n, m, l)
     .                   * (3.0**2) * (-2.0 * uperp(n))
                     df2_cql_uprl(n, m, l)
     .                   = f2_cql_cart(n, m, l)
     .                   * (3.0**2) * (-2.0 * upara(m))
                  end do
               end do
            end do

         end if



         if(myid .eq. 0)then
            WRITE (15,*)
            WRITE (15,*) "nuper = ",  nuper
            WRITE (15,*) "nupar = ",  nupar
            WRITE (15,*) "n_psi = ",  n_psi

            WRITE (15,*)
            WRITE (15,*) "vce_mks = ", vce_mks
            WRITE (15,*) "vc1_mks = ", vc1_mks
            WRITE (15,*) "vc2_mks = ", vc2_mks
            WRITE (15,*) "vc3_mks = ", vc3_mks

            write(15, *)
            write(15, *) "rho/a(i_psi)"
            write(15, 310) (rho_a(n), n = 1, n_psi)

            write(15, *)
            write(15, *) "uperp ="
            write(15, 310) (uperp(n), n = 1, nuper)

            write(15, *)
            write(15, *) "upara ="
            write(15, 310) (upara(n), n = 1, nupar)
         end if

      end if  ! ndist!=0


      if(myid .eq. 0)write(53, 309) nnodex, nnodey, nnodephi

      xwleft = rwleft - rt
      xwright = rwright - rt
      xmax = rwright - rwleft
      ymax = ytop - ybottom

*     ---------------------------------------
*     Define x mesh: x(i), xprime(i), capr(i)
*     ---------------------------------------
c--   xprime: 0 to xmax
c--   x(i) : -xmax / 2.0   to   xmax / 2.0
      dx = xmax / nnodex
      do i = 1, nnodex
         xprime(i) = (i-1) * dx  + dx / 2.0
c--   Note: the code gives slightly smoother results with dx/2.0 added
         x(i) = xprime(i) + xwleft
         capr(i) = rt + x(i)
         if(myid .eq. 0)write(53, 1312)i, capr(i)
      end do



      dxc = xmax / nnodecx
      do i = 1, nnodecx
         xprimec(i) = (i-1) * dxc + dxc / 2.0
         xcourse(i) = xprimec(i) + xwleft
         caprc(i) = rt + xcourse(i)

      end do

*     ---------------------------------------
*     Define y mesh: y(j), yprime(j)
*     ---------------------------------------
c--   yprime: 0 to ymax
c--   y(j) : -ymax / 2.0   to   ymax / 2.0
      dy = ymax / nnodey
      do j = 1, nnodey
         yprime(j) = (j-1) * dy  + dy / 2.0
c--      Note: the code gives slightly smoother results with dy/2.0 added
         y(j) = yprime(j) + ybottom
         if(myid .eq. 0)write(53, 1312)j, y(j)
      end do

      dyc = ymax / nnodecy
      do j = 1, nnodecy
         yprimec(j) = (j-1) * dyc + dyc / 2.0
         ycourse(j) = yprimec(j) + ybottom
      end do


      phimax = two_pi / nfp

      if(myid .eq. 0)then
         write(15, 920)
         write(15, 7164)
         write(15, 920)
      end if

*     ---------------------------------------
*     Define phi mesh: phi(k), phiprime(k)
*     ---------------------------------------
c--   phiprime: 0 to phimax
c--   phi(k) : -phimax / 2.0  to  phimax / 2.0
c--   phi0(k): 0 to phimax / 2.0 and 0 to -phimax / 2.0
      dphi = phimax / nnodephi


      do k = 1, nnodephi
         phiprime(k) = (k-1) * dphi
c     1       + dphi / 2.0
c--      Note: the code gives worse results with dz/2.0 added
         phi(k) = phiprime(k) - phimax / 2.0

         if(phi(k) .le. 0.0) then
             phi0(k) = phi(k) + phimax / 2.0
         else
            phi0(k) = phi(k) - phimax / 2.0
         endif

         if(myid .eq. 0)write(15, 1312)k, phi(k), phiprime(k), phi0(k)
         if(myid .eq. 0)write(53, 1312)k, phiprime(k)
      end do

      do i = 1, nnodex
         do j = 1, nnodey
            if(x(i) .ne. 0.0 .or. y(j) .ne. 0.0)then
               theta0(i,j) = atan2(y(j), x(i))
               if(theta0(i,j) .ge. 0.0) then
                   theta(i,j) = theta0(i,j)
               else
                   theta(i,j) = theta0(i,j) + two_pi
               endif
            end if
         end do
      end do


c--   Define rho mesh: rhon(n)
c--   rhon: 0 to rhomax
      rhomax = 1.0
      drho = rhomax / (nnoderho - 1)
      do n = 1, nnoderho
         rhon(n) = (n-1) * drho
c         write(6, 1312)n, rhon(n)
      end do


      xant = rant - rt
      iant=int((rant - rwleft) / dx) + 1
      if (rant .ne. 0.0) psiant = psi(iant, jdiag, kdiag)

c--note curden is in Amps per meter of toroidal length (2.*pi*rt).
      xjantx = curdnx / dx
      xjanty = curdny / dx
      xjantz = curdnz / dx
      xjant=sqrt(xjantx**2 + xjanty**2 + xjantz**2)

      rplasm = rt + aplasm
      rlim   = rt + alim



*     --------------------------------
*     Tokamak (Solovev) flux surfaces:
*     --------------------------------
      if(igeom .eq. 1)then
         np = mhel
         nfp = mhel
         q07qa = 0.0
         rhonorm = xwright * (1.0 + xwright / 2. / rt)
         if (myid .eq. 0) then
            write(15, 1022)rhonorm
         end if


         psiplasm = xiota0 * b0 / 2.0 *
     .       ( (rplasm**2 - rt**2)**2 / 4. / rt**2 )
         rhoplasm = 1.0/ rhonorm * sqrt(2.0 * psiplasm /(xiota0 * b0))


         psilim = xiota0 * b0 / 2.0 *
     .       ( (rlim**2 - rt**2)**2 / 4. / rt**2 )
         rholim = 1.0 / rhonorm * sqrt(2.0 * psilim / (xiota0 * b0))


         psiant = b0 * xiota0 / 2. *
     .       ( (rant**2 - rt**2)**2 / 4. / rt**2 )
         rhoant = 1.0 / rhonorm * sqrt(2.0 * psiant / (xiota0 * b0))


         psimax = b0 * xiota0 / 2. *
     1             ((capr(nnodex) * y(nnodey) / rt / ekappa)**2
     1             + (capr(nnodex)**2 - rt**2)**2 / 4. / rt**2 )
         rhomax = 1.0 / rhonorm *
     1                     sqrt(2.0 * psimax / (xiota0 * b0))

         do i = 1, nnodex
            do j = 1, nnodey
               do k = 1, nnodephi
                  psi(i, j, k) = b0 * xiota0 / 2. *
     1                ((capr(i) * y(j) / rt / ekappa)**2
     1                + (capr(i)**2 - rt**2)**2 / 4. / rt**2 )

                  denom = ekappa / 2. / capr(i) *(capr(i)**2 - rt**2)
                  if(y(j) .ne. 0.0 .or. denom .ne. 0.0) then
                        theta0(i,j) = atan2(y(j), denom)
                     if(theta0(i,j) .ge. 0.0) then
                         theta(i,j) = theta0(i,j)
                     else
                         theta(i,j) = theta0(i,j) + two_pi
                     endif
                  end if

                  rho(i,j,k) = 1.0 / rhonorm *
     1                     sqrt(2.0 * psi(i,j,k) / (xiota0 * b0))
                  if(rho(i,j,k) .eq. 0.0)rho(i,j,k) = 1.0e-08

                  bx = -b0 * xiota0 * capr(i) * y(j)
     1                / rt**2 / ekappa**2
                  by = xiota0 * b0 / 2.0 *
     1                   ( 2.0 * y(j)**2 / rt**2 / ekappa**2
     1                    + capr(i)**2 / rt**2 - 1.0)



                  gaussian =  exp(-rho(i,j,k)**2 / rhoplasm**2)
                  frho   = q07qa + (1.0 - q07qa) * gaussian**alphan


                  bx = bx * frho
                  by = by * frho


                  bz = b0 * rt / capr(i)
                  bmod(i, j, k) = sqrt(bx**2 + by**2 + bz**2)


                  bxn(i,j,k) = bx / bmod(i,j,k)
                  byn(i,j,k) = by / bmod(i,j,k)
                  bzn(i,j,k) = bz / bmod(i,j,k)
               end do
            end do
         end do


      end if


*     --------------------------
*     Stellarator flux surfaces
*     --------------------------
      if (igeom .eq. 2) then
         np = mhel
         nfp = mhel
         xlp = lhel * two_pi * rt / mhel
         xkh0 = two_pi / xlp
         capa = 0.5 * (ekappa**2 - 1.0) / (ekappa**2 + 1.0)
         epsl = - 2.0 * b0 * capa


         if(myid .eq. 0)then
            write(6, 1023)epsl
         end if


         btor = b0 * rt /(rt + aplasm)

         arg = mhel * aplasm / rt
         call besic(arg, lhel + 1, bl, ier)
         do l = 0, lhel
            xil(l) = bl(l+1)
         end do
         do l = 0, lhel
            if(l .eq. 0)xilp(0) = xil(1)
            if(l .ne. 0)xilp(l) = xil(l-1) - l / arg * xil(l)
         end do

         if(iexpnd .eq. 0)xilp(lhel) = arg / 4.0

         psiplasm = btor * xkh0 * aplasm**2 / 2.
     1                     - aplasm * epsl * xilp(lhel)


         btor = b0 * rt /(rt + alim)

         arg = mhel * alim / rt
         call besic(arg, lhel + 1, bl, ier)
         do l = 0, lhel
            xil(l) = bl(l+1)
         end do
         do l = 0, lhel
            if(l .eq. 0)xilp(0) = xil(1)
            if(l .ne. 0)xilp(l) = xil(l-1) - l / arg * xil(l)
         end do

         if(iexpnd .eq. 0)xilp(lhel) = arg / 4.0


         psilim = btor * xkh0 * alim**2 / 2.
     1                     - alim * epsl * xilp(lhel)



         btor = b0 * rt /(rt + xant)

         arg = mhel * xant / rt
         call besic(arg, lhel + 1, bl, ier)
         do l = 0, lhel
            xil(l) = bl(l+1)
         end do
         do l = 0, lhel
            if(l .eq. 0)xilp(0) = xil(1)
            if(l .ne. 0)xilp(l) = xil(l-1) - l / arg * xil(l)
         end do

         if(iexpnd .eq. 0)xilp(lhel) = arg / 4.0

         psiant = btor * xkh0 * xant**2 / 2.
     1                     - xant * epsl * xilp(lhel)



         btor = b0 * rt /(rt + xwright)

         arg = mhel * xwright / rt
         call besic(arg, lhel + 1, bl, ier)
         do l = 0, lhel
            xil(l) = bl(l+1)
         end do
         do l = 0, lhel
            if(l .eq. 0)xilp(0) = xil(1)
            if(l .ne. 0)xilp(l) = xil(l-1) - l / arg * xil(l)
         end do

         if(iexpnd .eq. 0)xilp(lhel) = arg / 4.0

         psiright = btor * xkh0 * xwright**2 / 2.
     1                     - xwright * epsl * xilp(lhel)


         rhonorm = sqrt(2.0 * psiright / (xkh0 * b0) )

         if(myid .eq. 0)then
            write(15, 1022)rhonorm
         end if



         do i = 1, nnodex
            do j = 1, nnodey
               r = sqrt (x(i)**2 + y(j)**2)

               btor = b0 * rt / capr(i)

               arg = lhel * xkh0 * r

               call besic(arg, lhel + 1, bl, ier)

               do l = 0, lhel
                  xil(l) = bl(l+1)
               end do

               do l = 0, lhel
                  if(l .eq. 0)xilp(0) = xil(1)
                  if(l .ne. 0)xilp(l) = xil(l-1) - l / arg * xil(l)
               end do

               if(iexpnd .eq. 0)then
                  xil(lhel) = arg**2 / 8.0
                  xilp(lhel) = arg / 4.0
               end if

               if(x(i) .ne. 0.0 .or. y(j) .ne. 0.0)then
                  theta0(i,j) = atan2(y(j), x(i))
                  if(theta0(i,j) .ge. 0.0) then
                      theta(i,j) = theta0(i,j)
                  else
                      theta(i,j) = theta0(i,j) + two_pi
                  endif
               end if
               costh = cos(theta(i,j))
               sinth = sin(theta(i,j))


               do k = 1, nnodephi

                  angle = lhel * theta(i,j) + mhel * phiprime(k)
                  cosang = cos(angle)
                  sinang = sin(angle)

                  psi(i, j, k) = btor * xkh0 * r**2 / 2.
     .                           - r * epsl * xilp(lhel) * cosang

                  rho(i, j, k) = sqrt(2.0 * psi(i,j,k) / (xkh0 * b0) )
     .               / rhonorm


                  br = lhel * epsl * xilp(lhel) * sinang
                  bth = 1. / (xkh0 * r) * lhel * epsl * xil(lhel)
     .                    * cosang
                  bphi = btor - lhel * epsl * xil(lhel) * cosang

                  bmod(i, j, k) = sqrt(br**2 + bth**2 + bphi**2)


                  bx = br * costh - bth * sinth
                  by = br * sinth + bth * costh
                  bz = bphi


                  bxn(i,j,k) = bx / bmod(i,j,k)
                  byn(i,j,k) = by / bmod(i,j,k)
                  bzn(i,j,k) = bz / bmod(i,j,k)

               end do
            end do
         end do


      end if


*     ------------------------------------
*     Stellarator (expanded) flux surfaces
*     ------------------------------------
      if (igeom .eq. 3 .or. igeom .eq. 4) then

         np = mhel
         nfp = mhel
         xlp = lhel * two_pi * rt / mhel
         xkh0 = two_pi / xlp
         capa = 0.5 * (ekappa**2 - 1.0) / (ekappa**2 + 1.0)


         if(myid .eq. 0)then
            write(15, 91023)capa
         end if


         btor = b0 * rt /(rt + alim)
         psi_lim = xkh0 * btor / 2.0 * alim**2 * (1.0 + 2.0 * capa)


         do i = 1, nnodex
            btor = b0 * rt / capr(i)
            xkh = mhel / (lhel * capr(i))
            do j = 1, nnodey
               do k = 1, nnodephi

                  z = rt * phiprime(k)
                  if(igeom .eq. 3) hz =   xkh0 * z
                  if(igeom .eq. 4) hz = - xkh0 * z
                  coshz = cos(hz)
                  sinhz = sin(hz)
                  capx = x(i) * coshz - y(j) * sinhz
                  capy = x(i) * sinhz + y(j) * coshz
                  r2 = capx**2 + capy**2


                  psi1 = xkh0 / 2.0 *  btor *
     .                          (r2 + 2.0 * capa * (capx**2 - capy**2))
                  psi(i, j, k) = psi1 / psi_lim
                  rho(i, j, k) = sqrt(psi(i, j, k))


                  denom = 1.0 + xkh0**2 * r2

                  bcapx = -2.0 * capa * btor / denom
     .                 * xkh0 * capy * (1.0 + 2.0 * xkh0**2 * capx**2)

                  bcapy = -2.0 * capa * btor / denom
     .                 * xkh0 * capx * (1.0 + 2.0 * xkh0**2 * capy**2)

                  bcapz = btor + 2.0 * capa * btor / denom
     .                 * xkh0**2 *(capx**2 - capy**2)

                  bx =   bcapx * coshz + bcapy * sinhz
                  by = - bcapx * sinhz + bcapy * coshz
                  bz =   bcapz


                  bmod(i, j, k) = sqrt(bx**2 + by**2 + bz**2)


                  bxn(i,j,k) = bx / bmod(i,j,k)
                  byn(i,j,k) = by / bmod(i,j,k)
                  bzn(i,j,k) = bz / bmod(i,j,k)

               end do
            end do
         end do


      end if

      rhomin = 1.0
      k0 = nnodephi / 2
      do i = 1, nnodex
         do j = 1, nnodey
            if (rho(i,j,k0) .lt. rhomin) then
               i0 = i
               j0 = j
               rhomin = rho(i,j,k0)
            end if
         end do
      end do

*     -----------------------
*     Calculate bmod_mid(i,j,k)
*     -----------------------
      bmod_min = 0.0
      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi
               n = int(rho(i,j,k) / drho) + 1
               if(n .ge. 1 .and. n .le. nnoderho)then
                  if (bmod_min(n) .eq. 0.0
     .                .or. bmod_min(n) .lt. bmod(i,j,k)) then
                     bmod_min(n) = bmod(i,j,k)
                  end if
               end if
            end do
         end do
      end do
      bmod_mid = bmod
      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi
               n = int(rho(i,j,k) / drho) + 1
               if(n .ge. 1 .and. n .le. nnoderho)then
                  bmod_mid(i,j,k) = bmod_min(n)
               end if
            end do
         end do
      end do


*     ------------------------------
*     Calculate the antenna current
*     ------------------------------

      dpsiant = dpsiant0
      yant_max = antlen / 2.0

      phiant = pi / nfp

      if(myid .eq. 0)then
         write(15, 920)
         write(15, 8165)
         write(15, 920)
      end if

      if (rant .ne. 0.0) then
         rant1 = rant
      else
         do i = nnodex, 1, -1
            if (psi(i, jdiag, kdiag) .lt. psiant) then
               rant1 = capr(i)
               exit
            end if
         end do
      end if
      if (myid .eq. 0) then
         write(6, *) "rant1 =", rant1
      end if

      do k = 1, nnodephi
         do i = 1, nstrap
            zstrap = (i - 1.0 - (nstrap - 1) / 2.0) * wd
            zn = (phiprime(k) - phiant) * rant1
            if (abs(zn - zstrap) .lt. (xlt / 2.0))
     .           shapephi(k) = amplt(i)
     .           * exp(zi * phase_deg(i) * pi / 180.)
         end do

         if(myid .eq. 0)write(15, 1312) k, shapephi(k)
      end do


        allocate( mask(1:nnodex,1:nnodey,1:nnodephi) )


      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi

               gausspsi = exp(-(psi(i,j,k) - psiant)**2 / dpsiant**2)

               shapey = 0.0
               deltay = y(j) - yant

               if(capr(i) .gt. rt) then
!                1) antenna current is cos(ky * y)  (DEFAULT)
                 if(i_antenna.eq.1 .and. abs(deltay).lt.yant_max) then
                    shapey = cos(xk0 * antlc * deltay)
!                0) antenna current is Gaussian
                 elseif (i_antenna .eq. 0) then
                    shapey = exp(-deltay**2 / yant_max**2)
                 endif
               endif

               gausspsi = gausspsi * shapey * shapephi(k)
               xjx(i,j,k) = xjantx * gausspsi
               xjy(i,j,k) = xjanty * gausspsi
               xjz(i,j,k) = xjantz * gausspsi

            end do
         end do
      end do


*    ----------------
*    plasma profiles:
*    ----------------
      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi

               if (iprofile .eq. 1) then
                  gaussian =  exp(-psi(i,j,k) / psiplasm)

                  xnea(i, j, k) = xnlim
     .                 + (xn0 - xnlim) * gaussian**alphan
                  xkte(i,j,k) = telimj
     .                 + (t0e - telimj) * gaussian**alphate
                  xkti(i,j,k) = tilimj
     .                 + (t0i - tilimj) * gaussian**alphati
                  xkti2(i,j,k) = ti2limj
     .                 + (t0i2 - ti2limj) * gaussian**alphati2
                  xkti3(i,j,k) = ti3limj
     .                 + (t0i3 - ti3limj) * gaussian**alphati3

               else if (iprofile .eq. 3) then
                  if (psi(i,j,k) .lt. psiplasm) then
                     xnea(i,j,k) = (psi(i,j,k) / psiplasm)**betan
                     xnea(i,j,k) = (1.0 - xnea(i,j,k))**alphan
                     xnea(i,j,k) = xnlim + (xn0 - xnlim) * xnea(i,j,k)

                     xkte(i,j,k) = (psi(i,j,k) / psiplasm)**betate
                     xkte(i,j,k) = (1.0 - xkte(i,j,k))**alphate
                     xkte(i,j,k) = telimj + (t0e - telimj) * xkte(i,j,k)
                     xkti(i,j,k) = (psi(i,j,k) / psiplasm)**betati
                     xkti(i,j,k) = (1.0 - xkti(i,j,k))**alphati
                     xkti(i,j,k) = tilimj + (t0i - tilimj) * xkti(i,j,k)
                     xkti2(i,j,k) = (psi(i,j,k) / psiplasm)**betati2
                     xkti2(i,j,k) = (1.0 - xkti2(i,j,k))**alphati2
                     xkti2(i,j,k) = ti2limj
     .                    + (t0i2 - ti2limj) * xkti2(i,j,k)
                     xkti3(i,j,k) = (psi(i,j,k) / psiplasm)**betati3
                     xkti3(i,j,k) = (1.0 - xkti3(i,j,k))**alphati3
                     xkti3(i,j,k) = ti3limj
     .                    + (t0i3 - ti3limj) * xkti3(i,j,k)
                  else
                     xnea(i,j,k) = xnlim
                     xkte(i,j,k) = telimj
                     xkti(i,j,k) = tilimj
                     xkti2(i,j,k) = ti2limj
                     xkti3(i,j,k) = ti3limj
                  end if
               end if

               if (iprofile .ne. 5) then
                  xn2a(i, j, k) = xnea(i, j, k) * eta2
                  xn3a(i, j, k) = xnea(i, j, k) * eta3
                  xn1a(i, j, k) = (xnea(i, j, k)
     .                 - xn2a(i, j, k) * z2
     .                 - xn3a(i, j, k) * z3) / z1
               end if

*               ------------------------------------
*               determine which points are in plasma
*               ------------------------------------
               mask(i,j,k) = 1
               if (psi(i,j,k).gt.psimask) mask(i,j,k) = 0

            end do
         end do
      end do

*     ---------------------------------------
*     read profile data from the plasma state:
*     ---------------------------------------
      if (iprofile .eq. 5) then
         rewind (63)
         read(63, state)

         s_t_s = q * s_t_s * 1.0e3    !input is kev TEMPS IN JOULES

         xnea = 1.0e-10
         call flux_to_rzphi(nnodex,nnodey,nnodephi,s_n_s(1:s_nrho_n,0),
     .        xnea(1:nnodex,1:nnodey,1:nnodephi),
     .        s_rho_n_grid(1:s_nrho_n),s_nrho_n,
     .        rho(1:nnodex,1:nnodey,1:nnodephi))

         xkte = 1.0e-10 * q
         call flux_to_rzphi(nnodex,nnodey,nnodephi,s_t_s(1:s_nrho_n,0),
     .        xkte(1:nnodex,1:nnodey,1:nnodephi),
     .        s_rho_n_grid(1:s_nrho_n),s_nrho_n,
     .        rho(1:nnodex,1:nnodey,1:nnodephi))

         xn1a = 1.0e-10
         call flux_to_rzphi(nnodex,nnodey,nnodephi,s_n_s(1:s_nrho_n,1),
     .        xn1a(1:nnodex,1:nnodey,1:nnodephi),
     .        s_rho_n_grid(1:s_nrho_n),s_nrho_n,
     .        rho(1:nnodex,1:nnodey,1:nnodephi))

         xkti = 1.0e-10 * q
         call flux_to_rzphi(nnodex,nnodey,nnodephi,s_t_s(1:s_nrho_n,1),
     .        xkti(1:nnodex,1:nnodey,1:nnodephi),
     .        s_rho_n_grid(1:s_nrho_n),s_nrho_n,
     .        rho(1:nnodex,1:nnodey,1:nnodephi))

         xn2a = 1.0e-10
         call flux_to_rzphi(nnodex,nnodey,nnodephi,s_n_s(1:s_nrho_n,2),
     .        xn2a(1:nnodex,1:nnodey,1:nnodephi),
     .        s_rho_n_grid(1:s_nrho_n),s_nrho_n,
     .        rho(1:nnodex,1:nnodey,1:nnodephi))

         xkti2 = 1.0e-10 * q
         call flux_to_rzphi(nnodex,nnodey,nnodephi,s_t_s(1:s_nrho_n,2),
     .        xkti2(1:nnodex,1:nnodey,1:nnodephi),
     .        s_rho_n_grid(1:s_nrho_n),s_nrho_n,
     .        rho(1:nnodex,1:nnodey,1:nnodephi))

         xn3a = 1.0e-10
         call flux_to_rzphi(nnodex,nnodey,nnodephi,s_n_s(1:s_nrho_n,3),
     .        xn3a(1:nnodex,1:nnodey,1:nnodephi),
     .        s_rho_n_grid(1:s_nrho_n),s_nrho_n,
     .        rho(1:nnodex,1:nnodey,1:nnodephi))

         xkti3 = 1.0e-10 * q
         call flux_to_rzphi(nnodex,nnodey,nnodephi,s_t_s(1:s_nrho_n,3),
     .        xkti3(1:nnodex,1:nnodey,1:nnodephi),
     .        s_rho_n_grid(1:s_nrho_n),s_nrho_n,
     .        rho(1:nnodex,1:nnodey,1:nnodephi))


*        ---------------------------------
*        Minority ion (species 2) profiles:
*        ---------------------------------
*        for now leave the species 2 minority option open
*        if eta comes in as zero or small, use namelist values
*        else create a minority

         if(eta .ge. 1.0e-6) then
            do i = 1, nnodex
               do j = 1, nnodey
                  do k = 1, nnodephi
                     xn2a(i,j,k) = eta * xnea(i,j,k)

*                 ----------------
*                 Cold minority H:
*                 ----------------
*                  xkti2(i,j) = xkti(i,j)

*                 ----------------
*                 Hot minority H:
*                 ----------------
c                 xkti2(i,j) = xkti3(i,j)

                  end do
               end do
            end do
         end if

*        --------------------------------------------------
*        Adjust ni1a for charge neutrality on all processors:
*        --------------------------------------------------
         do i = 1, nnodex
            do j = 1, nnodey
               do k = 1, nnodephi
                  xn1a(i, j, k) = (xnea(i, j, k) - z2 * xn2a(i, j, k)
     .                                  - z3 * xn3a(i, j, k)   ) / z1

               if (xn1a(i, j, k) .le. 0.0) xn1a(i, j, k) = 1.0e-10
               end do
            end do
         end do



         xn0 = xnea(i0, j0, k0)
         xn1 = xn1a(i0, j0, k0)
         xn2 = xn2a(i0, j0, k0)
         xn3 = xn3a(i0, j0, k0)
         eta1 = xn1 / xn0
         eta2 = xn2 / xn0
         eta3 = xn3 / xn0

         if(myid .eq. 0)then
            write(15, 6834)eta1
            write(15, 6835)eta2
            write(15, 6836)eta
         end if
      end if

*     ----------------------------------------------
*     Calculate the differential volume on half mesh:
*     ----------------------------------------------

      dvol  = 1.0e-05

      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi
               n = int(rho(i,j,k) / drho) + 1
               if(n .le. nnoderho)then
                  dvol(n)  =  dvol(n) + dx * dy * dphi * capr(i)
               end if
            end do
         end do
      end do

      dldbavg = 1.0

*     ------------------------
*     Do flux surface averages:
*     ------------------------

      call fluxavg(xnea, xnavg, rho, nxmx, nymx, nphimx, nrhomax,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)

      call fluxavg(xn1a, xn1avg, rho, nxmx, nymx, nphimx, nrhomax,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)

      call fluxavg(xn2a, xn2avg, rho, nxmx, nymx, nphimx, nrhomax,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)

      call fluxavg(xn3a, xn3avg, rho, nxmx, nymx, nphimx, nrhomax,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)

      call fluxavg(xkte, xkteavg, rho, nxmx, nymx, nphimx, nrhomax,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)

      call fluxavg(xkti, xktiavg, rho, nxmx, nymx, nphimx, nrhomax,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)

      call fluxavg(xkti2, xkti2avg, rho, nxmx, nymx, nphimx, nrhomax,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)

      call fluxavg(xkti3, xkti3avg, rho, nxmx, nymx, nphimx, nrhomax,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)


      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey
               signb = b0 / abs(b0)

               omgce(i, j, k)  = qe  * bmod(i, j, k) / xme  * signb
               omgci1(i, j, k) = qi1 * bmod(i, j, k) / xmi1 * signb
               omgci2(i, j, k) = qi2 * bmod(i, j, k) / xmi2 * signb
               omgci3(i, j, k) = qi3 * bmod(i, j, k) / xmi3 * signb

               omgpe2(i, j, k) = xnea(i ,j, k) * qe**2  / (eps0 * xme)
               omgp12(i, j, k) = xn1a(i, j, k) * qi1**2 / (eps0 * xmi1)
               omgp22(i, j, k) = xn2a(i, j, k) * qi2**2 / (eps0 * xmi2)
               omgp32(i, j, k) = xn3a(i, j, k) * qi3**2 / (eps0 * xmi3)
            end do
         end do
      end do


      if(myid .eq. 0)then
         write(15, 920)
         write(15, 8164)
         write(15, 920)


         do k = 1, nnodephi

            iant = 1
            xjymax = 0.0
            do i = 1, nnodex
               modxjy = sqrt(    conjg(xjy(i, jequat, k))
     .                               * xjy(i, jequat, k)   )
               if(modxjy .gt. xjymax)then
                  xjymax = modxjy
                  iant = i
               end if
            end do

            write(15, 1312) k, xjymax
         end do

         write(15, 920)
         write(15, 164)
         write(15, 920)
      end if


      k = kplot
      do i = 1, nnodex
         do j = 1, nnodey

            reomg1 = omgrf / omgci1(i,j,k)
            reomg2 = omgrf / omgci2(i,j,k)
            reomg3 = omgrf / omgci3(i,j,k)

            scold = 1.0 - omgpe2(i,j,k) / (omgrf**2 - omgce(i,j,k)**2)
     .                  - omgp12(i,j,k) / (omgrf**2 - omgci1(i,j,k)**2)
     .                  - omgp22(i,j,k) / (omgrf**2 - omgci2(i,j,k)**2)
     .                  - omgp32(i,j,k) / (omgrf**2 - omgci3(i,j,k)**2)

            acold(i,j) = xnphi0**2 - scold


            if(myid .eq. 0)then
               if(j .eq. jequat)then
                  write(15, 2163)i, capr(i), x(i), bmod(i, j, k),
     .               reomg1, reomg2, real(acold(i,j))
               end if
            end if

         end do
      end do



      if(myid .eq. 0)then
         write(15, 920)
         write(15, 9164)
         write(15, 920)
      end if

      k = kplot
      do i = 1, nnodex
         do j = 1, nnodey

            reomg1 = omgrf / omgci1(i,j,k)
            reomg2 = omgrf / omgci2(i,j,k)
            reomg3 = omgrf / omgci3(i,j,k)

            scold = 1.0 - omgpe2(i,j,k) / (omgrf**2 - omgce(i,j,k)**2)
     .                  - omgp12(i,j,k) / (omgrf**2 - omgci1(i,j,k)**2)
     .                  - omgp22(i,j,k) / (omgrf**2 - omgci2(i,j,k)**2)
     .                  - omgp32(i,j,k) / (omgrf**2 - omgci3(i,j,k)**2)

            acold(i,j) = xnphi0**2 - scold


            if(myid .eq. 0)then
               if (i .eq. icenter)then
                  write(15, 2163)j, capr(i), y(j), bmod(i, j, k),
     1               reomg1, reomg2, real(acold(i,j))
               end if
            end if

         end do
      end do

      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey

               reson1(i,j,k) = omgrf / omgci1(i,j,k)
               reson2(i,j,k) = omgrf / omgci2(i,j,k)

            end do
         end do
      end do


      do n = nkx1, nkx2
         xkxsav(n) = two_pi * n / xmax
         if(n .eq. 0)xkxsav(n) = epszet
      end do

      do m = nky1, nky2
         xkysav(m) = two_pi * m / ymax
         if(m .eq. 0)xkysav(m) = epszet
      end do

      if (myid .eq. 0)then
         write(15, 920)
         write(15, 6164)
         write(15, 920)
      end if

! -------------------------
!     MCAP LOOP
! -------------------------

      do 9000 imcap = 1, mcap_number  ! Start of mcap loop
      if(myid .eq. 0)then
         write(15, 7114)mcap(imcap)
      endif

      time0 = second1(dummy)

      do nphi = nphi1, nphi2
         xkzsav(nphi) = nfp * nphi + mcap(imcap)  ! POL: Slide 33
         if(nphi .eq. 0 .and. mcap(imcap) .eq. 0)xkzsav(nphi) = epszet
         if(myid .eq. 0)write(15, 1312)nphi, xkzsav(nphi)
      end do

      kperp_max_actual = sqrt(xkxsav(nkx2)**2 + xkysav(nky2)**2
     .     + (xkzsav(nphi2) / capr(1))**2)

      kperp_max = kperp_max_actual * 2.0

      xkx_cutoff = xkxsav(nkx2) * xkperp_cutoff
      xky_cutoff = xkx_cutoff
      xkz_cutoff = xkzsav(nphi2) * xkperp_cutoff

*     ---------------------------
*     Calculate rotation matrix U
*     ---------------------------

      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi

               sqx = sqrt(1.0 - bxn(i, j, k)**2)

               uxx(i, j, k) =   sqx
               uxy(i, j, k) = - bxn(i, j, k) * byn(i, j, k) / sqx
               uxz(i, j, k) = - bxn(i, j, k) * bzn(i, j, k) / sqx
               uyx(i, j, k) =   0.0
               uyy(i, j, k) =   bzn(i, j, k) / sqx
               uyz(i, j, k) = - byn(i, j, k) / sqx
               uzx(i, j, k) =   bxn(i, j, k)
               uzy(i, j, k) =   byn(i, j, k)
               uzz(i, j, k) =   bzn(i, j, k)

            end do
         end do
      end do


*     ------------------------------------------------
*     Calculate parallel gradient and derivatives of U
*     ------------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi

               call deriv_x(bmod, nxmx, nymx, nphimx, i, j, k,
     .               nnodex, nnodey, nnodephi, dx, dbdx, dum)
               call deriv_y(bmod, nxmx, nymx, nphimx, i, j, k,
     .               nnodex, nnodey, nnodephi, dy, dbdy, dum)
               call deriv_phi(bmod, nxmx, nymx, nphimx, i, j, k,
     .               nnodex, nnodey, nnodephi, dphi, dbdphi, dum)

               gradprlb(i,j,k) = 0.0

               if(iflag_gammab .ne. 0)
     .              gradprlb(i,j,k) = uzx(i,j,k) * dbdx
     .                         + uzy(i,j,k) * dbdy
     .                         + uzz(i,j,k) / capr(i) * dbdphi



               call deriv_x(uxx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dxuxx(i,j,k), dxxuxx(i,j,k))
               call deriv_x(uxy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dxuxy(i,j,k), dxxuxy(i,j,k))
               call deriv_x(uxz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dxuxz(i,j,k), dxxuxz(i,j,k))

               call deriv_x(uyx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dxuyx(i,j,k), dxxuyx(i,j,k))
               call deriv_x(uyy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dxuyy(i,j,k), dxxuyy(i,j,k))
               call deriv_x(uyz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dxuyz(i,j,k), dxxuyz(i,j,k))

               call deriv_x(uzx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dxuzx(i,j,k), dxxuzx(i,j,k))
               call deriv_x(uzy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dxuzy(i,j,k), dxxuzy(i,j,k))
               call deriv_x(uzz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dxuzz(i,j,k), dxxuzz(i,j,k))



               call deriv_y(uxx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dyuxx(i,j,k), dyyuxx(i,j,k))
               call deriv_y(uxy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dyuxy(i,j,k), dyyuxy(i,j,k))
               call deriv_y(uxz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dyuxz(i,j,k), dyyuxz(i,j,k))

               call deriv_y(uyx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dyuyx(i,j,k), dyyuyx(i,j,k))
               call deriv_y(uyy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dyuyy(i,j,k), dyyuyy(i,j,k))
               call deriv_y(uyz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dyuyz(i,j,k), dyyuyz(i,j,k))

               call deriv_y(uzx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dyuzx(i,j,k), dyyuzx(i,j,k))
               call deriv_y(uzy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dyuzy(i,j,k), dyyuzy(i,j,k))
               call deriv_y(uzz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dyuzz(i,j,k), dyyuzz(i,j,k))



               call deriv_xy(uxx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dy, dxyuxx(i,j,k))
               call deriv_xy(uxy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dy, dxyuxy(i,j,k))
               call deriv_xy(uxz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dy, dxyuxz(i,j,k))

               call deriv_xy(uyx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dy, dxyuyx(i,j,k))
               call deriv_xy(uyy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dy, dxyuyy(i,j,k))
               call deriv_xy(uyz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dy, dxyuyz(i,j,k))

               call deriv_xy(uzx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dy, dxyuzx(i,j,k))
               call deriv_xy(uzy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dy, dxyuzy(i,j,k))
               call deriv_xy(uzz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dy, dxyuzz(i,j,k))




               call deriv_phi(uxx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dphi,dphiuxx(i,j,k),d2phiuxx(i,j,k))
               call deriv_phi(uxy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dphi,dphiuxy(i,j,k),d2phiuxy(i,j,k))
               call deriv_phi(uxz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dphi,dphiuxz(i,j,k),d2phiuxz(i,j,k))

               call deriv_phi(uyx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dphi,dphiuyx(i,j,k),d2phiuyx(i,j,k))
               call deriv_phi(uyy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dphi,dphiuyy(i,j,k),d2phiuyy(i,j,k))
               call deriv_phi(uyz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dphi,dphiuyz(i,j,k),d2phiuyz(i,j,k))

               call deriv_phi(uzx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dphi,dphiuzx(i,j,k),d2phiuzx(i,j,k))
               call deriv_phi(uzy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dphi,dphiuzy(i,j,k),d2phiuzy(i,j,k))
               call deriv_phi(uzz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dphi,dphiuzz(i,j,k),d2phiuzz(i,j,k))



               call deriv_xphi(uxx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dphi, dxphiuxx(i,j,k))
               call deriv_xphi(uxy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dphi, dxphiuxy(i,j,k))
               call deriv_xphi(uxz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dphi, dxphiuxz(i,j,k))

               call deriv_xphi(uyx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dphi, dxphiuyx(i,j,k))
               call deriv_xphi(uyy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dphi, dxphiuyy(i,j,k))
               call deriv_xphi(uyz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dphi, dxphiuyz(i,j,k))

               call deriv_xphi(uzx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dphi, dxphiuzx(i,j,k))
               call deriv_xphi(uzy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dphi, dxphiuzy(i,j,k))
               call deriv_xphi(uzz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dx, dphi, dxphiuzz(i,j,k))



               call deriv_yphi(uxx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dphi, dyphiuxx(i,j,k))
               call deriv_yphi(uxy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dphi, dyphiuxy(i,j,k))
               call deriv_yphi(uxz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dphi, dyphiuxz(i,j,k))

               call deriv_yphi(uyx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dphi, dyphiuyx(i,j,k))
               call deriv_yphi(uyy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dphi, dyphiuyy(i,j,k))
               call deriv_yphi(uyz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dphi, dyphiuyz(i,j,k))

               call deriv_yphi(uzx, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dphi, dyphiuzx(i,j,k))
               call deriv_yphi(uzy, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dphi, dyphiuzy(i,j,k))
               call deriv_yphi(uzz, nxmx, nymx, nphimx, i, j, k, nnodex,
     .            nnodey, nnodephi, dy, dphi, dyphiuzz(i,j,k))


            end do
         end do
      end do


c     -------------------------------------
c     setup Z function table when nzfun = 3
c     -------------------------------------
      if (nzfun .eq. 3) then

         iflag = 0

         call ztable_setup(myid, iflag)

         call blacs_barrier(icontxt, 'All')

         if(iflag .eq. 1)call ztable_setup(myid, nproc)

      end if



c     ----------------------------------------
c     precompute xx(n,i), yy(m,j), zz(nphi, k)
c     ----------------------------------------

        allocate( inv_xx(nkx1:nkx2,1:nnodex) )
        allocate( inv_yy(nky1:nky2,1:nnodey) )
        allocate( inv_zz(nphi1:nphi2,1:nnodephi))

        allocate( inv_xx_t(1:nnodex,nkx1:nkx2) )
        allocate( inv_yy_t(1:nnodey,nky1:nky2) )
        allocate( inv_zz_t(1:nnodephi,nphi1:nphi2) )

      do i = 1, nnodex
         do n = nkx1, nkx2
            xx(n, i) = exp(zi * xkxsav(n) * xprime(i))
            inv_xx(n,i) = 1.0/xx(n,i)
            inv_xx_t(i,n) = inv_xx(n,i)
         end do
      end do

      do j = 1, nnodey
         do m = nky1, nky2
            yy(m, j) = exp(zi * xkysav(m) * yprime(j))
            inv_yy(m,j) = 1.0/yy(m,j)
            inv_yy_t(j,m) = inv_yy(m,j)
         end do
      end do

      do k = 1, nnodephi
         do nphi = nphi1, nphi2
            zz(nphi, k) = exp(zi * xkzsav(nphi) * phiprime(k))
            inv_zz(nphi,k) = 1.0/zz(nphi,k)
            inv_zz_t(k,nphi) = inv_zz(nphi,k)
         end do
      end do



*       ------------------------------------
*       determine reordering and permutation
*       ------------------------------------

*       -------------------
*       add layer of points
*       -------------------

        num_inside = 0
        do k=1,nnodephi
        do j=1,nnodey
        do i=1,nnodex
          if (mask(i,j,k).eq.1) then
                num_inside = num_inside + 1
          endif
        enddo
        enddo
        enddo

        allocate( new_to_org( 3*max(1,num_inside) ) )

        num_new = 0
        num_org = 0
        do k=1,nnodephi
        do j=1,nnodey
        do i=1,nnodex
                irnc = (k-1) * nnodex * nnodey * 3
     .             + (i-1) * nnodey * 3
     .             + (j-1) * 3
     .             + 1

           num_org = irnc

           if (mask(i,j,k).eq.1) then

                num_new = num_new + 1
                new_to_org( num_new ) = num_org
                num_new = num_new + 1
                new_to_org( num_new) = num_org + 1
                num_new = num_new + 1
                new_to_org( num_new) = num_org + 2
           endif
        enddo
        enddo
        enddo


*       ----------------------------------------------
*       note we wish to be solving a smaller
*       matrix and ignoring points outside the plasma
*       ----------------------------------------------

        org_nrow = 3*nnodex*nnodey*nnodephi
        org_ncol = org_nrow

        new_nrow = 3*num_inside
        new_ncol = new_nrow

        nrow = new_nrow
        ncol = nrow

        if ((myrow.eq.0).and.(mycol.eq.0)) then
          write(15,*) 'new_nrow,org_nrow ', new_nrow,org_nrow
        endif




      t1=second1(dummy)



*     ----------------------------------
*     define scalapack array descriptors
*     ----------------------------------
      mb = 3*33
      nb = mb


*       --------------------------------------------------------
*       setting mmb = mb is convenient but may take up too much
*       memory, especially for large problems
*
*       mmb and mb should be divisible by 3, mod(mmb,3).eq.0
*       The reduced matrix size is 39492 and there are 1936 cpus so on average,
*       there are about 21 rows to transform.
*       Since the time taken to invert Fourier transform is about 0.044min, just
*       the time  required to perform 21 transforms  is about 1min. However, if
*       the variable "mmb" is set too high, say mmb=mb=99, then 99 rows are
*       transformed at a time leaving some cpus  idle but some cpus need to
*       transform about 99 rows. The time for transforming 99 rows is about 4.4min.

*       A workaround is to set variable mmb to a smaller value (must be
*       divisible by 3) to say 21 (in this case) or even as small as mmb=3. This
*       should more evenly divide the work and most likely reduce the load time.

*       --------------------------------------------------------
        mmb = mb
        isok = (mmb.ge.3).and.(mod(mmb,3).eq.0)
        if (.not.isok) then
          mmb = 3
        endif

         mmb = 21

*       -------------------------------
*       setup datastructures for Bmat_all(0:(nprow-1),0:(npcol-1))
*       each can be considered a trivial scalapack array of size mb by ncol
*       will all data for Bmat_all(irow,icol) residing
*       on processor (irow,icol)
*       -------------------------------

        allocate(niabegin_all(0:(nprow-1),0:(npcol-1)) )
        allocate( isize_all(0:(nprow-1),0:(npcol-1)) )
        niabegin_all(:,:) = 0
        isize_all(:,:) = 0

        allocate( row(org_ncol) )
        allocate( rowk(org_ncol) )

        allocate(Btmp(mmb,org_ncol))
        allocate(brhs(mmb,1))

        Btmp(:,:) = 0.0
        brhs(:,:) = 0.0

        allocate( descBtmp_all(DLEN_,0:(nprow-1),0:(npcol-1)) )
        allocate( descbrhs_all(DLEN_,0:(nprow-1),0:(npcol-1)) )
        descBtmp_all(:,:,:) = 0
        descbrhs_all(:,:,:) = 0

           Locp = mmb
           Locq = ncol
           lld = max(1,Locp)

        do icol=0,npcol-1
        do irow=0,nprow-1
           call descinit( descBtmp,mmb, ncol, mb, ncol,
     &                  irow,icol,icontxt, lld, info)
           if (info.ne.0) then
            if ((myrow.eq.0).and.(mycol.eq.0)) then
              write(*,9710) irow,icol,info
 9710         format('descinit for descBtmp: irow,icol ',
     &                 2(1x,i7),' info ',i7 )
            endif
           endif

           call descinit( descbrhs,mmb,1, mmb,1,
     &                  irow,icol,icontxt, lld, info )
           descBtmp_all(:,irow,icol) = descBtmp(:)
           descbrhs_all(:,irow,icol) = descbrhs(:)
        enddo
        enddo



      ncol = nrow
      norder = nrow

      nrow_local = max(1,numroc( nrow,mb,myrow,0,nprow ))
      ncol_local = max(1,numroc( ncol,nb,mycol,0,npcol ))
      lld = max(1,nrow_local)

      p_amat_dim = nrow_local * ncol_local + 1
      p_amat_size = p_amat_dim * 16 / 1.0e+06




        if (myid.eq.0) then
         write (15,6839) p_amat_dim, p_amat_size
      end if
 6839 format('p_amat_dim = ', i10, ' words,   ',
     .       'p_amat_size = ', i5, ' MBytes')

*     ----------------------------------------
*     allocate storage for scalapack matrices
*     ----------------------------------------

      allocate( p_amat(p_amat_dim) )
      p_amat(:) = cmplx(0.0,0.0)

      p_brhs_dim = nrow
      p_ipiv_dim = nrow
      allocate( p_brhs(p_brhs_dim), p_ipiv(p_ipiv_dim) )
      p_brhs(:) = cmplx(0.0,0.0)
      p_ipiv(:) = 0


      call descinit( desc_amat, nrow,ncol, mb,nb, 0,0, icontxt,
     &         lld, info )
      if (info.ne.0) then
        write(6,*) '** descinit for amat returns info = ', info
        write(6,*) 'lld = ', lld
        write(15,*) '** descinit for amat returns info = ', info
        write(15,*) 'lld = ', lld
        stop '** error 1 ** '
      endif

      do irnc=1,nrow_local*ncol_local
          p_amat(irnc) = cmplx(0., 0.)
      enddo

      if (lld .gt. p_brhs_dim) then
        write(6,*) 'increase p_brhs_dim to ', lld+1
        write(15,*) 'increase p_brhs_dim to ', lld+1
        stop '** error 2 ** '
      endif


      call descinit( desc_brhs, nrow, 1, mb,nb, 0,0, icontxt,
     &         lld, info )
      if (info.ne.0) then
        write(6,*) '** descinit for brhs returns info ', info
        write(15,*) '** descinit for brhs returns info ', info
        stop '** error 3 ** '
      endif

! --------------------------------------
! precompute mapping from ia->(i,j,k)
! precompute mapping from ja->(n,m,nphi)
! --------------------------------------

       allocate( itable(3 * nnodex * nnodey * nnodephi) )
       allocate( jtable(3 * nnodex * nnodey * nnodephi) )
       allocate( ktable(3 * nnodex * nnodey * nnodephi) )

       allocate( mtable(3 * nnodex * nnodey * nnodephi) )
       allocate( ntable(3 * nnodex * nnodey * nnodephi) )
       allocate( nphitable(3 * nnodex * nnodey * nnodephi) )

       do i = 1, 3 * nnodex * nnodey * nnodephi
          itable(i) = undefined
          jtable(i) = undefined
          ktable(i) = undefined

          mtable(i) = undefined
          ntable(i) = undefined
          nphitable(i) = undefined
       enddo


       do i = 1, nnodex
          do j = 1, nnodey
             do k = 1, nnodephi
                irnc = (k-1) * nnodex * nnodey * 3
     .             + (i-1) * nnodey * 3
     .             + (j-1) * 3
     .             + 1
                itable(irnc) = i
                jtable(irnc) = j
                ktable(irnc) = k
             end do
          end do
       end do


       do nphi = nphi1, nphi2
          do n = nkx1, nkx2
             do m = nky1, nky2

                icnc = (nphi - nphi1) * 3 * nmodesx * nmodesy
     .             + (n - nkx1) * 3 * nmodesy
     .             + (m - nky1) * 3 + 1

                ntable(icnc) = n
                mtable(icnc) = m
                nphitable(icnc) = nphi
             enddo
          enddo
       end do


*     -----------------------------------------------------------------------
*     Load x, y, and z equations for spatial pint (i,j,k) and mode (n,m,nphi)
*     -----------------------------------------------------------------------


        nnb = nprow*npcol*mmb
        do niastart=1,desc_amat(M_),nnb
           niaend = min(desc_amat(M_), niastart+nnb-1)
           niasize = niaend-niastart+1
           if (niasize.le.0) exit


*       -------------------------------------------
*       determine what rows need to be computed for
*       (irow,icol) processor
*
*       important to preserve fortran column oriented
*       assignment of block rows to processor grid
*       -------------------------------------------
           nia = niastart
           do icol=0,(npcol-1)
           do irow=0,(nprow-1)
             niabegin_all(irow,icol) = nia
             isize_all(irow,icol) = max(0,min(mmb, niasize))

             nia = nia + isize_all(irow,icol)
             niasize = niasize - isize_all(irow,icol)
           enddo
           enddo



           do ni=1,isize_all(myrow,mycol),3

*           -------------------------------------------
*           construct the entire row of original matrix
*           -------------------------------------------
            nia = (ni-1) + niabegin_all(myrow,mycol)
            ia = new_to_org( nia )

            do ja=1,org_ncol,3



!           -----------------------------------------------
!           (ia,ja) are all local indices in a single block
!
!           note each ia and ja loop has increment by 3
!           to match the original code
!           -----------------------------------------------


!                 ------------------------------------
!                 obtain (i,j,k,m,n,nphi) from (ia,ja)
!                 ------------------------------------

                  i = itable(ia)
                  j = jtable(ia)
                  k = ktable(ia)
                  m = mtable(ja)
                  n = ntable(ja)
                  nphi = nphitable(ja)


                  isok = (1.le.i).and.(i.le.nnodex)  .and.
     &                   (1.le.j).and.(j.le.nnodey)  .and.
     &                   (1.le.k).and.(k.le.nnodephi).and.
     &                   (nkx1.le.n).and.(n.le.nkx2) .and.
     &                   (nky1.le.m).and.(m.le.nky2) .and.
     &                   (nphi1.le.nphi).and.(nphi.le.nphi2)

                  if (.not.isok) then
                      write(*,*) 'i,j,k ',i,j,k
                      write(*,*) 'm,n,nphi, ',m,n,nphi
                      write(*,*) 'ia,ja ', ia,ja
                      stop 'error in main '
                  endif

!                 ------------
!                 double check
!                 ------------

                  irnc = (k-1) * nnodex * nnodey * 3
     .                 + (i-1) * nnodey * 3
     .                 + (j-1) * 3
     .                 +  1


                  icnc = (nphi - nphi1) * 3 * nnodex * nnodey
     .                 +   (n - nkx1) * 3 * nnodey
     .                 +   (m - nky1) * 3
     .                 + 1


                  isok = (irnc.eq.ia).and.(icnc.eq.ja)
                  if (.not.isok) then
                      write(*,*) 'main: ia,ja ', ia,ja
                      write(*,*) 'i,j,k ',i,j,k
                      write(*,*) 'n,m,nphi ',n,m,nphi
                      write(*,*) 'main: irnc,icnc ',irnc,icnc
                      stop 'error in current '
                  endif



!                 -----------------------------
!                 copy variable to be compatible with previous code
!                 -----------------------------
                  irnc1 = irnc
                  icnc1 = icnc
                  rsrc1 = myrow
                  csrc1 = mycol


                     xb = -zi / omgrf / eps0 * xjx(i,j,k)
                     xc = -zi / omgrf / eps0 * xjy(i,j,k)
                     xd = -zi / omgrf / eps0 * xjz(i,j,k)

                     xkphi = xkzsav(nphi) / capr(i)

                     cexpkxkykz = xx(n,i) * yy(m,j) * zz(nphi,k)

                     sigxx = 0.0
                     sigxy = 0.0
                     sigxz = 0.0
                     sigyx = 0.0
                     sigyy = 0.0
                     sigyz = 0.0
                     sigzx = 0.0
                     sigzy = 0.0
                     sigzz = 0.0

*                    ----------------------
*                    interior plasma region:
*                    ----------------------
                     if(psi(i,j,k) .le. psilim) then


                        if (isigma .eq. 1)then

                           call sigmad_cql3d(i, j, k, n, m,
     .                        rho(i, j, k), rho_a,
     .                        gradprlb(i, j, k), bmod(i, j, k),
     .                        bmod_mid(i, j, k),
     .                        xme, qe, xnea(i, j, k), xnuomg,
     .                        xkte(i,j,k), omgce(i,j,k), omgpe2(i,j,k),
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkxsav(n), xkysav(m), xkzsav(nphi),
     .                        capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sigexx, sigexy, sigexz,
     .                        sigeyx, sigeyy, sigeyz,
     .                        sigezx, sigezy, sigezz,
     .                        delta0, ndiste, nupar, nuper, n_psi,
     .                        n_psi_dim, dfdupere, dfdupare,
     .                        UminPara, UmaxPara, UPERP, UPARA,
     .                        vce_mks, dfe_cql_uprp, dfe_cql_uprl,
     .                        nbessj, nkperp, zi, eps0, v0i, omgrf,
     .                        xk0, kperp_max,
     .                        i_sav, j_sav, k_sav, 1, damping,
     .                        xkx_cutoff, xky_cutoff, xkz_cutoff,
     .                        rt, nkx2, nky2)


                           call sigmad_cql3d(i, j, k, n, m,
     .                        rho(i, j, k), rho_a,
     .                        gradprlb(i, j, k), bmod(i, j, k),
     .                        bmod_mid(i, j, k),
     .                        xmi1, qi1, xn1a(i, j, k), xnuomg,
     .                        xkti(i,j,k), omgci1(i,j,k), omgp12(i,j,k),
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkxsav(n), xkysav(m), xkzsav(nphi),
     .                        capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sig1xx, sig1xy, sig1xz,
     .                        sig1yx, sig1yy, sig1yz,
     .                        sig1zx, sig1zy, sig1zz,
     .                        delta0, ndisti1, nupar, nuper, n_psi,
     .                        n_psi_dim, dfduper1, dfdupar1,
     .                        UminPara, UmaxPara, UPERP, UPARA,
     .                        vc1_mks, df1_cql_uprp, df1_cql_uprl,
     .                        nbessj, nkperp, zi, eps0, v0i, omgrf,
     .                        xk0, kperp_max,
     .                        i_sav, j_sav, k_sav, 1, damping,
     .                        xkx_cutoff, xky_cutoff, xkz_cutoff,
     .                        rt, nkx2, nky2)


                           if(eta2 .ne. 0.0)
     .                        call sigmad_cql3d(i, j, k, n, m,
     .                        rho(i, j, k), rho_a,
     .                        gradprlb(i, j, k), bmod(i, j, k),
     .                        bmod_mid(i, j, k),
     .                        xmi2, qi2, xn2a(i, j, k), xnuomg,
     .                        xkti2(i,j,k), omgci2(i,j,k),omgp22(i,j,k),
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkxsav(n), xkysav(m), xkzsav(nphi),
     .                        capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sig2xx, sig2xy, sig2xz,
     .                        sig2yx, sig2yy, sig2yz,
     .                        sig2zx, sig2zy, sig2zz,
     .                        delta0, ndisti2, nupar, nuper, n_psi,
     .                        n_psi_dim, dfduper2, dfdupar2,
     .                        UminPara, UmaxPara, UPERP, UPARA,
     .                        vc2_mks, df2_cql_uprp, df2_cql_uprl,
     .                        nbessj, nkperp, zi, eps0, v0i, omgrf,
     .                        xk0, kperp_max,
     .                        i_sav, j_sav, k_sav, 1, damping,
     .                        xkx_cutoff, xky_cutoff, xkz_cutoff,
     .                        rt, nkx2, nky2)


                           if(eta3 .ne. 0.0)
     .                        call sigmad_cql3d(i, j, k, n, m,
     .                        rho(i, j, k), rho_a,
     .                        gradprlb(i, j, k), bmod(i, j, k),
     .                        bmod_mid(i, j, k),
     .                        xmi3, qi3, xn3a(i, j, k), xnuomg,
     .                        xkti3(i,j,k), omgci3(i,j,k),omgp32(i,j,k),
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkxsav(n), xkysav(m), xkzsav(nphi),
     .                        capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sig3xx, sig3xy, sig3xz,
     .                        sig3yx, sig3yy, sig3yz,
     .                        sig3zx, sig3zy, sig3zz,
     .                        delta0, ndisti3, nupar, nuper, n_psi,
     .                        n_psi_dim, dfduper3, dfdupar3,
     .                        UminPara, UmaxPara, UPERP, UPARA,
     .                        vc3_mks, df3_cql_uprp, df3_cql_uprl,
     .                        nbessj, nkperp, zi, eps0, v0i, omgrf,
     .                        xk0, kperp_max,
     .                        i_sav, j_sav, k_sav, 1, damping,
     .                        xkx_cutoff, xky_cutoff, xkz_cutoff,
     .                        rt, nkx2, nky2)

                        end if ! (isigma.eq.1)


                        if (isigma .eq. 0)then

                           call sigmac_stix(i, j, n, m,
     .                        xme, qe, xnea(i, j, k), xnuomg,
     .                        xkte(i,j,k), omgce(i,j,k), omgpe2(i,j,k),
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkxsav(n), xkysav(m), xkzsav(nphi),
     .                        capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sigexx, sigexy, sigexz,
     .                        sigeyx, sigeyy, sigeyz,
     .                        sigezx, sigezy, sigezz)

                           call sigmac_stix(i, j, n, m,
     .                        xmi1, qi1, xn1a(i, j, k), xnuomg,
     .                        xkti(i,j,k), omgci1(i,j,k), omgp12(i,j,k),
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkxsav(n), xkysav(m), xkzsav(nphi),
     .                        capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sig1xx, sig1xy, sig1xz,
     .                        sig1yx, sig1yy, sig1yz,
     .                        sig1zx, sig1zy, sig1zz)

                           if(eta2 .ne. 0.0)
     .                        call sigmac_stix(i, j, n, m,
     .                        xmi2, qi2, xn2a(i, j, k), xnuomg,
     .                        xkti2(i,j,k), omgci2(i,j,k),omgp22(i,j,k),
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkxsav(n), xkysav(m), xkzsav(nphi),
     .                        capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sig2xx, sig2xy, sig2xz,
     .                        sig2yx, sig2yy, sig2yz,
     .                        sig2zx, sig2zy, sig2zz)

                           if(eta3 .ne. 0.0)
     .                        call sigmac_stix(i, j, n, m,
     .                        xmi3, qi3, xn3a(i, j, k), xnuomg,
     .                        xkti3(i,j,k), omgci3(i,j,k),omgp32(i,j,k),
     .                        -lmax, lmax, nzfun, ibessel,
     .                        xkxsav(n), xkysav(m), xkzsav(nphi),
     .                        capr(i),
     .                        bxn(i,j,k), byn(i,j,k), bzn(i,j,k),
     .                        uxx(i,j,k), uxy(i,j,k), uxz(i,j,k),
     .                        uyx(i,j,k), uyy(i,j,k), uyz(i,j,k),
     .                        uzx(i,j,k), uzy(i,j,k), uzz(i,j,k),
     .                        sig3xx, sig3xy, sig3xz,
     .                        sig3yx, sig3yy, sig3yz,
     .                        sig3zx, sig3zy, sig3zz)


                        end if ! (isigma.eq.0)


                        sigxx = sigexx + sig1xx + sig2xx + sig3xx
                        sigxy = sigexy + sig1xy + sig2xy + sig3xy
                        sigxz = sigexz + sig1xz + sig2xz + sig3xz

                        sigyx = sigeyx + sig1yx + sig2yx + sig3yx
                        sigyy = sigeyy + sig1yy + sig2yy + sig3yy
                        sigyz = sigeyz + sig1yz + sig2yz + sig3yz

                        sigzx = sigezx + sig1zx + sig2zx + sig3zx
                        sigzy = sigezy + sig1zy + sig2zy + sig3zy
                        sigzz = sigezz + sig1zz + sig2zz + sig3zz



                        xkxx = 1.0 + zi / (eps0 * omgrf) * sigxx
                        xkxy =       zi / (eps0 * omgrf) * sigxy
                        xkxz =       zi / (eps0 * omgrf) * sigxz

                        xkyx =       zi / (eps0 * omgrf) * sigyx
                        xkyy = 1.0 + zi / (eps0 * omgrf) * sigyy
                        xkyz =       zi / (eps0 * omgrf) * sigyz

                        xkzx =       zi / (eps0 * omgrf) * sigzx
                        xkzy =       zi / (eps0 * omgrf) * sigzy
                        xkzz = 1.0 + zi / (eps0 * omgrf) * sigzz


                        rnx   = xkxsav(n) / xk0
                        rny   = xkysav(m) / xk0
                        rnphi = xkzsav(nphi) / capr(i) / xk0


                        dxx = (xkxx - rny**2 - rnphi**2) * uxx(i,j,k)
     .                      + xkyx * uyx(i,j,k)
     .                      + xkzx * uzx(i,j,k)
     .                      + rnx * (rnphi * uxz(i,j,k)
     .                               + rny * uxy(i,j,k))
     .                      - zi * rnphi / xk0 *
     .                      (uxz(i,j,k) / capr(i) + dxuxz(i,j,k)
     .                         - 2. / capr(i) * dphiuxx(i,j,k))
     .                      - zi * rny / xk0 * (dxuxy(i,j,k)
     .                                           - 2. * dyuxx(i,j,k))
     .                      - zi * rnx / xk0 * (dphiuxz(i,j,k) / capr(i)
     .                                                + dyuxy(i,j,k))
     .                      + 1./ xk0**2 * (dyyuxx(i,j,k) -dxyuxy(i,j,k)
     .                                    - dxphiuxz(i,j,k) / capr(i))
     .                      - 1./ (xk0 * capr(i))**2 * (dphiuxz(i,j,k)
     .                                          - d2phiuxx(i,j,k))

                        dxy =  xkxy * uxx(i,j,k)
     .                      + (xkyy - rny**2 - rnphi**2) * uyx(i,j,k)
     .                      +  xkzy * uzx(i,j,k)
     .                      + rnx * (rnphi * uyz(i,j,k)
     .                               + rny * uyy(i,j,k))
     .                      - zi * rnphi / xk0 *
     .                        (uyz(i,j,k) / capr(i) + dxuyz(i,j,k)
     .                         - 2. / capr(i) * dphiuyx(i,j,k))
     .                      - zi * rny / xk0 * (dxuyy(i,j,k)
     .                                           - 2. * dyuyx(i,j,k))
     .                      - zi * rnx / xk0 * (dphiuyz(i,j,k) / capr(i)
     .                                                + dyuyy(i,j,k))
     .                      + 1./ xk0**2 * (dyyuyx(i,j,k) -dxyuyy(i,j,k)
     .                                    - dxphiuyz(i,j,k) / capr(i))
     .                      - 1./ (xk0 * capr(i))**2 * (dphiuyz(i,j,k)
     .                                          - d2phiuyx(i,j,k))

                        dxz =  xkxz * uxx(i,j,k)
     .                      + xkyz * uyx(i,j,k)
     .                      + (xkzz - rny**2 - rnphi**2) * uzx(i,j,k)
     .                      + rnx * (rnphi * uzz(i,j,k)
     .                               + rny * uzy(i,j,k))
     .                      - zi * rnphi / xk0 *
     .                         (uzz(i,j,k) / capr(i) + dxuzz(i,j,k)
     .                          - 2. / capr(i) * dphiuzx(i,j,k))
     .                      - zi * rny / xk0 * (dxuzy(i,j,k)
     .                                           - 2. * dyuzx(i,j,k))
     .                      - zi * rnx / xk0 * (dphiuzz(i,j,k) / capr(i)
     .                                                + dyuzy(i,j,k))
     .                      + 1./ xk0**2 * (dyyuzx(i,j,k) -dxyuzy(i,j,k)
     .                                   - dxphiuzz(i,j,k) / capr(i))
     .                      - 1./ (xk0 * capr(i))**2 * (dphiuzz(i,j,k)
     .                                          - d2phiuzx(i,j,k))


                        dyx = (xkxx - rnx**2 - rnphi**2) * uxy(i,j,k)
     .                      +  xkyx * uyy(i,j,k)
     .                      +  xkzx * uzy(i,j,k)
     .                      + rny * (rnphi * uxz(i,j,k)
     .                               + rnx * uxx(i,j,k))
     .                      - zi * rny / xk0 *
     .                         (dxuxx(i,j,k) + uxx(i,j,k) / capr(i)
     .                                       + dphiuxz(i,j,k) / capr(i))
     .                      - zi * rnphi / xk0 * (dyuxz(i,j,k)
     .                                  - 2. / capr(i) * dphiuxy(i,j,k))
     .                      - zi * rnx / xk0 *
     .                         (dyuxx(i,j,k) - uxy(i,j,k) / capr(i)
     .                                              - 2. * dxuxy(i,j,k))
     .                      - 1. / (xk0**2 * capr(i)) * (dyuxx(i,j,k)
     .                                 - dxuxy(i,j,k) + dyphiuxz(i,j,k))
     .                      + 1. / (xk0 * capr(i))**2 * d2phiuxy(i,j,k)
     .                      - 1. / xk0**2 * (dxyuxx(i,j,k)
     .                                     - dxxuxy(i,j,k))

                        dyy =  xkxy * uxy(i,j,k)
     .                      + (xkyy - rnx**2 - rnphi**2) * uyy(i,j,k)
     .                      +  xkzy * uzy(i,j,k)
     .                      + rny * (rnphi * uyz(i,j,k)
     .                               + rnx * uyx(i,j,k))
     .                      - zi * rny / xk0 *
     .                         (dxuyx(i,j,k) + uyx(i,j,k) / capr(i)
     .                                       + dphiuyz(i,j,k) / capr(i))
     .                      - zi * rnphi / xk0 * (dyuyz(i,j,k)
     .                                  - 2. / capr(i) * dphiuyy(i,j,k))
     .                      - zi * rnx / xk0 *
     .                         (dyuyx(i,j,k) - uyy(i,j,k) / capr(i)
     .                                              - 2. * dxuyy(i,j,k))
     .                      - 1. / (xk0**2 * capr(i)) * (dyuyx(i,j,k)
     .                                 - dxuyy(i,j,k) + dyphiuyz(i,j,k))
     .                      + 1. / (xk0 * capr(i))**2 * d2phiuyy(i,j,k)
     .                      - 1. / xk0**2 * (dxyuyx(i,j,k)
     .                                     - dxxuyy(i,j,k))

                        dyz =  xkxz * uxy(i,j,k)
     .                      +  xkyz * uyy(i,j,k)
     .                      + (xkzz - rnx**2 - rnphi**2) * uzy(i,j,k)
     .                      + rny * (rnphi * uzz(i,j,k)
     .                               + rnx * uzx(i,j,k))
     .                      - zi * rny / xk0 *
     .                         (dxuzx(i,j,k) + uzx(i,j,k) / capr(i)
     .                                       + dphiuzz(i,j,k) / capr(i))
     .                      - zi * rnphi / xk0 * (dyuzz(i,j,k)
     .                                  - 2. / capr(i) * dphiuzy(i,j,k))
     .                      - zi * rnx / xk0 *
     .                         (dyuzx(i,j,k) - uzy(i,j,k) / capr(i)
     .                                              - 2. * dxuzy(i,j,k))
     .                      - 1. / (xk0**2 * capr(i)) * (dyuzx(i,j,k)
     .                                 - dxuzy(i,j,k) + dyphiuzz(i,j,k))
     .                      + 1. / (xk0 * capr(i))**2 * d2phiuzy(i,j,k)
     .                      - 1. / xk0**2 * (dxyuzx(i,j,k)
     .                                     - dxxuzy(i,j,k))

                        dzx = (xkxx - rnx**2 - rny**2) * uxz(i,j,k)
     .                      +  xkyx * uyz(i,j,k)
     .                      +  xkzx * uzz(i,j,k)
     .                      + rnphi * (rny * uxy(i,j,k)
     .                               + rnx * uxx(i,j,k))
     .                      - zi * rny / xk0 * (dphiuxy(i,j,k) / capr(i)
     .                                            - 2. * dyuxz(i,j,k))
     .                      - zi * rnphi / xk0 *
     .                         (dyuxy(i,j,k) + dxuxx(i,j,k)
     .                                         - uxx(i,j,k) / capr(i))
     .                      - zi * rnx / xk0 * (dphiuxx(i,j,k) / capr(i)
     .                      - uxz(i,j,k) / capr(i) - 2.* dxuxz(i,j,k))
     .                      - 1./ (xk0 * capr(i))**2 * (uxz(i,j,k)
     .                                               - dphiuxx(i,j,k))
     .                      + 1./ xk0**2  * (dxxuxz(i,j,k)
     .                                      + dyyuxz(i,j,k))
     .                      - 1./ (xk0**2 * capr(i))  * (dxphiuxx(i,j,k)
     .                               + dyphiuxy(i,j,k) - dxuxz(i,j,k))

                        dzy =  xkxy * uxz(i,j,k)
     .                      + (xkyy - rnx**2 - rny**2) * uyz(i,j,k)
     .                      +  xkzy * uzz(i,j,k)
     .                      + rnphi * (rny * uyy(i,j,k)
     .                               + rnx * uyx(i,j,k))
     .                      - zi * rny / xk0 * (dphiuyy(i,j,k) / capr(i)
     .                                           - 2. * dyuyz(i,j,k))
     .                      - zi * rnphi / xk0 *
     .                         (dyuyy(i,j,k) + dxuyx(i,j,k)
     .                                         - uyx(i,j,k) / capr(i))
     .                      - zi * rnx / xk0 * (dphiuyx(i,j,k) / capr(i)
     .                      - uyz(i,j,k) / capr(i) - 2.* dxuyz(i,j,k))
     .                      - 1. / (xk0 * capr(i))**2 * (uyz(i,j,k)
     .                                               - dphiuyx(i,j,k))
     .                      + 1./ xk0**2  * (dxxuyz(i,j,k)
     .                                     + dyyuyz(i,j,k))
     .                      - 1./ (xk0**2 * capr(i))  * (dxphiuyx(i,j,k)
     .                               + dyphiuyy(i,j,k) - dxuyz(i,j,k))

                        dzz =  xkxz * uxz(i,j,k)
     .                      +  xkyz * uyz(i,j,k)
     .                      + (xkzz - rnx**2 - rny**2) * uzz(i,j,k)
     .                      + rnphi * (rny * uzy(i,j,k)
     .                               + rnx * uzx(i,j,k))
     .                      - zi * rny / xk0 * (dphiuzy(i,j,k) / capr(i)
     .                                            -2. * dyuzz(i,j,k))
     .                      - zi * rnphi / xk0 *
     .                      (dyuzy(i,j,k) + dxuzx(i,j,k)
     .                                         - uzx(i,j,k) / capr(i))
     .                      - zi * rnx / xk0 * (dphiuzx(i,j,k) / capr(i)
     .                      - uzz(i,j,k) / capr(i) - 2. * dxuzz(i,j,k))
     .                      - 1./ (xk0 * capr(i))**2 * (uzz(i,j,k)
     .                                               - dphiuzx(i,j,k))
     .                      + 1./ xk0**2  * (dxxuzz(i,j,k)
     .                                     + dyyuzz(i,j,k))
     .                      - 1./ (xk0**2 * capr(i))  * (dxphiuzx(i,j,k)
     .                               + dyphiuzy(i,j,k) - dxuzz(i,j,k))


                        fdk = dxx * cexpkxkykz
                        fek = dxy * cexpkxkykz
                        ffk = dxz * cexpkxkykz

                        fgk = dyx * cexpkxkykz
                        fak = dyy * cexpkxkykz
                        fpk = dyz * cexpkxkykz

                        frk = dzx * cexpkxkykz
                        fqk = dzy * cexpkxkykz
                        fsk = dzz * cexpkxkykz

                     end if ! (psi(i,j,k).le.psilim)

*                    -----------------------------
*                    metal boundary in edge region:
*                    -----------------------------
                     if(psi(i,j,k) .gt. psilim) then
                        fdk = cexpkxkykz
                        fek = 0.0
                        ffk = 0.0

                        fgk = 0.0
                        fak = cexpkxkykz
                        fpk = 0.0

                        frk = 0.0
                        fqk = 0.0
                        fsk = cexpkxkykz

                        xb = 0.0
                        xc = 0.0
                        xd = 0.0
                     end if


*       -------------------------------
*       copy the value into local array
*       -------------------------------
                        Btmp(ni,  ja) = fdk
                        Btmp(ni+1,ja) = fgk
                        Btmp(ni+2,ja) = frk

                        Btmp(ni,  ja+1) = fek
                        Btmp(ni+1,ja+1) = fak
                        Btmp(ni+2,ja+1) = fqk

                        Btmp(ni,  ja+2) = ffk
                        Btmp(ni+1,ja+2) = fpk
                        Btmp(ni+2,ja+2) = fsk

                        if (ja.eq.1) then
                                brhs(ni,1)   = xb
                                brhs(ni+1,1) = xc
                                brhs(ni+2,1) = xd
                        endif



               end do ! do ja
            end do ! do ni

*       -------------------------------------
*       perform transformation to real space
*       -------------------------------------

      do ni=1,isize_all(myrow,mycol)
        row(1:org_ncol) = Btmp(ni,1:org_ncol)

        call convert3d_row_kron( row, rowk,
     &    nnodex,nnodey,nnodephi,
     &    nkx1,nkx2, nky1,nky2, nphi1,nphi2,
     &    inv_xx_t, inv_yy_t,inv_zz_t,
     &    xmax, ymax, phimax, dx,dy, dphi)


*       -----------------------------------------
*       select only subset of variables in plasma
*       -----------------------------------------
        do nja=1,ncol
          ja = new_to_org(nja)
          Btmp(ni,nja) = rowk(ja)
        enddo

      enddo ! enddo ni

*       -------------------------
*       ready to copy into p_amat
*       -------------------------
        if (use_Btmp_all) then

        do icol=0,(npcol-1)
        do irow=0,(nprow-1)
          mm = isize_all(irow,icol)
          if (mm.le.0) cycle

          nn = ncol
           nia = niabegin_all(irow,icol)
           nja = 1
           descBtmp(:) = descBtmp_all(:,irow,icol)

           call pzgecopy( mm,nn,
     &          Btmp, 1,1, descBtmp,
     &          p_amat, nia,nja,desc_amat )

           descbrhs(:) = descbrhs_all(:,irow,icol)

           call pzgecopy( mm,1,
     &          brhs, 1,1, descbrhs,
     &          p_brhs, nia,1, desc_brhs )
        enddo
        enddo

        else
*       --------------------------------------------------------
*       perform communication by all nprow processors
*       (entire column on processor grid) at a time
*       --------------------------------------------------------
        mm = nprow*mmb
        nn = npcol*ncol
        Locp = mmb
        call descinit( descBtmp, mm,nn, mmb,ncol,
     &                  0,0,desc_amat(CTXT_), Locp, info )

        if (info.ne.0) then
          stop 'error with descinit(descBtmp)'
        endif

        mm = nprow*mmb
        nn = npcol
        Locp = mmb
        call descinit( descbrhs, mm,nn, mmb,1,
     &                  0,0,desc_brhs(CTXT_), Locp, info )
        if (info.ne.0) then
          stop 'error with descinit(descbrhs) '
        endif

        do icol=0,(npcol-1)
          mm = sum( isize_all(0:(nprow-1),icol) )
          if (mm.le.0) cycle

          nn = ncol
          ib = 1
          jb = 1 + (icol*ncol)
          nia = niabegin_all(0,icol)
          nja = 1

          call pzgecopy( mm,nn,
     &          Btmp, ib,jb, descBtmp,
     &          p_amat, nia,nja,desc_amat )

          ib = 1
          jb = 1+icol
          call pzgecopy( mm,1,
     &          brhs,ib,jb,descbrhs,
     &          p_brhs,nia,1,desc_brhs )
        enddo

        endif

        enddo ! do niastart



      time = second1(dummy) - t1
      tmin = time / 60.

      if (myid .eq. 0)then
           write(15,835) tmin
      end if
  835 format('time to load matrix =',f9.3,4h min)

      call blacs_barrier(icontxt, 'All')



*     --------------------------
*     scale each row by its norm
*     --------------------------
      global_col = 1
      incX = desc_amat(M_)

      do global_row = 1, nrow
        call pdznrm2( ncol, norm2, p_amat, global_row, global_col,
     .       desc_amat, incX )

        alpha = 1.0
        if (norm2 .ne. 0.0) alpha = 1.0 / norm2


        call pzscal ( ncol, alpha, p_amat, global_row, global_col,
     .       desc_amat, incX )

        call pzscal ( 1, alpha, p_brhs, global_row, 1,
     .       desc_brhs, desc_brhs(M_))


      enddo


      t1 = second1(dummy)

*       ---------------------
*       use scalapack solver
*       ---------------------
        call pzgetrf( nrow,ncol, p_amat, 1,1, desc_amat, p_ipiv, info )
        if (info.ne.0) then
           write(6,*) 'pzgetrf returns info = ', info
           write(15,*) 'pzgetrf returns info = ', info
           stop '** error 4 ** '
        endif

      time1 = second1(dummy) - t1
        tmin1 = time1 / 60.

      if (myid .eq. 0)then
           write(15,833) tmin1
      end if
  833 format('time taken by routine pzgetrf =',f9.3,4h min)

      call blacs_barrier(icontxt, 'All')



      t1 = second1(dummy)

        call pzgetrs( 'notrans', nrow, 1, p_amat, 1,1,desc_amat,
     &                p_ipiv, p_brhs, 1,1,desc_brhs, info )
        if (info.ne.0) then
           write(6,*) 'pzgetrs returns info = ', info
           write(15,*) 'pzgetrs returns info = ', info
           stop '** error 5 ** '
        endif

      time2 = second1(dummy) - t1
        tmin2 = time2 / 60.


      if (myid .eq. 0)then
           write(15,834) tmin2
      end if
  834 format('time taken by routine pzgetrs =',f9.3,4h min)

!       ---------------------------
!       release the bulk of storage
!       ---------------------------
        deallocate( p_amat )



      time = time1 + time2

c--   Operations for complex matrix factor and solve:
      ops = 8. / 3. * (real(nrow))**3 + 7. * (real(nrow))**2

      gflops  = ops / time / 1.0e+09
      gflopsp = gflops / nproc
        tmin = time / 60.


      if(myid .eq. 0)then
         write(15,839) nrow, nproc
           write(15,837) gflops
           write(15,838) gflopsp
      end if



  839 format('nrow =', i10, 5x, 'nproc =', i10)
  837 format('operations per sec by pzgetrf & pzgetrs ='
     1   ,f9.3,11h Gflops/sec)
  838 format('operations per sec per processor ='
     1   ,f9.3,21h Gflops/sec/processor)



      call blacs_barrier(icontxt, 'All')

*     ------------------
*     broadcast solution
*     ------------------
        allocate( brhs2(org_nrow) )
        allocate( brhs_tmp(new_nrow) )
        brhs2(:) = 0.0
        brhs_tmp(:) = 0.0

      icnc = 1
      do irnc=1,nrow
         call infog2l( irnc, icnc, desc_brhs, nprow,npcol,
     &                    myrow,mycol, lrindx,lcindx, rsrc,csrc )
         ismine = (rsrc .eq. myrow) .and. (csrc .eq. mycol)
         if (ismine) then
             ipos = lrindx + (lcindx-1)*desc_amat(LLD_)
             brhs_tmp(irnc) =  p_brhs(ipos)
         endif
      enddo
      call zgsum2d(icontxt, 'All', ' ', nrow, 1, brhs_tmp, nrow, -1, -1)


*     --------------------------------
*     expand solution to original size
*     --------------------------------
      do nia=1,nrow
        ia = new_to_org(nia)
        brhs2(ia) = brhs_tmp(nia)
      enddo



*     ---------------------------------------
*     adjust problem size
*     ---------------------------------------
      nrow = org_nrow
      ncol = nrow

*       ------------------------
*       convert to Fourier space
*       ------------------------
        do k=1,nnodephi
        do j=1,nnodey
        do i=1,nnodex
                  irnc = (k-1) * nnodex * nnodey * 3
     .                 + (i-1) * nnodey * 3
     .                 + (j-1) * 3
     .                 +  1
          ealpha(i,j,k) = brhs2(irnc)
          ebeta(i,j,k)  = brhs2(irnc+1)
          eb(i,j,k)     = brhs2(irnc+2)

        enddo
        enddo
        enddo


      call blacs_barrier(icontxt, 'All')

!     --------------------------------
!     write out solution in real space
!     --------------------------------

      if (myid .eq.0) then
         write(53,8310)(((ealpha(i,j,k),i=1,nnodex),
     &                 j=1,nnodey),k=1,nnodephi)
         write(53,8310)(((ebeta(i,j,k),i=1,nnodex),
     &                 j=1,nnodey),k=1,nnodephi)
         write(53,8310)(((eb(i,j,k),i=1,nnodex),
     &                 j=1,nnodey),k=1,nnodephi)

      end if

 8310 format(1p6e14.6)


      t1 = second1(dummy)

*     ------------------------------------------
*     Fourier transform the Stix electric fields:
*     ------------------------------------------

      call sft3d_kron(ealpha, ealphak,
     .   nxmx, nymx, nphimx,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz,
     .   xmax, ymax, phimax, dx,dy, dphi)

      call sft3d_kron(ebeta, ebetak,
     .   nxmx, nymx, nphimx,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz,
     .   xmax, ymax, phimax, dx,dy, dphi)

      call sft3d_kron(eb, ebk,
     .   nxmx, nymx, nphimx,
     .   nkdim1, nkdim2, mkdim1, mkdim2,  lkdim1, lkdim2,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz,
     .   xmax, ymax, phimax, dx,dy, dphi)

      if (anti_alias .eq. .true.) then
*        --------------------------------------
*        Anti-aliasing filter (two-thirds rule):
*        --------------------------------------
         n_upper = 2. / 3. * nkx2
         n_lower = 2. / 3. * nkx1
         m_upper = 2. / 3. * nky2
         m_lower = 2. / 3. * nky1
         l_upper = 2. / 3. * nphi2
         l_lower = 2. / 3. * nphi1

         if (myid .eq. 0) then
         write(6, *) "n_lower = ", n_lower, "   n_upper = ", n_upper
         write(6, *) "m_lower = ", m_lower, "   m_upper = ", m_upper
         write(6, *) "l_lower = ", l_lower, "   l_upper = ", l_upper
         end if

         do nphi = nphi1, nphi2
            do n = nkx1, nkx2
               do m = nky1, nky2

                  if (n .gt. n_upper .or. n .lt. n_lower .or.
     .                m .gt. m_upper .or. m .lt. m_lower .or.
     .                nphi .gt. l_upper .or. nphi .lt. l_lower)then

                     ealphak(n, m, nphi) = 0.0
                     ebetak(n, m, nphi) = 0.0
                     ebk(n, m, nphi) = 0.0
                  end if

               end do
            end do
         end do

         ealpha = 0.0
         ebeta  = 0.0
         eb     = 0.0

*        ---------------------------------------
*        Invert Fourier transform to real space:
*        ---------------------------------------

         call fftinv3d_kron(x, y, phi, xkxsav, xkysav, xkzsav,
     .        ealpha, ealphak, nxmx, nymx, nphimx,
     .        nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .        nnodex, nnodey, nnodephi,
     .        nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz)

         call fftinv3d_kron(x, y, phi, xkxsav, xkysav, xkzsav,
     .        ebeta, ebetak, nxmx, nymx, nphimx,
     .        nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .        nnodex, nnodey, nnodephi,
     .        nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz)

         call fftinv3d_kron(x, y, phi, xkxsav, xkysav, xkzsav,
     .        eb, ebk, nxmx, nymx, nphimx,
     .        nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .        nnodex, nnodey, nnodephi,
     .        nkx1, nkx2, nky1, nky2, nphi1, nphi2, xx, yy, zz)

      end if

      do nphi = nphi1, nphi2
         do n = nkx1, nkx2
            do m = nky1, nky2

               if(nphi .eq. nphiplot)then
                  ealphakmod(n, m) = sqrt(conjg(ealphak(n, m, nphi))
     .                                    * ealphak(n, m, nphi))
                  ebetakmod(n, m) = sqrt(conjg(ebetak(n, m, nphi))
     .                                    * ebetak(n, m, nphi))
                  ebkmod(n, m) = sqrt(conjg(ebk(n, m, nphi))
     .                                    * ebk(n, m, nphi))
               end if

               if(n .eq. nplot .and. m .eq. mplot)then
                  ealphakmod1(nphi) = sqrt(conjg(ealphak(n, m, nphi))
     .                                         * ealphak(n, m, nphi))
                  ebetakmod1(nphi) = sqrt(conjg(ebetak(n, m, nphi))
     .                                        * ebetak(n, m, nphi))
                  ebkmod1(nphi) = sqrt(conjg(ebk(n, m, nphi))
     .                                     * ebk(n, m, nphi))
               end if



            end do
         end do
      end do

!       -----------------------------------
!       write out solution in Fourier modes
!       -----------------------------------

      if (myid .eq.0) then
         write(54,309) nkx1,nkx2,nky1,nky2,nphi1,nphi2
         write(54,8310) (((ealphak(n,m,nphi),n=nkx1,nkx2),
     &             m=nky1,nky2),nphi=nphi1,nphi2)
         write(54,8310) (((ebetak(n,m,nphi),n=nkx1,nkx2),
     &             m=nky1,nky2),nphi=nphi1,nphi2)
         write(54,8310) (((ebk(n,m,nphi),n=nkx1,nkx2),
     &             m=nky1,nky2),nphi=nphi1,nphi2)

         write(54, 8310)(xkxsav(n), n = nkx1, nkx2)
         write(54, 8310)(xkysav(m), m = nky1, nky2)
         write(54, 8310)(xkzsav(nphi), nphi=nphi1,nphi2)
      endif

      time = second1(dummy) - t1
      tmin = time / 60.

      if(myid .eq. 0)then
         write(15, 1837) tmin
      end if

 1837 format('time taken to Fourier transform the fields =',f9.3,4h min)

c     ------------------------
c     Calculate E in Lab frame
c     ------------------------
      do k = 1, nnodephi
         do i = 1, nnodex
            do j= 1, nnodey

            ex(i,j,k)   = uxx(i,j,k) * ealpha(i,j,k)
     .                  + uyx(i,j,k) * ebeta(i,j,k)
     .                  + uzx(i,j,k) * eb(i,j,k)

            ey(i,j,k)   = uxy(i,j,k) * ealpha(i,j,k)
     .                  + uyy(i,j,k) * ebeta(i,j,k)
     .                  + uzy(i,j,k) * eb(i,j,k)

            ez(i,j,k)   = uxz(i,j,k) * ealpha(i,j,k)
     .                  + uyz(i,j,k) * ebeta(i,j,k)
     .                  + uzz(i,j,k) * eb(i,j,k)
            end do
         end do
      end do

      if (myid .eq. 0)then
         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 940)
         write(15, 920)


         i = idiag
         k = kdiag
         do j = 1, nnodey
            write(15, 900)i, j, ex(i,j,k), ey(i,j,k), ez(i,j,k),
     .                          eb(i,j,k)
         end do

         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 940)
         write(15, 920)

         j = jdiag
         k = kdiag
         do i = 1, nnodex
            write(15, 900)i, j, ex(i,j,k), ey(i,j,k), ez(i,j,k),
     .                          eb(i,j,k)
         end do

         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 950)
         write(15, 920)

         j = jdiag
         i = idiag
         do k = 1, nnodephi
            write(15, 900)k, j, ex(i,j,k), ey(i,j,k), ez(i,j,k),
     .                          eb(i,j,k)
         end do

      end if


  900 format(i5, i5, 1p9e12.3)
  901 format(i5, 1p9e12.3)
  940 format(1h , 4h   i, 5h    j,
     1                    4x,  6h re ex, 6x, 6h im ex,
     1                    6x,  6h re ey, 6x, 6h im ey,
     1                    6x,  6h re ez, 6x, 6h im ez,
     1                    6x,  6h re eb, 6x, 6h im eb)

  950 format(1h , 4h   k, 5h    j,
     1                    4x,  6h re ex, 6x, 6h im ex,
     1                    6x,  6h re ey, 6x, 6h im ey,
     1                    6x,  6h re ez, 6x, 6h im ez,
     1                    6x,  6h re eb, 6x, 6h im eb)




  910 format (1h1)
  920 format (1h0)
  930 format(3x, 'rf electric field')





 1940 format(1h , 4h   i, 5h    j,
     1                    4x,  7hre erho, 5x, 7him erho,
     1                    5x,  7hre eeta, 5x, 7him eeta,
     1                    6x,  7hre eb  , 5x, 7him eb  )
 1930 format(3x, 'rf electric field in Stix frame')


c
c-----calculate individual species currents in real space on solution mesh:
c

      t1 = second1(dummy)

      call ntilda_(ntilda_e, nxmx, nymx, nphimx,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   xme, qe, xnea, xkte, omgce, omgpe2, lmax,
     .   xkxsav, xkysav, xkzsav, nzfun, ibessel,
     .   ealphak, ebetak, ebk, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   icontxt,
     .   xx, yy, zz, xnuomg, psi, psilim,
     .   myid, nproc, gradprlb, bmod)

      call current_elect(xjpxe, xjpye, xjpze, nxmx, nymx, nphimx,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   xme, qe, xnea, xkte, omgce, omgpe2, lmax,
     .   xkxsav, xkysav, xkzsav, nzfun, ibessel,
     .   ealphak, ebetak, ebk, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   icontxt,
     .   xx, yy, zz, isigma, xnuomg, psi, psilim,
     .   myid, nproc, delta0, gradprlb, bmod,
     .   xk0, damping, xkx_cutoff, xky_cutoff, xkz_cutoff)


      call current(xjpx1, xjpy1, xjpz1, nxmx, nymx, nphimx,
     .   nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .   nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .   xmi1, qi1, xn1a, xkti, omgci1, omgp12, lmax,
     .   xkxsav, xkysav, xkzsav, nzfun, ibessel,
     .   ealphak, ebetak, ebk, capr,
     .   bxn, byn, bzn,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   icontxt,
     .   xx, yy, zz, isigma, xnuomg, psi, psilim,
     .   myid, nproc, delta0, gradprlb, bmod,
     .   xk0)

      if(eta2 .ne. 0.0)
     .   call current(xjpx2, xjpy2, xjpz2, nxmx, nymx, nphimx,
     .      nnodex, nnodey, nnodephi,
     .      nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .      nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .      xmi2, qi2, xn2a, xkti2, omgci2, omgp22, lmax,
     .      xkxsav, xkysav, xkzsav, nzfun, ibessel,
     .      ealphak, ebetak, ebk, capr,
     .      bxn, byn, bzn,
     .      uxx, uxy, uxz,
     .      uyx, uyy, uyz,
     .      uzx, uzy, uzz,
     .      icontxt,
     .      xx, yy, zz, isigma, xnuomg, psi, psilim,
     .      myid, nproc, delta0, gradprlb, bmod,
     .      xk0)


      if(eta3 .ne. 0.0)
     .   call current(xjpx3, xjpy3, xjpz3, nxmx, nymx, nphimx,
     .      nnodex, nnodey, nnodephi,
     .      nkx1, nkx2, nky1, nky2, nphi1, nphi2,
     .      nkdim1, nkdim2, mkdim1, mkdim2, lkdim1, lkdim2,
     .      xmi3, qi3, xn3a, xkti3, omgci3, omgp32, lmax,
     .      xkxsav, xkysav, xkzsav, nzfun, ibessel,
     .      ealphak, ebetak, ebk, capr,
     .      bxn, byn, bzn,
     .      uxx, uxy, uxz,
     .      uyx, uyy, uyz,
     .      uzx, uzy, uzz,
     .      icontxt,
     .      xx, yy, zz, isigma, xnuomg, psi, psilim,
     .      myid, nproc, delta0, gradprlb, bmod,
     .      xk0)


      time = second1(dummy) - t1
      tmin = time / 60.

      if (myid .eq. 0)then
         write(15,836) tmin
      end if




  836 format('time taken to calculate currents',f9.3,4h min)

      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey
               xjpx(i,j,k) = xjpxe(i,j,k) + xjpx1(i,j,k)
     .                     + xjpx2(i,j,k) + xjpx3(i,j,k)
               xjpy(i,j,k) = xjpye(i,j,k) + xjpy1(i,j,k)
     .                     + xjpy2(i,j,k) + xjpy3(i,j,k)
               xjpz(i,j,k) = xjpze(i,j,k) + xjpz1(i,j,k)
     .                     + xjpz2(i,j,k) + xjpz3(i,j,k)
            end do
         end do
      end do


      t1 = second1(dummy)

!     -----------------
!     allocate arrays
!     -----------------

      allocate( bqlavg_e(nuper, nupar, nnoderho) )
      allocate( cqlavg_e(nuper, nupar, nnoderho) )
      allocate( eqlavg_e(nuper, nupar, nnoderho) )
      allocate( fqlavg_e(nuper, nupar, nnoderho) )
      bqlavg_e = 0.0
      cqlavg_e = 0.0
      eqlavg_e = 0.0
      fqlavg_e = 0.0


      allocate( bqlavg_i1(nuper, nupar, nnoderho) )
      allocate( cqlavg_i1(nuper, nupar, nnoderho) )
      allocate( eqlavg_i1(nuper, nupar, nnoderho) )
      allocate( fqlavg_i1(nuper, nupar, nnoderho) )
      bqlavg_i1 = 0.0
      cqlavg_i1 = 0.0
      eqlavg_i1 = 0.0
      fqlavg_i1 = 0.0



      if(eta2 .ne. 0.0)then
         allocate( bqlavg_i2(nuper, nupar, nnoderho) )
         allocate( cqlavg_i2(nuper, nupar, nnoderho) )
         allocate( eqlavg_i2(nuper, nupar, nnoderho) )
         allocate( fqlavg_i2(nuper, nupar, nnoderho) )
         bqlavg_i2 = 0.0
         cqlavg_i2 = 0.0
         eqlavg_i2 = 0.0
         fqlavg_i2 = 0.0
      end if

      if(eta3 .ne. 0.0)then
         allocate( bqlavg_i3(nuper, nupar, nnoderho) )
         allocate( cqlavg_i3(nuper, nupar, nnoderho) )
         allocate( eqlavg_i3(nuper, nupar, nnoderho) )
         allocate( fqlavg_i3(nuper, nupar, nnoderho) )
         bqlavg_i3 = 0.0
         cqlavg_i3 = 0.0
         eqlavg_i3 = 0.0
         fqlavg_i3 = 0.0
      end if


      if (nzeta_wdot .gt. 0) then
        if (myid .eq. 0 .and. i_write .ne. 0) then
           open(unit=43,file='out_orbitrf.coef',  status='replace',
     .                                             form='formatted')
           write(43,309) nnodex, nnodey, nnodephi
           write(43,310) (capr(i), i = 1, nnodex)
           write(43,310) (y(j), j = 1, nnodey)
           write(43,310) (phi(k), k = 1, nnodephi)
           close(43)
        end if

      ! ------------------ !
      ! --  Electrons   -- !
      ! ------------------ !

      call ql_myra_write(bqlavg_e, cqlavg_e, eqlavg_e, fqlavg_e,
     .   wdote, fx0e, fy0e, fz0e, dvol, nrhomax, nnoderho, drho,
     .   dx, dy, dphi, nxmx, nymx, nphimx, nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, nkdim1, nkdim2,
     .   mkdim1, mkdim2, lkdim1, lkdim2,
     .   xme, qe, xkte, omgce, omgpe2, lmax,
     .   xkxsav, xkysav, xkzsav,
     .   ealphak, ebetak, ebk, capr,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   icontxt,
     .   xx, yy, zz, psi, psilim, nboundary,
     .   myid, nproc, gradprlb, bmod, ndiste, bmod_mid,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfdupere, dfdupare,
     .   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     .   vce_mks, dfe_cql_uprp, dfe_cql_uprl, rho, rho_a,
     .   rt, omgrf,
     .   lmaxdim, nzeta_wdot,
     .   dldbavg,
     .   1, i_write, xkx_cutoff, xky_cutoff, xkz_cutoff)



      ! ------------------ !
      ! --Majority ions -- !
      ! ------------------ !

      call ql_myra_write(bqlavg_i1, cqlavg_i1, eqlavg_i1, fqlavg_i1,
     .   wdoti1, fx0i1, fy0i1, fz0i1, dvol, nrhomax, nnoderho, drho,
     .   dx, dy, dphi, nxmx, nymx, nphimx, nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, nkdim1, nkdim2,
     .   mkdim1, mkdim2, lkdim1, lkdim2,
     .   xmi1, qi1, xkti, omgci1, omgp12, lmax,
     .   xkxsav, xkysav, xkzsav,
     .   ealphak, ebetak, ebk, capr,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   icontxt,
     .   xx, yy, zz, psi, psilim, nboundary,
     .   myid, nproc, gradprlb, bmod, ndisti1, bmod_mid,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper1, dfdupar1,
     .   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     .   vc1_mks, df1_cql_uprp, df1_cql_uprl, rho, rho_a,
     .   rt, omgrf,
     .   lmaxdim, nzeta_wdot,
     .   dldbavg,
     .   1, i_write, xkx_cutoff, xky_cutoff, xkz_cutoff)

      ! ------------------ !
      ! --Minority ions -- !
      ! ------------------ !

      if(eta2 .ne. 0.0)then
      call ql_myra_write(bqlavg_i2, cqlavg_i2, eqlavg_i2, fqlavg_i2,
     .   wdoti2, fx0i2, fy0i2, fz0i2, dvol, nrhomax, nnoderho, drho,
     .   dx, dy, dphi, nxmx, nymx, nphimx, nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, nkdim1, nkdim2,
     .   mkdim1, mkdim2, lkdim1, lkdim2,
     .   xmi2, qi2, xkti2, omgci2, omgp22, lmax,
     .   xkxsav, xkysav, xkzsav,
     .   ealphak, ebetak, ebk, capr,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   icontxt,
     .   xx, yy, zz, psi, psilim, nboundary,
     .   myid, nproc, gradprlb, bmod, ndisti2, bmod_mid,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper2, dfdupar2,
     .   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     .   vc2_mks, df2_cql_uprp, df2_cql_uprl, rho, rho_a,
     .   rt, omgrf,
     .   lmaxdim, nzeta_wdot,
     .   dldbavg,
     .   1, i_write, xkx_cutoff, xky_cutoff, xkz_cutoff)
      end if

      ! --------------------- !
      ! --Third ion species-- !
      ! --------------------- !

      if(eta3 .ne. 0.0)then
      call ql_myra_write(bqlavg_i3, cqlavg_i3, eqlavg_i3, fqlavg_i3,
     .   wdoti3, fx0i3, fy0i3, fz0i3, dvol, nrhomax, nnoderho, drho,
     .   dx, dy, dphi, nxmx, nymx, nphimx, nnodex, nnodey, nnodephi,
     .   nkx1, nkx2, nky1, nky2, nphi1, nphi2, nkdim1, nkdim2,
     .   mkdim1, mkdim2, lkdim1, lkdim2,
     .   xmi3, qi3, xkti3, omgci3, omgp32, lmax,
     .   xkxsav, xkysav, xkzsav,
     .   ealphak, ebetak, ebk, capr,
     .   uxx, uxy, uxz,
     .   uyx, uyy, uyz,
     .   uzx, uzy, uzz,
     .   icontxt,
     .   xx, yy, zz, psi, psilim, nboundary,
     .   myid, nproc, gradprlb, bmod, ndisti3, bmod_mid,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper3, dfdupar3,
     .   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     .   vc3_mks, df3_cql_uprp, df3_cql_uprl, rho, rho_a,
     .   rt, omgrf,
     .   lmaxdim, nzeta_wdot,
     .   dldbavg,
     .   1, i_write, xkx_cutoff, xky_cutoff, xkz_cutoff)
      end if


      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(15,2838) tmin
      endif

 2838 format('time to calculate quasilinear operator =', f9.3, ' min')

      t1 = second1(dummy)

      ! ------------------ !
      ! --  Electrons   -- !
      ! ------------------ !

      call wdot_qlcheck(wdote_ql,
     .   nnoderho, nrhomax,
     .   bqlavg_e, cqlavg_e, xme, omgrf, xkteavg, ndiste,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfdupere, dfdupare,
     .   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     .   vce_mks, dfe_cql_uprp, dfe_cql_uprl, rhon, rho_a,
     .   dldbavg)



      ! ------------------ !
      ! --Majority ions -- !
      ! ------------------ !

      call wdot_qlcheck(wdoti1_ql,
     .   nnoderho, nrhomax,
     .   bqlavg_i1, cqlavg_i1, xmi1, omgrf, xktiavg, ndisti1,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper1, dfdupar1,
     .   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     .   vc1_mks, df1_cql_uprp, df1_cql_uprl, rhon, rho_a,
     .   dldbavg)



      ! ------------------ !
      ! --Minority ions -- !
      ! ------------------ !

      if(eta2 .ne. 0.0)then

      call wdot_qlcheck(wdoti2_ql,
     .   nnoderho, nrhomax,
     .   bqlavg_i2, cqlavg_i2, xmi2, omgrf, xkti2avg, ndisti2,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper2, dfdupar2,
     .   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     .   vc2_mks, df2_cql_uprp, df2_cql_uprl, rhon, rho_a,
     .   dldbavg)

      end if



      ! --------------------- !
      ! --Third ion species-- !
      ! --------------------- !

      if(eta3 .ne. 0.0)then

      call wdot_qlcheck(wdoti3_ql,
     .   nnoderho, nrhomax,
     .   bqlavg_i3, cqlavg_i3, xmi3, omgrf, xkti3avg, ndisti3,
     .   nupar, nuper, n_psi,
     .   n_psi_dim, dfduper3, dfdupar3,
     .   UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work,
     .   vc3_mks, df3_cql_uprp, df3_cql_uprl, rhon, rho_a,
     .   dldbavg)

      end if

      end if

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(15, 2839) tmin
      endif

 2839 format('time to check the quasilinear operator =', f9.3, ' min')

      t1 = second1(dummy)

c
c--   Calculate E* dot J and Poynting flux (sp) in real space:
      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey
               wdot(i,j,k) = wdote(i,j,k)
     .                   + wdoti1(i,j,k) + wdoti2(i,j,k) + wdoti3(i,j,k)

               redotj(i,j,k) =
     .                     .5 * real(conjg(ealpha(i,j,k)) * xjpx(i,j,k)
     .                             + conjg(ebeta(i,j,k)) * xjpy(i,j,k)
     .                             + conjg(eb(i,j,k)) * xjpz(i,j,k))

               pc = -0.5 * (conjg(ex(i,j, k)) * xjx(i,j,k)
     1                    + conjg(ey(i,j, k)) * xjy(i,j,k)
     1                    + conjg(ez(i,j, k)) * xjz(i,j,k) )

               pcre(i,j,k) = real(pc)
               pcim(i,j,k) = aimag(pc)

            end do
         end do
      end do

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(15, 2840) tmin
      endif

 2840 format('time to calculate E* dot J and Poynting flux =', f9.3,
     . ' min')


      t1 = second1(dummy)


c
c--   Integrate E dot Jplasma and E dot Jantenna:
c

      call torint(x, y, phi, pcre, pcrto2,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

      call torint(x, y, phi, pcim, pcito2,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

      call torint(x, y, phi, redotj, ptot,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

      call torint(x, y, phi, wdot, ptot_wdot,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)


*     ------------------
*     Normalize to jdote:
*     ------------------
      powtot = ptot

      if (powtot .ne. 0.0) then
          pscale = prfin / powtot
      else
          pscale = 1.0
      endif

      if (prfin .eq. 0.0) pscale = 1.0
      if (pscale .le. 0.0) pscale = 1.0

      if (myid .eq. 0)then
         write(15, 1218) pscale
      end if

*     -----------------------
*     Scale powers to Prf,in:
*     -----------------------
      powtot = powtot * pscale
      ptot = ptot * pscale
      pcrto2 = pcrto2 * pscale
      pcito2 = pcito2 * pscale
      ptot_wdot = ptot_wdot * pscale

*     ------------------------------------
*     Scale fields and current to Prf, in:
*     ------------------------------------

      do n = nkx1, nkx2
         do m = nky1, nky2
            ealphakmod(n, m) = ealphakmod(n, m) * sqrt(pscale)
            ebetakmod(n, m) = ebetakmod(n, m) * sqrt(pscale)
            ebkmod(n, m) = ebkmod(n, m) * sqrt(pscale)
         end do
      end do


      do nphi = nphi1, nphi2
         ealphakmod1(nphi) = ealphakmod1(nphi) * sqrt(pscale)
         ebetakmod1(nphi) = ebetakmod1(nphi) * sqrt(pscale)
         ebkmod1(nphi) = ebkmod1(nphi) * sqrt(pscale)
      end do


      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey
               ex(i,j,k) = ex(i,j,k) * sqrt(pscale)
               ey(i,j,k) = ey(i,j,k) * sqrt(pscale)
               ez(i,j,k) = ez(i,j,k) * sqrt(pscale)


               ealpha(i,j,k) = ealpha(i,j,k) * sqrt(pscale)
               ebeta(i,j,k) = ebeta(i,j,k) * sqrt(pscale)
               eb(i,j,k) = eb(i,j,k) * sqrt(pscale)


               xjpxe(i,j,k) = xjpxe(i,j,k) * sqrt(pscale)
               xjpye(i,j,k) = xjpye(i,j,k) * sqrt(pscale)
               xjpze(i,j,k) = xjpze(i,j,k) * sqrt(pscale)

               xjpx1(i,j,k) = xjpx1(i,j,k) * sqrt(pscale)
               xjpy1(i,j,k) = xjpy1(i,j,k) * sqrt(pscale)
               xjpz1(i,j,k) = xjpz1(i,j,k) * sqrt(pscale)

               xjpx2(i,j,k) = xjpx2(i,j,k) * sqrt(pscale)
               xjpy2(i,j,k) = xjpy2(i,j,k) * sqrt(pscale)
               xjpz2(i,j,k) = xjpz2(i,j,k) * sqrt(pscale)

               xjpx3(i,j,k) = xjpx3(i,j,k) * sqrt(pscale)
               xjpy3(i,j,k) = xjpy3(i,j,k) * sqrt(pscale)
               xjpz3(i,j,k) = xjpz3(i,j,k) * sqrt(pscale)

               ntilda_e(i,j,k) = ntilda_e(i,j,k) * sqrt(pscale)

               redotj(i,j,k) = redotj(i,j,k) * pscale
               pcre(i,j,k) = pcre(i,j,k) * pscale
               pcim(i,j,k) = pcim(i,j,k) * pscale

               wdote(i,j,k)  = wdote(i,j,k)  * pscale
               wdoti1(i,j,k) = wdoti1(i,j,k) * pscale
               wdoti2(i,j,k) = wdoti2(i,j,k) * pscale
               wdoti3(i,j,k) = wdoti3(i,j,k) * pscale
               wdot(i,j,k)   = wdot(i,j,k)   * pscale

               fz0e(i,j,k)  = fz0e(i,j,k)  * pscale
               fz0i1(i,j,k) = fz0i1(i,j,k) * pscale
               fz0i2(i,j,k) = fz0i2(i,j,k) * pscale
               fz0i3(i,j,k) = fz0i3(i,j,k) * pscale
               fz0(i,j,k)   = fz0(i,j,k)   * pscale
            end do
         end do
      end do

      do n = 1, nnoderho
         wdote_ql(n)  = wdote_ql(n)  * pscale
         wdoti1_ql(n) = wdoti1_ql(n) * pscale
         wdoti2_ql(n) = wdoti2_ql(n) * pscale
         wdoti3_ql(n) = wdoti3_ql(n) * pscale

         do niu = 1, nuper
            do miu = 1, nupar

               bqlavg_e(niu, miu, n) = bqlavg_e(niu, miu, n) * pscale
               cqlavg_e(niu, miu, n) = cqlavg_e(niu, miu, n) * pscale
               eqlavg_e(niu, miu, n) = eqlavg_e(niu, miu, n) * pscale
               fqlavg_e(niu, miu, n) = fqlavg_e(niu, miu, n) * pscale


               bqlavg_i1(niu, miu, n) = bqlavg_i1(niu, miu, n) * pscale
               cqlavg_i1(niu, miu, n) = cqlavg_i1(niu, miu, n) * pscale
               eqlavg_i1(niu, miu, n) = eqlavg_i1(niu, miu, n) * pscale
               fqlavg_i1(niu, miu, n) = fqlavg_i1(niu, miu, n) * pscale

               if(eta2 .ne. 0.0)then
               bqlavg_i2(niu, miu, n) = bqlavg_i2(niu, miu, n) * pscale
               cqlavg_i2(niu, miu, n) = cqlavg_i2(niu, miu, n) * pscale
               eqlavg_i2(niu, miu, n) = eqlavg_i2(niu, miu, n) * pscale
               fqlavg_i2(niu, miu, n) = fqlavg_i2(niu, miu, n) * pscale
	       end if

	       if(eta3 .ne. 0.0)then
	       bqlavg_i3(niu, miu, n) = bqlavg_i3(niu, miu, n) * pscale
               cqlavg_i3(niu, miu, n) = cqlavg_i3(niu, miu, n) * pscale
               eqlavg_i3(niu, miu, n) = eqlavg_i3(niu, miu, n) * pscale
               fqlavg_i3(niu, miu, n) = fqlavg_i3(niu, miu, n) * pscale
	       end if


            end do
         end do

      end do

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(15, 2841) tmin
      endif

 2841 format('time to do power scaling =', f9.3, ' min')


      t1 = second1(dummy)
*     ------------------------------------------------------------------------
*     calculate individual species edotj`s and driven current on solution mesh
*     ------------------------------------------------------------------------

      do k = 1, nnodephi
         do i = 1, nnodex
            do j = 1, nnodey

               redotje(i,j,k)=
     .                     .5 * real(conjg(ealpha(i,j,k)) * xjpxe(i,j,k)
     .                             + conjg(ebeta(i,j,k)) * xjpye(i,j,k)
     .                             + conjg(eb(i,j,k)) * xjpze(i,j,k))

               redotj1(i,j,k)=
     .                     .5 * real(conjg(ealpha(i,j,k)) * xjpx1(i,j,k)
     .                             + conjg(ebeta(i,j,k)) * xjpy1(i,j,k)
     .                             + conjg(eb(i,j,k)) * xjpz1(i,j,k))

               redotj2(i,j,k)=
     .                     .5 * real(conjg(ealpha(i,j,k)) * xjpx2(i,j,k)
     .                             + conjg(ebeta(i,j,k)) * xjpy2(i,j,k)
     .                             + conjg(eb(i,j,k)) * xjpz2(i,j,k))

               redotj3(i,j,k)=
     .                     .5 * real(conjg(ealpha(i,j,k)) * xjpx3(i,j,k)
     .                             + conjg(ebeta(i,j,k)) * xjpy3(i,j,k)
     .                             + conjg(eb(i,j,k)) * xjpz3(i,j,k))

               redotjt(i,j,k) = redotje(i,j,k) +  redotj1(i,j,k)
     .                        + redotj2(i,j,k) +  redotj3(i,j,k)

            end do
         end do
      end do


*     ----------------------------------------------------------
*     Integrate individual species Estar dot J on solution mesh:
*     ----------------------------------------------------------

      call torint(x, y, phi, redotje, pedotje,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

      call torint(x, y, phi, redotj1, pedotj1,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

      call torint(x, y, phi, redotj2, pedotj2,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

      call torint(x, y, phi, redotj3, pedotj3,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

      call torint(x, y, phi, redotjt, pedotjt,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

      call torint(x, y, phi, wdote, pe,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)
      call torint(x, y, phi, wdoti1, pi1,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)
      call torint(x, y, phi, wdoti2, pi2,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)
      call torint(x, y, phi, wdoti3, pi3,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

      call torint(x, y, phi, wdot, pt,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)
      call torint(x, y, phi, fz0, fz0tot,
     .   1, nnodex, 1, nnodey, 1, nnodephi,
     .   capr, nxmx, nymx, nphimx)

*     -------------------------
*     Do flux surface averages:
*     -------------------------

      call fluxavg(redotje, redotjeavg, rho, nxmx, nymx, nphimx,
     .   nrhomax, nnodex, nnodey, nnodephi, nnoderho, drho,
     .   dx, dy, dphi, capr, fvol, vol)
      call fluxavg(redotj1, redotj1avg, rho, nxmx, nymx, nphimx,
     .   nrhomax, nnodex, nnodey, nnodephi, nnoderho, drho,
     .   dx, dy, dphi, capr, fvol, vol)
      call fluxavg(redotj2, redotj2avg, rho, nxmx, nymx, nphimx,
     .   nrhomax, nnodex, nnodey, nnodephi, nnoderho, drho,
     .   dx, dy, dphi, capr, fvol, vol)
      call fluxavg(redotj3, redotj3avg, rho, nxmx, nymx, nphimx,
     .   nrhomax, nnodex, nnodey, nnodephi, nnoderho, drho,
     .   dx, dy, dphi, capr, fvol, vol)

      call fluxavg(wdote, wdoteavg, rho, nxmx, nymx, nphimx,
     .   nrhomax, nnodex, nnodey, nnodephi, nnoderho, drho,
     .   dx, dy, dphi, capr, fvol, vol)
      call fluxavg(wdoti1, wdoti1avg, rho, nxmx, nymx, nphimx,
     .   nrhomax, nnodex, nnodey, nnodephi, nnoderho, drho,
     .   dx, dy, dphi, capr, fvol, vol)
      call fluxavg(wdoti2, wdoti2avg, rho, nxmx, nymx, nphimx,
     .   nrhomax, nnodex, nnodey, nnodephi, nnoderho, drho,
     .   dx, dy, dphi, capr, fvol, vol)
      call fluxavg(wdoti3, wdoti3avg, rho, nxmx, nymx, nphimx,
     .   nrhomax, nnodex, nnodey, nnodephi, nnoderho, drho,
     .   dx, dy, dphi, capr, fvol, vol)
      call fluxavg(wdot, wdotavg, rho, nxmx, nymx, nphimx,
     .   nrhomax, nnodex, nnodey, nnodephi, nnoderho, drho,
     .   dx, dy, dphi, capr, fvol, vol)



!     ----------------------------
!     Write to different outputs
!     ----------------------------

      if (myid .eq. 0)then
         write(15, 1309) ptot

         write(15, 1110) pcrto2, pcito2

         write(15, 1216)powtot
      end if

c--   Calculate species fractions from Estar dot J:
      pcedotje =  pedotje / pedotjt * 100.
      pcedotj1 =  pedotj1 / pedotjt * 100.
      pcedotj2 =  pedotj2 / pedotjt * 100.
      pcedotj3 =  pedotj3 / pedotjt * 100.

      pcedotjt = pcedotje + pcedotj1 + pcedotj2 + pcedotj3


71667 format(' Species absorption from Estar dot J:')


      if (myid .eq. 0)then
         write(15, 169)
         write(15, 71667)
         write(15, 1109)  pedotje, pcedotje,
     1                    pedotj1, pcedotj1,
     1                    pedotj2, pcedotj2,
     1                    pedotj3, pcedotj3,
     1                    pedotjt, pcedotjt


c--      Calculate species fractions from Wdot:
         if (pt .ne. 0.0) then
            pcte =   pe / pt  * 100.
            pcti1 =  pi1 / pt * 100.
            pcti2 =  pi2 / pt * 100.
            pcti3 =  pi3 / pt * 100.
         end if

         pctt = pcte + pcti1 + pcti2 + pcti3

81667 format(' Species absorption from Wdot:')

         write(15, 169)
         write(15, 81667)
         write(15, 71109)  pe, pcte,
     .                    pi1, pcti1,
     .                    pi2, pcti2,
     .                    pi3, pcti3,
     .                    pt,  pctt


*        --------------------------------
*        write psi and bmod data to out28
*        --------------------------------

         write(28, 309) mcap(imcap)
         write(28, 310) (((psi(i, j, k), i = 1, nnodex),j = 1,nnodey),
     .                                                  k = 1, nnodephi)
         write(28, 310) (((bmod(i, j, k), i = 1, nnodex),j = 1,nnodey),
     .                                                  k = 1, nnodephi)

*        --------------------------------
*        write out38
*        --------------------------------

         write(38, 309) nnodex, nnodey, nnodephi
         write(38, 310) psilim
         write(38, 310) (x(i), i = 1, nnodex)
         write(38, 310) (y(j), j = 1, nnodey)
         write(38, 310) (capr(i), i = 1, nnodex)

         write(38, 310) ((psi(i, j, kplot), i = 1, nnodex),j = 1,nnodey)
         write(38, 310) ((theta(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((bmod(i, j, kplot),i = 1, nnodex),j = 1,nnodey)

         write(38, 310) ((xnea(i, j, kplot),i = 1, nnodex),j = 1,nnodey)
         write(38, 310) ((xkte(i, j, kplot),i = 1, nnodex),j = 1,nnodey)
         write(38, 310) ((xkti(i, j, kplot), i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((xjy(i, j, kplot),  i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((xjx(i, j, kplot),  i = 1,nnodex),j = 1,nnodey)

         write(38, 309) nnoderho
         write(38, 310) (rhon(n), n = 1, nnoderho)
         write(38, 310) (xnavg(n), n = 1, nnoderho)


         write(38, 310) ((adisp(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((acold(i, j), i = 1, nnodex), j = 1, nnodey)


         write(38, 310) ((ex(i, j, kplot), i = 1, nnodex), j = 1,nnodey)
         write(38, 310) ((ey(i, j, kplot), i = 1, nnodex), j = 1,nnodey)
         write(38, 310) ((ez(i, j, kplot), i = 1, nnodex), j = 1,nnodey)

         write(38, 309) nkx1, nkx2
         write(38, 309) nky1, nky2

         write(38, 310) (xkxsav(n), n = nkx1, nkx2)
         write(38, 310) (xkysav(m), m = nky1, nky2)

         write(38, 310) ((ealphakmod(n, m), n = nkx1, nkx2),
     .      m = nky1, nky2)
         write(38, 310) ((ebetakmod(n, m), n = nkx1, nkx2),
     .      m = nky1, nky2)
         write(38, 310) ((ebkmod(n, m), n = nkx1, nkx2), m = nky1, nky2)

         write(38, 310) ((redotje(i,j,kplot),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((redotj1(i,j,kplot),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((redotj2(i,j,kplot),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((redotj3(i,j,kplot),i = 1,nnodex),j = 1,nnodey)


         kplot1 = 1
         kplot2 = kplot1 + ndkplot
         kplot3 = kplot2 + ndkplot
         kplot4 = kplot3 + ndkplot

         write(38, 310) ((redotj2(i,j,kplot1),i = 1,nnodex),j =1,nnodey)
         write(38, 310) ((redotj2(i,j,kplot2),i = 1,nnodex),j =1,nnodey)
         write(38, 310) ((redotj2(i,j,kplot3),i = 1,nnodex),j =1,nnodey)
         write(38, 310) ((redotj2(i,j,kplot4),i = 1,nnodex),j =1,nnodey)

         write(38, 310) ((ealpha(i, j, kplot1), i = 1, nnodex),
     .                                                     j = 1,nnodey)
         write(38, 310) ((ealpha(i, j, kplot2), i = 1, nnodex),
     .                                                     j = 1,nnodey)
         write(38, 310) ((ealpha(i, j, kplot3), i = 1, nnodex),
     .                                                     j = 1,nnodey)
         write(38, 310) ((ealpha(i, j, kplot4), i = 1, nnodex),
     .                                                     j = 1,nnodey)

         write(38, 310) ((psi(i, j, kplot1), i = 1, nnodex),j =1,nnodey)
         write(38, 310) ((psi(i, j, kplot2), i = 1, nnodex),j =1,nnodey)
         write(38, 310) ((psi(i, j, kplot3), i = 1, nnodex),j =1,nnodey)
         write(38, 310) ((psi(i, j, kplot4), i = 1, nnodex),j =1,nnodey)

         write(38, 309) nphi1, nphi2

         write(38, 310) (xkzsav(nphi),  nphi = nphi1, nphi2)
         write(38, 310) (ealphakmod1(nphi), nphi = nphi1, nphi2)
         write(38, 310) (ebetakmod1(nphi), nphi = nphi1, nphi2)
         write(38, 310) (ebkmod1(nphi), nphi = nphi1, nphi2)


         do nphi = nphi1, nphi2
            write(15, 901)nphi, xkzsav(nphi), ealphakmod1(nphi),
     .                  ebetakmod1(nphi), ebkmod1(nphi)
         end do


         write(38, 310) ((redotj1(i,j,kplot1),i = 1,nnodex),j =1,nnodey)
         write(38, 310) ((redotj1(i,j,kplot2),i = 1,nnodex),j =1,nnodey)
         write(38, 310) ((redotj1(i,j,kplot3),i = 1,nnodex),j =1,nnodey)
         write(38, 310) ((redotj1(i,j,kplot4),i = 1,nnodex),j =1,nnodey)

         write(38, 310) ((eb(i, j, kplot), i = 1, nnodex), j = 1,nnodey)

         write(38, 310) ((reson1(i,j,kplot1),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((reson2(i,j,kplot1),i = 1,nnodex),j = 1,nnodey)

         write(38, 310) ((reson1(i,j,kplot2),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((reson2(i,j,kplot2),i = 1,nnodex),j = 1,nnodey)

         write(38, 310) ((reson1(i,j,kplot3),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((reson2(i,j,kplot3),i = 1,nnodex),j = 1,nnodey)

         write(38, 310) ((reson1(i,j,kplot4),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((reson2(i,j,kplot4),i = 1,nnodex),j = 1,nnodey)

         write(38, 310) (((xjy(i,j,k),i = 1,nnodex), j = 1,nnodey),
     .                   k = 1, nnodephi)
         write(38, 309) kplot1, kplot2, kplot3, kplot4

         write(38, 310) ((bmod(i, j, kplot1), i = 1,nnodex),j =1,nnodey)
         write(38, 310) ((bmod(i, j, kplot2), i = 1,nnodex),j =1,nnodey)
         write(38, 310) ((bmod(i, j, kplot3), i = 1,nnodex),j =1,nnodey)
         write(38, 310) ((bmod(i, j, kplot4), i = 1,nnodex),j =1,nnodey)

         write(38, 310) ((ebeta(i, j, kplot1), i = 1, nnodex),
     .                                                     j = 1,nnodey)
         write(38, 310) ((ebeta(i, j, kplot2), i = 1, nnodex),
     .                                                     j = 1,nnodey)
         write(38, 310) ((ebeta(i, j, kplot3), i = 1, nnodex),
     .                                                     j = 1,nnodey)
         write(38, 310) ((ebeta(i, j, kplot4), i = 1, nnodex),
     .                                                     j = 1,nnodey)

         write(38, 310) ((eb(i, j, kplot1), i = 1, nnodex),j = 1,nnodey)
         write(38, 310) ((eb(i, j, kplot2), i = 1, nnodex),j = 1,nnodey)
         write(38, 310) ((eb(i, j, kplot3), i = 1, nnodex),j = 1,nnodey)
         write(38, 310) ((eb(i, j, kplot4), i = 1, nnodex),j = 1,nnodey)

         write(38, 310) ((wdote(i,j,kplot),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((wdoti1(i,j,kplot),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((wdoti2(i,j,kplot),i = 1,nnodex),j = 1,nnodey)
         write(38, 310) ((wdoti3(i,j,kplot),i = 1,nnodex),j = 1,nnodey)


*        -------------------------------
*        write wdot data to movie_wdot
*        -------------------------------
         call write_power_abs(69, nnodex, nnodey, nnodephi,
     .                        wdote, wdoti1, wdoti2, wdoti3)

*        -------------------------------
*        write edotj data to edotj.out
*        -------------------------------
         call write_power_abs(39, nnodex, nnodey, nnodephi,
     .                        redotje, redotj1, redotj2, redotj3)

*        -------------------
*        write fpm (34) file
*        -------------------
         if (imcap .eq. 1) then
            open(unit=34, file='fpm',status='unknown', form='formatted')
            write(34, 309) mcap_number

            write(34, 309) nnodex, nnodey, nnodephi
            write(34, 310) (x(i), i = 1, nnodex)
            write(34, 310) (y(j), j = 1, nnodey)
            write(34, 310) (capr(i), i = 1, nnodex)
            write(34, 310) (phiprime(k), k = 1, nnodephi)

            write(34, 310) psilim, rt, b0
            write(34, 310) (((psi(i, j, k), i = 1, nnodex),
     .           j = 1,nnodey), k = 1, nnodephi)
            write(34, 310) (((bmod(i, j, k), i = 1, nnodex),
     .           j = 1,nnodey), k = 1, nnodephi)

            write(34, 309) nnoderho
            write(34, 310) (rhon(n), n = 1, nnoderho)
            write(34, 310) (xnavg(n), n = 1, nnoderho)
         endif

         write(34, 309) mcap(imcap)

         write(34, 310)(((ealpha(i,j,k),i=1,nnodex),
     &                 j=1,nnodey),k=1,nnodephi)
         write(34, 310)(((ebeta(i,j,k),i=1,nnodex),
     &                 j=1,nnodey),k=1,nnodephi)
         write(34, 310)(((eb(i,j,k),i=1,nnodex),
     &                 j=1,nnodey),k=1,nnodephi)

         write(34, 310)(((ntilda_e(i,j,k),i=1,nnodex),
     &                  j=1,nnodey),k=1,nnodephi)

         write(34, 310) ptot, pcrto2, pcito2

         write(34, 310) (redotjeavg(n), n = 1, nnoderho)
         write(34, 310) (redotj1avg(n), n = 1, nnoderho)
         write(34, 310) (redotj2avg(n), n = 1, nnoderho)
         write(34, 310) (redotj3avg(n), n = 1, nnoderho)

         write(34, 310) (wdoteavg(n), n = 1, nnoderho)
         write(34, 310) (wdoti1avg(n), n = 1, nnoderho)
         write(34, 310) (wdoti2avg(n), n = 1, nnoderho)
         write(34, 310) (wdoti3avg(n), n = 1, nnoderho)

         write(34, 310) (wdote_ql(n),  n = 1, nnoderho)
         write(34, 310) (wdoti1_ql(n), n = 1, nnoderho)
         write(34, 310) (wdoti2_ql(n), n = 1, nnoderho)
         write(34, 310) (wdoti3_ql(n), n = 1, nnoderho)
      end if

      time = second1(dummy)-time0
      ttotal = time / 60.

      if(myid .eq. 0)then
         write(15, 162)
         write(15, 899) ttotal
      end if

!     --------------------
!     format declarations
!     --------------------

  310 format(1p6e12.4)
  309 format(10i10)
 9310 format(1p7e12.4)
  311 format(1p10e12.4)
 3117 format(1p11e12.4)
 1312 format(i10, 1p8e12.4)
 1314 format (2i10, 1p8e12.4)
 1310 format (5i10, 1p8e12.4)
13149 format (4i10, 1p8e12.4)
 1313 format(10i10)
 1311 format(1p9e12.4)
   10 format(i10,1p4e10.3,i10,1pe10.3)
 1010 format(1f4.0,4f8.3,3f7.3,1f8.3,1f9.3,2f8.3)
 1009 format(3x,21h        frequency  = ,1pe12.4,7h hertz )
 1012 format(3x,21h             omgrf = ,1pe12.4,7h hertz )
 2012 format(3x,51h2*pi*rt/Ly * real part of impedance (resistance) = ,
     1   1pe12.4,10h ohms     )
 2013 format(3x,
     1   55h2*pi*rt/Ly * imaginary part of impedance (reactance) = ,
     1   1pe12.4,10h ohms     )
 1014 format(3x,21h               xkz = ,1pe12.4,7h m-1   )
 1321 format(3x,21h         vph / vth = ,1pe12.4,7h       )
 1391 format(3x,21h critical shear(0) = ,1pe12.4,7h s-1   )
 1392 format(3x,21h critical shear(a) = ,1pe12.4,7h s-1   )
 1393 format(3x,21h         mu neo(a) = ,1pe12.4,7h s-1   )
 1322 format(3x,21h               vph = ,1pe12.4,7h       )
 1323 format(3x,21h               vth = ,1pe12.4,7h       )
 1714 format(3x,21h               xk0 = ,1pe12.4,7h m-1   )
 1021 format(3x,21h        n parallel = ,1pe12.4,7h       )
 1022 format(3x,21h           rhonorm = ,1pe12.4,7h m     )
 1023 format(3x,21h              epsl = ,1pe12.4,7h T     )
91023 format(3x,21h              capa = ,1pe12.4,7h       )
 1812 format(3x,21h                rt = ,1pe12.4,7h m     )
 1822 format(3x,21h            aplasm = ,1pe12.4,7h m     )
 1823 format(3x,21h              rant = ,1pe12.4,7h m     )
 1809 format(3x,21h                b0 = ,1pe12.4,7h T     )
 1813 format(3x,21h              xn10 = ,1pe12.4,7h m-3   )
 6813 format(3x,21h              xne0 = ,1pe12.4,7h m-3   )
 1814 format(3x,21h              xn20 = ,1pe12.4,7h m-3   )
 1834 format(3x,21h              xn30 = ,1pe12.4,7h m-3   )
 6834 format(3x,21h              eta1 = ,1pe12.4,7h       )
 6835 format(3x,21h              eta2 = ,1pe12.4,7h       )
 6836 format(3x,21h              eta3 = ,1pe12.4,7h       )
 1815 format(3x,21h               te0 = ,1pe12.4,7h eV    )
 1821 format(3x,21h               ti0 = ,1pe12.4,7h eV    )
 1016 format(3x,21h xnue/omgrf ad hoc = ,1pe12.4,7h       )
 1017 format(3x,21h xnu1/omgrf ad hoc = ,1pe12.4,7h       )
 1018 format(3x,21h xnu2/omgrf ad hoc = ,1pe12.4,7h       )
 1013 format(3x,21h              nphi = ,i12,7h       )
 7013 format(3x,21h           nmodesx = ,i12,7h       /
     1       3x,21h           nmodesy = ,i12,7h       /
     1       3x,21h         nmodesphi = ,i12,7h       )

 7113 format(3x,21h             nwdot = ,i12,7h       )
 7114 format(3x,21h              mcap = ,i12,7h       )
 7213 format(3x,21h           nnodecx = ,i12,7h       /
     1       3x,21h           nnodecy = ,i12,7h       )
 7014 format(3x,21h              lmax = ,i12,7h       )
 7015 format(3x,21h           ibessel = ,i12,7h       )
 7115 format(3x,21h             nzfun = ,i12,7h       )
 7016 format(3x,21h         rhoi1 / L = ,1pe12.4,7h       )
 7116 format(3x,21h         xnu / omg = ,1pe12.4,7h       )
 7017 format(3x,21h             rhoi1 = ,1pe12.4,7h m     )
 7217 format(3x,21h             qavg0 = ,1pe12.4,7h       )
 1020 format(3x,21h            nnodex = ,i12, 7h       /
     1       3x,21h            nnodey = ,i12, 7h       )

 3013 format(3x,21h                i0 = ,i12,7h       )
30131 format(3x,21h             ileft = ,i12,7h       )
30132 format(3x,21h            iright = ,i12,7h       )
 3017 format(3x,21h           xnup(0) = ,1pe12.4,7h s-1   )
 3018 format(3x,21h          eps**1.5 = ,1pe12.4,7h       )

 1309 format(3x,
     1   35h             total power absorbed = ,1pe12.4,9h watts   )
11091 format(3x,
     1   35h       total power absorbed (old) = ,1pe12.4,9h watts   )

 1109 format(
     .   3x, 35h      power absorbed by electrons = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       /
     .   3x, 35h  power absorbed by majority ions = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       /
     .   3x, 35h  power absorbed by minority ions = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       /
     .   3x, 35hpower absorbed by 3rd ion species = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       /
     .   3x, 35h             total power absorbed = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       )

71109 format(
     .   3x, 35h      power absorbed by electrons = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       /
     .   3x, 35h  power absorbed by majority ions = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       /
     .   3x, 35h  power absorbed by minority ions = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       /
     .   3x, 35hpower absorbed by 3rd ion species = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       /
     .   3x, 35h             total power absorbed = ,1e12.5,
     .   7h Watts , 3h = , f10.4, 9h %       )

81109 format(
     .   3x, 35h  power absorbed by majority ions = ,1e12.5,
     .   9h watts   /
     1   3x, 35h  power absorbed by minority ions = ,1e12.5,
     1   9h watts   /
     1   3x, 35h      power absorbed by electrons = ,1e12.5,
     1   9h watts   )

 1112 format(
     1   3x,35h      power absorbed by electrons = ,1f12.4,9h %       /
     1   3x,35h  power absorbed by majority ions = ,1f12.4,9h %       /
     1   3x,35h  power absorbed by minority ions = ,1f12.4,9h %       /
     1   3x,35hpower absorbed by 3rd ion species = ,1f12.4,9h %       /
     1   3x,35h             total power absorbed = ,1f12.4,9h %       )
 1289 format(
     1   3x,36h                         x(iedge) = ,1pe12.4,2h m/
     1   3x,36h                               bz = ,1pe12.4,2h T/
     1   3x,36h              real(eps parallel)  = ,1pe12.4,2h  /
     1   3x,36h    estimate: real(eps parallel)  = ,1pe12.4,2h  /
     1   3x,36h                 real (eps left)  = ,1pe12.4,2h  /
     1   3x,36h     estimate:   real (eps left)  = ,1pe12.4,2h  /
     1   3x,36h                real (eps right)  = ,1pe12.4,2h  /
     1   3x,36h    estimate:   real (eps right)  = ,1pe12.4,2h  /
     1   3x,36h                               xn = ,1pe12.4,4h m-3/
     1   3x,36h                               LN = ,1pe12.4,2h m/
     1   3x,36h                        mod2 (ez) = ,1pe12.4,9h (V/m)**2/
     1   3x,36h                            L rf  = ,1pe12.4,2h m/
     1   3x,36h                x=omgrf/omgci(i)  = ,1pe12.4,4h    /
     1   3x,36h                             xkti = ,1pe12.4,4h deg/
     1   3x,36h                             xkte = ,1pe12.4,4h deg)
 1389 format(
     1   3x,36h                           dlnndr = ,1pe12.4,4h m-1/
     1   3x,36h                            deprl = ,1pe12.4,
     1   11h (V/m)**2/m/
     1   3x,36h                           deleft = ,1pe12.4,
     1   11h (V/m)**2/m/
     1   3x,36h                           derght = ,1pe12.4,
     1   11h (V/m)**2/m/
     1   3x,36h                        real (ez) = ,1pe12.4,4h V/m/
     1   3x,36h                         imag(ez) = ,1pe12.4,4h V/m)
 1290 format(
     1   3x,36h             alpha rf parallel    = ,1pe12.4,2h  /
     1   3x,36h estimate of alpha rf parallel    = ,1pe12.4,2h  /
     1   3x,36h                 alpha rf left    = ,1pe12.4,2h  /
     1   3x,36h     estimate of alpha rf left    = ,1pe12.4,2h  /
     1   3x,36h                 alpha rf right   = ,1pe12.4,2h  /
     1   3x,36h     estimate of alpha rf right   = ,1pe12.4,2h  /
     1   3x,36h                 alpha rf total   = ,1pe12.4,2h  /
     1   3x,36h electron ponderomotive potential = ,1pe12.4,3h eV/
     1   3x,36h      ion ponderomotive potential = ,1pe12.4,3h eV)
 1110 format(
     1   3x,35h               real power emitted = ,1pe12.4,9h watts   /
     1   3x,35h          imaginary power emitted = ,1pe12.4,9h watts   )
 1113 format(3x,'driven current per meter = ',1pe12.4,' Amps/m')
 1114 format(3x,'current driven per watt = ',1pe12.4,' A/W')
 1215 format(3x,'current drive efficiency = ',1pe12.4,' A/W/m**2')
 1216 format(3x,'total RF power = ',1pe12.4,' Watts')
 1217 format(3x,'total driven current = ',1pe12.4,' Amps')
 1218 format(3x,'pscale = ',1pe12.4,' ')
 1219 format(3x,'gamma = ',1pe12.4,' ')
 1111 format(3x,8h At x = ,1pe12.4,2h m,
     1       3x,9h ifail = ,i5)
  162 format(1h1)
  169 format(1h )
  163 format(1h0)
 2162 format(1p8e12.4)
 2163 format(i5, 1p11e12.3)
 2165 format(3i10,5e10.3)
 1002 format(2i10,7e10.3)
11313 format(2i10, 1p8e12.4)
 1000 format(1i10,7e10.3)
 2000 format(1i10,1e10.3,1i10,1e10.3)
 1001 format(8e10.3)


  164 format(1h ,1x,4h  i ,  3x,6h R(x) ,
     1           6x,6h  x   ,6x,6h bmod ,6x,6hre om1,
     1           6x,6hre om2,6x,6h acold,6x,6h      , 6x,6h      )


 7164 format(1h ,3x,6h     k,3x,6h  phi ,
     1           6x,6hphiprm,6x,6h phi0 ,6x,6h      ,
     1           6x,6h      ,6x,6h      ,6x,6h      , 6x,6h      )

 8164 format(1h ,3x,6h     k,3x,6h  Jy  ,
     1           6x,6h      ,6x,6h      ,6x,6h      ,
     1           6x,6h      ,6x,6h      ,6x,6h      , 6x,6h      )

 8165 format(1h ,3x,6h     k,3x,8hshapephi,
     1           6x,6h      ,6x,6h      ,6x,6h      ,
     1           6x,6h      ,6x,6h      ,6x,6h      , 6x,6h      )


 9164 format(1h ,1x,4h  j ,  3x,6h R(x) ,
     1           6x,6h  y   ,6x,6h bmod ,6x,6hre om1,
     1           6x,6hre om2,6x,6h acold,6x,6h      , 6x,6h      )

 6164 format(1h ,3x,6h  nphi,3x,6hxkzsav,
     1           6x,6h      ,6x,6h      ,6x,6h      ,
     1           6x,6h      ,6x,6h      ,6x,6h      , 6x,6h      )



 1164 format(1h ,3x,6h R(x) ,

     1           6x,6h  x   , 6x,6hvymean,6x,6hmu*vyx,
     1           5x,7hmu*vzxx,6x,6hmu*vzx,
     1           6x,6h      ,6x,6h      ,6x,6h      ,
     1           6x,6h      )

 4000 format(1h ,3x,7h  x    ,5x,7hre k**2,5x,7hre k**2,
     1           5x,7hre k**2,5x,7hre k**2,5x,7hre k**2,
     1           5x,7hre k**2,5x,7hre k**2,5x,7hre k**2)

 7000 format(1h ,3x,6h  x   ,6x,6h amp1 ,6x,6h amp2 ,
     1           6x,6h amp3 ,6x,6h amp4 ,6x,6h amp5 ,
     1           6x,6h amp6 ,6x,6h amp7 ,6x,6h amp8 )

 7100 format(1h ,3x,6h  x   ,6x,6h  sx1 ,6x,6h  sx2 ,
     1           6x,6h  sx3 ,6x,6h  sx4 ,6x,6h  sx5 ,
     1           6x,6h  sx6 ,6x,6h  sx7 ,6x,6h  sx8 ,
     1           6x,6hsx sum)
 4002 format(1h ,3x,6h  x   ,6x,6hre dt1,6x,6hre dt2,
     1           6x,6hre dt3,6x,6hre dt4,6x,6hre dt5,
     1           6x,6hre dt6,6x,6hre dt7,6x,6hre dt8)

 4001 format(1h ,3x,6h  x   , 5x,7him k**2,5x,7him k**2,
     1           5x,7him k**2,5x,7him k**2,5x,7him k**2,
     1           5x,7him k**2,5x,7him k**2,5x,7him k**2)

 4003 format(1h ,3x,6h  x   ,6x,6him dt1,6x,6him dt2,
     1           6x,6him dt3,6x,6him dt4,6x,6him dt5,
     1           6x,6him dt6,6x,6him dt7,6x,6him dt8)


 1650 format(3x,32hdispersion relation             )
 7650 format(3x,32hamplitude of modes              )
 7750 format(3x,32hflux of power carried by modes  )
 1651 format(3x,32hdispersion relation check       )
 1115 format(3x,13h        kh = ,1pe12.4,7h m-1   )
 1116 format(3x,13h        lp = ,1pe12.4,7h m     )
 1117 format(3x,13h      eps2 = ,1pe12.4,7h tesla )
 1118 format(3x,13havg iotabr = ,1pe12.4,7h       )
 1119 format(2x,14hellipticity = ,1pe12.4,7h       )
 1120 format(3x,13h    psib   = ,1pe12.4,8h tesla-m)



  899 format('total cpu time used =',f9.3,4h min)
  898 format('total cpu time used for all mcaps =',f9.3,4h min)


 9930 format(3x, 'electron power and flow')
  933 format(3x, 'species #3 ion power and flow')
  932 format(3x, 'minority ion power and flow')
 1214 format(9x, 'i', 6x, 'R', 9x, 'wdoti2', 7x, 'fyi2')
  931 format(3x, 'majority ion power and flow')
 1319 format(9x, 'i', 6x, 'R', 9x, 'wdoti1', 7x, 'fyi1')

 9321 format(3x, 'Electric field:')
 1394 format(9x, 'i', 6x, 'R', 9x, 'Re Ex', 7x, 'Im Ex',
     1                         7x, 'Re Ey', 7x, 'Im Ey',
     1                         7x, 'Re Ez', 7x, 'Im Ez')

      call blacs_barrier(icontxt, 'All')

!     -----------------
!     deallocate arrays
!     -----------------


      deallocate( bqlavg_e )
      deallocate( cqlavg_e )
      deallocate( eqlavg_e )
      deallocate( fqlavg_e )

      deallocate( bqlavg_i1 )
      deallocate( cqlavg_i1 )
      deallocate( eqlavg_i1 )
      deallocate( fqlavg_i1 )

      if(eta2 .ne. 0.0) then
         deallocate( bqlavg_i2 )
         deallocate( cqlavg_i2 )
         deallocate( eqlavg_i2 )
         deallocate( fqlavg_i2 )
      end if

      if(eta3 .ne. 0.0) then
      deallocate( bqlavg_i3 )
         deallocate( cqlavg_i3 )
         deallocate( eqlavg_i3 )
         deallocate( fqlavg_i3 )
      end if


      deallocate( inv_xx )
      deallocate( inv_yy )
      deallocate( inv_zz )
      deallocate( inv_xx_t )
      deallocate( inv_yy_t )
      deallocate( inv_zz_t )
      deallocate( new_to_org )
      deallocate( niabegin_all )
      deallocate( isize_all )
      deallocate( row )
      deallocate( rowk )
      deallocate( Btmp )
      deallocate( brhs )
      deallocate( descBtmp_all )
      deallocate( descbrhs_all )
      deallocate( itable )
      deallocate( jtable )
      deallocate( ktable )
      deallocate( mtable )
      deallocate( ntable )
      deallocate( nphitable )
      deallocate( brhs2 )
      deallocate( brhs_tmp )

      deallocate( p_brhs )
      deallocate( p_ipiv )

 9000 continue  ! end of mcap loop

      deallocate( UPERP )
      deallocate( UPARA )
      deallocate( VPERP )
      deallocate( VPARA )
      deallocate( UPERP_work )
      deallocate( UPARA_work )

      deallocate( f )

      deallocate( dfdupere )
      deallocate( dfdupare )

      deallocate( dfduper1 )
      deallocate( dfdupar1 )
      deallocate( dfduper2 )
      deallocate( dfdupar2 )
      deallocate( dfduper3 )
      deallocate( dfdupar3 )

      if(ndisti1   .ne. 0 .or.
     .   ndisti2   .ne. 0 .or.
     .   ndisti3   .ne. 0 .or.
     .   ndiste    .ne. 0)   then

         deallocate( fe_cql_cart )
         deallocate( f1_cql_cart )
         deallocate( f2_cql_cart )
         deallocate( f3_cql_cart )

         deallocate( dfe_cql_uprp )
         deallocate( dfe_cql_uprl )

         deallocate( df1_cql_uprp )
         deallocate( df1_cql_uprl )

         deallocate( df2_cql_uprp )
         deallocate( df2_cql_uprl )

         deallocate( df3_cql_uprp )
         deallocate( df3_cql_uprl )
      end if

      deallocate( mask )

      if (myid .eq. 0) then
         close(28)
         close(34)
         close(38)
         close(39)
         close(53)
         close(54)
         close(69)
      endif

      deallocate( xjx, xjy, xjz )

      deallocate( psi, thetap, rho, xkte, xkti, xkti2, xkti3,
     .     xnea, xn1a, xn2a, xn3a,
     .     omgce, omgci1, omgci2, omgci3, omgpe2, omgp12, omgp22,omgp32,
     .     bmod, bmod_mid)


      deallocate( xjpxe, xjpye, xjpze, xjpx1, xjpy1, xjpz1, xjpx2,
     .     xjpy2, xjpz2, xjpx3, xjpy3, xjpz3, xjpx, xjpy, xjpz )

      deallocate( pcre, pcim, redotj, redotje, redotj1, redotj2,
     .     redotj3, redotjt )


      deallocate( bxn, byn, bzn)

      deallocate( uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz )

      deallocate( dxuxx, dxuxy, dxuxz,
     .     dxuyx, dxuyy, dxuyz,
     .     dxuzx, dxuzy, dxuzz )

      deallocate( dyuxx, dyuxy, dyuxz,
     .     dyuyx, dyuyy, dyuyz,
     .     dyuzx, dyuzy, dyuzz )


      deallocate( dyyuxx, dyyuxy, dyyuxz,
     .     dyyuyx, dyyuyy, dyyuyz,
     .     dyyuzx, dyyuzy, dyyuzz )

      deallocate( dxyuxx, dxyuxy, dxyuxz,
     .     dxyuyx, dxyuyy, dxyuyz,
     .     dxyuzx, dxyuzy, dxyuzz )

      deallocate( dxxuxx, dxxuxy, dxxuxz,
     .     dxxuyx, dxxuyy, dxxuyz,
     .     dxxuzx, dxxuzy, dxxuzz )

      deallocate( dphiuxx, dphiuxy, dphiuxz,
     .     dphiuyx, dphiuyy, dphiuyz,
     .     dphiuzx, dphiuzy, dphiuzz )


      deallocate( d2phiuxx,d2phiuxy, d2phiuxz,
     .     d2phiuyx, d2phiuyy, d2phiuyz,
     .     d2phiuzx, d2phiuzy, d2phiuzz )

      deallocate( dxphiuxx,dxphiuxy, dxphiuxz,
     .     dxphiuyx, dxphiuyy, dxphiuyz,
     .     dxphiuzx, dxphiuzy, dxphiuzz )

      deallocate( dyphiuxx,dyphiuxy, dyphiuxz,
     .     dyphiuyx, dyphiuyy, dyphiuyz,
     .     dyphiuzx, dyphiuzy, dyphiuzz )

      deallocate( gradprlb )


      deallocate( ealphak, ebetak, ebk )

      deallocate( ex, ey, ez )

      deallocate( ealpha, ebeta, eb )

      deallocate( ntilda_e )


      deallocate( xx, yy, zz )

      deallocate( reson1, reson2 )

      deallocate( wdote, wdoti1, wdoti2, wdoti3, wdot )

      deallocate( fx0e, fx0i1, fx0i2, fx0i3, fy0e, fy0i1, fy0i2, fy0i3,
     .     fz0e, fz0i1, fz0i2, fz0i3, fx0, fy0, fz0 )


      time = second1(dummy)-time00
      ttotal = time / 60.

      if(myid .eq. 0)then
         write(15, 162)
         write(15, 899) ttotal
         close(15)
      end if

*     -------------------------
*     stop parallel environment
*     -------------------------
      call blacs_gridexit( icontxt )
      call blacs_exit(0)

      end
c
c***************************************************************************
c
        subroutine pzgecopy( m,n, A,ia,ja,descA, B,ib,jb,descB)
        integer m,n,  ia,ja,descA(*), ib,jb,descB(*)
        complex A(*), B(*)

        complex alpha,beta

        alpha = 1.0
        beta = 0.0
*       ----------------------------------------------------
*       pzgeadd is new capability found in PBLAS V2 or later
*       otherwise, it can be implmented less efficiently by
*       repeated calls to pzcopy
*       ----------------------------------------------------
        call pzgeadd( 'N', m,n, alpha, A,ia,ja,descA,
     &        beta, B,ib,jb,descB )
        return
        end subroutine pzgecopy

c
c***************************************************************************
c
      subroutine fluxavg(f, favg, rho, nxdim, nydim, nphidim, nrhodim,
     .   nnodex, nnodey, nnodephi, nnoderho, drho, dx, dy, dphi, capr,
     .   fvol, vol)

      implicit none

      integer nxdim, nydim, nphidim, nrhodim, nnodex, nnodey, nnodephi,
     .   nnoderho
      integer n, i, j, k

      real f(nxdim, nydim, nphidim), favg(nrhodim),
     .    rho(nxdim, nydim, nphidim),
     .   drho, dx, dy, dphi, fvol(nrhodim), vol(nrhodim), capr(nxdim)


      do n = 1, nnoderho
          fvol(n) = 0.0
          vol(n) = 0.0
      end do


      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi
               n = int(rho(i,j,k) / drho) + 1
               if(n .ge. 1 .and. n .le. nnoderho)then

                  fvol(n) = fvol(n) + dx * dy * capr(i) * dphi *f(i,j,k)
                  vol(n) =   vol(n) + dx * dy * capr(i) * dphi

               end if
            end do
         end do
      end do


      do n = 1, nnoderho
         if(vol(n) .ne. 0.0)
     .     favg(n) = fvol(n) / vol(n)
c         write(6, 100)n, fvol(n), vol(n), favg(n)
      end do



  100 format (1i10, 1p8e12.4)
  101 format (10i10)
      return
      end
c
c***************************************************************************
c
      subroutine flux_to_rzphi(nnodex, nnodey, nnodephi, profile_in,
     .   profile_out, rho, nrho, rho_ijk)
      ! takes profile_in(rho) and returns profile_out(nx, ny, nphi)
      ! Used to transform the &state namelist arrays

      implicit none

      real :: profile_in(nrho), profile_out(nnodex,nnodey,nnodephi)
      real :: rho(nrho), rho_ijk(nnodex,nnodey,nnodephi)
      integer :: nnodex, nnodey, nnodephi, nrho, i, j, k, n

      do i = 1, nnodex
         do j = 1, nnodey
            do k = 1, nnodephi
               do n = 1, nrho
                  if (rho(n) .gt. rho_ijk(i,j,k)) exit
               end do

               if (n .eq. 1) then
                  profile_out(i,j,k) = profile_in(1)
               else if (n .gt. nrho) then
                  profile_out(i,j,k) = profile_in(nrho)
               else
                  profile_out(i,j,k) = profile_in(n-1)  +
     .                (profile_in(n) - profile_in(n-1)) /
     .                (rho(n) - rho(n-1)) * (rho_ijk(i,j,k) - rho(n-1))
               end if
           end do
        end do
      end do

      return
      end subroutine flux_to_rzphi


      subroutine write_power_abs(nunit, nx, ny, nphi, pe, p1, p2, p3)
      ! Write power absorved for different species
      implicit none

      integer :: i, j, k, nx, ny, nphi, nunit
      real, dimension(nx, ny, nphi) :: pe, p1, p2, p3

  310 format(1p6e12.4)
  309 format(10i10)

         write(nunit, 310) (((pe(i,j,k), i=1, nx), j=1, ny), k=1, nphi)
         write(nunit, 310) (((p1(i,j,k), i=1, nx), j=1, ny), k=1, nphi)
         write(nunit, 310) (((p2(i,j,k), i=1, nx), j=1, ny), k=1, nphi)
         write(nunit, 310) (((p3(i,j,k), i=1, nx), j=1, ny), k=1, nphi)

      end subroutine write_power_abs
