module common_grid
use configure

  implicit none

  ! Parameters
 ! integer, parameter :: FRSIZE_MAX = 5000
 ! integer, parameter :: FRSIZE_FREQ = 1000
 ! integer, parameter :: FRSIZE_MU_HALF = 300
 ! integer, parameter :: FRSIZE_PHI = 400
 ! integer, parameter :: FRSIZE_DRR = 10
 ! integer, parameter :: FRSIZE_X = 500
 ! integer, parameter :: FRSIZE_Y = 500

  ! >>>> Spatial grid <<<<
  integer :: irsi_place(0:3)
  integer :: irsi_frsizex, irsi_frsizey
  integer :: irsi_inp, irsi_ibnp

  double precision :: rsi_x_c(-1:FRSIZE_MAX+2, 1:2)
  double precision :: rsi_dx_c(-1:FRSIZE_MAX+2, 1:2)
  double precision :: rsi_x_i(-1:FRSIZE_MAX+2, 1:2)
  double precision :: rsi_dx_i(-1:FRSIZE_MAX+2, 1:2)
  double precision :: rsi_dx_cil(-1:FRSIZE_MAX+2, 1:2)
  double precision :: rsi_dx_cir(-1:FRSIZE_MAX+2, 1:2)

  ! >>>> Frequency grid <<<<
  double precision :: freq_nu(1:FRSIZE_FREQ)
  integer :: freq_nr, freq_grid_type

  ! >>>> Global spacegrid info <<<<
  character :: spacegrid_metric
  character :: spacegrid_type_x, spacegrid_type_y

  double precision :: spacegrid_dxi, spacegrid_xi, spacegrid_xo
  double precision :: spacegrid_equator_eps
  double precision :: spacegrid_refine_a, spacegrid_refine_thetar
  double precision :: spacegrid_dxixi

  integer :: spacegrid_refine_n, spacegrid_radius_read
  integer :: spacegrid_theta_read

  ! >>>> Angular grid <<<<
  double precision :: rmu(-FRSIZE_MU_HALF:FRSIZE_MU_HALF, 0:FRSIZE_DRR)
  double precision :: rphi(-1:FRSIZE_PHI+2, 0:FRSIZE_DRR)
  double precision :: dmuvdr
  double precision :: rmu_i(-FRSIZE_MU_HALF-1:FRSIZE_MU_HALF+1, 0:FRSIZE_DRR)
  double precision :: rphi_i(-1:FRSIZE_PHI+3, 0:FRSIZE_DRR)

  integer :: nrmu(0:FRSIZE_DRR), nrphi(0:FRSIZE_DRR)
  integer :: iextrmu(0:FRSIZE_DRR)
  integer :: iang_warn_mureslow

  double precision :: anggrid_drr, anggrid_muerr_max
  integer :: anggrid_frsizemu, anggrid_frsizephi
  integer :: anggrid_mu_type, anggrid_mu_zero

  ! >>>> Reindex arrays <<<<
  integer :: ridx_ip(-FRSIZE_PHI-4:2*FRSIZE_PHI+4, -4:FRSIZE_Y+4)
  integer :: ridx_it(-4:FRSIZE_Y+4)
  integer :: ridx_ready


end module common_grid

