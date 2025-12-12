#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - phytoplankton component

module rma_prey 
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_rma_prey
      ! Variable identifiers
      type (type_state_variable_id)      :: id_prey
      type (type_diagnostic_variable_id) :: id_logistic

      ! Model parameters
      real(rk) :: r,K
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_rma_prey), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_p

      self%dt = 86400._rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%r,    'r',    'd-1',        'growth rate ',            default=0.5_rk)
      call self%get_parameter(self%K,    'K',    'mg DW m-1',  'carrying capacity',       default=2.6_rk)

      ! Register state variables
      call self%register_state_variable(self%id_prey, 'DW', 'mg DW m-1', 'dryweight', 0.0_rk, minimum=0.0_rk)


      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_logistic,  'growth', 'mg DW m-1 d-1', 'growth rate')

      ! Register environmental dependencies

      ! Contribute to light attentuation
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_rma_prey), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: P
      real(rk)            :: logistic
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_prey,P)         ! phytoplankton

         ! Retrieve current environmental conditions.

         ! Light acclimation formulation based on surface light intensity.

         ! Loss rate of phytoplankton to detritus depends on local light intensity.

         ! Define some intermediate quantities that will be reused multiple times.
         logistic = self%r*P*(1.0-P/self%K)
         _SET_ODE_(self%id_prey, logistic)
         _SET_DIAGNOSTIC_(self%id_logistic, logistic)

         ! Set temporal derivatives

         ! If an externally maintained ...

         ! Export diagnostic variables

      ! Leave spatial loops (if any)
   _LOOP_END_
   end subroutine do



end module rma_prey

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
