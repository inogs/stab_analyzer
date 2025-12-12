#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - phytoplankton component

module rma_predator
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_rma_predator
      ! Variable identifiers
      type (type_state_variable_id)      :: id_predator
      type (type_state_variable_id)      :: id_prey
      type (type_diagnostic_variable_id) :: id_rosen
      type (type_diagnostic_variable_id) :: id_grazing

      ! Model parameters
      real(rk) :: m,e,g,H!,r,K
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_rma_predator), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_p

      self%dt = 86400._rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
!      call self%get_parameter(self%r,    'r',    'd-1',        'growth rate ',            default=0.5_rk)
      call self%get_parameter(self%m,    'm',    'd-1',        'mortality',               default=0.15_rk)
      call self%get_parameter(self%e,    'e',    '',           'efficency',               default=0.6_rk)
      call self%get_parameter(self%g,    'g',    'd-1',        'grazing rate',            default=0.4_rk)
      call self%get_parameter(self%H,    'H',    'mg DW m-1',  'half saturation',         default=0.6_rk)
!      call self%get_parameter(self%K,    'K',    'mg DW m-1',  'carrying capacity',       default=2.6_rk)

      ! Register state variables
      call self%register_state_variable(self%id_predator, 'DWz', 'mg DW m-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
    !  call self%register_state_variable(self%id_prey, 'DW', 'mg DW m-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
      call self%register_state_dependency(self%id_prey, 'DW', 'mg DW m-1', 'dryweight')

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_rosen,  'growthz', 'mg DW m-1 d-1', 'growth rate')
      call self%register_diagnostic_variable(self%id_grazing,  'grazing', 'mg DW m-1 d-1', 'grazing')

      ! Register environmental dependencies

      ! Contribute to light attentuation
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_rma_predator), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: Z,P
      real(rk)            :: rosen,grazing
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_predator,Z)         ! zooplankton
         _GET_(self%id_prey,P)             ! phytoplankton

         ! Retrieve current environmental conditions.

         ! Light acclimation formulation based on surface light intensity.

         ! Loss rate of phytoplankton to detritus depends on local light intensity.

         ! Define some intermediate quantities that will be reused multiple times.
         rosen = self%e*self%g*P*Z/(P+self%H)-self%m*Z
         grazing = -self%g*P*Z/(P+self%H)
         _SET_ODE_(self%id_predator, rosen)
         _SET_ODE_(self%id_prey, grazing)
         _SET_DIAGNOSTIC_(self%id_rosen, rosen)
         _SET_DIAGNOSTIC_(self%id_grazing, grazing)
         
         ! Set temporal derivatives

         ! If an externally maintained ...

         ! Export diagnostic variables

      ! Leave spatial loops (if any)
   _LOOP_END_
   end subroutine do



end module rma_predator

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
