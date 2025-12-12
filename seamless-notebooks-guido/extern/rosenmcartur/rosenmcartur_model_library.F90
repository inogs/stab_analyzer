module rosenmcartur_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use rma_prey 
   use rma_predator

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: rosenmcartur_model_factory

contains

   subroutine create(self, name, model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('prey'); allocate(type_rma_prey::model)
         case ('predator'); allocate(type_rma_predator::model)
      end select
   end subroutine create

end module rosenmcartur_model_library
