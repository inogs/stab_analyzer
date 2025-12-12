add_library(fabm_models_rosenmcartur OBJECT
            rosenmcartur_model_library.F90
	    predator.F90
	    prey.F90
	    )


add_dependencies(fabm_models_rosenmcartur fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------
