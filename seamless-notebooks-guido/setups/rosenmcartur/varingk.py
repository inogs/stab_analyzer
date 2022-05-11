import numpy
import scipy.integrate
import pyfabm
from mpl_toolkits.mplot3d import Axes3D
# Axes3D import has side effects, it enables using projection='3d' in add_subplot
import matplotlib.pyplot as plt
# Regular imports
from copy import deepcopy


# Yaml commentary
from ruamel.yaml.comments import \
    CommentedMap as OrderedDict, \
    CommentedSeq as OrderedList

# For manual creation of tokens
from ruamel.yaml.tokens import CommentToken
from ruamel.yaml.error import CommentMark

import ruamel.yaml as yaml

# Number of spaces for an indent 
INDENTATION = 2 
# Used to reset comment objects
tsRESET_COMMENT_LIST = [None, [], None, None]


K = numpy.arange(0.1,4.0,0.08)

fig0 = plt.figure()
h = fig0.add_subplot(1,1,1, projection='3d')
Npoints = 300

for k in K :
    k = float(k)
    fabm_list = {
        "check_conservation": False,
        "require_initialization": True,
        "instances": {
            "P1": {
                "long_name": "prey",
                "model": "rosenmcartur/prey",
                "parameters": {
                    "r": 0.5,
                    "K": k
                },
                "initialization": {
                    "DW": 3.0},
                "coupling": None 
                },
            "Z1": {
                "long_name": "predator",
                "model": "rosenmcartur/predator",
                "parameters": {
                    "r": 0.5,
                    "m": 0.15,
                    "e": 0.6,
                    "g": 0.4,
                    "H": 0.6,
                    "K": k
                },
                "initialization": {
                    "DWz": 3.0},
                "coupling": {
                    "DW": "P1/DW"
                }
                }
            }
    }
    
    with open(r'fabm0.yaml', 'w') as outfile:
        documents = yaml.round_trip_dump(fabm_list, outfile)
    # Create model (loads fabm.yaml)
    model = pyfabm.Model('fabm0.yaml')

    # Configure the environment
    # Note: the set of environmental dependencies depends on the loaded biogeochemical model.
    #model.dependencies['surface_downwelling_photosynthetic_radiative_flux'].value = 50.
#model.dependencies['downwelling_photosynthetic_radiative_flux'].value = 25.

    # Verify the model is ready to be used
    model.cell_thickness=1.

    assert model.checkReady(), 'One or more model dependencies have not been fulfilled.'

    # Time derivative
    def dy(y, t0):
        model.state[:] = y
        return model.getRates()

    # Time-integrate over 200 days (note: FABM's internal time unit is seconds!)
    t = numpy.linspace(0, 2000., 1000)
    y = scipy.integrate.odeint(dy, model.state, t*86400)
    vk = []
    for i in range(Npoints) :
        vk.append(k)
    h.scatter(vk,y[-Npoints:,0],y[-Npoints:,1], c='black')

plt.xlabel('K[mg DW L^-1]')
plt.ylabel('P[mg DW L^-1]')
h.set_zlabel('Z[mg DW L^-1]')
plt.show()
# Plot results
#import pylab
#pylab.plot(t, y)
#pylab.legend([variable.path for variable in model.state_variables])
#pylab.show()
