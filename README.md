# stab_analyzer
This repository contains the code which has been used to study the temporal stability of ecosystems in (Occhipinti et al. 2023, Occhipinti et al. 2024):
- **lyapunov.py** A python class which computes the maximum Lyapunov exponent of a timeseries following Wolf et al. (1985).
- **lyapunovV.py** A faster version of the class lyapunov.py, thanks to a vectorization of procedures.
- **seamless-notebooks-guido** The code used in the studies (Occhipinti et al. 2023, Occhipinti et al. 2024)
  - setups: the setup files and the python scripts to launch models
  - parsac: contains scripts for sensitivity analysis
  - extern: contain a model to reproduce the results of Fussmann et al. (2002) and a modification to Parsac sensitivity analysis tool.
- **ExternalForcings** Contains the script for a paper under revision 



## References
Fussmann, G. F., & Heber, G. (2002). Food web complexity and chaotic population dynamics. Ecology Letters, 5(3), 394-401.

Occhipinti, G., Solidoro, C., Grimaudo, R., Valenti, D., & Lazzari, P. (2023). Marine ecosystem models of realistic complexity rarely exhibits significant endogenous non-stationary dynamics. Chaos, Solitons & Fractals, 175, 113961.

Occhipinti, G., Piani, S., & Lazzari, P. (2024). Stochastic effects on plankton dynamics: Insights from a realistic 0-dimensional marine biogeochemical model. Ecological Informatics, 83, 102778.

Wolf, A., Swift, J. B., Swinney, H. L., & Vastano, J. A. (1985). Determining Lyapunov exponents from a time series. Physica D: nonlinear phenomena, 16(3), 285-317
