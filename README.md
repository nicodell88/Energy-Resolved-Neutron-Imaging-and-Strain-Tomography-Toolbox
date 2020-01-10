# Bragg-Edge-Fitting

Toolbox for processing neutron time-of-flight data for strain tomography applications.

This toolbox implements three Bragg-edge fitting methods.

1. The method presented in [[1]](#1) uses a Bragg-edge attenuation model in conjunction with a 3 parameter sigmoidal function modelling the transition.

2. The method used in [[2]](#2) uses a similar 3 parameter sigmoidal function to [[1]](#1), modified with a height and offset to model the entire Bragg-edge.

3. A newly developed method, presented in [[3]](#3), which uses the Bragg-edge attenuation model, as in [[1]](#1), but instead of modelling the transition with a parametric function, we model the residual with a Gaussian Process (GP) and utilise the properties of the GP and a Monte-Carlo step to obtain an estimate of the Bragg-edge location and confidence interval.

## References
<a id="1">[1]</a>  Santisteban, J., Edwards, L., Steuwer, A., Withers, P., 2001. Time-of-flight neutron transmission diffraction. Journal of applied crystallography 34 (3), 289-297. [https://onlinelibrary.wiley.com/doi/pdf/10.1107/S0021889801003260](https://onlinelibrary.wiley.com/doi/pdf/10.1107/S0021889801003260)

<a id="2">[2]</a> 
Tremsin, A. S., Gao, Y., Dial, L. C., Grazzi, F., Shinohara, T., 2016. Investigation of microstructure in additive manufactured inconel 625 by spatially resolved neutron transmission spectroscopy. Science and Technology of advanced MaTerialS 17 (1), 324?336. [https://doi.org/10.1080/14686996.2016.1190261](https://doi.org/10.1080/14686996.2016.1190261)

<a id="3">[3]</a> 
TODO