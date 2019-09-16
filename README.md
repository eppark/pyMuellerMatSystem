# pyMuellerMatSystem
Three Python programs representing a Mueller Matrix system for a dual-channel polarimeter.
These are its possible functions:
1) Given some Stokes Parameters, the two beams of a Wollaston Prism are computed.
2) Given a set of measurements from the two beams of the Wollaston prism and knowledge of the HWP angle and parallactic angle, the corresponding on-sky polarization is retrieved.
3) Given a set of targets, plot the tracks over hour angle and parallactic angle from the Keck telescope.
4) Test the model fit to find the diattenuation and retardance values from normal randomized Stokes values.

## manual_mm.py
This uses numpy arrays to represent Mueller matrices.
Its three main components are a Wollaston Prism, a rotatable HWP and a rotation matrix out front that compensates for parallactic rotation.

It can perform functions **1 and 2**.

## astropy_mm.py
This uses the pyMuellerMat library for its Mueller matrices and the Astropy package for celestial target calculations from the Keck telescope.
Its components are a Wollaston Prism, a diattenuating retarder representing the derotator, a rotatable HWP, a diattenuating retarder representing the third mirror, and a rotation matrix that compensates for parallactic rotation.

It can perform functions **1-3**.

## final_offline_mm.py
This uses the pyMuellerMat library for its Mueller matrices, but _not_ the Astropy package.
Its components are a Wollaston Prism, a diattenuating retarder representing the derotator, a rotatable HWP, a diattenuating retarder representing the third mirror, and a rotation matrix that compensates for parallactic rotation.

It can perform functions **1-4**.
