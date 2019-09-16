# Tutorial
This tutorial will go through the final_offline_mm.py program. The other two programs are intuitive if you know this one.

This program will always loop back to the main menu prompt, unless specified to quit (`e`):

```
This offline program represents a Mueller matrix system for a dual channel polarimeter using the pyMuellerMat library.

What would you like to do?
```
[a) compute the two beams of a Wollaston prism from Stokes parameters](#computing-the-two-beams-of-a-wollaston-prism-from-stokes-parameters)

[b) find the corresponding on-sky polarization with a set of measurements from the two beams of the Wollaston prism and HWP/parallactic angles data](#finding-the-corresponding-on-sky-polarization-with-the-wollaston-prism-beams)

[c) plot the tracks of a set of targets over time and parallactic angle from the Keck telescope](#plotting-the-tracks-of-a-set-of-targets-over-time-and-parallactic-angle-from-the-keck-telescope)

[d) test the model fit to find the diattenuation and retardance values from randomized Stokes values](#testing-the-model-fit-to-find-the-diattenuation-and-retardance-values-from-randomized-stokes-values)

e) quit the program
```
(a/b/c/d/e):
```

## Computing the two beams of a Wollaston prism from Stokes parameters
Entering `a` gives us
```
Enter the Stokes parameters, separated by a space:
```

As an example, let's say the Stokes I = 0, Q = 1, U = 0, and V = 0. So we type `0 1 0 0`, giving us
```
The two beams from the Wollaston prism are  0.5  and  -0.5

Would you like to see a plot of the intensities? (y/n): 
```
`y` gives us a graph of the difference in intensities of the Wollaston prism beams over the HWP angle.
![woll_prism_plot](https://user-images.githubusercontent.com/12475630/64978821-398aaa80-d86b-11e9-9b69-af73fb550cf9.png)
`n` would skip this step.

## Finding the corresponding on-sky polarization with the Wollaston prism beams
Entering `b` at the main menu prompts us for some information, each separated by a newline. You will need **at least 2** pairs of data for this function to work, otherwise it invokes a `Singular Matrix` error.
```
Enter the first I parameter in the pair (positive Wollaston):
Enter the second I parameter in the pair (negative Wollaston): 
Enter the HWP angle (deg):
Enter the parallactic angle (deg):

Do you have another pair of data to add? (y/n): 
```
It then calculates and outputs the corresponding on-sky polarization Stokes vector.

**Note: This least-squares method is explained in Wictorowicz (2012).**

## Plotting the tracks of a set of targets over time and parallactic angle from the Keck telescope
Entering `c` at the prompt gives us
```
Enter the name of the target to track:
Enter the hour angle:
Enter the declination:
Add another target? (y/n):
```
**Note: If you use the astropy_mm.py program, it just asks for the target's name, since it looks up the target in the SIMBAD astronomical database.**

You can change what telescope to use by editing the initialization at the top of the program:
```python
# Initialize the telescope's latitude
keck = np.radians(19.8260)
```

We then get two graphs, one over parallactic angle and one over time, when the HWP is 0 degrees:
![parallactic_angle_hwp0](https://user-images.githubusercontent.com/12475630/64979168-05fc5000-d86c-11e9-949f-2effaa888d1e.png)
![time_hwp0](https://user-images.githubusercontent.com/12475630/64979173-085eaa00-d86c-11e9-9843-ef946505adb3.png)
And two similar graphs for when the HWP is 22.5 degrees:
![parallactic_angle_hwp225](https://user-images.githubusercontent.com/12475630/64979174-085eaa00-d86c-11e9-99cc-d11f8317fae4.png)
![time_hwp225](https://user-images.githubusercontent.com/12475630/64979175-08f74080-d86c-11e9-954d-6d212190a918.png)

## Testing the model fit to find the diattenuation and retardance values from randomized Stokes values
Entering `d` at the main menu finds the diattenuation and retardance values by first taking known polarization standards to find initial Stokes vectors, using the Mueller matrix system to find the final Stokes vectors, and then using a given signal-to-noise ratio to generate some noisy Stokes values. The input and noisy Stokes vectors are then used to estimate the diattenuation and retardance values of the original derotator and third mirror.

First the program prompts for what signal-to-noise ratio to use.
```
Enter the signal-to-noise ratio:
```
The higher the signal-to-noise ratio given, the more accurate the estimations will be.

Sample diattenuation and retardance values were used as the "original" values. These were taken from van Holstein (2016). These values can be changed in the corresponding function:
```python
def fit_model(noise):
    # Initialize the original diattenuation and retardance values
    derotator_d_i, derotator_r_i = 0.9662, 186.6
    m3_d_i, m3_r_i = 0.9761, 186.6
```
Thus we have
```
Original derotator diattenuation: 0.9662 
Original derotator retardance: 186.6 
Original mirror 3 diattenuation: 0.9761 
Original mirror 3 retardance: 186.6
```
Then, the program uses three polarization standards to find the initial Stokes values. **At least two, preferrably three or more,** polarization standards should be used. These can be changed by editing the array of arrays at the top of the program. The ones currently in use are from the [United Kingdom Infrared Telescope](http://www.ukirt.hawaii.edu/instruments/irpol/irpol_stds.html) site.
```python
# Initialize the standards and convert to degrees
# HDE 279652: RA = 04 14 50.2, dec = +37 35 54, P = 0.61
# HDE 279658: RA = 04 13 47.3, dec = +37 09 32, P = 1.42
# HDE 283637: RA = 04 22 53.3, dec = +27 30 18, P = 1.57
ra = np.array([[4, 14, 50.2], [4, 13, 47.3], [4, 22, 53.3]])
dec = np.array([[37, 35, 54], [37, 9, 32], [27, 30, 18]])
p = np.array([0.0061, 0.0142, 0.0157])
```
The program then finds the final Stokes values and generates noisy Stokes vectors by drawing random samples from a normal distribution based on the signal-to-noise ratio. Finally, the noisy Stokes and the input Stokes values are used to estimate the original diattenuation and retardance values for the derotator and third mirror.
```python
curve_fit(system, np.repeat(stokes_i[str(theta)], 2, axis=0), noisy_f, bounds=(0, (1, 360, 1, 360)))
```
**Note: The input Stokes vector array is repeated once to accomodate for the fact that each of the Stokes vectors generated from the Wollaston prism beam came from the same input Stokes vector. As such, the input needs to be repeated once to have the same number of inputs as outputs.**

The bounds represent the bounds for the derotator diattenuation, the derotator retardance, the third mirror diattenuation, and the third mirror retardance respectively.

This entire process is then repeated for a total of four times, with the HWP angle changing from 0 degrees, 22.5 degrees, 45 degrees, and 60 degrees. It then reports the estimated values with their percent errors:
```
Estimated derotator diattenuation: 0.9702 with 0.411 % error 
Estimated derotator retardance: 183.6684 with 1.596 % error 
Estimated mirror 3 diattenuation: 0.9875 with 1.169 % error 
Estimated mirror 3 retardance: 184.0129 with 1.386 % error
```
This indicates a relatively low error rate.
