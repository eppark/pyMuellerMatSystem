# This program represents a Mueller matrix system for a Dual Channel Polarimeter, which includes 3 main components:
# a Wollaston Prism, a rotatable HWP and a rotation matrix out front that compensates for parallactic rotation
# This uses the pyMuellerMat library to represent the Mueller matrices.
#
# Main functions:
# a) Given some Stokes Parameters, the two beams of a Wollaston Prism are computed.
# b) Given a set of measurements from the two beams of the Wollaston prism and knowledge of the HWP angle and
# parallactic angle, the corresponding on-sky polarization is retrieved.
# c) Given a set of targets, plot the tracks over hour angle and parallactic angle from the Keck telescope.
# d) Test the model fit to find the diattenuation and retardance values from randomized Stokes values.


from pyMuellerMat import common_mms as cmm
from pyMuellerMat import MuellerMat

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Initialize the telescope's latitude
keck = np.radians(19.8260)

# Initialize the finalized system Mueller Matrix
derotator = cmm.DiattenuatorRetarder()
derotator.name = 'derotator'
m3 = cmm.DiattenuatorRetarder()
m3.name = 'm3'
master_sys_mm = MuellerMat.SystemMuellerMatrix([cmm.WollastonPrism(), derotator, cmm.HWP(), m3, cmm.Rotator()])

# Initialize the standards and convert to degrees
# http://www.ukirt.hawaii.edu/instruments/irpol/irpol_stds.html
# HDE 279652: RA 04 14 50.2, dec +37 35 54, P = 0.61
# HDE 279658: RA 04 13 47.3, dec +37 09 32, P = 1.42
# HDE 283637: RA 04 22 53.3, dec +27 30 18, P = 1.57
ra = np.array([[4, 14, 50.2], [4, 13, 47.3], [4, 22, 53.3]])
dec = np.array([[37, 35, 54], [37, 9, 32], [27, 30, 18]])
rad = np.radians(15 * (ra[:, 0] + ra[:, 1] / 60 + ra[:, 2] / 3600))
decd = np.radians(dec[:, 0] + dec[:, 1] / 60 + dec[:, 2] / 3600)
p = np.array([0.0061, 0.0142, 0.0157])

# matplotlib formatting
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 22
plt.rcParams['figure.figsize'] = (20, 20)


# Function to find the two beams of a Wollaston prism based on the Stokes parameters
# Input:
#       stokes: a Stokes vector (an array of 4 single-item arrays), ie. [[I], [Q], [U], [V]]
# Output:
#       [pos, neg]: an array with the positive Wollaston beam and the negative Wollaston beam
def wollaston(stokes):
    sys_mm = MuellerMat.SystemMuellerMatrix([cmm.WollastonPrism(), cmm.HWP()])
    sys_mm.master_property_dict['WollastonPrism']['beam'] = 'o'
    pos = sys_mm.evaluate() @ stokes
    sys_mm.master_property_dict['WollastonPrism']['beam'] = 'e'
    neg = sys_mm.evaluate() @ stokes

    return [pos, neg]


# Function that plots the difference of two beams of a Wollaston prism with a half-wave plate over the half-wave
# plate angle
# Input:
#       stokes: a Stokes vector (an array of 4 single-item arrays)
def plot_wollaston(stokes):
    data = np.empty(shape=[0, 2])
    sys_mm = MuellerMat.SystemMuellerMatrix([cmm.WollastonPrism(), cmm.HWP()])

    # Find data points from 0 to 2 * pi
    for angle in np.arange(0, 2 * np.pi, 0.001):
        sys_mm.master_property_dict['HalfwaveRetarder']['theta'] = np.degrees(angle)
        sys_mm.master_property_dict['WollastonPrism']['beam'] = 'o'
        I1 = sys_mm.evaluate() @ stokes
        sys_mm.master_property_dict['WollastonPrism']['beam'] = 'e'
        I2 = sys_mm.evaluate() @ stokes
        data = np.append(data, [[np.degrees(angle), (I1[0] - I2[0])]], axis=0)

    # Plot the data points
    plt.scatter(*data.T, s=1)
    plt.title('Difference between Wollaston prism beams over HWP angle')
    plt.ylabel('Difference between $\mathdefault{I^+}$ and $\mathdefault{I^-}$')
    plt.xlabel('HWP angle (deg)')
    ax = plt.gca()
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    plt.show()


# Function to find the corresponding on-sky polarization based on data of the intensities of the two beams
# of the Wollaston prism, the HWP angle, and the parallactic angle
# Input:
#       values: an array of arrays with the positive Wollaston beam, the negative Wollaston beam, the HWP angle,
#               and the parallactic angle, ie. [[pos Woll 1, neg Woll 1, HWP 1, parallactic 1], [pos Woll 2,
#               neg Woll 2, HWP 2, parallactic 2], ... [pos Woll n, neg Woll n, HWP n, parallactic n]]
# Output:
#       stokes: a Stokes vector of the calculated on-sky polarization (an array of 4 single-item arrays), ie.
#               [[I], [Q], [U], [V]]
def on_sky(values):
    i = np.empty(shape=[0, 1])
    m_system = np.empty(shape=[0, 4])
    sys_mm = MuellerMat.SystemMuellerMatrix([cmm.WollastonPrism(), cmm.HWP(), cmm.Rotator()])

    # Calculate the Mueller matrices
    for j in range(len(values)):
        sys_mm.master_property_dict['HalfwaveRetarder']['theta'] = values[j][2]
        sys_mm.master_property_dict['Rotator']['pa'] = values[j][3]

        sys_mm.master_property_dict['WollastonPrism']['beam'] = 'o'
        row1 = sys_mm.evaluate()
        sys_mm.master_property_dict['WollastonPrism']['beam'] = 'e'
        row2 = sys_mm.evaluate()

        i = np.append(i, [[values[j][0]]], axis=0)
        m_system = np.append(m_system, [row1[0]], axis=0)
        i = np.append(i, [[values[j][1]]], axis=0)
        m_system = np.append(m_system, [row2[0]], axis=0)

    # Return a least-squares solution
    return inv(np.transpose(m_system) @ m_system) @ np.transpose(m_system) @ i


# Function that plots the difference of two beams of a Wollaston prism with a half-wave plate of fixed targets
# over the parallactic angle and hour angle
# Input:
#       targets: an array of arrays with the hour angles and the declinations of the targets, ie. [[ha 1, dec 1],
#               [ha 2, dec2], ... [ha n, dec n]]
def track_plot(targets):
    # Initialize the initial stokes vector
    stokes = [[0], [1], [0], [0]]
    hwp_angles = [0, 22.5]

    for hwp in hwp_angles:
        angle_plot = []
        time_plot = []
        master_sys_mm.master_property_dict['HalfwaveRetarder']['theta'] = hwp

        for j in range(len(targets)):
            wollaston_data = []
            start_ha = targets[j][0]
            dec = np.radians(targets[j][1])
            ha = np.radians(np.arange(start=start_ha, stop=start_ha + 45, step=0.001, dtype='float'))

            # Calculate the parallactic angles and the altitudes
            angles = np.degrees(np.arctan(np.sin(ha) / (np.cos(dec) * np.tan(keck) - np.sin(dec) * np.cos(ha))))
            altitudes = np.degrees(np.arcsin(np.sin(keck) * np.sin(dec) + np.cos(keck) * np.cos(dec) * np.cos(ha)))

            # Calculate the Wollaston beams and parallactic angle as time passes
            for pa, alt in zip(angles, altitudes):
                master_sys_mm.master_property_dict['Rotator']['pa'] = pa
                m3.properties['theta'] = alt
                master_sys_mm.master_property_dict['WollastonPrism']['beam'] = 'o'
                I1 = master_sys_mm.evaluate() @ stokes
                master_sys_mm.master_property_dict['WollastonPrism']['beam'] = 'e'
                I2 = master_sys_mm.evaluate() @ stokes
                wollaston_data.append((I1[0] - I2[0]).item())

            angle_plot.append(np.array([angles, wollaston_data]).T)
            time_plot.append(np.array([np.degrees(ha), wollaston_data]).T)

        # Plot the angle data points
        for k in range(len(targets)):
            x, y = angle_plot[k].T
            plt.scatter(x, y, s=1, label=targets[k][2])
        plt.title('Difference between Wollaston prism beams over parallactic angle with HWP at %.1f degrees' % hwp)
        plt.ylabel('Difference between $\mathdefault{I^+}$ and $\mathdefault{I^-}$')
        plt.xlabel('Parallactic angle (deg)')
        plt.legend(loc="upper left")
        plt.show()

        # Plot the time data points
        for k in range(len(targets)):
            x, y = time_plot[k].T
            plt.scatter(x, y, s=1, label=targets[k][2])
        plt.title('Difference between Wollaston prism beams over hour angle with HWP at %.1f degrees' % hwp)
        plt.ylabel('Difference between $\mathdefault{I^+}$ and $\mathdefault{I^-}$')
        plt.xlabel('Coordinate (hour angle)')
        plt.legend(loc="upper left")
        plt.show()


# Function representing a system with two diattenuating retarders that returns the output Stokes value given input
# on-sky values
# Input:
#       x: on-sky Stokes values (an array of arrays of 4 floats), ie. [[I 1, Q 1, U 1, V 1],
#               [I 2, Q 2, U 2, V 2], ... [I n, Q n, U n, V n]]
#       dd: a float representing the derotator diattenuation
#       dr: a float representing the derotator retardance
#       md: a float representing the mirror 3 diattenuation
#       mr: a float representing the mirror 3 retardance
# Output:
#       I: an array of floats representing the output Stokes vector intensity values, ie. [I1, I2, ... I n]
def system(x, dd, dr, md, mr):
    # Calculate the parallactic angles and the altitudes. If these values are known, comment out the lines
    # calculating them and set the variables instead.
    angles = np.degrees(np.arctan(np.sin(rad) / (np.cos(decd) * np.tan(keck) - np.sin(decd) * np.cos(rad))))
    altitudes = np.degrees(np.arcsin(np.sin(keck) * np.sin(decd) + np.cos(keck) * np.cos(decd) * np.cos(rad)))
    # angles = []
    # altitudes = []

    derotator.properties['epsilon'] = dd
    derotator.properties['phi'] = dr
    m3.properties['epsilon'] = md
    m3.properties['phi'] = mr
    angles = np.repeat(angles, 2, axis=0)
    altitudes = np.repeat(altitudes, 2, axis = 0)

    I = []
    for (n, stokes), pa, alt in zip(enumerate(x), angles, altitudes):
        master_sys_mm.master_property_dict['Rotator']['pa'] = pa
        m3.properties['theta'] = alt
        if n % 2 == 0:
            master_sys_mm.master_property_dict['WollastonPrism']['beam'] = 'o'
        else:
            master_sys_mm.master_property_dict['WollastonPrism']['beam'] = 'e'
        I.append((master_sys_mm.evaluate() @ np.reshape(stokes, (4, 1)))[0].item())
    return I


# Function that tests the model fit to find the diattenuation and retardance values from randomized Stokes values
# Input:
#       noise: a float representing the signal-to-noise ratio
def fit_model(noise):
    # Initialize the original diattenuation and retardance values
    # Van holstein pg 36
    derotator_d_i, derotator_r_i = 0.9662, 186.6
    m3_d_i, m3_r_i = 0.9761, 186.6

    print("\nOriginal derotator diattenuation:", derotator_d_i, "\nOriginal derotator retardance:", derotator_r_i,
          "\nOriginal mirror 3 diattenuation:", m3_d_i, "\nOriginal mirror 3 retardance:", m3_r_i)

    derotator.properties['theta'] = 0
    derotator.properties['epsilon'] = derotator_d_i
    derotator.properties['phi'] = derotator_r_i
    m3.properties['theta'] = 0
    m3.properties['epsilon'] = m3_d_i
    m3.properties['phi'] = m3_r_i

    # Calculate the parallactic angles and the altitudes
    angles = np.degrees(np.arctan(np.sin(rad) / (np.cos(decd) * np.tan(keck) - np.sin(decd) * np.cos(rad))))
    altitudes = np.degrees(np.arcsin(np.sin(keck) * np.sin(decd) + np.cos(keck) * np.cos(decd) * np.cos(rad)))

    # Set up variables to store the original Stokes values and the final Stokes values
    stokes_i, stokes_f = {}, {}

    estimate = []

    hwp = [0, 22.5, 45, 60]
    for theta in hwp:
        # Initialize the system
        derotator.properties['epsilon'] = derotator_d_i
        derotator.properties['phi'] = derotator_r_i
        m3.properties['epsilon'] = m3_d_i
        m3.properties['phi'] = m3_r_i

        stokes_i[str(theta)] = []
        stokes_f[str(theta)] = []
        noisy_f = []

        # Find the initial Stokes values
        q = p * np.cos(2 * theta)
        u = p * np.sin(2 * theta)
        for x, y in zip(q, u):
            stokes_i[str(theta)].append([1, x, y, 0])

        master_sys_mm.master_property_dict['HalfwaveRetarder']['theta'] = theta

        # Calculate the final Stokes values
        for pa, alt, stokes in zip(angles, altitudes, stokes_i[str(theta)]):
            master_sys_mm.master_property_dict['Rotator']['pa'] = pa
            m3.properties['theta'] = alt
            master_sys_mm.master_property_dict['WollastonPrism']['beam'] = 'o'
            I1 = master_sys_mm.evaluate() @ np.vstack(stokes)
            master_sys_mm.master_property_dict['WollastonPrism']['beam'] = 'e'
            I2 = master_sys_mm.evaluate() @ np.vstack(stokes)
            stokes_f[str(theta)].extend([I1[0].item(), I2[0].item()])

        # Generate noisy values with the given noise we have
        for stokes in stokes_f[str(theta)]:
            noisy_f.append(np.random.normal(stokes, stokes / noise, 1)[0])

        # Estimate the diattenuation and the retardance from this noisy data
        estimate.append(curve_fit(system, np.repeat(stokes_i[str(theta)], 2, axis=0), noisy_f,
                                  bounds=(0, (1, 360, 1, 360)))[0])

    # Report the estimated values found
    estimate = np.array(estimate)
    derotator_d_f = np.mean(estimate[:, 0])
    derotator_r_f = np.mean(estimate[:, 1])
    m3_d_f = np.mean(estimate[:, 2])
    m3_r_f = np.mean(estimate[:, 3])
    error = np.abs(np.array([(derotator_d_f - derotator_d_i) / derotator_d_i * 100,
                             (derotator_r_f - derotator_r_i) / derotator_r_f * 100, (m3_d_f - m3_d_i) / m3_d_i * 100,
                             (m3_r_f - m3_r_i) / m3_r_i * 100]))

    print("\nEstimated derotator diattenuation: %.4f" % derotator_d_f, "with %.3f" % error[0], "% error",
          "\nEstimated derotator retardance: %.4f" % derotator_r_f, "with %.3f" % error[1], "% error",
          "\nEstimated mirror 3 diattenuation: %.4f" % m3_d_f, "with %.3f" % error[2], "% error",
          "\nEstimated mirror 3 retardance: %.4f" % m3_r_f, "with %.3f" % error[3], "% error")


# Main function that prompts user for input
def main():
    print("This offline program represents a Mueller matrix system for a dual channel polarimeter using the "
          "pyMuellerMat library.")
    find = ""

    while find != "e":
        # Prompt
        find = input("\nWhat would you like to do?\na) compute the two beams of a Wollaston prism from Stokes "
                     "parameters\nb) find the corresponding on-sky polarization with a set of measurements from "
                     "the two beams of the Wollaston prism and HWP/parallactic angles data\nc) plot the tracks of "
                     "a set of targets over time and parallactic angle from the Keck telescope\n"
                     "d) test the model fit to find the diattenuation and retardance values from randomized Stokes "
                     "values\ne) quit the program\n(a/b/c/d/e): ").lower()

        if find == "a":
            stokes = []

            # Get the Stokes parameters input and store into a list
            while len(stokes) != 4:
                stokes = input("\nEnter the Stokes parameters, separated by a space: ")
                stokes = stokes.split()

                if len(stokes) != 4:
                    print("Enter all four parameters!")

            stokes = np.array([[float(stokes[0])], [float(stokes[1])], [float(stokes[2])], [float(stokes[3])]])
            woll = wollaston(stokes)

            print("\nThe two beams from the Wollaston prism are ", woll[0][0][0], " and ", woll[1][0][0])

            input("\n---------------------------------------------------------------------------------------------")

            # Plot the intensities if the user chooses to
            plot = input("\nWould you like to see a plot of the intensities? (y/n): ").lower()
            if plot == 'y':
                plot_wollaston(stokes)

            print("\n---------------------------------------------------------------------------------------------")

        elif find == "b":
            # Get the intensities and angle data
            values = []
            cont = "y"

            while cont == "y":
                I_1 = float(input("\nEnter the first I parameter in the pair (positive Wollaston): "))
                I_2 = float(input("Enter the second I parameter in the pair (negative Wollaston): "))
                hwp = np.radians(float(input("Enter the HWP angle (deg): ")))
                sky = np.radians(float(input("Enter the parallactic angle (deg): ")))
                values.append([I_1, I_2, hwp, sky])
                cont = input("\nDo you have another pair of data to add? (y/n): ").lower()

            print("\nThe corresponding on-sky polarization is:\n")
            print(on_sky(values))

            input("\n---------------------------------------------------------------------------------------------")

        elif find == "c":
            # Get the targets to track
            targets = []
            cont = "y"

            while cont == "y":
                name = input("\nEnter the name of the target to track: ")
                hour_angle = float(input("Enter the hour angle: "))
                declination = float(input("Enter the declination: "))
                targets.append([hour_angle, declination, name])
                cont = input("Add another target? (y/n): ").lower()
            print("\nTracking...\n")
            track_plot(targets)
            input("\n---------------------------------------------------------------------------------------------")

        elif find == "d":
            fit_model(float(input("\nEnter the signal-to-noise ratio: ")))
            input("\n---------------------------------------------------------------------------------------------")

    print("\nProgram ended.\n")


if __name__ == "__main__":
    main()
