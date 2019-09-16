# This program represents a Mueller matrix system for a Dual Channel Polarimeter, which includes 3 main components:
# a Wollaston Prism, a rotatable HWP and a rotation matrix out front that compensates for parallactic rotation
# This uses the pyMuellerMat library to represent the Mueller matrices and the astropy library.
#
# Main functions:
# a) Given some Stokes Parameters, the two beams of the Wollaston Prism are computed.
# b) Given a set of measurements from the two beams of the Wollaston prism and knowledge of the HWP angle and
# parallactic angle, the corresponding on-sky polarization is retrieved.
# c) Given a set of targets, plot the tracks over time and parallactic angle from the Keck telescope.


from pyMuellerMat import common_mms as cmm
from pyMuellerMat import MuellerMat

import math
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from astroplan import Observer, FixedTarget, download_IERS_A
from astropy.time import Time
import datetime

# Initialize the telescope
keck = Observer.at_site("Keck Observatory", timezone="US/Hawaii")
fig, ax = plt.subplots()
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 22
plt.rcParams['figure.figsize'] = (20, 20)


# Function to find the two beams of the Wollaston prism based on the Stokes parameters
def wollaston(stokes):
    sys_mm = MuellerMat.SystemMuellerMatrix([cmm.WollastonPrism(), cmm.HWP()])
    sys_mm.master_property_dict['WollastonPrism']['beam'] = 'o'
    pos = sys_mm.evaluate() @ stokes
    sys_mm.master_property_dict['WollastonPrism']['beam'] = 'e'
    neg = sys_mm.evaluate() @ stokes

    return [pos, neg]


# Function that plots the difference of two beams of a Wollaston prism with a half-wave plate over the half-wave
# plate angle
def plot_wollaston(stokes):
    data = np.empty(shape=[0, 2])
    sys_mm = MuellerMat.SystemMuellerMatrix([cmm.WollastonPrism(), cmm.HWP()])

    # Find data points from 0 to 2 * pi
    for angle in np.arange(0, 2 * math.pi, 0.001):
        sys_mm.master_property_dict['HalfwaveRetarder']['theta'] = math.degrees(angle)
        sys_mm.master_property_dict['WollastonPrism']['beam'] = 'o'
        I1 = sys_mm.evaluate() @ stokes
        sys_mm.master_property_dict['WollastonPrism']['beam'] = 'e'
        I2 = sys_mm.evaluate() @ stokes
        data = np.append(data, [[math.degrees(angle), (I1[0] - I2[0])]], axis=0)

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
# over the parallactic angle and time
def track_plot(targets):
    # Initialize the start time, the targets, and the initial stokes vector
    time = Time("2015-09-13")
    step = np.arange(0, 1, 1 / 86400)
    stokes = [[0], [1], [0], [0]]
    hwp_angles = [0, 22.5]

    derotator = cmm.DiattenuatorRetarder()
    m3 = cmm.DiattenuatorRetarder()
    sys_mm = MuellerMat.SystemMuellerMatrix([cmm.WollastonPrism(), derotator, cmm.HWP(), m3, cmm.Rotator()])

    # Put in M3 - use astropy for altitude - diattenuating rotator - as a perfect mirror with an angle
    # "perfect" - no retardance and no diattenuation
    # Derotator - diattenuating retarder at a given parallactic angle
    # Check diattenuating retarder form with goldstein and witzel
    # Can calculate coefficients from material parameters - Fresnel reflection - index of refraction
    # Fresnel coefficients - how to get the r values and possibly retardance

    # use hour angle and dec to find the parallactic angle
    # find the altitude given an hour angle and a target

    for hwp in hwp_angles:
        angle_plot = []
        time_plot = []
        sys_mm.master_property_dict['HalfwaveRetarder']['theta'] = hwp

        for j in range(len(targets)):
            wollaston_data = []
            target = FixedTarget.from_name(targets[j])

            # Calculate the parallactic angles and the altitudes
            angles = np.degrees((keck.parallactic_angle(time + step, target)).to_value())
            altitudes = (keck.altaz(time + step, target)).alt.to_value()

            # Calculate the Wollaston beams and parallactic angle as time passes
            for pa, alt in zip(angles, altitudes):
                sys_mm.master_property_dict['Rotator']['pa'] = pa
                m3.properties['theta'] = alt
                sys_mm.master_property_dict['WollastonPrism']['beam'] = 'o'
                I1 = sys_mm.evaluate() @ stokes
                sys_mm.master_property_dict['WollastonPrism']['beam'] = 'e'
                I2 = sys_mm.evaluate() @ stokes
                wollaston_data.append(np.asscalar(I1[0] - I2[0]))

            angle_plot.append(np.array([angles, wollaston_data]).T)
            time_plot.append(np.array([((time + step).to_datetime()), wollaston_data]).T)

        # Plot the angle data points
        for k in range(len(targets)):
            x, y = angle_plot[k].T
            plt.scatter(x, y, s=1, label=targets[k])
        plt.title('Difference between Wollaston prism beams over parallactic angle with HWP at %.1f degrees' % hwp)
        plt.ylabel('Difference between $\mathdefault{I^+}$ and $\mathdefault{I^-}$')
        plt.xlabel('Parallactic angle (deg)')
        plt.legend(loc="upper left")
        plt.show()

        # Plot the time data points
        for k in range(len(targets)):
            x, y = time_plot[k].T
            plt.scatter(x, y, s=1, label=targets[k])
        plt.title('Difference between Wollaston prism beams over time with HWP at %.1f degrees' % hwp)
        plt.ylabel('Difference between $\mathdefault{I^+}$ and $\mathdefault{I^-}$')
        plt.xlabel('Time (hour of day)')
        plt.legend(loc="upper left")
        ax = plt.gca()
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_xlim(datetime.date(2015, 9, 13), datetime.date(2015, 9, 14))
        plt.show()


# Main function that prompts user for input
def main():
    download_IERS_A()

    print("This program represents a Mueller matrix system for a dual channel polarimeter using the pyMuellerMat"
          "library.")
    find = ""

    while find != "d":
        # Prompt
        find = input("\nWhat would you like to do?\na) compute the two beams of the Wollaston prism from Stokes "
                     "parameters\nb) find the corresponding on-sky polarization with a set of measurements from "
                     "the two beams of the Wollaston prism and HWP/parallactic angles data\nc) plot the tracks of "
                     "a set of targets over time and parallactic angle from the Keck telescope\n"
                     "d) quit the program\n(a/b/c/d): ").lower()

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
                hwp = math.radians(float(input("Enter the HWP angle (deg): ")))
                sky = math.radians(float(input("Enter the parallactic angle (deg): ")))
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
                target = input("\nEnter the name of the target to track: ")
                targets.append(target)
                cont = input("Add another target? (y/n): ").lower()
            print("\nTracking", ', '.join(targets), "...\n")
            track_plot(targets)
            input("\n---------------------------------------------------------------------------------------------")

    print("\nProgram ended.\n")


if __name__ == "__main__":
    main()
