# This program represents a Mueller matrix system for a Dual Channel Polarimeter, which includes 3 main components:
# a Wollaston Prism, a rotatable HWP and a rotation matrix out front that compensates for parallactic rotation
# This uses numpy arrays to represent Mueller matrices.
#
# Main functions:
# a) Given some Stokes Parameters, the two beams of the Wollaston Prism are computed.
# b) Given a set of measurements from the two beams of the Wollaston prism and knowledge of the HWP angle and
# parallactic angle, the corresponding on-sky polarization is retrieved.

import math
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt


# Represents a Mueller matrix
def m(theta, phi):
    return np.array([[1, -math.cos(2 * theta), 0, 0], [-math.cos(2 * theta), 1, 0, 0],
                     [0, 0, math.sin(2 * theta) * math.cos(phi), math.sin(2 * theta) * math.sin(phi)],
                     [0, 0, -math.sin(2 * theta) * math.sin(theta), math.sin(2 * theta) * math.cos(phi)]])


# Represents a rotation matrix
def t(angle):
    return np.array([[1, 0, 0, 0], [0, math.cos(4 * angle), math.sin(4 * angle), 0],
                     [0, -math.sin(4 * angle), math.cos(4 * angle), 0], [0, 0, 0, 1]])


# Wollaston prism Mueller matrices, representing opposite polarization states
m_woll_pos = m((math.pi / 2), 0)
m_woll_neg = m(math.pi, 0)


# Function to find the two beams of the Wollaston prism based on the Stokes parameters
def wollaston(stokes):
    pos = 0.5 * (m_woll_pos @ stokes)
    neg = 0.5 * (m_woll_neg @ stokes)

    return [pos, neg]


# Function that plots the difference of two beams of a Wollaston prism over the half-wave plate angle
def plot_wollaston(stokes):
    data = np.empty(shape=[0, 2])

    # Find data points from 0 to 2 * pi
    for angle in np.arange(0, 2 * math.pi, 0.001):
        I1 = 0.5 * (m_woll_pos @ t(angle) @ stokes)
        I2 = 0.5 * (m_woll_neg @ t(angle) @ stokes)
        data = np.append(data, [[math.degrees(angle), (I1[0] - I2[0])]], axis=0)

    # Plot the data points
    plt.rcParams['font.family'] = 'Times New Roman'
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

    # Calculate the Mueller matrices
    for j in range(len(values)):
        hwp = values[j][2]
        sky = values[j][3]
        row1 = 0.5 * (m_woll_pos @ t(hwp) @ t(sky))
        row2 = 0.5 * (m_woll_neg @ t(hwp) @ t(sky))
        i = np.append(i, [[values[j][0]]], axis=0)
        m_system = np.append(m_system, [row1[0]], axis=0)
        i = np.append(i, [[values[j][1]]], axis=0)
        m_system = np.append(m_system, [row2[0]], axis=0)

    # Return a least-squares solution
    return inv(np.transpose(m_system) @ m_system) @ np.transpose(m_system) @ i


# Main function that prompts user for input
def main():
    print("This program represents a Mueller matrix system for a dual channel polarimeter.")
    find = ""

    while find != "c":
        # Prompt
        find = input("\nWhat would you like to do?\na) compute the two beams of the Wollaston prism from Stokes "
                     "parameters\nb) find the corresponding on-sky polarization with a set of measurements from "
                     "the two beams of the Wollaston prism and HWP/parallactic angles data\n"
                     "c) quit the program\n(a/b/c): ").lower()

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

    print("\nProgram ended.\n")


if __name__ == "__main__":
    main()
