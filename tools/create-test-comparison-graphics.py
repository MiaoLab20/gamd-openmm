#!/usr/bin/env python

import sys
import os
import matplotlib.pyplot as plt
from pprint import pprint


def main():
    if len(sys.argv) < 5:
        usage()
        sys.exit(1)
    cmd_directory = sys.argv[1]
    gamd_1_directory = sys.argv[2]
    gamd_2_directory = sys.argv[3]
    gamd_3_directory = sys.argv[4]
    phi_filepath = "/graphics-out/pmf-phi-reweight-CE2.xvg"
    psi_filepath = "/graphics-out/pmf-psi-reweight-CE2.xvg"

    graph_title = os.getcwd()
    create_graphic(cmd_directory, gamd_1_directory, gamd_2_directory, gamd_3_directory, phi_filepath, "Phi",
                   graph_title)
    create_graphic(cmd_directory, gamd_1_directory, gamd_2_directory, gamd_3_directory, psi_filepath, "Psi",
                   graph_title)


#
# Driver method for creating the graphics from the coordinates in the directory for a given file type.
#
def create_graphic(cmd_directory, gamd_1_directory, gamd_2_directory, gamd_3_directory, filepath, xlabel, graph_title):
    coordinates = gather_coordinates(filepath, cmd_directory, gamd_1_directory, gamd_2_directory, gamd_3_directory)
    graph_values = create_averages_and_errors(coordinates)
    generate_graphic_file(graph_values, xlabel, graph_title)


def usage():
    print("This program is meant to take a conventional MD run and three gamd runs to generate a graph to compare the conventional MD run to the average of the three short gamd runs.\n")
    print("Usage:")
    print("create_graphic.py cmd-directory gamd-directory gamd-directory gamd-directory")


def is_number(candidate):
    try:
        float(candidate)
    except ValueError:
        return False
    return True

#
#  Function that actually generates the image file.
#


def generate_graphic_file(coordinates, xlabel, graph_title):
    cmd_x = coordinates[0]["x"]
    cmd_y = coordinates[0]["y"]
    gamd_x = coordinates[1]
    gamd_y = coordinates[2][0]
    gamd_y_errors = [coordinates[2][1], coordinates[2][2]]
    plt.plot(cmd_x, cmd_y, "", label="Conventional MD")
    plt.axis([-180,180,0,8])

    plt.ylabel("PMF")
    plt.xlabel(xlabel)
    plt.errorbar(gamd_x, gamd_y, yerr=gamd_y_errors, capsize=2, errorevery=1, alpha=0.5, label="GaMD")
    plt.legend()
    plt.title(graph_title)

    plt.savefig("1D-" + xlabel + ".png")
    plt.close()



#
# These are the routines used to calculate the averages, errors, and put things in an easy to access
# way for doing our graphing.
#

def create_averages_and_errors(coordinates):
    cmd = coordinates[0]
    gamd_1 = coordinates[1]
    gamd_2 = coordinates[2]
    gamd_3 = coordinates[3]

    x_values = []
    y_averages = []
    y_min_errors = []
    y_max_errors = []

    if len(gamd_1) == len(gamd_2) == len(gamd_3):
        for entry in range(len(gamd_1['x'])):
            x_values.append(gamd_1["x"][entry])
            y_values = calculate_average_and_errors(gamd_1, gamd_2, gamd_3, "y", entry)

            y_averages.append(y_values[0])
            y_min_errors.append(y_values[1])
            y_max_errors.append(y_values[2])

    else:
        print("Error: uneven number of coordinates between files.")
        sys.exit(2)
    
    return [cmd, x_values, [y_averages, y_min_errors, y_max_errors]]


def get_average(gamd_1, gamd_2, gamd_3, coordinate_type, entry):
    return (gamd_1[coordinate_type][entry] + gamd_2[coordinate_type][entry] + gamd_3[coordinate_type][entry]) / 3


def get_minimum_error(gamd_1, gamd_2, gamd_3, coordinate_type, entry, average):
    return average - min(gamd_1[coordinate_type][entry], gamd_2[coordinate_type][entry], gamd_3[coordinate_type][entry])


def get_maximum_error(gamd_1, gamd_2, gamd_3, coordinate_type, entry, average):
    return max(gamd_1[coordinate_type][entry], gamd_2[coordinate_type][entry], gamd_3[coordinate_type][entry]) - average


def calculate_average_and_errors(gamd_1, gamd_2, gamd_3, coordinate_type, entry):
    average = get_average(gamd_1, gamd_2, gamd_3, coordinate_type, entry)
    minimum_error = get_minimum_error(gamd_1, gamd_2, gamd_3, coordinate_type, entry, average)
    maximum_error = get_maximum_error(gamd_1, gamd_2, gamd_3, coordinate_type, entry, average)

    return [average, minimum_error, maximum_error]


#
# These are the routines to gather the information from the different files.
#
#

def gather_coordinates(filepath, cmd_directory, gamd_1_directory, gamd_2_directory, gamd_3_directory):
    cmd = get_coordinates_from_file(cmd_directory + filepath)
    gamd_1 = get_coordinates_from_file(gamd_1_directory + filepath)
    gamd_2 = get_coordinates_from_file(gamd_2_directory + filepath)
    gamd_3 = get_coordinates_from_file(gamd_3_directory + filepath)
    return [cmd, gamd_1, gamd_2, gamd_3]


def get_coordinates_from_file(filepath):
    with open(filepath, "r") as input_file:
        lines = input_file.readlines()
    result = {}
    x_values = []
    y_values = []
    for line in lines:
        fields = line.split()
        if len(fields) > 0 and is_number(fields[0]):
            x_values.append(float(fields[0]))
            y_values.append(float(fields[1]))
    result["x"] = x_values
    result["y"] = y_values

    return result


if __name__ == '__main__':
    main()
