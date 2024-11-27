import numpy as np
import sys
sys.path.append("C:/Users/ivana/BackUp/PyKinematics/PyKinematics/src")
# print(sys.path)
from Kinematics import *
import time
import pathlib
import csv
import math
import copy
import tqdm
import os

os.makedirs('paraview', exist_ok=True)
import pandas as pd


### This f gives you the lat value from load_multiplier.txt

def check_last_value():
    global load_multiplier
    try:
        with open("load_multiplier.txt", "r") as file:
            # Read all lines in the file
            all_lines = file.readlines()  # This returns a list of all lines
            # Strip whitespace and newline characters from each line and convert to float
            all_lines = [float(line.strip()) for line in all_lines if line.strip()]
            if all_lines:  # Check if there are any lines
                load_multiplier = all_lines[-1]  # Get the last float value
            else:
                load_multiplier = None  # Set to None if the file is empty
        return load_multiplier
    except FileNotFoundError:
        return None
    except ValueError:
        # Handle cases where conversion to float fails
        return None


def read_last_displacement_values():  # !!!!! Find more general way to extract nodes values
    # Define the path to the file in the current directory
    file_path = os.path.join(os.getcwd(),
                             'C:/Users/ivana/BackUp/OpenSees-master/OpenSees-master/Win64/bin/output_from_OS/AllDispl_th.out')

    # Read the file and process the last row
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract the last row and split it into columns
    last_row = lines[-1].strip().split()

    # Assuming the file has at least three columns, and we skip the first column
    value_node_3 = float(last_row[1])  # Second column (index 1)
    value_node_6 = float(last_row[2])  # Third column (index 2)

    # Print the results
    print(f"Node 3: {value_node_3}")
    print(f"Node 6: {value_node_6}")

    # Calculate the relative displacement
    relative_displacement = value_node_6 - value_node_3

    return value_node_3, value_node_6, relative_displacement


def solve_pushover_rigid_friction_associative_2d_new(elems, contps, node_control, disp_incre, max_iteration,  current_iteration=0):
    component_index = np.nonzero(disp_incre)
    if len(component_index) > 1:
        raise Exception("Only one component of displacement can be controlled")

    contacts = []
    forces = []
    displacements = []

    start_iteration = current_iteration
    end_iteration = current_iteration + 1

    for i in tqdm.tqdm(range(start_iteration, end_iteration)):
        cal_gap_2d(contps)
        solution, solsta, contact_forces = solve_infinitefc_associative(elems, contps)
       # print("Inside pushover f, contact forces", contact_forces)

        # # Check if contact forces are all zeros
        # if np.all(np.isclose(contact_forces, 0)):
        #     print("Inside pushover f: All contact forces are zero. Stopping the analysis.")
        #     return None  # Returning None when forces are zero

        # if solution['limit_force'] <= 0:
        #     break
        # if solsta in [solsta.dual_infeas_cer, solsta.prim_infeas_cer]:
        #     print("Solution is infeasible (either primal or dual). Stopping analysis.")
        #     #return
        #     sys.exit("inside pushover f TH AN infeasible solution encountered. Exiting program.")
        #_control_disp_2d(contps, elems, node_control, disp_incre, component_index[0][0])
        forces.append(solution['limit_force'])
        #print("inside pushover solve, forces are after an", forces)
        contacts.append(solution["contact_forces"])
        # displacements.append(contps[node_control].displacement[component_index[0][0]])

        # for contp_id, contp in contps.items():
        #         # Print the ID of the contact point and its gap value
        #         print(f"Inside pushover f, Contact point {contp_id} displacement: {contp.displacement}")
        #
        # for elem_id, elem in elems.items():
        #     # Print the ID of the element and its displacement
        #     print(f"Inside pushover function, Element {elem_id} displacement: {elem.displacement}")


        displacements.append(contps[node_control].displacement[component_index[0][0]])
        _displace_model_2d(elems, contps)

        print("here")

    return forces, displacements, solsta, solution["contact_forces"]




def solve_pushover_elastic_associative_2d_new(elems, contps, node_control, unit_horizontal_load, max_iteration, Aglobal, current_iteration=0):
    print("Push over starts:")
    i = 0
    contacts = []
    forces = []
    displacements = []

    start_iteration = current_iteration
    end_iteration = current_iteration + 1

    for i in tqdm.tqdm(range(start_iteration, end_iteration)):
        #i=1
        live_load_of_this_step = 1 * unit_horizontal_load
        #live_load_of_this_step = i * unit_horizontal_load / max_iteration
       # A_matrix = cal_A_global_2d(elems, contps, sparse=True)
        cal_gap_2d(contps)

        for e in elems.values():
            if 'brick' in e.type:  # Check if 'brick' is a substring of e.type
               # e.ll = [live_load_of_this_step,0,0]
                e.dl[0] = live_load_of_this_step
                #e.dl[1] = live_load_of_this_step
                print("Live load live_load_of_this_step", live_load_of_this_step)
        solution_cf, solution_d, solution_conv = solve_elastic_infinitefc_associative_2d(
            elems, contps, Aglobal=Aglobal)
        print("Print sol [disp],", solution_d)
        print("Print sol [contf],", solution_cf)
        if solution_conv == False:
            print("Solution does not converge")
            return forces, displacements
        contacts.append(solution_cf)
        print("Contact forces are in contacts", contacts)
        forces.append(live_load_of_this_step)
        displacements.append(
            contps[node_control].displacement[0])
        print("inside pushover displ are", displacements)

        _displace_model_2d(elems, contps)

    return forces, displacements, solution_cf




def _control_disp_2d(contps, elems, node_control, disp_incre, component_index):
    print("component index is,", component_index)
    print(
        f"Converting {contps[node_control].displacement[component_index]} to {disp_incre[component_index]} for node {node_control} disp {component_index}")
    factor = disp_incre[component_index] / \
             contps[node_control].displacement[component_index]
    print("factor in displacement control", factor)
    for element in elems.values():
        element.displacement = (np.asarray(
            element.displacement) * factor).tolist()

    # ! Need to update contact point information, to record displacement
    for k, value in contps.items():
        elem_disp = np.asarray(elems[value.cand].displacement)
        # print(f"element displacement {elem_disp}")
        elem_center = elems[value.cand].center
        # print(f"element center {elem_center}")
        node_x = value.coor[0] - elem_center[0]
        node_y = value.coor[1] - elem_center[1]
        new_x = node_x * \
                math.cos(elem_disp[2]) + node_y * \
                math.sin(elem_disp[2]) + elem_disp[0] + elem_center[0]
        new_y = -node_x * \
                math.sin(elem_disp[2]) + node_y * \
                math.cos(elem_disp[2]) + elem_disp[1] + elem_center[1]
        value.displacement = [new_x - value.coor[0], new_y - value.coor[1]]
    print("Control point displacement", contps[node_control].displacement)


def _displace_model_2d(elems, contps):
    # ! NEED to update vertices information because the next step could fail
    # element center
    for key, value in elems.items():
        for pt in value.vertices:
            node_x = pt[0] - value.center[0]
            node_y = pt[1] - value.center[1]
            pt[0] = node_x * \
                    math.cos(value.displacement[2]) + node_y * \
                    math.sin(value.displacement[2]) + \
                    value.displacement[0] + value.center[0]
            pt[1] = -node_x * \
                    math.sin(value.displacement[2]) + node_y * \
                    math.cos(value.displacement[2]) + \
                    value.displacement[1] + value.center[1]

    for k, value in contps.items():
        elem_disp = np.asarray(elems[value.cand].displacement)
        # print(f"element displacement {elem_disp}")
        elem_center = elems[value.cand].center
        # print(f"element center {elem_center}")
        node_x = value.coor[0] - elem_center[0]
        node_y = value.coor[1] - elem_center[1]
        value.coor[0] = node_x * \
                        math.cos(elem_disp[2]) + node_y * \
                        math.sin(elem_disp[2]) + elem_disp[0] + elem_center[0]
        value.coor[1] = -node_x * \
                        math.sin(elem_disp[2]) + node_y * \
                        math.cos(elem_disp[2]) + elem_disp[1] + elem_center[1]
        # value.displacement = [0, 0]

    for key, value in elems.items():
        value.center[0] = value.center[0] + value.displacement[0]
        value.center[1] = value.center[1] + value.displacement[1]
        # value.displacement = [0, 0, 0]



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Functions necessary to calculate and pass reactions to the OS.
# # # FOR STATIC ANALYSIS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Function that calculates resultant moment on the central point of an interface (shorter edge of voussoir).
def calculate_moment(x1, y1, x2, y2, force_n1, force_n2):
    # Calculate the full length between two points
    full_length = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    # Lever arm is half the full length
    lever_arm = full_length / 2
    # Calculate the resultant moment
    return (force_n1 * lever_arm) - (force_n2 * lever_arm)


# # # Calculate resultant of contact forces on one interface; then on side of the arch (ground el. included);
# # # then total reactions after self weight analysis to be passed to OS

def sum_forces_to_reactions_static():
    global reactions_to_pass_static_final
    forces_filename = "reactions_of_blocks_static.csv"
    points_filename = "data/point.csv"

    # Read forces and points data
    forces = {}
    point_data = {}
    candidate_0_points = []  # List to store points belonging to candidate ID 0

    # Read the forces file
    with open(forces_filename, mode='r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            contact_id = int(row[0])
            forces[contact_id] = {
                't': float(row[1]),
                'n': float(row[2])
            }

    # Read the points file and collect points for candidate ID 0
    with open(points_filename, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            point_id = int(row['id'])
            candidate_id = int(row['candidate_id'])
            point_data[point_id] = {
                'x': float(row['x']),
                'y': float(row['y']),
                'candidate_id': candidate_id,
                'counter_point': int(row['counter_point'])
            }
            if candidate_id == 0:
                candidate_0_points.append(point_id)

    # Initialize moments dictionary and include special groups for candidate ID 0
    moments = {"0a": {'t': 0, 'n': 0, 'moment': 0}, "0b": {'t': 0, 'n': 0, 'moment': 0}}

    # Calculate moments for each contact point
    for point_id, data in point_data.items():
        if point_id in forces:
            # Find the next consecutive point_id for the same candidate_id
            next_point_id = point_id + 1 if (point_id + 1) in point_data and \
                                            point_data[point_id + 1]['candidate_id'] == data['candidate_id'] else None

            if next_point_id and next_point_id in forces:
                force_n1, force_t1 = forces[point_id]['n'], forces[point_id]['t']
                force_n2, force_t2 = forces[next_point_id]['n'], forces[next_point_id]['t']
                x1, y1 = data['x'], data['y']
                x2, y2 = point_data[next_point_id]['x'], point_data[next_point_id]['y']
                moment = calculate_moment(x1, y1, x2, y2, force_n1, force_n2)
                candidate_id = data['candidate_id']

                if candidate_id == 0:
                    # Check if point_id is the first of the pair for 0a or 0b
                    if point_id == candidate_0_points[0]:
                        group_key = "0a"
                    elif point_id == candidate_0_points[2]:
                        group_key = "0b"
                    else:
                        continue  # Skip if it's not the first point of the pair

                    moments[group_key]['t'] += force_t1 + force_t2
                    moments[group_key]['n'] += force_n1 + force_n2
                    moments[group_key]['moment'] += moment
                elif candidate_id != 0:
                    if candidate_id not in moments:
                        moments[candidate_id] = {'t': 0, 'n': 0, 'moment': 0}
                    moments[candidate_id]['t'] += force_t1 + force_t2
                    moments[candidate_id]['n'] += force_n1 + force_n2
                    moments[candidate_id]['moment'] += moment
            elif next_point_id:
                print(f"Warning: No force data for point_id {next_point_id}")

    # Find the highest numeric candidate ID for combination with group 0b
    numeric_candidate_ids = [cid for cid in moments.keys() if isinstance(cid, int)]
    candidate_n_id = max(numeric_candidate_ids) if numeric_candidate_ids else None

    # Retrieve values for candidate 1 and candidate n
    candidate_1_values = moments.get(1, {'t': 0, 'n': 0, 'moment': 0})  # Corrected to use integer key
    candidate_n_values = moments.get(candidate_n_id, {'t': 0, 'n': 0, 'moment': 0}) if candidate_n_id is not None else {
        't': 0, 'n': 0, 'moment': 0}

    # Considering that we have 4 contpts on 'ground' element, we divide them in two groups, first two will be summed
    # up with element 1 and second two with element n

    reactions_to_pass_static_final = [
        [
            candidate_1_values['t'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_1_values['n'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_1_values['moment'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ],
        [
            candidate_n_values['t'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_n_values['n'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_n_values['moment'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ]
    ]

    # print(reactions_to_pass_static_final)
    return reactions_to_pass_static_final


# This code chooses points of cand_id 0 (ground el) and els 1 and n.
def read_points():
    points_file = "./data/point.csv"
    points_bottom = []
    points_with_y_1 = []  # To store points with y-coordinates for candidate_id 1
    points_with_y_n = []  # To store points with y-coordinates for candidate_id n
    max_candidate_id = 0

    # Read points
    with open(points_file, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            candidate_id = int(row['candidate_id'])
            point_id = int(row['id'])
            y_coord = float(row['y'])

            max_candidate_id = max(max_candidate_id, candidate_id)
            if candidate_id == 0:
                points_bottom.append(point_id)
            elif candidate_id == 1:
                points_with_y_1.append((point_id, y_coord))
            elif candidate_id == max_candidate_id:
                points_with_y_n.append((point_id, y_coord))

    # Select the two points with the smallest y-coordinates for candidate_id 1
    points_with_y_1.sort(key=lambda x: x[1])
    points_bottom.extend([point[0] for point in points_with_y_1[:2]])

    # Select the two points with the smallest y-coordinates for candidate_id n
    points_with_y_n.sort(key=lambda x: x[1])
    points_bottom.extend([point[0] for point in points_with_y_n[:2]])

    return points_bottom


# # # f that saves all contact forces after static analysis # # #
def save_contacts_to_csv_static(contacts):
    filename = f"contactForces_static.csv"
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Assuming contacts is a list or similar structure with 72 force values
        for i in range(0, len(contacts), 2):  # Step through contacts two at a time
            contact_number = i // 2 + 1  # Calculate contact number (1-36)
            force1 = contacts[i]  # First force for this contact
            force2 = contacts[i + 1]  # Second force for this contact
            writer.writerow([contact_number, force1, force2])


# # # f that find 4 contps with lowest y coor not on el_0 and el_100, and saves them to csv after static an. # # #
def save_specific_contacts_to_csv_static(contacts):
    global model
    # Use the points obtained from read_points
    selected_points = read_points()

    # Filename for the specific contacts
    filename = f"reactions_of_blocks_static.csv"

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Process each specific contact ID from the selected points
        for contact_id in selected_points:
            # Calculate the index in the contacts list
            index = (contact_id - 1) * 2

            # Extract the forces for this contact
            if index < len(contacts) - 1:  # Check to avoid index out of range
                force1 = contacts[index]
                force2 = contacts[index + 1]
                writer.writerow([contact_id, force1, force2])

    # print(f"Contacts saved to {filename}")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Functions necessary to calculate and pass reactions to the OS.
# # # FOR TIME HISTORY ANALYSIS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Function that calculates resultant moment on the central point of an interface (shorter edge of voussoir).
# # # f that saves all contact forces after static analysis # # #
# # # f that saves all contact forces after one step of TH # # #
def sum_forces_to_reactions():
    global reactions_to_pass_final
    forces_filename = f"reactions_of_blocks_th.csv"
    points_filename = "data/point.csv"

    # Read forces and points data
    forces = {}
    point_data = {}
    candidate_0_points = []  # List to store points belonging to candidate ID 0

    # Read the forces file
    with open(forces_filename, mode='r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            contact_id = int(row[0])
            forces[contact_id] = {
                't': float(row[1]),
                'n': float(row[2])
            }

    # Read the points file and collect points for candidate ID 0
    with open(points_filename, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            point_id = int(row['id'])
            candidate_id = int(row['candidate_id'])
            point_data[point_id] = {
                'x': float(row['x']),
                'y': float(row['y']),
                'candidate_id': candidate_id,
                'counter_point': int(row['counter_point'])
            }
            if candidate_id == 0:
                candidate_0_points.append(point_id)

    # Initialize moments dictionary and include special groups for candidate ID 0
    moments = {"0a": {'t': 0, 'n': 0, 'moment': 0}, "0b": {'t': 0, 'n': 0, 'moment': 0}}

    # Calculate moments for each contact point
    for point_id, data in point_data.items():
        if point_id in forces:
            # Find the next consecutive point_id for the same candidate_id
            next_point_id = point_id + 1 if (point_id + 1) in point_data and \
                                            point_data[point_id + 1]['candidate_id'] == data['candidate_id'] else None

            if next_point_id and next_point_id in forces:
                force_n1, force_t1 = forces[point_id]['n'], forces[point_id]['t']
                force_n2, force_t2 = forces[next_point_id]['n'], forces[next_point_id]['t']
                x1, y1 = data['x'], data['y']
                x2, y2 = point_data[next_point_id]['x'], point_data[next_point_id]['y']
                moment = calculate_moment(x1, y1, x2, y2, force_n1, force_n2)
                candidate_id = data['candidate_id']

                if candidate_id == 0:
                    # Check if point_id is the first of the pair for 0a or 0b
                    if point_id == candidate_0_points[0]:
                        group_key = "0a"
                    elif point_id == candidate_0_points[2]:
                        group_key = "0b"
                    else:
                        continue  # Skip if it's not the first point of the pair

                    moments[group_key]['t'] += force_t1 + force_t2
                    moments[group_key]['n'] += force_n1 + force_n2
                    moments[group_key]['moment'] += moment
                elif candidate_id != 0:
                    if candidate_id not in moments:
                        moments[candidate_id] = {'t': 0, 'n': 0, 'moment': 0}
                    moments[candidate_id]['t'] += force_t1 + force_t2
                    moments[candidate_id]['n'] += force_n1 + force_n2
                    moments[candidate_id]['moment'] += moment
            elif next_point_id:
                print(f"Warning: No force data for point_id {next_point_id}")

    # Find the highest numeric candidate ID for combination with group 0b
    numeric_candidate_ids = [cid for cid in moments.keys() if isinstance(cid, int)]
    candidate_n_id = max(numeric_candidate_ids) if numeric_candidate_ids else None

    # Retrieve values for candidate 1 and candidate n
    candidate_1_values = moments.get(1, {'t': 0, 'n': 0, 'moment': 0})  # Corrected to use integer key
    candidate_n_values = moments.get(candidate_n_id, {'t': 0, 'n': 0, 'moment': 0}) if candidate_n_id is not None else {
        't': 0, 'n': 0, 'moment': 0}

    # Considering that we have 4 contpts on 'ground' element, we divide them in two groups, first two will be summed
    # up with element 1 and second two with element n

    reactions_to_pass_final = [
        [
            candidate_1_values['t'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_1_values['n'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_1_values['moment'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ],
        [
            candidate_n_values['t'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_n_values['n'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_n_values['moment'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ]
    ]

    # print(reactions_to_pass_final)
    return reactions_to_pass_final


# def save_contacts_to_csv(contacts, iteration):
def save_contacts_to_csv(contacts):
    # filename = f"contactForces_i{iteration}.csv" # Promenila si jer nemas vise te i od akcelerograma
    filename = f"contactForces_th.csv"
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Assuming contacts is a list or similar structure with 72 force values
        for i in range(0, len(contacts), 2):  # Step through contacts two at a time
            contact_number = i // 2 + 1  # Calculate contact number (1-36)
            force1 = contacts[i]  # First force for this contact
            force2 = contacts[i + 1]  # Second force for this contact
            writer.writerow([contact_number, force1, force2])


# def save_specific_contacts_to_csv(contacts, iteration):
def save_specific_contacts_to_csv(contacts):
    global model

    # Use the points obtained from read_points
    selected_points = read_points()

    # Filename for the specific contacts
    # filename = f"reactions_of_blocks_i{iteration}.csv"
    filename = f"reactions_of_blocks_th.csv"

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Process each specific contact ID from the selected points
        for contact_id in selected_points:
            # Calculate the index in the contacts list
            index = (contact_id - 1) * 2

            # Extract the forces for this contact
            if index < len(contacts) - 1:  # Check to avoid index out of range
                force1 = contacts[index]
                force2 = contacts[index + 1]
                writer.writerow([contact_id, force1, force2])

    # print(f"Contacts saved to {filename}")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Functions necessary to calculate and pass reactions to the OS.
# # # FOR TIME HISTORY ANALYSIS FAKE ACC IN FIRST CALL # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Function that calculates resultant moment on the central point of an interface (shorter edge of voussoir).
# # # f that saves all contact forces after static analysis # # #
# # # f that saves all contact forces after one step of TH and after TH with fake acc in first step # # #
def sum_forces_to_reactions_fakeacc():
    global reactions_to_pass_final_fakeacc
    forces_filename = f"reactions_of_blocks_th_fakeacc.csv"
    points_filename = "data/point.csv"

    # Read forces and points data
    forces = {}
    point_data = {}
    candidate_0_points = []  # List to store points belonging to candidate ID 0

    # Read the forces file
    with open(forces_filename, mode='r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            contact_id = int(row[0])
            forces[contact_id] = {
                't': float(row[1]),
                'n': float(row[2])
            }

    # Read the points file and collect points for candidate ID 0
    with open(points_filename, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            point_id = int(row['id'])
            candidate_id = int(row['candidate_id'])
            point_data[point_id] = {
                'x': float(row['x']),
                'y': float(row['y']),
                'candidate_id': candidate_id,
                'counter_point': int(row['counter_point'])
            }
            if candidate_id == 0:
                candidate_0_points.append(point_id)

    # Initialize moments dictionary and include special groups for candidate ID 0
    moments = {"0a": {'t': 0, 'n': 0, 'moment': 0}, "0b": {'t': 0, 'n': 0, 'moment': 0}}

    # Calculate moments for each contact point
    for point_id, data in point_data.items():
        if point_id in forces:
            # Find the next consecutive point_id for the same candidate_id
            next_point_id = point_id + 1 if (point_id + 1) in point_data and \
                                            point_data[point_id + 1]['candidate_id'] == data['candidate_id'] else None

            if next_point_id and next_point_id in forces:
                force_n1, force_t1 = forces[point_id]['n'], forces[point_id]['t']
                force_n2, force_t2 = forces[next_point_id]['n'], forces[next_point_id]['t']
                x1, y1 = data['x'], data['y']
                x2, y2 = point_data[next_point_id]['x'], point_data[next_point_id]['y']
                moment = calculate_moment(x1, y1, x2, y2, force_n1, force_n2)
                candidate_id = data['candidate_id']

                if candidate_id == 0:
                    # Check if point_id is the first of the pair for 0a or 0b
                    if point_id == candidate_0_points[0]:
                        group_key = "0a"
                    elif point_id == candidate_0_points[2]:
                        group_key = "0b"
                    else:
                        continue  # Skip if it's not the first point of the pair

                    moments[group_key]['t'] += force_t1 + force_t2
                    moments[group_key]['n'] += force_n1 + force_n2
                    moments[group_key]['moment'] += moment
                elif candidate_id != 0:
                    if candidate_id not in moments:
                        moments[candidate_id] = {'t': 0, 'n': 0, 'moment': 0}
                    moments[candidate_id]['t'] += force_t1 + force_t2
                    moments[candidate_id]['n'] += force_n1 + force_n2
                    moments[candidate_id]['moment'] += moment
            elif next_point_id:
                print(f"Warning: No force data for point_id {next_point_id}")

    # Find the highest numeric candidate ID for combination with group 0b
    numeric_candidate_ids = [cid for cid in moments.keys() if isinstance(cid, int)]
    candidate_n_id = max(numeric_candidate_ids) if numeric_candidate_ids else None

    # Retrieve values for candidate 1 and candidate n
    candidate_1_values = moments.get(1, {'t': 0, 'n': 0, 'moment': 0})  # Corrected to use integer key
    candidate_n_values = moments.get(candidate_n_id, {'t': 0, 'n': 0, 'moment': 0}) if candidate_n_id is not None else {
        't': 0, 'n': 0, 'moment': 0}

    # Considering that we have 4 contpts on 'ground' element, we divide them in two groups, first two will be summed
    # up with element 1 and second two with element n

    reactions_to_pass_final_fakeacc = [
        [
            candidate_1_values['t'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_1_values['n'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_1_values['moment'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ],
        [
            candidate_n_values['t'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_n_values['n'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_n_values['moment'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ]
    ]

    # print(reactions_to_pass_final)
    return reactions_to_pass_final_fakeacc


# def save_contacts_to_csv(contacts, iteration):
def save_contacts_to_csv_fakeacc(contacts):
    # filename = f"contactForces_i{iteration}.csv" # Promenila si jer nemas vise te i od akcelerograma
    filename = f"contactForces_th_fakeacc.csv"
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Assuming contacts is a list or similar structure with 72 force values
        for i in range(0, len(contacts), 2):  # Step through contacts two at a time
            contact_number = i // 2 + 1  # Calculate contact number (1-36)
            force1 = contacts[i]  # First force for this contact
            force2 = contacts[i + 1]  # Second force for this contact
            writer.writerow([contact_number, force1, force2])


# def save_specific_contacts_to_csv(contacts, iteration):
def save_specific_contacts_to_csv_fakeacc(contacts):
    global model

    # Use the points obtained from read_points
    selected_points = read_points()

    # Filename for the specific contacts
    # filename = f"reactions_of_blocks_i{iteration}.csv"
    filename = f"reactions_of_blocks_th_fakeacc.csv"

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Process each specific contact ID from the selected points
        for contact_id in selected_points:
            # Calculate the index in the contacts list
            index = (contact_id - 1) * 2

            # Extract the forces for this contact
            if index < len(contacts) - 1:  # Check to avoid index out of range
                force1 = contacts[index]
                force2 = contacts[index + 1]
                writer.writerow([contact_id, force1, force2])

    print(f"Contacts saved to {filename}")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Functions necessary to calculate and pass reactions to the OS.
# # # FOR PUSHOVER ANALYSIS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Function that calculates resultant moment on the central point of an interface (shorter edge of voussoir).
# # # f that saves all contact forces after static analysis # # #
# # # f that saves all contact forces after one step of TH # # #
# # # f that saves all contact forces after one step of pushover # # #
def sum_forces_to_reactions_pushover():
    global reactions_to_pass_pushover
    forces_filename = f"reactions_of_blocks_pushover.csv"
    points_filename = "data/point.csv"

    # Read forces and points data
    forces = {}
    point_data = {}
    candidate_0_points = []  # List to store points belonging to candidate ID 0

    # Read the forces file
    with open(forces_filename, mode='r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            contact_id = int(row[0])
            forces[contact_id] = {
                't': float(row[1]),
                'n': float(row[2])
            }

    # Read the points file and collect points for candidate ID 0
    with open(points_filename, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            point_id = int(row['id'])
            candidate_id = int(row['candidate_id'])
            point_data[point_id] = {
                'x': float(row['x']),
                'y': float(row['y']),
                'candidate_id': candidate_id,
                'counter_point': int(row['counter_point'])
            }
            if candidate_id == 0:
                candidate_0_points.append(point_id)

    # Initialize moments dictionary and include special groups for candidate ID 0
    moments = {"0a": {'t': 0, 'n': 0, 'moment': 0}, "0b": {'t': 0, 'n': 0, 'moment': 0}}

    # Calculate moments for each contact point
    for point_id, data in point_data.items():
        if point_id in forces:
            # Find the next consecutive point_id for the same candidate_id
            next_point_id = point_id + 1 if (point_id + 1) in point_data and \
                                            point_data[point_id + 1]['candidate_id'] == data['candidate_id'] else None

            if next_point_id and next_point_id in forces:
                force_n1, force_t1 = forces[point_id]['n'], forces[point_id]['t']
                force_n2, force_t2 = forces[next_point_id]['n'], forces[next_point_id]['t']
                x1, y1 = data['x'], data['y']
                x2, y2 = point_data[next_point_id]['x'], point_data[next_point_id]['y']
                moment = calculate_moment(x1, y1, x2, y2, force_n1, force_n2)
                candidate_id = data['candidate_id']

                if candidate_id == 0:
                    # Check if point_id is the first of the pair for 0a or 0b
                    if point_id == candidate_0_points[0]:
                        group_key = "0a"
                    elif point_id == candidate_0_points[2]:
                        group_key = "0b"
                    else:
                        continue  # Skip if it's not the first point of the pair

                    moments[group_key]['t'] += force_t1 + force_t2
                    moments[group_key]['n'] += force_n1 + force_n2
                    moments[group_key]['moment'] += moment
                elif candidate_id != 0:
                    if candidate_id not in moments:
                        moments[candidate_id] = {'t': 0, 'n': 0, 'moment': 0}
                    moments[candidate_id]['t'] += force_t1 + force_t2
                    moments[candidate_id]['n'] += force_n1 + force_n2
                    moments[candidate_id]['moment'] += moment
            elif next_point_id:
                print(f"Warning: No force data for point_id {next_point_id}")

    # Find the highest numeric candidate ID for combination with group 0b
    numeric_candidate_ids = [cid for cid in moments.keys() if isinstance(cid, int)]
    candidate_n_id = max(numeric_candidate_ids) if numeric_candidate_ids else None

    # Retrieve values for candidate 1 and candidate n
    candidate_1_values = moments.get(1, {'t': 0, 'n': 0, 'moment': 0})  # Corrected to use integer key
    candidate_n_values = moments.get(candidate_n_id, {'t': 0, 'n': 0, 'moment': 0}) if candidate_n_id is not None else {
        't': 0, 'n': 0, 'moment': 0}

    # Considering that we have 4 contpts on 'ground' element, we divide them in two groups, first two will be summed
    # up with element 1 and second two with element n

    reactions_to_pass_pushover = [
        [
            candidate_1_values['t'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_1_values['n'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_1_values['moment'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ],
        [
            candidate_n_values['t'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_n_values['n'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_n_values['moment'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ]
    ]

    # print(reactions_to_pass_final)
    return reactions_to_pass_pushover


# def save_contacts_to_csv(contacts, iteration):
def save_contacts_to_csv_pushover(contacts):
    # filename = f"contactForces_i{iteration}.csv" # Promenila si jer nemas vise te i od akcelerograma
    filename = f"contactForces_pushover.csv"
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Assuming contacts is a list or similar structure with 72 force values
        for i in range(0, len(contacts), 2):  # Step through contacts two at a time
            contact_number = i // 2 + 1  # Calculate contact number (1-36)
            force1 = contacts[i]  # First force for this contact
            force2 = contacts[i + 1]  # Second force for this contact
            writer.writerow([contact_number, force1, force2])


# def save_specific_contacts_to_csv(contacts, iteration):
def save_specific_contacts_to_csv_pushover(contacts):
    global model

    # Use the points obtained from read_points
    selected_points = read_points()

    # Filename for the specific contacts
    # filename = f"reactions_of_blocks_i{iteration}.csv"
    filename = f"reactions_of_blocks_pushover.csv"

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Process each specific contact ID from the selected points
        for contact_id in selected_points:
            # Calculate the index in the contacts list
            index = (contact_id - 1) * 2

            # Extract the forces for this contact
            if index < len(contacts) - 1:  # Check to avoid index out of range
                force1 = contacts[index]
                force2 = contacts[index + 1]
                writer.writerow([contact_id, force1, force2])

    # print(f"Contacts saved to {filename}")


#####################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Functions necessary to calculate and pass reactions to the OS.
# # # FOR LIMIT ANALYSIS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Function that calculates resultant moment on the central point of an interface (shorter edge of voussoir).
# # # f that saves all contact forces after static analysis # # #
# # # f that saves all contact forces after one step of TH # # #
# # # f that saves all contact forces after one step of pushover # # #
def sum_forces_to_reactions_LA():
    global reactions_to_pass_LA
    forces_filename = f"reactions_of_blocks_LA.csv"
    points_filename = "data/point.csv"

    # Read forces and points data
    forces = {}
    point_data = {}
    candidate_0_points = []  # List to store points belonging to candidate ID 0

    # Read the forces file
    with open(forces_filename, mode='r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            contact_id = int(row[0])
            forces[contact_id] = {
                't': float(row[1]),
                'n': float(row[2])
            }

    # Read the points file and collect points for candidate ID 0
    with open(points_filename, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            point_id = int(row['id'])
            candidate_id = int(row['candidate_id'])
            point_data[point_id] = {
                'x': float(row['x']),
                'y': float(row['y']),
                'candidate_id': candidate_id,
                'counter_point': int(row['counter_point'])
            }
            if candidate_id == 0:
                candidate_0_points.append(point_id)

    # Initialize moments dictionary and include special groups for candidate ID 0
    moments = {"0a": {'t': 0, 'n': 0, 'moment': 0}, "0b": {'t': 0, 'n': 0, 'moment': 0}}

    # Calculate moments for each contact point
    for point_id, data in point_data.items():
        if point_id in forces:
            # Find the next consecutive point_id for the same candidate_id
            next_point_id = point_id + 1 if (point_id + 1) in point_data and \
                                            point_data[point_id + 1]['candidate_id'] == data['candidate_id'] else None

            if next_point_id and next_point_id in forces:
                force_n1, force_t1 = forces[point_id]['n'], forces[point_id]['t']
                force_n2, force_t2 = forces[next_point_id]['n'], forces[next_point_id]['t']
                x1, y1 = data['x'], data['y']
                x2, y2 = point_data[next_point_id]['x'], point_data[next_point_id]['y']
                moment = calculate_moment(x1, y1, x2, y2, force_n1, force_n2)
                candidate_id = data['candidate_id']

                if candidate_id == 0:
                    # Check if point_id is the first of the pair for 0a or 0b
                    if point_id == candidate_0_points[0]:
                        group_key = "0a"
                    elif point_id == candidate_0_points[2]:
                        group_key = "0b"
                    else:
                        continue  # Skip if it's not the first point of the pair

                    moments[group_key]['t'] += force_t1 + force_t2
                    moments[group_key]['n'] += force_n1 + force_n2
                    moments[group_key]['moment'] += moment
                elif candidate_id != 0:
                    if candidate_id not in moments:
                        moments[candidate_id] = {'t': 0, 'n': 0, 'moment': 0}
                    moments[candidate_id]['t'] += force_t1 + force_t2
                    moments[candidate_id]['n'] += force_n1 + force_n2
                    moments[candidate_id]['moment'] += moment
            elif next_point_id:
                print(f"Warning: No force data for point_id {next_point_id}")

    # Find the highest numeric candidate ID for combination with group 0b
    numeric_candidate_ids = [cid for cid in moments.keys() if isinstance(cid, int)]
    candidate_n_id = max(numeric_candidate_ids) if numeric_candidate_ids else None

    # Retrieve values for candidate 1 and candidate n
    candidate_1_values = moments.get(1, {'t': 0, 'n': 0, 'moment': 0})  # Corrected to use integer key
    candidate_n_values = moments.get(candidate_n_id, {'t': 0, 'n': 0, 'moment': 0}) if candidate_n_id is not None else {
        't': 0, 'n': 0, 'moment': 0}

    # Considering that we have 4 contpts on 'ground' element, we divide them in two groups, first two will be summed
    # up with element 1 and second two with element n

    reactions_to_pass_LA = [
        [
            candidate_1_values['t'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_1_values['n'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_1_values['moment'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ],
        [
            candidate_n_values['t'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_n_values['n'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_n_values['moment'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ]
    ]

    # print(reactions_to_pass_final)
    return reactions_to_pass_LA


# def save_contacts_to_csv(contacts, iteration):
def save_contacts_to_csv_LA(contacts):
    # filename = f"contactForces_i{iteration}.csv" # Promenila si jer nemas vise te i od akcelerograma
    filename = f"contactForces_LA.csv"
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Assuming contacts is a list or similar structure with 72 force values
        for i in range(0, len(contacts), 2):  # Step through contacts two at a time
            contact_number = i // 2 + 1  # Calculate contact number (1-36)
            force1 = contacts[i]  # First force for this contact
            force2 = contacts[i + 1]  # Second force for this contact
            writer.writerow([contact_number, force1, force2])


# def save_specific_contacts_to_csv(contacts, iteration):
def save_specific_contacts_to_csv_LA(contacts):
    global model

    # Use the points obtained from read_points
    selected_points = read_points()

    # Filename for the specific contacts
    # filename = f"reactions_of_blocks_i{iteration}.csv"
    filename = f"reactions_of_blocks_LA.csv"

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Process each specific contact ID from the selected points
        for contact_id in selected_points:
            # Calculate the index in the contacts list
            index = (contact_id - 1) * 2

            # Extract the forces for this contact
            if index < len(contacts) - 1:  # Check to avoid index out of range
                force1 = contacts[index]
                force2 = contacts[index + 1]
                writer.writerow([contact_id, force1, force2])

    # print(f"Contacts saved to {filename}")


#####################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Functions necessary to calculate and pass reactions to the OS.
# # # FOR SETTLEMENT ANALYSIS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # Function that calculates resultant moment on the central point of an interface (shorter edge of voussoir).
# # # f that saves all contact forces after static analysis # # #
# # # f that saves all contact forces after one step of TH # # #
# # # f that saves all contact forces after one step of pushover # # #
# # # f that saves all contact forces after one step of settlements # # #
def sum_forces_to_reactions_settle():
    global reactions_to_pass_settle_r
    forces_filename = f"reactions_of_blocks_settle.csv"
    points_filename = "data/point.csv"

    # Read forces and points data
    forces = {}
    point_data = {}
    candidate_0_points = []  # List to store points belonging to candidate ID 0

    # Read the forces file
    with open(forces_filename, mode='r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            contact_id = int(row[0])
            forces[contact_id] = {
                't': float(row[1]),
                'n': float(row[2])
            }

    # Read the points file and collect points for candidate ID 0
    with open(points_filename, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            point_id = int(row['id'])
            candidate_id = int(row['candidate_id'])
            point_data[point_id] = {
                'x': float(row['x']),
                'y': float(row['y']),
                'candidate_id': candidate_id,
                'counter_point': int(row['counter_point'])
            }
            if candidate_id == 0:
                candidate_0_points.append(point_id)

    # Initialize moments dictionary and include special groups for candidate ID 0
    moments = {"0a": {'t': 0, 'n': 0, 'moment': 0}, "0b": {'t': 0, 'n': 0, 'moment': 0}}

    # Calculate moments for each contact point
    for point_id, data in point_data.items():
        if point_id in forces:
            # Find the next consecutive point_id for the same candidate_id
            next_point_id = point_id + 1 if (point_id + 1) in point_data and \
                                            point_data[point_id + 1]['candidate_id'] == data['candidate_id'] else None

            if next_point_id and next_point_id in forces:
                force_n1, force_t1 = forces[point_id]['n'], forces[point_id]['t']
                force_n2, force_t2 = forces[next_point_id]['n'], forces[next_point_id]['t']
                x1, y1 = data['x'], data['y']
                x2, y2 = point_data[next_point_id]['x'], point_data[next_point_id]['y']
                moment = calculate_moment(x1, y1, x2, y2, force_n1, force_n2)
                candidate_id = data['candidate_id']

                if candidate_id == 0:
                    # Check if point_id is the first of the pair for 0a or 0b
                    if point_id == candidate_0_points[0]:
                        group_key = "0a"
                    elif point_id == candidate_0_points[2]:
                        group_key = "0b"
                    else:
                        continue  # Skip if it's not the first point of the pair

                    moments[group_key]['t'] += force_t1 + force_t2
                    moments[group_key]['n'] += force_n1 + force_n2
                    moments[group_key]['moment'] += moment
                elif candidate_id != 0:
                    if candidate_id not in moments:
                        moments[candidate_id] = {'t': 0, 'n': 0, 'moment': 0}
                    moments[candidate_id]['t'] += force_t1 + force_t2
                    moments[candidate_id]['n'] += force_n1 + force_n2
                    moments[candidate_id]['moment'] += moment
            elif next_point_id:
                print(f"Warning: No force data for point_id {next_point_id}")

    # Find the highest numeric candidate ID for combination with group 0b
    numeric_candidate_ids = [cid for cid in moments.keys() if isinstance(cid, int)]
    candidate_n_id = max(numeric_candidate_ids) if numeric_candidate_ids else None

    # Retrieve values for candidate 1 and candidate n
    candidate_1_values = moments.get(1, {'t': 0, 'n': 0, 'moment': 0})  # Corrected to use integer key
    candidate_n_values = moments.get(candidate_n_id, {'t': 0, 'n': 0, 'moment': 0}) if candidate_n_id is not None else {
        't': 0, 'n': 0, 'moment': 0}

    # Considering that we have 4 contpts on 'ground' element, we divide them in two groups, first two will be summed
    # up with element 1 and second two with element n

    reactions_to_pass_settle_r = [
        [
            candidate_1_values['t'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_1_values['n'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_1_values['moment'] + moments.get("0a", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ],
        [
            candidate_n_values['t'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['t'],
            candidate_n_values['n'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['n'],
            candidate_n_values['moment'] + moments.get("0b", {'t': 0, 'n': 0, 'moment': 0})['moment']
        ]
    ]

    # print(reactions_to_pass_final)
    return reactions_to_pass_settle_r


# def save_contacts_to_csv(contacts, iteration):
def save_contacts_to_csv_settle(contacts):
    # filename = f"contactForces_i{iteration}.csv" # Promenila si jer nemas vise te i od akcelerograma
    filename = f"contactForces_settle.csv"
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Assuming contacts is a list or similar structure with 72 force values
        for i in range(0, len(contacts), 2):  # Step through contacts two at a time
            contact_number = i // 2 + 1  # Calculate contact number (1-36)
            force1 = contacts[i]  # First force for this contact
            force2 = contacts[i + 1]  # Second force for this contact
            writer.writerow([contact_number, force1, force2])


# def save_specific_contacts_to_csv(contacts, iteration):
def save_specific_contacts_to_csv_settle(contacts):
    global model

    # Use the points obtained from read_points
    selected_points = read_points()

    # Filename for the specific contacts
    # filename = f"reactions_of_blocks_i{iteration}.csv"
    filename = f"reactions_of_blocks_settle.csv"

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['ContactNumber', 't', 'n'])  # Column headers

        # Process each specific contact ID from the selected points
        for contact_id in selected_points:
            # Calculate the index in the contacts list
            index = (contact_id - 1) * 2

            # Extract the forces for this contact
            if index < len(contacts) - 1:  # Check to avoid index out of range
                force1 = contacts[index]
                force2 = contacts[index + 1]
                writer.writerow([contact_id, force1, force2])

    # print(f"Contacts saved to {filename}")


############# Global var for different model states #################################
model = None
# global saved_model_state
# saved_model_state = None
global init_state
init_state = None
global last_committed_state
last_committed_state = None
global random_state
random_state = None
global temp_saved_model_state
temp_saved_model_state = None
global last_commited_random_state
last_commited_random_state = None

# ###### GLobal var for a copy of a model for parallel application of displ and acc ####
# global rel_displaced_state
# rel_displaced_state = None
# global init_model_copy
# init_model_copy = None
global cycle_number
cycle_number = 0

global dt
dt = None

global theta
theta = None
# dt = 0.01
# theta = 0.7
############# Global var for reactions to be passed to OS ###########################
current_contacts = []
global reactions_to_pass_final  # Reactions that are updated after each i; and passed to OS as reactions after TH
reactions_to_pass_final = []
global reactions_to_pass_static_final
reactions_to_pass_static_final = []
global reactions_to_pass_final_fakeacc
reactions_to_pass_final_fakeacc = []  # Global variable to store the final static reactions
global reactions_to_pass_pushover
reactions_to_pass_pushover = []
global reactions_to_pass_LA
reactions_to_pass_LA = []
global reactions_to_pass_settle
reactions_to_pass_settle = []
global reactions_to_pass_settle_r
reactions_to_pass_settle_r = []
############# Global var to check result.txt - which analysis will be done ##########
global result
result = "0"

global load_multiplier
load_multiplier = 0
############ Global var to check if we are in the first call ##################
global first_call
first_call = True

############ Global var to take acceleration from OS ######################################
global global_accel_val
############ Global var to store the load multiplier from previous step of settlements analysis
############ needed for the scaling of displacements
global old_force
old_force = None
global new_force
new_force = None
########### Global var for pushover ##########################
global Dmax
Dmax = None
global incr
incr = None
global nSteps
nSteps = None
global controlled_node
controlled_node = None
global controlled_dof
controlled_dof = None
global next_iteration
next_iteration = None
global current_iteration
current_iteration = None


########## Global var for stopping the analysis on arch


############ Function to check the value in analysis_result.txt ###########################
def check_analysis_result():
    global result
    try:
        with open("analysis_result.txt", "r") as file:
            result = file.read().strip()
            return result
    except FileNotFoundError:
        return None

########### Function to initialize the model ##############################################
def init():
    global model
    print("Init py starts")
    if model is None:
        # print("pre fje u init, dt i theta", dt, theta)
        # get_parameters_TH_analysis([dt, theta])
        # print("posle fje u init, dt i theta", dt, theta)
        model = Model_TH_OS()  # Initialize the model
        print("Model initialized.")
        # set_dimension(2)
        # _data_dir = pathlib.Path(__file__).resolve().parent
        # # dt = dt
        # # theta = theta
        # model.from_csv(_data_dir / "data", dt, theta)  # Import geometry here
        # # model.from_csv(_data_dir / "data")  # Import geometry here
        # cal_anta_id(model.contps)


def pass_dt_theta(new_dt, new_theta):
    # Use the global keyword to indicate which variables you are referring to
    global dt, theta

    # Check if the vectors are not empty and update the global variables with the first element
    if new_dt and isinstance(new_dt, (list, tuple)):  # Check if it's list-like and not empty
        dt = new_dt[0]  # Assign the first element of new_dt to the global dt
    else:
        print("Warning: new_dt is empty or not a list/tuple.")

    if new_theta and isinstance(new_theta, (list, tuple)):  # Check if it's list-like and not empty
        theta = new_theta[0]  # Assign the first element of new_theta to the global theta
    else:
        print("Warning: new_theta is empty or not a list/tuple.")

    # Optionally print the new values to confirm they are received correctly
    print("Updated global dt and theta in Python:", dt, theta)


def set_domain(type):
    global model, dt, theta
    print(f"Inside set_domain with type (path): {type}")

    if model is None:
        print("Error: Model is not initialized. Please initialize before setting the domain.")
        return

    set_dimension(2)
    # Convert type to a Path object for more robust path handling
    # data_dir = pathlib.Path(type)
    # print("in set domain f in py data dir is: ", data_dir)
    data_dir = pathlib.Path(type) / "data"
    print(f"In set_domain, full path to data directory is: {data_dir}")
    if not data_dir.exists():
        print(f"Error: Data directory does not exist: {data_dir}")
        return

    model.from_csv(data_dir, dt, theta)

    # for e_id, e in model.elems.items():
    #     print(f"Restored dl and ll for element {e_id}: dl={e.dl}, ll={e.ll}")
    cal_anta_id(model.contps)


def get_parameters_pushover():
    # Path to the file
    filename = "pushover_parameters.txt"

    # Initialize the variables to None
    Dmax = None
    incr = None
    nSteps = None
    controlled_node = None
    controlled_dof = None

    # Dictionary to map keywords to variables
    parameters = {}

    # Open the file and read lines
    try:
        with open(filename, 'r') as file:
            for line in file:
                # Strip whitespace and split each line by space
                parts = line.strip().split()
                if len(parts) == 2:  # Ensure there are exactly two parts
                    key, value = parts[0], parts[1]
                    # Convert the value based on its nature
                    try:
                        # Try to convert to float first
                        float_value = float(value)
                        # If it's an integer, use int, otherwise keep float
                        if float_value.is_integer():
                            value = int(float_value)
                        else:
                            value = float_value
                    except ValueError:
                        pass  # Keep as string if it's not a number
                    # Assign the converted value to the correct variable
                    parameters[key] = value

        # Map values from the dictionary to variables, if present
        Dmax = parameters.get('Dmax', None)
        incr = parameters.get('incr', None)
        nSteps = parameters.get('nSteps', None)
        controlled_node = parameters.get('controlled_node', None)
        controlled_dof = parameters.get('controlled_dof', None)

    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")

    # Return the values
    return Dmax, incr, nSteps, controlled_node, controlled_dof


# Example usage
# get_parameters_pushover('pushover_parameters.txt')
#
# # Now you can access the variables globally, for example:
# print(Dmax, incr, nSteps, controlled_node, controlled_dof)


def get_parameters_settlements():
    # Path to the file
    filename = "settlements_parameters.txt"

    # Initialize the variables to None
    Dmax = None
    incr = None
    nSteps = None
    controlled_node = None
    controlled_dof = None
    displaced_node = None
    displaced_dof = None

    # Dictionary to map keywords to variables
    parameters = {}

    # Open the file and read lines
    try:
        with open(filename, 'r') as file:
            for line in file:
                # Strip whitespace and split each line by space
                parts = line.strip().split()
                if len(parts) == 2:  # Ensure there are exactly two parts
                    key, value = parts[0], parts[1]
                    # Convert the value based on its nature
                    try:
                        # Try to convert to float first
                        float_value = float(value)
                        # If it's an integer, use int, otherwise keep float
                        if float_value.is_integer():
                            value = int(float_value)
                        else:
                            value = float_value
                    except ValueError:
                        pass  # Keep as string if it's not a number
                    # Assign the converted value to the correct variable
                    parameters[key] = value

        # Map values from the dictionary to variables, if present
        Dmax = parameters.get('Dmax', None)
        incr = parameters.get('incr', None)
        nSteps = parameters.get('nSteps', None)
        controlled_node = parameters.get('controlled_node', None)
        controlled_dof = parameters.get('controlled_dof', None)
        displaced_node = parameters.get('displaced_node', None)
        displaced_dof = parameters.get('displaced_dof', None)

    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")

    # Return the values
    return Dmax, incr, nSteps, controlled_node, controlled_dof, displaced_node, displaced_dof


init()


# ########### Function to save rel displaced model and restore rel displaced model #######################################
# # Function to save the state
# def save_rel_displaced_state():
#     global rel_displaced_state, init_model_copy
#     rel_displaced_state = copy.deepcopy(init_model_copy)
#
# # Function to restore the state
# def restore_rel_displaced_state():
#     global rel_displaced_state, init_model_copy
#     init_model_copy = copy.deepcopy(rel_displaced_state)

########### Function to save committed model and restore committed model #######################################
# def save_committed_state(from_temp=False):
#     global saved_model_state, temp_saved_model_state
#     if from_temp:
#         saved_model_state = copy.deepcopy(temp_saved_model_state)
#     else:
#         saved_model_state = copy.deepcopy(model)


########## Function to save scaled state when not first_call and before TH ############################################


def restore_committed_state():
    global model, last_committed_state

    model = copy.deepcopy(last_committed_state)


########### Function to save temporary model and restore last temporary model #######################################

def save_temp_state():
    global temp_saved_model_state, model

    temp_saved_model_state = copy.deepcopy(model)


def restore_temp_state():
    global model, temp_saved_model_state

    if temp_saved_model_state is not None:
        model = copy.deepcopy(temp_saved_model_state)
    else:
        print("No temporary saved state to restore.")


########### Function to save initial model and restore initial model #######################################

def save_init_state():
    global init_state, model
    init_state = copy.deepcopy(model)


def restore_init_state():
    global init_state, model
    model = copy.deepcopy(init_state)


########### Function to restore random state for DS model ##################################################

def restore_random_state():
    global random_state, init_state, first_call, last_saved_random_state

    if first_call:
        random_state = copy.deepcopy(init_state)
    else:
        random_state = copy.deepcopy(last_saved_random_state)
        model = copy.deepcopy(random_state)


def save_random_state():
    global model, random_state

    random_state = copy.deepcopy(model)


def recover_last_commited_random():
    global last_commited_random_state, random_state

    random_state = copy.deepcopy(last_commited_random_state)


########## Update function

def update(disp, accel):
    # Importing and Initializing Global Variables

    global model, current_contacts, result, first_call, global_accel_val, temp_saved_model_state, random_state, last_commited_random_state, cycle_number
    global init_state, last_committed_state, dt, theta, next_iteration, current_iteration, old_force, new_force,  reactions_to_pass_static_final
    global reactions_to_pass_settle_r, reactions_to_pass_final
    # Convert disp and accel to numpy arrays
    disp = np.array(disp)
    accel = np.array(accel)

    print("In the beginning of update function, first call is:", first_call)

    # Update the global acceleration value with the first component of the accel array
    global_accel_val = accel[0]

    # Check the value in analysis_result.txt and set result
    result = check_analysis_result()
    print('Result inside update function:', result)

    # Exit if result is "0" (No analysis to be performed)
    if result == "0":
        return


    elif result == "1":  # Static Analysis

        print("Performing static analysis.")

        Aglobal = cal_A_global_2d(model.elems, model.contps)

        # Step 1: Save the initial dl and ll values for elements of type 's'
        saved_dl_ll = {}  # Dictionary to store the initial dl and ll values for each element

        # for e_id, e in model.elems.items():
        #     print(f"res 1, element {e_id}: dl={e.dl}, ll={e.ll}")


        for e_id, e in model.elems.items():
            if e.type == 's':
                saved_dl_ll[e_id] = {'dl': e.dl.copy(), 'll': e.ll.copy()}  # Store the current dl and ll values
                e.dl = [0, 0, 0]  # Set to zero as per current code
                e.ll = [0, 0, 0]
                print(f"Ako je blok s tipa, element {e_id}: dl={e.dl}, ll={e.ll}")

        static_reactions, static_disp, static_conv = solve_elastic_infinitefc_associative_2d(model.elems, model.contps,
                                                                                             Aglobal=Aglobal)
        # Save reactions and contacts
        save_contacts_to_csv_static(static_reactions)
        save_specific_contacts_to_csv_static(static_reactions)

        # Sum forces to reactions
        reactions_to_pass_static_final = sum_forces_to_reactions_static()
        print("Reactions to pass (static final):", reactions_to_pass_static_final)

        # Step 2: Restore the initial dl and ll values before saving the state
        for e_id, e in model.elems.items():
            if e_id in saved_dl_ll:  # Check if the element was of type 's'
                e.dl = saved_dl_ll[e_id]['dl']  # Restore the original dl
                e.ll = saved_dl_ll[e_id]['ll']  # Restore the original ll
                print(f"Restored dl and ll for element {e_id}: dl={e.dl}, ll={e.ll}")

        for elem_id, elem in model.elems.items():
            # Print the ID of the element and its displacement
            print(f"Static an, model, Element {elem_id} displacement: {elem.displacement}")

        # Save the initial state after static analysis
        save_temp_state()
        save_init_state()  # Da li je ok ovo?


    elif result == "2":  # Dynamic Analysis

        print("Performing dynamic analysis.")
        acc_val = accel[0]
        print("Acc_vel is", acc_val)

        # Preserve first_call status specifically for calc_vel_after_OS
        is_first_call_for_calc_vel = first_call and not np.isclose(accel[0], 0.0)
        print(f"Before if first_call: first_call = {first_call}")

        # Handle the first call with zero acceleration
        # Check if acceleration vector is zero
        if np.all(np.isclose(accel[0], 0.0)):

            print("Acceleration vector is zero.")
            # If it's the first call and accel is zero, do not proceed with the analysis
            if first_call:
                print("First call detected but acceleration is zero, skipping real analysis. Fake acc = 0.0001")
                acc_val = 0.0001
                model.add_gr_mot([acc_val, 0.0, 0.0])
                model.add_init_vel()
                model.add_init_vel_part(dt)
                model.tot_load_step0_acc(1)
               # print('Acceleration from OS is - before cal_gap:', acc_val)
                cal_gap_2d(model.contps)

                r, contact_forces = solve_TH_rigid_associative_2d_FB(model.elems, model.contps)

                save_contacts_to_csv_fakeacc(contact_forces)
                save_specific_contacts_to_csv_fakeacc(contact_forces)
                reactions_to_pass_final_fakeacc = sum_forces_to_reactions_fakeacc()
                # print('Reactions to be passed, final, fake acc:', reactions_to_pass_final_fakeacc())


        else:

            print("Theta is", theta)
            print("Dt is", dt)
            # If here, accel is not zero
            # Read from tcl which node and dof you're controlling for the DS analysis

            Dmax, incr, nSteps, controlled_node, controlled_dof, displaced_node, displaced_dof = get_parameters_settlements()

            # U om pravcu ide incr ? u py je incr u minusu ako ide na levo
            print("increment is,:", incr)
            print("Else, controlled node for DS is", controlled_node)

            if first_call:
                print("else,1st, before DS, find reactions of static an.", reactions_to_pass_static_final)
                print("First dynamic analysis call with non-zero acceleration.")

                # ----- Run one random DS analysis -----
                print("First call, random state")
                restore_random_state()

                for elem_id, elem in random_state.elems.items():
                    # Print the ID of the element and its displacement
                    print(f"1st call, before random DS an, random state, Element {elem_id} displacement: {elem.displacement}")

                # !!!!! Check acc_val on which side it is? Depending then  on the sign, the reaction should be in s.ll.

                for e in random_state.elems.values():
                    if e.type == 's':
                        print("Ako je blok s tipa")
                        e.dl = [0, -reactions_to_pass_static_final[0][1], 0]
                        #e.dl = [0, -3.6, 0]
                        print("Dead load on block s is", e.dl)
                        e.ll = [reactions_to_pass_static_final[0][1], 0, 0]
                        #e.ll = [3.6, 0, 0]
                        print("Live load on block s is", e.ll)

                    elif 'brick' in e.type:  # If e.type contains 'brick'
                        print("Element type contains 'brick'")
                        e.ll = [0, 0, 0]  # Live load set to [0,0,0]
                        print("Live load on brick element is", e.ll)
                        print("Dead load on brick element is", e.dl)

                random_state.pre_check_push_over()
                current_iteration_r = 0

                Aglobal_copy = cal_A_global_2d(random_state.elems, random_state.contps)

                # Perform the analysis on the copied model

                forces, displacements_incre, solsta, contact_forces_settle_r = solve_pushover_rigid_friction_associative_2d_new(
                    random_state.elems, random_state.contps, controlled_node, (incr, 0, 0), 1,
                    current_iteration=current_iteration_r)

                if solsta in [solsta.dual_infeas_cer, solsta.prim_infeas_cer]:
                    print("Solution is infeasible (either primal or dual). Stopping analysis.")
                  # return
                #
                save_contacts_to_csv_settle(contact_forces_settle_r)
                save_specific_contacts_to_csv_settle(contact_forces_settle_r)  # Save specific contacts
                reactions_to_pass_settle_r = sum_forces_to_reactions_settle()
                print("1st call; reactions of random DS:",
                      reactions_to_pass_settle_r)  # Random state se vec sam updateovao, ne treba save_random_state jer to sajvuje model u random

                for elem_id, elem in random_state.elems.items():
                    # Print the ID of the element and its displacement

                    print(f"1st call, random DS an, random state, Element {elem_id} displacement: {elem.displacement}")

                # ----- Run the iteration of TH -----

                # Make a deep copy of the initial model before modifying it
                print("Inside first_call block with non-zero acceleration TH an.")
                print('Acceleration from OS is - inside block firstcall:', accel)

                restore_init_state()

                for elem_id, elem in model.elems.items():
                    # Print the ID of the element and its displacement
                    print(
                        f"1st call, 1st iter of TH ,after restore init, model, Element {elem_id} displacement: {elem.displacement}")

                # # This block will only execute for the first time with non-zero acceleration


                model.add_gr_mot([accel[0], 0.0, 0.0])
                model.add_init_vel()
                model.add_init_vel_part(dt)
                model.tot_load_step0_acc()

                for e_id, e in model.elems.items():
                    if e.type == 's':
                        #e.dl = [0, 0, 0]  # Set to zero as per current code
                        e.ll = [0, 0, 0]
                        e.totload = [0, e.totload[1], 0]
                    #print(f"TH 1st call: Loads on block are element {e_id}: dl={e.dl}, ll={e.ll}, totload={e.totload}")


                for e_id, e in model.elems.items():
                    # Assuming all elements have the attributes dl, ll, and totload
                    print(f"Element {e_id}: dl={e.dl}, ll={e.ll}, sl={e.sl}, totload={e.totload}")

                # Save the state before committing any changes
                save_temp_state()
                print('Acceleration from OS is - before cal_gap:', accel)
                cal_gap_2d(model.contps)

                r, contact_forces, solsta, y = solve_TH_rigid_associative_2d_FB(model.elems, model.contps)

                if solsta in [solsta.dual_infeas_cer, solsta.prim_infeas_cer]:
                    print("Firstcall - Solution is infeasible (either primal or dual). Stopping analysis. Solsta is", solsta)
                    return

                    # analysis_should_stop = True

                    # sys.exit("FORCE TH Infeasible solution encountered. Exiting program.")

                # # ----- Code you need for the visualization -----
                #
                # result_data = result  # No need to convert to string just yet
                #
                # # Group every three elements together
                #
                # grouped_data = [result_data[i:i + 3] for i in range(0, len(result_data), 3)]
                #
                # # Convert list of groups into a string for writing to file
                #
                # result_string = '\n'.join(['\t'.join(map(str, group)) for group in grouped_data])
                #
                # # Filename in the current working directory
                #
                # file_path = 'center_disp_i1.txt'
                #
                # # Open the file in write mode and save the data
                #
                # with open(file_path, 'w') as file:
                #
                #     file.write(result_string)

                # ----------------------------------------------------

                save_contacts_to_csv(contact_forces)
                save_specific_contacts_to_csv(contact_forces)
                reactions_to_pass_final = sum_forces_to_reactions()
                model.calc_vel_after_OS(theta, dt,
                                        is_first_call_for_calc_vel)  # Assuming theta = 0.7 and dt = 0.001 as per your setup

                for elem_id, elem in model.elems.items():
                    # Print the ID of the element and its displacement
                    print(f"1st call, TH an, model state, Element {elem_id} displacement: {elem.displacement}")

                # ------------------------------------------------

                save_temp_state()


            else:

                # Handle subsequent calls with non-zero acceleration

                # Increment cycle number each time this block is executed

                # ----- Restore last_commited_state -----

                restore_committed_state()
                for elem_id, elem in model.elems.items():
                    # Print the ID of the element and its displacement
                    print(f"Else, restore com state, model, Element {elem_id} displacement: {elem.displacement}")

                # ----- Adjust the model based on the previous DS analysis before the TH step -----

                print("Read out files")
                disp_l_node_out, disp_r_node_out, relative_disp_out = read_last_displacement_values()
                print("Values from out files, left, right", disp_l_node_out, disp_r_node_out)
               #
                print("recover last commited random")
                recover_last_commited_random()
                for elem_id, elem in random_state.elems.items():  # u ovom trenutku isto kao model.msM DA NIJE VEC KAO ds MODEL
                    # Print the ID of the element and its displacement
                    print(f"Else, random_state, check before scaling Element {elem_id} displacement: {elem.displacement}")

                last_TH_model_disp1 = random_state.elems[1].displacement[
                    0]  # !!!! Check if last_committed_state or last random. S markom last_committed_state
                print("Last TH model disp1", last_TH_model_disp1)

                last_TH_model_disp21 = random_state.elems[21].displacement[0]
                print("Last TH model disp21", last_TH_model_disp21)

                # last_TH_model_relative = last_TH_model_disp21 - last_TH_model_disp1 #!!!!! If cond, koji je veci displ
                last_TH_model_relative = last_TH_model_disp1 - last_TH_model_disp21

                if last_TH_model_relative != 0:
                    print(
                        f"Disp r.out {disp_r_node_out} - disp l.out{disp_l_node_out} / DS elem1{last_TH_model_disp1} - DS elem21 {last_TH_model_disp21}")
                    scale_factor = (disp_r_node_out - disp_l_node_out) / last_TH_model_relative  # !!!!! Da li treba abs?
                    print("Scale factor is",scale_factor)

                else:
                    scale_factor = 0.0

                DS_reaction = 1/(disp_r_node_out - disp_l_node_out)
                print("DS_reactions is", DS_reaction)
                #print("Scale factor is", scale_factor)
                # OVde je bilo pre recover_last_commited_random()
                number_of_items = len(model.elems)

                for idx in range(number_of_items):

                    if scale_factor != 0.0:
                        # print("Skaliranje u konacni model")

                        # ----- Sums both the DS and TH displacements -----

                        model.elems[idx].displacement[0] = model.elems[idx].displacement[0] + scale_factor * \
                                                           random_state.elems[idx].displacement[0]

                        model.elems[idx].displacement[1] = model.elems[idx].displacement[1] + scale_factor * \
                                                           random_state.elems[idx].displacement[
                                                               1]  # !!!!! Da li skalirati rotaciju

                #print("Ovde si sada")
                #
                for elem_id, elem in model.elems.items():
                    # Print the ID of the element and its displacement

                    print(f"Else, after scaling, final model, Element {elem_id} displacement: {elem.displacement}")
                print("ovde smo")
            ######## DA LI OVDE TREBA save_random_state()
                save_random_state()  # !!!!! Proveri da li da skaliramo sile ili ti solver vec izbaci prave vrednosti jer jer na skairanom modelu

                # save_random_state saveuje model trenutni, tj ovaj finalni.

                # ----- Create DS analysis for the next iteration -----


                random_state.pre_check_push_over()


                for e in random_state.elems.values():
                    if e.type == 's':
                       # print("Ako je blok s tipa")
                        e.dl = [0, -DS_reaction, 0]
                       # print("Dead load on block s is", e.dl)
                        e.ll = [DS_reaction, 0, 0]
                       # print("Live load on block s is", e.ll)

                    elif 'brick' in e.type:  # If e.type contains 'brick'
                        # print("Element type contains 'brick'")
                        e.ll = [0, 0, 0]  # Live load set to [0,0,0]
                       # print("Live load on brick element is", e.ll)
                       # print("Dead load on brick element is", e.dl)

                current_iteration_r = 0

                Aglobal_copy = cal_A_global_2d(random_state.elems, random_state.contps)

                for elem_id, elem in random_state.elems.items():
                    # Print the ID of the element and its displacement

                    print(f"After total model, before second DS an, Element {elem_id} displacement: {elem.displacement}")

                forces, displacements_incre, solsta, contact_forces_settle_r = solve_pushover_rigid_friction_associative_2d_new(
                    random_state.elems, random_state.contps, controlled_node, (incr, 0, 0), 1,
                    current_iteration=current_iteration_r)


                # # First, check for infeasibility in solsta as you already did

                # if solsta in [solsta.dual_infeas_cer, solsta.prim_infeas_cer]:

                #     print("ELSE DS Solution is infeasible (either primal or dual). Stopping analysis. Solsta is", solsta)

                #     return None

                save_contacts_to_csv_settle(contact_forces_settle_r)
                save_specific_contacts_to_csv_settle(contact_forces_settle_r)  # Save specific contacts
                reactions_to_pass_settle_r = sum_forces_to_reactions_settle()
                #print("ELSE DS reactions_to_pass_settle_r", reactions_to_pass_settle_r)

                for elem_id, elem in random_state.elems.items():
                    # Print the ID of the element and its displacement
                    print(f"After total model, after second DS an, Element {elem_id} displacement: {elem.displacement}")


                ##### Ovde se sacuva total model ove iteracije i nastavlja se dalje, ali sledeca iteracija ovo ima za svoj random_state.

                # ----- Run one TH iteration -----

                #print("Subsequent call with non-zero acceleration.")
                model.add_gr_mot([accel[0], 0.0, 0.0])

                # Any specific operations for subsequent updates go here

                # For example, if you have operations that should only be done after the first non-zero acceleration update

                model.tot_load_steps_acc_OS(dt)
                for e_id, e in model.elems.items():
                    if e.type == 's':
                        # e.dl = [0, 0, 0]  # Set to zero as per current code
                        e.ll = [0, 0, 0]
                        e.totload = [0, e.totload[1], 0]
                    print(f"Ako je blok s tipa, element {e_id}: dl={e.dl}, ll={e.ll}, totload={e.totload}")

                print('Acceleration from OS is - before cal_gap:', accel)
                cal_gap_2d(model.contps)

                r, contact_forces, solsta, y = solve_TH_rigid_associative_2d_FB(model.elems, model.contps)

                if solsta in [solsta.dual_infeas_cer, solsta.prim_infeas_cer]:
                    print("ELSE TH Solution is infeasible (either primal or dual). Stopping analysis.")
                    return

                    # analysis_should_stop = True

                    # sys.exit("ELSE FORCE TH AN infeasible solution encountered. Exiting program.")

                # Next, check if all values in contact_forces_settle_r are zero

                if all(force == 0.0 for force in forces):
                    print("ELSE TH All contact forces are zero. Stopping analysis.")

                    # sys.exit("ELSE TH All contact forces are zero. Exiting program.")

                # ----- Visualisation -----
                # result_data = y # No need to convert to string just yet
                # # Group every three elements together
                # grouped_data = [result_data[i:i + 3] for i in range(0, len(result_data), 3)]
                # # Convert list of groups into a string for writing to file
                # result_string = '\n'.join(['\t'.join(map(str, group)) for group in grouped_data])
                # # Filename in the current working directory
                # # Use cycle_number in the file name
                # file_path = f'center_disp_i{cycle_number}.txt'  # Dynamic file name based on the cycle number
                # # Open the file in write mode and save the data
                # with open(file_path, 'w') as file:
                #     file.write(result_string)

                # ----------------------------------------------------------------------------------------

                save_contacts_to_csv(contact_forces)
                save_specific_contacts_to_csv(contact_forces)
                reactions_to_pass_final = sum_forces_to_reactions()
                model.calc_vel_after_OS(theta, dt,
                                        is_first_call_for_calc_vel)  # Assuming theta = 0.7 and dt = 0.001 as per your setup

                print("Else, TH, reactions", reactions_to_pass_final)

                for elem_id, elem in model.elems.items():
                    # Print the ID of the element and its displacement
                    print(f"After total model, after second DS and after TH cycle *last step*, Element {elem_id} displacement: {elem.displacement}")

                # -------------------------

                save_temp_state()  # This now acts as saving the state after analysis

    #     elif result == "3":  # Perform the pushover analysis code here
    #         Dmax, incr, nSteps, controlled_node, controlled_dof = get_parameters_pushover()
    #         # if controlled_dof == 1:
    #         #     disp_incre = (incr, 0, 0)
    #         # elif controlled_dof == 2:
    #         #     disp_incre = (0, incr, 0)
    #         # elif controlled_dof == 3:
    #         #     disp_incre = (0, 0, incr)
    #         # else:
    #         #     disp_incre = (0, 0, 0)
    #         load_incr = incr
    #         print("Load increment is", load_incr)
    #         # component_index = np.nonzero(disp_incre)
    #         print("Starting pushover analysis")
    #         #
    #         # contacts = []
    #         # forces = []
    #         # displacements = []
    #         if first_call:
    # ############################################################ With pushover code
    #             restore_init_state()
    #             model.pre_check_push_over()
    #             save_temp_state()
    #             current_iteration=0
    #             print("Inside first call, before solver, next_iteration_cont is", next_iteration)
    #             Aglobal = cal_A_global_2d(model.elems, model.contps)
    #             # forces, displacements_incre, solsta, contact_forces = solve_pushover_rigid_friction_associative_2d_new(
    #             #     model.elems, model.contps, controlled_node, (incr, 0, 0), 1 , current_iteration=current_iteration)
    #             # print("Inside first call, before solver, next_iteration_cont is", next_iteration_cont)
    #             forces_PO1, displacement_PO1, contact_forces_P01 = solve_pushover_elastic_associative_2d_new(model.elems,
    #                         model.contps, controlled_node, load_incr, 1, Aglobal=Aglobal, current_iteration=0)
    #             print("Update 3 after opt, print contact_forces", contact_forces_P01)
    #             # if solsta in [solsta.dual_infeas_cer, solsta.prim_infeas_cer]:
    #             #     print("Solution is infeasible (either primal or dual). Stopping analysis.")
    #             #     return
    #
    #             save_contacts_to_csv_LA( contact_forces_P01)
    #             save_specific_contacts_to_csv_LA( contact_forces_P01)  # Save specific contacts
    #             reactions_to_pass_LA = sum_forces_to_reactions_LA()
    #             print("reactions to pass pushover:", reactions_to_pass_LA)
    #             # print('Reactions to be passed, final:', reactions_to_pass_final)
    #             # for contp_id, contp in model.contps.items():
    #             #     # Print the ID of the contact point and its gap value
    #             #     print(f"Contact point {contp_id} displacement: {contp.displacement}")
    #             save_temp_state()
    #         # # **************************************************************plot
    #         else:
    #             #cycle_number += 1  # Increment cycle number each time this block is executed
    #             restore_committed_state()
    #             print("Inside else pushover")
    #             for contp_id, contp in model.contps.items():
    #                 # Print the ID of the contact point and its gap value
    #                 print(f"Contact point {contp_id} displacement: {contp.displacement}")
    #
    #             model.pre_check_push_over()
    #             print("In else doing pushover")
    #             print("Next iteration is",next_iteration)
    #             print("Current iteration is",current_iteration)
    #             print("NSteps is", nSteps)
    #             # next_iteration = 1
    #             print("Disp incre je", )
    #             #save_temp_state()disp_incre
    #             # forces, displacements_incre, solsta, contact_forces_n = solve_pushover_rigid_friction_associative_2d_new(
    #             #     model.elems, model.contps, controlled_node, (incr, 0, 0), 1, current_iteration=current_iteration)
    #
    #             forces, displacement,  contact_forces_n = solve_pushover_elastic_associative_2d_new(model.elems,
    #                                                                 model.contps, controlled_node, load_incr, 1, Aglobal,
    #                                                                 current_iteration=current_iteration)
    #             print("IN ELSE AFTER SOLVER")
    #             print("contact forces are in else:,", contact_forces_n)
    #             save_contacts_to_csv_LA(contact_forces_n)
    #             save_specific_contacts_to_csv_LA(contact_forces_n)  # Save specific contacts
    #             reactions_to_pass_LA = sum_forces_to_reactions_LA()
    #             print("reactions to pass pushover:", reactions_to_pass_LA)
    #             # if solsta in [solsta.dual_infeas_cer, solsta.prim_infeas_cer]:
    #             #     print("Solution is infeasible (either primal or dual). Stopping analysis.")
    #             #     reactions_to_pass_LA = [0,0,0,0,0,0, 0,0,0,0,0,0]
    #             #     return
    #                         # for contp_id, contp in model.contps.items():
    #             #     # Print the ID of the contact point and its gap value
    #             #     print(f"Contact point {contp_id} displacement: {contp.displacement}")
    #
    #             # contacts= [item for sublist in contacts for item in sublist]
    #             # save_contacts_to_csv_pushover(contact_forces)
    #             # save_specific_contacts_to_csv_pushover(contact_forces)  # Save specific contacts
    #             # reactions_to_pass_pushover = sum_forces_to_reactions_LA()
    #             # print("reactions to pass pushover in else :", reactions_to_pass_pushover)
    #             save_temp_state()

    # elif result == "4":  # Perform the settlements analysis code here
    #     # Simple Python script to print the contents of the file
    #
    #
    #     Dmax, incr, nSteps, controlled_node, controlled_dof, displaced_node, displaced_dof = get_parameters_settlements()
    #     # U om pravcu ide incr ? u py je incr u minusu ako ide na levo
    #     print("increment is,:",incr)
    #     if first_call:
    #         ############################################################ With pushover code
    #         restore_init_state()
    #
    #         model.pre_check_push_over()
    #         save_temp_state()
    #         current_iteration = 0
    #         print("Inside 4 first call , before solver, next_iteration_cont is", next_iteration)
    #         Aglobal = cal_A_global_2d(model.elems, model.contps)
    #         #new_force = get_last_value_from_file("load_multiplier.txt")
    #         print("In 4 first call, new force read is", new_force)
    #         forces, displacements_incre, solsta, contact_forces_settle = solve_pushover_rigid_friction_associative_2d_new(
    #                      model.elems, model.contps, controlled_node, (incr, 0, 0), 1 , current_iteration=current_iteration)
    #         # print("Inside first call, before solver, next_iteration_cont is", next_iteration_cont)
    #         print("Update 4, first call, forces (load mult)is", forces)
    #         old_force = forces[0]
    #         print("Update 4 first call, old force is", old_force)
    #         print("Update 4 after opt, print contact_forces", contact_forces_settle)
    #         # if solsta in [solsta.dual_infeas_cer, solsta.prim_infeas_cer]:
    #         #     print("Solution is infeasible (either primal or dual). Stopping analysis.")
    #         #     return
    #
    #         save_contacts_to_csv_settle(contact_forces_settle)
    #         save_specific_contacts_to_csv_settle(contact_forces_settle)  # Save specific contacts
    #         reactions_to_pass_settle= sum_forces_to_reactions_settle()
    #         print("reactions to pass settlements:", reactions_to_pass_settle)
    #         # print('Reactions to be passed, final:', reactions_to_pass_final)
    #         # for contp_id, contp in model.contps.items():
    #         #     # Print the ID of the contact point and its gap value
    #         #     print(f"Contact point {contp_id} displacement: {contp.displacement}")
    #         save_temp_state()
    #         # # **************************************************************plot
    #     else:
    #         # cycle_number += 1  # Increment cycle number each time this block is executed
    #         restore_committed_state()
    #         print("Inside else settlements")
    #         # for contp_id, contp in model.contps.items():
    #         #     # Print the ID of the contact point and its gap value
    #         #     print(f"Contact point {contp_id} displacement: {contp.displacement}")
    #
    #         model.pre_check_push_over()
    #         print("In else doing settlement")
    #         print("Next iteration is", next_iteration)
    #         print("Current iteration is", current_iteration)
    #         print("NSteps is", nSteps)
    #         # next_iteration = 1
    #         print("Disp incre je", incr )
    #         load_multiplier = check_last_value()
    #         new_force = load_multiplier
    #         print("In 4 else, new force read is", new_force)
    #         # save_temp_state()disp_incre
    #         forces, displacements_incre, solsta, contact_forces_s = solve_pushover_rigid_friction_associative_2d_new(
    #             model.elems, model.contps, controlled_node, (incr, 0, 0), 1, old_force, new_force,  current_iteration=current_iteration)
    #         #new_force = forces[0]
    #         # if solsta in [solsta.dual_infeas_cer, solsta.prim_infeas_cer]:
    #         #     print("Solution is infeasible (either primal or dual). Stopping analysis.")
    #         #     return
    #
    #         print("IN ELSE AFTER SOLVER settlements")
    #         print("contact forces are in else:,", contact_forces_s)
    #         save_contacts_to_csv_settle(contact_forces_s)
    #         save_specific_contacts_to_csv_settle(contact_forces_s)  # Save specific contacts
    #         reactions_to_pass_settle = sum_forces_to_reactions_settle()
    #         print("reactions to pass settlements:", reactions_to_pass_settle)
    #
    #         save_temp_state()

    else:
        print("Invalid value in analysis_result.txt")
    # first_call = False


# Function that gets called by OpenSees when the state is committed

def commit_state_in_python(successVec):
    global first_call, last_committed_state, temp_saved_model_state, result, next_iteration, current_iteration, old_force, \
        new_force, last_commited_random_state, cycle_number

    cycle_number += 1

    if first_call and result == "1":

        print("Res 1 - Success vector is:", successVec)
        last_committed_state = temp_saved_model_state

        # ----- Visualisation -----
        result_data = []  # Initialize an empty list to store the element data
        # Loop through each element and format its displacements
        for elem_id, elem in last_committed_state.elems.items():
            # Convert displacement list to a string with tab-separated values
            displacement_str = '\t'.join([str(val) for val in elem.displacement])
            result_data.append(displacement_str)  # Add the formatted displacement to result_data

        # Join the displacements, one per line
        result_string = '\n'.join(result_data)

        file_path = f'plotArch/Arch_Calderini/plot_elemDispl/center_disp_i{cycle_number}.txt'

        with open(file_path, 'w') as file:
            file.write(result_string)

        first_call = True

    if first_call and result == "2":

        print("Res 2 - Success vector is:", successVec)

        last_committed_state = temp_saved_model_state
        last_commited_random_state = random_state

        # ----- Visualisation -----
        result_data = []  # Initialize an empty list to store the element data

        # Loop through each element and format its displacements
        for elem_id, elem in last_committed_state.elems.items():
            # Convert displacement list to a string with tab-separated values
            displacement_str = '\t'.join([str(val) for val in elem.displacement])
            result_data.append(displacement_str)  # Add the formatted displacement to result_data

        # Join the displacements, one per line
        result_string = '\n'.join(result_data)

        file_path = f'plotArch/Arch_Calderini/plot_elemDispl/center_disp_i{cycle_number}.txt'

        with open(file_path, 'w') as file:
            file.write(result_string)

        first_call = False

    elif (not first_call) and result == "2":

        last_committed_state = temp_saved_model_state
        last_commited_random_state = random_state

        # ----- Visualisation -----
        result_data = []  # Initialize an empty list to store the element data

        # Loop through each element and format its displacements
        for elem_id, elem in last_committed_state.elems.items():
            # Convert displacement list to a string with tab-separated values
            displacement_str = '\t'.join([str(val) for val in elem.displacement])
            result_data.append(displacement_str)  # Add the formatted displacement to result_data

        # Join the displacements, one per line
        result_string = '\n'.join(result_data)
        # result_string = '\n'.join(['\t'.join(group) for group in grouped_data])

        file_path = f'plotArch/Arch_Calderini/plot_elemDispl/center_disp_i{cycle_number}.txt'

        with open(file_path, 'w') as file:
            file.write(result_string)

        first_call = False

    if first_call and result == "3":
        print("Res 3 - Success vector is:", successVec)
        last_committed_state = temp_saved_model_state
        first_call = False
        current_iteration += 1
        print("After commit, current iteration is:", current_iteration)
        print("After commit, next iteration is:", next_iteration)

    elif (not first_call) and result == "3":
        last_committed_state = temp_saved_model_state
        current_iteration += 1
        print("Commit in general case, current iteration is:", current_iteration)
        print("Commit in general case, next iteration is:", next_iteration)

    if first_call and result == "4":
        print("Res 4 - Success vector is:", successVec)
        last_committed_state = temp_saved_model_state
        # old_force = forces[-1] if forces else old_force
        first_call = False
        current_iteration += 1
        # new_force = get_last_value_from_file("load_multiplier.txt")
        print("Inside commit, first call, new force is", new_force)
        print("After commit, current iteration is:", current_iteration)
        print("After commit, next iteration is:", next_iteration)

    elif (not first_call) and result == "4":
        last_committed_state = temp_saved_model_state
        load_multiplier = check_last_value()
        new_force = load_multiplier
        print('load mult inside update f is,', load_multiplier)
        print("Inside commit, else, new force is", new_force)
        current_iteration += 1
        print("Commit in general case, current iteration is:", current_iteration)
        print("Commit in general case, next iteration is:", next_iteration)

    return


def get_resisting_force():
    global result, global_accel_val, first_call, reactions_to_pass_final, reactions_to_pass_static_final, \
        reactions_to_pass_LA, reactions_to_pass_pushover, reactions_to_pass_settle, reactions_to_pass_settle_r

    print('result isnide get res force', result)
    # print("before if in get res force , reactions to pass LA", reactions_to_pass_LA)
    # Initialize force with a default value to ensure it always has a value
    # force = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    if result == "0":
        # When the result is "0" it is before Model has been created, use the existing list
        print('we are inside if res ==0; get res force', result)
        force = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    elif result == "1":
        # When the result is "1" (static analysis), use the new list
        # force = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        print("Inside getResForce 1")
        initial_force = [item for sublist in reactions_to_pass_static_final for item in sublist]
        # print("initial force:",initial_force)
        # print("reactions to pass static final", reactions_to_pass_static_final)
        # print("Reactions to pass static final are:", reactions_to_pass_static_final)
        # # After adding zeroes for the 3d case:

        force = [

            initial_force[0],
            0,
            initial_force[1],
            0,
            initial_force[2],
            0,  # Mz_1

            initial_force[3],
            0,
            initial_force[4],  # Pz_2
            0,
            initial_force[5],
            # My_2
            0

        ]


    elif result == '2':

        # Check if acceleration vector is zero

        if np.all(np.isclose(global_accel_val, 0.0)):

            # If it's the first call and accel is zero

            if first_call:
                print("First call detected but acceleration is zero, skipping analysis. Return with fake acc")

                # force = [-896.7, 8151.568848719811, -597.6734472479266, 4622.498833490921, 6680.047769971034, -643.9079868353318]  # Return zeros as specified

                # = [Nx1, Vx1, 0, 0, 0, Mz]

                initial_force = [item for sublist in reactions_to_pass_static_final for item in sublist]

                force = [

                    initial_force[0],

                    0,

                    initial_force[1],

                    0,

                    initial_force[2],

                    0,  # Mz_1

                    initial_force[3],

                    0,

                    initial_force[4],  # Pz_2

                    0,

                    initial_force[5],

                    # My_2

                    0

                ]

                # force = [item for sublist in reactions_to_pass_final_fakeacc for item in sublist]

        #     else:

        #         # If accel is zero but not the first call, proceed with normal reactions logic

        #         print("Acceleration is zero, using normal reactions.")

        #         force = [item for sublist in reactions_to_pass_final for item in sublist]  # Use existing logic for normal reactions

        else:

            # Handling for non-zero acceleration or subsequent calls.

            # if all(force == 0.0 for force in contact_forces_settle_r):

            #     print("get_res_force: All contact forces are zero. Stopping analysis.")

            #     sys.exit("All contact forces are zero. Exiting program.")

            # print("In grf f, check reactions_to_pass_settle_r", reactions_to_pass_settle_r)
            #
            # # Check if reactions_to_pass_settle_r contains all zeros
            #
            # if reactions_to_pass_settle_r is not None and np.all(np.isclose(reactions_to_pass_settle_r, 0.0)):
            #     print("reactions_to_pass_settle_r contains all zeros. Stopping the analysis.")
            #
            #     return [np.nan] * 12  # Return NaNs or some special value to signal the analysis should stop

            print("Handling non-zero acceleration or subsequent calls.")

            initial_force = [item for sublist in reactions_to_pass_final for item in
                             sublist]  # Apply dynamic analysis results.

            # Check if all forces are zeros

            # if np.all(np.isclose(initial_force, 0.0)):
            #     print("All forces are zero. Returning from function.")
            #
            #     return  # Exit the function  if all forces are zero

            force = [

                initial_force[0],

                0,

                initial_force[1],

                0,

                initial_force[2],

                0,  # Mz_1

                initial_force[3],

                0,

                initial_force[4],  # Pz_2

                0,

                initial_force[5],

                # My_2

                0

            ]

    elif result == '3':
        if first_call:
            print('getResForce first call', result)
            initial_force = [item for sublist in reactions_to_pass_LA for item in sublist]
            print("Initial force in getResForce  pfirst call", initial_force)

            force = [
                initial_force[0],
                0,
                initial_force[1],
                0,
                initial_force[2],
                0,  # Mz_1

                initial_force[3],
                0,
                initial_force[4],  # Pz_2
                0,
                initial_force[5],
                # My_2
                0]

            print(force)

        else:
            print('getResForce in else', result)
            initial_force = [item for sublist in reactions_to_pass_LA for item in sublist]
            print("Initial force in getResForce  else", initial_force)

            force = [
                initial_force[0],
                0,
                initial_force[1],
                0,
                initial_force[2],
                0,  # Mz_1

                initial_force[3],
                0,
                initial_force[4],  # Pz_2
                0,
                initial_force[5],
                # My_2
                0]

            print(force)

    elif result == '4':
        if first_call:
            print('getResForce first call settlements', result)
            print('getResForce first call settle reactions', reactions_to_pass_settle)
            initial_force = [item for sublist in reactions_to_pass_settle for item in sublist]
            print("Initial force in getResForce  first call", initial_force)

            force = [
                initial_force[0],
                0,
                initial_force[1],
                0,
                initial_force[2],
                0,  # Mz_1

                initial_force[3],
                0,
                initial_force[4],  # Pz_2
                0,
                initial_force[5],
                # My_2
                0]

            print(force)

        else:
            print('getResForce in else', result)
            print('getResForce else settle reactions', reactions_to_pass_settle)
            initial_force = [item for sublist in reactions_to_pass_settle for item in sublist]
            print("Initial force in getResForce  else", initial_force)

            force = [
                initial_force[0],
                0,
                initial_force[1],
                0,
                initial_force[2],
                0,  # Mz_1

                initial_force[3],
                0,
                initial_force[4],  # Pz_2
                0,
                initial_force[5],
                # My_2
                0]

            print(force)
        # force = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    else:
        print("Invalid value in analysis_result.txt")
        # force = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]

    return force


def get_initial_stiff():
    # stiff = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    # return stiff.flatten().tolist() # flatten() returns a 1D array, easier to process in C++
    # Absolute path to the .txt file
    #file_path = "C:/Users/ivana/BackUp/Arch_OpenSees/Arch3D_12/Basic_Stiffness_3d_ElBeam.txt"
    file_path = "C:/Users/ivana/BackUp/Arch_OpenSees/Arch3D_12/Basic_Stiffness_3d_newE.txt"
    matrix = np.loadtxt(file_path, comments="#")
    # Load the matrix from the file using numpy's loadtxt function
    # matrix = np.loadtxt(file_path)
    # Flatten the matrix into a 1D array
    # If the matrix is not 1D (e.g., 3x3), flatten it into a 1D array
    if matrix.ndim > 1:
        flattened_array = matrix.flatten()
    else:
        flattened_array = matrix

    return flattened_array.tolist()


mat = get_initial_stiff()  # This will always read stiff_mat_10.txt
print(mat)

def get_mass_mom_inertia():

    # Absolute path to the .txt file
    file_path = "C:/Users/ivana/BackUp/python_dev_TH/kinematic/src/Tests/Arch_onButtresses_Dyn_Block2d_Fig5/massMatrix_3D.txt"
    # Read the data from the file
    matrix = np.loadtxt(file_path)

    # Extract the diagonal values from the matrix
    diagonal_values = np.diag(matrix)

    # Reshape the diagonal values to a 2D array with one row
    reshaped_values = diagonal_values.reshape(1, -1)

    # Flatten the 2D array to a 1D array and convert it to a list
    return reshaped_values.flatten().tolist()


# # Example usage
inertia_values = get_mass_mom_inertia()
print(inertia_values)
# # # #
# disp = [0, 0, 0]
# # # # accel = [0.0, 0.0, 0.0]
# # # # # # # # #
# # # # update(disp, accel)
# # # # a = get_resisting_force()
# # # # print(a)
# # # # # accel=[0.3, 0.0, 0.0]
# accel=[0.200057, 0.166348, -0.438867]
# result = 3
#
# update(disp, accel)
# accel = [-0.56,-0.6,0.6]
# update(disp, accel)
# accel =[0.3, 0.1, -0.867]
# update(disp, accel)
# # # print("van updatea, ", reactions_to_pass_final)
# # result = check_analysis_result()
# # # print('na kraju result:',result)
# # a = get_resisting_force()
# # print(a)
#
#





















































