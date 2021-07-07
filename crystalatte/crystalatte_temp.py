#!/usr/bin/env python

#
# @BEGIN LICENSE
#
# CrystaLattE: The tool for the automated calculation of crystal lattice 
# energies.
#
# Copyright (c) 2017-2020 
# Carlos H. Borca
# Brandon W. Bakr
# Lori A. Burns
# C. David Sherrill
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of CrystaLattE.
#
# CrystaLattE is free software; you can redistribute it and/or modify
# it under the tesms of the GNU Lesser General Public License as 
# published by the Free Software Foundation, version 3.
#
# CrystaLattE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public 
# License along with CrystaLattE; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
# 02110-1301 USA.
#
# @END LICENSE
#

# Import standard Python modules.
import itertools
import math
import multiprocessing
import numpy as np
import os
import sys
import shutil
import time

# Import parts of Psi4.
import psi4

# Import QCEelemental and QCEngine
import qcelemental as qcel
import qcengine as qcng

# Import parts of PyCIFRW
from CifFile import CifFile

# Debugging module
from pprint import pprint
from scipy.spatial import distance_matrix

# ======================================================================
def input_parser(in_f_name):
    """Reads a CrystaLattE input file, sets the arguments for the main()
    function, and calls main() to start the execution of the program.
    
    Arguments:
    <str> in_f_name
        Input filename.
    """

    keywords = {}
    keywords["nmers_up_to"] = 2 
    keywords["cif_a"] = 5
    keywords["cif_b"] = 5
    keywords["cif_c"] = 5
    keywords["bfs_thresh"] = 1.2
    keywords["uniq_filter"] = "ChSEV"
    keywords["r_cut_com"] = 10.0 
    keywords["r_cut_monomer"] = 12.0 
    keywords["r_cut_dimer"] = 10.0 
    keywords["r_cut_trimer"] = 8.0
    keywords["r_cut_tetramer"] = 6.0 
    keywords["r_cut_pentamer"] = 4.0 
    keywords["cle_run_type"] = ["test"]
    keywords["psi4_method"] = "HF/STO-3G"
    keywords["psi4_bsse"] = "cp"
    keywords["psi4_memory"] = "500 MB"
    keywords["qca_server_uri"] = None
    keywords["qca_mode"] = 'run'
    keywords["verbose"] = 1

    with open(in_f_name, "r") as input_file:
        
        non_blank_lines = []

        for line_inp_f in input_file:

            if line_inp_f == "\n":
                continue

            elif line_inp_f.startswith("#"):
                continue
            
            else:
                non_blank_lines.append(line_inp_f)

        for non_blank_line in non_blank_lines:

            split_line = non_blank_line.split("=")
            keyword_name = split_line[0].lower().strip()

            try:
                keyword_value = split_line[1].split('\n')[0]
            
            except IndexError:
                print("\nERROR: Invalid input file. Check that keywords and values are separated by an equal sign.\n")
                sys.exit()
           
            if keyword_name == "cle_run_type":
                keyword_value = keyword_value.lower()
                keyword_value = keyword_value.replace(" ", "").split("+")
                   
            elif keyword_name in ["uniq_filter", "psi4_bsse", "psi4_memory", "psi4_method", "cif_input", "cif_output", "qca_server_uri","qca_mode"]:
                keyword_value = keyword_value.strip()

            elif keyword_name in ["nmers_up_to", "verbose"]:
                
                try:
                    keyword_value = int(keyword_value)

                except ValueError:
                    print("\nERROR: Invalid input file. Check that given {} is an integer.\n".format(keyword_name))
                    sys.exit()
            
            elif keyword_name in ["cif_a", "cif_b", "cif_c"]:
                
                try:
                    keyword_value = int(keyword_value)
                    
                    if (keyword_value%2) == 0:
                        keyword_value += 1

                except ValueError:
                    print("\nERROR: Invalid input file. Check that given {} is an integer.\n".format(keyword_name))
                    sys.exit()

            else:
                try:
                    keyword_value = float(keyword_value)

                except ValueError:
                    print("\nERROR: Invalid input file. Check that given {} is an float.\n".format(keyword_name))
                    sys.exit()

            keywords[keyword_name] = keyword_value

    # Check CIF is provided.
    if "cif_input" not in keywords.keys():
        print("\nERROR: No CIF given as input.\n")
        sys.exit()

    # Check CIF extension.
    else:
        if not keywords["cif_input"].endswith(".cif"):
            print("\nERROR: Invalid CIF file name. Check that given file ends with a .cif extension.\n")
            sys.exit()

    # Attempt to create a filename for the supercell XYZ file.
    if "cif_output" not in keywords.keys():
    
        # Check proper CIF input filename.
        if keywords["cif_input"].endswith(".cif"):
            # Create supercell filename changing the extension of the input.
            keywords["cif_output"] = keywords["cif_input"][:-4] + ".xyz"

        else:
            print("\nERROR: Invalid CIF file name. Check that given file ends with a .cif extension.\n")
            sys.exit()

    # Print program header.
    if keywords["verbose"] >= 1:

        print_header()

    if keywords["verbose"] >= 2:

        print("CrystaLattE execution setup:\n")
        
        # Get the keys of the keywords dictionary, and put them on a list.
        kw_keys = list(keywords.keys())

        # Sort the list in decreasing order.
        kw_keys.sort()
        
        for kw_key in kw_keys:
            print("  {:15} = {}".format(kw_key, str(keywords[kw_key])))

    if "makefp" in keywords["cle_run_type"]:

        if len(keywords["cle_run_type"]) > 1:

            if "timings" not in keywords["cle_run_type"]: # TODO: What if there are more than two methods?!
                print("\nERROR: makefp mode cannot be run at the same time with any other mode, except timings.\n")
                sys.exit()
    
    elif "quiet" in keywords["cle_run_type"]:
        
        if (len(keywords["cle_run_type"]) < 2):
            print("\nERROR: quiet mode be must run together with psi4api or test mode.\n")
            sys.exit()
        
        elif "psithon" in keywords["cle_run_type"]:
            print("\nERROR: quiet and psithon modes cannot be run at the same time.\n")
            sys.exit()

        elif "libefpmbe" in keywords["cle_run_type"]:
            print("\nERROR: quiet and libefpmbe modes cannot be run at the same time.\n")
            sys.exit()

        elif "timings" in keywords["cle_run_type"]:
            print("\nERROR: if running quiet and timings modes together, psi4api mode must be included too.\n")
            sys.exit()

    elif "psi4api" in keywords["cle_run_type"]:

        if "psithon" in keywords["cle_run_type"]:
            print("\nERROR: psi4api and psithon modes cannot be run at the same time.\n")
            sys.exit()

        elif "libefpmbe" in keywords["cle_run_type"]:
            print("\nERROR: psi4api and libefpmbe modes cannot be run at the same time.\n")
            sys.exit()

        else:
            print("HELP!! I AM TRAPPED! elif psi4api in keywords[cle_run_type]:")

    if "timings" in keywords["cle_run_type"]:

        if (len(keywords["cle_run_type"]) < 2):
            print("\nERROR: timings mode be must run together with another mode.\n")
            sys.exit()

        else:
        
            func_str = "main("

            for kw_key in kw_keys:
                
                if type(keywords[kw_key]) == str:
                    func_str = func_str + "{}='{}', ".format(kw_key, str(keywords[kw_key]))
                
                else:
                    func_str = func_str + "{}={}, ".format(kw_key, str(keywords[kw_key])) 
            
            func_str = func_str[:-2] + ")"
            
            import cProfile as profile
            profile.run(func_str)

    else:
        main(
            keywords["cif_input"], 
            keywords["cif_output"],
            keywords["cif_a"], 
            keywords["cif_b"],
            keywords["cif_c"], 
            keywords["bfs_thresh"],
            keywords["uniq_filter"], 
            keywords["nmers_up_to"],
            keywords["r_cut_com"], 
            keywords["r_cut_monomer"], 
            keywords["r_cut_dimer"], 
            keywords["r_cut_trimer"], 
            keywords["r_cut_tetramer"], 
            keywords["r_cut_pentamer"], 
            keywords["cle_run_type"],
            keywords["psi4_method"],
            keywords["psi4_bsse"],
            keywords["psi4_memory"],
            keywords["qca_server_uri"],
            keywords["qca_mode"],
            keywords["verbose"])

    return keywords
    
# ======================================================================


# ======================================================================
def extract_element(label):
    """Converts an "_atom_type_label" into an element name.
    """

    elem2 = ['He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar','Ca','Sc','Ti',
             'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
             'Rb','Sr','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
             'Sb','Te','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',
             'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','Re','Os','Ir','Pt',
             'Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',
             'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']

    if (label[0:2] in elem2):
        return label[0:2]

    elem1 = ['H','B','C','N','O','F','P','S','K','V','Y','I','W','U']

    if (label[0] in elem1):
        return label[0]

    print("{}".format("~"*(shutil.get_terminal_size().columns)))
    print('WARNING: could not convert "%s" into element name!' % label)
    print("{}".format("~"*(shutil.get_terminal_size().columns)))
    return label
# ======================================================================


# ======================================================================
def write_xyz(atoms, box, f):
    """Writes a basic XYZ format file.
    """

    # Write the number of atoms.
    N = len(atoms)
    f.write('%d\n' % N)

    # Write a comment line with the box size.
    if (len(box) == 6):
        # box = (ax,bx,by,cx,cy,cz)
        ax = box[0]
        bx = box[1]
        by = box[2]
        cx = box[3]
        cy = box[4]
        cz = box[5]
        f.write('Crystal created from CIF file. Box vectors:')
        f.write(' a= %10.5f %10.5f %10.5f' % (ax, 0.0, 0.0))
        f.write(' b= %10.5f %10.5f %10.5f' % (bx, by, 0.0))
        f.write(' c= %10.5f %10.5f %10.5f\n' % (cx, cy, cz))
    else:
        # box = (ax,by,cz)
        f.write('Crystal created from CIF file. Box size:') 
        f.write(' %10.5f %10.5f %10.5f\n' % box)

    # Write atom data (units are Angstroms).
    # The argument "atoms" has format ('Si', x, y, z) for example
    for i in range(N):
        f.write('%-10s %10.6f %10.6f %10.6f\n' % atoms[i])
# =============================================================================


# =============================================================================
def read_cif(fNameIn):
    """Read CIF file, and extract the necessary info in the form of a 
    dictionary. E.g., the value of "_cell_volume" can be found with 
    data['_cell_volume'].
    """
    data = {}

    # Open the CIF file and read all the lines into a list of strings.
    try:
        f = open(fNameIn, 'r')
        lines = []

        for line in f:
            stripped = line.strip()

            if (len(stripped) > 0):  lines.append(stripped)

    except TypeError:
        print("\nERROR: A type error occured when opening or reading CIF file '{0}'".format(fNameIn))
        sys.exit()

    except NameError:
        print("\nERROR: Failed to open CIF file '{0}'".format(fNameIn))
        sys.exit()

    except FileNotFoundError:
        print("\nERROR: Failed to find CIF file '{0}'".format(fNameIn))
        sys.exit()

    except:
        print("\nERROR: An unknown error occured when opening or reading CIF file '{0}'".format(fNameIn))
        sys.exit()

    # Use the CifFile parser to extract the data. Although there might
    # be multiple data blocks, we'll only use the first one.
    
    #TODO: Known bug: Sometimes not all the blocks in the CIF files are
    #      read, and the following extraction fails.
    cif_file = CifFile(fNameIn)

    for db in cif_file:
        data_block = db
        break

    try:
        # Extract some parameters, and convert them to floats.
        data['_cell_length_a']    = float(data_block['_cell_length_a'].replace("(", "").replace(")", ""))
        data['_cell_length_b']    = float(data_block['_cell_length_b'].replace("(", "").replace(")", ""))
        data['_cell_length_c']    = float(data_block['_cell_length_c'].replace("(", "").replace(")", ""))
        data['_cell_angle_alpha'] = float(data_block['_cell_angle_alpha'].replace("(", "").replace(")", ""))
        data['_cell_angle_beta']  = float(data_block['_cell_angle_beta'].replace("(", "").replace(")", ""))
        data['_cell_angle_gamma'] = float(data_block['_cell_angle_gamma'].replace("(", "").replace(")", ""))
        data['_cell_volume']      = float(data_block['_cell_volume'].replace("(", "").replace(")", ""))

        # Get the symbolic operations that define the space group. In
        # a CIF file that's the part that looks like:
        
        # loop_
        # _symmetry_equiv_pos_as_xyz
        #   'x,y,z'
        #   'y,x,2/3-z'
        #   '-y,x-y,2/3+z'
        #   '-x,-x+y,1/3-z'
        #   '-x+y,-x,1/3+z'
        #   'x-y,-y,-z'
        
        # In some cases it is called:
        # "_space_group_symop_operation_xyz".
        data['_symmetry_equiv_pos_as_xyz'] = []

        try:
            xyz = data_block["_symmetry_equiv_pos_as_xyz"]

        except KeyError:
        
            try:
                xyz = data_block["_space_group_symop_operation_xyz"]
            
            except KeyError:
                print('\nERROR: Missing item in CIF file: need either \'_symmetry_equiv_pos_as_xyz\' or \'_space_group_symop_operation_xyz\'.')
                sys.exit()

        # Copy the x,y,z symmetry group operations. Remove the quotes 
        # if there are any.
        for op_xyz in xyz:

            if (op_xyz[0] == '\''):
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz[1:-1])
            
            else:
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz)


        # Add x,y,z of the atoms to "data", but make sure to convert
        # e.g. "0.1549(8)" to "0.1549".
        data['_atom_site_label'] = data_block['_atom_site_label']

        data['_atom_site_fract_x'] = []
        
        for str_x in data_block['_atom_site_fract_x']:
            data['_atom_site_fract_x'].append( float(str_x.split('(')[0]) )

        data['_atom_site_fract_y'] = []
        
        for str_y in data_block['_atom_site_fract_y']:
            data['_atom_site_fract_y'].append( float(str_y.split('(')[0]) )

        data['_atom_site_fract_z'] = []
        
        for str_z in data_block['_atom_site_fract_z']:
            data['_atom_site_fract_z'].append( float(str_z.split('(')[0]) )

    
    except KeyError as e:
        print('\nERROR: Missing item in file.')
        print(e)
        sys.exit()

    # Return the extracted data.
    return data
# =============================================================================


# ======================================================================
def cif_driver(cif_input, cif_output, cif_a, cif_b, cif_c, verbose=1):
    """Takes the name of a CIF input file and the name of a .xyz output
    file, as well as the number of replicas of the rectangular cell in
    each direction (A, B, and C). It then calls Read_CIF() and passes
    that information as arguments to generate an .xyz file of the
    supercell.

    Arguments:
    <str> cif_input
        CIF input filename.
    <str> cif_output
        XYZ output filename.
    <int> cif_a
        Number of replicas of the cartesian unit cell in `a` direction.
    <int> cif_b
        Number of replicas of the cartesian unit cell in `b` direction.
    <int> cif_c
        Number of replicas of the cartesian unit cell in `c` direction.
    <int> verbose
        Adjusts the level of detail of the printouts.
    """

    cif_arguments = ["", "-i", cif_input, "-o", cif_output, "-b", cif_a, cif_b, cif_c, "-r"]

    if verbose >= 2:
        print("\nGenerating the supercell .xyz file.")
        print("\nThe following arguments will be passed to the CIF reader script:")
        print("./Read_CIF.py" + " ".join(str(cif_argument) for cif_argument in cif_arguments) + "\n")

    return cif_arguments
# ======================================================================


# =============================================================================
def cif_main(args):

    # Default settings.
    fNameIn = ''
    fNameOut = ''
    Nx = 1
    Ny = 1
    Nz = 1
    make_rect_box = False
    
    
    # Read the arguments. We expect at least 4.
    if (len(args) <= 4):
        print_usage()
    
    i = 1
    while (i < len(args)):
    
        # Check if the name of the input file was given.
        if (args[i] == '-i'):
            
            # Make sure a file name is given.
            if (i+1 == len(args)):
                print('\nERROR: no input file name given')
            
            fNameIn = args[i+1]
            i = i + 2
    
        # Check if the name of the output file was given.
        elif (args[i] == '-o'):
    
            # Make sure a file name is given.
            if (i+1 == len(args)):
                print('\nERROR: no output file name given')
    
            # Check we have a valid file extension.
            fNameOut = args[i+1]
            unknown = True
    
            for ext in ['.xyz', '.lammpstrj', '.gro', '.cif']:
                if (fNameOut.endswith(ext)):
                    unknown = False
    
            if (unknown):
                print('\nERROR: unknown file extension of output file')
    
            i = i + 2
    
        # Check if the box size was given.
        elif (args[i] == '-b'):
    
            # Make sure 3 integers are given.
            if (i+3 >= len(args)):
                print('\nERROR: need 3 integers to indicate box size')
    
            Nx = int(args[i+1])
            Ny = int(args[i+2])
            Nz = int(args[i+3])
    
            if (Nx == 0  or  Ny == 0  or  Nz == 0):
                print('\nERROR: box size integers need to be larger than zero')
    
            i = i + 4
    
        # Check if the final configuration should be in a rectangular
        # shape, or in the same shape as the unit cell.
        elif (args[i] == '-r'):
    
            make_rect_box = True
            i = i + 1
    
    
        # Anything else is wrong.
        else:
            print('\nERROR: invalid argument "%s"' % args[i])
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Read input file.
    
    # Make sure an input file was given.
    if (fNameIn == ''):
        print('\nERROR: no input file given.  Use:  -i filename')
    
    # Open the CIF file and read the data.
    data = read_cif(fNameIn)
    
    # Extract lengths and angles from the CIF file.
    La = float(data['_cell_length_a'])
    Lb = float(data['_cell_length_b'])
    Lc = float(data['_cell_length_c'])
    alpha = math.radians(float(data['_cell_angle_alpha']))
    beta = math.radians(float(data['_cell_angle_beta']))
    gamma = math.radians(float(data['_cell_angle_gamma']))
    volume = float(data['_cell_volume'])
    
    # Extract the symmetry operations.  This will be a list of strings
    # such as:
    #    ['x,y,z', 'y,x,2/3-z', '-y,x-y,2/3+z', '-x,-x+y,1/3-z', ... ]
    ops = data['_symmetry_equiv_pos_as_xyz']
    
    # For proper evaluation, we need to convert "2/3" into "2./3", etc.
    # to prevent integer division which would turn e.g. 2/3 into 0.
    for i in range(len(ops)):
        ops[i] = ops[i].replace("0/", "0./") # also for e.g. 10/9
        ops[i] = ops[i].replace("1/", "1./")
        ops[i] = ops[i].replace("2/", "2./")
        ops[i] = ops[i].replace("3/", "3./")
        ops[i] = ops[i].replace("4/", "4./")
        ops[i] = ops[i].replace("5/", "5./")
        ops[i] = ops[i].replace("6/", "6./")
        ops[i] = ops[i].replace("7/", "7./")
        ops[i] = ops[i].replace("8/", "8./")
        ops[i] = ops[i].replace("9/", "9./")
    #    ops[i] = ops[i].replace("/", "./")
    
    # Get the atom labels and coordinates.
    labels = data['_atom_site_label']
    fX = [ float(s) for s in data['_atom_site_fract_x'] ]
    fY = [ float(s) for s in data['_atom_site_fract_y'] ]
    fZ = [ float(s) for s in data['_atom_site_fract_z'] ]
    
    # Create a list of 4-tuples, where each tuple is an atom:
    #   [ ('Si', 0.4697, 0.0, 0.0),  ('O', 0.4135, 0.2669, 0.1191),  ... ]
    atoms = [ (labels[i], fX[i], fY[i], fZ[i]) for i in range(len(labels)) ]
    
    # Make sure that all atoms lie within the unit cell.  Also convert 
    # names such as 'Oa1' into 'O'.
    for i in range(len(atoms)):
        (name,xn,yn,zn) = atoms[i]
        xn = (xn + 10.0) % 1.0
        yn = (yn + 10.0) % 1.0
        zn = (zn + 10.0) % 1.0
        name = extract_element(name)
        atoms[i] = (name,xn,yn,zn)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Use symmetry operations to create the unit cell.
    
    # The CIF file consists of a few atom positions plus several "symmetry
    # operations" that indicate the other atom positions within the unit cell.  So
    # using these operations, create copies of the atoms until no new copies can be
    # made.
    
    # Two atoms are on top of each other if they are less than "eps" away.
    eps = 0.01  # in Angstrom
    
    # For each atom, apply each symmetry operation to create a new atom.
    imax = len(atoms)
    i=0
    while (i < imax):
    
        label,x,y,z = atoms[i]
    
        for op in ops:
    
            # Python is awesome: calling e.g. eval('x,y,1./2+z') will convert the
            # string into a 3-tuple using the current values for x,y,z!
            xn,yn,zn = eval(op)
    
            # Make sure that the new atom lies within the unit cell.
            xn = (xn + 10.0) % 1.0
            yn = (yn + 10.0) % 1.0
            zn = (zn + 10.0) % 1.0
    
            # Check if the new position is actually new, or the same as a previous
            # atom.
            new_atom = True
            for at in atoms:
                if (abs(at[1]-xn) < eps  and  abs(at[2]-yn) < eps  and  abs(at[3]-zn) < eps):
                    new_atom = False
    
                    # Check that this is the same atom type.
                    if (at[0] != label):
                        print('\nERROR: invalid CIF file: atom of type %s overlaps with atom of type %s' % (at[0],label))
    
            # If the atom is new, add it to the list!
            if (new_atom):
                atoms.append( (label,xn,yn,zn) )  # add a 4-tuple
    
        # Update the loop iterator.
        i = i + 1
        imax = len(atoms)
    
    # Sort the atoms according to type alphabetically.
    atoms = sorted(atoms, key=lambda at: at[0])
    atoms.reverse()
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Create a larger box made of several unit cells: the super cell.
    
    atomlist = []
    
    for atom in atoms:
    
        # Get label and fractional coordinates.
        label,xf,yf,zf = atom
    
        for i in range(Nx):
                x = i+xf
    
                for j in range(Ny):
                    y = j+yf
    
                    for k in range(Nz):
                        z = k+zf
                        atomlist.append( (label,x,y,z) ) # add 4-tuple
    
    atoms = atomlist
    
    # If the user wants us to create a copy of the current CIF file, with
    # additional atoms, then do that.  Note that the atoms here have *fractional*
    # coordinates!
    if (fNameOut.endswith('.cif')):
        
        write_cif(fNameIn, atoms, fNameOut)
    
        print('\nERROR: Done writing extended CIF file (%d atoms in total).' % len(atoms))
        exit(0)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Convert the fractional coordinates into real coordinates.
    # The primitive vectors a,b,c are such that 
    #
    #   cos(alpha) = b.c / |b||c|
    #   cos(beta)  = a.c / |a||c|
    #   cos(gamma) = a.b / |a||b|
    #
    # with the convention
    #
    #   a = La*xhat
    #   b = bx*xhat + by*yhat
    #   c = cx*xhat + cy*yhat + cz*zhat
    #
    cosa = math.cos(alpha)
    #sina = math.sin(alpha)
    cosb = math.cos(beta)
    #sinb = math.sin(beta)
    cosg = math.cos(gamma)
    sing = math.sin(gamma)
    
    cosa2 = cosa * cosa
    cosb2 = cosb * cosb
    sing2 = sing * sing
    
    ax = La
    
    bx = Lb * cosg
    by = Lb * sing
    
    cx = Lc * cosb
    cy = Lc * (cosa - cosg*cosb) / sing
    cz = Lc * math.sqrt( 1 - (cosa2 + cosb2 - 2*cosg*cosb*cosa) / sing2 )
    
    # Use the volume to check if we did the vectors right.
    V = ax*by*cz
    
    if ( abs(V - volume) > 0.1):
        print("{}".format("~"*(shutil.get_terminal_size().columns)))
        print('WARNING: Volume of the unit cell declared in CIF ({:.2f} A^3) is different than the calculated from primitive vectors ({:.2f} A^3).\n'.format(volume, V))
        print("{}".format("~"*(shutil.get_terminal_size().columns)))
    
    # Check if we have a rectangular box.
    if (bx < eps  and  cx < eps  and cy < eps):
        make_rect_box = True
    
    # Determine the box size.
    Lx = Nx * La
    Ly = Ny * Lb
    Lz = Nz * Lc
    
    for i in range(len(atoms)):
    
        # Get label and fractional coordinates.
        label,xf,yf,zf = atoms[i]
    
        xa = xf * ax  # contribution of a-vector to the x-coordinate of this atom
        #ya = 0       # a-vector has no y-component, so does not affect y of atom
        #za = 0       # a-vector has no z-component, so does not affect z of atom
        
        xb = yf * bx  # contribution of b-vector to the x-coordinate of this atom
        yb = yf * by  # contribution of b-vector to the y-coordinate of this atom
        #zb = 0       # b-vector has no z-component, so does not affect z of atom
    
        xc = zf * cx  # contribution of c-vector to the x-coordinate of this atom
        yc = zf * cy  # contribution of c-vector to the y-coordinate of this atom
        zc = zf * cz  # contribution of c-vector to the z-coordinate of this atom
    
        # Add all contributions.
        xn = xa + xb + xc
        yn = yb + yc
        zn = zc
    
        if (make_rect_box):
            xn = (xn + Lx) % Lx
            yn = (yn + Ly) % Ly
            zn = (zn + Lz) % Lz
    
        atoms[i] = (label, xn, yn, zn)
    
    # Determine the box-vector.
    if (make_rect_box):
        box = (Lx, Ly, Lz)
    
    else:
        box = (Lx, Ly, Lz, Nx*cx, Ny*cy, Nz*cz)
    
    try:
        fOut = open(fNameOut, 'w')
    
        if (fNameOut.endswith('.xyz')):
            write_xyz(atoms, box, fOut)

    except TypeError:
        print("\nERROR: A type error occured when opening or writting XYZ file '{0}'".format(fNameOut))
        sys.exit()

    except NameError:
        print("\nERROR: Failed to open XYZ file.")
        sys.exit()

    except:
        print('\nERROR: Failed to write to XYZ output file')

    fOut.close()
# ======================================================================


# ======================================================================
def center_supercell(cif_output, verbose=0):
    """Takes the supercell file produced by Read_CIF and computes the
    center of the supercell coordinates to translate the supercell to
    the origin.

    Arguments:
    <str> cif_output
        Name of the file with the cartesian coordinates of the
        supercell.
    <int> verbose
        Adjusts the level of detail of the printouts.

    Returns:
    <numpy.ndarray> scell_geom_max_coords
        Array of 3 numbers with the maximum value of the coordinates x,
        y, and z of the centered supercell.
    <numpy.ndarray> scell_geom
        3-columns array with the x, y, and z coordinates of the 
        centered supercell.
    <numpy.ndarray> scell_elem
        1-column array with element symbols of the centered supercell.
    """
    
    # Find out if there are duplicate coordinates in the XYZ of the
    # supercell. Generate the indexes of unique atoms to exclude
    # superimposed atoms.
    scell_dupl = np.loadtxt(cif_output, skiprows=2, usecols=(0, 1, 2, 3), dtype="str")
    dirty_cell = [tuple(row) for row in scell_dupl]
    clean_cell, unique_idx = np.unique(dirty_cell, return_index=True, axis=0)

    # Compare the size of the clean (all unique coordinates) supercell
    # with that of the dirty supercell. If different, need to clean up
    # the supercell to avoid superimposed atoms.
    if len(scell_dupl) != len(clean_cell):
        print("{}".format("~"*(shutil.get_terminal_size().columns)))
        print("WARNING: The supercell contains {} duplicate coordinates. Skipping duplicates.".format(len(scell_dupl) - len(clean_cell)))
    
        # Creates two NumPy arrays: one with the coordinates of atoms in the
        # supercell and other with the element symbols of the atoms in it.
        scell_geom_dupl = np.loadtxt(cif_output, skiprows=2, usecols=(1, 2, 3), dtype=np.float64)
        scell_elem_dupl = np.loadtxt(cif_output, skiprows=2, usecols=(0), dtype="str")

        # Clean up duplicate coordiantes in the new NumPy arrays.
        scell_geom = scell_geom_dupl[unique_idx]
        scell_elem = scell_elem_dupl[unique_idx]

        print("         The number of unique coordinates in the supercell is now: {}\n".format(len(clean_cell)))
        print("{}".format("~"*(shutil.get_terminal_size().columns)))

    else:
        # Creates two NumPy arrays: one with the coordinates of atoms in the
        # supercell and other with the element symbols of the atoms in it.
        scell_geom = np.loadtxt(cif_output, skiprows=2, usecols=(1, 2, 3), dtype=np.float64)
        scell_elem = np.loadtxt(cif_output, skiprows=2, usecols=(0), dtype="str")

    if verbose >= 2:
        print("Generating monomers for all complete molecules in the supercell:")

    # Distances will be handled in Bohr.
    scell_geom = scell_geom / qcel.constants.bohr2angstroms

    # Calculation of the supercell center as the midpoint of all
    # coordinates.
    scell_cntr = (np.max(scell_geom, axis=0) - np.min(scell_geom, axis=0))/2.0
    
    if verbose >= 2:
        print("\nCurrently, the center of the supercell is located at:")
        print("x = %10.5f" % (scell_cntr[0]))
        print("y = %10.5f" % (scell_cntr[1]))
        print("z = %10.5f" % (scell_cntr[2]))
        print("\nThe supercell coordinates will be translated and centered on the origin.")

    # Translate the supercell to the origin.
    scell_geom -= scell_cntr
    
    # Return a NumPy array with 3 numbers: the maximum of each x, y, and
    # z.
    scell_geom_max_coords = np.max(scell_geom, axis=0)

    return scell_geom_max_coords, scell_geom, scell_elem
# ======================================================================


# ======================================================================
def supercell2monomers(cif_output, r_cut_monomer, bfs_thresh, verbose=1):
    """Takes the supercell cartesian coordinates file produced by
    Read_CIF, and passes it to the `center_supercell()` function which
    translates the supercell to the origin.

    The centered supercell geometries and elements arrays are passed to
    the Breadth-First Search which returns all fragments found in the
    supercell.

    This function also takes a cutoff distance, measured from the
    origin, which is then used to decide if a monomer should be included
    or not based on its proximity to the origin.

    Returns a dictionary with all the fragments (monomers) that are
    located within the cutoff region.
    
    Arguments:
    <str> cif_output
        Name of the file with the cartesian coordinates of the
        supercell.
    <float> r_cut_monomer
        Cutoff value to include fragments found in the supercell in the
        dictionary of N-mers.
    <int> verbose
        Adjusts the level of detail of the printouts.

    Returns:
    <dict> nmers
        A dictionary populated with N-mers (monomers at this time) and
        their corresponding atributes.
    """

    # Centering the supercell.
    scell_geom_max_coords, scell_geom, scell_elem = center_supercell(cif_output, verbose)

    # Check if each of dimensions of the supercell satisfies the 
    # condition of being the twice as long as the cutoff. This helps
    # filtering out fragments that contain incomplete molecules located
    # at the edges of the supercell.
    if (r_cut_monomer / qcel.constants.bohr2angstroms) > np.min(scell_geom_max_coords):
        print("\nERROR: Cutoff (%3.2f A) longer than half the smallest dimension of the supercell (%3.2f A)." \
              % (r_cut_monomer, np.min(scell_geom_max_coords)*qcel.constants.bohr2angstroms))
        print("       Please increase the dimensions of the supercell to at least twice r_cut_monomer or reduce the lenght of the cutoff.\n")
        sys.exit()

    # Start the BFS timer.
    bfs_start_time = time.time()
    
    # Passes the supercell geometry and elements to the breadth-first
    # search algorithm of QCDB to obtain fragments.
    #fragments = BFS(scell_geom, scell_elem, None, bfs_thresh)
    fragments = psi4.driver.qcdb.bfs.BFS(scell_geom, scell_elem, None, bfs_thresh)

    # Stop the BFS timer.
    bfs_stop_time = time.time() - bfs_start_time
    
    # Two lists containing geometries and elements of a fragment.
    frag_geoms = [scell_geom[fr] for fr in fragments]
    frag_elems = [scell_elem[fr] for fr in fragments]

    # List containing the magnitude of the shortest position vector of
    # each fragment.
    frag_r_mins = []

    for frag in frag_geoms:
        # Get the minium of the norm of the position vectors of all
        # atoms in the fragment.
        r_min_atom = np.min(np.linalg.norm(frag, axis=1))
        frag_r_mins.append(r_min_atom)
    
    # Convert the list of minimums to a NumPy array
    frag_r_mins = np.array(frag_r_mins)

    # Array mapping the indices of the fragments sorted by r_min_atom to
    # the order in the frag_r_mins list.
    mapper = frag_r_mins.argsort()

    # A dictionary for storing N-mers will be created. 

    # Then, for each fragment in the frag_geoms list, an index will be 
    # defined based on the mapper list.
    
    # If the shortest position vector is shorter that the monomers 
    # cutoff, a new dictionary will be created to store all the 
    # information corresponding to such N-mer.
    nmers = {}

    for i in range(len(frag_geoms)):
        index = mapper[i]

        # Discard edges of the supercell and keeps the monomers within a
        # sphere of the cutoff radius.
        if frag_r_mins[index] <= (r_cut_monomer / qcel.constants.bohr2angstroms):
            name = "1mer-" + str(i)
            
            # The dictionary of one N-mer.
            nmers[name] = {}
            
            nmers[name]["monomers"] = [i]
            nmers[name]["elem"] = frag_elems[index]
            nmers[name]["coords"] = frag_geoms[index]
            nmers[name]["rmin"] = frag_r_mins[index]
            nmers[name]["delimiters"] = []
            nmers[name]["com"] = center_of_mass(frag_elems[index], frag_geoms[index])
            nmers[name]["priority_min"] = 0.0
            nmers[name]["priority_com"] = 0.0
            nmers[name]["priority_cutoff"] = 0.0

    total_number_of_monomers = len(nmers.keys())

    if verbose >= 2:
        print("\nThe BFS algorithm found {} monomers in the supercell in {:.2f} s".format(total_number_of_monomers, bfs_stop_time))
    
    return nmers
# ======================================================================


# ======================================================================
def create_nmer(nmers, ref_monomer, other_monomers, verbose=1):
    """Takes a `nmers` dictionary, and two strings with the keys of a
    refrence monomer and other monomer.

    This function will merge both monomers and create a new N-mer.

    It returns a string `nm_new_name` and a `nm_new` dictionary
    containing the information of the newly generated N-mer.

    Arguments:

    Returns:
    <str> nm_new_name
        Name of the new N-mer, typically `n`mer-`i`+`j`+`k`+... where
        `n` is the order of the N-mer and `i`, `j`, ... the indices of
        the monomer that compose it.
    <dict> nm_new
        Dictionary of the new N-mer containing its corresponding
        information.
    """
    
    nm_new = {}

    # Monomers included in the new N-mer.
    nm_new_monomers = []
    nm_new_monomers.append(ref_monomer["monomers"][0])
    nm_new_monomers.extend(other_monomers)
    nm_new["monomers"] = nm_new_monomers
    
    # Elements and coordinates of the atoms in the new N-mer.
    nm_new_elem_arrays = []
    nm_new_coords_arrays = []
    nm_new["atoms_per_monomer"] = []

    for monomer in nm_new_monomers:
        name = "1mer-" + str(monomer)
        nm_new["atoms_per_monomer"].append(len(nmers[name]["elem"]))
        nm_new_elem_arrays.append(nmers[name]["elem"])
        nm_new_coords_arrays.append(nmers[name]["coords"])

    nm_new["elem"] = np.concatenate(nm_new_elem_arrays, axis=0)
    nm_new["coords"] = np.concatenate(nm_new_coords_arrays, axis=0)
    
    # Indices of each monomer in the array of elements and coordinates 
    # (For spliting N-mers into monomers).
    nm_new["delimiters"] = np.cumsum(nm_new["atoms_per_monomer"])
    
    # Norm of the shortest position vector of an atom in the new N-mer.
    # The reference is always closest to the origin so it is always the
    # minimum rmin
    nm_new["rmin"] = ref_monomer["rmin"]

    # Nuclear repulsion energy of the new N-mer.
    nm_new["nre"] = nre(nm_new["elem"], nm_new["coords"])

    # Chemical space eigenvalues of the new N-mer.
    nm_new["chsev"] = chemical_space(nm_new["elem"], nm_new["coords"])

    # Non additive N-body energy of the N-mer.
    nm_new["nambe"] = 0.0

    # Counter of the number of replicas of this N-mer.
    nm_new["replicas"] = 1

    # Contribution of unique structure to crystal lattice energy.
    nm_new["contrib"] = 0.0

    # Shortest separation vector between the atoms of the monomers in
    # the N-mer.
    nm_new["min_monomer_separations"] = []

    # Separation vector between center of mass of the monomers in the
    # N-mer.
    nm_new["com_monomer_separations"] = []

    for a in nm_new_monomers:
        a_name = "1mer-" + str(a)
        
        for b in nm_new_monomers:
            
            if b <= a:
                continue
            
            b_name = "1mer-" + str(b)

            distm, r_min = distance_matrix(nmers[a_name]["coords"], nmers[b_name]["coords"]) 
            nm_new["min_monomer_separations"].append(r_min)
            nm_new["com_monomer_separations"].append(np.linalg.norm(nmers[a_name]["com"] - nmers[b_name]["com"]))
    
    # Distance to compare against cutoff values.
    nm_new["max_mon_sep"] = max(nm_new["min_monomer_separations"])
    nm_new["max_com_sep"] = max(nm_new["com_monomer_separations"])

    #TODO: Create a keyword to select which criterion to use.

    # Criterion to launch energy calculations.
    nm_new["priority_min"] = 0.0

    priority_min = 1.0
    
    for rmin in nm_new["min_monomer_separations"]:
        one_over_rmin3 = 1.0/rmin**3
        priority_min *= one_over_rmin3
    
    nm_new["priority_min"] = priority_min

    # Alternative criterion to launch energy calculations.
    nm_new["priority_com"] = 0.0

    priority_com = 1.0
    
    for rcom in nm_new["com_monomer_separations"]:
        one_over_rcom3 = 1.0/rcom**3
        priority_com *= one_over_rcom3
    
    nm_new["priority_com"] = priority_com

    # Another alternative criterion to launc energy calculations based
    # on the maximum separation between monomers in an N-mer.
    nm_new["priority_cutoff"] = 0.0

    max_sep = max(nm_new["min_monomer_separations"])
    main_contrib = 1.0/(max_sep**(len(nm_new_monomers)**2.0))

    idx = 0
    add = 0.0
    for rmin in nm_new["min_monomer_separations"]:

        if idx != 0:
            one_over_rmin3 = 1.0e-10/(rmin**(len(nm_new_monomers)**2.0))
            add += one_over_rmin3

        idx += 1

    priority_cutoff = main_contrib + add

    nm_new["priority_cutoff"] = priority_cutoff

    # Key of the new N-mer in the nmers dictionary.
    nm_new_name = str(len(nm_new_monomers)) + "mer-" + "+".join(map(str, nm_new_monomers))

    return nm_new_name, nm_new
# ======================================================================


# ======================================================================
def center_of_mass(elems, geoms):
    """Takes the element symbols and coordinates of a set of atoms and
    computes the center of mass of the molecule.
    
    Arguments:
    <numpy.ndarray> elems
        Array of 1 column with the atomic symbols of each atom.
    <numpy.ndarray> geoms
        Array of 3 columns with the coordinates of the system.

    Returns:
    <numpy.ndarray> com
        Array of 3 numbers with the center of mass of the system.
    """
    
    com = np.array([0.0, 0.0, 0.0])
    total_mass = 0.0

    for at in range(len(geoms)):
        m = qcel.periodictable.to_mass(elems[at])
        com += m*geoms[at]
        total_mass += m
    
    com /= total_mass
    
    return com
# ======================================================================


# ======================================================================
def distance_matrix(a, b):
    """Euclidean distance matrix between rows of arrays `a` and `b`.
    Equivalent to `scipy.spatial.distance.cdist(a, b, 'euclidean')`.
    Returns a.shape[0] x b.shape[0] array. It also returns the shortest
    distance between the atoms of `a` and of `b`.
    """

    assert(a.shape[1] == b.shape[1])
    distm = np.zeros([a.shape[0], b.shape[0]])
    
    for i in range(a.shape[0]):
    
        for j in range(b.shape[0]):
            distm[i, j] = np.linalg.norm(a[i] - b[j])

    r_min = np.min(distm)

    return distm, r_min
# ======================================================================


# ======================================================================
def nre(elem, geom):
    """Takes two Numpy arrays, one with the element symbols and one with
    coordinates of a set of atoms and returns a float number with the 
    computed nuclear repulsion energy.
    """
    
    nre = 0.
    for at1 in range(geom.shape[0]):

        for at2 in range(at1):
            dist = np.linalg.norm(geom[at1] - geom[at2])
            nre += qcel.periodictable.to_Z(elem[at1]) * qcel.periodictable.to_Z(elem[at2]) / dist

    return nre
# ======================================================================


# ======================================================================
def chemical_space(elem, geom):
    """Takes the element symbols and coordinates of a set of atoms and
    computes the chemical space matrix, eigenvalues, and eigenvectors
    of the chemical system and returns a list with the sorted
    eigenvalues.

    Arguments:
    <numpy.ndarray> elems
        Array of 1 column with the atomic symbols of each atom.
    <numpy.ndarray> geoms
        Array of 3 columns with the coordinates of the system.

    Returns:
    <numpy.ndarray> chem_spc_eigen_values
        Array of x numbers with the sorted eigenvalues of the molecular
        system, where x is the number of atoms in the system.
    """

    # Get the number of atoms in the N-mer
    natoms = geom.shape[0]
    
    # Create a NumPy matrix
    M = np.zeros((natoms,natoms))

    # Iterate over atoms
    for i in range(natoms):
   
        # Get the charge of atom i
        charge_i = qcel.periodictable.to_Z(elem[i])
        
        # Fill the diagonal with the special polynomial from of:
        # DOI: 10.1103/PhysRevLett.108.058301     
        M[i,i] = 0.5 * np.power(charge_i, 2.4)
    
        for j in range(i):
            
            # Get the charge of atom j
            charge_j = qcel.periodictable.to_Z(elem[j])
            
            # Compute distance between i and j
            dist = np.linalg.norm(geom[i] - geom[j])
            
            # Compute Coulomb interaction between i and j
            ij_elem = charge_i * charge_j / dist

            M[i,j] = ij_elem
            # Symmetric Matrix 
            M[j,i] = ij_elem

    # Solve the eigenvalue problem
    eigenvalues, eigenvectors = np.linalg.eig(M)

    # Eigenvalues must be in list to use 'sort'
    sorted_eigenvalues = list(eigenvalues)
    
    # Sort the eigenvalues in order of "decreasing absolute value".
    # This first sort is done to guarantee the same sort in the case 
    # that two eigenvalues are the same magnitude but different sign.
    sorted_eigenvalues.sort(key = lambda x: -x)
    sorted_eigenvalues.sort(key = lambda x: -abs(x))
    
    # Cast back to a NumPy array.
    chem_spc_eigen_values = np.array(sorted_eigenvalues)

    return chem_spc_eigen_values
# ======================================================================


# ======================================================================
def build_nmer(nmers, total_monomers, nmer_type, nmer_separation_cutoff, coms_separation_cutoff, uniq_filter, verbose=1):
    """Takes the nmers dictionary and the filters criteria and builds a
    dictionary with generated N-mers of a given order, i.e. dimers. Such
    N-mers are generated after filtering out all other combinations that
    do not fit the filters criteria.

    Arguments:
    <dict> nmers
        Dictionary containing dictionaries for each N-mer in the system.
    <int> total_monomers
        Number of total monomers in the supercell, as determined by the
        BFS algorithm.
    <str> nmer_type
        The order of the new N-mers to be created in string format.
    <float> nmer_separation_cutoff
        Maximum allowed separation between closest atoms of different
        monomers in the new N-mers.
    <float> coms_separation_cutoff
        Maximum allowed separation between closest COMs of different 
        monomers in the new N-mers.
    <string> uniq_filter
        Type of filter to be used to determine chemical equivalency.
        Either the chemical-space or the dreamaligner filter.
    <int> verbose
        Adjusts the level of detail of the printouts.

    Returns:
    <dict> nmers
        Dictionary containing dictionaries for each N-mer in the system.
        Updated to include any newly generated N-mers.
    """


    # Function reused for different types of N-mers.
    build_nmers_start_time = time.time()

    if nmer_type == "dimers":
        num_monomers = 2

    elif nmer_type == "trimers":
        num_monomers = 3

    elif nmer_type == "tetramers":
        num_monomers = 4

    elif nmer_type == "pentamers":
        num_monomers = 5
    
    else:
        print("\nERROR: The N-mer type must be defined as 'dimers', 'trimers', 'tetramers', or 'pentamers'.")

    if verbose >= 2:
        print("")

    counter_new_nmers = 0 # Number of new N-mers generated
    counter_dscrd_sep = 0 # N-mers filtered out by atomic separation
    counter_dscrd_com = 0 # N-mers filtered out by COM separation
    counter_dscrd_rep = 0 # N-mers filtered out as a replicas

    new_nmers = {}

    # TODO: Support for crystals with more than one molecule in the
    #       primitive unit cell.
    num_ref_monomers = 1

    # Take the first monomer as the reference, aka 1mer-0.
    for ref_monomer_idx in range(num_ref_monomers):
        monomer_key = "1mer-" + str(ref_monomer_idx)
        ref_monomer = nmers[monomer_key]

        # Iterate over combinations that include all other monomers in unrepeated fashion.
        for other_monomers in itertools.combinations(range(ref_monomer_idx+1, total_monomers), num_monomers-1):

            # Start N-mer building timer.
            nmer_start_time = time.time()

            # Combine reference monomer with other monomers to create a new N-mer.
            new_nmer_name, new_nmer = create_nmer(nmers, ref_monomer, other_monomers, verbose)

            # Start applying N-mer filters to weed out unnecessary calculations.
            # First, use the corresponding N-mer cutoff (dimer cutoff, trimer cutoff, ...)
            if new_nmer["max_mon_sep"] > (nmer_separation_cutoff / qcel.constants.bohr2angstroms):
                
                if verbose >= 2:
                    
                    # Stop N-mer building timer.
                    nmer_stop_time = time.time() - nmer_start_time

                    print(
                        "{} discarded in {:.2f} s: Maximum separation between closest atoms of different monomers is {:3.2f} A, longer than cutoff {:3.2f} A.".format(
                            new_nmer_name, nmer_stop_time, new_nmer["max_mon_sep"]*qcel.constants.bohr2angstroms, nmer_separation_cutoff))
                        
                    counter_dscrd_sep += 1
            
            # Second, use the global COM cutoff. Same distance applies to all N-mer orders.
            elif new_nmer["max_com_sep"] > (coms_separation_cutoff / qcel.constants.bohr2angstroms):
                
                if verbose >= 2:
                    
                    # Stop N-mer building timer.
                    nmer_stop_time = time.time() - nmer_start_time
                    print(
                        "{} discarded in {:.2f} s: Maximum separation between closest COMs of different monomers is {:3.2f} A, longer than cutoff {:3.2f} A.".format(
                            new_nmer_name, nmer_stop_time, new_nmer["max_com_sep"]*qcel.constants.bohr2angstroms, coms_separation_cutoff))

                    counter_dscrd_com += 1

            # Third, employ more sophisticated filters: NRE and ChSEV or Dreamaligner.
            else:
                found_duplicate = False
                
                #if False:
                if uniq_filter == "Dreamaligner":
                    chemical_space = False
                    dreamliner = True
                
                else:
                    chemical_space = True
                    dreamliner = False

                nre_filter_ran = False
                chsev_filter_ran = False
                rmsd_filter_ran = False

                for kexisting, existing in new_nmers.items():

                    # Fourth, the nuclear repulsion energy filter.
                    nre_diff = abs(existing["nre"] - new_nmer["nre"])
                    nre_filter_ran = True

                    # If NRE difference is large, this is a new N-mer.
                    # Threfore reset posterior filters ran flags, there
                    # is no need to run further filters than NRE.
                    if nre_diff > 1.e-5:

                        chsev_filter_ran = False
                        rmsd_filter_ran = False                        
                    
                    # Fifth, if the NRE is small, proceed to next filters.
                    else:
                        # Chemical space eigenvalues filter.
                        if chemical_space == True:

                            chsev_diff = np.linalg.norm(existing["chsev"] - new_nmer["chsev"])
                            chsev_filter_ran = True

                            if chsev_diff < 1.e-3:
                                found_duplicate = True

                                if verbose >= 2:

                                    # Stop N-mer building timer.
                                    nmer_stop_time = time.time() - nmer_start_time

                                    print("{} discarded in {:.2f} s: Replica of {}. NRE difference is {:.1e} and ChSEV difference is {:.1e}.".format(
                                        new_nmer_name, nmer_stop_time, kexisting, nre_diff, chsev_diff))

                                existing["replicas"] += 1
                                counter_dscrd_rep += 1

                                break

                        # RMSD filter.
                        if dreamliner == True:
                        
                            # Block B787 printout
                            sys.stdout = open(os.devnull, 'w') #NOTE: Check if this is still needed (?)
                            
                            # Call the dreamliner from QCDB.
                            rmsd, mill = qcel.molutil.B787(rgeom=existing["coords"], cgeom=new_nmer["coords"], runiq=existing["elem"], cuniq=new_nmer["elem"], run_mirror=True, verbose=2)
                            rmsd_filter_ran = True

                            # Reanable printout
                            sys.stdout = sys.__stdout__ #NOTE: Check if this is still needed (?)

                            if rmsd < 1.e-3:
                                found_duplicate = True

                                if verbose >= 2:
                                    
                                    # Stop N-mer building timer.
                                    nmer_stop_time = time.time() - nmer_start_time

                                    print("{} discarded in {:.2f} s: Replica of {}. NRE difference is {:.1e} and RMSD is {:.1e}.".format(
                                        new_nmer_name, nmer_stop_time, kexisting, nre_diff, rmsd))

                                existing["replicas"] += 1
                                counter_dscrd_rep += 1

                                break

                # If the structure was not filtered out, then store the new N-mer.
                if not found_duplicate:
                    new_nmers[new_nmer_name] = new_nmer

                    if verbose >= 2:

                        # Stop N-mer building timer.
                        nmer_stop_time = time.time() - nmer_start_time

                        if nre_filter_ran == True:
                            
                            if chsev_filter_ran == True:
                                print("{} generated in {:.2f} s: New N-mer NRE difference is {:.1e}, lowest ChSEV difference found is {:.1e}.".format(
                                    new_nmer_name, nmer_stop_time, nre_diff, chsev_diff))
                        
                            if rmsd_filter_ran == True:
                                print("{} generated in {:.2f} s: New N-mer NRE difference is {:.1e}, lowest RMSD found is {:.1e}.".format(
                                    new_nmer_name, nmer_stop_time, nre_diff, rmsd))

                            if chsev_filter_ran == False and rmsd_filter_ran == False:
                                print("{} generated in {:.2f} s: New N-mer NRE difference is {:.1e}.".format(new_nmer_name, nmer_stop_time, nre_diff))

                        else:
                            print("{} generated in {:.2f} s: New N-mer NRE is {:.12f}.".format(new_nmer_name, nmer_stop_time, new_nmer["nre"]))

                    counter_new_nmers += 1
    
    if verbose >= 2:
        build_nmers_stop_time = time.time() - build_nmers_start_time
        print("\n{} unique {} were found and generated in {:.1f} s.".format(counter_new_nmers, nmer_type, build_nmers_stop_time))

    if verbose >= 2:
        print("\n{} {} did not meet the atomic separation cutoff and were discarded.".format(counter_dscrd_sep, nmer_type))
        print("{} {} did not meet the center of mass separation cutoff and were discarded.".format(counter_dscrd_com, nmer_type))
        print("{} {} were duplicates of another dimer and were discarded.".format(counter_dscrd_rep, nmer_type))
    
    nmers.update(new_nmers)

    return nmers
# ======================================================================


# ======================================================================
def nmer2psiapimol(nmers, keynmer, nmer, verbose=0):
    """Takes the `nmers` dictionary; `keynmer`, the key of a given N-mer of
    the N-mers dictionary; and `nmer`, its corresponding dictionary.
    Returns a string `psi_api_molecule` that defines a molecule in PSI4
    API mode. This string starts with three lines specifying the use of
    atomic units, avoidance of center of mass translation, and
    avoidance of reorientation. The next lines contain the element
    symbol and coordinates, for each atom in the N-mer, separating each
    monomer with lines containing a double hyphen.
    """

    psi_api_molecule = "\nunits = au\n"
    psi_api_molecule += "no_com\n"
    psi_api_molecule += "no_reorient\n"

    for at in range(nmer["coords"].shape[0]):

        if at in nmer["delimiters"]:
            psi_api_molecule += "--\n"

        psi_api_molecule += "  {:6} {:16.8f} {:16.8f} {:16.8f} \n".format(nmer["elem"][at], nmer["coords"][at][0], nmer["coords"][at][1], nmer["coords"][at][2])

    return psi_api_molecule
# ======================================================================


# ======================================================================
def monomer2makefp(cif_output, nmer, verbose=0):
    """.
    """

    makefp_input =  "! CrystaLattE: GAMESS input for generation of EFP parameters\n"
    makefp_input += " $contrl units=angs local=boys runtyp=makefp\n"
    makefp_input += "   mult=1 icharg=0 coord=cart icut=11 $end\n"
    makefp_input += " $system timlim=99999 mwords=200 $end\n"
    makefp_input += " $scf dirscf=.t. soscf=.f. diis=.t. conv=1.0d-06 $end\n"
    makefp_input += " $basis gbasis=n311 ngauss=6 npfunc=2 ndfunc=3 nffunc=1\n"
    makefp_input += "   diffs=.t. diffsp=.t. $end\n"
    makefp_input += " $makefp pol=.t. disp=.t. exrep=.t. chtr=.f. $end\n"
    makefp_input += " $stone\n"
    makefp_input += "   bigexp=0.0\n"
    makefp_input += " $end\n"
    makefp_input += " $data\n"
    makefp_input += "  {}\n".format(cif_output.split(".")[0])
    makefp_input += "c1\n"

    for at in range(nmer["coords"].shape[0]):

        if at in nmer["delimiters"]:
            makefp_input += "--\n"
        
        makefp_input += " {:6} {:4.1f} {:16.8f} {:16.8f} {:16.8f} \n".format(nmer["elem"][at] + str(at + 1), qcel.periodictable.to_Z(nmer["elem"][at]), nmer["coords"][at][0], nmer["coords"][at][1], nmer["coords"][at][2])

    makefp_input += " $end"

    owd = os.getcwd()
    makefp_folder = cif_output[:-4]

    try:
        os.mkdir(makefp_folder)
    
    except FileExistsError:
        pass

    os.chdir(makefp_folder)
    makefp_filename = cif_output.split(".")[0] + "-makefp-lb.inp"

    with open(makefp_filename, "w") as makefp_f:
        
        for line in makefp_input:
            makefp_f.write(line)

    os.chdir(owd)
# ======================================================================


# ======================================================================
def nmer2psithon(cif_output, nmers, keynmer, nmer, rminseps, rcomseps, psi4_method, psi4_bsse, psi4_memory, verbose=0):
    """.
    """

    # This strings are used by the psithonyzer script. Any changes
    # applied to them here must be synchronized there The order in
    # which they are printed is also important, especially for the
    # different types of priorities.
    psithon_input =  "# PSI4 file produced by CrystaLattE\n\n"
    psithon_input += "# Generated from:               {}\n".format(cif_output)
    psithon_input += "# Psithon input for N-mer:      {}\n".format(keynmer)
    psithon_input += "# Number of atoms per monomer:  {}\n".format(nmer["atoms_per_monomer"])
    psithon_input += "# Number of replicas:           {}\n".format(nmer["replicas"])
    psithon_input += "# COM priority:                 {:12.12e}\n".format(nmer["priority_com"])
    psithon_input += "# Minimum COM separations:      {}\n".format(rcomseps.lstrip(" "))
    psithon_input += "# Separation priority:          {:12.12e}\n".format(nmer["priority_min"])
    psithon_input += "# Minimum monomer separations:  {}\n".format(rminseps.lstrip(" "))
    psithon_input += "# Cutoff priority:              {:12.12e}\n".format(nmer["priority_cutoff"])
    psithon_input += "# Nuclear repulsion energy:     {}\n".format(nmer["nre"])
    
    psithon_input += "\nmemory {}\n".format(psi4_memory)
    
    mymol = keynmer.replace("2mer", "Dimer").replace("3mer", "Trimer").replace("4mer", "Tetramer").replace("5mer", "Pentamer").replace("-", "_").replace("+", "_")
    psithon_input += "\nmolecule {} {{\n".format(mymol)
    
    for at in range(nmer["coords"].shape[0]):

        if at in nmer["delimiters"]:
            psithon_input += "--\n"
        
        psithon_input += "  {:6} {:16.8f} {:16.8f} {:16.8f} \n".format(nmer["elem"][at], nmer["coords"][at][0], nmer["coords"][at][1], nmer["coords"][at][2])

    psithon_input += "units = au\n"
    psithon_input += "no_com\n"
    psithon_input += "no_reorient\n"
    
    psithon_input += "}\n"

    psithon_input += "\nset {\n"
    psithon_input += "  e_convergence 10\n"
    psithon_input += "  d_convergence 10\n"
    psithon_input += "  scf_type df\n"
    psithon_input += "  mp2_type df\n" 
    psithon_input += "  cc_type df\n"
    psithon_input += "  freeze_core true\n"
    psithon_input += "}\n"

    # Hartree-Fock is called with the 'scf' string in Psithon mode.
    if psi4_method.lower().startswith("hf"):
        psithon_method = "scf" + psi4_method[2:]

    else:
        psithon_method = psi4_method

    psithon_input += "\nenergy('{}', bsse_type = '{}')\n".format(psithon_method, psi4_bsse)
    psithon_input += "\n"

    owd = os.getcwd()
    psithon_folder = cif_output[:-4]

    try:
        os.mkdir(psithon_folder)
    
    except FileExistsError:
        pass

    os.chdir(psithon_folder)
    psithon_filename = keynmer + ".in"

    with open(psithon_filename, "w") as psithon_f:
        
        for line in psithon_input:
            psithon_f.write(line)

    os.chdir(owd)
# ======================================================================

def gamessefp_mbe(fragment_potential, fragment_geometries, indices, keynmer, nmer):

    fragment_geometry = """ $efrag\n"""
    fragment_geometry += """coord=cart\n"""
    for i in indices:
        fragment_geometry += fragment_geometries[i]
    fragment_geometry += """$end"""

    efp_command = fragment_geometry.rstrip() + '\n ' + fragment_potential.lstrip()

    # now we'll make a dummy molecule to hold the efp command
    # the symbols and geometry don't matter at all 
    dummy_mol = qcel.models.Molecule(
        symbols=nmer['elem'],
        geometry=nmer['coords'],
        extras={"efp" : efp_command}
    )

    # evaluate the efp energy of the H2O dimer using our calculated EFP potential
    inp = qcel.models.AtomicInput(
        molecule=dummy_mol,
        driver="energy",
        model={"method": "efp"},
        keywords={},
    )
    res = qcng.compute(inp, program="gamess", raise_error=True, return_dict=True)
    assert res["success"] is True

    efp_ene = res["return_result"]
    #pprint(res["extras"]["qcvars"])

    #print(f'computed interaction energy is {627.509*efp_ene:0.5f} kcal/mol')
    return efp_ene

# ======================================================================
def gamessefp_nambe(fragment_potential, fragment_orientation, nmer, keynmer):

    tstart = time.time()

    ## the nmer for calculating the non-additive nbody energy
    #nmer = nmers[keynmer]

    # the numer of monomers in this nmer
    count = len(nmer['monomers'])

    coords = nmer['coords']
    elems = nmer['elem']
    delimiters = nmer['delimiters']
    atoms_per_monomer = nmer['atoms_per_monomer']
    monomers = nmer['monomers']

    fragment_geometries = []
    for monomer_ind in range(count):
        atom_start = delimiters[monomer_ind] - atoms_per_monomer[monomer_ind]
        atom_end = delimiters[monomer_ind]
        coords_monomer = coords[atom_start:atom_end]
        elems_monomer = elems[atom_start:atom_end]
        monomer_distances = distance_matrix(coords_monomer, coords_monomer)[0]

        # TODO: unhardcode orientation
        found_match = False
        for a0, a1, a2 in itertools.product(range(6), range(6), range(6)):
            temp_coords = np.array([coords_monomer[a0], coords_monomer[a1], coords_monomer[a2]])
            temp_orientation = distance_matrix(temp_coords, temp_coords)[0]
            match = np.max(np.abs(temp_orientation - fragment_orientation))
            if match < 0.00001:
                found_match = True
                break
        assert found_match

        fragment_geometry = """FRAGNAME=FRAGNAME\n"""
        fragment_geometry += f"""A01{elems_monomer[a0]:2s}   {coords_monomer[a0,0]}   {coords_monomer[a0,1]}   {coords_monomer[a0,2]}\n"""
        fragment_geometry += f"""A02{elems_monomer[a1]:2s}   {coords_monomer[a1,0]}   {coords_monomer[a1,1]}   {coords_monomer[a1,2]}\n"""
        fragment_geometry += f"""A03{elems_monomer[a2]:2s}   {coords_monomer[a2,0]}   {coords_monomer[a2,1]}   {coords_monomer[a2,2]}\n"""
        fragment_geometries.append(fragment_geometry)

    if count == 2:
        nambe = gamessefp_mbe(fragment_potential, fragment_geometries, [0, 1], keynmer, nmer)
    elif count == 3:
        e_abc = gamessefp_mbe(fragment_potential, fragment_geometries, [0, 1, 2], keynmer, nmer)
        e_ab = gamessefp_mbe(fragment_potential, fragment_geometries, [0, 1], keynmer, nmer)
        e_ac = gamessefp_mbe(fragment_potential, fragment_geometries, [0, 2], keynmer, nmer)
        e_bc = gamessefp_mbe(fragment_potential, fragment_geometries, [1, 2], keynmer, nmer)

        #key_ab, key_ac, key_bc = f'2mer-{monomers[0]}+{monomers[1]}', f'2mer-{monomers[0]}+{monomers[2]}', f'2mer-{monomers[1]}+{monomers[2]}'
        #if key_ab in nmers:
        #    e_ab = nmers[key_ab]['nambe']
        #else:
        #    e_ab = gamessefp_mbe(fragment_potential, fragment_geometries, [0, 1], keynmer, nmer)

        #if key_ac in nmers:
        #    e_ac = nmers[key_ac]['nambe']
        #else:
        #    e_ac = gamessefp_mbe(fragment_potential, fragment_geometries, [0, 2], keynmer, nmer)

        #if key_bc in nmers:
        #    e_bc = nmers[key_bc]['nambe']
        #else:
        #    e_bc = gamessefp_mbe(fragment_potential, fragment_geometries, [1, 2], keynmer, nmer)
        nambe = e_abc - e_ab - e_ac - e_bc
    else:
        exit()

    dt = time.time() - tstart
    print(f'{keynmer:20}   {dt:5.1f} (s)            {2625.5*nambe:9.5f} kJ/mol')

    return nambe

    

# ======================================================================


# ======================================================================
def nmer2libefpmbe(cif_output, nmers, keynmer, nmer, rminseps, rcomseps, verbose=0):
    """This function will write to disk a LibEFP input for the passed
    N-mer. To be able to compute the non-additive many-body energy, a
    trick has been implemented. N calculations will be performed for an
    N-mer, those will contain all possible dimers. At the moment, this
    should not be used for tetramers and pentamers.

    Arguments:
    <str> cif_output
        Name of the file with the cartesian coordinates of the 
        supercell.
    <dict> nmers
        Dictionary containing dictionaries for each N-mer in the system.
    <str> keynmer
        Name of the N-mer as it appears in the key of its corresponding
        entry in the nmers dictionary.
    <dict> nmer
        Dictionary containing all the information of one N-mer of the
        system.
    <str> rminseps
        List of separations between closest atoms of different monomers
        in the N-mer.
    <str> rcomseps
        List of separations between COMs of different monomers in the
        N-mer.
    <int> verbose
        Adjusts the level of detail of the printouts.

    Returns:
    None
    """

    # List to use pairwise energy substractions
    ligands = ["a","b","c","d","e"]

    # Name of the EFP fragments. Can be found in the .efp pothential
    # file preceeded by a $, i.e. $BENZENE_L or $water
    frgname = cif_output.split(".")[0] #TODO: Implement input keyword for fragment names.

    # Generate as many inputs as monomers exist in the N-mer.
    # This is an intent to circumvent the problem of computing the
    # substraction of (N-1)-body energies from N-mers.
    for i in range(len(nmer['monomers'])):

        libefpmbe_input =  "# LibEFP input file produced by CrystaLattE\n\n"
        
        libefpmbe_input += "# WARNING: This part of the CrystaLattE code is designed to work with\n"
        libefpmbe_input += "#          Prof. Lyudmila V. Slipchenko's fork of LibEFP 1.5.0 which can\n"
        libefpmbe_input += "#          be found at:\n"
        libefpmbe_input += "#          https://github.com/libefp2/libefp.git\n\n"

        libefpmbe_input += "# WARNING: This implementation should be limited to the computation of\n"
        libefpmbe_input += "#          dimers and trimers. This is due to the limitation introduced\n"
        libefpmbe_input += "#          by the usage of LibEFP, which at the moment is not able to\n"
        libefpmbe_input += "#          decompose the many-body energy into lower-order energies.\n\n"

        libefpmbe_input += "# WARNING: Name of the fragment are assigned based on the CIF name.\n"
        libefpmbe_input += "#          Check that the name of the fragments corresponds to the\n" 
        libefpmbe_input += "#          fragment name field in the .efp potential file.\n\n"
        
        # This strings are used by the analysis script. Any changes
        # applied to them here must be synchronized there The order in
        # which they are printed is also important, especially for the
        # different types of priorities.

        libefpmbe_input += "# Generated from:               {}\n".format(cif_output)
        libefpmbe_input += "# LibEFP input for N-mer:       {}\n".format(keynmer)
        libefpmbe_input += "# Number of atoms per monomer:  {}\n".format(nmer["atoms_per_monomer"])
        libefpmbe_input += "# Number of replicas:           {}\n".format(nmer["replicas"])
        libefpmbe_input += "# COM priority:                 {:12.12e}\n".format(nmer["priority_com"])
        libefpmbe_input += "# Minimum COM separations:      {}\n".format(rcomseps.lstrip(" "))
        libefpmbe_input += "# Separation priority:          {:12.12e}\n".format(nmer["priority_min"])
        libefpmbe_input += "# Minimum monomer separations:  {}\n".format(rminseps.lstrip(" "))
        libefpmbe_input += "# Cutoff priority:              {:12.12e}\n".format(nmer["priority_cutoff"])
        libefpmbe_input += "# Nuclear repulsion energy:     {}\n".format(nmer["nre"])

        libefpmbe_input += "\n"
        libefpmbe_input += " run_type sp\n"
        libefpmbe_input += " coord points\n"
        libefpmbe_input += " terms elec pol disp xr\n"
        libefpmbe_input += " elec_damp overlap\n" #TODO: Implement input keywords to support for other EFP damping methods.
        libefpmbe_input += " disp_damp overlap\n" #TODO: Implement input keywords to support for other EFP damping methods.
        libefpmbe_input += " pol_damp tt\n"       #TODO: Implement input keywords to support for other EFP damping methods.

        # To compute the non-addditive many-body energy using LibEFP
        # we are going to use the pairwise experimental feature of the
        # code. However, for dimers there is no need to use such 
        # feature. Ideally, in the future, this will be fix on LibEFP.
        if len(nmer["monomers"]) > 2:
            libefpmbe_input += " enable_pairwise true\n" #NOTE: Requires fork from https://github.com/libefp2/libefp.git
            libefpmbe_input += " ligand {}\n".format(i)  #NOTE: Requires fork from https://github.com/libefp2/libefp.git
        
        libefpmbe_input += " userlib_path .\n"

        libefpmbe_input += "\nfragment {}\n".format(frgname)
        
        # Because LibEFP only requests the first three atoms of a 
        # fragment, we skip writting the coordinates of atoms with
        # indexes greater than 2 for each fragment.
        line_idx = 0

        for at in range(nmer["coords"].shape[0]):
            
            if at in nmer["delimiters"]:
                libefpmbe_input += "\nfragment {}\n".format(frgname)
                line_idx = 0

            if line_idx < 3:
                libefpmbe_input += "{:16.8f} {:16.8f} {:16.8f} \n".format(nmer["coords"][at][0], nmer["coords"][at][1], nmer["coords"][at][2])
            
            line_idx += 1

        # Is this blank line needed?
        libefpmbe_input += "\n"

        # Create the inputs folder.
        # First, get the current working directory.
        owd = os.getcwd()
        libefpmbe_folder = cif_output[:-4]

        # Check if the directory exist, and if not, create it.
        try:
            os.mkdir(libefpmbe_folder)

        except FileExistsError:
            pass

        # Go to the directory where the file is supposed to be written.
        os.chdir(libefpmbe_folder)
        
        # Ligands trick is only applied to (N>2)-mers.
        if len(nmer["monomers"]) > 2:
            libefpmbe_filename = "{}-{}.in".format(keynmer, ligands[i])

        # This is the normal way of naming files.
        else:
            libefpmbe_filename = "{}.in".format(keynmer)

        # Create the new input file.
        with open(libefpmbe_filename, "w") as libefpmbe_f:

            for line in libefpmbe_input:
                libefpmbe_f.write(line)

        os.chdir(owd)
def qcarchive_energies(cif_output, nmers, keynmer, nmer, cpus, cle_run_type, psi4_method, psi4_bsse, psi4_memory, qca_client=None, qca_mode='run', verbose=0):
    """
    Arguments:
    
    """

    # If the output is going to be kept, setup the filename.
    if "quiet" or "test" in cle_run_type:
        psi4.core.be_quiet()
        
    # If the output is not kept, do not print to screen.
    else:
        owd = os.getcwd()

        p4folder = cif_output[:-4]
        
        try:
            os.mkdir(p4folder)
        except FileExistsError:
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("WARNING: A folder with the same name as the CIF file already exists.")
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
            pass

        os.chdir(p4folder)
        
        p4out = keynmer + ".dat"
        psi4.core.set_output_file(p4out)
        
        os.chdir(owd)

    psi_api_molecule = nmer2psiapimol(nmers, keynmer, nmer, verbose)
    mymol = psi4.geometry(psi_api_molecule)
    
    # Set the number of threads to run Psi4.
    psi4.core.set_num_threads(cpus)
    
    # Set Psi4 memory.
    psi4.set_memory(psi4_memory)
    
    # Set Psi4 block of options.
    psi4.set_options({'scf_type': 'df', 'mp2_type': 'df', 'cc_type': 'df', 'freeze_core': 'true', 'e_convergence': '8'})
    print('keynmer', keynmer)
    print()
    
    # Execute Psi4 energy calculations, unless running on test mode.

    if "test" not in cle_run_type:
        # Build a plan
        plan = psi4.energy(psi4_method, molecule=mymol, bsse_type=[psi4_bsse], return_plan=True, return_total_data=True)
        # Sumbit jobs in the plan to the QCA jobqueue
        #print(type(plan)) # the type is `ManyBodyComputer'
        plan.compute(qca_client)

        # Query results separate from running them
        if qca_mode == 'query':
            print("Querying results from the QC manager")
            # NOTE: All results are stored in this nested dict,
            #       could dump to a json file if you want

            qca_results = plan.get_results(qca_client)

            # Get the non-additive n-body contribution, exclusive of all
            # previous-body interactions.
            varstring = "{}-CORRECTED {}-BODY INTERACTION ENERGY".format(psi4_bsse.upper(), str(len(nmer["monomers"])))

            
            n_body_energy = qca_results.extras['qcvars']['nbody'][varstring]
            
            if len(nmer["monomers"]) > 2:
                varstring = "{}-CORRECTED {}-BODY INTERACTION ENERGY".format(psi4_bsse.upper(), str(len(nmer["monomers"]) - 1))
                n_minus_1_body_energy = qca_results.extras['qcvars']['nbody'][varstring]
                nmer["nambe"] = n_body_energy - n_minus_1_body_energy

            else:
                nmer["nambe"] = n_body_energy

        else:
            print("Submitted jobs to QC manager")
            exit()
# ======================================================================


# ======================================================================
def psi4api_energies(cif_output, nmers, keynmer, nmer, cpus, cle_run_type, psi4_method, psi4_bsse, psi4_memory, verbose=0):
    """
    Arguments:
    
    """

    # If the output is going to be kept, setup the filename.
    if ("quiet" in cle_run_type) or ("test" in cle_run_type):
        psi4.core.be_quiet()
        
    # If the output is not kept, do not print to screen.
    else:
        owd = os.getcwd()

        p4folder = cif_output[:-4]
        
        try:
            os.mkdir(p4folder)

        except FileExistsError:
            pass

        os.chdir(p4folder)
        
        p4out = keynmer + ".dat"
        psi4.core.set_output_file(p4out, True)
        
        os.chdir(owd)

    psi_api_molecule = nmer2psiapimol(nmers, keynmer, nmer, verbose)
    mymol = psi4.geometry(psi_api_molecule)
    
    # Set the number of threads to run Psi4.
    psi4.core.set_num_threads(cpus)
    
    # Set Psi4 memory.
    psi4.set_memory(psi4_memory)
    
    # Set Psi4 block of options.
    psi4.set_options({'scf_type': 'df', 'mp2_type': 'df', 'cc_type': 'df', 'freeze_core': 'true', 'e_convergence': '10', 'd_convergence': '10'})
    
    # Execute Psi4 energy calculations, unless running on test mode.

    # Example:  psi4.energy('MP2/aug-cc-pV[D,T]Z', molecule=mymol, bsse_type=['cp', 'nocp', 'vmfc'])
    #           psi4.energy('HF/STO-3G', molecule=mymol, bsse_type=['nocp'])
    #           psi4.energy('MP2/aug-cc-pVDZ', molecule=mymol, bsse_type=['vmfc'])
    
    if "test" not in cle_run_type:
        psi4.energy(psi4_method, molecule=mymol, bsse_type=[psi4_bsse])
    
        # This is how you print out all the variables to figure out errors.
        #import pprint
        #pprint.pprint(psi4.core.variables())

        # Get the non-additive n-body contribution, exclusive of all
        # previous-body interactions.
        varstring = "{}-CORRECTED {}-BODY INTERACTION ENERGY".format(psi4_bsse.upper(), str(len(nmer["monomers"])))
    
        n_body_energy = psi4.core.variable(varstring)
    
        if len(nmer["monomers"]) > 2:
            varstring = "{}-CORRECTED {}-BODY INTERACTION ENERGY".format(psi4_bsse.upper(), str(len(nmer["monomers"]) - 1))
            n_minus_1_body_energy = psi4.core.variable(varstring)
            nmer["nambe"] = n_body_energy - n_minus_1_body_energy
        
        else:
            nmer["nambe"] = n_body_energy
# ======================================================================

def run_makefp(monomer):

    # TODO: charge, multiplicity, etc.
    mol = qcel.models.Molecule(symbols=monomer['elem'], geometry=monomer['coords'])
    mol.pretty_print()
    fragment_orientation = distance_matrix(monomer['coords'][:3], monomer['coords'][:3])[0]

    # small basis
    inp_small = qcel.models.AtomicInput(
        molecule=mol,
        driver="energy",
        model={"method": "hf-makefp", "basis": "N31"},
        keywords={'basis__ngauss' : '6', 'basis__ndfunc' : '1', 'damp__ifttyp(1)' : '2,0', 'damp__thrsh' : '500.0'},
    )

    # large basis
    inp_large = qcel.models.AtomicInput(
        molecule=mol,
        driver="energy",
        model={"method": "hf-makefp", "basis": "N311"},
        keywords={'basis__ngauss' : '6', 'basis__npfunc': '2', 'basis__ndfunc' : '3', 'basis__nffunc' : '1', 'basis__diffs' : '.t.', 'basis__diffsp' : '.t.'},
    )

    res_small = qcng.compute(inp_small, program="gamess", raise_error=True, return_dict=True, local_options={"ncores":6})
    assert res_small["success"] is True
    fragment_potential_small_lines = res_small["extras"]["outfiles"]["gamess.efp"].split('\n')

    res_large = qcng.compute(inp_large, program="gamess", raise_error=True, return_dict=True, local_options={"ncores":6})
    assert res_large["success"] is True
    fragment_potential_large_lines = res_large["extras"]["outfiles"]["gamess.efp"].split('\n')

    print(len(fragment_potential_small_lines))
    print(len(fragment_potential_large_lines))

    ind1_small, ind1_large = 0, 0

    for i, line in enumerate(fragment_potential_small_lines):
        if line.strip().startswith('POLARIZABLE POINTS'):
            ind1_small = i
    for i, line in enumerate(fragment_potential_large_lines):
        if line.strip().startswith('POLARIZABLE POINTS'):
            ind1_large = i
    print(ind1_small, ind1_large)

    ind2_small, ind2_large = 0, 0

    for i, line in enumerate(fragment_potential_small_lines):
        if line.strip().startswith('SCREEN2'):
            ind2_small = i
    for i, line in enumerate(fragment_potential_large_lines):
        if line.strip().startswith('SCREEN2'):
            ind2_large = i
    print(ind2_small, ind2_large)

    fragment_potential = fragment_potential_small_lines[2:ind1_small] + fragment_potential_large_lines[ind1_large:ind2_large] + fragment_potential_small_lines[ind2_small:]

    # remove the first two lines of the efp file
    fragment_potential = '\n'.join(fragment_potential)

    # TODO remove this
    with open('zach_gamess.efp', 'w') as fp:
        fp.write(fragment_potential)

    return fragment_potential, fragment_orientation


# ======================================================================
def cle_manager(cif_output, nmers, cle_run_type, psi4_method, psi4_bsse, psi4_memory, qca_server_uri=None, qca_mode='run',verbose=0):
    """Manages which mode of CrystaLattE calculation will be employed.
    
    Arguments:
    <dict> nmers
        Dictionary containing dictionaries for each N-mer in the 
        system.
    <list> cle_run_type
        List of keywords indicating the modes in which CrystaLattE is
        going to be run.
    <str> psi4_method
        Method and basis set for the energy calculation, separated by a
        slash.
    <str> psi4_bsse
        Method for correction of the basis set superposition error.
    <str> psi4_memory
        Memory allocation for PSI4, written as `500 MB`, or `60 GB`.
    <int> verbose
        Adjusts the level of detail of the printouts.

    Returns:
    <float> crystal_lattice_energy
        The value of the accumulated crystal lattice energy in atomic
        units.
    <list> results
        A summary of the results for all N-mers including energies,
        replicas, contributions, cumulative lattice energy, priority,
        and minimum atomic separations; for the print_results function.
    """

    assert "gamessefp" in cle_run_type

    # Find out the number of CPUs in the local system.
    cpus = multiprocessing.cpu_count()

    crystal_lattice_energy = 0.0
    results = []
    csv_lines = []
    csv_header = "N-mer Name,"\
            + "Non-Additive MB Energy (kJ/mol),"\
            + "Num. Rep. (#), N-mer Contribution (kJ/mol),"\
            + "Partial Crys. Lattice Ener. (kJ/mol),"\
            + "Calculation Priority (Arb. Units),"\
            + "Minimum Monomer Separations (A)"
    csv_lines.append(csv_header)

    print('\nRunning MAKEFP on monomer ...')
    tstart_makefp = time.time()
    fragment_potential, fragment_orientation = run_makefp(nmers['1mer-0'])
    dt_makefp = time.time() - tstart_makefp
    print(f'... Done in {dt_makefp:.1f} seconds\n')

    # Get the keys of the N-mers dictionary, and put them on a list.
    nmer_keys = list(nmers.keys())
    nmer_keys = [nmer_key for nmer_key in nmer_keys if len(nmers[nmer_key]['monomers']) > 1]

    # Sort the list in decreasing priority order.
    nmer_keys.sort(key = lambda x: -nmers[x]['priority_cutoff'])

    # The next line was replaced to trigger the calculations in
    # priority order, and not in the order the N-mers were created.
    #for keynmer, nmer in nmers.items():

    gamess_inps = [(fragment_potential, fragment_orientation, nmers[keynmer], keynmer) for keynmer in nmer_keys]
    gamess_outs = multiprocessing.Pool(24).starmap(gamessefp_nambe, gamess_inps)

    #for gamess_inp in gamess_inps:
    #    #nmer = nmers[keynmer]
    #    #nmer["nambe"] = gamessefp_nambe(fragment_potential, fragment_orientation, nmers[keynmer], keynmer)
    #    gamess_outs.append(gamessefp_nambe(*gamess_inp))

    for i, keynmer in enumerate(nmer_keys):
        nmers[keynmer]["nambe"] = gamess_outs[i]

    for keynmer in nmer_keys:

        nmer = nmers[keynmer]
        nmer["contrib"] = nmer["nambe"] * nmer["replicas"] / float(len(nmer["monomers"]))
        crystal_lattice_energy += nmer["contrib"]

        # Generate a string with an ordered list of minimum separations
        # between atoms belonging to different monomers.
        rminseps = ""
        nmer_min_monomer_separations = nmer["min_monomer_separations"]
        nmer_min_monomer_separations.sort()
        for r in nmer_min_monomer_separations:
            rminseps += "{:6.3f} ".format(r * qcel.constants.bohr2angstroms)

        nmer_result = "{:26} | {:>12.8f} | {:>4} | {:>12.8f} | {:>13.8f} | {:12.6e} | {}".format(
                keynmer, 
                nmer["nambe"] * qcel.constants.hartree2kcalmol * qcel.constants.cal2J, 
                nmer["replicas"], 
                nmer["contrib"] * qcel.constants.hartree2kcalmol * qcel.constants.cal2J,
                crystal_lattice_energy * qcel.constants.hartree2kcalmol * qcel.constants.cal2J,
                nmer["priority_cutoff"],
                rminseps)
        
        results.append(nmer_result)

        nmer_csv = "{:},{:.8f},{:},{:.8f},{:.8f},{:.6e},{}".format(
                keynmer,
                nmers[keynmer]["nambe"] * qcel.constants.hartree2kcalmol * qcel.constants.cal2J,
                nmers[keynmer]["replicas"],
                nmers[keynmer]["contrib"] * qcel.constants.hartree2kcalmol * qcel.constants.cal2J,
                crystal_lattice_energy * qcel.constants.hartree2kcalmol * qcel.constants.cal2J,
                nmers[keynmer]["priority_min"],
                rminseps)
        
        csv_lines.append(nmer_csv)
 
    with open('Results.csv', 'w') as csvf:
        for line in csv_lines:
            csvf.write(line + "\n")



    return crystal_lattice_energy, results
# ======================================================================


# ======================================================================
def print_header():
    """.
    """

    print("")
    print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
    print("                              CrystaLattE                              \n")
    print("  The tool for the automated calculation of crystal lattice energies.  \n")
    print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
# ======================================================================


# ======================================================================
def print_results(results, crystal_lattice_energy, verbose=0):
    """Prints a summary of the energy results at the end of the
    execution.
    
    Arguments:
    <int> verbose
        Adjusts the level of detail of the printouts.
    """

    if verbose >= 1:
        print("Summary of results:")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}".format("-"*(shutil.get_terminal_size().columns - 112)))
        print("                           | Non-Additive | Num. |        N-mer | Partial Crys. |  Calculation | Minimum Monomer")
        print("N-mer Name                 |    MB Energy | Rep. | Contribution | Lattice Ener. |     Priority | Separations")
        print("                           |     (kJ/mol) |  (#) |     (kJ/mol) |      (kJ/mol) | (Arb. Units) | (A)")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}".format("-"*(shutil.get_terminal_size().columns - 112)))
        for result in results:
            print(result)
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}\n".format("-"*(shutil.get_terminal_size().columns - 112)))
        #print("Crystal Lattice Energy (Eh)       = {:5.8f}".format(crystal_lattice_energy))
        print("Crystal Lattice Energy (kJ/mol)   = {:9.8f}".format(crystal_lattice_energy * qcel.constants.hartree2kcalmol * qcel.constants.cal2J))
        print("Crystal Lattice Energy (kcal/mol) = {:9.8f}\n".format(crystal_lattice_energy * qcel.constants.hartree2kcalmol))
# ======================================================================


# ======================================================================
def print_end_msg(start, verbose=0):
    """Prints a success message and timing information at the end of the
    execution.
       
    Arguments:
    <int> verbose
        Adjusts the level of detail of the printouts.
    """
    
    if verbose >= 1:
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
        print("Execution terminated succesfully.")
        print("Total elapsed wall-clock time: {:.2f} s\n".format(time.time() - start))
        print("Thank you for using CrystaLattE.\n")
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

# ======================================================================


# ======================================================================
def main(cif_input, cif_output="sc.xyz", cif_a=5, cif_b=5, cif_c=5, bfs_thresh=1.2, uniq_filter="ChSEV", nmers_up_to=2, r_cut_com=10.0, r_cut_monomer=12.0, r_cut_dimer=10.0, r_cut_trimer=8.0, r_cut_tetramer=6.0, r_cut_pentamer=4.0, cle_run_type=["test"], psi4_method="HF/STO-3G", psi4_bsse="cp", psi4_memory="500 MB", qca_server_uri=None, qca_mode='run', verbose=1):
    """Takes a CIF file and computes the crystal lattice energy using a
    many-body expansion approach.
    """
   
    # Start counting time.
    start = time.time()

    # Read a CIF file and generate the unit cell.
    cif_arguments = cif_driver(cif_input, cif_output, cif_a, cif_b, cif_c, verbose)
    cif_main(cif_arguments)
    
    # Read the output of the automatic fragmentation.
    nmers = supercell2monomers(cif_output, r_cut_monomer, bfs_thresh, verbose)
    total_monomers = len(nmers)

    # If makefp mode requested, produce a MAKEFP input file for gamess
    # and exit.
    if "makefp" in cle_run_type:

        monomer2makefp(cif_output, nmers["1mer-0"], verbose)
        
        print("\nThe makefp file for {} has been created.\n".format(cif_output[:-4]))
        
        print_end_msg(start, verbose)
        sys.exit(0)

    # Loop through all monomers and monomers in the central unit cell
    # to generate dimers with at least one monomer in the central cell.
    # Then loop through all existing N-mers and generate higher-order 
    # N-mers with all monomers.

    if nmers_up_to < 2:
        print("\nERROR: CrystaLattE is designed to use at least dimers.")
        print("       Please use 2 <= nmer_up_to < 5.")
        sys.exit()
    
    if nmers_up_to >= 2:
        
        if verbose >= 2:
            print("\nMerging monomers with monomers to obtain dimers.")
        build_nmer(nmers, total_monomers, "dimers", r_cut_dimer, r_cut_com, uniq_filter, verbose)

    if nmers_up_to >= 3:
        
        if verbose >= 2:
            print("\nMerging dimers with monomers to obtain trimers.")
        build_nmer(nmers, total_monomers, "trimers", r_cut_trimer, r_cut_com, uniq_filter, verbose)

    if nmers_up_to >= 4:
        
        if verbose >= 2:
            print("\nMerging trimers with monomers to obtain tetramers.")
        build_nmer(nmers, total_monomers, "tetramers", r_cut_tetramer, r_cut_com, uniq_filter, verbose)

    if nmers_up_to == 5:
        
        if verbose >= 2:
            print("\nMerging tetramers with monomers to obtain pentamers.")
        build_nmer(nmers, total_monomers, "pentamers", r_cut_pentamer, r_cut_com, uniq_filter, verbose)

    if nmers_up_to > 5:
        print("\nERROR: The current implementation of CrystaLattE is limited to pentamers.")
        print("       Please use 2 <= nmer_up_to < 5.")
        sys.exit()
   
    if verbose >= 2:
        
        if "psi4api" in cle_run_type:
            print ("\nComputing interaction energies of N-mers:")
        
        if "psithon" in cle_run_type:
            print ("\nWriting N-mers coordinates to Psithon input files:")

        if "libefpmbe" in cle_run_type:
            print ("\nWriting N-mers coordinates as a non-embedded framgents into LibEFP input files:")

    crystal_lattice_energy, results = cle_manager(cif_output, nmers, cle_run_type, psi4_method, psi4_bsse, psi4_memory, qca_server_uri, qca_mode, verbose)
    # ------------------------------------------------------------------

    if verbose >= 2:
        print("")

    # Print the final results.
    print_results(results, crystal_lattice_energy, verbose)
    
    # Print exit message and timings information.
    print_end_msg(start, verbose)

    # Debug only.
    #import pprint
    #pprint.pprint(nmers)

    return nmers, crystal_lattice_energy

# ======================================================================


if __name__ == "__main__":

    # Hard-coded Test
    if "crystalatte.py" in sys.argv[-1]:
        
        print_header()

        main(   cif_input="../Tests/Ammonia/Ammonia.cif",
                cif_output="../Tests/Ammonia/Ammonia.xyz",
                cif_a=3,
                cif_b=3,
                cif_c=3,
                bfs_thresh=1.2,
                uniq_filter="ChSEV",
                nmers_up_to=2,
                r_cut_com=6.5,
                r_cut_monomer=3.5,
                r_cut_dimer=2.6,
                r_cut_trimer=3.7,
                r_cut_tetramer=3.7,
                r_cut_pentamer=6.1,
                cle_run_type=["psi4api", "quiet"],
                psi4_method="HF/STO-3G",
                psi4_bsse="nocp",
                psi4_memory="500 MB",
                qca_server_uri=None,
                qca_mode='run',
                verbose=2)

        print("{}".format("~"*(shutil.get_terminal_size().columns)))
        print("WARNING: No input was provided. The previous execution was just a test.")
        print("{}".format("~"*(shutil.get_terminal_size().columns)))
        print("\nCrystaLattE execution command:  ./crystalatte.py YourInput.cle\n")

    # Normal execution using an input file.
    else:
        input_parser(sys.argv[-1])
