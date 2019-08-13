#!/usr/bin/env python

import numpy as np
import qcelemental as qcel
from psi4.driver import qcdb

def atm_psithon2qmol(fname):
    """
    ."""
    
    with open(fname, "r") as fp:
        lines = fp.readlines()
    
    molecule = []
    save = False
    
    for line in lines:        
        line = line.strip()
        
        if line.endswith("}"):
            break

        if save:
            molecule.append(line)

        if line.startswith("molecule") and line.endswith("{"):
            save = True

    molecule = "\n".join(molecule)
    print(molecule)

    qmol = qcdb.Molecule(molecule)
    
    return qmol

def atm_molecular(qmol, c9_abc):

    mol_a = qmol.extract_subsets(1)
    mol_b = qmol.extract_subsets(2)
    mol_c = qmol.extract_subsets(3)

    print("mol_a\n")
    print(mol_a)
    print("mol_b\n")
    print(mol_b)
    print("mol_c\n")
    print(mol_c)

    com_a = np.array(mol_a.center_of_mass())
    com_b = np.array(mol_b.center_of_mass())
    com_c = np.array(mol_c.center_of_mass())

    E_atm = atm_energy(com_a, com_b, com_c, c9_abc)
    
    return E_atm

def atm_energy(coord_a, coord_b, coord_c, c9_abc):
    """
    Calculate the Axilrod-Teller-Muto triple dipole dispersion energy from coordinates.
    
    Parameters
    ----------
    coord_a: ndarray
        (3,) First set of Cartesian coordinates
    coord_b: ndarray
        (3,) Second set of Cartesian coordinates
    coord_c: ndarray
        (3,) Third set of Cartesian coordinates
    c9_abc:  float
        Three-body dispersion parameter, must be in the same units as the coordinates

    Returns
    -------
    E_abc:   float
        Three-body dispersion energy

    Notes
    -----
    B.M. Axilrod, and E. Teller, J. Chem. Phys. 11, 299 (1943)
    DOI: 10.1063/1.1723844
    Y. Muto, J. Phys. Soc. Jpn. 17, 629, (1943)
    DOI: 10.11429/subutsukaishi1927.17.10-11-12_629 
    """

    v_ab = coord_b - coord_a
    v_bc = coord_c - coord_b
    v_ca = coord_a - coord_c

    r_ab = np.linalg.norm(v_ab)
    r_bc = np.linalg.norm(v_bc)
    r_ca = np.linalg.norm(v_ca)

    cos_theta_abc = np.dot(-v_ab, v_bc)/(r_ab*r_bc)
    cos_theta_bca = np.dot(-v_bc, v_ca)/(r_bc*r_ca)
    cos_theta_cab = np.dot(-v_ca, v_ab)/(r_ca*r_ab)
    
    E_abc = c9_abc*(3*cos_theta_abc*cos_theta_bca*cos_theta_cab + 1)/((r_ab*r_bc*r_ca)**3)
    
    return E_abc

fname = "3mer-0+1+2.in"

qmol_abc = atm_psithon2qmol(fname)
print(atm_molecular(qmol_abc, 82657.65)) # c9_abc in a.u. from the JCP (2014).

c9_abc = 8.0/11.0

coord_a = np.array([0.000, 0.000, 0.000])
coord_b = np.array([1.000, 0.000, 0.000])
coord_c = np.array([0.500, np.sqrt(3.0)/2.0, 0.000])

E_abc = atm_energy(coord_a, coord_b, coord_c, c9_abc)
