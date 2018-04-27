#!/usr/bin/env python2

import sys

import numpy as np
import numpy.linalg
import qml

from copy import deepcopy

GAUSSIAN_BOHR_TO_ANGS = 0.52917721092
BOHR_TO_ANGS = GAUSSIAN_BOHR_TO_ANGS

HARTREE_TO_KCAL_MOL = 627.5095


def get_energy(nuclear_charges, coordinates):

    e = 0.0

    for i, qi in enumerate(nuclear_charges):
        for j, qj in enumerate(nuclear_charges):

            if (j <= i): continue

            e += qi * qj / (np.linalg.norm(coordinates[i]/BOHR_TO_ANGS \
                    - coordinates[j]/BOHR_TO_ANGS))

    return e



def get_engrad(nuclear_charges, coordinates):

    energy = 0.0
    grad = np.zeros(coordinates.shape)

    for i, qi in enumerate(nuclear_charges):
        for j, qj in enumerate(nuclear_charges):
            
            
            if (j >= i): continue
            # if (j <= i): continue

            dX = (coordinates[i, 0] - coordinates[j, 0])/BOHR_TO_ANGS
            dY = (coordinates[i, 1] - coordinates[j, 1])/BOHR_TO_ANGS
            dZ = (coordinates[i, 2] - coordinates[j, 2])/BOHR_TO_ANGS

            R = np.sqrt(dX*dX + dY*dY + dZ*dZ)

            dF= -1.0 * (qi * qj )/ (R*R*R)

            grad[i, 0] += dF * dX
            grad[i, 1] += dF * dY
            grad[i, 2] += dF * dZ

            grad[j, 0] -= dF * dX
            grad[j, 1] -= dF * dY
            grad[j, 2] -= dF * dZ 

            energy += qi * qj / R

    return energy, grad


def test_engrad():

    dH = 1e-6

    comp = qml.Compound(xyz="water.xyz")

    energy, grad = get_engrad(comp.nuclear_charges, comp.coordinates)

    grad_numm = np.zeros(grad.shape)

    for i in range(len(comp.nuclear_charges)):
        for j in range(3):

            coords_displaced = deepcopy(comp.coordinates)
            coords_displaced[i,j] += dH
            e_plus = get_energy(comp.nuclear_charges, coords_displaced)
            
            coords_displaced = deepcopy(comp.coordinates)
            coords_displaced[i,j] -= dH
            e_minus = get_energy(comp.nuclear_charges, coords_displaced)

            grad_numm[i,j] = (e_plus - e_minus) / (2*dH/BOHR_TO_ANGS)

    assert np.allclose(grad, grad_numm)


if __name__ == "__main__":

    test_engrad()

    filename = sys.argv[1]

    comp = qml.Compound(xyz=filename)

    e = get_energy(comp.nuclear_charges, comp.coordinates)
    print e

    energy, grad = get_engrad(comp.nuclear_charges, comp.coordinates)

    print energy
    print grad
