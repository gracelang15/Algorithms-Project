import numpy as numpy
import pandas as pd
import sys

THRESHOLD = 3  # Avoid sharp turn.

pairdict = {
    'a': 'u',
    'u': 'a',
    'g': 'c',
    'c': 'g'
}


# solve the pairing problem

def check_pairing(k, j, sequence):
    if abs(j - k) > THRESHOLD:
        return pairdict.get(sequence[k].lower()) == sequence[j].lower()  # actuall we could perprocess the dna
    return 0


def zuker(sequence):
    # Calculate the V and W matrix of the Zuker algorithm
    # Scan through subsequences with the length from THRESHOLD + 1 to total RNA length
    v_energy = numpy.zeros(shape=(len(sequence), len(sequence)))  # Matrix V in the Zuker et al. 1981 paper
    w_energy = numpy.zeros(shape=(len(sequence), len(sequence)))
    for dist in range(THRESHOLD + 1, len(sequence)):
        for base_start in range(0, len(sequence) - dist):
            base_end = base_start + dist
            min_energy(base_start, base_end, sequence, v_energy, w_energy)

    container = traceback(0, len(sequence), bp_m)
    rna_structure = encode_output(container, sequence)
    return rna_structure
    # return total_bp_m, bp_m


def min_energy(bi, bj, sequence, v_energy, w_energy):
    # Calculate matrix V
    def v_cal(i, j, seq, v, w):
        if not check_pairing(i, j, seq):
            return float('inf')
        else:
            energy_1 = FH(i, j)  # E1 in the Zuker et al. 1981 paper

            # Calculate E2
            energy_2 = float('inf')  # Initialize E2
            for ii in range(i + 1, j - 1):
                for jj in range(ii + 1, j):
                    energy_2_temp = FL(i, j, ii, jj) + v[ii, jj]
                    if energy_2 > energy_2_temp:
                        energy_2 = energy_2_temp

            # Calculate E3
            energy_3 = float('inf')  # Initialize E3
            for ii in range(i + 2, j - 2):
                if energy_3 > w[i + 1, ii] + w[ii + 1, j - 1]:
                    energy_3 = w[i + 1, ii] + w[ii + 1, j - 1]
            return min(energy_1, energy_2, energy_3)

    v_energy[bi, bj] = v_cal(bi, bj, sequence, v_energy, w_energy)

    # Calculate matrix W
    def w_cal(i, j, seq, v, w):
        # Calculate E4
        energy_4 = float('inf')  # Initialize E4
        for ii in range(i + 1, j - 1):
            if energy_4 > w[i, ii] + w[ii + 1, j]:
                energy_4 = w[i, ii] + w[ii + 1, j]
        return min(w(i + 1, j), w(i, j - 1), v(i, j), energy_4)

    w_energy[bi, bj] = w_cal(bi, bj, sequence, v_energy, w_energy)

    return v_energy, w_energy
