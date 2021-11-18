import numpy as numpy
import math as math
import pandas as pd
import sys


sequence = "GUCUACGGCCAUACC"

THRESHOLD = 3  # Avoid sharp turn.
GAS_CONSTANT = 1.9872e-3
TEMP = 310.15  # Body temperature

pairdict = {
    'a': 'u',
    'u': 'a',
    'g': 'c',
    'c': 'g'
}
'''
Theoretically in this algorithm, GU pair is also allowed if it is not a terminal base pair. 
Here to simplify we ignore this case.
'''

# solve the pairing problem


def check_pairing(k, j, sequence):
    if abs(j-k) > THRESHOLD:
        return pairdict.get(sequence[k].lower()) == sequence[j].lower()  # actuall we could perprocess the dna
    return 0


def fh(i, j, sequence):  # Hairpin loop free energy
    # sequence[i] and sequence[j] must base-pair with each other, otherwise this function should not be called
    # if abs(j-i) <= THRESHOLD:
    #    return float('inf')  # If the hairpin loop has fewer than four exterior edges, then E = inf.
    """
    Hairpin loop energy number from Salser (1978). "Globin mRNA Sequences: Analysis of Base Pairing and
    Evolutionary Implications".

    For energy of loops longer than 10, modified extrapolation equation from
    https://github.com/Lattice-Automation/seqfold/blob/master/seqfold/fold.py
    to fit the numbers in Salser (1978):

    gas_constant = 1.9872e-3
    query_energy = known_energy + 1.31 * gas_constant * temp * math.log(query_len / float(known_len))

    See Jacobson and Stockmayer (2004), SantaLucia and Hicks (2004) for more information.
    """
    # Hairpin loop closed by a GC pair:
    hg = {
        1: float('inf'),
        2: float('inf'),
        3: 8.40,
        4: 5.90,
        5: 4.10,
        6: 4.30,
        7: 4.50,
        8: 4.60,
        9: 4.80,
        10: 4.89
    }

    # Hairpin loop closed by an AU pair:
    ha = {
        1: float('inf'),
        2: float('inf'),
        3: 8.00,
        4: 7.50,
        5: 6.90,
        6: 6.40,
        7: 6.60,
        8: 6.80,
        9: 6.90,
        10: 7.00
    }

    if sequence[i].lower() == 'g' or sequence[i].lower() == 'c':
        if abs(j-i-1) > 10:  # Loop longer than 10 nt
            return float(hg.get(10)) + 1.31 * GAS_CONSTANT * TEMP * math.log((j-i-1)/10)
        else:
            return float(hg.get(j-i-1))
    else:
        if abs(j-i-1) > 10:  # Loop longer than 10 nt
            return float(ha.get(10)) + 1.31 * GAS_CONSTANT * TEMP * math.log((j-i-1)/10)
        else:
            return float(ha.get(j-i-1))


def stacking(i, ii, sequence):
    # Call this function only if i < ii < jj < j and i pairs with j and ii pairs with jj
    fs = {
        'aa': -1.2,
        'ac': -2.1,
        'ag': -2.1,
        'au': -1.8,
        'ca': -2.1,
        'cc': -4.8,
        'cg': -3.0,
        'cu': -2.1,
        'ga': -2.1,
        'gc': -4.3,
        'gg': -4.8,
        'gu': -2.1,
        'ua': -1.8,
        'uc': -2.1,
        'ug': -2.1,
        'uu': -1.2
    }
    stack_nt = ''.join([sequence[i].lower(), sequence[ii].lower()])
    return fs.get(stack_nt)


def bulge(i, j, sequence):  # Make sure that j > i
    """
       Bulge loop energy number from Salser (1978). "Globin mRNA Sequences: Analysis of Base Pairing and
       Evolutionary Implications".
       For energy of loops longer than 10, modified extrapolation equation from
       https://github.com/Lattice-Automation/seqfold/blob/master/seqfold/fold.py
       to fit the numbers in Salser (1978).
       See Jacobson and Stockmayer (2004), SantaLucia and Hicks (2004) for more information.
       """
    b = {
        1: 2.80,
        2: 3.90,
        3: 4.45,
        4: 5.00,
        5: 5.15,
        6: 5.30,
        7: 5.45,
        8: 5.60,
        9: 5.69,
        10: 5.78
    }
    if abs(j-i-1) > 10:  # Loop longer than 10 nt
        energy = float(b.get(10)) + 1.31 * GAS_CONSTANT * TEMP * math.log((j-i-1)/10)
    else:
        energy = float(b.get(j-i-1))
    energy = energy + stacking(i, j, sequence)  # Including the stacking energy
    return energy


def int_loop(i, j, sequence):  # Interior loops: j-i>1
    il = numpy.array([[0.10, 0.90, 1.60, 2.10, 2.50, 2.62, 2.72, 2.82, 2.90],
                      [0.95, 1.75, 2.45, 2.95, 3.35, 3.47, 3.57, 3.67, 3.75],
                      [1.80, 2.60, 3.30, 3.80, 4.20, 4.32, 4.42, 4.51, 4.60]])
    '''
    Line 0: energy when both stems on two sides of the interior loop are closed by GC pairs
    Line 1: by one GC and one AU pair
    Line 2: by both AU pairs
    '''
    il_closure = {  # Judging the base pairs on the two sides and choose the correct line
        'aa': 2,
        'ac': 1,
        'ag': 1,
        'au': 2,
        'ca': 1,
        'cc': 0,
        'cg': 0,
        'cu': 1,
        'ga': 1,
        'gc': 0,
        'gg': 0,
        'gu': 1,
        'ua': 2,
        'uc': 1,
        'ug': 1,
        'uu': 2
    }
    line = il_closure.get(''.join([sequence[i].lower(), sequence[j].lower()]))

    if abs(j-i-1) > 10:  # Loop longer than 10 nt
        return il[line, 8] + 1.31 * GAS_CONSTANT * TEMP * math.log((j-i-1)/10)
    else:
        return il[line, j-i-1-2]


def fl(i, j, ii, jj, sequence):
    """
    Calculate the energy of the face within i, ii, jj, and j, where i < ii < jj < j.
    here ii and jj must pair with each other, otherwise the function should not be called.

    The face can be:
    1. stacking region, where ii=i+1, jj=j-1
    2. bulge loop, where only ii=i+1 or only jj=j-1
    3. interior loop, where ii>i+1 and jj<j-1
    """

    if (ii == i+1) and (jj == j-1):  # Stacking region
        return stacking(i, ii, sequence)
    elif (ii > i+1) and (jj < j-1):  # two interior loops: i--ii, jj--j
        return int_loop(i, ii, sequence) + int_loop(jj, j, sequence)
    else:  # either ii=i+1 or jj=j+1: bulge loop
        if ii == i+1:
            return bulge(jj, j, sequence)
        else:
            return bulge(i, ii, sequence)


def zuker(sequence):
    # Calculate the V and W matrix of the Zuker algorithm
    # Scan through subsequences with the length from THRESHOLD + 1 to total RNA length
    v_energy = numpy.zeros(shape=(len(sequence), len(sequence)))  # Matrix V in the Zuker et al. 1981 paper
    w_energy = numpy.zeros(shape=(len(sequence), len(sequence)))
    for dist in range(THRESHOLD + 1, len(sequence)):
        for base_start in range(0, len(sequence) - dist):
            base_end = base_start + dist
            (v_energy, w_energy) = min_energy(base_start, base_end, sequence, v_energy, w_energy)

    #container = traceback(0, len(sequence), bp_m)
    #rna_structure = encode_output(container, sequence)
    #return rna_structure
    return v_energy, w_energy


def min_energy(bi, bj, sequence, v_energy, w_energy):
    # Calculate matrix V
    def v_cal(i, j, seq, v, w):
        if not check_pairing(i, j, seq):
            return float('inf')
        else:
            energy_1 = fh(i, j, sequence)  # E1 in the Zuker et al. 1981 paper

            # Calculate E2
            energy_2 = float('inf')  # Initialize E2
            for ii in range(i + 1, j - 1):
                for jj in range(ii + 1, j):
                    if check_pairing(ii, jj, sequence):
                        energy_2_temp = fl(i, j, ii, jj, sequence) + v[ii, jj]
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
        return min(w[i + 1, j], w[i, j - 1], v[i, j], energy_4)

    w_energy[bi, bj] = w_cal(bi, bj, sequence, v_energy, w_energy)

    return v_energy, w_energy

if __name__ == '__main__':
    print(zuker(sequence))