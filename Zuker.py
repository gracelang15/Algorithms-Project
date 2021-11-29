import numpy as numpy
import math as math
import scoring as sc
import formating as fo
import sys
import pandas as pd
#from main import convert_input, rna_base_pairing, encode_output, calculate_distances, calculate_RBP_score


#sequence = "GUCUACGGCCAGAAA"

THRESHOLD = 3  # Avoid sharp turn.
GAS_CONSTANT = 1.9872e-3
TEMP = 310.15  # Body temperature

pairdict = {  # Base pairing rules
    'a': 'u',
    'u': 'a',
    'g': 'c',
    'c': 'g'
}
'''
Theoretically in this algorithm, GU pair is also allowed if it is not a terminal base pair. 
Here to simplify we ignore this case.
'''


# Base pairing judgement
def check_pairing(k, j, sequence):
    if abs(j - k) > THRESHOLD:
        return pairdict.get(sequence[k].lower()) == sequence[j].lower()  # actuall we could perprocess the dna
    return 0


def fh(i, j, sequence):  # Hairpin loop free energy
    # sequence[i] and sequence[j] must base pair with each other, otherwise this function should not be called

    """
    Hairpin loop energy from Salser (1978). "Globin mRNA Sequences: Analysis of Base Pairing and
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
        # These two lines are equivalent to the thresholding process, prohibiting hairpin loops with less than 3 nt
        # in the loop
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
        # These two lines are equivalent to the thresholding process, prohibiting hairpin loops with less than 3 nt
        # in the loop
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
        if abs(j - i - 1) > 10:  # Loop longer than 10 nt: extrapolation
            return float(hg.get(10)) + 1.31 * GAS_CONSTANT * TEMP * math.log((j - i - 1) / 10)
        else:
            return float(hg.get(j - i - 1))
    else:
        if abs(j - i - 1) > 10:  # Loop longer than 10 nt: extrapolation
            return float(ha.get(10)) + 1.31 * GAS_CONSTANT * TEMP * math.log((j - i - 1) / 10)
        else:
            return float(ha.get(j - i - 1))


def stacking(i, ii, sequence):
    """
        Calculating the free energy of two continuous base pairs.
        Call this function only if i < ii < jj < j and i pairs with j and ii pairs with jj
        Stacking energy from Salser (1978). "Globin mRNA Sequences: Analysis of Base Pairing and
        Evolutionary Implications".
    """
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
       Bulge loop energy from Salser (1978). "Globin mRNA Sequences: Analysis of Base Pairing and
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
    if abs(j - i - 1) > 10:  # Loop longer than 10 nt: extrapolation
        energy = float(b.get(10)) + 1.31 * GAS_CONSTANT * TEMP * math.log((j - i - 1) / 10)
    else:
        energy = float(b.get(j - i - 1))
    energy = energy + stacking(i, j, sequence)  # Including the stacking energy
    return energy


def int_loop(i, j, sequence):  # Interior loops: j-i>1
    """
       Interior loop energy from Salser (1978). "Globin mRNA Sequences: Analysis of Base Pairing and
       Evolutionary Implications".
       For energy of loops longer than 10, modified extrapolation equation from
       https://github.com/Lattice-Automation/seqfold/blob/master/seqfold/fold.py
       to fit the numbers in Salser (1978).
       See Jacobson and Stockmayer (2004), SantaLucia and Hicks (2004) for more information.
    """
    il = numpy.array([[0.10, 0.90, 1.60, 2.10, 2.50, 2.62, 2.72, 2.82, 2.90],
                      [0.95, 1.75, 2.45, 2.95, 3.35, 3.47, 3.57, 3.67, 3.75],
                      [1.80, 2.60, 3.30, 3.80, 4.20, 4.32, 4.42, 4.51, 4.60]])

    """ 
        Line 0: energy when both stems on two sides of the interior loop are closed by GC pairs
        Line 1: by one GC and one AU pair
        Line 2: by both AU pairs
    """

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

    if abs(j - i - 1) > 10:  # Loop longer than 10 nt: extrapolation
        return il[line, 8] + 1.31 * GAS_CONSTANT * TEMP * math.log((j - i - 1) / 10)
    else:
        return il[line, j - i - 1 - 2]


def fl(i, j, ii, jj, sequence):
    """
    Calculate the energy of the face within i, ii, jj, and j, where i < ii < jj < j.
    here ii and jj must pair with each other, otherwise the function should not be called.

    The face can be:
    1. stacking region, where ii=i+1, jj=j-1
    2. bulge loop, where only ii=i+1 or only jj=j-1
    3. interior loop, where ii>i+1 and jj<j-1
    """

    if (ii == i + 1) and (jj == j - 1):  # Stacking region
        return stacking(i, ii, sequence)
    elif (ii > i + 1) and (jj < j - 1):  # two interior loops: i--ii, jj--j
        return int_loop(i, ii, sequence) + int_loop(jj, j, sequence)
    else:  # either ii=i+1 or jj=j+1: bulge loop
        if ii == i + 1:
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
            (v_energy, w_energy) = min_energy(base_start, base_end, sequence, v_energy, w_energy)  # Get matrices v, w

    container = traceback_z(0, len(sequence), v_energy, w_energy)  # Get the base pairs
    #rna_structure = fo.encode_output(container, sequence)  # Get the dot-bracket structure
    return container
    # return v_energy, w_energy


def min_energy(bi, bj, sequence, v_energy, w_energy):
    # Calculate matrix V
    def v_cal(i, j, seq, v, w):
        if not check_pairing(i, j, seq):  # i and j do not base pair
            # For energy, return inf
            # For type, return inf as "do not pair"
            return float('inf'), float('inf')
        else:
            # Calculate E1 for hairpin loop in the Zuker et al. 1981 paper.
            # If a hairpin loop has the smallest free energy,
            # record that there is no more base pairs between i and j.
            energy_1 = fh(i, j, sequence)

            # Calculate E2
            # If E2 is the smallest, record that some ii pairs with jj between i and j.
            energy_2 = float('inf')  # Initialize E2
            for ii in range(i + 1, j - 1):
                for jj in range(ii + 1, j):
                    if check_pairing(ii, jj, sequence):
                        energy_2_temp = fl(i, j, ii, jj, sequence) + v[ii, jj]
                        if energy_2 > energy_2_temp:
                            energy_2 = energy_2_temp
                            # turn [ii,jj] into an 1-D integer. pos_2 > 0
                            pos_2 = ii * len(sequence) + jj

            # Calculate E3
            # If E3 is the smallest, no further base pairs can be confirmed.
            # Record the internal position ii
            energy_3 = float('inf')  # Initialize E3
            for ii in range(i + 2, j - 2):
                if energy_3 > w[i + 1, ii] + w[ii + 1, j - 1]:
                    energy_3 = w[i + 1, ii] + w[ii + 1, j - 1]
                    pos_3 = -ii  # ii > 0, and -ii < 0

            v_cache = [energy_1, energy_2, energy_3]
            v_min = min(v_cache)
            """ v_min:
                0: E1, v_type := 0, hairpin loop: no more base pairs. Finish
                1: E2, v_type := pos_2 = ii * len(sequence) + jj > 0
                    The next base pair is at ii--jj, and then go check v[ii, jj]. 
                2: E3, v_type := pos_3 = -ii < 0
                    An internal point for further w search. 
            """
            v_type = v_cache.index(v_min)
            if v_type == 2:  # E3
                v_type = pos_3
            elif v_type == 1:  # E2
                v_type = pos_2

            return v_min, v_type
            # return v_type

    # Calculate matrix W
    def w_cal(i, j, v, w):
        # Calculate E4
        energy_4 = float('inf')  # Initialize E4
        for ii in range(i + 1, j - 1):
            if energy_4 > w[i, ii] + w[ii + 1, j]:
                energy_4 = w[i, ii] + w[ii + 1, j]
                pos = ii

        w_cache = [w[i + 1, j], w[i, j - 1], v[i, j], energy_4]
        w_min = min(w_cache)  # Find the smallest w energy
        """
        For recording the w energy info, we define:
        -1: w[i, j] = w[i + 1, j]
        -2: w[i, j] = w[i, j - 1]
        -3: w[i, j] = v[i, j], and i--j pair. Then go check:
            --v[i, j] = E1: hairpin loop: no more base pairs. Finish
            --v[i, j] = E2: the next base pair is at ii--jj, and then go check v[ii, jj]
            --v[i, j] = E3: a internal point for further w search
        ii: w[i, j] = E4 and record where the internal point point ii is.
        
        This info is recorded in the other half of the w_energy matrix
        """
        w_pos = -w_cache.index(w_min) - 1  # Which type of energy is the smallest?
        if w_pos == -4:
            w_pos = pos

        return w_min, w_pos

    # Save the energy and the corresponding info at the opposite positions of the matrix
    # v_energy[bi, bj] = v_cal(bi, bj, sequence, v_energy, w_energy)
    (v_energy[bi, bj], v_energy[bj, bi]) = v_cal(bi, bj, sequence, v_energy, w_energy)
    (w_energy[bi, bj], w_energy[bj, bi]) = w_cal(bi, bj, v_energy, w_energy)

    return v_energy, w_energy


def traceback_z(i, j, v, w):
    """
    traceback function for getting the base pair information
    Arguments:
    :param i: start
    :param j: the length of the sequence
    :param v: matrix V
    :param w: matrix W
    :return: the container with (i,j) pairs
    """
    container = []
    LENGTH = j

    def tracebackstep_v(i, j, v):
        # print("tracebackstep_v Current i j: " + str(i) + "," + str(j))
        if v[j, i] == 0:
            return
        elif v[j, i] < 0:
            tracebackstep_w(i + 1, -int(v[j, i]))
            tracebackstep_w(-int(v[j, i]) + 1, j - 1)
        else:
            container.append((int(v[j, i] // LENGTH), int(v[j, i] % LENGTH)))
            tracebackstep_v(int(v[j, i] // LENGTH), int(v[j, i] % LENGTH), v)

    def tracebackstep_w(i, j):
        # print("tracebackstep_w Current i j: " + str(i) + "," + str(j))
        if w[j, i] == 0:
            return  # w[j, i] == 0 means the j-i < = THRESHOLD + 1, cannot form a structure
        elif w[j, i] == -1:
            tracebackstep_w(i + 1, j)
        elif w[j, i] == -2:
            tracebackstep_w(i, j - 1)
        elif w[j, i] == -3:
            container.append((i, j))  # w[i,j] == v[i,j] means that i and j base pair
            tracebackstep_v(i, j, v)
            return
        else:
            tracebackstep_w(i, int(w[j, i]))
            tracebackstep_w(int(w[j, i]) + 1, j)

    tracebackstep_w(i, j-1)
    return container

if __name__ == '__main__':
    #use below line for testing with one row of data
    #rna_data = pd.read_excel("./RNAData2.xlsx", usecols="A:D")
    #print(rna_data)
    rna_data = pd.read_excel("D:/Courses/CS5112_Algorithm/data/RNAData_10-410.xlsx")
    rna_data["RNA_true_base_pairing"] = rna_data["RNA_structure"].apply(lambda x: fo.convert_input(x))
    print('done rna true base pairing')
    rna_data["RNA_predicted_base_pairing"] = rna_data["RNA_sequence"].apply(lambda x: zuker(x))
    print('done rna predict base pairing')
    rna_data["RNA_predicted_structure"] = rna_data.apply(
        lambda x: fo.encode_output(x.RNA_predicted_base_pairing, x.RNA_sequence),
        axis=1)
    print('done rna predicted structure')
    rna_data["RNA_distances"] = rna_data.apply(
        lambda x: sc.calculate_distances(x.RNA_true_base_pairing, x.RNA_predicted_base_pairing),
        axis=1)
    print('done rna distances')
    rna_data["RBP_score"] = rna_data.apply(
       lambda x: sc.calculate_RBP_score(1, x.RNA_distances), axis=1)
    print('done rna rbp score')
    print(rna_data)
    rna_data.to_excel("D:/Courses/CS5112_Algorithm/data/results_zuker2.xlsx", sheet_name="zuker")

