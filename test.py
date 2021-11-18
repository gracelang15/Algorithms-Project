import numpy as numpy


def fh(i, j, sequence):  # Hairpin loop free energy
    # if abs(j-i) <= THRESHOLD:
    #    return float('inf')  # If the hairpin loop has fewer than four exterior edges, then E = inf.
    '''
        Hairpin loop energy number from Salser (1978). "Globin mRNA Sequences: Analysis of Base Pairing and
        Evolutionary Implications".

        For energy of loops longer than 10, modified extrapolation equation from
        https://github.com/Lattice-Automation/seqfold/blob/master/seqfold/fold.py
        to fit the numbers in Salser (1978):

        gas_constant = 1.9872e-3
        query_energy = known_energy + 1.31 * gas_constant * temp * math.log(query_len / float(known_len))

        See Jacobson and Stockmayer (2004), SantaLucia and Hicks (2004) for more information.
        '''
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
        if abs(j - i - 1) > 10:  # Loop longer than 10 nt
            return float(hg.get(10)) + 1.31 * GAS_CONSTANT * TEMP * math.log((j - i - 1) / 10)
        else:
            return float(hg.get(j - i - 1))
    else:
        if abs(j - i - 1) > 10:  # Loop longer than 10 nt
            return float(ha.get(10)) + 1.31 * GAS_CONSTANT * TEMP * math.log((j - i - 1) / 10)
        else:
            return float(ha.get(j - i - 1))


fs = {
    'AA': -1.2,
    'AC': -2.1,
    'AG': -2.1,
    'AU': -1.8,
    'CA': -2.1,
    'CC': -4.8,
    'CG': -3.0,
    'CU': -2.1,
    'GA': -2.1,
    'GC': -4.3,
    'GG': -4.8,
    'GU': -2.1,
    'UA': -1.8,
    'UC': -2.1,
    'UG': -2.1

}

il = numpy.array([[0.10, 0.90, 1.60, 2.10, 2.50, 2.62, 2.72, 2.82, 2.90],
                  [0.95, 1.75, 2.45, 2.95, 3.35, 3.47, 3.57, 3.67, 3.75],
                  [1.80, 2.60, 3.30, 3.80, 4.20, 4.32, 4.42, 4.51, 4.60]])
if __name__ == '__main__':
    print(il)
