import numpy as numpy
from base_pairing import base_pairing

RNA_length = 10
THRESHOLD = 3  # Avoid sharp turn.


# THRESHOLD = 3 means that two bases can pair only if they are at least 3 bases away from each other.

def rna_base_pairing(sequence):
    # For an RNA sequence, calculate a matrix of the maximum base pairing number
    total_bp_m = numpy.zeros(shape=(RNA_length, RNA_length))  # Matrix M in the Nussinov et al. 1980 paper
    # Scan through subsequences with the length from THRESHOLD + 1 to total RNA length
    for dist in range(THRESHOLD + 1, RNA_length):
        for base_start in range(0, RNA_length - dist):
            base_end = base_start + dist
            total_bp_m = max_bp(sequence, base_start, base_end, total_bp_m)
    return total_bp_m


def max_bp(sequence, bi, bj, total_bp_m):
    # This function is used to make the matrix M(i,j) in the Nussinov et al. 1980 paper
    # bi, bj: integers, represent the Bi and Bj in the Nussinov et al. 1980 paper
    for bk in range(bi, bj):
        bp = base_pairing(sequence[bk], sequence[bj])
        if bp == 0:  # bk cannot pair with bj
            continue  # Move to the next bk
        elif (bp == 1) and \
                (total_bp_m[bi, bj] < 1 + total_bp_m[bi, bk - 1] + total_bp_m[bk + 1, bj - 1]):
            # Current bk is a better choice
            total_bp_m[bi, bj] = 1 + total_bp_m[bi, bk - 1] + total_bp_m[bk + 1, bj - 1]
            # Saving the base pairing: bk pairs with bj
            total_bp_m[bj, bi] = bk  # This is the Kij (= Mji). bk pairs with bj
        else:  # Input error for base_pairing
            print(f'Input error for base_pairing()! {bp}')
            break  # Report the error and break the loop

    if total_bp_m[bi, bj - 1] > total_bp_m[bi, bj]:  # It is better that bj does not pair with any base
        total_bp_m[bi, bj] = total_bp_m[bi, bj - 1]
        total_bp_m[bj, bi] = 0
    return total_bp_m

    # Press the green button in the gutter to run the script.


if __name__ == '__main__':
    print("hi")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
