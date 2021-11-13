import numpy as numpy
import pandas as pd
import sys
from base_pairing import base_pairing
numpy.set_printoptions(threshold=sys.maxsize)

#RNA_length = 15
THRESHOLD = 3  # Avoid sharp turn.
## actually we could use something like len(RNA) to get RNA's length 
## have some advice to this part, we could simply use a map structue to store those basic pairing pairs 
## and use tolowercase to preprocess the character 

#sequence = "GUCUACGGCCAUACC"
# THRESHOLD = 3 means that two bases can pair only if they are at least 3 bases away from each other.

pairdict = {
  'a': 'u',
  'u': 'a',
  'g': 'c',
  'c': 'g'
}
## solve the pairing problem 

## create structure used to store M(i,j) and K(i,j)
#total_bp_m = numpy.zeros(shape=(RNA_length, RNA_length))  # Matrix M in the Nussinov et al. 1980 paper
#bp_m = numpy.zeros(shape=(RNA_length, RNA_length))    # Matrix K
#bp_m.fill(-1)

##currently not used 
def check_pairing(k, j, sequence):
    if abs(j-k)>THRESHOLD:
        return pairdict.get(sequence[k].lower())==sequence[j].lower() ## actuall we could perprocess the dna
    return 0


def rna_base_pairing(sequence):
    # For an RNA sequence, calculate a matrix of the maximum base pairing number
    # Scan through subsequences with the length from THRESHOLD + 1 to total RNA length
    total_bp_m = numpy.zeros(shape=(len(sequence), len(sequence)))  # Matrix M in the Nussinov et al. 1980 paper
    bp_m = numpy.zeros(shape=(len(sequence), len(sequence)))  # Matrix K
    bp_m.fill(-1)
    for dist in range(THRESHOLD + 1, len(sequence)):
        for base_start in range(0, len(sequence) - dist):
            base_end = base_start + dist
            max_bp(base_start, base_end, sequence, total_bp_m, bp_m)
    container = traceback(0, len(sequence), bp_m)
    rna_structure = encode_output(container, sequence)
    return rna_structure
    #return total_bp_m, bp_m


def max_bp(bi, bj, sequence, total_bp_m, bp_m):
    # This function is used to make the matrix M(i,j) in the Nussinov et al. 1980 paper
    # bi, bj: integers, represent the Bi and Bj in the Nussinov et al. 1980 paper
    for bk in range(bi, bj):
        #bp = base_pairing(sequence[bk], sequence[bj])
        bp = check_pairing(bk, bj, sequence)
        # instead of passing value, passing the index of characters inside the array, and let the fuction to judge whether they distance overpass the therold 

        ##print(bp)
        if bp == 0:  # bk cannot pair with bj
            continue  # Move to the next bk
        elif (bp == 1) and \
                (total_bp_m[bi, bj] < 1 + total_bp_m[bi, bk - 1] + total_bp_m[bk + 1, bj - 1]):
            # Current bk is a better choice
            total_bp_m[bi, bj] = 1 + total_bp_m[bi, bk - 1] + total_bp_m[bk + 1, bj - 1]
            # Saving the base pairing: bk pairs with bj
            bp_m[bi, bj] = bk  # This is the Kij (= Mji). bk pairs with bj
        #else:  # Input error for base_pairing
        #    print(f'Input error for base_pairing()! {bp}')
        #    break  # Report the error and break the loop

    if total_bp_m[bi, bj - 1] > total_bp_m[bi, bj]:  # It is better that bj does not pair with any base
        total_bp_m[bi, bj] = total_bp_m[bi, bj - 1]
        bp_m[bi, bj] = -1
    #return total_bp_m, bp_m
    # Press the green button in the gutter to run the script.


def traceback(i, j, bp_m):
    container=[]
    length=j-1
    start=i
    def tracebackstep(i, j):
        #print("call trace step ")
        while j>i and bp_m[i][j]==-1:
            j=j-1
        if i>=j:
            return
        pairValue= int(bp_m[i][j])
        #print("current pair is: " + str(pairValue))
        container.append((pairValue,j))
        tracebackstep(i,pairValue-1)
        tracebackstep(pairValue+1,j-1)

    tracebackstep(start,length)
    #print("called")
    return container

def encode_output(container, sequence):
    rna_structure = ["."] * len(sequence)
    for pair in container:
        rna_structure[pair[0]] = "("
        rna_structure[pair[1]] = ")"
    rna_structure = ''.join(rna_structure)
    return rna_structure


if __name__ == '__main__':
    rna_data = pd.read_excel("./RNAData.xlsx")
    rna_data["Calculated Structure"] = rna_data["RNA_sequence"].apply(lambda x: rna_base_pairing(x))
    print(rna_data["Calculated Structure"])
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
