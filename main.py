import numpy as numpy
import pandas as pd
import sys
numpy.set_printoptions(threshold=sys.maxsize)
import time

#RNA_length = 15
THRESHOLD = 3  # Avoid sharp turn.
## actually we could use something like len(RNA) to get RNA's length 
## have some advice to this part, we could simply use a map structue to store those basic pairing pairs 
## and use tolowercase to preprocess the character 

#sequence = "GUCUACGGCCAGAAA"
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
    #rna_structure = encode_output(container, sequence)
    return container

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


def traceback(i, j, bp_m):
    container=[]
    length=j-1
    start=i

    def tracebackstep(i, j):
        while j>i and bp_m[i][j]==-1:
            j=j-1
        if i>=j:
            return
        pairValue= int(bp_m[i][j])
        container.append((pairValue,j))
        tracebackstep(i,pairValue-1)
        tracebackstep(pairValue+1,j-1)

    tracebackstep(start,length)
    return container


def encode_output(container, sequence):
    rna_structure = ["."] * len(sequence)
    for pair in container:
        rna_structure[pair[0]] = "("
        rna_structure[pair[1]] = ")"
    rna_structure = ''.join(rna_structure)
    return rna_structure


def convert_input2(rna_true_structure):
    save_element = []
    pair_index = []
    for i in range(len(rna_true_structure)):
        if rna_true_structure[i] == ".":
            continue
        elif rna_true_structure[i] == "(":
            save_element.append(i)
        else:
            # When you firstly meet a ")", pair it with the closest unpaired "(" on the left side of
            # the current ")". Then remove this "(" from the list of the unpaired "(".
            pair_index.append((save_element[len(save_element)-1], i))
            save_element.pop()
    return pair_index


def calculate_distances(rna_true_pairing, rna_predicted_pairing):
    s1 = []
    s2 = []
    for pair1 in rna_predicted_pairing:
        # print(pair1)
        pair1_distances = []
        for true_pair in rna_true_pairing:
            pair1_distances.append(max(abs(pair1[0] - true_pair[0]), abs(pair1[1] - true_pair[1])))
        s1.append(min(pair1_distances))
    for pair2 in rna_true_pairing:
        pair2_distances = []
        for pred_pair in rna_predicted_pairing:
            pair2_distances.append(max(abs(pair2[0] - pred_pair[0]), abs(pair2[1] - pred_pair[1])))
        s2.append(min(pair2_distances))
    delta = s2 + s1
    delta.sort()
    delta = delta[::-1]
    return delta


def calculate_RBP_score(t, rna_distances):
    #check if there are no cases where delta <= t*m:
    if rna_distances[-1] > (len(rna_distances)-1)*t:
        print("found_case")
        return len(rna_distances)-1
    else:
        m=0
        delta = rna_distances[m]
        while delta > t*m:
            m=m+1
            delta = rna_distances[m]
        return m

if __name__ == '__main__':
    rna_data = pd.read_excel("./RNADataSubset.xlsx")
    print(rna_data)
    #use line below when wanting to do a quick test on one row of data
    #rna_data = pd.read_excel("./RNAData2.xlsx", usecols="A:D", sheet_name=1)
    rna_data["RNA_true_base_pairing"] = rna_data["RNA_structure"].apply(lambda x: convert_input(x))
    start = time.perf_counter()
    rna_data["RNA_predicted_base_pairing"] = rna_data["RNA_sequence"].apply(lambda x: rna_base_pairing(x))
    end = time.perf_counter()
    print("time_elapsed:", end - start)
    print("starting scoring")
    rna_data["RNA_predicted_structure"] = rna_data.apply(lambda x: encode_output(x.RNA_predicted_base_pairing, x.RNA_sequence),
                                                         axis=1)
    rna_data["RNA_distances"] = rna_data.apply(lambda x: calculate_distances(x.RNA_true_base_pairing, x.RNA_predicted_base_pairing),
                                                         axis=1)
    rna_data["RBP_score"] = rna_data.apply(
        lambda x: calculate_RBP_score(0, x.RNA_distances),
        axis=1)
    print("done!")
    print(rna_data)

    #rna_data.to_excel("results_nuss_BP.xlsx", sheet_name="nussinov")
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
