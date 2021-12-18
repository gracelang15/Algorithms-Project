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
