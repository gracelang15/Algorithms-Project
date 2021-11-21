def encode_output(container, sequence):
    rna_structure = ["."] * len(sequence)
    for pair in container:
        rna_structure[pair[0]] = "("
        rna_structure[pair[1]] = ")"
    rna_structure = ''.join(rna_structure)
    return rna_structure


def convert_input(rna_true_structure):
    first_pair_element = []
    second_pair_element = []
    for i in range(len(rna_true_structure)):
        if rna_true_structure[i] == "(":
            first_pair_element.append(i)
        elif rna_true_structure[i] == ")":
            second_pair_element.append(i)
    pair_index = list(zip(first_pair_element, second_pair_element))
    return pair_index
