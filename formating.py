def encode_output(container, sequence):
    rna_structure = ["."] * len(sequence)
    for pair in container:
        rna_structure[pair[0]] = "("
        rna_structure[pair[1]] = ")"
    rna_structure = ''.join(rna_structure)
    return rna_structure


def convert_input(rna_true_structure):
    # This is not correct for getting base pairs: just keep for a record
    first_pair_element = []
    second_pair_element = []
    for i in range(len(rna_true_structure)):
        if rna_true_structure[i] == "(":
            first_pair_element.append(i)
        if rna_true_structure[-i-1] == ")":
            second_pair_element.append(len(rna_true_structure)+(-1-i))
    pair_index = list(zip(first_pair_element, second_pair_element))
    return pair_index


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

if __name__ == '__main__':
    structure = "....((((......))))...(((())))."
    print(convert_input2(structure))