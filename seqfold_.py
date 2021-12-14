import seqfold as sf
import formating as fo
import pandas as pd
import numpy as numpy
import scoring as sc

'''
Installation of seqfold package:
pip install seqfold
(Available on Linux, Mac, and Windows)
https://github.com/Lattice-Automation/seqfold
'''


def seqfold_struct(seq):
    """
    :param seq: a string, the input RNA sequence
    :return: the container of all base pairs
    """
    structs = sf.fold(seq)
    db_output = []
    for struct in structs:
        db_output.append(struct.ij[0])
    return db_output


if __name__ == "__main__":
    rna_data = pd.read_excel("D:/Courses/CS5112_Algorithm/data/RNAData_10-410-for_seqfold.xlsx")
    # Use below line for testing
    # rna_data = pd.read_excel("D:/Courses/CS5112_Algorithm/data/testdata.xlsx")
    rna_data["RNA_true_base_pairing"] = rna_data["RNA_structure"].apply(lambda x: fo.convert_input(x))
    print('done rna true base pairing')
    rna_data["RNA_predicted_base_pairing"] = rna_data["RNA_sequence"].apply(lambda x: seqfold_struct(x))
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
    rna_data.to_excel("D:/Courses/CS5112_Algorithm/data/results_seqfold_full.xlsx", sheet_name="seqfold")
    # Use below line for testing
    # rna_data.to_excel("D:/Courses/CS5112_Algorithm/data/results_testdata.xlsx", sheet_name="seqfold")
