import numpy as numpy
import math as math
import scoring as sc
import formating as fo
import sys
import pandas as pd


if __name__ == '__main__':
    #bps = main.rna_base_pairing(sequence)
    #print(bps)
    #print(fo.encode_output(bps, sequence))
    data_nuss = pd.read_excel("D:/Courses/CS5112_Algorithm/data/results_correct_nuss.xlsx")
    data_zuker = pd.read_excel("D:/Courses/CS5112_Algorithm/data/results_correct_zuker.xlsx")
    data_seqfold = pd.read_excel("D:/Courses/CS5112_Algorithm/data/results_correct_seqfold.xlsx")
    data_output = []
    pd.DataFrame(data_output)
    data_output = data_nuss[["RNA_name", "RNA_length", "RNA_sequence", "RNA_structure", "RNA_true_base_pairing"]]
    data_output["RNA_true_base_pairing"] = data_output["RNA_structure"].apply(lambda x: fo.convert_input2(x))

    data_output["Nuss_predicted_structure"] = data_nuss["RNA_predicted_structure"]
    data_output["Nuss_predicted_base_pairing"] = data_output["Nuss_predicted_structure"].apply(
        lambda x: fo.convert_input2(x))
    data_output["Nuss_distances"] = data_output.apply(
        lambda x: sc.calculate_distances(x.RNA_true_base_pairing, x.Nuss_predicted_base_pairing),
        axis=1)
    data_output["Nuss_BP_score"] = data_output.apply(lambda x: sc.calculate_RBP_score(0, x.Nuss_distances), axis=1)

    data_output["Zuker_predicted_structure"] = data_zuker["RNA_predicted_structure"]
    data_output["Zuker_predicted_base_pairing"] = data_output["Zuker_predicted_structure"].apply(
        lambda x: fo.convert_input2(x))
    data_output["Zuker_distances"] = data_output.apply(
        lambda x: sc.calculate_distances(x.RNA_true_base_pairing, x.Zuker_predicted_base_pairing),
        axis=1)
    data_output["Zuker_BP_score"] = data_output.apply(lambda x: sc.calculate_RBP_score(0, x.Zuker_distances), axis=1)

    data_output["Seqfold_predicted_structure"] = data_seqfold["RNA_predicted_structure"]
    data_output["Seqfold_predicted_base_pairing"] = data_output["Seqfold_predicted_structure"].apply(
        lambda x: fo.convert_input2(x))
    data_output["Seqfold_distances"] = data_output.apply(
        lambda x: sc.calculate_distances(x.RNA_true_base_pairing, x.Seqfold_predicted_base_pairing),
        axis=1)
    data_output["Seqfold_BP_score"] = data_output.apply(lambda x: sc.calculate_RBP_score(0, x.Seqfold_distances), axis=1)
    # data_output["Nuss_predicted_base_pairing_number"] = data_output["Nuss_predicted_base_pairing"].apply(lambda x: x.count('('))
    # data_output["Zuker_predicted_base_pairing_number"] = data_output["Zuker_predicted_base_pairing"].apply(lambda x: x.count('('))
    # data_output["Seqfold_predicted_base_pairing_number"] = data_output["Seqfold_predicted_base_pairing"].apply(lambda x: x.count('(')
    data_output.to_excel("D:/Courses/CS5112_Algorithm/data/BP_scores_full.xlsx", sheet_name="BP_scores")

