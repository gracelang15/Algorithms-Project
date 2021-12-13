import RNA
import pandas as pd


def cal_tree(rnadata, struct):
    '''

    :param rnadata: the input dataframe
    :param struct: a string, the name of the column containing the dot-bracket sequences that you want to convert to tree
    :return: an appended dataframe containing converted trees
    '''
    rnadata["HIT_tree"] = rna_data[struct].apply(lambda x: RNA.db_to_tree_string(x, 1))
    rnadata["Short_tree"] = rna_data[struct].apply(lambda x: RNA.db_to_tree_string(x, 2))
    rnadata["Full_tree"] = rna_data[struct].apply(lambda x: RNA.db_to_tree_string(x, 3))
    rnadata["Extended_tree"] = rna_data[struct].apply(lambda x: RNA.db_to_tree_string(x, 4))
    rnadata["Weighted_tree"] = rna_data[struct].apply(lambda x: RNA.db_to_tree_string(x, 5))
    rnadata["Expanded_tree"] = rna_data[struct].apply(lambda x: RNA.db_to_tree_string(x, 6))
    return rnadata

if __name__ == '__main__':
    rna_data = pd.read_excel("~/RNAData_10-410.xlsx")
    rna_data = cal_tree(rna_data, "RNA_structure")
    rna_data.to_excel("~/RNAData_10-410_tree.xlsx", sheet_name="RNAData")

    rna_data = pd.read_excel("~/results_nuss_full.xlsx")
    rna_data = cal_tree(rna_data, "RNA_predicted_structure")
    rna_data.to_excel("~/results_nuss_full_tree.xlsx", sheet_name="nussinov")

    rna_data = pd.read_excel("~/results_zuker_full.xlsx")
    rna_data = cal_tree(rna_data, "RNA_predicted_structure")
    rna_data.to_excel("~/results_zuker_full_tree.xlsx", sheet_name="zuker")
