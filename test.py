import numpy as numpy
import main
import formating as fo
import pandas as pd

sequence = "GUCUACGGCCAGAAA"
if __name__ == '__main__':
    #bps = main.rna_base_pairing(sequence)
    #print(bps)
    #print(fo.encode_output(bps, sequence))
    rna_data = pd.read_excel("~/RNAData_10-410.xlsx")
