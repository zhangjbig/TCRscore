import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser




parser = ArgumentParser(description="Specifying Input Parameters")
parser.add_argument("-i", "--input", required=True, help="Path to the input data file (data_test.csv).")
parser.add_argument("-tcr","--tcrratio", default="./output1/pmhc_counts.csv", help="tcrratio result")
parser.add_argument("-hed","--hlahed",default="./output1/testResult.txt", help="hla hed result")
parser.add_argument("-o","--output",required=True, help="Output file.")
args = parser.parse_args()


current_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_path)
####  Input file
datainput = pd.read_csv(args.input)
HLA_hed = pd.read_csv(args.hlahed, sep="\t")
TCRratio = pd.read_csv(args.tcrratio)

### Data process
datainput['sampleHLAHed'] = datainput['sample'].apply(lambda x: HLA_hed.loc[HLA_hed['Sample'] == x, 'Mean_HE'].values[0] if x in HLA_hed['Sample'].values else np.nan)
datainput['HLAtype'] = datainput['hla'].apply(lambda x: x.split('*')[0])
def calculate_hlaheduse(row):
    if row['sample'] in HLA_hed['Sample'].values:
        hla_type = row['HLAtype']
        return HLA_hed.loc[HLA_hed['Sample'] == row['sample'], f'HED_{hla_type}'].values[0]
    else:
        return np.nan

datainput['HLAheduse'] = datainput.apply(calculate_hlaheduse, axis=1)
datainput['pmhc'] = datainput['hla'] + "_" + datainput['antigen']
datainput['tcrratio'] = datainput['pmhc'].apply(lambda x: TCRratio.loc[TCRratio['pmhc'] == x, 'count'].values[0] if x in TCRratio['pmhc'].values else np.nan)
datainput['score'] = datainput['tcrratio'] * (1 + 0.5 * np.log10((1 + datainput['HLAheduse']) * (1 + datainput['sampleHLAHed'])))
#print(datainput.head())

columns_to_output = ['sample','hla','antigen','pmhc','score']
output_file = os.path.join(args.output, "processed_data.csv")
datainput[columns_to_output].to_csv(output_file, index=False)


