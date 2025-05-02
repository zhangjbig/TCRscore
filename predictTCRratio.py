from sklearn import preprocessing
import pandas as pd
from keras.models import Model, load_model
import numpy as np
import os
from argparse import ArgumentParser


# Args parser
parser = ArgumentParser(description="Specifying Input Parameters")
parser.add_argument("-i", "--inputdata", help="Specify the file with standard format you want to predict")
parser.add_argument("-tcr", "--tcr_encoder", default='./Model/tcr_encoder.h5', help="Specify the encode tcr model")
parser.add_argument("-pmhc", "--pmhc_encoder", default='./Model/hla_encoder.h5', help="Specify the encode pmhc model")
parser.add_argument("-model", "--model", default='./Model/model.h5', help="Specify the model from deep learning")
parser.add_argument("-o", "--output", help="The output folder")
parser.add_argument("-tcrseq", "--tcr_allseq", default='./tcr_data/tcr.csv', help="all tcr sequence")
parser.add_argument("-tcrfea", "--tcr_allseqfea", default='./tcr_data/tcr.npy', help="all tcr sequence features")
args = parser.parse_args()

Z = pd.read_csv('Descriptorfiles/Z-descriptors.csv')
des_aa = Z['Amino acid'].tolist()
hla_antigen_encoder_path = args.pmhc_encoder
tcr_encoder_path = args.tcr_encoder
final_model_path = args.model
tcr_seq = pd.read_csv(args.tcr_allseq)
tcr_array = np.load(args.tcr_allseqfea)

hla_db_dir = 'hla_library'
HLA_ABC = [hla_db_dir + '/A_prot.fasta', hla_db_dir + '/B_prot.fasta', hla_db_dir + '/C_prot.fasta',
           hla_db_dir + '/E_prot.fasta']
HLA_seq_lib = {}
for one_class in HLA_ABC:
    prot = open(one_class)
    # pseudo_seq from netMHCpan:https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000796
    pseudo_seq_pos = [7, 9, 24, 45, 59, 62, 63, 66, 67, 79, 70, 73, 74, 76, 77, 80, 81, 84, 95, 97, 99, 114, 116, 118,
                      143, 147, 150, 152, 156, 158, 159, 163, 167, 171]
    # write HLA sequences into a library
    # class I alles
    name = ''
    sequence = ''
    for line in prot:
        if len(name) != 0:
            if line.startswith('>HLA'):
                pseudo = ''
                for i in range(0, 33):
                    if len(sequence) > pseudo_seq_pos[i]:
                        pseudo = pseudo + sequence[pseudo_seq_pos[i]]
                HLA_seq_lib[name] = pseudo
                name = line.split(' ')[1]
                sequence = ''
            else:
                sequence = sequence + line.strip()
        else:
            name = line.split(' ')[1]


def hla_sequence(HLA_names):
    hla_list = []
    print('Start getting HLA sequence!')
    for index, HLA_name in enumerate(HLA_names):
        print("iterate", index)
        if HLA_name not in HLA_seq_lib.keys():
            try:
                HLA_name = [hla_allele for hla_allele in HLA_seq_lib.keys() if hla_allele.startswith(str(HLA_name))][0]

                if HLA_name not in HLA_seq_lib.keys():
                    print('Not proper HLA allele:' + HLA_name)
                HLA_sequence = HLA_seq_lib[HLA_name]
                hla_list.append(HLA_sequence)
                primary_file.loc[index, "hla_sequence"] = HLA_sequence
            except:
                HLA_sequence = None
                hla_list.append(HLA_sequence)
                primary_file.loc[index, "hla_sequence"] = HLA_sequence
                print('FAIL')
        else:
            HLA_sequence = HLA_seq_lib[HLA_name]
            hla_list.append(HLA_sequence)
            primary_file.loc[index,"hla_sequence"] = HLA_sequence
 #   file_path = save_path + '/predict_res.csv'
 #   primary_file.to_csv(file_path)
    return primary_file


def get_descriptor(protein, maxlen):
    aas = list(protein)
    proteinDES = []
    n = 3
    DES = Z

    if len(aas) <= maxlen:
        for aa in aas:
            if aa not in des_aa:
                return (0, None)

            aaDES = np.array(DES[DES['Amino acid'] == aa]).tolist()[0]
            aaDES = aaDES[1:]
            proteinDES.append(aaDES)

        for i in range(maxlen - len(aas)):
            proteinDES.append([0] * n)

        non_normalized = np.asarray(proteinDES)
        normalized = preprocessing.MaxAbsScaler().fit_transform(non_normalized)
        return (1, normalized)


def get_descriptor_hla(protein, maxlen):
    aas = list(protein)
    proteinDES = []
    DES = Z

    if len(aas) <= maxlen:
        for aa in aas:
            if aa not in des_aa:
                return (0, None)

            aaDES = np.array(DES[DES['Amino acid'] == aa]).tolist()[0]
            aaDES = aaDES[1:]
            proteinDES.append(aaDES)

        non_normalized = np.asarray(proteinDES)
        normalized = preprocessing.MaxAbsScaler().fit_transform(non_normalized)
        return (1, normalized)

    elif len(aas) > maxlen:
        print('exceed maxlen!')
        return (0, None)


def get_descriptor_antigen(protein, maxlen):
    aas = list(protein)
    proteinDES = []
    n = 3
    DES = Z

    if len(aas) <= maxlen:
        for aa in aas:
            if aa not in des_aa:
                return (0, None)

            aaDES = np.array(DES[DES['Amino acid'] == aa]).tolist()[0]
            aaDES = aaDES[1:]
            proteinDES.append(aaDES)

        zeros = [[0] * n for i in range(maxlen - len(aas))]
        proteinDES = proteinDES[:len(aas) // 2] + zeros + proteinDES[len(aas) // 2:]

        non_normalized = np.asarray(proteinDES)
        normalized = preprocessing.MaxAbsScaler().fit_transform(non_normalized)
        return (1, normalized)

    elif len(aas) > maxlen:
        print('exceed maxlen!')
        return (0, None)


def encode(file):
    n = 3
    pos = 0
    maxlen_hla = 33
    maxlen_antigen = 15
    maxlen_tcr = 50
    error_num = 0

    hlas = file['hla_sequence'].tolist()
    antigens = file['antigen'].tolist()
    proteins = file['tcr'].tolist()

    tcr_array = np.zeros((len(proteins), maxlen_tcr, n), dtype=np.float64)
    hla_array = np.zeros((len(proteins), maxlen_hla, n), dtype=np.float64)
    antigen_array = np.zeros((len(proteins), maxlen_antigen, n), dtype=np.float64)

    print('Start Encoding!')
    for index, protein in enumerate(proteins):
        hla = hlas[index]
        antigen = antigens[index]
        tcr = proteins[index]

        res_hla, hlaDES = get_descriptor_hla(hla, maxlen_hla)
        res_antigen, antigenDES = get_descriptor_antigen(antigen, maxlen_antigen)
        res_tcr, tcrDES = get_descriptor(tcr, maxlen_tcr)

        if res_hla and res_antigen and res_tcr:
            hla_array[pos] = hlaDES.reshape(1, maxlen_hla, n)
            antigen_array[pos] = antigenDES.reshape(1, maxlen_antigen, n)
            tcr_array[pos] = tcrDES.reshape(1, maxlen_tcr, n)
            pos += 1
        else:
            file.loc[index, "res"] = "FAIL"
            print(file.loc[index, "res"])
            error_num += 1
        # if index%10000 == 0 and index != 0:
        #     name_hla = 'data/training/res_hla_n' + str(index // 20000)
        #     name_antigen = 'data/training/res_antigen_n'+ str(index//20000)
        #     name_tcr = 'data/training/res_tcr_n' + str(index // 20000)
        #
        #     np.save(name_tcr, tcr_array)
        #     np.save(name_hla, hla_array)
        #     np.save(name_antigen, antigen_array)
        print("iterate", index)
    # antigen_hla_file.to_csv('data/training/post_descriptor_file_n.csv')
    hla_array_res = hla_array[:(len(proteins) - error_num)]
    antigen_array_res = antigen_array[:(len(proteins) - error_num)]
    tcr_array_res = tcr_array[:(len(proteins) - error_num)]
    print('Conversion Finished!')
    return hla_array_res, antigen_array_res, tcr_array_res


def encode_notcr(file,tcr_array):
    n = 3
    pos = 0
    maxlen_hla = 33
    maxlen_antigen = 15
    error_num = 0
    
    hlas = file['hla_sequence'].tolist()  
    antigens = file['antigen'].tolist()
    
    hla_array = np.zeros((len(file), maxlen_hla, n), dtype=np.float64)  
    antigen_array = np.zeros((len(file), maxlen_antigen, n), dtype=np.float64)
    
    print('Start Encoding HLA-Antigen!')
    for index, row in file.iterrows():
        hla = hlas[index]
        antigen = antigens[index]
            
        res_hla, hlaDES = get_descriptor_hla(hla, maxlen_hla)
        res_antigen, antigenDES = get_descriptor_antigen(antigen, maxlen_antigen)
            
        if res_hla and res_antigen:
            hla_array[pos] = hlaDES.reshape(1, maxlen_hla, n)
            antigen_array[pos] = antigenDES.reshape(1, maxlen_antigen, n)
            pos += 1
        else:
            file.loc[index, "res"] = "FAIL"
            print(file.loc[index, "res"])
            error_num += 1
            print("iterate", index)
    hla_array_res = hla_array[:(len(file) - error_num)]
    antigen_array_res = antigen_array[:(len(file) - error_num)]
    # Repeat hla_array and antigen_array for each TCR in tcr_array
    hla_array_res = np.repeat(hla_array_res, tcr_array.shape[0], axis=0)
    antigen_array_res = np.repeat(antigen_array_res, tcr_array.shape[0], axis=0)
    tcr_array_res = np.tile(tcr_array, (hla_array_res.shape[0] // tcr_array.shape[0], 1, 1)) 
    # Expanded file 
    expanded_file = file.loc[file.index.repeat(tcr_array.shape[0])].reset_index(drop=True)
    expanded_file['tcr'] = np.repeat(tcr_seq['tcr'].values, len(file))
    print('Conversion Finished!')
    return hla_array_res, antigen_array_res, tcr_array_res,expanded_file



def extract_features(hla_array, antigen_array, tcr_array):
    hla_antigen_encoder = load_model(hla_antigen_encoder_path)
    pmhc_features = hla_antigen_encoder.predict([hla_array, antigen_array])

    tcr_encoder = load_model(tcr_encoder_path)
    tcr_features = tcr_encoder.predict(tcr_array)
    tcr_features = tcr_features.reshape((tcr_features.shape[0], tcr_features.shape[1]))

    return pmhc_features, tcr_features


def get_model(pmhc, tcr, file):
    model = load_model(final_model_path)
    result = model.predict([tcr, pmhc]).tolist()
    pred_res = []
    for i in range(len(result)):
        #pred = round(result[i][0])
        pred = result[i][0]
        pred_res.append(pred)
    file.insert(loc=0, column='result', value=pred_res)
    # file.to_csv(save_path + 'predict_res.csv')
    return file


if __name__ == "__main__":
    save_path = args.output
   
    primary_file = pd.read_csv(args.inputdata)
    hla_names = primary_file['hla'].tolist()
    primary_file = hla_sequence(hla_names)
    primary_file = primary_file.dropna(subset=['hla_sequence'])


    hla_names = primary_file['hla'].tolist()
    file = hla_sequence(hla_names)

    if 'tcr' in primary_file.columns:
        hla_array, antigen_array, tcr_array = encode(file)
    else:
        hla_array, antigen_array, tcr_array,file = encode_notcr(file, tcr_array)
    print("hla_array shape:", hla_array.shape)
    print("antigen_array shape:", antigen_array.shape)
    print("tcr_array shape:", tcr_array.shape)


    pmhc_features, tcr_features = extract_features(hla_array, antigen_array, tcr_array)
    file = get_model(pmhc_features, tcr_features, file)
    result = pd.DataFrame(file, columns=['tcr', 'hla', 'antigen', 'result'])
 #   result.to_csv(save_path + '/predict_res.csv')



    # Filter rows where result > 0.99
    result['result'] = pd.to_numeric(result['result'], errors='coerce')
    filtered_result = result[result['result'] > 0.99].copy()
    filtered_result['pmhc'] = filtered_result['hla'] + '_' + filtered_result['antigen']
 #   filtered_result.to_csv(save_path + '/filtered_predict_res.csv')

    # Count the number of rows for each pmhc
    pmhc_counts = filtered_result['pmhc'].value_counts().reset_index()
    pmhc_counts.columns = ['pmhc', 'count']
    pmhc_counts['count'] = pmhc_counts['count'] / 67300
    pmhc_counts.to_csv(save_path + '/pmhc_counts.csv')
    print("Filtered results and pmhc counts have been saved to the output folder.")

  # Save the pmhc counts with sample information
  # pmhc_sample_info = filtered_result[['pmhc', 'sample']].drop_duplicates()
  # pmhc_counts = pmhc_counts.merge(pmhc_sample_info, on='pmhc', how='left')
  #  pmhc_counts.to_csv(save_path + '/pmhc_counts.csv')
  #  print("Filtered results and pmhc counts have been saved to the output folder.")
    pass