# TCRscore
##### TCRscore is a deep learning-based model to predict the immunogenicity of neoantigen peptides by integrating the intrinsic features of TCRs, neoantigen peptides, HLA class I molecules, and the HLA evolutionary divergence metric (HED). The TCRscore model comprises five modules: the protein sequence coding, TCR representation, pMHC feature extraction, TCR–pMHC, and HED modules. The protein sequence coding module numerically encoded HLA sequences, peptide sequences, and TCR CDR3 regions using amino acid Z-descriptors. Considering that the fourth and fifth descriptors are not clearly derivable, we utilized the first three Z-descriptors, which correspond to lipophilicity (z1 or A), volume (z2 or B), and polarity (z3 or C). 

## I. Environmental installation
### Method 1：  
#### Install Anaconda and Create a conda virtual environment
##### conda create -n tcr python=3.7

#### Configure the environment - Download the source:：
##### conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2/
##### conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
##### conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
##### conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/
##### conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
##### conda config --set show_channel_urls yes

#### cudnn and cudatoolkit are required only if installing the GPU-enabled version of TensorFlow
##### conda install cudatoolkit=11.1.3
##### conda install cudnn=8.2.1 
##### conda install tensorflow-gpu=2.6

#### Runtime Environment Configuration
##### conda install python=3.7
##### conda install tensorflow=2.6
##### conda install numpy=1.19.5
##### conda install keras=2.6
##### conda install pandas=1.3.5
##### conda install scikit-learn=1.0.2
##### conda install scipy=1.7.3
##### conda install imbalanced-learn=0.9.0
##### conda install matplotlib


### Method 2：  
##### pip intsall tensorflow==2.6
##### pip install numpy
##### pip install keras==2.6
##### pip install pandas==1.3.5
##### pip install scikit-learn==1.0.2
##### pip install scipy==1.7.3
##### pip install imbalanced-learn

## II. Code running
### Step1: Setting the path
###### cd E:\TCRscore\TCRscore   

### Step2: predictTCRratio.py
#### Function: A Python script for predicting TCR ratios for antigen-HLA combinations.
#### This tool predicts TCR ratios using pre-trained encoders for TCR and pMHC, along with a trained classifier model. It takes input data containing TCR, HLA, and antigen information, and outputs prediction results.
##### Usage: python predictTCRratio.py -i <input> -tcr <tcr_encoder> -pmhc <pmhc_encoder> -model <model> -o <output>
#### Parameters:
##### -i: Path to input CSV file(required columns: hla, antigen, sample, optional:tcr).
##### -tcr: Path to the pre-trained TCR encoder (default:./Model/tcr_encoder.h5).
##### -pmhc:Path to the pre-trained pMHC encoder (default:./Model/hla_encoder.h5).
##### -model:Path to the trained classifier model (default:./Model/model.h5).
##### -o:Output directory to save prediction results (default:./output1/pmhc_counts.csv).
##### Example:  python predictTCRratio.py -i data_test.csv -o output1


### Step3:predict_hlahed.py
#### Function: A Python script for predicting HLA epitope dissimilarity using Grantham matrix and protein sequences.
##### Usage: python predict_hlahed.py -d <grantham_matrix> -f <protein_fasta> -i <input> -o <output>
#### Parameters:
##### -d, --grantham_matrix file Path (default: ./hlahed/grantham_matrix.txt)
##### -f, --HLA protein_fasta	Path (default: ./hlahed/ABC_prot.fa)
##### -i, --input Path to input file
##### -o, --output Path to output file for prediction results (default:./output1/testResult.txt)
##### Example:  python predict_hlahed.py -d ./hlahed/grantham_matrix.txt -f ./hlahed/ABC_prot.fa -i test.txt -o ./output1/testResult.txt

### Step4:predict_TCRscore.py
#### Function:A Python script that combines TCRratio and HLA epitope dissimilarity (HLAHED) predictions to calculate comprehensive TCR scores.
##### Usage:  python predictTCRscore.py -i <input> -o <output> [options]
#### Parameters:
##### -i, --input	   Path to input CSV file (same format as TCRratio)
##### -o, --output  Output directory for results(default: processed_data.csv)
#### Options
##### -tcr , --tcrratio	: Path to TCRratio predictions (default:"./output1/pmhc_counts.csv")
##### -hed, --hlahed	Path to HLAHED predictions (default: "./output1/testResult.txt")
##### Example: python predictTCRscore.py -i data_test.csv -o output1


