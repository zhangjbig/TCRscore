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
##### conda install biopython



### Method 2：  
##### pip intsall tensorflow==2.6
##### pip install numpy
##### pip install keras==2.6
##### pip install pandas==1.3.5
##### pip install scikit-learn==1.0.2
##### pip install scipy==1.7.3
##### pip install imbalanced-learn
##### pip install biopython







## II. Code running
### Step1: Setting the path
###### cd E:\TCRscore\TCRscore   


### Step2: predictTCRscore.py
#### Function: A Python script for predicting TCR score for antigen-HLA combinations.
#### This tool predicts TCR score using pre-trained encoders for TCR and pMHC, along with a trained classifier model. It takes input data containing TCR, HLA, and antigen information, and outputs prediction results.
##### Usage:  python predictTCRscore.py -i <input> -ihla <input2> -o <output> [options]
#### Parameters:
##### -i, --inputdata, Specify the csv file with standard format you want to predict. (required columns: hla, antigen, sample, optional: tcr).
##### -ihla, --inputdatahla, Sample HLA file with standard format you want to predict.(txt file,columns: Sample,A1,A2,B1,B2,C1,C2)
##### -o, --output, help=The output folder.
##### -tcr, --tcr_encoder, Path to the pre-trained TCR encoder (default:./Model/tcr_encoder.h5).
##### -pmhc,--pmhc_encoder, Path to the pre-trained pMHC encoder (default:./Model/hla_encoder.h5).
##### -model,--model,Path to the trained classifier model.(default='./Model/model.h5')
##### -tcrseq, --tcr_allseq, default='./tcr_data/tcr.csv', help=all tcr sequence csv.
##### -tcrfea, --tcr_allseqfea, default='./tcr_data/tcr.npy', help=all tcr sequence features npy.
##### -hlahed_d,--hlahed_d, grantham_matrix file Path for HLA Hed.(default: ./hlahed/grantham_matrix.txt)
##### -hlahed_f,--hlahed_f, HLA protein_fasta Path for HLA Hed.(default: ./hlahed/ABC_prot.fa)

##### Example:  python predictTCRscore.py -i data_test.csv -ihla Sample_HLA.txt -o output
