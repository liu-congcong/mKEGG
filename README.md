# mKEGG

KEGG for metagenome.

## Dependencies

* Numpy
* Samtools
* Hmmer
* Prodigal

## Installation

### Download and install mKEGG

```shell
# You may need to install pip3 before. #
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3 get-pip.py

# You may need to install or upgrade setuptools and wheel using pip3 before. #
pip3 install --upgrade setuptools wheel

# Install mKEGG #
pip3 install https://github.com/liu-congcong/mKEGG/releases/download/v1.0.0/mkegg.1.0.0.tar.gz
```

## Usage

### Preparations

**Before running mKEGG, you may need to prepare some files by yourself.**

```shell
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
gunzip ko_list.gz
echo "export KO_LIST=/path/ko_list" >> ~/.bashrc
```

```shell
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
tar xvf profiles.tar.gz
cd profiles
cat K*.hmm > ../k.hmm
rm -rf profiles
echo "export K_HMM=/path/k.hmm" >> ~/.bashrc
```

```shell
wget -O ko00001.json "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json"
echo "export KO_JSON=/path/ko00001.json" >> ~/.bashrc
```

**Make sure that samtools (Samtools), hmmsearch (Hmmer) and prodigal (Prodigal) has been added to the environment variables. **

### Run mKEGG

#### Map K to metagenome assembly

```shell
mkegg mapping -threads CPUs -assembly ASSEMBLY.FASTA -output ASSEMBLY.KMAP
```

#### Estimate abundance of K assignments

```shell
mkegg abundance -threads CPUs -bam SORTED&INDEXED.BAM -kmap ASSEMBLY.KMAP -output ASSEMBLY.AKMAP
```

#### Estimate abundance of Pathway

```shell
mkegg pathway -kmap ASSEMBLY.AKMAP -clusters [CLSUTER1.FASTA CLUSTER2.FASTA ... | ASSEMBLY.FASTA] -output file PATHWAY.MATRIX
```
