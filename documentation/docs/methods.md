# Methods

## Code

### Set up

Begin by cloning this repository

```
git clone https://github.com/ira-zibbu/picocell.git
```

This additional submodule is also required

```
cd picocell/src/
git submodule add https://github.com/ira-zibbu/EvolvePro
```

Create the necessary environments

```
cd EvolvePro
conda env create -f environment.yml
sh setup_plm.sh
```

### Training the regressor

The first part requires us to a build a regression model that maps the embeddings of the proteins to their actvitiy or fitness. We will use a deep mutational scanning dataset of cas12a proteins to train the regressor.

First, we need to generate a FASTA file of all possible single amino acid substitution variants of the wildtype cas12a sequence. 

```
conda activate evolvepro
./trial_run.py --file ../../run/01_wildtype/cas12f_WT.fasta --output ../../run/02_mutant_fasta/cas12f_aa_mutants.fasta
conda deactivate evolvepro
```

With all the mutants generated, we can now use a PLM to generated embeddings. We will use ESM1 (650M) model, and average the embeddings across all of the amino acids.

```
conda activate plm
python3 evolvepro/plm/esm/extract.py  esm1b_t33_650M_UR50S ../../run/02_mutant_fasta/cas12f_aa_mutants.fasta  ../../run/04_embeddings/casf21f_esm1b_t33_650M_UR50S --toks_per_batch 512 --include mean --concatenate_dir ../../run/04_embeddings
conda deactivate plm
```

This step can take a few hours on a personal computer.


