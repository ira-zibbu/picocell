# Code

### Set up

Begin by cloning this repository

```
git clone https://github.com/ira-zibbu/picocell.git
```

Create the picocell conda environment

```
cd picocell
conda env create -f picocell.yaml
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

The first part requires us to a build a regression model that maps the embeddings of the proteins to their actvitiy or fitness. We will use a deep mutational scanning dataset of cas12f protein to train the regressor.

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

Using these embeddings, and the deep mutational scanning dataset, we can train a random forest regressor to learn to map the embeddings to the activity of the protein (which in this case is also the fitness of the protein.

```
conda activate picocell
src/scripts/train_RF.py --dms run/03_dms_data/DMS_AsCas12f_preprocessed.xlsx --embeddings run/04_embeddings/cas12f_aa_mutants_esm1b_t33_650M_UR50S.csv --model_path run/05_rf_model/rf_model.pkl 
conda deactivate picocell
```

This step can take a few hours on a personal computer. This outputs a pickle object of the model in `run/05_rf_model`. 

### Running picocell

Picocell needs its own conda environment

```
conda activate picocell
```

Here is the general usage for picocell:

```
usage: picocell.py [-h] [--output OUTPUT] [--no_selection] [--drift DRIFT] [--generations GENERATIONS] [--mutation_rate MUTATION_RATE]

Picocell.py, a model of evolution

options:
  -h, --help            show this help message and exit
  --output OUTPUT       Path to output the results of the run
  --no_selection        Disable natural selection
  --drift DRIFT         Set genetic drift level (0-1).
  --generations GENERATIONS
                        Number of generations to run the simulation for
  --mutation_rate MUTATION_RATE
                        Per base per generation mutation rate (float)
```

The `output` is a csv file, with the following fields:

- Cell ID
- DNA sequence
- Amino acid sequence
- Fitness
- Mutation rate
- Generation (aka when was this cell born)
- Death (aka when did this cell die)

Here's is an example usage and output:

```
src/scripts/picocell.py --output run/06_output/run_2.csv --drift 0 --generations 20  --mutation_rate 0.001
```

This run has no genetic drift at play, goes for 20 generations and uses a fairly high mutation rate. Remember, cells divide exponentially, so make sure you limit the total number of generations.

```
================================
       _                    _ _ 
      (_)                  | | |
_ __  _  ___ ___   __ _ ___| | |
| '_ || |/ __/ _ | / __/ _ | | |
| |_) | | (_| (_) | (_|  __/ | |
| .__/|_||___|___/ |___|___|_|_|
| |                             
|_|                             
================================
Picocell is a PLM-based model of cellular evolution
Initial parameters (loading)...
Natural selection present? True
Genetic drift level? 0
Number of starting cells: 1
Information about starting cells:
Cell #1, with fitness 1.0 and mutation rate 0.01
Starting simulation...
generation count is 1
cell #1 fails to reproduce
generation count is 2
cell #1 fails to reproduce
generation count is 3
Time for cell #1 to reproduce!
mutant child with non-synonymous substitution
['F48Y', 'K96R', 'D100V', 'I109N', 'I166S', 'L337R']
Subprocess extraction completed successfully.
Size of embeddings is (1, 1280)
generation count is 4
cell #1 fails to reproduce
cell #2 fails to reproduce
generation count is 5
cell #1 fails to reproduce
cell #2 fails to reproduce
generation count is 6
cell #1 fails to reproduce
Time for cell #2 to reproduce!
Lethal mutation! This child cell is killed
generation count is 7
Time for cell #1 to reproduce!
mutant child with non-synonymous substitution
['A172S', 'N296D', 'R335H', 'A400S', 'I416S']
Subprocess extraction completed successfully.
Size of embeddings is (1, 1280)
cell #2 fails to reproduce
cell #3 fails to reproduce
```





