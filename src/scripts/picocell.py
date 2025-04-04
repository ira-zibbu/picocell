#!/usr/bin/env python3
"""
picocell: an evolution simulator!
"""

__author__ = "Ira Zibbu"
__version__ = "0.1.0"

""" imports """
import re
from Bio.Seq import Seq
from Bio import SeqIO
import random
from train_RF import load_DMS_data, load_embeddings
import copy
import joblib
import subprocess
import numpy as np
import os
import shutil

class cell:

    """ A class that represents a cell  """

    def __init__(self, id=0, dna_sequence="",aa_sequence="",fitness=0.0,mutation_rate=0.0):
        self.id = id # each cell has an integer ID number to keep track of it
        self.dna_sequence = dna_sequence # a string that represents the DNA sequence
        self.aa_sequence = aa_sequence # a string that represents the DNA sequence
        self.fitness = fitness # float that represents the fitness of the cell, based on its sequence
        self.mutation_rate = mutation_rate # a float that is the per-base per-generation mutation rate

    def __str__(self):
        return "Cell ID: {0}, with dna sequence {1}, aa sequence {2} fitness {3} and mutation rate {4}".format(self.id,self.dna_sequence,self.aa_sequence,self.fitness, self.mutation_rate)

def test_reproduction(fitness, kappa):

    """
    Tests to see if a cell will divide, based on its fitness value and the constant kappa
    args:
        fitness (float): float representing fitness of the cell
        kappa (float): integer constant to scale fitness
    returns:
        bool: Returns true with a probability fitness*kappa
    """
    reproduction_probability = fitness * kappa
    random_value = random.random()
    # if the random value is less than the reproduction probability, the cell divides
    return random_value < reproduction_probability

def mutate_sequence(wt_sequence, mutation_rate):

    """
    Accept a wt dna sequence, and a per base mutation rate. Randomly choose a base to mutate based on the mutation rate. Return sequence.
    args:
        wt_sequence (string): dna sequence for the gene
        mutation_rate (float): a mutation rate
    returns:
        mutant_sequence (string): mutated sequence. It will be the same as wt_sequence if no mutation occurs.
    """

    valid_bases = ['A','T','G','C']
    mutated_sequence = ""
    for base in wt_sequence:
        random_value = random.random()
        if random_value > mutation_rate:
            mutated_sequence+=base
            continue
        else:
            base_choices = [b for b in valid_bases if b != base] # list of possible bases to mutate to
            random_base = random.choice(base_choices) # pick a random base to replace it
            mutated_sequence+=random_base

    mutated_sequence = Seq(mutated_sequence)
    mutated_aa_sequence=mutated_sequence.translate(to_stop=True)

    return str(mutated_sequence), str(mutated_aa_sequence)

def fetch_fitness(wt_aa_sequence,aa_sequence,path_to_dms,path_to_rf_model):

    """ Accept an aa_sequence and return fitness score. If this a single aa variant that is present in the DMS data, read the value from the table. Else, use the trained random forest regressor to predict values
    args:
        wt_aa_sequence (string): amino acid sequence
        aa_sequence (string): amino acid sequence
        path_to_dms (string): path to the dms csv file
        path_to_rf_model (string): path to the rf model
    returns:
        fitness (float): fitness value corresponding to the
    """

    # load dms data as a dataframe
    df_dms = load_DMS_data(path_to_dms)

    # find number of aa that are mutated wrt the wildtype sequence. Generate strings representing the mutations eg" N2K"

    mutant_sites = [] # list to store mutants
    for index, aa in enumerate(wt_aa_sequence):
        if aa != aa_sequence[index]:
            mutation_variant = f"{aa}{index+1}{aa_sequence[index]}"
            mutant_sites.append(mutation_variant)

    print(mutant_sites)

    # if only one mutation is present wrt to wildtype sequence, attempt to find it in the dms data and fetch fitness score from them
    if len(mutant_sites) == 1:
        variant = mutant_sites[0]
        fitness = df_dms[df_dms["variant"] == variant]["avg_activity"].iloc[0]
        return fitness

    # else load rf model
    rf_model = None

    if len(mutant_sites) > 1:
        rf_model = joblib.load(path_to_rf_model)

    # extract embeddings for the new amino acid sequence

    temp_fasta_file="tmp/temp.fasta"

    with open(temp_fasta_file,'w') as file:
        file.write(">amino_acid_sequence\n")
        file.write(aa_sequence)

    temp_dir = "tmp/"
    plm_extraction_command=["python3", "src/EvolvePro/evolvepro/plm/esm/extract.py", "esm1b_t33_650M_UR50S", temp_fasta_file,  "tmp/casf21f_esm1b_t33_650M_UR50S", "--toks_per_batch", "512", "--include", "mean","--concatenate_dir", temp_dir]

    try:
        # Run the subprocess
        subprocess.run(plm_extraction_command, check=True)
        print(f"Subprocess extraction completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess: {e}")

    # pass embeddings to the rf model to predict fitness
    path_to_embeddings = "tmp/temp_esm1b_t33_650M_UR50S.csv"

    df_embeddings = load_embeddings(path_to_embeddings)

    fitness = rf_model.predict(df_embeddings) # use embeddings to predict fitness

    # clean up temp files
    os.remove("tmp/temp.fasta")
    os.remove(path_to_embeddings)
    shutil.rmtree("tmp/casf21f_esm1b_t33_650M_UR50S")

    return fitness[0]

def reproduce(cell_counter,parent,path_to_dms,path_to_rf_model):

    """
    Makes a copy of the parent cell to produce a child cell. The sequence of the child cell is mutated with a certain probability, and its fitness is updated with the random forest model
    args:
        parent (object): parent cell
    returns:
        child (object): child cell object produced by parent
    """

    # mutate the dna sequence with a certain proability
    mutated_sequence, mutated_aa_sequence = mutate_sequence(parent.dna_sequence, parent.mutation_rate)

    if len(mutated_aa_sequence) < len(parent.aa_sequence): # did we truncate the protein with a nonsense codon?
        print("Lethal mutation!")
        child = cell(cell_counter,mutated_sequence,mutated_aa_sequence,0, parent.mutation_rate)
        return child # this child cell is inviable and is killed

    if mutated_aa_sequence == parent.aa_sequence: # no mutation occured, and the child cell is a copy of the parent
        print("Unmutated child")
        child = copy.copy(parent)
        child.id = cell_counter
        return child

    if mutated_aa_sequence != parent.aa_sequence:
        print("mutant child")
        fitness = fetch_fitness(parent.aa_sequence,mutated_aa_sequence,path_to_dms,path_to_rf_model)
        child = cell(cell_counter,mutated_sequence,mutated_aa_sequence,fitness, parent.mutation_rate)
        return child

def main():

    path_to_fasta = "run/01_wildtype/cas12f_WT_dna.fasta"
    wt_dna = SeqIO.read(path_to_fasta,"fasta")
    wt_aa = "MIKVYRYEIVKPLDLDWKEFGTILRQLQQETRFALNKATQLAWEWMGFSSDYKDNHGEYPKSKDILGYTNVHGYAYHTIKTKAYRLNSGNLSQTIKRATDRFKAYQKEILRGDMSIPSYKRDIPLDLIKENISVNRMNHGDYIASLSLLSNPAKQEMNVKRKISVIIIVRGAGKTIMDRILSGEYQVSASQIIHDDRKNKWYLNISYDFEPQTRVLDLNKIMGIDLGVAVAVYMAFQHTPARYKLEGGEIENFRRQVESRRISMLRQGKYAGGARGGHGRDKRIKPIEQLRDKIANFRDTTNHRYSRYIVDMAIKEGCGTIQMEDLTNIRDIGSRFLQNWTYYDLQQKIIYKAEEAGIKVIKIDPQYTSQRCSECGNIDSGNRIGQAIFKCRACGYEANADYNAARNIAIPNIDKIIAESIK"
    kappa = 0.5
    list_of_parents=[]
    list_of_parents.append(cell(1,wt_dna.seq,wt_aa,1.0,0.001))

    path_to_dms = "run/03_dms_data/DMS_AsCas12f_preprocessed.xlsx"
    path_to_rf_model = "run/05_rf_model/rf_model.pkl"
    number_of_generations = 3
    cell_counter = 1 # keep a count of how many cells are present
    for generation in range(number_of_generations):
        print(f"generation count is {generation}")
        list_of_children = [] #make a fresh list of children cells born each generation
        for parent in list_of_parents:
            reproduction_bool = test_reproduction(parent.fitness,kappa) # trigger cell division with a probability proportional to fitness
            if reproduction_bool:
                print(f"Time for cell #{parent.id} to reproduce!")
                cell_counter+=1
                child = reproduce(cell_counter,parent, path_to_dms, path_to_rf_model)
                list_of_children.append(child)
            if not(reproduction_bool):
                print(f"cell #{parent.id} fails to reproduce")

        list_of_parents = list_of_parents+list_of_children # all children born this generation become children of the next generation

    for parent in list_of_parents:
        print(str(parent))

main()
