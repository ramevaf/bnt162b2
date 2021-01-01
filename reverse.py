
import pandas as pd
from dnachisel import *
import python_codon_tables

# codon table
codon_usage_table = python_codon_tables.get_codons_table("h_sapiens_9606")

def load_sequence():
    rawInput = pd.read_csv("side-by-side.csv")
    rawInput.head()
    seqVirus = "".join(rawInput.codonOrig.tolist())
    seqVaccine = "".join(rawInput.codonVaccine.tolist())
    # proline substitution
    seqVirus.replace("CCTCCT", "AAAGTT")

    return seqVirus, seqVaccine

def compute_match(one, two):
    num_mismatches = 0
    for base1, base2 in zip(one, two):
        if base1 != base2:
            num_mismatches = num_mismatches + 1
    return (1 - float(num_mismatches) / float(len(one)))

def optimize_virus_sequence(sequence: str) -> str:
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            EnforceGCContent(mini=0.4, maxi=0.7, window=50),
            EnforceTranslation()
        ],
        objectives=[
            CodonOptimize(codon_usage_table=codon_usage_table)
        ]
    )
    problem.resolve_constraints()
    problem.optimize()
    return problem.sequence



if __name__ == '__main__':
    seqVirus, seqVaccine = load_sequence()

    optimized_sequence = optimize_virus_sequence(sequence=seqVirus)
    
    score = compute_match(seqVaccine, optimized_sequence)
    print(f"{score:.2%}")