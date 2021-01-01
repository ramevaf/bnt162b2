
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

def compose_sequence_series(sequence: str, *, as_codons=None) -> pd.Series:
    chunk_len = 3 if as_codons else 1
    chunks = (sequence[0+i:chunk_len+i] for i in range(0, len(sequence), chunk_len))
    return pd.Series(chunks)

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

def print_sequence_side_by_side(seqVirus: str, seqVaccine: str, seqOptimized: str):
    f = open("side-by-side_optimized.csv",'w')

    codonsVirus = compose_sequence_series(seqVirus, as_codons=True)
    codonsVaccine = compose_sequence_series(seqVaccine, as_codons=True)
    codonsOptimized = compose_sequence_series(seqOptimized, as_codons=True)
    codonsVirus.rename("Virus")
    codonsVaccine.rename("Vaccine")
    codonsOptimized.rename("Optimized")
    
    out = pd.concat([codonsVirus, codonsVaccine, codonsOptimized], axis=1)
    out.to_csv(f, sep=";")
    f.close()

if __name__ == '__main__':
    seqVirus, seqVaccine = load_sequence()

    seqOptimized = optimize_virus_sequence(sequence=seqVirus)
    print_sequence_side_by_side(seqVirus, seqVaccine, seqOptimized)
    

    score = compute_match(seqVaccine, seqOptimized)
    print(f"{score:.2%}")