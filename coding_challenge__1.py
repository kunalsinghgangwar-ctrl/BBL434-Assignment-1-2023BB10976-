#!/usr/bin/env python3

import argparse

def read_fasta(file_path):
    seq = []
    with open(file_path) as fh:
        for ln in fh:
            if not ln.startswith(">"):
                seq.append(ln.strip().upper())
    return "".join(seq)

def find_ori_by_gc(seq):
    score = 0
    values = []
    for nt in seq:
        if nt == "G":
            score += 1
        elif nt == "C":
            score -= 1
        values.append(score)
    return values.index(min(values))

def read_design(file_path):
    cut_sites = []
    resistances = []
    with open(file_path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln:
                continue
            key, val = [x.strip() for x in ln.split(",")]
            if "site" in key.lower():
                cut_sites.append(val)
            else:
                resistances.append(val)
    return cut_sites, resistances

def read_marker_db(file_path):
    data = {}
    with open(file_path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or "," not in ln:
                continue
            k, v = ln.split(",", 1)
            data[k.strip()] = v.strip()
    return data

RESTRICTIONS = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "SalI": "GTCGAC",
    "XbaI": "TCTAGA",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG",
    "NotI": "GCGGCCGC"
}

ANTIBIOTIC_GENES = {
    "Ampicillin": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTG",
    "Kanamycin": "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCC",
    "Chloramphenicol": "ATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCC",
    "Blue_White_Selection": "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGAC"
}

SCAR = "TACTAGAG"

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--markers", required=True)
    parser.add_argument("--output", default="Output.fa")
    args = parser.parse_args()

    genome_seq = read_fasta(args.input)
    ori_index = find_ori_by_gc(genome_seq)
    ori_block = genome_seq[ori_index-300:ori_index+300]

    sites, marker_list = read_design(args.design)
    marker_table = read_marker_db(args.markers)

    plasmid_seq = ori_block

    for m in marker_list:
        if m in ANTIBIOTIC_GENES:
            plasmid_seq += SCAR + ANTIBIOTIC_GENES[m]
        else:
            print(f"WARNING: marker {m} missing")

    for e in sites:
        if e in RESTRICTIONS:
            plasmid_seq += SCAR + RESTRICTIONS[e]
        else:
            print(f"WARNING: enzyme {e} missing")

    with open(args.output, "w") as out:
        out.write(">Universal_Plasmid\n")
        for i in range(0, len(plasmid_seq), 70):
            out.write(plasmid_seq[i:i+70] + "\n")

    print("Plasmid generated")
    print("ORI index:", ori_index)
    print("Total size:", len(plasmid_seq))

if __name__ == "__main__":
    run()
