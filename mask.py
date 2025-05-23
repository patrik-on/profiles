#!/usr/bin/env python
from __future__ import annotations
import os
import sys
import csv
import re
import argparse
from collections import defaultdict
from typing import TextIO, Tuple, List
import numpy as np

FASTA_LABEL_SYMBOL = ">"
SIMPLE_ALPHABET = "ACGT"
FULL_ALPHABET = "ACGTN-"


def load_fasta(f: TextIO):
    label, buffer = None, []
    for line in f:
        if line.startswith(FASTA_LABEL_SYMBOL):
            if buffer:
                yield label, "".join(buffer)
            buffer = []
            label = line.strip()[1:]
        else:
            buffer.append(line.strip())
    if buffer:
        yield label, "".join(buffer)


def is_genomic_alphabet(s: str):
    return re.fullmatch(rf"[{SIMPLE_ALPHABET}]*", s) is not None


def load_virus_genome(filename: str):
    try:
        with open(filename, "r", encoding="utf-8") as f:
            full_genome = list(load_fasta(f))
    except Exception as e:
        print(f"Error loading genome '{filename}': {e}")
        return None
    if not full_genome:
        print(f"No sequence found in '{filename}'.")
        return None
    genome = full_genome[0][1].upper()
    if not is_genomic_alphabet(genome):
        print(f"Invalid characters in genome: {genome[:20]}...")
        return None
    return genome


def _load_variants_raw(f: TextIO, alphabet: str) -> Tuple[List[str], np.ndarray]:
    reader = csv.reader(f, delimiter='\t')
    raw_data = defaultdict(lambda: defaultdict(lambda: [0.0] * len(alphabet)))
    variants_order, seen, max_pos = [], set(), -1
    for row in reader:
        if len(row) < 4:
            continue
        try:
            var, pos, letter, cnt = row[0], int(row[1]), row[2].upper(), float(row[3])
        except ValueError:
            continue
        try:
            idx = alphabet.index(letter)
        except ValueError:
            continue
        if var not in seen:
            variants_order.append(var)
            seen.add(var)
        raw_data[var][pos][idx] = cnt
        max_pos = max(max_pos, pos)
    if not raw_data:
        raise ValueError("No valid variant data.")
    length = max_pos + 1
    table = np.zeros((len(variants_order), length, len(alphabet)))
    for i, var in enumerate(variants_order):
        for p, counts in raw_data[var].items():
            table[i, p, :] = counts
    return variants_order, table


def _positions_of_missing_letters(a: str, b: str) -> List[int]:
    bset = set(b)
    return [i for i, l in enumerate(a) if l not in bset]


def _remove_letters_from_variants(table: np.ndarray, positions: List[int]) -> np.ndarray:
    return np.delete(table, positions, axis=2)


def _normalized_variants(table: np.ndarray) -> np.ndarray:
    if table.size == 0:
        return table
    result = np.zeros_like(table)
    uniform = np.full(table.shape[2], 1.0 / table.shape[2])
    for v in range(table.shape[0]):
        coverage = table[v].sum(axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            normed = np.divide(table[v], coverage[:, None])
        zero_cov = coverage == 0
        normed[zero_cov] = uniform
        result[v] = normed
    return result


def load_variants(filename: str) -> Tuple[List[str], np.ndarray]:
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            names, table = _load_variants_raw(f, FULL_ALPHABET)
    except Exception as e:
        print(f"Error loading variants '{filename}': {e}")
        return None, None
    if table is None or table.size == 0:
        print("No variant data loaded.")
        return None, None
    to_remove = _positions_of_missing_letters(FULL_ALPHABET, SIMPLE_ALPHABET)
    filtered = _remove_letters_from_variants(table, to_remove)
    if filtered.shape[2] == 0:
        print("All bases removed by alphabet filter.")
        return names, None
    normalized = _normalized_variants(filtered)
    return names, normalized


def generate_uniformity_mask_vcf(variant_table, variant_names, genome_seq, output_vcf, chrom: str, epsilon=1e-9):
    if variant_table is None or variant_table.size == 0 or not genome_seq:
        print("Cannot generate mask.")
        return
    nv, gl, nb = variant_table.shape
    uniform_p = 1.0 / nb
    is_uniform = np.abs(variant_table - uniform_p) < epsilon
    all_uniform = is_uniform.all(axis=2).all(axis=0)
    mask_pos = set(np.nonzero(all_uniform)[0] + 1)
    with open(output_vcf, "w", encoding="utf-8", newline='') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FILTER=<ID=mask,Description="Uniform across all variants">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, base in enumerate(genome_seq):
            flt = "mask" if (i + 1) in mask_pos else "."
            f.write(f"{chrom}\t{i+1}\t.\t{base}\t.\t.\t{flt}\t.\n")


def parse_vcf_pos(line: str) -> int | None:
    if line.startswith("#"):
        return None
    parts = line.split("\t")
    try:
        pos = int(parts[1])
        return pos if pos > 0 else None
    except (IndexError, ValueError):
        return None


def load_mask(filename: str) -> dict[int, str]:
    data = {}
    with open(filename, encoding='utf-8') as f:
        for L in f:
            Ls = L.strip()
            if not Ls or Ls.startswith("#"):
                continue
            p = parse_vcf_pos(Ls)
            if p:
                data[p] = Ls
    return data


def process_first_block(lines: list[str], mask: dict[int, str]) -> list[str]:
    out, exp, filled, skip = [], 1, 0, 0
    for L in lines:
        p = parse_vcf_pos(L)
        if p is None or p < exp:
            skip += 1
            continue
        if p > exp:
            for x in range(exp, p):
                if x in mask:
                    out.append(mask[x])
                    filled += 1
            exp = p
        out.append(L)
        exp += 1
    print(f"Filled gaps: {filled}, skipped: {skip}")
    return out


def finalize_vcf(primary: str, mask: dict[int, str], output: str) -> bool:
    hdr, first, second, sec = [], [], [], False
    with open(primary, encoding='utf-8') as f:
        for L in f:
            if L.startswith("#"):
                hdr.append(L)
            else:
                if not sec and L.startswith("*"):
                    sec = True
                (second if sec else first).append(L.strip())
    out1 = process_first_block(first, mask)
    with open(output, "w", encoding='utf-8') as f:
        f.writelines(hdr)
        f.writelines(l + "\n" for l in out1 + second)
    return True


def main():
    parser = argparse.ArgumentParser(description="Generate and apply uniformity mask for VCF files.")
    parser.add_argument("-g", "--genome-file", default="genome.fa", help="Cesta k FASTA súboru s referenčným genomom.")
    parser.add_argument("-p", "--profiles-file", default="calculated_profiles.tsv", help="Cesta k TSV súboru s profilmi variantov.")
    parser.add_argument("--primary-mask-file", default="ont-short.masking.vcf", help="Cesta k pôvodnému maskovaciemu VCF súboru.")
    parser.add_argument("--uniform-mask-file", default="uniformity_mask.vcf", help="Výstupný VCF pre uniformitu.")
    parser.add_argument("-o", "--output-file", default="final.vcf", help="Cesta k finálnemu VCF súboru.")
    parser.add_argument("--chrom", default="MN908947.3", help="Názov chromozómu alebo kontigu.")
    args = parser.parse_args()

    # Check input files
    for path in (args.profiles_file, args.genome_file, args.primary_mask_file):
        if not os.path.exists(path):
            print(f"Missing file: {path}")
            sys.exit(1)

    genome = load_virus_genome(args.genome_file)
    names, table = load_variants(args.profiles_file)
    generate_uniformity_mask_vcf(table, names, genome, args.uniform_mask_file, args.chrom)
    mask_data = load_mask(args.uniform_mask_file)
    if finalize_vcf(args.primary_mask_file, mask_data, args.output_file):
        print(f"Final VCF written to {args.output_file}")


if __name__ == "__main__":
    main()
