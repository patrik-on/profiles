import sys
import csv
import yaml
import argparse


def load_reference_genome(file_path):
    sequence_lines = []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                continue
            sequence_lines.append(line.strip())
    return "".join(sequence_lines)


def load_barcode(file_path):
    barcodes = {}
    with open(file_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            clade = row[0]
            mutations = {
                header[i]: 1
                for i in range(1, len(row))
                if float(row[i]) == 1.0
            }
            barcodes[clade] = mutations
    return barcodes


def load_lineages(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data


def load_clade_counts(file_path):
    clade_counts = {}
    with open(file_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            clade = row["clade"]
            clade_counts[clade] = {
                "inclusive_count": int(row["inclusive_count"]),
                "exclusive_count": int(row["exclusive_count"])
            }
    return clade_counts


def build_parent_map(lineages):
    parent_map = {}
    for entry in lineages:
        clade = entry.get("name")
        parents = []
        p = entry.get("parent")
        if p:
            parents.extend(p if isinstance(p, list) else [p])
        rp = entry.get("recombinant_parents")
        if rp:
            if isinstance(rp, list):
                parents.extend(rp)
            else:
                parents.extend([x.strip() for x in rp.split(",") if x.strip()])
        if parents:
            parent_map[clade] = parents
    return parent_map


def build_children_map(parent_map):
    children_map = {}
    for child, parents in parent_map.items():
        for p in parents:
            children_map.setdefault(p, []).append(child)
    return children_map


def get_all_descendants(children_map, variant):
    descendants = set()
    stack = [variant]
    while stack:
        cur = stack.pop()
        for c in children_map.get(cur, []):
            if c not in descendants:
                descendants.add(c)
                stack.append(c)
    return descendants


def get_all_ancestors(parent_map, variant):
    ancestors = set()
    stack = [variant]
    while stack:
        cur = stack.pop()
        for p in parent_map.get(cur, []):
            if p not in ancestors:
                ancestors.add(p)
                stack.append(p)
    return ancestors


def parse_mutation_id(mutation_id):
    original_nuc = mutation_id[0]
    new_nuc = mutation_id[-1]
    position = mutation_id[1:-1]
    return position, original_nuc, new_nuc


def compute_weighted_profile(variant, barcodes, clade_counts, children_map):
    group = {variant} | get_all_descendants(children_map, variant)
    total_weighted = sum(
        clade_counts.get(v, {}).get("exclusive_count", 0)
        for v in group
    )
    mutation_weighted_sum = {}
    if variant in barcodes:
        for mut in barcodes[variant]:
            mutation_weighted_sum[mut] = total_weighted
    return group, total_weighted, mutation_weighted_sum


def write_combined_profiles(file_path, profiles, genome_length):
    with open(file_path, "w", encoding="utf-8") as f:
        for variant, (mutation_weighted_sum, total_weighted) in profiles.items():
            default_value = total_weighted / 4 if total_weighted else 0
            profile_by_pos = {}
            for mut, wt in mutation_weighted_sum.items():
                pos_str, orig, new = parse_mutation_id(mut)
                try:
                    pos = int(pos_str) - 1
                except ValueError:
                    continue
                profile_by_pos.setdefault(pos, {"A": 0, "C": 0, "G": 0, "T": 0})
                profile_by_pos[pos][new] += wt
            for pos in range(genome_length):
                if pos not in profile_by_pos:
                    profile_by_pos[pos] = {b: default_value for b in ["A", "C", "G", "T"]}
            for pos in range(genome_length):
                for b in ["A", "C", "G", "T"]:
                    f.write("\t".join([variant, str(pos), b, str(profile_by_pos[pos][b])]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Vypočíta vážené profily mutácií pre vybrané varianty a 'other' skupinu.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-b", "--barcode-file", default="usher_barcodes.csv",
                        help="Cesta k CSV súboru s barkódmi (mutáciami).")
    parser.add_argument("-l", "--lineage-file", default="lineages.yml.txt",
                        help="Cesta k YAML súboru s definíciami línií a rodičov.")
    parser.add_argument("-c", "--clades-file", default="clades.tsv",
                        help="Cesta k TSV súboru s počtami sekvencií pre clade_counts.")
    parser.add_argument("-g", "--genome-file", default="genome.fa",
                        help="Cesta k FASTA súboru s referenčným genómom.")
    parser.add_argument("-o", "--output-file", default="calculated_profiles.tsv",
                        help="Cesta k výstupnému súboru s kombinovanými profilmi (formát TSV).")
    parser.add_argument("-v", "--variants", nargs='+',
                        default=["B.1.617.2", "B.1.1.7"],
                        help="Zoznam variantov, pre ktoré sa majú počítať profily.")
    args = parser.parse_args()

    barcodes = load_barcode(args.barcode_file)
    lineage_data = load_lineages(args.lineage_file)
    clade_counts = load_clade_counts(args.clades_file)

    lineage_names = {entry["name"] for entry in lineage_data}

    missing = False
    for var in args.variants:
        if var not in barcodes:
            print(f"ERROR: Variant '{var}' not found in barcode file.")
            missing = True
        if var not in lineage_names:
            print(f"ERROR: Variant '{var}' not found in lineage file.")
            missing = True
        if var not in clade_counts:
            print(f"ERROR: Variant '{var}' not found in clades file.")
            missing = True
    if missing:
        sys.exit(1)

    print("--- Konfigurácia ---")
    print(f"Súbor barkódov:     {args.barcode_file}")
    print(f"Súbor línií:        {args.lineage_file}")
    print(f"Súbor clade_counts: {args.clades_file}")
    print(f"Súbor genómu:       {args.genome_file}")
    print(f"Výstupný súbor:     {args.output_file}")
    print(f"Vybrané varianty:   {', '.join(args.variants)}")
    print("--------------------")

    print("Načítavam dáta...")
    ref_genome = load_reference_genome(args.genome_file)
    GENOME_LENGTH = len(ref_genome)
    print(f"Načítaný referenčný genóm má dĺžku: {GENOME_LENGTH}")

    print(f"Načítaných {len(barcodes)} barkódov.")
    print("Načítané dáta o líniách.")
    print(f"Načítaných {len(clade_counts)} záznamov o počtoch.")

    parent_map = build_parent_map(lineage_data)
    children_map = build_children_map(parent_map)
    print(f"Vytvorená mapa rodičov s {len(parent_map)} záznamami.")

    profiles = {}
    for var in args.variants:
        _, total_wt, mut_sum = compute_weighted_profile(
            var, barcodes, clade_counts, children_map
        )
        profiles[var] = (mut_sum, total_wt)

    # build 'other'
    exclude_set = set(args.variants)
    for var in args.variants:
        exclude_set |= get_all_ancestors(parent_map, var)

    all_vars = set(clade_counts.keys())
    other = all_vars - exclude_set

    total_wt_o = 0
    mut_sum_o = {}
    for anc in other:
        group = {anc} | get_all_descendants(children_map, anc)
        wt = sum(clade_counts.get(v, {}).get("exclusive_count", 0)
                 for v in group)
        total_wt_o += wt
        if anc in barcodes:
            for mut in barcodes[anc]:
                mut_sum_o[mut] = mut_sum_o.get(mut, 0) + wt

    profiles["other"] = (mut_sum_o, total_wt_o)

    print("Zapisujem profily do súboru...")
    write_combined_profiles(args.output_file, profiles, GENOME_LENGTH)
    print(f"Profily zapísané do {args.output_file}")
