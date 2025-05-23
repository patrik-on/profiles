## Stručná dokumentácia:

- 1 Vstupné súbory a ich zdroje Obsah
   - 1.1 usher_barcodes.csv
   - 1.2 lineages.yml.txt
   - 1.3 clades.tsv
   - 1.4 genome.fa
- 2 Spustenie skriptu
- 3 Nastavenia skriptu
- 4 Zhrnutie postupu


## 1 Vstupné súbory a ich zdroje Obsah

### 1.1 usher_barcodes.csv

- Popis:Definuje mutácie charakteristické pre jednotlivé SARS-CoV-2 línie.
- Formát:CSV (Comma Separated Values).
- Zdroj:Repozitár Freyja.
- Download:https://raw.githubusercontent.com/andersen-lab/Freyja/main/freyja/data/usher_barcodes.csv
- Argument:-balebo–barcode-file

### 1.2 lineages.yml.txt

- Popis:Obsahuje hierarchiu (strom) SARS-CoV-2 línií a ich rodičovské vzťahy.
- Formát:YAML (ak sťahujete priamo, možno uložiť ako .txt).
- Zdroj:Repozitár outbreak.info.
- Download:https://raw.githubusercontent.com/outbreak-info/outbreak.info/master/curated_reports_prep/lineages.yml
- Argument:-lalebo–lineage-file

### 1.3 clades.tsv

- Popis:Obsahuje počty sekvencií priradené ku každému variantu, pričom kľúčový je stĺpec
    exclusive_count.
- Formát:TSV s hlavičkou (stĺpce:clade,exclusive_count,inclusive_count).
- Zdroj:UCSC Genome Browser – UShER SARS-CoV-2 strom.
- Download:https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/
- Argument:-calebo–clades-file

### 1.4 genome.fa

- Popis:Referenčný genóm SARS-CoV-2.
- Formát:FASTA.
- Zdroj:NCBI.
- Download:https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta
- Argument:-galebo–genome-file


## 2 Spustenie skriptu

- Skript spúšťajte z príkazového riadku:  
    ```bash
    python virpool_profile_builder.py
    ```  
  alebo jednorazovo so všetkými argumentmi:  
    ```bash
    python ./virpool_profile_builder.py -b usher_barcodes.csv -l lineages.yml.txt -c clades.tsv -g genome.fa -o calculated_profiles.tsv -v B.1.617.2 B.1.1.7
    ```

- Ak vstupné súbory majú predvolené názvy a nachádzajú sa v rovnakom adresári, postačí spustenie prvého príkazu.
- Pre zmenu predvolených nastavení použite uvedené argumenty (`-b`, `-l`, `-c`, `-g`, `-o`, `-v`).

## Príkazové argumenty skriptu



| Skratka            | Dlhá forma           | Typ               | Predvolené                    | Popis                                                                                  |
|--------------------|----------------------|-------------------|-------------------------------|----------------------------------------------------------------------------------------|
| `-b`               | `--barcode-file`     | `str`             | `usher_barcodes.csv`          | Cesta k CSV súboru s barkódmi (mutáciami).                                             |
| `-l`               | `--lineage-file`     | `str`             | `lineages.yml.txt`            | Cesta k YAML súboru s definíciami línií a rodičovských vzťahov.                         |
| `-c`               | `--clades-file`      | `str`             | `clades.tsv`                  | Cesta k TSV súboru s počtami sekvencií pre jednotlivé počty (`inclusive_count`, `exclusive_count`). |
| `-g`               | `--genome-file`      | `str`             | `genome.fa`                   | Cesta k FASTA súboru s referenčným genómom.                                            |
| `-o`               | `--output-file`      | `str`             | `calculated_profiles.tsv`     | Cesta k výstupnému TSV súboru, kde sa zapíšu vypočítané profily.                        |
| `-v`               | `--variants`         | `list[str]`       | `["B.1.617.2", "B.1.1.7"]`     | Zoznam variantov, pre ktoré sa majú počítať vážené profily. Zadajte ľubovoľný počet.    |


## 3 Nastavenia skriptu

- Vstupné súbory:
    - usher_barcodes.csv
    - lineages.yml.txt
    - clades.tsv
    - genome.fa
- Výstupný súbor:Predvolenecalculated_profiles.tsv.

## 4 Zhrnutie postupu

1. Stiahnite potrebné vstupné súbory z uvedených odkazov.
2. Uistite sa, že máte nainštalovaný Python 3 a knižnicu PyYAML.
3. Skript upravuje cesty k súborom a zoznam variantov pomocou príkazových argumentov.
4. Spustite skript z príkazového riadku a skontrolujte výsledný výstupný súbor.

## Príkazové argumenty skriptu

| Skratka               | Dlhá forma                  | Typ     | Predvolené                | Popis                                                                                  |
|-----------------------|-----------------------------|---------|---------------------------|----------------------------------------------------------------------------------------|
| `-g`                  | `--genome-file`             | `str`   | `genome.fa`               | Cesta k FASTA súboru s referenčným genómom.                                            |
| `-p`                  | `--profiles-file`           | `str`   | `calculated_profiles.tsv` | Cesta k TSV súboru s vypočítanými profilmi variantov (výstup prvého skriptu).          |
| `--primary-mask-file` |                             | `str`   | `ont-short.masking.vcf`   | Cesta k pôvodnému VCF so základnou maskou (primárna maska, vstup do finalizácie).     |
| `--uniform-mask-file` |                             | `str`   | `uniformity_mask.vcf`     | Cesta k dočasnému VCF, kde sú pozície uniformné naprieč všetkými variantmi.            |
| `-o`                  | `--output-file`             | `str`   | `final.vcf`           | Cesta k výslednému VCF, ktorý kombinuje pôvodnú a uniformitu masku.                    |
| `--chrom`             |                             | `str`   | `MN908947.3`              | Názov chromozómu/kontigu použitý v generovaných VCF riadkoch.                          |

---

## Výstupné súbory

- **Uniform mask VCF** (`--uniform-mask-file`):  
  Vygenerovaný pomocou funkcie `generate_uniformity_mask_vcf()`. Obsahuje hlavičky VCF a záznamy so stavom filtra `mask` pre pozície, kde je rozloženie báz uniformné naprieč všetkými variantmi.

- **Finálny VCF** (`-o` / `--output-file`):  
  Kombinuje pôvodný maskovací VCF (`--primary-mask-file`) so uniform maskou. Skript doplňuje chýbajúce riadky (gap‐fill) z uniform masky do sekvencie prvej časti a zachováva zvyšok originálneho súboru.

---

## Spustenie

```bash
python ./mask.py -g genome.fa -p calculated_profiles.tsv --primary-mask-file ont-short.masking.vcf --uniform-mask-file uniformity_mask.vcf -o final.vcf --chrom MN908947.3


