#!/usr/bin/env bash
# run_ldsc.sh -- official S-LDSC (+ MAGMA) production recipe (issue #1).
#
# These runs need the 1000 Genomes EUR reference panel and the official
# tools, which are too large/heavy for the analysis sandbox, so they were NOT
# run here. The committed results in reports/ come from the in-repo
# MAGMA-style pipeline (scripts/run_magma_enrichment.py) run on the same
# programs. This script is the exact, runnable recipe for the production run.
#
# Prereqs:
#   conda create -n ldsc -c bioconda ldsc magma
#   (or the community python3 port: pip install git+https://github.com/bulik/ldsc.git)
set -euo pipefail

DATA=${DATA:-data}          # staged inputs (see docs/DATA_SOURCES.md)
OUT=${OUT:-results_ldsc}    # outputs stay out of git
mkdir -p "$OUT"

# ------------------------------------------------- 0. 1000G EUR reference
# Download once (~2 GB total):
#   wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz
#   wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
#   wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
#   wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
REF=${REF:-ref}             # directory where the tarballs were unpacked

# -------------------------------------- 1. harmonize + munge summary stats
for TRAIT in AD_Bellenguez2022 PD_Nalls2019; do
  # produces SNP/A1/A2/N/Z columns from the staged harmonised files
  python -m celltype_heritability.sumstats harmonize \
      --src "$DATA/gwas/${TRAIT}.raw.tsv.gz" \
      --dest "$OUT/$TRAIT.harmonised.tsv.gz"
  ldsc munge_sumstats \
      --sumstats "$OUT/$TRAIT.harmonised.tsv.gz" \
      --snp rsid --N-col n --a1 effect_allele --a2 other_allele \
      --out "$OUT/$TRAIT" \
      --merge-alleles "$REF/w_hm3.snplist"
done

# ------------------ 2. marker programs -> per-SNP annotations -> LD scores
# scripts/programs_to_annot.py converts reports/gene_programs/*.json into
# per-chromosome .annot.gz files against the 1000G EUR BIM positions with a
# 10 kb flank (same window as the MAGMA-style run):
#   python scripts/programs_to_annot.py \
#       --programs reports/gene_programs/programs_cluster_top100.json \
#       --bim-prefix "$REF/1000G_EUR_Phase3_plink/1000G.EUR.QC." \
#       --window 10000 --out "$OUT/annots/cluster"
for CHR in {1..22}; do
  for ANNOT in "$OUT"/annots/cluster/*.${CHR}.annot.gz; do
    ldsc --l2 --bfile "$REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR}" \
        --ld-wind-cm 1 --annot "$ANNOT" --thin-annot \
        --out "${ANNOT%.annot.gz}" \
        --print-snps "$REF/hm3_no_hla.snps.${CHR}"
  done
done
# one .ldcts file per hierarchy level, then:
#   printf '%s\t%s\n' "$PROGRAM" "$OUT/annots/cluster/$PROGRAM." >> "$OUT/whb_cluster.ldcts"

# --------------------------------------------- 3. cell-type partitioned h2
for TRAIT in AD_Bellenguez2022 PD_Nalls2019; do
  ldsc --h2-cts "$OUT/$TRAIT.sumstats.gz" \
      --ref-ld-chr "$REF/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD." \
      --ref-ld-chr-cts "$OUT/whb_cluster.ldcts" \
      --w-ld-chr "$REF/weights_hm3_no_hla/weights." \
      --out "$OUT/${TRAIT}_whb_cluster"
done

# ----------------------------------------- 4. official MAGMA cross-check
# magma --annotate window=10,10 --snp-loc "$REF/g1000_eur.bim" \
#       --gene-loc "$DATA/genes/ensembl_grch38_protein_coding.bed" \
#       --out "$OUT/magma_annot"
# for TRAIT in AD_Bellenguez2022 PD_Nalls2019; do
#   magma --bfile "$REF/g1000_eur" --gene-annot "$OUT/magma_annot.genes.annot" \
#         --pval "$OUT/$TRAIT.magma_pval.txt" ncol=n --out "$OUT/magma_$TRAIT"
#   magma --gene-results "$OUT/magma_$TRAIT.genes.raw" \
#         --set-annot reports/gene_programs/programs_supercluster_top100.genesets.txt \
#         --out "$OUT/magma_${TRAIT}_supercluster"
# done

echo "Inputs and checksums: docs/DATA_SOURCES.md"
