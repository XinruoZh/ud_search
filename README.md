# ud_search — Downstream Gene Neighborhood Analysis

Analyzes the genomic neighborhood downstream of a target gene across a collection of bacterial genomes. Runs after [allele_search](https://github.com/...) produces per-query BLAST hit FASTAs.

---

## Prerequisites

- `allele_search` must be run first to produce `standard_dna_seq/{query}.fasta` hit files
- Conda environment `ud_search` with: `prokka`, `emapper.py`, `parallel`, `pyyaml`
- Databases: Prokka HMM, EggNOG-mapper (`emapperdb-5.0.2`)

---

## Workflow

### Step 1 — Genome annotation + neighborhood extraction

For each query gene × genome set:

| Sub-step | What it does | Script |
|---|---|---|
| 1A | Run Prokka on all genomes; extract neighborhood FAA around target gene | `run_annotate_prokka.sh` |
| 1B | Pool unique neighborhood proteins; run EggNOG-mapper | `run_annotate_eggnog.sh` |

**Run all queries at once** (rgg144, rgg939, rgg1518, TprA × S_mitis, S_pneumo, S_oralis):

```bash
bash scripts/run_all_sequence_annotation.sh          # full run
bash scripts/run_all_sequence_annotation.sh --test   # test mode: 10 genomes per set
```

Per-query YAML configs are auto-generated to `scripts/configs/{query}_prokka.yaml` and `{query}_eggnog.yaml`.

**Run a single query manually:**

```bash
# edit scripts/annotate_prokka.yaml first
bash scripts/run_annotate_prokka.sh [scripts/annotate_prokka.yaml] [--test]
bash scripts/run_annotate_eggnog.sh [scripts/annotate_eggnog.yaml]
```

### Step 2 — Gene chain analysis and visualization

```bash
# edit scripts/analysis.yaml first
bash scripts/run_analysis.sh [scripts/analysis.yaml]
```

---

## Configuration

| File | Used by | Purpose |
|---|---|---|
| `scripts/annotate_prokka.yaml` | `run_annotate_prokka.sh` | Manual single-query Prokka config |
| `scripts/annotate_eggnog.yaml` | `run_annotate_eggnog.sh` | Manual single-query EggNOG config |
| `scripts/configs/{query}_prokka.yaml` | auto-generated | Per-query Prokka config (created by `run_all_sequence_annotation.sh`) |
| `scripts/configs/{query}_eggnog.yaml` | auto-generated | Per-query EggNOG config (created by `run_all_sequence_annotation.sh`) |
| `scripts/analysis.yaml` | `run_analysis.sh` | Post-annotation analysis config |

---

## Output structure

```
annotation_results/
  {genome_name}_{query}/
    annotation/          # Prokka GFFs and FAAs (one subdir per strain)
    eggnog_subset_faa/   # Per-strain neighborhood FAAs
    eggnog/
      unique_proteins.faa          # Deduplicated neighborhood proteins
      protein_origins.tsv          # Maps unique_id → strain, locus, position
      {genome_name}_{query}.emapper.annotations   # EggNOG results
    log/
```

---

## Key scripts in `code/`

| Script | Purpose |
|---|---|
| `extract_neighborhood_faa.py` | Extracts ≤N CDS on each side of the target gene from a Prokka GFF+FAA pair |
| `pool_unique_proteins.py` | Pools and deduplicates neighborhood FAAs across all strains |
