# ISBT V15 → V16 migration delta

Release: **V16 (June 2026)**, applied to `data/config/` on 2026-06-10.
Engine: unchanged (data-only). Method: full per-allele re-fetch + `variants[]` diff.

## Variant table (`variation_annotation.dat`)

4 new variants; 5 old single-SNV / placeholder rows removed as they were
consolidated. Net row count unchanged (1767 data rows).

| Action | System | HGVS (transcript) | hg19 | hg38 | Notes |
|---|---|---|---|---|---|
| + | AUG  | `c.589+1G>C`        | chr6:44198219  | 44230482 | replaces placeholder `c.589+1=` |
| + | CTL2 | `c.455A>G`          | chr19:10742170 | 10631494 | new SLC44A2 reference SNP |
| + | DO   | `c.431_432delinsAA` | chr12:14993800 | 14840866 | consolidates `c.431C>A`+`c.432C>A` |
| + | GYPB | `c.71_72delinsGT`   | chr4:144922402 | 144001249 | consolidates `c.71A>G`+`c.72G>T` |
| + | KN   | `c.4828A>T`         | chr1:207782916 | —        | added across the KN system |
| − | AUG  | `c.589+1=`          | | | placeholder removed |
| − | DO   | `c.431C>A`, `c.432C>A` | | | superseded by delins |
| − | GYPB | `c.71A>G`, `c.72G>T`   | | | superseded by delins |

All new variants carry both hg19 and hg38 coordinates → **no liftOver required**.

## Alleles with changed variant composition (24)

Authoritative diff of per-allele `variants[]` sets, V15 vs V16:

| Allele | System | change |
|---|---|---|
| `KN*01`, `KN*01.-05`, `KN*01.-13`, `KN*01.06`, `KN*01.07`, `KN*01.10`, `KN*01.12`, `KN*01W`, `KN*02`, `KN*02.10` | KN | **+ `c.4828A>T`** (systematic; 10 alleles) |
| `GYPB*03N.02`, `GYPB*03N.04`, `GYPB*03N.05`, `GYPB*06.01`, `GYPB*06.02` | MNS | `71A>G`+`72G>T` → `71_72delinsGT` |
| `CTL2*01`, `CTL2*01.-02`, `CTL2*01.-05`, `CTL2*01N.01` | SLC44A2 | **+ `c.455A>G`** |
| `DO*02.-07` | ART4 | `431C>A`+`432C>A` → `431_432delinsAA` |
| `AUG*01N` | SLC29A1 | `589+1=` → `589+1G>C` |
| `ABCB6*01N.24`, `ABCB6*01N.27` | LAN | dropped spurious `c.459del` |

## Alleles retired in V16 (−3)

`DI*02.04`, `JK*01W.11`, `ABCB6*01N.27` — removed from the in-release allele set.

## Schema note

V16 adds a `show_comment` boolean to the variant API schema. It is irrelevant to
bloodAGENT's `.dat` extraction (the build scripts read named fields only) and
requires no handling.

## Validation results

- `validate_dat.py`: **0 issues** on both `variation_annotation.v16.dat` and
  `genotype_to_phenotype_annotation.v16.dat` (header order OK, 28/10 cols, all
  `stoi()` integer columns parse, strand ∈ {+,−}, unique index keys).
- `diff_against_current.py` vs live HGDP config: 0 added, 3 retired, 0 phenotype
  regressions.
- Cross-file join integrity (every gt2pt `base_change` token resolvable in the
  variation table): **0 new** unresolved tokens vs V15 (the 15 pre-existing ones
  are coverage-detected/legacy alleles, unchanged).
- exonic_annotation files: untouched (no new systems).

## Applied files

`data/config/variation_annotation.dat` + 7 per-pipeline
`genotype_to_phenotype_annotation_*.dat` (CMR, Dragen, HGDP, Microarray, ONT,
PacBio×2). Each original backed up as `<file>.pre-v16.bak`.
