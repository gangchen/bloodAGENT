# bloodAGENT × ISBT V16 — data integration workspace

This directory aligns bloodAGENT's `.dat` knowledge base with the ISBT Blood Group
Database **V16** snapshot (release **June 2026**, `releaseVersion: 16`). It mirrors
the `data/source/v15/` toolchain — the C++ engine is unchanged; this is a
**data-only** refresh on top of V15.

> **TL;DR** V16 is a *small, contained* delta over V15: no new genes, systems, or
> antigens, no liftOver needed. 4 new variants, 24 alleles re-pointed, 3 alleles
> retired. See [`MIGRATION_V16.md`](MIGRATION_V16.md) for the exact change list.

```
data/source/v16/
├── README.md                   # this file
├── SYSTEM_MAPPING.md           # ISBT 001-048 ↔ bloodAGENT system/gene mapping (unchanged vs V15)
├── MIGRATION_V16.md            # exact V15→V16 delta + validation results
├── raw/                        # raw API dumps (gitignored when large)
│   ├── system.json   gene.json   antigen.json   allele.json
│   ├── variant.json  phenotype.json  publication.json  genbank.json
│   ├── release.json            # release index (now includes v16)
│   ├── release_16.json         # V16 release object + updatedAlleleIds
│   └── alleles/<id>.json       # per-allele detail WITH variants[] (2053 files, full V16 snapshot)
├── tools/                      # same scripts as v15, reparametrized for v16
└── derived/                    # build outputs (variation_annotation.v16.dat, gt2pt, pipeline/, reports)
```

## ISBT V16 counts (from `/api/release/16`)

| | V15 | V16 | Δ |
|---|---|---|---|
| Blood group systems | 48 | 48 | 0 |
| Genes | 57 | 57 | 0 |
| Antigens | 397 | 397 | 0 |
| Alleles (in-release) | 2004 | 2001 | −3 |
| Variants | 1828 | 1832 | +4 |

Bulk endpoint superset is 2053 alleles / 1869 variants in both releases; the
release counts above reflect the approved/in-release subset.

## ⚠️ Lesson learned: `updatedAlleleIds` is NOT a complete delta

The V16 release object lists 75 `updatedAlleleIds`. **This is incomplete** — 20 of
the 24 alleles whose *variant composition* actually changed (e.g. the entire KN
system gaining `c.4828A>T`, and the CTL2 alleles gaining `c.455A>G`) are **not**
in that list, and their `updatedAt`/`version` fields in the bulk `allele.json`
were not bumped either. **You must do a full per-allele re-fetch and diff the
`variants[]` sets** to get the true delta. The incremental shortcut is unsafe.

## End-to-end pipeline

```bash
# 1) Fetch full V16 snapshot from the public API (~15-30 min; per-IP throttled)
data/source/v16/tools/fetch_isbt_v16.sh

# 2) Build bloodAGENT-format master files
python3 data/source/v16/tools/build_variation_annotation.py \
  --raw-dir data/source/v16/raw --allele-detail-dir data/source/v16/raw/alleles \
  --out-dir data/source/v16/derived --old-master data/config/variation_annotation.dat
python3 data/source/v16/tools/build_gt2pt.py \
  --raw-dir data/source/v16/raw --allele-detail-dir data/source/v16/raw/alleles \
  --out-dir data/source/v16/derived
python3 data/source/v16/tools/build_pipeline_gt2pt.py \
  --master data/source/v16/derived/genotype_to_phenotype_annotation.v16.dat \
  --out-root data/source/v16/derived/pipeline

# 3) Validate (mimics CISBTAnno::readAnnotation + generateIndex — catches stoi() crashes)
python3 data/source/v16/tools/validate_dat.py \
  data/source/v16/derived/variation_annotation.v16.dat \
  data/source/v16/derived/genotype_to_phenotype_annotation.v16.dat

# 4) Diff vs current config (allele-level)
python3 data/source/v16/tools/diff_against_current.py \
  --current data/config/HGDP/genotype_to_phenotype_annotation_HGDP.dat \
  --v16 data/source/v16/derived/pipeline/HGDP/genotype_to_phenotype_annotation_HGDP.v16.dat \
  --out-dir data/source/v16/derived

# 5) Apply to data/config (creates <file>.pre-v16.bak backups; idempotent)
data/source/v16/tools/apply_to_config.sh --apply
#   undo: for f in $(find data/config -name '*.pre-v16.bak'); do mv "$f" "${f%.pre-v16.bak}"; done
```

## What stayed manual (same as V15)

- **No exonic_annotation changes.** V16 adds no new systems/genes, so
  `exonic_annotation.{hg19,hg38}.BGStarget.txt` is untouched. All 4 new V16
  variants fall in genes already in the target panel (SLC29A1/AUG, GYPB/MNS,
  ART4/DO, SLC44A2/CTL2) and ship with both hg19+hg38 coords (no liftOver).
- **Pipeline-specific overlay files** (RHCE 109bp insertion per VCF caller) remain
  hand-curated and are unaffected by V16.
- **Container regression** is the authoritative final check. It needs a built
  `bloodagent:v16` image (`docker build -f Dockerfile.cn -t bloodagent:v16 .`)
  and then `data/source/v16/tools/run_regression.sh bloodagent:v16`. Because the
  change is data-only with no new join breakage (see MIGRATION_V16.md), this was
  deferred during data prep, exactly as in V15.
