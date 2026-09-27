# TNBC Response Residuals

A reproducible analysis of longitudinal treatment-response transcriptomics in triple-negative breast cancer using public NeoTRIP data.

## Study concept

Observed early molecular response is decomposed into an expected component predicted from pretreatment state and a patient-specific residual component. Residuals are generated with patient-level cross-fitting without using pathological complete response (pCR).

## Primary dataset

GSE319641 / NeoTRIP comprises 401 RNA-seq samples from 251 patients, including 150 verified baseline-D1C2 pairs: 63 with pCR and 87 with residual disease. All analysis uses public data requiring no controlled access.

## Analysis

- Sample reconstruction and QC
- InstaPrism composition estimation
- Six frozen Hallmark pathways
- singscore pathway scoring
- Cross-fitted ridge expected-response models
- Composition-adjusted pCR association
- Sensitivity and negative-control analyses
- Nested prediction analysis
- Orthogonal single-cell interpretation

## Main findings

No pathway passed the prespecified multiplicity-corrected confirmatory threshold. Immune/interferon and apoptosis showed the most coherent hypothesis-generating residual associations. Composition adjustment materially affected several associations. Residual features did not show consistent incremental pCR prediction.

GSE122630 transportability was not performed because no pathway passed the prespecified eligibility gate. GSE263995 was used only for descriptive orthogonal cell-type interpretation.

## Repository structure

- `R/`: analysis functions
- `scripts/`: analysis stages and audits
- `config/`: analysis configuration
- `metadata/`: dataset provenance, frozen specifications, and checksums
- `results/`: statistical outputs and validation records
- `tables/`: analysis tables
- `figures/`: analysis figures
- `logs/`: execution records, methodological decisions, and software versions

## Reproducibility

The R environment is captured in `renv.lock`. Raw public data are not stored in Git. Dataset provenance and checksums are recorded in `metadata/data_inventory.tsv`. Analysis stages and methodological decisions are recorded in `logs/decision_log.tsv`.

## Data availability

The completed analysis used public accessions GSE319641 and GSE263995. Other datasets described in the research planning documents were not necessarily analysed.

## License

License to be selected before public release.
