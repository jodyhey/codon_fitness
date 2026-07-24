# Drosophila Synonymous-Codon SFS and Selection Analyses

This repository archives the scripts used for the resubmission analyses of the
manuscript on polymorphism-based estimates of selection on synonymous codons in
*Drosophila melanogaster*.  The archive is organized by analysis stage and is
intended to be copied into the public `codon_fitness` GitHub repository.

The code is provided as a reproducibility archive.  Some scripts preserve
absolute paths from the analysis workstation and may need path edits, or
equivalent command-line arguments, before being rerun elsewhere.  Generated
large data files, manuscript documents, figures, caches, and Python bytecode are
not included.

## License

Code in this archive is released under the MIT License; see `LICENSE`.

## Repository Layout

### `SFRatios_pipeline/`

Core SFS and SFRatios scripts:

- count paired synonymous and short-intron SNPs from annotated VCF files;
- convert paired SNP counts into unfolded site-frequency spectra;
- run SFRatios over codon-pair SFS files;
- summarize multiple SFRatios runs;
- estimate directional codon fitness values by least squares;
- set up and run bootstrap SFRatios analyses;
- apply optional polarization-error deconvolution corrections.

Important entry points include:

- `get_short_intron_paired_SNP_allele_counts_with_ids_fixed.py`
- `make_codon_pair_SFS_from_SNP_paired_allele_counts.py`
- `SFRatios.py`
- `run_multiple_SFRatios_jobs.py`
- `summarize_multiple_SFRatios_runs.py`
- `Leastsquares_2Ns_estimates_with_masking_v2.py`
- `bootstrap/setup_bootstrap.py`
- `bootstrap/run_SFRatios_and_LeastSquares_on_bootstrap_samples.py`

### `rooting/`

Scripts used to polarize variants and prepare rooted VCFs.  This includes both
the SINGER/ARG-rooting workflow used for the primary resubmission analyses and
the est-sfs preparation scripts retained for comparison analyses.

Important subdirectories and scripts:

- `rooting/singer_rooting/`
  - `convert_to_tskit_compatible.py`
  - `fix_tskit_backmutation_parents.py`
  - `convert_and_fix_zi_singer_trees.py`
  - `root_vcf_with_singer_trees.py`
- `rooting/est-sfs_prep/`
  - `prepare_est_sfs_dataset.py`
  - `prepare_est_sfs_from_existing_downsample_mafs.py`
  - `make_est_sfs_rooted_vcf.py`
- additional ancestral-state and VCF-polarization helpers:
  - `estimate_mel_ancestral_bases.py`
  - `make_rooted_EFF_vcfs_manuscript_archive.py`
  - `make_rooted_EFF_vcfs_DmelDsim.py`
  - `polarize_vcf_with_ancestors_Dsim_and_Dmel.py`
  - `polarize_vcf_with_ancestors_selective_EFF_or_ANN.py`
  - `compare_rooting_quality.py`

### `gene_expression_work/`

Scripts for gene-expression and expression-scaled mutation-selection-drift
analyses.

- `analyze_gene_expression_codon_selection.py` fits gene-level multinomial
  models of synonymous codon usage as a function of standardized
  `ln(FPKM + 1)`.
- `fit_expression_scaled_mutation_selection_model.py` fits the model in which
  gene expression changes the codon-selection scale:
  `lambda_g = exp(alpha + beta * standardized_log_expression_g)`.
- `run_revision2_gene_expression_work.py` is the workflow wrapper used for the
  resubmission analysis.

### `figure_generation/`

Scripts used to generate manuscript figure panels and helper tables.  Filenames
begin with the figure panel they generate or support.

Examples:

- `Figure_1A_phenylalanine_example_SFSs.py`
- `Figure_1B_forward_reverse_2Ns_dotplot.py`
- `Figure_1C_g_2Ns_histogram.py`
- `Figure_1D_compare_g_values_across_methods_dotplot.py`
- `Figure_2A_observed_RSCU_g_dotplot.py`
- `Figure_2B_predicted_codon_frequencies_RSCU_dotplot.py`
- `Figure_3A_codon_fitness_vs_expression_slope.py`
- `Figure_3B_observed_predicted_frequency_change_by_ghat.py`
- `Figure_4A_factor1_loadings_by_g_dotplot.py`
- `Figure_4B_codon_stability_by_g_dotplot.py`
- `Figure_4C_RNA_stem_fold_change_by_2Ns.py`

The `figure_generation/revision2/` directory contains the most recent wrapper
scripts used to assemble the final resubmission figures.

### `supplementary_information/`

Scripts for extracting codon-fitness tables and rebuilding the Supplementary
Information workbook from processed analysis outputs.

### `workflow/`

High-level wrappers used during the resubmission analysis.  These scripts
coordinate bootstrap sampling, deconvolution analyses, figure generation,
supplement rebuilding, and expression-model refitting.

## Main Analysis Inputs

The scripts assume processed inputs from the manuscript analysis directory,
including:

- ARG-rooted, SnpEff-annotated *D. melanogaster* VCF files;
- SINGER tree files converted to tskit-compatible tree sequences;
- est-sfs rooted VCFs used for comparison analyses;
- short-intron BED files and reference genome FASTA files;
- codon-pair count files;
- SFRatios output summaries and least-squares codon-fitness estimates;
- *D. melanogaster* coding-sequence and gene annotation files;
- gene-expression data summarized across life stages as FPKM;
- RNA stability and RNA secondary-structure summary tables used in Figure 4.

Large input and output datasets are not stored in this code archive.  They
should be obtained from the manuscript data archive or regenerated with the
workflow scripts after updating paths for the local environment.

## Typical Workflow

1. Root/polarize the annotated VCF using SINGER tree files or est-sfs outputs.
2. Count paired synonymous and short-intron SNPs from the rooted VCF.
3. Convert paired SNP counts into codon-pair SFS files.
4. Run SFRatios for all ordered codon pairs.
5. Estimate directional codon fitness values by least squares.
6. Bootstrap paired SNP counts to obtain confidence intervals.
7. Generate figure panels and rebuild the supplementary workbook.
8. Fit gene-expression and expression-scaled mutation-selection-drift models.

## Notes For GitHub Release

This archive is designed to replace the contents of the public
`jodyhey/codon_fitness` code repository for the resubmission.  Before release,
check `ARCHIVE_MANIFEST.txt` and confirm that no large generated data files have
been added accidentally.
