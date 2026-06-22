# Drosophila synonymous-codon SFS and selection analyses

This repository archives the analysis scripts used for the manuscript on
polymorphism-based estimates of selection on synonymous codons in
*Drosophila melanogaster*.  The scripts are organized by analysis stage.  They
are provided as a reproducibility archive: several scripts retain absolute paths
from the original working environment and may need path edits before rerunning
elsewhere.

## License

Code in this archive is released under the MIT License; see `LICENSE`.

## Repository layout

### `SFRatios_pipeline/`

Scripts for building site-frequency spectra, running SFRatios repeatedly, and
summarizing codon-pair selection estimates.

- `get_short_intron_paired_SNP_allele_counts_with_ids.py`
- `make_codon_pair_SFS_from_SNP_paired_allele_counts.py`
- `run_multiple_SFRatios_jobs.py`
- `summarize_multiple_SFRatios_runs.py`
- `setup_bootstrap.py`
- `run_SFRatios_and_LeastSquares_on_bootstrap_samples.py`
- `Leastsquares_2Ns_estimates_with_masking_v2.py`
- `SFRatios.py`
- `SFRatios_functions.py`

### `rooting/`

Scripts used to estimate ancestral states, polarize VCF records, and check
rooting quality.  These scripts expect ancestral-base probability tables,
MAF/FASTA/GTF resources, and annotated VCF files prepared in the original
pipeline.

- `estimate_mel_ancestral_bases.py`
- `make_rooted_EFF_vcfs_manuscript_archive.py`
- `make_rooted_EFF_vcfs_DmelDsim.py`
- `polarize_vcf_with_ancestors_Dsim_and_Dmel.py`
- `polarize_vcf_with_ancestors_selective_EFF_or_ANN.py`
- `compare_rooting_quality.py`

### `gene_expression_work/`

Scripts for the gene-expression analyses added during revision.

- `analyze_gene_expression_codon_selection.py` fits gene-level multinomial
  models of synonymous codon usage as a function of standardized
  `log(FPKM + 1)`.
- `fit_expression_scaled_mutation_selection_model.py` fits the model in which
  gene expression changes the codon-selection scale.

### `figure_generation/`

Scripts for manuscript figures.  Filenames begin with the figure panel they
generate or support.

- `Figure_1A_phenylalanine_example_SFSs.py`
- `Figure_1B_forward_reverse_2Ns_dotplot.py`
- `Figure_1C_g_2Ns_histogram.py`
- `Figure_1D_compare_g_values_across_methods_dotplot.py`
- `Figure_2A_predicted_codon_frequencies_RSCU_dotplot.py`
- `Figure_2B_observed_RSCU_g_dotplot.py`
- `Figure_3A_codon_fitness_vs_expression_slope.py`
- `Figure_3B_helper_observed_frequency_change.py`
- `Figure_3B_helper_predicted_frequency_change.py`
- `Figure_3B_observed_predicted_frequency_change_by_ghat.py`
- `Figure_3B_diagnostic_lambda_vs_expression.py`
- `Figure_3_model_implied_expression_slope_diagnostic.py`
- `Figure_4A_factor1_loadings_by_g_dotplot.py`
- `Figure_4B_codon_stability_by_g_dotplot.py`
- `Figure_4C_RNA_stem_fold_change_by_2Ns.py`
- `Figure_4_helper_stem_loop_SFSs.py`

## Input data

The scripts ingest processed and external data files that are not all included
in this code-only archive.  The main required inputs are:

- annotated *D. melanogaster* VCF files with SnpEff `EFF` or `ANN` annotations;
- ancestral-base estimates or MAF files used to polarize variants;
- short-intron polymorphism data used as the neutral reference class;
- codon-pair SFS files generated from polarized synonymous SNPs;
- SFRatios output files and least-squares summaries of ordered codon-pair
  estimates;
- *D. melanogaster* coding sequences from Ensembl;
- *D. melanogaster* GTF annotation files for transcript/CDS mapping;
- gene-expression data from DGET, summarized across life stages as FPKM;
- processed manuscript tables such as codon fitness estimates, codon counts,
  gene expression summaries, and RNA secondary-structure summaries.

For the revision analyses, the working manuscript directory contained processed
inputs including:

- `Dmel_gene_codon_freqs_and_expression_rank.tsv`
- `codon_fitnesses.tsv`
- `gene_expression_multinomial_log_expression/`
- `expression_scaled_mutation_selection_model_bootstrap200/`
- `figwork/`

Large source data files should be obtained from their original sources or from
the manuscript data archive.  Before rerunning scripts on another system, update
absolute paths near the top of each script or pass equivalent command-line
arguments where available.

## Suggested archival workflow

For public release, create a tagged GitHub release from this repository and
archive that release with Zenodo to mint a DOI.  The DOI can then be cited in
the manuscript code-availability statement.
