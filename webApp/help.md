## How to use PACNet

PACNet scores how closely your engineered cells resemble their target cell or tissue type, compares them with published engineered-cell protocols, and ranks transcription factors whose perturbation is predicted to improve the result.

Everything runs in your browser: your data never leaves your computer.

### Quick start

1. Prepare two CSV files, a **sample metadata table** and an **expression matrix** (formats below). To try the app first, download the example files: [example_sample_metadata_table.csv](https://github.com/CahanLab/PACNet/raw/main/example_data/example_sample_metadata_table.csv) and [example_counts_matrix.csv](https://github.com/CahanLab/PACNet/raw/main/example_data/example_counts_matrix.csv).
2. On the **Analyze** tab, keep *CellNet analysis* selected and upload both files.
3. Choose the **target cell/tissue type**, the type your protocol aims to make.
4. Click **Submit**. Results appear on the right; a progress bar shows each step.

### Input files

**Sample metadata table** (CSV, one row per sample)

| Column | Required | Meaning |
|---|---|---|
| `sample_name` | yes | Sample identifier; must match a column name of the expression matrix |
| `description1` | yes | Group label (e.g. protocol or time point) shown in the plots |
| anything else | no | Ignored |

**Expression matrix** (CSV)

- First column: human **gene symbols** (e.g. `ALB`, `NKX2-5`). Other identifiers such as Ensembl IDs will not match.
- One column per sample, named exactly as in `sample_name`; column order does not matter.
- Counts, TPM or FPKM are all accepted.
- Genes listed more than once: the row with the highest median is kept. Genes with missing values are dropped.

### Target cell/tissue types

| Choice | Classification | Engineered reference comparison | Network influence scores |
|---|---|---|---|
| heart, hspc, intestine_colon, liver, lung, neuron, skeletal_muscle | yes | yes | yes |
| esc, endothelial_cell, fibroblast, kidney, b_cell, t_cell, monocyte_macrophage | yes | no | yes |
| none_of_the_above | yes | no | no |

### Results

**Classification.** A heatmap with one column per sample and one row per cell/tissue type in the training data. Each score runs from 0 (distinct from that type) to 1 (indistinguishable from it). Three randomized profiles are included as a baseline. For the seven types with reference panels, a second heatmap shows how published engineered samples score, so you can compare your protocol with others.

**NIS (network influence score).** For each sample, transcription factors in the target type's gene regulatory network. Negative scores suggest the factor should be *up*-regulated to move closer to the target type; positive scores suggest *down*-regulation. The plots show selected factors; *Download Network Influence Scores* gives every factor.

**New Training Validation.** Appears only when retraining was needed (next section). It shows how well the retrained classifier separates held-out training samples.

Both result tables can be downloaded as CSV.

### Retraining

The classifier uses a fixed set of genes. If your matrix is missing some of them, PACNet retrains a classifier on the genes you do have. This takes several minutes in the browser; the progress bar says when it is happening. If fewer than 1,000 genes are shared, you will be warned that the classifier may be unreliable, which usually means the row names are not gene symbols.

### Browsing the reference data

Choose *Browse reference data* in **What do you want to do?**, pick a panel, and click **Submit** to see how published engineered samples classify. No upload is needed.

### Loading times

- **First visit:** the browser downloads the R runtime and analysis packages, which takes about a minute. Later visits are faster.
- **Each target type:** its reference data (a few MB up to about 15 MB) is downloaded the first time you use it in a session.

### Troubleshooting

- **"does not share any genes"** or **"fewer than 1000 genes"**: make sure the first column holds human gene symbols.
- **"Column names of expression matrix DO NOT match"** or **"samples … are not columns"**: every `sample_name` must be a column of the expression matrix and vice versa.
- **The page becomes unresponsive with very large matrices:** the browser has limited memory. [Run PACNet in R on your own computer](https://github.com/CahanLab/PACNet#webapp-local) instead.
