#!/usr/bin/env Rscript
# Writes webApp/cnlib.R: the verbatim source of every cancerCellNet function
# the PACNet web app calls, directly or indirectly. Web-specific patches are
# applied by hand afterwards and marked with "# PACNet-web:" comments.
#
# Usage (from repo root): Rscript tools/vendor_cnlib.R <cancerCellNet checkout> [out]

args <- commandArgs(trailingOnly = TRUE)
ccn_dir <- args[1]
out <- if (length(args) > 1) args[2] else "webApp/cnlib.R"

fns <- sort(c(
  "addRandToSampTab", "broadClass_predict", "broadClass_train", "ccn_classAssess",
  "cn_clPerf", "cn_eval", "downSampleW", "findBestPairs", "findClassyGenes",
  "getClassGenes", "gnrAll", "logRank", "makeClassifier", "makePairTab",
  "plot_class_PRs", "ptGetTop", "ptSmall", "query_transform", "randomize",
  "rf_classPredict", "sc_compAlpha", "sc_compMu", "sc_cov", "sc_fano",
  "sc_filterGenes", "sc_sampR_to_pattern", "sc_statTab", "sc_testPattern",
  "splitCommon_proportion", "trans_prop", "weighted_down"))

sha <- system2("git", c("-C", shQuote(ccn_dir), "rev-parse", "HEAD"), stdout = TRUE)

# Collate order matches R package default (alphabetical); later definitions win.
chunks <- list()
for (f in sort(list.files(file.path(ccn_dir, "R"), pattern = "\\.R$", full.names = TRUE))) {
  exprs <- parse(f, keep.source = TRUE)
  refs <- attr(exprs, "srcref")
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]])) {
      name <- as.character(e[[2]])
      if (name %in% fns) chunks[[name]] <- paste(as.character(refs[[i]]), collapse = "\n")
    }
  }
}
missing <- setdiff(fns, names(chunks))
if (length(missing)) stop("Not found in ", ccn_dir, ": ", paste(missing, collapse = ", "))

header <- c(
  "# Vendored from https://github.com/pcahan1/cancerCellNet (MIT License),",
  paste0("# commit ", sha, ", by tools/vendor_cnlib.R."),
  "# Only the functions the PACNet web app calls are included. Local changes",
  "# are marked with \"# PACNet-web:\" comments.",
  "")
writeLines(c(header, unlist(lapply(fns, function(n) c(chunks[[n]], "")))), out)
cat("wrote", out, "with", length(fns), "functions\n")
