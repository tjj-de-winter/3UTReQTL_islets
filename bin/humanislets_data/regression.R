#!/usr/bin/env Rscript

# Description: Limma regression of variant allele dosage vs GSIS metrics and protein abundance. 

library(limma)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
csv_file <- args[1]

outfile <- gsub(".csv", ".regression.csv", csv_file)

df <- read.csv(csv_file, check.names = FALSE)
names(df) <- trimws(names(df))

GSIS_cols    <- names(df)[startsWith(names(df), "GSIS-")]
variant_cols <- names(df)[startsWith(names(df), "eQTL-")]
protein_cols <- names(df)[startsWith(names(df), "Proteinexp-")]
meta_cols    <- names(df)[startsWith(names(df), "meta-")]

# Covariates
covars <- c("meta-donorage", "meta-donorsex", "meta-puritypercentage",
            "meta-beta_end", "meta-alpha_end", "meta-bodymassindex", "meta-hba1c")

make_pheno <- function(var_col) {
  ph <- df[, c(var_col, covars), drop=FALSE]
  if ("meta-donorsex" %in% names(ph) && !is.factor(ph[["meta-donorsex"]])) {
    ph[["meta-donorsex"]] <- factor(ph[["meta-donorsex"]])
  }
  ph[[var_col]] <- as.integer(ph[[var_col]])
  ph
}

si_gsis <- t(as.matrix(df[, GSIS_cols, drop=FALSE]))
rownames(si_gsis) <- GSIS_cols
expr_prot <- t(as.matrix(df[, protein_cols, drop=FALSE]))
rownames(expr_prot) <- protein_cols

.drop_invariant_covars <- function(ph, var_col) {
  keep <- rep(TRUE, ncol(ph))
  for (j in seq_len(ncol(ph))) {
    nm <- colnames(ph)[j]; cj <- ph[[j]]
    if (nm == var_col) next
    if (is.factor(cj) && nlevels(droplevels(cj)) < 2) keep[j] <- FALSE
    if ((is.numeric(cj) || is.integer(cj)) && stats::var(cj, na.rm=TRUE) == 0) keep[j] <- FALSE
  }
  ph[, keep, drop=FALSE]
}

make_design <- function(ph, predictor_orig) {
  safe <- make.names(colnames(ph), unique=TRUE)
  idx <- match(predictor_orig, colnames(ph))
  if (is.na(idx)) return(list(X=NULL, coef_name=NA_character_))
  coef_name <- safe[idx]
  ph2 <- ph; colnames(ph2) <- safe
  X <- model.matrix(reformulate(safe), data=ph2)
  if ((nrow(X) - qr(X)$rank) <= 0) return(list(X=NULL, coef_name=NA_character_))
  list(X=X, coef_name=coef_name)
}

safe_limma_tab <- function(Y, X, coef_name) {
  if (is.null(X) || is.na(coef_name)) return(NULL)
  if (!(coef_name %in% colnames(X))) return(NULL)
  if (is.null(Y) || ncol(Y) != nrow(X) || nrow(Y) < 1) return(NULL)
  if ((nrow(X) - qr(X)$rank) <= 0) return(NULL)
  fit <- tryCatch(lmFit(Y, X), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  fit2 <- tryCatch(eBayes(fit), error=function(e) NULL)
  if (is.null(fit2)) return(NULL)
  coef_idx <- match(coef_name, colnames(X))
  if (is.na(coef_idx)) return(NULL)
  tab <- tryCatch(topTable(fit2, coef=coef_idx, number=Inf, sort.by="none"), error=function(e) NULL)
  tab
}

parse_variant <- function(v) {
  no_pref <- sub("^eQTL-", "", v)
  parts <- strsplit(no_pref, "_")[[1]]
  list(chrom_pos=parts[1], ref=parts[2], alt=parts[3], gene=parts[4])
}

valid_snp <- function(x) {
  x <- x[is.finite(x)]
  length(x) >= 5 && length(unique(x)) >= 2 && stats::var(x) > 0
}

gsis_out <- list(); prot_out <- list()

for (var_col in variant_cols) {
  ph <- make_pheno(var_col)
  ok <- complete.cases(ph)
  ph <- ph[ok,,drop=FALSE]
  if (!valid_snp(ph[[var_col]])) next
  ph <- .drop_invariant_covars(ph, var_col)
  des <- make_design(ph, var_col)
  Y <- si_gsis[, ok, drop=FALSE]

  tab <- safe_limma_tab(Y, des$X, des$coef_name)
  if (is.null(tab)) next

  parts <- parse_variant(var_col)
  res <- data.frame(analysis="GSIS", predictor=var_col, response=rownames(tab), n=ncol(Y),
                    estimate=tab$logFC, t_value=tab$t, p_value=tab$P.Value,
                    adj_P_within=tab$adj.P.Val, AveExpr=tab$AveExpr, B=tab$B,
                    chrom_pos=parts$chrom_pos, ref=parts$ref, alt=parts$alt, gene=parts$gene,
                    stringsAsFactors=FALSE)
  gsis_out[[length(gsis_out)+1]] <- res
}

for (var_col in variant_cols) {
  parts <- parse_variant(var_col)
  target <- paste0("Proteinexp-", parts$gene)
  if (!(target %in% rownames(expr_prot))) next

  ph <- make_pheno(var_col)
  ok <- complete.cases(ph)
  ph <- ph[ok,,drop=FALSE]
  if (!valid_snp(ph[[var_col]])) next
  ph <- .drop_invariant_covars(ph, var_col)
  des <- make_design(ph, var_col)
  Y <- expr_prot[target, ok, drop=FALSE]

  tab <- safe_limma_tab(Y, des$X, des$coef_name)
  if (is.null(tab)) next

  row <- tab[1, , drop=FALSE]
  rownames(row) <- rownames(Y)[1]
  if (nrow(row) != 1) next

  res <- data.frame(analysis="Proteomics", predictor=var_col, response=target, n=ncol(Y),
                    estimate=row$logFC, t_value=row$t, p_value=row$P.Value,
                    adj_P_within=row$adj.P.Val, AveExpr=row$AveExpr, B=row$B,
                    chrom_pos=parts$chrom_pos, ref=parts$ref, alt=parts$alt, gene=parts$gene,
                    stringsAsFactors=FALSE)
  prot_out[[length(prot_out)+1]] <- res
}

all_results <- bind_rows(if (length(gsis_out)) bind_rows(gsis_out) else NULL,
                         if (length(prot_out)) bind_rows(prot_out) else NULL)

write.csv(all_results, outfile, row.names=FALSE)
message(gsub("<OUTFILE>", outfile, "Done with limma. Results written to <OUTFILE>"))