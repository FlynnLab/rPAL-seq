#!/usr/bin/env Rscript
# FINAL replacement panels for Fig 2e and Fig 4b: LOOCV ROC with the dimension count set
# by a declared variance rule -- retain every MDS axis carrying >= 5 % of the variance of
# its own embedding. One rule, stated in Methods, applied identically to all four matrices.
#
# Selected k: cell-line glyco 5, cell-line expr 4, EV glyco 3, EV expr 4.
#
# Note the knife edge: axis 5 of the cell-line EXPRESSION embedding carries 4.8 %, so it
# falls just outside the cut. Had it been retained, expression would read 0.69 instead of
# 0.82. This does not change any conclusion -- DeLong is non-significant either way
# (P = 0.855 as computed, P = 0.353 at the alternative pairing) -- but it should be stated
# rather than discovered. See lineage_stats_reanalysis.md 3e.
#
# Outputs -> revision_figs/roc_var5_{cellline,ev}.{pdf,png} + roc_var5_summary.csv
#         -> EV panel also copied into ../counts_analysis_DiFi/ where its data lives
suppressMessages({library(ggplot2); library(e1071); library(pROC); library(paletteer)
                  library(ragg); library(vegan)})

tp <- path.expand(c("../../theme_nature.R"))
source(tp[file.exists(tp)][1])
luc <- as.character(paletteer_d("PrettyCols::Lucent"))
READOUTS <- c("bold(Enrichment~intensity)~bold((italic(z)))", "bold(Expression)~bold((VST))")
COL <- setNames(c(luc[1], luc[4]), READOUTS)
dir.create("revision_figs", showWarnings = FALSE)

PANELS <- list(
  cellline = list(key="cellline", title="Cell lineage",
    sub="14 human cell lines (7 hematopoietic, 7 non-hematopoietic)",
    lab="Lineage", lev=c("Non-hematopoietic","Hematopoietic"),
    z="ml_input_mds_z_delta.csv", b="ml_input_mds_baseexpr_vst.csv",
    dz="origdist_z.rds", db="origdist_b.rds"),
  ev = list(key="ev", title="Cellular vs secreted",
    sub="8 groups (5 epithelial, 3 EV/NVEP)",
    lab="Category", lev=c("Epithelial","EV/NVEP"),
    z="../counts_analysis_DiFi/ml_input_mds_ev_z_delta.csv",
    b="../counts_analysis_DiFi/ml_input_mds_ev_baseexpr_vst.csv",
    dz="../counts_analysis_DiFi/origdist_ev_z.rds",
    db="../counts_analysis_DiFi/origdist_ev_b.rds")
)

load1 <- function(p, w) {
  d <- read.csv(p[[w]], stringsAsFactors = FALSE)
  y <- factor(d[[p$lab]], levels = p$lev)   # never leave the order to factor() defaults
  stopifnot(!any(is.na(y)))
  list(y = y, X = as.matrix(d[, 3:ncol(d)]))
}
cosine_pairwise <- function(A) {
  n <- nrow(A); D <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in i:n) {
    ok <- is.finite(A[i,]) & is.finite(A[j,])
    d <- if (!any(ok)) 1 else { nu <- sum(A[i,ok]*A[j,ok])
      de <- sqrt(sum(A[i,ok]^2))*sqrt(sum(A[j,ok]^2)); if (de==0) 1 else 1-nu/de }
    D[i,j] <- D[j,i] <- d }
  as.dist(D)
}

VAR_CUT <- 5   # per cent; the rule that goes in Methods

# k = number of MDS axes carrying >= VAR_CUT % of that embedding's variance
select_k <- function(X) { v <- apply(X, 2, var); sum(100 * v / sum(v) >= VAR_CUT) }

fixed_roc <- function(y, X) {
  k <- select_k(X); Xk <- X[, seq_len(k), drop = FALSE]
  s <- numeric(length(y)); set.seed(42)
  for (i in seq_along(y)) {
    m <- svm(x = Xk[-i, , drop = FALSE], y = y[-i],
             type = "C-classification", kernel = "linear")
    # decision values are oriented level1-vs-level2; negate so higher = level 2 (case)
    s[i] <- -as.numeric(attr(predict(m, Xk[i, , drop = FALSE],
                                     decision.values = TRUE), "decision.values")[1,1])
  }
  list(roc = roc(y, s, levels = levels(y), direction = "<", quiet = TRUE), k = k)
}
roc_pts <- function(r, nm) {
  p <- as.data.frame(coords(r, x="all", transpose=FALSE))
  p$fpr <- 1 - p$specificity; p$metric <- nm
  p[rev(seq_len(nrow(p))), ]
}

rows <- list()
for (p in PANELS) {
  fits <- lapply(c("z","b"), function(w) { d <- load1(p, w)
    r <- fixed_roc(d$y, d$X); r$n <- length(d$y); r$dims <- ncol(d$X); r$d <- d; r })
  aucs <- sapply(fits, function(f) as.numeric(auc(f$roc)))
  nms  <- sprintf('%s~~~AUC == "%.2f"', READOUTS, aucs)

  rdf <- rbind(roc_pts(fits[[1]]$roc, nms[1]), roc_pts(fits[[2]]$roc, nms[2]))
  rdf$metric <- factor(rdf$metric, levels = nms)

  ks <- sapply(fits, function(f) f$k)
  g <- ggplot(rdf, aes(fpr, sensitivity, colour = metric)) +
    geom_abline(slope=1, intercept=0, colour="grey55", linewidth=0.35) +
    geom_step(data = subset(rdf, metric==nms[1]), linewidth=1.3, direction="vh") +
    geom_step(data = subset(rdf, metric==nms[2]), linewidth=0.7, direction="vh") +
    scale_colour_manual(values = setNames(unname(COL), nms), name = NULL,
                        labels = function(l) parse(text = l)) +
    scale_x_continuous(limits=c(0,1), breaks=seq(0,1,.25), expand=c(.01,.01)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.25), expand=c(.01,.01)) +
    coord_equal() +
    labs(title = p$title,
         subtitle = sprintf("%s\nLOOCV linear SVM; MDS axes above 5%% variance (k = %d, %d)",
                            p$sub, ks[1], ks[2]),
         x = "1 - Specificity", y = "Sensitivity") +
    theme_nature +
    theme(legend.position = c(0.53, 0.15), legend.background = element_blank(),
          legend.key.height = unit(9, "pt"), legend.text = element_text(size=6.5),
          plot.subtitle = element_text(size=6.5, colour="grey35"))

  ggsave(sprintf("revision_figs/roc_var5_%s.pdf", p$key), g, width=3.2, height=3.3)
  ggsave(sprintf("revision_figs/roc_var5_%s.png", p$key), g, width=3.2, height=3.3,
         dpi=200, device=ragg::agg_png)

  for (i in 1:2) {
    d <- fits[[i]]$d; set.seed(42)
    a <- adonis2(readRDS(p[[c("dz","db")[i]]])$D ~ d$y, permutations = 4999)
    rows[[paste(p$key,i)]] <- data.frame(
      panel = p$title, readout = c("Glycosylation (delta-z)","Expression (VST)")[i],
      n = fits[[i]]$n, mds_dims_available = fits[[i]]$dims,
      AUC = aucs[i],
      k_selected = fits[[i]]$k,
      permanova_F = a$F[1], permanova_R2 = a$R2[1], permanova_P = a$`Pr(>F)`[1])
  }
}
for (f in list.files("revision_figs", "roc_var5_.*\\.pdf$", full.names = TRUE))
  try(grDevices::embedFonts(f), silent = TRUE)

res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, "revision_figs/roc_var5_summary.csv", row.names = FALSE)
print(res, digits = 3)

# the EV panel belongs with its own data as well
file.copy("revision_figs/roc_var5_ev.pdf", "../counts_analysis_DiFi/roc_var5_ev.pdf", overwrite = TRUE)
file.copy("revision_figs/roc_var5_ev.png", "../counts_analysis_DiFi/roc_var5_ev.png", overwrite = TRUE)
cat("copied EV panel into ../counts_analysis_DiFi/\n")
