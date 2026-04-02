# Consolidated reproducibility script for the GPNMB CNS myeloid analyses.
#
# Purpose:
# 1. Document the provenance of the myeloid meta-atlas construction workflow.
# 2. Reproduce the atlas-level summaries from the released integrated Seurat object.
# 3. Re-derive the cross-sample robust NMF programs from the archived sample-wise
#    NMF result object.
# 4. Optionally rebuild the integrated atlas and rerun sample-wise NMF from the
#    cohort-level component objects when those inputs are available.
#
# Source analyses consolidated here:
# - Analyses/ZMyeloid_01_myeloid_meta_analysis_A.Rmd
# - Analyses/ZMyeloid_02_myeloid_integration_NMF.Rmd
# - Analyses/GPNMB_submission/02_export_myeloid_meta_atlas_publishable.Rmd
# - Analyses/GPNMB_submission/03_gpnmb_scRNA_methods_reproducibility.Rmd
#
# Default behavior uses the archived release objects already present in Analyses/.
# That mode is the one intended for the GPNMB submission bundle. The optional
# rebuild/rerun modes are off by default because the full cohort-level inputs are
# not all shipped in this repository snapshot. The robust-program rederivation
# step is also off by default because it is computationally heavy; enable it with
# GPNMB_REDERIVE_MYELOID_PROGRAMS=TRUE when a full rerun is desired.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
})

timestamp_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

as_flag <- function(value, default = FALSE) {
  if (is.null(value) || !nzchar(value)) {
    return(default)
  }
  toupper(value) %in% c("1", "TRUE", "T", "YES", "Y")
}

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
  }

  frame_file <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(frame_file)) {
    return(normalizePath(frame_file, winslash = "/", mustWork = TRUE))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
submission_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
analysis_dir <- normalizePath(file.path(submission_dir, ".."), winslash = "/", mustWork = TRUE)
output_dir <- file.path(submission_dir, "outputs")
stats_dir <- file.path(output_dir, "stats")
fig_dir <- file.path(output_dir, "figures")
repro_object_dir <- file.path(output_dir, "reproducibility_objects")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repro_object_dir, recursive = TRUE, showWarnings = FALSE)

config <- list(
  rebuild_atlas_from_components = as_flag(Sys.getenv("GPNMB_REBUILD_MYELOID_ATLAS"), default = FALSE),
  rerun_sample_nmf = as_flag(Sys.getenv("GPNMB_RERUN_MYELOID_NMF"), default = FALSE),
  rederive_robust_programs = as_flag(Sys.getenv("GPNMB_REDERIVE_MYELOID_PROGRAMS"), default = FALSE),
  top_n_genes = as.integer(Sys.getenv("GPNMB_TOP_N_GENES", "50")),
  intra_threshold = as.numeric(Sys.getenv("GPNMB_INTRA_THRESHOLD", "0.5")),
  inter_threshold = as.numeric(Sys.getenv("GPNMB_INTER_THRESHOLD", "0.3")),
  min_cross_sample_hits = as.integer(Sys.getenv("GPNMB_MIN_CROSS_SAMPLE_HITS", "3")),
  n_robust_clusters = as.integer(Sys.getenv("GPNMB_N_ROBUST_CLUSTERS", "17")),
  prevalence_threshold = as.numeric(Sys.getenv("GPNMB_PREVALENCE_THRESHOLD", "0.3")),
  min_cells_per_sample = as.integer(Sys.getenv("GPNMB_MIN_CELLS_PER_SAMPLE", "50")),
  k_ranks = 2:15
)

config_tbl <- tibble::tibble(
  setting = names(config),
  value = vapply(config, function(x) paste(x, collapse = ","), character(1))
)
readr::write_tsv(config_tbl, file.path(stats_dir, "gpnmb_myeloid_reproducibility_config.tsv"))

pick_first_match <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[[1]]
}

pick_feature_name <- function(object, candidates) {
  feat <- rownames(object)
  hit <- feat[toupper(feat) %in% toupper(candidates)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[[1]]
}

pick_umap_reduction <- function(object, preferred = character()) {
  reduction_names <- names(object@reductions)

  preferred_hit <- preferred[preferred %in% reduction_names]
  if (length(preferred_hit) > 0) {
    return(preferred_hit[[1]])
  }

  umap_hit <- reduction_names[grepl("umap", reduction_names, ignore.case = TRUE)]
  if (length(umap_hit) > 0) {
    return(umap_hit[[1]])
  }

  two_dim_hit <- reduction_names[vapply(
    reduction_names,
    function(reduction_name) {
      embedding <- tryCatch(SeuratObject::Embeddings(object[[reduction_name]]), error = function(e) NULL)
      !is.null(embedding) && ncol(embedding) >= 2
    },
    logical(1)
  )]
  if (length(two_dim_hit) > 0) {
    return(two_dim_hit[[1]])
  }

  stop("No UMAP-like or other 2D reduction found in the Seurat object.")
}

save_plot_pair <- function(plot_handle, stem, width = 8, height = 6) {
  ggsave(filename = paste0(stem, ".pdf"), plot = plot_handle, width = width, height = height, units = "in")
  ggsave(filename = paste0(stem, ".png"), plot = plot_handle, width = width, height = height, units = "in", dpi = 300)
}

named_list_to_padded_df <- function(x) {
  max_len <- max(lengths(x))
  out <- lapply(x, function(values) {
    values <- as.character(values)
    length(values) <- max_len
    values
  })
  tibble::as_tibble(out)
}

flatten_programs <- function(programs_by_rank) {
  flat <- list()
  for (program_set in programs_by_rank) {
    flat <- c(flat, program_set)
  }
  flat
}

unique_length <- function(x) {
  length(unique(x))
}

jaccard_similarity_matrix <- function(gene_sets) {
  gene_sets <- lapply(gene_sets, function(x) unique(as.character(stats::na.omit(x))))
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  stopifnot(length(gene_sets) > 0)

  if (length(gene_sets) > 1 && requireNamespace("scMiko", quietly = TRUE)) {
    return(scMiko::jaccardSimilarityMatrix(gene_sets))
  }

  n_sets <- length(gene_sets)
  set_names <- names(gene_sets)
  jmat <- matrix(0, nrow = n_sets, ncol = n_sets, dimnames = list(set_names, set_names))

  for (i in seq_len(n_sets)) {
    jmat[i, i] <- 1
    if (i == 1) {
      next
    }
    for (j in seq_len(i - 1L)) {
      inter <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
      uni <- length(union(gene_sets[[i]], gene_sets[[j]]))
      score <- if (uni == 0) 0 else inter / uni
      jmat[i, j] <- score
      jmat[j, i] <- score
    }
  }

  jmat
}

pearson_profile_distance <- function(mat) {
  corr <- suppressWarnings(stats::cor(t(mat), method = "pearson", use = "pairwise.complete.obs"))
  corr[is.na(corr)] <- 0
  stats::as.dist(1 - corr)
}

resolve_component_file <- function(file_name, search_roots) {
  search_roots <- unique(search_roots[!is.na(search_roots) & nzchar(search_roots)])
  for (root in search_roots) {
    if (!dir.exists(root)) {
      next
    }

    direct_path <- file.path(root, file_name)
    if (file.exists(direct_path)) {
      return(normalizePath(direct_path, winslash = "/", mustWork = TRUE))
    }

    hits <- list.files(
      path = root,
      pattern = paste0("^", gsub("\\.", "\\\\.", file_name), "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(hits) > 0) {
      return(normalizePath(hits[[1]], winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}

ensure_percent_mt <- function(object) {
  if ("percent.mt" %in% colnames(object@meta.data)) {
    return(object)
  }

  assay_name <- if ("RNA" %in% names(object@assays)) "RNA" else DefaultAssay(object)
  gene_names <- rownames(object[[assay_name]])
  mito_genes <- gene_names[grepl("^MT-|^mt-", gene_names)]

  if (length(mito_genes) == 0) {
    object$percent.mt <- 0
  } else {
    object$percent.mt <- Seurat::PercentageFeatureSet(object, features = mito_genes, assay = assay_name)
  }

  object
}

filter_seurat_object <- function(object) {
  object <- ensure_percent_mt(object)
  subset(
    x = object,
    subset = percent.mt < 10 & nFeature_RNA < 9000 & nFeature_RNA > 200
  )
}

run_sctransform <- function(object, variable_features_n = 2000, conserve_memory = FALSE) {
  object <- filter_seurat_object(object)

  if (requireNamespace("glmGamPoi", quietly = TRUE)) {
    try_result <- tryCatch(
      Seurat::SCTransform(
        object = object,
        method = "glmGamPoi",
        verbose = FALSE,
        vst.flavor = "v2",
        vars.to.regress = "percent.mt",
        variable.features.n = variable_features_n,
        conserve.memory = conserve_memory
      ),
      error = function(e) NULL
    )
    if (!is.null(try_result)) {
      return(try_result)
    }
  }

  Seurat::SCTransform(
    object = object,
    verbose = FALSE,
    vars.to.regress = "percent.mt",
    variable.features.n = variable_features_n,
    conserve.memory = conserve_memory
  )
}

normalize_sample_list <- function(object_list) {
  object_list <- lapply(object_list, filter_seurat_object)
  object_list <- Filter(function(x) ncol(x) > config$min_cells_per_sample, object_list)
  lapply(object_list, run_sctransform)
}

rename_sample_list <- function(object_list) {
  stopifnot(length(object_list) > 0)

  new_names <- vapply(
    seq_along(object_list),
    function(i) {
      object <- object_list[[i]]
      type_value <- unique(as.character(object@meta.data$type))
      study_value <- unique(as.character(object@meta.data$study))
      type_value <- if (length(type_value) == 0) "UNSPECIFIED" else type_value[[1]]
      study_value <- if (length(study_value) == 0) "UNSPECIFIED" else gsub("_", "", study_value[[1]])
      paste0(type_value, "-", i, "_", study_value)
    },
    character(1)
  )

  names(object_list) <- new_names
  object_list
}

build_cell_to_sample_mapping <- function(object_list) {
  mapping <- unlist(lapply(names(object_list), function(sample_name) {
    cells <- colnames(object_list[[sample_name]])
    out <- rep(sample_name, length(cells))
    names(out) <- cells
    out
  }))
  mapping
}

extract_counts_matrix <- function(object) {
  if ("RNA" %in% names(object@assays)) {
    return(GetAssayData(object, assay = "RNA", slot = "counts"))
  }
  if ("SCT" %in% names(object@assays)) {
    return(GetAssayData(object, assay = "SCT", slot = "counts"))
  }
  GetAssayData(object, assay = DefaultAssay(object), slot = "counts")
}

prepare_current_cohort <- function(path) {
  object <- readRDS(path)
  object@meta.data$sample <- object@meta.data[["Barcode"]]
  object@meta.data$type <- paste0(tolower(object@meta.data[["PR"]]), "GBM")
  so_list <- Seurat::SplitObject(object = object, split.by = "sample")
  so_list <- normalize_sample_list(so_list)
  lapply(so_list, function(x) {
    x@meta.data$study <- "Mikolajewicz_2023"
    x
  })
}

prepare_abdel_cohort <- function(path) {
  object <- readRDS(path)
  counts_mat <- extract_counts_matrix(object)
  object <- CreateSeuratObject(counts = counts_mat, meta.data = object@meta.data)
  object@meta.data$type <- paste0(tolower(object@meta.data[["PR"]]), "GBM")
  so_list <- Seurat::SplitObject(object = object, split.by = "sample")
  so_list <- normalize_sample_list(so_list)
  lapply(so_list, function(x) {
    x@meta.data$study <- "Abdelfattah_2022"
    x
  })
}

prepare_wang_cohort <- function(path) {
  object <- readRDS(path)
  object@meta.data$PR2 <- object@meta.data[["PR"]]
  object@meta.data$PR2[grepl("Primary", object@meta.data[["PR"]])] <- "P"
  object@meta.data$PR2[grepl("Recurrent", object@meta.data[["PR"]])] <- "R"
  object@meta.data$type <- paste0(tolower(object@meta.data[["PR"]]), "GBM")
  object@meta.data$sample <- object@meta.data[["Barcode"]]
  so_list <- Seurat::SplitObject(object = object, split.by = "sample")
  so_list <- normalize_sample_list(so_list)
  lapply(so_list, function(x) {
    x@meta.data$study <- "Wang_2022"
    x@meta.data$type <- paste0(tolower(x@meta.data[["PR2"]]), "GBM")
    x
  })
}

prepare_legacy_cohort <- function(path) {
  object_list <- readRDS(path)
  if (inherits(object_list, "Seurat")) {
    stop("Legacy cohort file is a Seurat object, but a sample-split list was expected: ", path)
  }
  object_list
}

build_atlas_from_components <- function(component_lists) {
  if (!requireNamespace("scMiko", quietly = TRUE)) {
    stop("Package 'scMiko' is required to rebuild the BBKNN myeloid atlas.")
  }

  run_bbknn <- get("runBBKNN", envir = asNamespace("scMiko"))

  all_samples <- unlist(component_lists, recursive = FALSE)
  if (length(all_samples) == 0) {
    stop("No component sample lists were provided for atlas rebuilding.")
  }

  timestamp_message("Merging cohort-level sample objects for atlas rebuilding.")
  so_merge <- merge(all_samples[[1]], y = all_samples[-1])
  so_merge <- UpdateSeuratObject(so_merge)
  DefaultAssay(so_merge) <- "RNA"

  timestamp_message("Running SCTransform / PCA / BBKNN integration.")
  so_merge <- run_sctransform(so_merge, variable_features_n = 3000, conserve_memory = TRUE)
  DefaultAssay(so_merge) <- "SCT"
  so_merge <- RunPCA(so_merge, features = VariableFeatures(so_merge), verbose = FALSE)
  so_merge <- run_bbknn(so_merge, batch = "sample")
  so_merge <- FindClusters(so_merge, graph.name = "bbknn", resolution = 0.5, verbose = FALSE)
  so_merge <- FindClusters(so_merge, graph.name = "bbknn", resolution = 1.0, verbose = FALSE)

  saveRDS(so_merge, file.path(repro_object_dir, "seurat_all_myeloid_160523_rebuilt_for_submission.rds"))
  so_merge
}

get_expressed_genes <- function(object, min_pct = 0.005) {
  counts_mat <- GetAssayData(object, assay = "SCT", slot = "counts")
  detected_fraction <- Matrix::rowSums(counts_mat > 0) / ncol(counts_mat)
  names(detected_fraction)[detected_fraction >= min_pct]
}

get_nmf_genes <- function(feature_loading, norm_cutoff = 0.5, n_cutoff = NA_integer_) {
  stopifnot(is.matrix(feature_loading))

  nmf_kme <- t(feature_loading)
  nmf_kme <- apply(nmf_kme, 2, function(x) (x^2) / sum(x^2))

  if (!is.na(n_cutoff)) {
    module_genes <- apply(nmf_kme, 1, function(x) colnames(nmf_kme)[order(-x)][seq_len(min(n_cutoff, length(x)))])
    module_genes <- as.list(as.data.frame(module_genes, stringsAsFactors = FALSE))
    module_genes <- lapply(module_genes, function(x) unique(as.character(stats::na.omit(x))))
  } else {
    module_genes <- apply(nmf_kme, 1, function(x) colnames(nmf_kme)[x > norm_cutoff], simplify = FALSE)
  }

  module_genes[lengths(module_genes) > 0]
}

run_samplewise_nmf <- function(sample_list, k_ranks) {
  if (!requireNamespace("NNLM", quietly = TRUE) || !requireNamespace("NMF", quietly = TRUE)) {
    stop("Packages 'NNLM' and 'NMF' are required to rerun the sample-wise NMF models.")
  }

  results <- list()

  for (sample_name in names(sample_list)) {
    timestamp_message("Running sample-wise NMF for ", sample_name)
    object <- sample_list[[sample_name]]

    expr_genes <- get_expressed_genes(object, min_pct = 0.005)
    expr_genes <- expr_genes[!grepl("^mt-|^rps|^rpl", tolower(expr_genes))]
    object <- Seurat::GetResidual(object, features = expr_genes, assay = "SCT")

    expr_mat <- GetAssayData(object, assay = "SCT", slot = "scale.data")
    expr_mat <- expr_mat[rownames(expr_mat) %in% expr_genes, , drop = FALSE]
    expr_mat[expr_mat < 0] <- 0

    nmf_model_list <- lapply(k_ranks, function(k_rank) {
      fit <- NNLM::nnmf(expr_mat, k = k_rank, verbose = FALSE, check.k = FALSE)
      NMF::nmfModel(H = fit$H, W = fit$W)
    })
    names(nmf_model_list) <- paste0(sample_name, "_k", k_ranks)

    gene_programs <- lapply(names(nmf_model_list), function(model_name) {
      current_model <- nmf_model_list[[model_name]]
      current_genes <- get_nmf_genes(feature_loading = current_model@W, norm_cutoff = NA, n_cutoff = config$top_n_genes)
      model_suffix <- sub(paste0("^", sample_name, "_"), "", model_name)
      names(current_genes) <- paste0(model_suffix, "_", seq_along(current_genes))
      current_genes
    })
    names(gene_programs) <- names(nmf_model_list)

    gene_programs_flat <- flatten_programs(gene_programs)

    results[[sample_name]] <- list(
      sample = sample_name,
      nmf = nmf_model_list,
      programs = gene_programs,
      programs.flat = gene_programs_flat,
      jmat = jaccard_similarity_matrix(gene_programs_flat)
    )
  }

  saveRDS(results, file.path(repro_object_dir, "NMF_myeloid_v1_160523_rebuilt_for_submission.rds"))
  results
}

derive_robust_programs <- function(nmf_results) {
  nmf_results_refreshed <- nmf_results

  for (sample_name in names(nmf_results_refreshed)) {
    nmf_model_list <- nmf_results_refreshed[[sample_name]]$nmf
    gene_programs <- lapply(names(nmf_model_list), function(model_name) {
      current_model <- nmf_model_list[[model_name]]
      current_genes <- get_nmf_genes(
        feature_loading = current_model@W,
        norm_cutoff = NA,
        n_cutoff = config$top_n_genes
      )
      model_suffix <- sub(paste0("^", sample_name, "_"), "", model_name)
      names(current_genes) <- paste0(model_suffix, "_", seq_along(current_genes))
      current_genes
    })
    names(gene_programs) <- names(nmf_model_list)
    nmf_results_refreshed[[sample_name]]$programs <- gene_programs
    nmf_results_refreshed[[sample_name]]$programs.flat <- flatten_programs(gene_programs)
  }

  nmf_gene_all <- list()

  for (sample_name in names(nmf_results_refreshed)) {
    sample_programs <- nmf_results_refreshed[[sample_name]]$programs.flat
    sample_jmat <- jaccard_similarity_matrix(sample_programs)
    robust_rows <- rownames(sample_jmat)[rowSums(sample_jmat > config$intra_threshold) > 1]

    if (length(robust_rows) == 0) {
      next
    }

    sample_programs <- sample_programs[robust_rows]
    names(sample_programs) <- paste0(sample_name, "_", names(sample_programs))
    nmf_gene_all <- c(nmf_gene_all, sample_programs)
  }

  strip_component_suffix <- function(x) {
    stringr::str_remove(x, "_k[0-9]+_[0-9]+$")
  }

  jmat_all <- jaccard_similarity_matrix(nmf_gene_all)
  robust_inter <- apply(jmat_all, 1, function(x) {
    partners <- names(x)[x > config$inter_threshold]
    unique_length(strip_component_suffix(partners))
  })
  robust_component_names <- names(robust_inter)[robust_inter >= config$min_cross_sample_hits]

  jmat_robust <- jmat_all[robust_component_names, robust_component_names, drop = FALSE]
  robust_tree <- hclust(pearson_profile_distance(jmat_robust))
  robust_clusters <- cutree(robust_tree, k = config$n_robust_clusters)

  cluster_annotation <- tibble::tibble(
    component = names(robust_clusters),
    cluster = paste0("G", robust_clusters),
    cross_sample_support = robust_inter[names(robust_clusters)]
  )

  consensus_programs <- list()
  for (cluster_id in sort(unique(robust_clusters))) {
    cluster_name <- paste0("G", cluster_id)
    cluster_members <- names(robust_clusters)[robust_clusters == cluster_id]
    cluster_gene_sets <- nmf_gene_all[cluster_members]
    n_members <- length(cluster_gene_sets)
    tally_df <- as.data.frame(table(unlist(cluster_gene_sets)), stringsAsFactors = FALSE)
    tally_df$prop <- tally_df$Freq / n_members
    consensus_programs[[cluster_name]] <- as.character(tally_df$Var1[tally_df$prop > config$prevalence_threshold])
  }

  list(
    nmf_results_refreshed = nmf_results_refreshed,
    nmf_gene_all = nmf_gene_all,
    jmat_all = jmat_all,
    jmat_robust = jmat_robust,
    cluster_annotation = cluster_annotation,
    consensus_programs = consensus_programs
  )
}

match_programs_by_jaccard <- function(derived_programs, archived_programs) {
  overlap_rows <- lapply(names(derived_programs), function(derived_name) {
    best_name <- NA_character_
    best_score <- NA_real_

    for (archived_name in names(archived_programs)) {
      derived_genes <- unique(as.character(derived_programs[[derived_name]]))
      archived_genes <- unique(as.character(archived_programs[[archived_name]]))
      inter <- length(intersect(derived_genes, archived_genes))
      uni <- length(union(derived_genes, archived_genes))
      score <- if (uni == 0) 0 else inter / uni

      if (is.na(best_score) || score > best_score) {
        best_name <- archived_name
        best_score <- score
      }
    }

    tibble::tibble(
      derived_program = derived_name,
      archived_program = best_name,
      best_jaccard = best_score
    )
  })

  dplyr::bind_rows(overlap_rows)
}

component_root <- Sys.getenv("GPNMB_MYELOID_COMPONENT_DIR", unset = "")
component_root <- if (nzchar(component_root)) normalizePath(component_root, winslash = "/", mustWork = TRUE) else NA_character_

source_manifest <- tibble::tribble(
  ~category, ~label, ~path, ~role,
  "legacy_script", "Per-cohort myeloid extraction", file.path(analysis_dir, "ZMyeloid_01_myeloid_meta_analysis_A.Rmd"), "Raw-data to cohort-level myeloid sample lists",
  "legacy_script", "Atlas integration and myeloid NMF", file.path(analysis_dir, "ZMyeloid_02_myeloid_integration_NMF.Rmd"), "Cohort harmonization, integration, and sample-wise NMF",
  "release_script", "Publishable atlas export", file.path(submission_dir, "02_export_myeloid_meta_atlas_publishable.Rmd"), "Submission-ready atlas export",
  "release_script", "Methods reproducibility notebook", file.path(submission_dir, "03_gpnmb_scRNA_methods_reproducibility.Rmd"), "Submission-ready figure/stat reproduction"
)
readr::write_tsv(source_manifest, file.path(stats_dir, "gpnmb_myeloid_source_manifest.tsv"))

object_manifest <- tibble::tribble(
  ~label, ~path, ~role, ~status,
  "Integrated myeloid object", file.path(analysis_dir, "seurat_all_myeloid_160523.rds"), "Primary release object", "required",
  "Myeloid NMF object", file.path(analysis_dir, "NMF_myeloid_v1_160523.rds"), "Primary archived sample-wise NMF result", "required",
  "Archived myeloid program gene sets", file.path(analysis_dir, "ZMyeloid_MyCat_NMF_genesets_050524.rds"), "Archived release gene-set object", "optional",
  "Archived inferred states", file.path(analysis_dir, "ZM_myeloid_inferred_states.rds"), "Archived state calls", "optional",
  "Current cohort object", file.path(analysis_dir, "seurat_object_myeloid_current_150523.rds"), "Local cohort object used upstream", "optional",
  "Abdelfattah cohort object", file.path(analysis_dir, "seurat_object_myeloid_abdel_150523.rds"), "Local cohort object used upstream", "optional",
  "Wang cohort object", file.path(analysis_dir, "seurat_object_myeloid_wang_150523.rds"), "Local cohort object used upstream", "optional"
)
object_manifest$exists <- file.exists(object_manifest$path)
readr::write_tsv(object_manifest, file.path(stats_dir, "gpnmb_myeloid_object_manifest.tsv"))

cohort_catalog <- tibble::tribble(
  ~cohort_key, ~input_kind, ~expected_file,
  "current", "local_seurat_object", "seurat_object_myeloid_current_150523.rds",
  "abdel", "local_seurat_object", "seurat_object_myeloid_abdel_150523.rds",
  "wang", "local_seurat_object", "seurat_object_myeloid_wang_150523.rds",
  "franjic", "legacy_sample_list", "seurat_object_myeloid_franjic.rds",
  "kanton", "legacy_sample_list", "seurat_object_myeloid_kanton.rds",
  "khrameeva", "legacy_sample_list", "seurat_object_myeloid_khrameeva.rds",
  "jakel", "legacy_sample_list", "seurat_object_myeloid_jakel.rds",
  "schirmer", "legacy_sample_list", "seurat_object_myeloid_schirmer.rds",
  "sun", "legacy_sample_list", "seurat_object_myeloid_sun.rds",
  "kim", "legacy_sample_list", "seurat_object_myeloid_kim.rds",
  "biermann", "legacy_sample_list", "seurat_object_myeloid_biermann.rds",
  "heming", "legacy_sample_list", "seurat_object_myeloid_heming.rds",
  "wei", "legacy_sample_list", "seurat_object_myeloid_wei.rds",
  "riemond", "legacy_sample_list", "seurat_object_myeloid_riemond.rds",
  "lau", "legacy_sample_list", "seurat_object_myeloid_lau.rds",
  "smajic", "legacy_sample_list", "seurat_object_myeloid_smajic.rds",
  "yu", "legacy_sample_list", "seurat_object_myeloid_yu.rds",
  "fan", "legacy_sample_list", "seurat_object_myeloid_fan.rds",
  "cao", "legacy_sample_list", "seurat_object_myeloid_cao.rds",
  "bhaduri", "legacy_sample_list", "seurat_object_myeloid_bhaduri.rds",
  "aldinger", "legacy_sample_list", "seurat_object_myeloid_aldinger.rds"
)

search_roots <- unique(c(analysis_dir, component_root))
cohort_catalog$resolved_path <- vapply(cohort_catalog$expected_file, resolve_component_file, character(1), search_roots = search_roots)
cohort_catalog$available <- !is.na(cohort_catalog$resolved_path)
readr::write_tsv(cohort_catalog, file.path(stats_dir, "gpnmb_myeloid_component_availability.tsv"))

component_lists <- list()

if (config$rebuild_atlas_from_components || config$rerun_sample_nmf) {
  missing_components <- cohort_catalog$expected_file[!cohort_catalog$available]
  if (length(missing_components) > 0) {
    stop(
      "Component-level rebuild requested but the following cohort objects were not found: ",
      paste(missing_components, collapse = ", "),
      ". Set GPNMB_MYELOID_COMPONENT_DIR to the directory containing the legacy cohort files."
    )
  }

  timestamp_message("Preparing cohort-level sample lists for optional rebuild/rerun mode.")
  component_lists$current <- rename_sample_list(prepare_current_cohort(cohort_catalog$resolved_path[cohort_catalog$cohort_key == "current"]))
  component_lists$abdel <- rename_sample_list(prepare_abdel_cohort(cohort_catalog$resolved_path[cohort_catalog$cohort_key == "abdel"]))
  component_lists$wang <- rename_sample_list(prepare_wang_cohort(cohort_catalog$resolved_path[cohort_catalog$cohort_key == "wang"]))

  legacy_keys <- cohort_catalog$cohort_key[cohort_catalog$input_kind == "legacy_sample_list"]
  for (legacy_key in legacy_keys) {
    component_lists[[legacy_key]] <- rename_sample_list(
      prepare_legacy_cohort(cohort_catalog$resolved_path[cohort_catalog$cohort_key == legacy_key])
    )
  }

  c2n_mapping <- build_cell_to_sample_mapping(unlist(component_lists, recursive = FALSE))
  saveRDS(c2n_mapping, file.path(repro_object_dir, "PR11_sample2cell_mapping_submission.rds"))
}

if (config$rebuild_atlas_from_components) {
  so_myeloid <- build_atlas_from_components(component_lists)
  atlas_object_path <- file.path(repro_object_dir, "seurat_all_myeloid_160523_rebuilt_for_submission.rds")
} else {
  atlas_object_path <- file.path(analysis_dir, "seurat_all_myeloid_160523.rds")
  stopifnot(file.exists(atlas_object_path))
  so_myeloid <- readRDS(atlas_object_path)
}

timestamp_message("Loaded myeloid atlas object: ", basename(atlas_object_path))

myeloid_umap <- pick_umap_reduction(so_myeloid, preferred = "b")
myeloid_meta <- so_myeloid@meta.data %>%
  tibble::rownames_to_column("cell_id")

study_col <- pick_first_match(colnames(myeloid_meta), c("study", "Study"))
etiology_col <- pick_first_match(colnames(myeloid_meta), c("etiology", "type", "diagnosis", "PR", "PR2"))
class_col <- pick_first_match(colnames(myeloid_meta), c("class_consensus", "class", "Class", "subcluster", "cell_type3"))
sample_col <- pick_first_match(colnames(myeloid_meta), c("sample", "clean.id2", "clean.id", "Barcode", "orig.ident"))

myeloid_overview <- tibble::tribble(
  ~field, ~value,
  "Atlas object", basename(atlas_object_path),
  "UMAP reduction", myeloid_umap,
  "Number of cells", as.character(ncol(so_myeloid)),
  "Number of genes", as.character(nrow(so_myeloid)),
  "Study column", ifelse(is.na(study_col), "not present", study_col),
  "Etiology column", ifelse(is.na(etiology_col), "not present", etiology_col),
  "Class column", ifelse(is.na(class_col), "not present", class_col),
  "Sample column", ifelse(is.na(sample_col), "not present", sample_col)
)
readr::write_tsv(myeloid_overview, file.path(stats_dir, "gpnmb_myeloid_atlas_overview.tsv"))

counts_by_group <- list()
if (!is.na(study_col)) {
  counts_by_group$study <- myeloid_meta %>% count(.data[[study_col]], name = "n_cells")
}
if (!is.na(etiology_col)) {
  counts_by_group$etiology <- myeloid_meta %>% count(.data[[etiology_col]], name = "n_cells")
}
if (!is.na(class_col)) {
  counts_by_group$class <- myeloid_meta %>% count(.data[[class_col]], name = "n_cells")
}
if (!is.na(sample_col)) {
  counts_by_group$sample <- myeloid_meta %>% count(.data[[sample_col]], name = "n_cells")
}
for (name in names(counts_by_group)) {
  readr::write_tsv(counts_by_group[[name]], file.path(stats_dir, paste0("gpnmb_myeloid_counts_by_", name, ".tsv")))
}

if (!is.na(study_col)) {
  p_study <- DimPlot(so_myeloid, reduction = myeloid_umap, group.by = study_col, raster = TRUE) +
    ggtitle("GPNMB myeloid meta-atlas", subtitle = "Study / cohort")
  save_plot_pair(p_study, file.path(fig_dir, "gpnmb_myeloid_umap_by_study"), width = 10, height = 8)
}

if (!is.na(etiology_col)) {
  p_etiology <- DimPlot(so_myeloid, reduction = myeloid_umap, group.by = etiology_col, raster = TRUE) +
    ggtitle("GPNMB myeloid meta-atlas", subtitle = "Etiology / diagnosis")
  save_plot_pair(p_etiology, file.path(fig_dir, "gpnmb_myeloid_umap_by_etiology"), width = 10, height = 8)
}

if (!is.na(class_col)) {
  p_class <- DimPlot(so_myeloid, reduction = myeloid_umap, group.by = class_col, raster = TRUE) +
    ggtitle("GPNMB myeloid meta-atlas", subtitle = "Class / subtype")
  save_plot_pair(p_class, file.path(fig_dir, "gpnmb_myeloid_umap_by_class"), width = 10, height = 8)
}

gpnmb_feature <- pick_feature_name(so_myeloid, c("GPNMB", "Gpnmb"))
if (!is.na(gpnmb_feature)) {
  p_gpnmb <- FeaturePlot(so_myeloid, reduction = myeloid_umap, features = gpnmb_feature, raster = TRUE) +
    ggtitle("GPNMB myeloid meta-atlas", subtitle = paste0(gpnmb_feature, " expression"))
  save_plot_pair(p_gpnmb, file.path(fig_dir, "gpnmb_myeloid_featureplot_gpnmb"), width = 8, height = 6)
}

archived_gene_set_path <- file.path(analysis_dir, "ZMyeloid_MyCat_NMF_genesets_050524.rds")
archived_gene_sets_available <- file.exists(archived_gene_set_path)

robust_nmf <- NULL
nmf_object_path <- NA_character_

if (config$rederive_robust_programs || !archived_gene_sets_available) {
  if (config$rerun_sample_nmf) {
    nmf_results <- run_samplewise_nmf(unlist(component_lists, recursive = FALSE), k_ranks = config$k_ranks)
    nmf_object_path <- file.path(repro_object_dir, "NMF_myeloid_v1_160523_rebuilt_for_submission.rds")
  } else {
    nmf_object_path <- file.path(analysis_dir, "NMF_myeloid_v1_160523.rds")
    stopifnot(file.exists(nmf_object_path))
    nmf_results <- readRDS(nmf_object_path)
  }

  timestamp_message("Loaded myeloid NMF object: ", basename(nmf_object_path))
  robust_nmf <- derive_robust_programs(nmf_results)
  consensus_programs <- robust_nmf$consensus_programs
  program_source_tbl <- tibble::tibble(
    source_mode = "rederived_from_nmf",
    atlas_object = basename(atlas_object_path),
    nmf_object = basename(nmf_object_path)
  )
} else {
  archived_gene_sets <- readRDS(archived_gene_set_path)
  consensus_programs <- lapply(archived_gene_sets, as.character)
  names(consensus_programs) <- sub("^MyCat_", "", names(consensus_programs))
  program_source_tbl <- tibble::tibble(
    source_mode = "archived_release_gene_sets",
    atlas_object = basename(atlas_object_path),
    nmf_object = NA_character_
  )
}

readr::write_tsv(program_source_tbl, file.path(stats_dir, "gpnmb_myeloid_program_source.tsv"))

programs_long <- bind_rows(lapply(names(consensus_programs), function(program_name) {
  tibble::tibble(
    program = program_name,
    rank = seq_along(consensus_programs[[program_name]]),
    gene = consensus_programs[[program_name]]
  )
}))
programs_wide <- named_list_to_padded_df(consensus_programs)
readr::write_tsv(programs_long, file.path(stats_dir, "gpnmb_myeloid_consensus_nmf_programs_long.tsv"))
readr::write_tsv(programs_wide, file.path(stats_dir, "gpnmb_myeloid_consensus_nmf_programs_wide.tsv"))

if (!is.null(robust_nmf)) {
  cluster_annotation <- robust_nmf$cluster_annotation %>%
    dplyr::arrange(cluster, dplyr::desc(cross_sample_support), component)
  readr::write_tsv(cluster_annotation, file.path(stats_dir, "gpnmb_myeloid_nmf_component_cluster_annotation.tsv"))
}

program_summary <- programs_long %>%
  count(program, name = "n_genes") %>%
  arrange(program)
readr::write_tsv(program_summary, file.path(stats_dir, "gpnmb_myeloid_nmf_program_summary.tsv"))

if (!is.null(robust_nmf) && archived_gene_sets_available) {
  archived_gene_sets <- readRDS(archived_gene_set_path)
  archived_gene_sets <- lapply(archived_gene_sets, as.character)
  overlap_tbl <- match_programs_by_jaccard(consensus_programs, archived_gene_sets)
  overlap_tbl <- overlap_tbl %>%
    left_join(program_summary, by = c("derived_program" = "program"))
  readr::write_tsv(overlap_tbl, file.path(stats_dir, "gpnmb_myeloid_nmf_program_overlap_vs_archived_release.tsv"))

  overlap_heatmap <- expand.grid(
    derived_program = names(consensus_programs),
    archived_program = names(archived_gene_sets),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(
      jaccard = {
        derived_genes <- consensus_programs[[derived_program]]
        archived_genes <- archived_gene_sets[[archived_program]]
        inter <- length(intersect(derived_genes, archived_genes))
        uni <- length(union(derived_genes, archived_genes))
        if (uni == 0) 0 else inter / uni
      }
    ) %>%
    ungroup()

  p_overlap <- overlap_heatmap %>%
    ggplot(aes(x = archived_program, y = derived_program, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "grey95", high = "firebrick", limits = c(0, 1)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Derived versus archived myeloid NMF programs",
      x = "Archived release program",
      y = "Derived program",
      fill = "Jaccard"
    )
  save_plot_pair(p_overlap, file.path(fig_dir, "gpnmb_myeloid_nmf_overlap_heatmap"), width = 9, height = 7)
}

session_info_lines <- capture.output(sessionInfo())
writeLines(session_info_lines, con = file.path(stats_dir, "gpnmb_myeloid_reproducibility_session_info.txt"))

timestamp_message("Finished GPNMB myeloid meta-atlas and NMF reproducibility workflow.")
