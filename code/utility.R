
normCounts <-function(x, log=FALSE, prior.count=0.5)
  # Function to normalise to median library size instead of counts per million
  # Input is DGEList object
  # Belinda Phipson
  # 30 November 2015
{
  lib.size <- x$samples$lib.size*x$samples$norm.factors
  M <- median(lib.size)
  if(log){
    prior.count.scaled <- lib.size/mean(lib.size)*prior.count
    lib.size <- lib.size + 2*prior.count.scaled
    log2(t((t(x$counts)+prior.count.scaled)/lib.size*M))
  }
  else t(t(x$counts)/lib.size*M)
}

intDat <- function(seu, split = "donor", type = c("RNA","ADT"),
                   reference = NULL, k.weight = 100, adt.norm = NULL,
                   int.assay.name = "integrated"){

  type <- match.arg(type)
  seuLst <- SplitObject(seu, split.by = split)
  if (!is.null(reference)) reference <- which(names(seuLst) %in% reference)

  if(type == "RNA"){

    seuLst <- lapply(X = seuLst, FUN = SCTransform, method = "glmGamPoi")

    features <- SelectIntegrationFeatures(object.list = seuLst,
                                          nfeatures = 3000)
    seuLst <- PrepSCTIntegration(object.list = seuLst, anchor.features = features)
    seuLst <- lapply(X = seuLst, FUN = RunPCA, features = features)
    anchors <- FindIntegrationAnchors(object.list = seuLst,
                                      normalization.method = "SCT",
                                      anchor.features = features,
                                      reference = reference,
                                      dims = 1:30, reduction = "rpca",
                                      k.anchor = 20)
    int <- IntegrateData(anchorset = anchors, k.weight = k.weight,
                         normalization.method = "SCT",
                         dims = 1:30, new.assay.name = int.assay.name)

  } else if (type == "ADT"){

    seuLst <- lapply(X = seuLst, FUN = function(x) {
      VariableFeatures(x) <- rownames(x)
      if(is.null(adt.norm)) {
        x <- NormalizeData(x, verbose = FALSE, normalization.method = "CLR",
                           margin = 2)
      }
      x
    })

    features <- SelectIntegrationFeatures(object.list = seuLst)
    seuLst <- lapply(X = seuLst, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
      x
    })

    anchors <- FindIntegrationAnchors(object.list = seuLst, reduction = "rpca",
                                      dims = 1:30, reference = reference)
    int <- IntegrateData(anchorset = anchors, dims = 1:30,
                         new.assay.name = int.assay.name)
  }

  int
}

test_my_model_fixed <-
  function (model, testing.version = "hsa.latest", custom.dataset = NULL,
            target = NULL, plot = T, maxRank = 1500)
  {
    performance.computation <- ifelse(is.null(target), F, T)
    if (class(custom.dataset) == "Seurat") {
      testing.datasets <- list()
      testing.datasets$user.dataset <- custom.dataset
      custom = TRUE
    } else {
      custom = FALSE
    }
    if (!custom) {
      targets <- c("immune", "Lymphoid", "Myeloid", "Tcell",
                   "Bcell", "CD8T", "CD4T", "NK", "MoMacDC", "Plasma_cell",
                   "PanBcell")
      if (is.null(target)) {
        message("warning: target cell_type not provided. Avoiding performance computation")
        performace.computation = F
      }
      else if (!target %in% targets) {
        stop(sprintf("target must be one of %s; or NULL for avoiding performance computation",
                     paste(targets, collapse = "';'")))
      }
      available.datasets = c("hsa.latest")
      if (!testing.version %in% available.datasets) {
        stop("please, provide a valid testing.version paramter or provide a custom.dataset in seurat format")
      }
      if (testing.version == "hsa.latest") {
        testing.datasets <- get_testing_data(version = testing.version)
      }
    }
    if (custom) {
      if (!"cell_type" %in% colnames(custom.dataset@meta.data)) {
        stop("please, provide a 'cell_type' column to be used as reference cell type")
      }
      if (is.null(target)) {
        message("warning: target cell_type not provided. Avoiding performance computation")
        performace.computation = F
      }
      else if (any(!target %in% custom.dataset$cell_type)) {
        stop("all target celltypes must be included in cell_type metadata field. Otherwise, set target = NULL for avoiding performance computation")
      }
    }
    plt.out <- list()
    perf.out <- list()
    output <- list()
    for (dset in names(testing.datasets)) {
      obj <- testing.datasets[[dset]]
      plt <- list()
      dropcols = obj@meta.data %>% colnames() %>% grep("^is.pure",
                                                       ., value = T) %>% unique()
      if (length(dropcols) > 0) {
        for (col in dropcols) {
          obj[[col]] <- NULL
        }
      }
      obj <- scGate(obj, model = model, assay = DefaultAssay(obj),
                    maxRank = maxRank)
      nname <- sprintf("%s's manual annot", dset)
      plt[[nname]] <- DimPlot(obj, group.by = "cell_type",
                              label = T, repel = T, label.size = 3) + ggtitle(nname) +
        NoLegend() + theme(aspect.ratio = 1)
      level.plots <- plot_levels(obj)
      names(level.plots) <- paste(dset, names(level.plots),
                                  sep = "_")
      plt <- c(plt, level.plots)
      plt.out[[dset]] <- patchwork::wrap_plots(plt, ncol = length(plt))
      if (performance.computation) {
        if (!custom) {
          performance = scGate::performance.metrics(actual = obj@meta.data[,
                                                                           target], pred = obj$is.pure == "Pure")
        }
        else {
          performance = scGate::performance.metrics(actual = obj$cell_type %in%
                                                      target, pred = obj$is.pure == "Pure")
        }
        perf.out[[dset]] <- performance
      }
      output[[dset]] <- obj
    }
    if (performance.computation) {
      perf <- bind_rows(perf.out) %>% data.frame
      rownames(perf) <- names(perf.out)
    }
    if (performance.computation) {
      return(list(performance = perf, plots = plt.out, objects = output))
    }
    else {
      return(list(plots = plt.out, objects = output))
    }
  }
