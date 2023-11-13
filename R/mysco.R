setClass("mysco",slots= c(meta="data.frame",
                          counts="matrix",
                          normalised="matrix",
                          hvg="data.frame",
                          data.scale="matrix",
                          principal_components="matrix",
                          louvain.clustering = "numeric"
                          ))

setMethod("show", signature = c("mysco"),
          definition = function(object)
          {
            cat(paste("An object of class", class(object)),"\n")
            cat(paste("with", nrow(object@counts), "genes and", ncol(object@counts),"cells"),"\n")
          })

setValidity("mysco", function(object)
{
  if (ncol(object@counts) != nrow(object@meta))
  {
    "@counts and @meta must be the same length"
  }
  else if (any(colnames(object@counts) != rownames(object@meta)))
  {
    "@counts columns and @meta rows must have the same names"
  }
  else
  {
    TRUE
  }
})

CreateMySCO = function(matrix)
{
  n_count_rna = apply(pbmc.data, MARGIN = 2, sum)
  n_feature_rna = apply(pbmc.data > 0, MARGIN = 2, sum)
  meta = data.frame(n_count = n_count_rna, n_feature = n_feature_rna)
  pbmc <- new("mysco", meta=meta, counts=matrix)
  return(pbmc)
}

CalcMitoPct = function(mysco, pattern)
{
  n_mt_count_rna = apply(mysco@counts[grep(pattern, row.names((mysco@counts))),], MARGIN = 2, sum)
  mysco@meta$mt_percent = n_mt_count_rna / mysco@meta$n_count
  return(mysco)
}

MakeQCPlots = function(mysco)
{
  library(ggplot2)
  library(gridExtra)
  # make violin plots
  plots = lapply(names(mysco@meta), function(category)
  {
    ggplot(mysco@meta, aes(x = 1, y = mysco@meta[[category]])) +
      geom_violin() +
      geom_jitter(shape = ".", position = position_jitter(0.2)) +
      labs(x = category, y = "the rest of the plot")
  })

  # make scatter plots of n_features and mt_percent vs n_count
  plots = c(plots, lapply(c("n_feature", "mt_percent"), function(category)
  {
    ggplot(mysco@meta, aes(x=mysco@meta[["n_count"]], y=mysco@meta[[category]])) +
      geom_point() +
      labs(x = "n_count", y = category)
  }))

  grid.arrange(grobs = plots, ncol = 3, nrow = 2)

}

FilterData = function(mysco, min.features, max.features, max.mt_percent, min.cells)
{
  library(dplyr)
  meta.filtered = mysco@meta %>% filter(n_feature > min.features & n_feature < max.features & mt_percent < max.mt_percent)
  counts.filtered = mysco@counts[, rownames(meta.filtered)]
  n_cells_per_feature = apply(counts.filtered > 0, 1, sum)
  counts.filtered = counts.filtered[names(which(n_cells_per_feature >= min.cells)), ]
  mysco@meta = meta.filtered
  mysco@counts = counts.filtered
  return(mysco)
}

NormaliseData = function(mysco, scaling_factor)
{
  #mysco@normalised = t((t(mysco@counts) / mysco@meta$n_count) * scaling_factor)
  mysco@normalised = log1p(apply(mysco@counts, 2, function(v) v / sum(v)) * scaling_factor)
  return(mysco)
}

FindHVGs = function(mysco, n_top)
{
  # calculate the logarithm of mean and variance
  feature_mean = apply(mysco@counts, 1, mean)
  feature_var = apply(mysco@counts, 1, var)
  feature_mean_log = log10(feature_mean)
  feature_var_log = log10(feature_var)
  # produce loess fitting
  fitting = loess(feature_var_log ~ feature_mean_log)
  # anti-log loess fitting
  fitting_anti_logged = 10 ** predict(fitting, feature_mean_log)
  # produce zscored matrix
  zscore = (mysco@counts - feature_mean) / (sqrt(fitting_anti_logged))
  # cap zscore to the square root of the number of cells
  sqrt_n_cells = sqrt(ncol(mysco@counts))
  zscore[which(zscore > sqrt_n_cells)] = sqrt_n_cells
  # calculate zscore variances
  zscore_var = (apply(zscore, 1, var))

  zscore_var_sorted = sort(zscore_var, decreasing = T)
  n_features = length(feature_mean)
  zscore_var_top = zscore_var_sorted[1:n_top]
  # plot feature means vs zscore variances
  plot(log(feature_mean), zscore_var)
  mysco@hvg = data.frame(gene_means = feature_mean[names(zscore_var_top)], hvgs_zscore_vars = zscore_var_top)
  return(mysco)
}

PlotHVGs = function(mysco)
{
  plot(log(mysco@hvg$gene_means), log(mysco@hvg$hvgs_zscore_vars))
}

ScaleData = function(mysco)
{
  mysco@data.scale = t(apply(mysco@normalised, 1, function(v) (v - mean(v)) / sd(v)))
  mysco@data.scale[mysco@data.scale > 10] = 10
  return(mysco)
}

CalcPC = function(mysco)
{
  library(RANN)
  # get the scaled counts of only the hvgs
  req.data = mysco@data.scale[rownames(mysco@hvg),]
  # calculate principal components
  mysco@principal_components = prcomp(t(req.data), center = F, rank. = 50)$x
  return(mysco)
}

ClusterData = function(mysco, n_pc, n_neighbours)
{
  req.data = mysco@data.scale[rownames(mysco@hvg),]
  # create nearest neighbour graph
  snn = nn2(mysco@principal_components[,1:n_pc], k=n_neighbours)$nn.idx
  # create adjacency matrix
  n_cells = nrow(mysco@principal_components)
  adjacency_matrix = matrix(0L, n_cells, n_cells)
  rownames(adjacency_matrix) = colnames(req.data)
  colnames(adjacency_matrix) = colnames(req.data)
  # populate adjacency matrix
  for (i in 1:n_cells)
  {
    adjacency_matrix[i, rownames(mysco@principal_components)[snn[i, ]]] = 1L
  }
  #sum(adjacency_matrix[1, ]) == n_neighbours
  #table(apply(adjacency_matrix, 1, sum))

  # perform leiden clustering
  library(igraph)
  mysco@louvain.clustering = cluster_louvain(graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected"))$membership
  return(mysco)
}

MakeUMAP = function(mysco)
{
  library(umap)
  umap.config = umap.defaults
  umap.config$n_components = 2

  umap.out = umap(mysco@principal_components[,1:10], config = umap.config)

  library(ggplot2)
  layout_data_frame = as.data.frame(umap.out$layout)
  layout_data_frame$group = mysco@louvain.clustering
  ggplot(layout_data_frame,
         aes(x = layout_data_frame[,2],
             y = layout_data_frame[,1])) +
    geom_point(aes(colour = factor(layout_data_frame$group)))
}
