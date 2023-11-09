setClass("mysco",slots= c(meta="data.frame",counts="matrix"))

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
  n_cells_per_feature = apply(mysco@counts > 0, 1, sum)
  counts.filtered = mysco@counts[names(which(n_cells_per_feature >= min.cells)), rownames(meta.filtered)]
  mysco@meta = meta.filtered
  mysco@counts = counts.filtered
  return(mysco)
}
