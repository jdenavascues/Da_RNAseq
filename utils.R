# +-------------------------------------------------------------+
# |  1 Differential gene expression.                            |
# +-------------------------------------------------------------+




# +-------------------------------------------------------------+
# |  2 Visualisation ...                   |
# +-------------------------------------------------------------+




# +-------------------------------------------------------------+
# |  3 Visualisation ...                   |
# +-------------------------------------------------------------+




### +-------------------------------------------------------------+
### |  4 Visualisation: Gene Sets                                 |
### +-------------------------------------------------------------+

# extract genes with abs(log2FC)
extract_regulated_sets <- function(list_of_degs, names_degs, fc_thresh=1.5) {
  # fc_thresh must be a positive number
  # list_of_degs is a list of dataframes
  # names_degs is a list of strings
  # they must have the same length
  reg_lvl <- paste0(' reg@log~2~FC≥', fc_thresh)
  breaks <- brk_manual(c(-fc_thresh, fc_thresh), left_vec = c(FALSE, TRUE))
  regulated_sets <- NULL
  for (l in 1:length(list_of_degs)) {
    degs_na <- dplyr::select(list_of_degs[[l]], c(log2FoldChange, padj, ensemblGeneID))
    # better not pass NAs to kiru
    degs <- na.omit(degs_na)
    # filter by log2FC threshold
    degs$reg <- 
      kiru(
        degs$log2FoldChange,
        breaks = breaks,
        extend = TRUE,
        labels=c("down", "non", "up")
      )
    # filter by p-val
    degs$reg <- as.character(degs$reg)
    x <- 1:nrow(degs)
    degs$reg <- ifelse(degs[x,'padj']<0.5, degs[x,'reg'], 'non')
    # recover NAs as non-regulataed
    added_nas <- degs_na[!(rownames(degs_na) %in% rownames(na.omit(degs_na))),]
    added_nas$reg <- 'non'
    degs <- rbind(degs, added_nas)
    # store reg status
    #rownames(degs) <- degs$ensemblGeneID
    degs[ , names_degs[[l]] ] <- degs$reg
    regulated_sets[[l]] <- dplyr::select(degs, c(ensemblGeneID, names_degs[[l]]) )
  }
  regulated_sets <- purrr::reduce(regulated_sets, full_join, by='ensemblGeneID')
  rownames(regulated_sets) <- regulated_sets$ensemblGeneID
  regulated_sets <- dplyr::select(regulated_sets, -ensemblGeneID)
  return(regulated_sets)
}

get_deg_logical <- function(regs, direction){
  # direction must be 'up' or 'down'
  # regs is the output of `extract_regulated_sets`
  deg_logical_set <- (regs == direction)[
    ( rowSums(regs==direction)!=0 ), 
  ]
  return(deg_logical_set)
}



### +-------------------------------------------------------------+
### |  5 Visualisation with scatter plots (MA).                   |
### +-------------------------------------------------------------+


# classify lists of genes as up/non/down-regulated
extract_regulated_sets2 <- function(list_of_degs, names_degs, fc_thresh=1.5, cols) {
  # fc_thresh must be a positive number
  # list_of_degs is a list of dataframes
  # names_degs is a list of strings
  # they must have the same length
  reg_lvl <- paste0(' reg@log~2~FC≥', fc_thresh)
  breaks <- brk_manual(c(-fc_thresh, fc_thresh), left_vec = c(FALSE, TRUE))
  regulated_sets <- NULL
  for (l in 1:length(list_of_degs)) {
    degs_na <- dplyr::select(list_of_degs[[l]], all_of(cols))
    # better not pass NAs to kiru
    degs <- na.omit(degs_na)
    # filter by log2FC threshold, create `reg`(ulated) col with up/down/non
    degs$reg <- 
      kiru(
        degs$log2FoldChange,
        breaks = breaks,
        extend = TRUE,
        labels=c("down", "non", "up")
      )
    # filter by p-val
    degs$reg <- as.character(degs$reg) # remove factor
    x <- 1:nrow(degs)
    degs$reg <- ifelse(degs[x,'padj']<0.05, degs[x,'reg'], 'non')
    # recover NAs as non-regulataed
    added_nas <- degs_na[!(rownames(degs_na) %in% rownames(na.omit(degs_na))),]
    added_nas$reg <- 'non'
    degs <- rbind(degs, added_nas)
    # store reg status
    degs[ , names_degs[[l]] ] <- degs$reg
    regulated_sets[[l]] <- dplyr::select(degs, c(ensemblGeneID, names_degs[[l]]) )
  }
  regulated_sets <- purrr::reduce(regulated_sets, full_join, by='ensemblGeneID')
  rownames(regulated_sets) <- regulated_sets$ensemblGeneID
  regulated_sets <- dplyr::select(regulated_sets, -ensemblGeneID)
  return(regulated_sets)
}

# custom MA plot
ggmaplot2 <- function(deg, markers, fc_thresh = 1.5,
                      md_label='*genotype^—^*', repulsion) {
  # deg is the output of DESeq2 with an added column of gene symbols
  # gene.symbols is that column
  
  make_goilist_ggmaplot2 <- function(deg, fc_thresh, markerlist) {
    goi_list <- deg %>% dplyr::filter(gene_symbol %in% markerlist &
                                        abs(log2FoldChange) > fc_thresh &
                                        padj<0.05)
    return(goi_list)
  }
  
  # standard MA plot with `ggpubr`, using the CBD1 palette from `cetcolor` 
  ma <- ggmaplot(
    deg,
    fdr = 0.05,
    fc = fc_thresh,
    genenames = deg$gene_symbols,
    size = 2,
    alpha = 0.5,
    seed = NA,
    font.label = c(16, "bold", "black"),
    label.rectangle = FALSE,
    palette = c(cet_pal(n = 3, name = "cbd1", alpha = 1)[3], # "#A89008"
                cet_pal(n = 3, name = "cbd1", alpha = 1)[1], # "#3A90FE"
                "#AAAAAA"),
    top = 0,
    main = NULL,
    xlab = "log~2~(mean expression)",
    ylab = "log~2~(fold change)",
    ggtheme = theme_linedraw(),
    legend = 'top'
  )
  
  annotations <- data.frame(
    xpos = floor(min(ggplot_build(ma)$layout$panel_params[[1]]$x.range)),
    ypos = ceiling(max(abs(ggplot_build(ma)$layout$panel_params[[1]]$y.range))),
    annotateText = md_label,
    hjustvar = 0,
    vjustvar = 1)
  
  goi_list <- make_goilist_ggmaplot2(deg, fc_thresh, unlist(markers$xmarkers))
  
  ma <- ma +
     # mark special genes if upreg
     geom_text_repel(data = deg %>% filter(gene_symbol %in% goi_list$gene_symbol &
                                             log2FoldChange > 0),
                     aes(x = log2(baseMean),
                         y = log2FoldChange,
                         label = gene_symbol,
                         segment.square  = FALSE,
                         segment.inflect = TRUE),
                     color         = '#6E5C03', # ==hue, <ligthness than than "#A89008"
                     fontface = 'bold',
                     segment.alpha = 0.8,
                     segment.linetype = 3,
                     hjust = 0,
                     # positioning
                     box.padding   = repulsion$box.padding,
                     point.padding = repulsion$point.padding,
                     nudge_x       = repulsion$nudge_x.up,
                     nudge_y       = repulsion$nudge_y.up,
                     force         = repulsion$force,
                     force_pull    = repulsion$force_pull,
                     max.overlaps  = Inf,
                     xlim          = repulsion$xlims.up,    # NA repels from edges
                     ylim          = repulsion$ylims.up,
                     seed          = repulsion$seed.up) +
     # mark special genes if downreg
     geom_text_repel(data = deg %>% filter(gene_symbol %in% goi_list$gene_symbol &
                                             log2FoldChange < 0),
                     aes(x = log2(baseMean),
                         y = log2FoldChange,
                         label = gene_symbol,
                         segment.square  = FALSE,
                         segment.inflect = TRUE),
                     color         = '#0565AF', # ==hue, <ligthness than "#3A90FE"
                     fontface = 'bold',
                     segment.alpha = 0.8,
                     segment.linetype = 3,
                     hjust = 0,
                     # positioning
                     box.padding   = repulsion$box.padding,
                     point.padding = repulsion$point.padding,
                     nudge_x       = repulsion$nudge_x.dn,
                     nudge_y       = repulsion$nudge_y.dn,
                     force         = repulsion$force,
                     force_pull    = repulsion$force_pull,
                     max.overlaps  = Inf,
                     xlim          = repulsion$xlims.dn,             # NA repels from edges
                     ylim          = repulsion$ylims.dn,
                     seed          = repulsion$seed.dn) +
     # genotype label
     geom_richtext(
       data = annotations,
       aes(x = xpos,y = ypos,
           hjust = hjustvar, vjust = vjustvar,
           label = annotateText,
           fontface = 'bold',
           size = 9),
       fill = NA,
       label.color = NA,
       label.padding = grid::unit(rep(0, 4), "pt"),
       label.margin = grid::unit(rep(0, 4), "pt")
     ) +
     # legend
     guides(size = 'none') +
     # apply markdown formatting, etc
     theme(# markdown
       axis.title.x = element_markdown(size = 14, face = 'bold'),
       axis.title.y = element_markdown(size = 14, face = 'bold'),
       axis.text = element_text(size = 10),
       # grid and panel
       panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
       panel.grid = element_blank(),
       # legend
       legend.margin = margin(b = -10),
       legend.key.size = unit(1.2, 'cm'),
       legend.text = element_text(size = 12, face = 'bold')
       )
  return(ma)
}

# custom MA plot with labels coloured by cell type
ggmaplot3 <- function(deg, markers, fc_thresh = 1.5,
                      md_label='*genotype^—^*', repulsion) {
  # deg is the output of DESeq2 with an added column of gene symbols
  # gene.symbols is that column
  
  goi_list  <- tidyr::unnest_longer(markers, xmarkers) %>%
    rename(gene_symbol = xmarkers)
  
  # standard MA plot with `ggpubr`, using the CBD1 palette from `cetcolor` 
  ma <- ggmaplot(
    deg,
    fdr = 0.05,
    fc = fc_thresh,
    genenames = deg$gene_symbols,
    size = 2,
    alpha = 0.5,
    seed = NA,
    font.label = c(16, "bold", "black"),
    label.rectangle = FALSE,
    palette = c(cet_pal(n = 3, name = "cbd1", alpha = 1)[3], # "#A89008"
                cet_pal(n = 3, name = "cbd1", alpha = 1)[1], # "#3A90FE"
                "#AAAAAA"),
    top = 0,
    main = NULL,
    xlab = "log~2~(mean expression)",
    ylab = "log~2~(fold change)",
    ggtheme = theme_linedraw(),
    legend = 'top'
  )
  
  annotations <- data.frame(
    xpos = floor(min(ggplot_build(ma)$layout$panel_params[[1]]$x.range)),
    ypos = ceiling(max(abs(ggplot_build(ma)$layout$panel_params[[1]]$y.range))),
    annotateText = md_label,
    hjustvar = 0,
    vjustvar = 1)
  
  genlabsup <- deg %>%
    filter(gene_symbol %in% goi_list$gene_symbol &
             log2FoldChange > fc_thresh &
             padj < 0.05) %>%
    dplyr::select(c(gene_symbol, baseMean, log2FoldChange)) %>%
    dplyr::left_join(dplyr::select(goi_list, !'celltype'), by = 'gene_symbol')
  
  genlabsdn <- deg %>%
    filter(gene_symbol %in% goi_list$gene_symbol &
             log2FoldChange < -fc_thresh &
             padj < 0.05) %>%
    dplyr::select(c(gene_symbol, baseMean, log2FoldChange)) %>%
    dplyr::left_join(dplyr::select(goi_list, !'celltype'), by = 'gene_symbol')
  
  ma <- ma +
    # mark special genes if upreg
    geom_text_repel(data = genlabsup,
                    aes(x = log2(baseMean),
                        y = log2FoldChange,
                        label = gene_symbol,
                        segment.square  = FALSE,
                        segment.inflect = TRUE),
                    colour = genlabsup$cellcolour,
                    segment.color = '#6E5C03', # ==hue, <ligthness than than "#A89008"
                    fontface = 'bold',
                    segment.alpha = 0.8,
                    segment.linetype = 3,
                    hjust = 0,
                    # positioning
                    box.padding   = repulsion$box.padding,
                    point.padding = repulsion$point.padding,
                    nudge_x       = repulsion$nudge_x.up,
                    nudge_y       = repulsion$nudge_y.up,
                    force         = repulsion$force,
                    force_pull    = repulsion$force_pull,
                    max.overlaps  = Inf,
                    xlim          = repulsion$xlims.up,    # NA repels from edges
                    ylim          = repulsion$ylims.up,
                    seed = repulsion$seed.up) +
    # mark special genes if downreg
    geom_text_repel(data = genlabsdn,
                    aes(x = log2(baseMean),
                        y = log2FoldChange,
                        label = gene_symbol,
                        segment.square  = FALSE,
                        segment.inflect = TRUE),
                    colour = genlabsdn$cellcolour,
                    segment.color = '#0565AF', # ==hue, <ligthness than "#3A90FE"
                    fontface = 'bold',
                    segment.alpha = 0.8,
                    segment.linetype = 3,
                    hjust = 0,
                    # positioning
                    box.padding   = repulsion$box.padding,
                    point.padding = repulsion$point.padding,
                    nudge_x       = repulsion$nudge_x.dn,
                    nudge_y       = repulsion$nudge_y.dn,
                    force         = repulsion$force,
                    force_pull    = repulsion$force_pull,
                    max.overlaps  = Inf,
                    xlim          = repulsion$xlims.dn,             # NA repels from edges
                    ylim          = repulsion$ylims.dn,
                    seed = repulsion$seed.dn) +
    # genotype label
    geom_richtext(
      data = annotations,
      aes(x = xpos,y = ypos,
          hjust = hjustvar, vjust = vjustvar,
          label = annotateText,
          fontface = 'bold',
          size = 9),
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt"),
      label.margin = grid::unit(rep(0, 4), "pt")
    ) +
    # legend
    guides(size = 'none') +
    # apply markdown formatting, etc
    theme(# markdown
      axis.title.x = element_markdown(size = 14, face = 'bold'),
      axis.title.y = element_markdown(size = 14, face = 'bold'),
      axis.text = element_text(size = 10),
      # grid and panel
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      # legend
      legend.margin = margin(b = -10),
      legend.key.size = unit(1.2, 'cm'),
      legend.text = element_text(size = 12, face = 'bold')
    )
  return(ma)
}


### +-------------------------------------------------------------+
### |  6-1 Gene Set over-representation/enrichment.               |
### +-------------------------------------------------------------+


# To get a differentially expressed gene set for Over-Representation Analysis
make_degset <- function(deg, up, fc_thresh) {
  converter <- ifelse(up,1,-1)
  degset <- deg %>%
    filter(padj<0.05) %>%
    filter(abs(log2FoldChange) > fc_thresh) %>%
    filter(log2FoldChange*converter > 0)
  return(degset)
}

# To get a ranked gene list for Gene Set Enrichment Analysis
make_degrank <- function(deg, mode='log2fc', key='gene_symbol') {
  if (!(key %in% c('gene_symbol', 'ensemblGeneID'))) {
    stop("The `key` parameter needs to be one of 'gene_symbol', 'ensemblGeneID'")
  }
  # mode: whether to use padj or log2FC
  if (mode=='log2fc') {
    deg_sorted <- deg %>% arrange(desc(log2FoldChange))
    degrank <- deg_sorted$log2FoldChange
    names(degrank) <- deg_sorted[,key]
  } else if (mode=='padj') {
    deg$p.rank <- sign(deg$log2FoldChange) * -log10(deg$padj)
    deg_sorted <- deg %>% arrange(desc(p.rank))
    degrank <- deg_sorted$p.rank
    names(degrank) <- deg_sorted[,key]
  } else {
    stop('the `mode` parameter can only be either `log2fc` or `padj`')
  }
  return(degrank)
}

# To import the gene sets in the gmx files
import_from_gmx <- function(gmxfile) {
  # read a GMX file and turn it into a df input appropriate for `clusterProfiler`
  df <- read.csv(gmxfile, header=TRUE, sep='\t')
  # remove 'na's in 1st row
  df <- df[ !df[1]=='na', ]
  # pivot
  df <- df %>% pivot_longer(cols = everything(),
                            names_to = 'term',
                            values_to = 'gene')
  # remove rows with no genes
  df <- df %>% subset(gene != '')
  df <- df %>% mutate(gene = str_replace(gene, "FBGN", "FBgn"),
                      term = str_replace(term, ".only", "-only"))
  df <- df %>% mutate(term = str_replace(term, "\\.", " "))
  df <- df %>% mutate(term = str_replace(term, " genes", ""))
  return(df)
}

gseCP_summarise <- function(gmx, gseCP_list, conditions, sets.as.factors, cluster=FALSE, nsig.out=FALSE) {
  # purr::map to convert the gse_list from S4 objects to their @result slots 
  gseCP_list <- map( gseCP_list, \(x) dplyr::select(x@result, NES, p.adjust, ID) )
  # name them to associate conditions with the data
  # add condition as an extra column
  gseCP_list <- lapply( 1:length(gseCP_list), \(x) cbind(gseCP_list[[x]],
                                                         condition = conditions[[x]]) )
  df <- bind_rows(gseCP_list)
  df$condition <- factor(df$condition, levels = conditions)
  df$ID <- factor(df$ID, levels = sets.as.factors)
  # apply clustering
  if (cluster & length(sets.as.factors)>2) { # `hclust` must have n >= 2 objects to cluster
    vectors <- df %>%
      dplyr::select(ID, NES, condition) %>%
      pivot_wider(names_from = condition, values_from = NES) %>%
      mutate_all(replace_na, 100)
    clust <- hclust(dist(as.matrix(vectors[2:length(gseCP_list)])))
    df$ID <- factor(df$ID, levels=vectors$ID[clust$order])
  } else if (cluster & length(sets.as.factors < 3)) {
    warning("`hclust` must have n >= 2 objects to cluster. NES columns will not be clustered.")
  }
  # apply filtering
  if (nsig.out) {
      # get the IDs for which there is at least one condition with significant enrichment
      sign_lgl <- lapply( lapply(unique(df$ID),
                                 \(x) filter(df, ID==x)),
                          \(x) any(x$p.adjust<0.05) )
      df_bycondition <- lapply(unique(df$ID),
                               \(x) filter(df, ID==x) )[unlist(sign_lgl)]
      df <- bind_rows(df_bycondition)
      # to avoid passing on non-filtered terms
      df$ID <- as.character(df$ID)
      levels(df$ID) <- factor(unique(df$ID))
  }
  return(df)
}

## Extract data from a list of gse objects from clusterProfiler into a single df
## gmx is a df containing the data of a GMX file
## gseCP_list is a list of clusterProfiler GSEA objects
## conditions is a list with the names of the experimental conditions for those GSEA objects
## dvar (data variable) for "name injection" with <data-masking> tidyverse functions
# gseCP_summarise_old <- function(gmx, gseCP_list, conditions, sets.as.factors, dvar, cluster=FALSE, nsig.out=FALSE) {
#   if ( !(dvar %in% c('NES', 'p.adjust')) ) {
#     stop("The argument `dvar` must be one of `NES` and `p.adjust`")
#   }
#   # create template with all term IDs (clusterProfiler trims those that do not give results)
#   empty <- data.frame(rep(0, length(unique(gmx$term))),
#                       unique(gmx$term))
#   names(empty) <- c(dvar, 'ID')
#   # create function to merge `empty` with GSE object data for each term ID
#   gsemerge <- function(gse, empty, dvar) {
#     return(
#       # full join selecting injecting dvar with
#       # embracing syntax {{}}
#       # bangbang !! operator with as.name function
#       # glue syntax "{}"
#       # dynamic assignment :=
#       full_join(dplyr::select(gse, all_of({{dvar}}), ID), empty, by='ID') %>%
#         mutate( "{dvar}" := !!as.name(paste0(dvar,".x")) ) %>%
#         dplyr::select( {{dvar}}, ID )
#     )
#   }
#   # purr::map to convert the gse_list from S4 objects to their @result slots 
#   gseCP_list <- map(gseCP_list, \(x) x@result)
#   # apply gsemerge
#   gseCP_list <- lapply( gseCP_list, \(x) gsemerge(x, empty, dvar) )
#   # name them to associate conditions with the data
#   names(gseCP_list) <- conditions
#   df <- gseCP_list %>%
#     # change the colnames to condition names
#     imap(.x = ., ~ set_names(.x, c(.y, "ID"))) %>%
#     # merge them all together
#     purrr::reduce(full_join, by='ID') %>%
#     # tidy up all dvar values in one col
#     pivot_longer(cols=-c('ID'), names_to = "condition", values_to = dvar)
#   df$condition <- factor(df$condition, levels = conditions)
#   df$ID <- factor(df$ID, levels = sets.as.factors)
#   # apply clustering
#   if (cluster & length(sets.as.factors)>2) { # `hclust` must have n >= 2 objects to cluster
#     if (dvar != 'NES') warning("Clustering is intended for enrichment scores, not p-values")
#     vectors <- df %>%
#       pivot_wider(names_from = condition, values_from = {{dvar}}) %>%
#       mutate_all(replace_na, 100)
#     clust <- hclust(dist(as.matrix(vectors[2:length(gseCP_list)])))
#     df$ID <- factor(df$ID, levels=vectors$ID[clust$order])
#   } else if (cluster & length(sets.as.factors < 3)) {
#     warning("`hclust` must have n >= 2 objects to cluster. NES columns will not be clustered.")
#   }
#   # apply filtering
#   if (nsig.out) {
#     if (dvar=='p.adjust') {
#       # get the IDs for which there is at least one condition with significant enrichment
#       sign_lgl <- lapply( lapply(unique(df$ID),
#                                  \(x) filter(df, ID==x)),
#                           \(x) any(x$p.adjust<0.05) )
#       df_bycondition <- lapply(unique(df$ID),
#                                \(x) filter(df, ID==x) )[unlist(sign_lgl)]
#       df <- bind_rows(df_bycondition)
#     }
#     if (dvar=='NES') {
#       
#     }
#     
#   }
#   return(df)
# }

# create a heatmap with NES in colour and coloured point as padj.
layer.heatmap <- function(layerhm.df, subt) {
  NES.df <- layerhm.df %>% dplyr::select(-p.adjust)
  padj.df <- layerhm.df %>% dplyr::select(-NES)
  p <- ggplot(NES.df, aes(x=ID, y=condition)) +
    # plot statistic (NES)
    geom_tile(aes(fill=NES), width=1) +
    scale_fill_gradient2(low = cet_pal(3, name='cbd1')[1],
                         mid = cet_pal(3, name='cbd1')[2],
                         high = cet_pal(3, name='cbd1')[3],
                         midpoint = 0) +
    # plot p-value
    # statistically significant with colour and shape `*`
    geom_point(data=subset(padj.df, p.adjust<0.05),
               aes(x=ID, y=condition, colour=-log10(p.adjust)),
               size=3, shape=8, stroke=1.5, alpha=1) +
    #size=5, shape=23, alpha=1) +
    # non-significant in dark gray and shape `x`
    # geom_point(data=subset(padj.df, p.adjust>0.05),
    #            aes(x=ID, y=condition),
    #            size=3, shape=4, stroke=2, alpha=1, colour='gray50') +
    scale_colour_gradient(low = cet_pal(3, name='d2')[2],
                          high = cet_pal(3, name='d2')[3]) +
    coord_equal() +
    theme_bw() +
    ggtitle("Normalised Enrichment Scores",
            subtitle = subt) +
    scale_x_discrete(labels=toupper(levels(NES.df$ID)),
                     position = "top",
                     expand = expansion(mult = 0, add = 0)) +
    scale_y_discrete(expand = expansion(mult = 0, add = 0)) +
    labs(fill = 'NES',
         colour = '-log~10~(_p.adj_)<br><span style = "font-size:8pt;">_p.adj_<0.05</span>') +
    theme(plot.subtitle = element_markdown(),
          axis.text.x = element_markdown(angle=25, hjust=0,
                                         face='bold', size=10.5),
          axis.text.y = element_markdown(hjust=1, face='bold', size=12),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_markdown(hjust=0.5, vjust=0.75),
          legend.direction = 'horizontal',
          legend.position = 'bottom',
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 1))
  return(p)
}

# For GLAD sub.(sub.)group terms
refine_glad_by <- function(dataset, col.name) {
  # refined by col.name
  subgrouping <- dataset %>%
    # select terms with subgroups
    filter(!is.na(.data[[col.name]])) %>%
    # select the term and subgroup cols
    dplyr::select(term, {{col.name}}) %>%
    # get unique pairs of terms/subgroup
    distinct()
  # re-make term/gene cols as per subgroups now
  gmx <- lapply(
    unique(subgrouping$term),
    \(x) dataset %>%
      filter( !is.na(.data[[col.name]]) & term=={{x}} ) %>%
      dplyr::select({{col.name}}, FBgn) %>%
      rename(gene = FBgn, term = {{col.name}})
  )
  names(gmx) <- unique(subgrouping$term)
  # remove terms that only have one subgroup
  gmx <- gmx[ sapply(gmx, \(x) length(unique(x$term)))>1 ]
  return(gmx)
}

subglad_gsea <- function(deg, gmx, perms=1000) {
  g <- lapply(
    1:length(gmx), \(x) do.call(
      GSEA,
      c(list(geneList=make_degrank(deg, mode='log2fc', key='ensemblGeneID'),
             nPermSimple = perms,
             TERM2GENE = gmx[[x]]),
        GSEAparams)
    )
  )
  names(g) <- names(gmx)
  return(g)
}


###
### 6-2
###




###
### 6-3
###

handled_pathview <- function(params) {
  withCallingHandlers(
    do.call ( pathview, c(params, kegg.native = FALSE) ),
    message = function(m) {
      if( length(grep('Try \"kegg.native=T\" instead!', m$message))==1) {
        cat('OK, I will try that\n')
        do.call ( pathview, c(params, kegg.native = TRUE) )
        cat('now it is done with native KEGG PNG\n>>>>> Ignore what comes below: <<<<<\n')
      }
    }
  )
}


###
### 7
###