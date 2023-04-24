###
### 1
###




###
### 2
###




###
### 3
###




###
### 4
###




###
### 5
###




###
### 6-1 Gene Set over-representation/enrichment
###

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
  }
  return(df)
}

## Extract data from a list of gse objects from clusterProfiler into a single df
## gmx is a df containing the data of a GMX file
## gseCP_list is a list of clusterProfiler GSEA objects
## conditions is a list with the names of the experimental conditions for those GSEA objects
## dvar (data variable) for "name injection" with <data-masking> tidyverse functions
# gseCP_summarise <- function(gmx, gseCP_list, conditions, sets.as.factors, dvar, cluster=FALSE, nsig.out=FALSE) {
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