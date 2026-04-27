################################################################################
# Import BIOM tables and merge samples into a single phyloseq object for each 
# taxonomic rank (P, C, O, F, G, S), and produce relative abundance barplots 
# per taxonomic ranks
################################################################################

library(dplyr)
library(ggplot2)
library(phyloseq)

####################### Parameters used to run MetNano #########################
# Run name
run = "example_run"

# Read mininum length
min_len = 300

# Read minimum quality
min_qual = 15

# Kraken2 database
database = "k2_core_nt_20251015"

# Kraken2 confidence score
confidence = 0.1

filter_tag = sprintf("min_%sbp_%sQ", min_len, min_qual)

############################# Formatting function ##############################
format_tax_level = function(tax_rank) {
  taxonomic_level = substr(tax_rank, start = 1, stop = 1)
  path = file.path(
    run, filter_tag, database, confidence, taxonomic_level
  )

  # Retrieve list of biom tables for the specified taxonomic level
  biom_tables = list.files(path)

  # Import biom tables into a single phyloseq object ---------------------------
  for (i in 1:length(bioms)){
    biom = bioms[i]
    biom_path = file.path(
      path, biom
    )
    metagenome = import_biom(biom_path)
    colnames(metagenome@otu_table) = paste0(colnames(metagenome@otu_table),
                                            "_", i)
    sample_data(metagenome) = data.frame(row.names = sample_names(metagenome),
                                         ID = sample_names(metagenome),
                                         Run = run)
    if (i == 1){
      ps = metagenome
    } else {
      ps = merge_phyloseq(ps, metagenome)
    }
  }
  metagenome = ps

  ## Clean taxonomy ------------------------------------------------------------
  tax_table(metagenome) = tax_table(metagenome) %>%
    as.matrix() %>%
    substring(4) %>%
    `colnames<-`(c("Kingdom","Phylum","Class","Order",
                   "Family","Genus","Species"))

  ## Remove unwanted taxa
  metagenome = subset_taxa(metagenome, Species != "uncultured bacterium")

  ## Filter data ---------------------------------------------------------------
  # Remove unwanted taxa
  metagenome = subset_taxa(metagenome, Phylum != "Chordata")
  
  taxa_tab = as(tax_table(metagenome), "matrix")
  keep = taxa_tab[, tax_rank] != ""
  ps_filt = prune_taxa(keep, metagenome)
  ps_percent = transform_sample_counts(
    ps_filt,
    function(x) 100 * x / sum(x)
  )
  ps_glom = tax_glom(ps_percent, taxrank = tax_rank)
  df = psmelt(ps_glom)
  df$Sample = unlist(lapply(strsplit(df$Sample, "_"), "[", 1))
  df$ID = unlist(lapply(strsplit(df$ID, "_"), "[", 1))
  return(df)
}

#################### Generate relative abundance dataframes ####################
df_P = format_tax_level("Phylum")

df_C = format_tax_level("Class")

df_O = format_tax_level("Order")

df_F = format_tax_level("Family")

df_G = format_tax_level("Genus")

df_S = format_tax_level("Species")

########################## General plotting function ###########################
plot_tax_level = function(ps_obj,
                          tax_rank,
                          low_abundance_threshold,
                          output_name) {
  taxa_tab <- as(tax_table(ps_obj), "matrix")
  keep <- taxa_tab[, tax_rank] != ""
  ps_filt <- prune_taxa(keep, ps_obj)
  ps_percent <- transform_sample_counts(
    ps_filt,
    function(x) 100 * x / sum(x)
  )
  ps_glom = tax_glom(ps_percent, taxrank = tax_rank)
  df = psmelt(ps_glom)
  df$Run = run

  # Group low abundance taxa ---------------------------------------------------
  df[[tax_rank]] = as.character(df[[tax_rank]])
  df[[tax_rank]][df$Abundance < low_abundance_threshold] =
    paste0(tax_rank, " < ", low_abundance_threshold, "%", " abundance")

  df[[tax_rank]] = factor(df[[tax_rank]])

  # Create color palette -------------------------------------------------------
  colors = colorRampPalette(
    brewer.pal(8, "Dark2")
  )(length(levels(df[[tax_rank]])))

  # Create plot ----------------------------------------------------------------
  p = ggplot(df, aes(x = Sample,
                      y = Abundance,
                      fill = .data[[tax_rank]])) +
    geom_bar(stat = "identity") +
    facet_wrap(~Run, scales = "free_x") +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      plot.title = element_text(size = 10, hjust = 0.5)
    ) +
    labs(x = "Samples",
         y = "Relative abundance (%)",
         fill = tax_rank,
         title = plot_title)
  
  # Save plot ------------------------------------------------------------------
  ggsave(
    filename = file.path("results", run, output_name),
    plot = p,
    width = 12,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

################################ Generate plots ################################
# plot_tax_level() takes 4 arguments: 
#   - a phyloseq object
#   - the studied taxonomic level
#   - the relative abundance threshold for grouping low abundance taxa
#   - a name for the output JPEG file

plot_P = plot_tax_level(df_P, "Phylum", 0.5, "relative_abundance_barplot_P.jpg")

plot_C = plot_tax_level(df_C, "Class", 1, "relative_abundance_barplot_C.jpg")

plot_O = plot_tax_level(df_O, "Order", 2, "relative_abundance_barplot_O.jpg")

plot_F = plot_tax_level(df_F, "Family", 3, "relative_abundance_barplot_F.jpg")

plot_G = plot_tax_level(df_G, "Genus", 4, "relative_abundance_barplot_G.jpg")

plot_S = plot_tax_level(df_S, "Species", 5, "relative_abundance_barplot_S.jpg")
