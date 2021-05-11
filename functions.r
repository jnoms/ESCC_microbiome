#!/bin/R 

# This script contains shared functions, and is not meant to be executed.

# For use as a hashmap replacement...
hash_lookup <- function(in_hash, in_vector){
  # This function allows replacement of the hashmap object with similar behavior from the 
  # hash object. in_hash should be generated as such:
  #     in_hash <- hash(vector_of_keys, vector_of_values)
  # In this function, in_vector is a vector of the keys that you desire to be looked up.
  # The output is a vector of same length as the in_vector, but containing the values.
  # If an input is not found in the hash, NA will be returned as value.
  
  result <- c()
  for (key in in_vector) {
    value = in_hash[[key]]
    
    # If item can't be found, value is NA
    if (is.null(value)) {
      value = NA
    }
    
    result = c(result, value)
  }
  result
}



# Sort tibble for a superkingdom and level
extract_data <- function(in_tbl, superkingdom, level, remove_cols) {
  
  # Filter the unneeded cols for only those that are present
  remove_cols <- remove_cols[remove_cols %in% colnames(in_tbl)]
  
  # This function sorts tbl for a superkingdom, level, and removes 
  # specified columns.
  out_tbl <- in_tbl %>%
    filter(kingdom == superkingdom, type == level) %>%
    dplyr::select(-remove_cols)
  out_tbl
}

format_metrics <- function(metrics_infile, hash) {
  # Takes in the initial metrics input tibble and makes a tidy tibble that is transposed so it's easier to 
  # deal with. If a hash is input, will also convert sampleID's to patientIDs.
  
  # Format the colnames
  colnames(metrics_infile) <- str_replace_all(colnames(metrics_infile), ".pathseq.filter_metrics", "")
  
  # Keep only the columns where the sampleID was successfully translated to patient ID
  metrics_infile <- subset(metrics_infile, select=colnames(metrics_infile)[!is.na(colnames(metrics_infile))])
  
  # Transpose and make colnames better
  metrics_infile_transposed <- as.data.frame(t(metrics_infile))
  colnames(metrics_infile_transposed) <- t(metrics_infile)[1,]
  metrics_infile_transposed <- metrics_infile_transposed[-1,]
  
  # Move the rownames to first row and make tibble
  metrics_infile_transposed$name <- row.names(metrics_infile_transposed)
  output <- as.tibble(metrics_infile_transposed) %>%
    dplyr::select(name, everything())
  
  output
}

# Get relative abundance
relative_abundance <- function(in_tbl) {
  # The input to this functino is a tbbl with the column "name"
  # detailing taxa, and all the other columns sample counts.
  out <- as_tibble(apply(in_tbl[,-1], 2, function(i) i/sum(i))) %>%
    add_column(name = in_tbl$name) %>%
    dplyr::select(name, everything())
  
  # return output
  out
}

get_jaccard <- function(tbl_1, tbl_2, sample_name, threshold=0.01) {
  # This function find the jaccard index between a given, common
  # sample_name in two tables. The threshold is a numeric value 
  # n, where only taxa with abundance greater than n are counted.
  
  # Get the column with the name and counts
  tbl1_sample_col <- tbl_1 %>% dplyr::select(name, sample_name)
  tbl2_sample_col <- tbl_2 %>% dplyr::select(name, sample_name)
  
  # get names of taxa with non-zero counts
  tbl_1_names <- tbl1_sample_col[tbl1_sample_col[2] > threshold,]$name
  tbl_2_names <- tbl2_sample_col[tbl2_sample_col[2] > threshold,]$name
  
  # Get shared names
  shared_names <- intersect(tbl_1_names, tbl_2_names)
  
  # Find total names in either
  total_names <- c(tbl_1_names, tbl_2_names)
  total_names <- total_names[!duplicated(total_names)]
  
  # Divide shared over total - resultant value is jaccard index
  jaccard_index <- length(shared_names)/length(total_names)
  jaccard_index
}

# Get the top n rows, by rowMean
get_top_by_rowMeans <- function(in_tbl, n_rows){
  # Assumes the in_tbl is properly formatted - The first
  # column should the taxa names, or whatever. ALL OTHER
  # COLUMNS should be numeric values.
  # 
  # This function works by calculating the rowMean for each
  # row and taking the top n_rows (specified). It then returns
  # the top rows. 
  
  # Make a new column that has row means
  in_tbl$means <- rowMeans(in_tbl[2:length(in_tbl)])
  
  # Sort, keep top N rows, and remove means column
  out_tbl <- in_tbl %>% arrange(desc(means))%>%
    top_n(n_rows) %>%
    dplyr::select(-means)
  
  # return output
  out_tbl
}

name_col_to_rownames <- function(in_tbl){
  # The purpose of this function is to take in a dataframe
  # and make the column titled "name" into the row names, 
  # and then delete that column... Helpful to use pheatmap
  # on a tibble.
  
  # Conver to dataframe
  out_tbl <- as.data.frame(in_tbl)
  
  # Reset rownames and dump name col
  rownames(out_tbl) <- out_tbl$name
  out_tbl$name <- NULL
  
  # output...
  out_tbl
}

isolate_taxa <- function(in_tbl, taxa) {
  # This function takes in an input tibble whose
  # values are RELATIVE ABUNDANCE. IT also takes
  # in taxa, which is a vector detailing a list
  # of taxa which are wanted. These taxa are
  # specifically isolated, and an "Other Taxa"
  # row is added that makes the columns sum to 1.
  
  # Filter for the desired taxa
  out_tbl <- in_tbl %>%
    filter(name %in% taxa)
  
  # Get the remainder
  tmp <- as.vector(1-colSums(out_tbl[2:length(out_tbl)]))
  tmp <- c("Other Taxa", tmp)
  
  # Add remainder as a new row
  out_tbl <- rbind(out_tbl, tmp)
  
  # name numeric
  out_tbl[2:length(out_tbl)] <- sapply(out_tbl[2:length(out_tbl)], as.numeric)
  
  # output
  out_tbl
}

make_taxon_vs_other <- function(in_tbl, taxon){
  # This is similar to the function 'isolate_taxa', except
  # it doesn't assume relative abundance. It takes one taxa,
  # then reports that as well as the sum of all the others.
  
  # Start output tibble
  out <- in_tbl %>%
    filter(name == taxon)
  
  # Remove the taxa frorm the initial tibble
  in_tbl_modified <- in_tbl %>%
    filter(!name == taxon)
  
  # Get colSums of the modified tbl
  sums <- c("Other Taxa", colSums(in_tbl_modified[2:length(in_tbl_modified)]))
  
  # Add colsums to the new tibble
  out <- rbind(out, sums)
  
  # Make it numeric
  out[2:length(out)] <- sapply(out[2:length(out)], as.numeric)
  
  return(out)
}

get_order <- function(in_tbl, in_tbl.m, taxon) {
  # This function is for deciding the ORDER OF THE COLUMNS
  # for plotting. in_tbl is the tibble before melting, while
  # in_tbl.m is the tibble after melting. The taxon is the
  # one you want the order for. For example, if taxon="Fusobacterium"
  # the output order will specify the orer of in_tbl.m such that
  # the bars are structured from more to less fusobacterium.
  
  # Transpose the tibble
  in_tbl.t <- as_tibble(cbind(nms = names(in_tbl), t(in_tbl)))
  colnames(in_tbl.t) <- in_tbl.t[1,]
  in_tbl.t <- in_tbl.t[-1,]
  
  # Take only the desired taxon, and make sure
  # its numeric.
  in_tbl.t.taxon <- in_tbl.t %>%
    dplyr::select(name, matches(taxon))
  in_tbl.t.taxon[,2] <- as.numeric(unlist(in_tbl.t.taxon[,2]))
  
  # Order based on the taxon, and select just the
  # name
  name_list_tbl <- in_tbl.t.taxon %>%
    arrange_at(vars(contains(taxon))) %>%
    dplyr::select(name)
  
  # Need to reverse the list to make it go from more-->less
  name_list <- rev(name_list_tbl$name)
  
  # In the input tbl, find position of each sample. 
  # This is the output order.
  order <- match(in_tbl.m$Sample, name_list)
  order
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_svg <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  svglite(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

normalize_to_human_reads <- function(in_tbl, hu_reads_hashmap, number_of_nonname_columns=1, div_constant=1000000){
  # Normalized input tibble to # human reads
  
  # Inputs:
  # - in_tbl:                    Input tibble (no rownames) that has sample values as columns. 
  #                              There can be non-sample-count columns at the front, the
  #                              number of which are specified by number_of_nonname_columns.
  #
  # - hu_reads_hashmap:          HASH that converts a sample name (present in the
  #                              columns of the in_tbl) to # human reads. Not a hashmap, a hash object.
  #
  # - number_of_nonname_columns: Specifies the # of columns, at the start of 
  #                              the in_tbl, do NOT contain sample values.
  #
  # - div_constant:              Constant # of reads to multiply by - i.e. if 1,000,000, 
  #                              the resultant values from this function are # reads 
  #                              per million human reads.
  
  result <- in_tbl
  names_start_at <- number_of_nonname_columns + 1
  names <- colnames(in_tbl)[names_start_at:length(colnames(in_tbl))]
  
  for (name in names) {
    
    # Get # human reads
    hu_reads <- hash_lookup(hu_reads_hashmap, name)
    
    # Extract the values for the input name
    current_counts <- in_tbl[name]
    
    # Divide by human reads
    new_counts <- (current_counts * div_constant)/hu_reads
    
    # Reset the input column in the result
    result[name] <- new_counts
  }
  
  # return result
  result
}

make_factor_numeric <- function(in_vector) {
  as.numeric(as.character(in_vector))
}

metrics_extract_wanted_cols_as_numeric <- function(metrics_tibble) {
  metrics_tibble  %>%
    dplyr::select(name, HOST_READS_FILTERED, MAPPED_READS, UNMAPPED_READS) %>%
    mutate(HOST_READS_FILTERED = as.numeric(as.character(HOST_READS_FILTERED))) %>%
    mutate(MAPPED_READS = as.numeric(as.character(MAPPED_READS))) %>%
    mutate(UNMAPPED_READS = as.numeric(as.character(UNMAPPED_READS)))
}


#--------------------------------------------------------------------------------------------------------#
# PLOTTING BACTERIA  ####
#--------------------------------------------------------------------------------------------------------#
bacteral_taxa_barchart <- function(in_tbl,
                                   sort_by_taxon="Fusobacterium",
                                   taxa_of_interest= c("Selenomonas", "Prevotella", "Campylobacter", "Porphyromonas", "Streptococcus", "Fusobacterium"),
                                   specified_sample_order=FALSE,
                                   unneeded_cols = c("tax_id", "taxonomy", "type", "kingdom", "reference_length", "mean", "median", "max"),
                                   color_codes = c("#F4F6F6", "#C39BD3", "#E59866", "#F7DC6F", "#7DCEA0", "#85C1E9", "#225ea8")
) {
  
  # PLOTS a bacterial relative abundance barchart. The output barchart will
  # have the columns sorted based on the abundance of the sort_by_taxon, OR
  # based on the order specified by specified_sample_order. 
  
  # Determine which of the unneeded cols are in the input tibble
  unneeded_cols <- unneeded_cols[unneeded_cols %in% colnames(in_tbl)]
  
  # Sort for bacterial genera, make relative abundance
  bact_genera <- extract_data(in_tbl, 'Bacteria', 'genus', unneeded_cols)
  bact_genera <- relative_abundance(bact_genera)
  
  # Sort for taxa of interest
  bact_genera.subset <- isolate_taxa(bact_genera, taxa_of_interest)
  
  # If specified_sample_order input, remove samples that aren't present in the desired order
  if (specified_sample_order != FALSE) {
    bact_genera.subset <- bact_genera.subset %>%
      dplyr::select(name, specified_sample_order)
  }
  
  # Melt for plotting
  bact_genera.subset.m <- melt(bact_genera.subset, id="name")
  colnames(bact_genera.subset.m) <- c("Genera", "Sample", "rel_abundance")
  
  # Plotting
  #-----------------------------------------------------------------------------#
  
  # Define order of the columns
  if (specified_sample_order == FALSE) {
    
    # If specified_sample_order isn't specified, order based on taxa_of_interest
    order <- get_order(bact_genera.subset, bact_genera.subset.m, sort_by_taxon)
  } else {
    
    # Otherwise, order by the specified sample order
    order <- match(bact_genera.subset.m$Sample, specified_sample_order)
  }
  
  # Add order
  bact_genera.subset.m$order <- order
  
  # HERE - DEFINING THE ORDER OF THE STACKED BAR - THE LAST ON THE LIST IS TOWARDS THE BOTTOM
  levels  <- c("Other Taxa", taxa_of_interest)
  
  # Plotting
  ggplot(bact_genera.subset.m,
         aes(x=reorder(Sample, order),
             y=rel_abundance,
             fill=factor(Genera, levels=levels, ordered=T)
         )
  ) +
    geom_bar(stat='identity', color="white", width=1) + #color=white sets line between bars
    theme_classic(base_size = 13) + # Sets theme, and makes all text 22pt
    theme(axis.title.x=element_blank(), # These element_blank()s remove x axis labels
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(hjust = 0.5)
    ) +
    labs(fill = "Genus") + #rename the legend title
    ylab("Relative Abundance") + #rename y axis title
    scale_y_continuous(expand = c(0, 0)) + #Removes the white space bewteen graph and the bottom black line
    scale_fill_manual(values = color_codes)
}

rel_abundance_pheatmap <- function(in_tbl,
                                   top_n,
                                   unneeded_cols = c("tax_id", "taxonomy", "type", "kingdom", "reference_length", "mean", "median", "max"),
                                   output_plot_ready_data = FALSE
) {
  
  # This function takes in a standard pathseq tibble and outputs either 
  # the bacterial genera pheatmap, with the top_n genera by rowmean
  # labeled, or the processed dataframe ready for plotting (if 
  # output_plot_ready_data is set to TRUE). This is helpful if you
  # want to do more customized pheatmap plotting.
  
  # Determine which of the unneeded cols are in the input tibble
  unneeded_cols <- unneeded_cols[unneeded_cols %in% colnames(in_tbl)]
  
  bact_genera <- extract_data(in_tbl, 'Bacteria', 'genus', unneeded_cols)
  bact_genera <- relative_abundance(bact_genera)  %>%
    get_top_by_rowMeans(30)
  bact_genera <- name_col_to_rownames(bact_genera)
  
  # Either return the pheatmap or the data for more customized plotting
  if (output_plot_ready_data == FALSE) {
    pheatmap(bact_genera)
  } else {
    bact_genera
  }
}


#--------------------------------------------------------------------------------------------------------#
# PLOTTING FUNGI ####
#--------------------------------------------------------------------------------------------------------#

get_fungi_genera <- function(in_tbl,
                             unneeded_cols = c("tax_id", "taxonomy", "type", "kingdom", "reference_length", "mean", "median", "max")) {
  
  unneeded_cols <- unneeded_cols[unneeded_cols %in% colnames(in_tbl)]
  
  in_tbl %>%
    # There are two candida genera. This one is the minor genera, will re-name to it's full name.
    mutate(name = ifelse(tax_id == 5475, "Saccharomycetales_incertae_sedis_Candida", name)) %>%
    extract_data('Fungi', 'genus', unneeded_cols)
}

plot_candida_vs_other_genera <- function(
  in_tbl,
  sample_order = FALSE,
  colors = c("#bdcebe", "#e06377"),
  hu_hashmap = FALSE,
  min_candida_abundance = 100
) {
  
  # in_tbl: Just a regular, old pathseq tibble. This will take care of
  # removing the unneccessary columns and so on!
  
  # If sample_order==FALSE, will sort in order of Candida abundance.
  # If hu_hashmap==FALSE (not entered), will not normalize to # human reads.
  
  # Convert to fungal genera dataframe
  fungal_genera <- get_fungi_genera(in_tbl)
  
  # Normalize to # human reads if specified
  if (!isFALSE(hu_hashmap)) {
    fungal_genera <- normalize_to_human_reads(fungal_genera, hu_hashmap)
  }
  
  # Make it taxon vs other
  fungal_genera <- make_taxon_vs_other(fungal_genera, "Candida")
  
  # Melt for plotting
  fungal_genera.m <- melt(fungal_genera, id="name")
  colnames(fungal_genera.m) <- c("Taxa", "Sample", "Abundance")
  
  # Need to define all levels
  levels = c("Other Taxa", "Candida")
  
  # Filter samples with less than specified candida abundance
  if (min_candida_abundance > 0) {
    samples_to_remove <- fungal_genera.m %>% 
      filter(Taxa == "Candida", Abundance < min_candida_abundance) %>%
      pull(Sample)
    
    fungal_genera.m <- fungal_genera.m %>%
      filter(!Sample %in% samples_to_remove)
  }
  
  # Find sample order if it wasn't input
  if (sample_order == FALSE) {
    sample_order <- fungal_genera.m %>%
      filter(Taxa=="Candida") %>%
      arrange(-Abundance) %>% pull(Sample) %>% as.character()
  }
  
  order <- match(fungal_genera.m$Sample, sample_order)
  fungal_genera.m$order <- order
  
  plt <- ggplot(fungal_genera.m,
         aes(x=reorder(Sample, order),
             y=Abundance,
             fill=factor(Taxa, levels=levels, ordered=T)
         )
  ) +
    geom_bar(stat='identity', color="white") + #color=white sets line between bars
    theme_classic(base_size = 10) +  #Sets theme, and makes all text 22pt
    theme(axis.title.x=element_blank(), # These element_blank()s remove x axis labels
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = colors) + 
    labs(fill = "Taxa") + #rename the legend title
    ylab("Abundance") + #rename y axis title
    scale_y_continuous(expand = c(0, 0), labels=comma) #Removes the white space bewteen graph and the bottom black line
  plt
}
