library(readr)

load_gtf <- function(filepath) {
  # Load a gtf or gff file, at the given filepath, into a dataframe
  gtf_cols <- c("chromosome", "source", "feature", "start", "end", "score",
                "strand", "frame", "attributes")
  genes <- data.frame(read_tsv(filepath, col_names=gtf_cols, comment='#'))
  genes$frame[genes$frame == '.'] <- '0'
  return(genes)
}

merge_gtfs <- function(gtf_list) {
  # Given some gts files as dataframes in a list, rbind and drop duplicates.
  # Among duplicates (ignoring the attributes column), the first entry is kept.
  all_genes <- Reduce(rbind, gtf_list)
  essentials <- cbind(all_genes[1], all_genes[4], all_genes[5], all_genes[7])
  return(all_genes[!duplicated(essentials), ])
}

clean_gerstein_2014 <- function(filepath) {
  # Given a filepath to the Gersetin et al 2014 data, produce a merged GTF file
  files <- c(paste(filepath, "comparable/lncRNA_exons.gtf", sep="/"),
             paste(filepath, "comparable/miRNA.gtf", sep="/"),
             paste(filepath, "comparable/pri_worm_intergenic_WS220.gff",
                   sep="/"),
             paste(filepath, "comparable/snoRNA.gtf", sep="/"),
             paste(filepath, "comparable/snRNA.gtf", sep="/"),
             paste(filepath, "comparable/tRNA.gtf", sep="/"),
             paste(filepath, "noncomparable/worm_non_consensus_ncRNAs.gtf",
                   sep="/"))
  genes <- lapply(files, load_gtf)
  unique_genes <- merge_gtfs(genes)
  novelty <- nrow(unique_genes) / sum(sapply(genes, nrow))
  write.table(unique_genes, file=paste(filepath, "gerstein_2014.gtf", sep='/'),
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, na=".")
  return(novelty)
}

# clean_gerstein_2014("ncRNA_Literature/Gerstein_2014")


clean_lu_2011 <- function(filepath) {
  # Given a filepath to the Lu et al 2011 data, convert (header-removed)
  # arff files to gtf files
  files <- c(paste(filepath, "full_data", sep="/"),
             paste(filepath, "gold", sep="/"))
  for (file in files) {
    genes <- read.csv(paste(file, "arff", sep="."))
    gtf_genes <- data.frame(chromosome = genes$chromosome,
                            source = rep("Lu2011", nrow(genes)),
                            feature = genes$class,
                            start = genes$start,
                            end = genes$end,
                            score = rep(".", nrow(genes)),
                            strand = genes$strand,
                            frame = integer(nrow(genes)),
                            attributes = rep("", nrow(genes)))
    write.table(gtf_genes, file=paste(file, "gtf", sep="."), quote=FALSE,
                sep="\t", col.names=FALSE, row.names=FALSE, na=".")
  }
}

# clean_lu_2011("ncRNA_Literature/Lu_2011")


strand_sign <- function(strand) {
  # Convert an integer representation of strand (-1 / 1) into a character
  # representation (+ / -)
  if (sign(strand) == 1) {
    return("+")
  } else {
    return("-")
  }
}


clean_wormbase <- function(filepath) {
  # Given a file to the Wormbase dataset file, comvert it from tsv to gtf.
  genes <- data.frame(read_tsv(paste(filepath, "C_elegans_ncRNA.tsv", sep="/"),
                               col_names=TRUE))
  gtf_genes <- data.frame("chromosome" = genes$Gene.chromosome.primaryIdentifier,
                          "source" = rep("Wormbase", nrow(genes)),
                          "feature" = genes$Gene.transcripts.method,
                          "start" = genes$Gene.chromosomeLocation.start,
                          "end" = genes$Gene.chromosomeLocation.end,
                          "score" = genes$Gene.score,
                          "strand" = sapply(genes$Gene.chromosomeLocation.strand,
                                            strand_sign),
                          "frame" = integer(nrow(genes)),
                          "attributes" = paste("gene_symbol \"",
                                               genes$Gene.symbol,
                                               "\"; ",
                                               "transcript_symbol \"",
                                               genes$Gene.transcripts.symbol,
                                               "\"",
                                               sep=""))
  write.table(gtf_genes, file=paste(filepath, "wormbase.gtf", sep="/"),
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, na=".")
}

# clean_wormbase("ncRNA_Literature/Wormbase")


clean_WBcel235 <- function(filepath) {
  # Given a filepath to the NCBI data folder, select genes labelled as "ncRNA"
  # or "Non-coding" from the wbcel235 dataset
  genes <- load_gtf(paste(filepath, "WBcel235.104.gtf", sep="/"))
  print(nrow(genes))
  nc_genes <- genes[(grepl("gene_biotype \"tRNA\"", genes$attributes) |
                    grepl("gene_biotype \"ncRNA\"", genes$attributes) |
                    grepl("gene_biotype \"rRNA\"", genes$attributes) |
                    grepl("gene_biotype \"snRNA\"", genes$attributes) |
                    grepl("gene_biotype \"snoRNA\"", genes$attributes) |
                    grepl("gene_biotype \"lincRNA\"", genes$attributes) |
                    grepl("gene_biotype \"miRNA\"", genes$attributes)), ]
  write.table(nc_genes, file=paste(filepath, "WBcel235.104_nc.gtf", sep="/"),
              quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, na=".")
  return(nrow(nc_genes) / nrow(genes))
}

clean_WBcel235("ncRNA_Literature/wbcel235.104")


make_combined_dataset <- function(in_filepaths, out_filepath) {
  # Given a vector of filepaths to GTF files (in_filepaths), combine them into a
  # single dataset, remove duplicates, and save at out_filepath. Returns the
  # novelty of the combined dataset (fraction of total genes which are unique)
  genes <- lapply(in_filepaths, load_gtf)
  unique_genes <- merge_gtfs(genes)
  novelty <- nrow(unique_genes) / sum(sapply(genes, nrow))
  write.table(unique_genes, file=out_filepath, quote=FALSE, sep="\t",
              col.names=FALSE, row.names=FALSE, na='.')
  return(novelty)
}

# make_combined_dataset(c('ncRNA_Literature/wbcel235.104/wbcel235.104_nc.gtf',
#                         'ncRNA_Literature/Gerstein_2014/gerstein_2014_wbc235.gtf'),
#                       'ncRNA_Literature/Combined/wbc235.104nc_Gerstein.gtf')
# make_combined_dataset(c('ncRNA_Literature/wbcel235.104/wbcel235.104_nc.gtf',
#                         'ncRNA_Literature/Wormbase/wormbase.gtf'),
#                       'ncRNA_Literature/Combined/wbc235.104nc_Wormbase.gtf')


make_test_dataset <- function(nsample, in_filepath, out_filepath) {
  # Given in_filepath and a path to a GTF file, create a dataset containing
  # nsample randomly selected genes from the full dataset. Intended for testing
  # other algorithms before running on the full set.
  all_genes <- load_gtf(in_filepath)
  sampled_genes <- all_genes[floor(runif(nsample, min=1, max=nrow(all_genes))), ]
  write.table(sampled_genes, file=out_filepath, quote=FALSE, sep="\t",
              col.names=FALSE, row.names=FALSE, na='.')
}

# make_test_dataset(1000,
#                   'ncRNA_Literature/wbcel235.104/wbcel235.104_nc.gtf',
#                   'ncRNA_Literature/test/wbcel235.104_nc_1000.gtf')
# make_test_dataset(1000,
#                   'ncRNA_Literature/Combined/wbc235.104nc_Gerstein.gtf',
#                   'ncRNA_Literature/test/wbcel235.104nc_Gerstein_1000.gtf')
# make_test_dataset(1000,
#                   'ncRNA_Literature/Combined/wbc235.104nc_Wormbase.gtf',
#                   'ncRNA_Literature/test/wbcel235.104nc_Wormbase_1000.gtf')


make_chr_subset_dataset <- function(nsample, in_filepath, out_filepath) {
  # Given in_filepath and a path to a GTF file, create a dataset containing
  # nsample randomly selected genes from the full dataset which are on
  # Chromosomes I or II. Intended for testing
  # other algorithms before running on the full set.
  all_genes <- load_gtf(in_filepath)
  chr_genes <- all_genes[all_genes$chromosome %in% c('I', 'II'), ]
  sampled_genes <- chr_genes[floor(runif(nsample, min=1, max=nrow(chr_genes))), ]
  write.table(sampled_genes, file=out_filepath, quote=FALSE, sep="\t",
              col.names=FALSE, row.names=FALSE, na='.')
}

# make_chr_subset_dataset(1000,
#                         'ncRNA_Literature/wbcel235.104/wbcel235.104_nc.gtf',
#                         'ncRNA_Literature/test/wbcel235.104_nc_chr_500.gtf')
