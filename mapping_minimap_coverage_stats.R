library(ggplot2)

data_dir = '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/minimap/'

depth_files = list.files(path = data_dir, 
                         pattern = ".depth",
                         full.names = TRUE)


## either select the files *with* or *without* filtering on the mapping quality

## WITHOUT mq20
# depth_files = depth_files[! grepl("mq20", depth_files)]

## WITH mq20
# depth_files = depth_files[grepl("mq20", depth_files)]

depth_files





covered_data = data.frame()
covered_data_perchrom = data.frame()


for (depth_file in depth_files) {
  
  if (grepl("barcode01",depth_file) | grepl("barcode02",depth_file)) {
    next
  }
  
  depth_data = read.table(depth_file)
  colnames(depth_data) = c('chrom', 'pos', 'depth')
  head(depth_data)
  
  ## extract the sample name
  sample = gsub(data_dir, "", depth_file)
  sample = gsub("/", "", sample)
  sample = gsub(".depth", "", sample)
  print(sample)
  

  ## extract the sample and the reference genome
  barcode <- sub("_.*", "", sample)
  ref_strain  <- sub("^[^_]+_", "", sub("\\.depth$", "", sample))
  ref_strain = gsub(".mq20", "", ref_strain)
  print(barcode)
  print(ref_strain)

  ## step 1: overall coverage stats
  genome_covered = sum(depth_data$depth > 0)
  genome_length = nrow(depth_data)
  perc_cov_genome = genome_covered/genome_length
  new_row = c(sample, barcode, ref_strain, perc_cov_genome, genome_covered, genome_length)
  covered_data = rbind(covered_data, new_row)
    
  ## step 2: per chromosome stats
  # 
  # for (chrom in unique(depth_data$chrom)) {
  #   depth_data_chrom = depth_data[depth_data$chrom == chrom, ]
  #   chrom_covered = sum(depth_data_chrom$depth > 0)
  #   chrom_length = nrow(depth_data_chrom)
  #   perc_cov_chrom = chrom_covered/chrom_length
  #   new_row = c(sample,chrom,perc_cov_chrom, chrom_length, chrom_covered)
  #   covered_data_perchrom = rbind(covered_data_perchrom, new_row)
  # }
  # 
}

colnames(covered_data) = c('sample', 'barcode', 'ref_strain', 'perc_cov', 'covered_bases', 'length')
# colnames(covered_data_perchrom) = c('sample', 'chrom', 'perc_cov', 'covered_bases', 'length')

covered_data$filter = ifelse(grepl("mq20", covered_data$sample), 'mq20', 'raw')
# covered_data_perchrom$filter = ifelse(grepl("mq20", covered_data_perchrom$sample), 'mq20', 'raw')

covered_data$sample_name = gsub(".mq20", "", covered_data$sample)
# covered_data_perchrom$sample_name = gsub(".mq20", "", covered_data_perchrom$sample)


covered_data$perc_cov = as.numeric(covered_data$perc_cov)
# covered_data_perchrom$perc_cov = as.numeric(covered_data_perchrom$perc_cov)
covered_data$length = as.numeric(covered_data$length)
# covered_data_perchrom$perc_cov = as.numeric(covered_data_perchrom$length)

## overall coverage
ggplot(covered_data, aes(x = sample_name, y = perc_cov, fill = filter)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Sample", y = "Percentage coverage", fill = "Filter") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # + 
  # facet_wrap( ~ barcode, ncol = 1)





## coverage per chromosome

for (barcode_nr in seq(3,6,1)) {
  barcode = paste0('barcode0', barcode_nr)
  covered_data_perchrom_filter = covered_data_perchrom[covered_data_perchrom$length > 500000,]
  covered_data_perchrom_filter = covered_data_perchrom_filter[covered_data_perchrom_filter$sample_name == barcode,]
  
  ggplot(covered_data_perchrom_filter, aes(x = chrom, y = perc_cov, fill = filter)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Sample", y = "Percentage coverage", fill = "Filter") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(barcode)
}
