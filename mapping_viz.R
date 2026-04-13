library(ggplot2)
library(patchwork)

data_dir = '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/minimap'
depth_files = list.files(data_dir, 
                         pattern = "*.depth", 
                         full.names = TRUE)

## check percentage of the genome covered by at least X times
genome_covered = data.frame()
for (depth_file in depth_files) {
  if (grepl("barcode01",depth_file) | grepl("barcode02",depth_file)) {
    next
  }
  
  sample_name = sub(paste0(data_dir, "/"), "", depth_file)
  sample_name = sub(".depth", "", sample_name)
  print(sample_name)
  
  ## extract the sample and the reference genome
  barcode <- sub("_.*", "", sample_name)
  ref_strain  <- sub("^[^_]+_", "", sub("\\.depth$", "", sample_name))
  ref_strain = gsub(".mq20", "", ref_strain)
  print(barcode)
  print(ref_strain)
  
  ## read in the depth data
  depth_data = read.table(depth_file, header=FALSE)
  colnames(depth_data) = c('chrom', 'pos', 'depth')
  genome_size = dim(depth_data)[1]
  
  for (min_cov in 1:50) {
    print(min_cov)
    covered = sum(depth_data$depth >= min_cov)/genome_size
    print(covered)
    genome_covered = rbind.data.frame(genome_covered, c(sample_name, min_cov, covered, barcode, ref_strain))
  }
  
}

colnames(genome_covered) = c('sample_name', 'min_cov', 'frac_genome_covered', 'barcode', 'ref_strain')
genome_covered$filter = ifelse(grepl("mq20", genome_covered$sample_name), 'mq20', 'raw')
genome_covered$min_cov <- as.numeric(as.character(genome_covered$min_cov))
genome_covered$frac_genome_covered <- as.numeric(as.character(genome_covered$frac_genome_covered))


ggplot(data = genome_covered, aes(x=min_cov, y=frac_genome_covered)) + 
  geom_line(aes(color = sample_name, linetype = filter)) +
  theme_minimal() + 
  facet_wrap( ~ barcode)



## make boxplot of per chromosome
boxplots = list()
for (depth_file in depth_files) {
  sample_name = sub(paste0(data_dir, "/"), "", depth_file)
  sample_name = sub(".depth", "", sample_name)
  
  ## read in the depth data
  depth_data = read.table(depth_file, header=FALSE)
  colnames(depth_data) = c('chrom', 'pos', 'depth')
  
  ## remove small contigs
  if (sample_name == "barcode03" | sample_name == "barcode04") {
    depth_data = depth_data[grepl("^Tb927_", depth_data$chrom),]
  }  
  
  if (sample_name == "barcode05" | sample_name == "barcode06") {
    depth_data = depth_data[! grepl("bin_contig_", depth_data$chrom),]
    depth_data = depth_data[! grepl("BAC", depth_data$chrom),]
  }  
  
  
  
  unique(depth_data$chrom)  

  boxplot = ggplot(depth_data, aes(x=chrom, y=depth)) +
    geom_boxplot() +
    theme_minimal() + 
    ggtitle(sample_name)
  
  boxplots[[sample_name]] = boxplot

}

wrap_plots(boxplots, ncol=3) & 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) &
  coord_cartesian(ylim=c(0,400))




mmstats_file = '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/minimap/minimap_summary_stats.xlsx'
stats_data = read.xlsx(mmstats_file)
plot_data = melt(stats_data)
colnames(plot_data) = c('barcode', 'type', 'count')

ggplot(plot_data, aes(x = barcode, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +  # 'dodge' makes grouped bars
  labs(x = "Barcode", y = "Count", fill = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # rotate x labels


stats_data2 = stats_data[,c('barcode', 'reads', 'primary_mapped')]
stats_data2$perc_unique_mapped = stats_data2$primary_mapped / stats_data2$reads  
ggplot(stats_data2, aes(x = barcode, y = perc_unique_mapped)) +
  geom_bar(stat = "identity", position = "dodge") +  # 'dodge' makes grouped bars
  labs(x = "Barcode", y = "percentage unique mapped") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +   # rotate x labels 
  coord_cartesian(ylim = c(0,1))
