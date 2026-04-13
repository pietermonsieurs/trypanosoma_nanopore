library(ggplot2)
library(patchwork)

data_dir = '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/minimap_vsgs/'

depth_files = list.files(path = data_dir, 
                         pattern = ".depth",
                         full.names = TRUE)


## either select the files *with* or *without* filtering on the mapping quality

## WITHOUT mq20
# depth_files = depth_files[! grepl("mq20", depth_files)]

## WITH mq20
# depth_files = depth_files[grepl("mq20", depth_files)]

depth_files

plots = list()
for (depth_file in depth_files) {
  depth_data = read.table(depth_file)
  
  ## extract the sample name
  sample = gsub(data_dir, "", depth_file)
  sample = gsub("/", "", sample)
  sample = gsub(".depth", "", sample)
  print(sample)
  
  colnames(depth_data) = c('vsg', 'pos', 'depth')
  head(depth_data)

  covered_data = data.frame()
  
  ## per VSG stats
  for (vsg in unique(depth_data$vsg)) {
    depth_data_vsg = depth_data[depth_data$vsg == vsg, ]
    vsg_covered = sum(depth_data_vsg$depth > 0)
    vsg_length = nrow(depth_data_vsg)
    perc_cov_vsg = vsg_covered/vsg_length
    new_row = c(sample,vsg,perc_cov_vsg, vsg_length, vsg_covered)
    covered_data = rbind(covered_data, new_row)
  }
  colnames(covered_data) = c('sample', 'vsg', 'perc_cov', 'vsg_length', 'vsg_covered')
  
  covered_data$filter = ifelse(grepl("mq20", covered_data$sample), 'mq20', 'raw')
  covered_data$sample_name = gsub(".mq20", "", covered_data$sample)
  covered_data$perc_cov = as.numeric(covered_data$perc_cov)
  covered_data$vsg_length = as.numeric(covered_data$vsg_length)
  
  head(covered_data)
  p = ggplot(covered_data, aes(x = vsg, y = perc_cov, fill = filter)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Sample", y = "Percentage coverage", fill = "Filter") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(sample)
  
  plots[[sample]] = p
}

wrap_plots(plots[['barcode03']], plots[['barcode04']], ncol = 1)
wrap_plots(plots[['barcode03.mq20']], plots[['barcode04.mq20']], ncol = 1)

wrap_plots(plots[['barcode05.mq20']] + theme(axis.title.x = element_blank()), 
           plots[['barcode06.mq20']] + theme(axis.title.x = element_blank()),
           ncol = 1)



