library(ggplot2)


data_dir = '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/merging/'
length_files = list.files(path = data_dir, pattern = "lengths.txt", full.names = TRUE)

length_data = data.frame()
for (length_file in length_files) {
  sample = gsub("_lengths.txt", "", length_file)
  sample = gsub(".lengths.txt", "", sample)
  sample = gsub(data_dir, "", sample)
  sample = gsub("\\/", "", sample)
  print(sample)
  
  length_data_sample = read.csv(length_file)
  length_data_sample$sample = sample
  colnames(length_data_sample) = c('length', 'sample')
  length_data = rbind.data.frame(length_data, length_data_sample)
  
  print(sum(length_data_sample$length))
}

length_data_barcode = length_data[grep("barcode", length_data$sample),]
ggplot(length_data_barcode, aes(x=length)) + 
  geom_density(aes(color=sample)) + 
  theme_bw() + 
  coord_cartesian(xlim = c(0,100000))

length_data_merged = length_data[grep("merged", length_data$sample),]
ggplot(length_data_merged, aes(x=length)) + 
  geom_density(aes(color=sample)) + 
  theme_bw() + 
  coord_cartesian(xlim = c(0,100000))




