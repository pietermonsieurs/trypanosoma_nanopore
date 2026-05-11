library(circlize)
library(ggplot2)

paf_dir = "/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_QC/"
cytoscape_dir = paste0(paf_dir, 'cytoscape/')
paf_files = list.files(paf_dir, pattern = "filtered.paf", full.names = TRUE)

## read in the text file with the repetitive regions
rep_regions_file = '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_QC/repeatfinder/consensus_TbgI_merged.fasta.2.5.5.80.15.40.100.parsed.dat'
rep_regions = read.table(rep_regions_file, header = TRUE, stringsAsFactors = FALSE)
head(rep_regions)

paf_all = data.frame()
for (paf_file in paf_files) {
  
  #### 1. preparation of data ####
  
  ## read in the PAF file
  # paf_file = '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_QC/consensus_bc03_vs_TREU927_filtered.paf'
  paf <- read.table(paf_file, header = FALSE, stringsAsFactors = FALSE)
  sample_name = basename(paf_file)
  sample_name = gsub(".paf", "", sample_name)
  
  ## give column names to the paf file
  colnames(paf)[1:9] <- c(
    "qname", "qlen", "qstart", "qend",
    "strand",
    "tname", "tlen", "tstart", "tend"
  )
  head(paf)
  
  ## do some replacing of chromosome names to make it more readable in the 
  ## circos plots
  paf$qname_original = paf$qname
  paf$tname = gsub("Tb927_", "", paf$tname)
  paf$tname = gsub("Tbg972_", "", paf$tname)
  paf$tname = gsub("_v5.1", "", paf$tname)
  paf$qname = gsub("contig_", "C", paf$qname)
  paf$qname = gsub("scaffold_", "S", paf$qname)
  
  ## do some additional filtering, only those alignments which are bigger 
  ## than 100,000 nt for example
  paf = paf[paf[,10] > 50000,]
  
  ## create an additional column to highlight the contig and from which 
  ## assembly it comes from. This can it a next loop be used to loop 
  ## over the different chromosomes. Also, extract the name of the reference
  ## as we have to run per reference genome too (DAL972 and TREU927)
  sample_name_short = strsplit(sample_name, "_")[[1]][2]
  paf$qname_full = paste0(sample_name_short, "_", paf$qname)
  
  ## extract ref genome
  if (grepl("merged", sample_name)) {
    ref_name = strsplit(sample_name, "_")[[1]][5]
  }else{
    ref_name = strsplit(sample_name, "_")[[1]][4]
  }
  paf$ref_genome = ref_name
  paf$strain = sample_name_short

  ## add to the overall data frame
  paf_all = rbind.data.frame(paf_all, paf)
  
  
  #### 2. define sectors ####
  
  # Unique contigs and chromosomes
  chroms  <- sort(unique(paf$tname))
  contigs <- rev(unique(paf[order(paf$tname), ]$qname))
  

  # Named vector of lengths
  contig_lengths <- tapply(paf$qlen, paf$qname, max)
  contig_lengths_sorted = contig_lengths[contigs]
  chrom_lengths  <- tapply(paf$tlen, paf$tname, max)
  
  all_sectors <- c(contigs, chroms)
  all_lengths <- c(contig_lengths_sorted, chrom_lengths)
  
  
  #### 3. enforce left contigs - right reference genome ####
  n_contigs <- length(contigs)
  n_chroms  <- length(chroms)
  
  ## add additional gaps between the contigs and the chromosomes
  gap.after <- c(
    rep(1, n_contigs - 1),
    10,  # big gap between contigs and chromosomes
    rep(1, n_chroms - 1),
    10   # close the circle
  )
  
  
  
  #### 4. creating the circos plot ####
  ## start creating the circos picture
  pdf_file = paste0(paf_dir, sample_name, '.pdf')
  pdf(pdf_file, width=8, height=8)
  
  circos.clear()
  circos.par(cell.padding = c(0.001, 0, 0.001, 0))
  
  circos.par(gap.after = gap.after)
  
  circos.initialize(
    factors = all_sectors,
    xlim = cbind(rep(0, length(all_lengths)), all_lengths)
  )
  
  
  
  circos.trackPlotRegion(
    ylim = c(0,1),
    track.height = 0.20,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      circos.text(
        x = mean(get.cell.meta.data("xlim")),
        y = 0.5,
        labels = sector.name,
        cex = .4,
        facing = "inside"
      )
    }
  )
  
  for (i in seq_len(nrow(paf))) {
    circos.link(
      paf$qname[i], c(paf$qstart[i], paf$qend[i]),
      paf$tname[i], c(paf$tstart[i], paf$tend[i]),
      col = ifelse(paf$strand[i] == "+", "#1f77b4AA", "#d62728AA")
    )
  }
  dev.off()
  
}  


## loop over the different chromosomes, and make one plot per chromosome but
## combined for all 4 or 5 strains

for (chrom in unique(paf_all$tname)) {
  for (ref in unique(paf_all$ref_genome)) {
    paf = paf_all[paf_all$tname == chrom & paf_all$ref_genome == ref, c(colnames(paf_all)[1:10], 'qname_full')]
    
    # Unique contigs and chromosomes
    chroms  <- sort(unique(paf$tname))
    contigs <- rev(unique(paf[order(paf$tname), ]$qname_full))
    
    chroms  <- unique(paf$tname)
    contigs <- unique(paf$qname_full)
    
    # Named vector of lengths
    contig_lengths <- tapply(paf$qlen, paf$qname_full, max)
    contig_lengths_sorted = contig_lengths[contigs]

    chrom_lengths  <- tapply(paf$tlen, paf$tname, max)
    
    all_sectors <- c(contigs, chroms)
    all_lengths <- c(contig_lengths, chrom_lengths)
    
    
    #### 3. enforce left contigs - right reference genome ####
    n_contigs <- length(contigs)
    n_chroms  <- length(chroms)
    
    ## add additional gaps between the contigs and the chromosomes
    gap.after <- c(
      rep(1, n_contigs - 1),
      10,  # big gap between contigs and chromosomes
      rep(1, n_chroms - 1),
      10   # close the circle
    )
    
    
    
    #### 4. creating the circos plot ####
    ## start creating the circos picture
    pdf_file = paste0(paf_dir, chrom, '_', ref, '.pdf')
    pdf(pdf_file, width=8, height=8)

    circos.clear()
    circos.par(cell.padding = c(0.001, 0, 0.001, 0))
    
    circos.par(gap.after = gap.after)
    
    circos.initialize(
      factors = all_sectors,
      xlim = cbind(rep(0, length(all_lengths)), all_lengths)
    )
    
    
    
    circos.trackPlotRegion(
      ylim = c(0,1),
      track.height = 0.20,
      panel.fun = function(x, y) {
        sector.name <- get.cell.meta.data("sector.index")
        circos.text(
          x = mean(get.cell.meta.data("xlim")),
          y = 0.5,
          labels = sector.name,
          cex = .7,
          facing = "inside"
        )
      }
    )
    
    for (i in seq_len(nrow(paf))) {
      print(i)
      circos.link(
        paf$qname_full[i], c(paf$qstart[i], paf$qend[i]),
        paf$tname[i], c(paf$tstart[i], paf$tend[i]),
        col = ifelse(paf$strand[i] == "+", "#1f77b440", "#d6272840")
      )
    }
    dev.off()
  }  
}


## loop over the different chromosomes and loop over the different strains, and 
## creqte one plot per chromosome and per reference genome, for each strain
## separately

for (chrom in unique(paf_all$tname)) {
  for (ref in unique(paf_all$ref_genome)) {
    for (strain in unique(paf_all$strain)) {
      paf = paf_all[paf_all$tname == chrom & 
                      paf_all$ref_genome == ref & 
                      paf_all$strain == strain, 
                    c(colnames(paf_all)[1:10],'qname_full')]
      
      # Unique contigs and chromosomes
      chroms  <- sort(unique(paf$tname))
      contigs <- rev(unique(paf[order(paf$tname), ]$qname_full))
      
      chroms  <- unique(paf$tname)
      contigs <- unique(paf$qname_full)
      
      # Named vector of lengths
      contig_lengths <- tapply(paf$qlen, paf$qname_full, max)
      contig_lengths_sorted = contig_lengths[contigs]
      
      chrom_lengths  <- tapply(paf$tlen, paf$tname, max)
      
      all_sectors <- c(contigs, chroms)
      all_lengths <- c(contig_lengths, chrom_lengths)
      
      
      #### 3. enforce left contigs - right reference genome ####
      n_contigs <- length(contigs)
      n_chroms  <- length(chroms)
      
      ## add additional gaps between the contigs and the chromosomes
      gap.after <- c(
        rep(1, n_contigs - 1),
        10,  # big gap between contigs and chromosomes
        rep(1, n_chroms - 1),
        10   # close the circle
      )
      
      
      
      #### 4. creating the circos plot ####
      ## start creating the circos picture
      pdf_file = paste0(paf_dir, 'per_strain/', strain, '_', chrom, '_', ref, '.pdf')
      pdf(pdf_file, width=8, height=8)
      
      circos.clear()
      circos.par(cell.padding = c(0.001, 0, 0.001, 0))
      
      circos.par(gap.after = gap.after)
      
      circos.initialize(
        factors = all_sectors,
        xlim = cbind(rep(0, length(all_lengths)), all_lengths)
      )
      
      
      
      # circos.trackPlotRegion(
      #   ylim = c(0,1),
      #   track.height = 0.20,
      #   panel.fun = function(x, y) {
      #     sector.name <- get.cell.meta.data("sector.index")
      #     circos.text(
      #       x = mean(get.cell.meta.data("xlim")),
      #       y = 0.5,
      #       labels = sector.name,
      #       cex = .7,
      #       facing = "inside"
      #     )
      #   }
      # )
      
      circos.trackPlotRegion(
        ylim = c(0,1),
        track.height = 0.20,
        panel.fun = function(x, y) {
          sector.name <- get.cell.meta.data("sector.index")
          
          # draw repeat regions
          df <- rep_regions[rep_regions$sequence_id == sector.name, ]
          
          if(nrow(df) > 0) {
            circos.rect(
              xleft  = df$start,
              ybottom = 0,
              xright = df$end,
              ytop   = 0.2,   # small band so it doesn't cover text
              col = rgb(1, 0, 0, 0.5),
              border = NA
            )
          }
          
          # existing labels
          circos.text(
            x = mean(get.cell.meta.data("xlim")),
            y = 0.5,
            labels = sector.name,
            cex = .7,
            facing = "inside"
          )
        }
      )
      
      
      for (i in seq_len(nrow(paf))) {
        print(i)
        circos.link(
          paf$qname_full[i], c(paf$qstart[i], paf$qend[i]),
          paf$tname[i], c(paf$tstart[i], paf$tend[i]),
          col = ifelse(paf$strand[i] == "+", "#1f77b440", "#d6272840")
        )
      }
      dev.off()
    }
  }  
}


## for TbgI we also have repetitive regions which we can try to map on the 
## outer circle. Same code as above, but now adding additional part to 
## visualize the repetitve regions
paf_all_tbgI = paf_all[paf_all$strain == "TbgI",]
rep_regions$sequence_id_to_match = paf_all_tbgI[match(rep_regions$sequence_id, paf_all_tbgI$qname_original),]$qname_full

strains_with_rep_regions = c('TbgI')
for (chrom in unique(paf_all$tname)) {
  for (ref in unique(paf_all$ref_genome)) {
    for (strain in strains_with_rep_regions) {
      paf = paf_all[paf_all$tname == chrom & 
                      paf_all$ref_genome == ref & 
                      paf_all$strain == strain, 
                    c(colnames(paf_all)[1:10],'qname_full')]
      
      # Unique contigs and chromosomes
      chroms  <- sort(unique(paf$tname))
      contigs <- rev(unique(paf[order(paf$tname), ]$qname_full))
      
      chroms  <- unique(paf$tname)
      contigs <- unique(paf$qname_full)
      
      # Named vector of lengths
      contig_lengths <- tapply(paf$qlen, paf$qname_full, max)
      contig_lengths_sorted = contig_lengths[contigs]
      
      chrom_lengths  <- tapply(paf$tlen, paf$tname, max)
      
      all_sectors <- c(contigs, chroms)
      all_lengths <- c(contig_lengths, chrom_lengths)
      
      
      #### 3. enforce left contigs - right reference genome ####
      n_contigs <- length(contigs)
      n_chroms  <- length(chroms)
      
      ## add additional gaps between the contigs and the chromosomes
      gap.after <- c(
        rep(1, n_contigs - 1),
        10,  # big gap between contigs and chromosomes
        rep(1, n_chroms - 1),
        10   # close the circle
      )
      
      
      
      #### 4. creating the circos plot ####
      ## start creating the circos picture
      pdf_file = paste0(paf_dir, 'per_strain/', strain, '_', chrom, '_', ref, '_with_rep_regions.pdf')
      pdf(pdf_file, width=8, height=8)
       
      circos.clear()
      circos.par(cell.padding = c(0.001, 0, 0.001, 0))
      
      circos.par(gap.after = gap.after)
      
      circos.initialize(
        factors = all_sectors,
        xlim = cbind(rep(0, length(all_lengths)), all_lengths)
      )
      
      
      circos.trackPlotRegion(
        ylim = c(0,1),
        track.height = 0.20,
        panel.fun = function(x, y) {
          sector.name <- get.cell.meta.data("sector.index")
          
          # draw repeat regions
          df <- rep_regions[rep_regions$sequence_id_to_match == sector.name, ]
          
          if(nrow(df) > 0) {
            circos.rect(
              xleft  = df$start,
              ybottom = 0,
              xright = df$end,
              ytop   = 1,   # small band so it doesn't cover text
              col = "darkgray",
              border = NA
            )
          }
          
          # existing labels
          circos.text(
            x = mean(get.cell.meta.data("xlim")),
            y = 0.5,
            labels = sector.name,
            cex = .7,
            facing = "inside"
          )
        }
      )
      
      
      for (i in seq_len(nrow(paf))) {
        print(i)
        circos.link(
          paf$qname_full[i], c(paf$qstart[i], paf$qend[i]),
          paf$tname[i], c(paf$tstart[i], paf$tend[i]),
          col = ifelse(paf$strand[i] == "+", "#1f77b440", "#d6272840")
        )
      }
      dev.off()
    }
  }  
}



## export to cytoscape. For this we need to have edges and nodes files that 
## can be loaded into Cytoscape



paf_all = data.frame()
for (paf_file in paf_files) {
  
  #### 1. preparation of data ####
  
  ## read in the PAF file
  # paf_file = '/Users/pmonsieurs/programming/trypanosoma_nanopore/results/assembly_QC/consensus_bc03_vs_TREU927_filtered.paf'
  paf <- read.table(paf_file, header = FALSE, stringsAsFactors = FALSE)
  sample_name = basename(paf_file)
  sample_name = gsub(".paf", "", sample_name)
  
  ## give column names to the paf file
  colnames(paf)[1:9] <- c(
    "qname", "qlen", "qstart", "qend",
    "strand",
    "tname", "tlen", "tstart", "tend"
  )
  head(paf)
  
  ## do some replacing of chromosome names to make it more readable in the 
  ## circos plots
  paf$tname = gsub("Tb927_", "", paf$tname)
  paf$tname = gsub("Tbg972_", "", paf$tname)
  paf$tname = gsub("_v5.1", "", paf$tname)
  paf$qname = gsub("contig_", "C", paf$qname)
  paf$qname = gsub("scaffold_", "S", paf$qname)
  
  ## do some additional filtering, only those alignments which are bigger 
  ## than 100,000 nt for example
  paf = paf[paf[,10] > 50000,]
  
  ## create an additional column to highlight the contig and from which 
  ## assembly it comes from. This can it a next loop be used to loop 
  ## over the different chromosomes. Also, extract the name of the reference
  ## as we have to run per reference genome too (DAL972 and TREU927)
  sample_name_short = strsplit(sample_name, "_")[[1]][2]
  paf$qname_full = paste0(sample_name_short, "_", paf$qname)
  
  ## extract ref genome
  if (grepl("merged", sample_name)) {
    ref_name = strsplit(sample_name, "_")[[1]][5]
  }else{
    ref_name = strsplit(sample_name, "_")[[1]][4]
  }
  paf$ref_genome = ref_name
  
  ## add to the overall data frame
  paf_all = rbind.data.frame(paf_all, paf)
  
  
  ## create the edages file

  edges <- data.frame()
  singletons <- data.frame()
  
  for (chr in unique(paf$tname)) {
    
    sub <- paf[paf$tname == chr, ]
    sub <- sub[order(sub$tstart), ]
    
    if (nrow(sub) > 1) {
      df <- data.frame(
        source = sub$qname[-nrow(sub)],
        target = sub$qname[-1],
        chromosome = chr
      )
      
      edges <- rbind(edges, df)
      
    } else {
      # store BOTH contig and chromosome
      singletons <- rbind(singletons, data.frame(
        contig = sub$qname,
        chromosome = chr
      ))
    }
  }
  
  ## for the singletons, create a "self-edge" to make clear that they are 
  ## the only contig mapping to a chromosome, so (potentially) perfect assembly
  singleton_edges <- data.frame(
    source = singletons$contig,
    target = singletons$contig,
    chromosome = singletons$chromosome
  )
  
  edges <- rbind(edges, singleton_edges)
  
  
  ## create the nodes files
  nodes <- data.frame(
    id = unique(c(edges$source, edges$target))
  )
  nodes$length <- tapply(paf$qlen, paf$qname, max)[nodes$id]
  nodes
  
  ## export to text files to visualize in cytoscape
  edge_file = paste0(cytoscape_dir, sample_name, '_edges.txt')
  node_file = paste0(cytoscape_dir, sample_name, '_nodes.txt')
  
  write.table(edges[, c("source", "target", "chromosome")],
              edge_file,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  
  write.table(nodes,
              node_file,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  
}
