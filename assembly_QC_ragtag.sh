## installation of ragtag: 
# conda create -n ragtag-env python=3.11
# conda activate ragtag-env
# conda install -c bioconda ragtag

## run ragtag to scaffold the assembly
ragtag.py scaffold -t 4 \
    /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-67_TbruceiTREU927_Genome_nuclear.fasta \
    /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/consensus_TbgI_merged.fasta \
    -o /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/ragtag/ragtag_TREU927/ 
    #-m 5000

awk '$5 == "W"' /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/ragtag/ragtag_TREU927/ragtag.scaffold.agp


## do the same but do now for another reference genome DAL972
ragtag.py scaffold -t 4 \
    /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-68_TbruceigambienseDAL972_Genome.fasta \
    /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/consensus_TbgI_merged.fasta \
    -o /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/ragtag/ragtag_DAL972/ 
    #-m 5000

awk '$5 == "W"' /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/ragtag/ragtag_DAL972/ragtag.scaffold.agp

