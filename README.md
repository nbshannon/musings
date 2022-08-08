# musings
MusiNGS (Microbiome using Next Generation Sequencing): computational analysis of nonhuman DNA sequences

MusiNGS: Microbiome using Next Generation Sequencing
Based on the hypothesis that a portion of the unmapped reads present in human sequencing data may represent sequences with microbial origin, a pipeline was constructed to align unmapped reads again a microbial sequence database.
After initial alignment of reads against human reference, the unmapped reads are extracted and fed into MusiNGS. For paired end data we use the first read in an unmapped pair (samtools view -f 76). We also exclude reads failing QC, or marked as duplicates (samtools view -F 1536)
Optional filters include additional human reference genomes (for those used see sup x), which are aligned using BWA and the exclusion of reads containing repeat sequences of 2 or more nucleotides (using repeatmasker).
Initial microbial alignment is performed against a database of bacterial, viral and fungal genomes obtained from the NCBI database using blastn (n=695285 sequences). At the same time a final human alignment is performed against RNA sequence data from NCBI also using blastn. Reads aligning to microbial sequence, but not human sequence are summarised.
To increase specificity contigs are assembled from reads xxx using velvet and realigned to microbial and human sequence in a repeat of the blastn step.


Fig x, Microbial sequence alignment pipeline (MusiNGS)
Supp x Human databases used in MusiNGS
NCBI Assembled chromosomes:
Each Downloaded from ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/
hs_ref_GRCh37.p10
hs_alt_CHM1_1.0 
hs_alt_HuRef
1000 genomes project reference genome:
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
Human genome and transcriptome BLAST database:
ftp://ftp.ncbi.nih.gov/blast/db/FASTA/human_genomic.gz
Ensembl Homo sapiens  cDNA database:
ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.70.cdna.all.fa.gz

Supp x2 Tools and parameters used in MusiNGS
Extraction of unmapped reads from aligned bams
Samtools (-f 76 -F 1536)
Require: Read unmapped, mate unmapped, first in pair
Exclude: read fails platform/vendor quality checks, read is PCR or optical duplicate

Alignment against human genomes
Bwa 

Repeat masking
Repeatmasker  (-species vertebrates –pa 4)

Blastn alignment
Blastn (-task megablast, -evalue 0.0000001, -word_size 16, -num_alignments 0, -num_descriptiontions 5 –dust no)

Contig assembly
velveth  (19 -fasta -short)
velvetg (-min_contig_lgth 75)
