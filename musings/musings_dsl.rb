# MusiNGS DSL
# This is the DSL wrapper for musings functionality
# Providing alignment and filtering of microbial sequence from bam files
# Using a compute cluster (LSF system)

# Author:: Nicholas Shannon (mailto:nbshannon@gmail.com)
# ---

# require the core musings code
require '<path>/musings.rb'

# Initialise muse (initialised with a default clusterSize of 150
MUSE = Musings::Muse.new

# = DSL functionality

# == Setup of run

# Set the output directory
def set_output_dir(outputDir)
	MUSE.set_output_dir(outputDir)
end

# construct directory structure in output directory
def make_directories
	MUSE.make_directories
end

# Specify input unmapped bam
def set_unmapped_bam(unmappedBam)
	MUSE.set_unmapped_bam(unmappedBam)
end

# Extract unmapped reads from another bam (as someMappedBam.unmapped.bam, ssets unmappedBam)
def	extract_unmapped_bam(someMappedBam)
	MUSE.extract_unmapped_bam(someMappedBam)
end

# Set number of fragments to split reads into for batch alignment and repeat masking (150 used as default)
def set_batch_size(clusterSize)
	MUSE.set_batch_size(clusterSize)
end

# Split reads for batch alignment and repeat masking
# Splits reads into a number of fasta files, as indicated by clusterSize
def split_unmapped_bam
	MUSE.split_unmapped_bam
end

# Add additional references for bwa alignment step (alignments will run on each of these)
def add_bwa(bwaDb)
	MUSE.add_bwa(bwaDb)
end

# Set microbial blast database used for microbial sequence identification
def set_microbial_blast_db(microbialDb)
	MUSE.set_microbial_blast_db(microbialDb)
end

# Set human blast database
def set_human_blast_db(humanBlastDb)
	MUSE.set_human_blast_db(humanBlastDb)
end

# == Start alignment

# Start bwa alignment jobs
def start_bwa_jobs
	MUSE.start_bwa_jobs
end

# Start repeat masking
def start_repeat_masking
	MUSE.start_repeat_masking
end

# Start microbial and human blast jobs
def start_blast_jobs
	MUSE.start_blast_jobs
end

# == Import results

# Read human alignmed reads from bwa alignments
def import_bwa_results
	MUSE.import_bwa_results
end

# Read human aligned reads from blast
def import_human_blast_results
	MUSE.import_human_blast_results
end

# Read microbial aligned reads from blast (without applying filtering)
def import_microbial_blast_results
	MUSE.import_microbial_blast_results
end

# Read repeat masked reads from RepeatMasker
def import_repeat_masking
	MUSE.import_repeat_masking
end

# == Filter and export microbial reads

# Process microbial blast results to remove human and masked reads
def strip_blast_results
	MUSE.strip_blast_results
end

# Extract microbial sequences from split fasta files
def extract_microbial_reads
	MUSE.extract_microbial_reads
end

# == Contig assembly and alignment

# Assemble contigs form filtered microbial reads and align against human and microbial blast databases
def assemble_and_align_contigs
	MUSE.assemble_and_align_contigs
end

# Build contigs from microbial sequences # deprecated in favour of assembling and aligning contigs in one step
#def build_contigs
#	MUSE.build_contigs
#end

# Start microbial and human blast jobs for contig alignment 
#def start_contig_blast_jobs # deprecated in favour of assembling and aligning contigs in one step
#	MUSE.start_contig_blast_jobs
#end

# Process contig blast results to remove human reads
def strip_contig_blast_results
	MUSE.strip_contig_blast_results
end

# Delete temporary files from contig directory
def cleanup_contig_directory
	MUSE.cleanup_contig_directory
end

# == Additional functionality

# Seeds list of queries from unmapped bam (non-destructive to current query list)
def seed_queries
	MUSE.seed_queries
end

# Write query information (human query and masking status)
def write_queries
	MUSE.write_queries
end

# Check for existence of queries file
def does_queries_exist
	MUSE.does_queries_exist
end

# Read query information (human query and masking status)
def read_queries
	MUSE.read_queries
end

# Summarise query counts (counts for each combination of query filtering status)
def summarise_query_counts
	MUSE.summarise_query_counts
end

# Export muse for direct import in subsequent usage (query writing recommended instead)
def export_muse
	MUSE.export_muse
end

# Import prior muse run
def load_muse(museFilename=nil)
	MUSE.load_muse(museFilename)
end

# summaries mapping, not current implemented
def summarise_mapping
	MUSE.summarise_mapping
end

# check successful completion of jobs, not currently implemented
def check_job_success
	MUSE.check_job_success
end

# zip results in output folder, not currently implemented
def zip_results
	MUSE.zip_results
end

# split and print query summary file
def show_query_summary
	MUSE.show_query_summary
end

# deletes bwa, repeat and split directories
def cleanup
	MUSE.cleanup
end

# Reads in queries for cross tabbing (per database mapped status # not implemented in preference for rewriting @queries hash
#def crosstab_queries
#	MUSE.crosstab_queries
#end

# Reruns blast alignment for specified cluster number (by default 1 to 150)
def rerun_blast(clusterNumber,type='microbial')
	MUSE.rerun_blast(clusterNumber,type)
end

# Reruns filtering for specified cluster number (by default 1 to 150)
def rerun_blast_filter(clusterNumber,type='microbial')
	MUSE.rerun_blast_filter(clusterNumber,type)
end

# Reruns microbial read extraction for specified cluster number (by default 1 to 150)
def rerun_extract_microbial_reads(clusterNumber)
	MUSE.rerun_extract_microbial_reads(clusterNumber)
end


list_commands = <<eos
# INITIAL ALIGNMENT
set_output_dir('/lustre/stlab/shanno01/unmapped_bam/muse3318')
set_unmapped_bam('/lustre/stlab/shanno01/unmapped_bam/3318.bam')
make_directories
add_bwa("/lustre/reference_data/stlab/shanno01/hs/hs_ref_GRCh37.p10.fa")
add_bwa("/lustre/reference_data/stlab/shanno01/hs/hs_alt_CHM1_1.0.fa")
add_bwa("/lustre/reference_data/stlab/shanno01/hs/hs_alt_HuRef.fa")
add_bwa("/lustre/reference_data/stlab/shanno01/other_hs/hs_GRCh37.70.cdna.fa")
add_bwa("/lustre/reference_data/stlab/shanno01/other_hs/hs_g1k_v37.fa")
add_bwa("/lustre/reference_data/stlab/shanno01/other_hs/hs_genomic.fa")
set_microbial_blast_db("/lustre/reference_data/stlab/shanno01/microbial/microbial")
set_human_blast_db("/lustre/reference_data/stlab/shanno01/blast/rna.fa")
start_bwa_jobs
start_repeat_masking
split_unmapped_bam
start_blast_jobs

# PROCESSING ALIGNMENT AND FILTERING MICROBIAL READS
import_bwa_results
import_repeat_masking
import_human_blast_results
write_queries
strip_blast_results
extract_microbial_reads

# ASSEMBLING AND ALIGNING CONTIGS
### Needs combining into one step
assemble_and_align_contigs
strip_contig_blast_results

# ADDITIONAL COMMANDS
[check_job_success] - not implemented
zip_results
show_query_summary
cleanup
eos

# = Simplified command line interface using defaults if run with arguments
if ARGV.length > 0
	require 'optparse'

	options = {}

	OptionParser.new do |opts|
		opts.banner = "Usage: musings_dsl.rb [options]"

		opts.on("-s", "--sample [SAMPLEID]", "Run musings on SAMPLEID (e.g. 3316)") do |sampleid|
			options[:sampleid] = sampleid
		end
	  
		opts.on("-r", "-run [SEGMENT]", "Run musings SEGMENT (align, queries, filter, contig)") do |segment|
			options[:segment] = segment
		end
		
		opts.on("-i", "--input [BAMFILE]", "Use BAMFILE as input") do |bamfile|
			options[:bamfile] = bamfile
		end
		
		opts.on("-o", "--output [FOLDER]", "Work and export results in FOLDER") do |outputFolder|
			options[:outputFolder] = outputFolder
		end
		
		opts.on("-c", "--command [COMMAND,COMMAND2,COMMAND3]", Array, "Call COMMANDs to muse about (run after any segment)") do |commands|
			options[:commands] = commands
		end
	  
		opts.on("-h", "--help", "Display help message") do |h|
			puts "Example usage: "
			puts "Setup and start alignments:           musings_dsl.rb -s 3316 -r \"align\""
			puts "Summarise queries:                    musings_dsl.rb -s 3316 -r \"queries\""
			puts "Read alignments and filter microbial: musings_dsl.rb -s 3316 -r \"filter\""
			puts "Assemble and align contigs:           musings_dsl.rb -s 3316 -r \"contig\""
			puts "Compress results to output.zip:		musings_dsl.rb -s 3316 -r \"zip\""
			puts ""
			puts "For DSL usage, use the following commands"
			puts "NB: These can be used via terminal as follwos musings_dsl.rb -s 3316 -c \"show_query_summary\""
			puts list_commands
			puts "" 
			exit
		end
	end.parse!

	# Use default settings to setup muse
	options[:outputFolder] ||= "/lustre/stlab/shanno01/unmapped_bam/muse#{options[:sampleid]}"
	#options[:outputFolder] ||= "/lustre/projects/stlab-icgc/dev/shanno01/unmapped_bams/muse#{options[:sampleid]}"
	set_output_dir(options[:outputFolder])
	options[:bamfile] ||= "/lustre/stlab/shanno01/unmapped_bam/#{options[:sampleid]}.bam"
	#options[:bamfile] ||= "/lustre/projects/stlab-icgc/dev/shanno01/unmapped_bams/#{options[:sampleid]}_unmapped.bam"
	set_unmapped_bam(options[:bamfile])
	make_directories
	add_bwa("/lustre/reference_data/stlab/shanno01/hs/hs_ref_GRCh37.p10.fa")
	add_bwa("/lustre/reference_data/stlab/shanno01/hs/hs_alt_CHM1_1.0.fa")
	add_bwa("/lustre/reference_data/stlab/shanno01/hs/hs_alt_HuRef.fa")
	add_bwa("/lustre/reference_data/stlab/shanno01/other_hs/hs_GRCh37.70.cdna.fa")
	add_bwa("/lustre/reference_data/stlab/shanno01/other_hs/hs_g1k_v37.fa")
	add_bwa("/lustre/reference_data/stlab/shanno01/other_hs/hs_genomic.fa")
	set_microbial_blast_db("/lustre/reference_data/stlab/shanno01/microbial/microbial")
	set_human_blast_db("/lustre/reference_data/stlab/shanno01/blast/rna.fa")

	
	# If called for alignments run alignment and masking jobs
	if options[:segment] == 'align'
		puts "Running segment: alignment"; $stdout.flush
		extract_unmapped_bam(options[:bamfile])
		start_bwa_jobs	
		split_unmapped_bam
		start_repeat_masking
		start_blast_jobs
	end
	
	# If called for query summarising, import alignments without filtering
	if options[:segment] == 'queries'
		puts "Running query summary"; $stdout.flush
		import_bwa_results
		import_repeat_masking
		import_human_blast_results
		import_microbial_blast_results
		write_queries
		seed_queries
		summarise_query_counts
	end

	# If called for processing alignments, parse alignments and masking and filter microbial reads
	if options[:segment] == 'filter'
		puts "Running segment: filtering"; $stdout.flush
		if does_queries_exist
			read_queries	
			strip_blast_results
		else
			import_bwa_results
			import_repeat_masking
			import_human_blast_results
			# import_microbial_blast_results # we don't need this step here as we get the data in strip_blast_results
			strip_blast_results
			write_queries
		end
		summarise_query_counts
		extract_microbial_reads
	end

	# If called for contig assembly and alignment, run contig commands, not currently implemented as one segment
	if options[:segment] == 'contig'
		assemble_and_align_contigs
		strip_contig_blast_results
		cleanup_contig_directory
	end
	
	# If called for compressing results, compress output folder
	if options[:segment] == 'zip'
		zip_results
	end
	
	if options[:segment] == 'dev'
		read_queries
		write_queries
		seed_queries
		summarise_query_counts
	end
	
	if options[:commands]
		options[:commands].each do |command|
			if command.match(/\(/)
				command, args = command.match(/(\w+)\((.*)\)/)[1..2]
				args = args.split(",")
				send(command,*args)
			else
				send(command)
			end
		end
	end
	
end
