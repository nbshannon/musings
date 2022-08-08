# MusiNGS module for ruby
# Run alignments and process results using the computer cluster
# Run via the musings_dsl

NUMREPEATS = 2
REPEATFILTER = "N" * NUMREPEATS
CLUSTERDEFAULT = 150

# -f 76, read and mate unmapped, first in pair
# -F 1536, exclude reads that fail platform/vendor quality checks, or are PCR or optical duplicates
SAMTOOLSUNMAPPED = "samtools view -b -f 76 -F 1536 "

# -F, exclude unmapped reads
SAMTOOLSMAPPED = "samtools view -S -F4 "
  
# Test for directory existence and ability to make
def testDir(dir)
	if FileTest::directory?(dir)
	else
		begin
			Dir::mkdir(dir)
		rescue
			abort "Cannot make directory: #{dir}"
		end
	end
end



# Wrapper for muse functionality
module Musings

	# container for run information, for the purpose of functions see musings_dsl
	class Muse
	attr_accessor :queries, :outputDir, :unmappedBam, :batchSize, :bwaDbs, :microbialDb, :humanBlastDb
	
	# Initialise muse with default batchSize of 150
	def initialize
		@batchSize = CLUSTERDEFAULT
	end
	
	def set_output_dir(outputDir)
		testDir(outputDir)
		@outputDir = outputDir
	end
	
	def set_unmapped_bam(unmappedBam)
		@unmappedBam = unmappedBam
	end
	
	# Use command in SAMTOOLSUNMAPPED to extract unmapped bam
	# samtools view -b -f 76 -F 1536
	def extract_unmapped_bam(someMappedBam)
		%x[#{SAMTOOLSUNMAPPED} #{someMappedBam} > #{someMappedBam}.unmapped.bam]; result=$?.success?
		self.set_unmapped_bam("#{someMappedBam}.unmapped.bam")
	end
	
	def seed_queries
		@queries ||= {}
		puts "Seeding queries from #{@unmappedBam}.."
		bam = IO.popen("samtools view #{@unmappedBam}")
			totalReadCount = 0
			until bam.eof?
				totalReadCount += 1
				query = bam.readline.chomp.split("\t")[0]
				@queries[query] ||= {:human=>[],:masked=>false,:microbial=>false}
			end
		bam.close
		puts ".. seeded #{totalReadCount} reads"
	end
	
	def set_batch_size(clusterSize)
		@batchSize = clusterSize.to_i
	end
	
	def add_bwa(bwaDb)
		@bwaDbs ||= []
		@bwaDbs << bwaDb
	end
	
	# Check that the blast database exists (.nhr, or .00.nhr if split)
	def set_microbial_blast_db(microbialDb)
		if File.exist?("#{microbialDb}.nhr") or File.exist?("#{microbialDb}.00.nhr")
			@microbialDb = microbialDb
		else
			abort "cannot find microbial blast db file: #{microbialDb}.nhr or #{microbialDb}.00.nhr"
		end
	end
	
	# Check that the blast database exists (.nhr, or .00.nhr if split)
	def set_human_blast_db(humanBlastDb)
		if File.exist?("#{humanBlastDb}.nhr") or File.exist?("#{humanBlastDb}.00.nhr")
			@humanBlastDb = humanBlastDb
		else
			abort "cannot find human blast db: #{humanBlastDb}.nhr or #{humanBlastDb}.00.nhr"
		end
	end
	
	# Directory paths used inside output directory
	def path
		{
		:root => @outputDir + "/",
		:bwa => @outputDir + "/bwa/",
		#:bwaf => @outputDir + "/bwa/bwa_",
		:blast => @outputDir + "/split/",
		#:blastf => @outputDir + "/split/split",
		:repeat => @outputDir + "/masked/",
		:output => @outputDir + "/output/",
		:report => @outputDir + "/output/report/",
		:filteredReads => @outputDir + "/output/fasta/",
		:filteredHits => @outputDir + "/output/hits/",
		#:filteredHitsf => @outputDir + "/output/hits/split",
		:contig => @outputDir + "/output/contig/",
		:contigHits => @outputDir + "/output/contighits/",
		}
	end

	def make_directories
		[self.path[:root],self.path[:bwa],self.path[:blast],self.path[:repeat],self.path[:output],self.path[:report],self.path[:filteredReads],self.path[:filteredHits],self.path[:contig],self.path[:contigHits]].each do |dir|
			testDir(dir)
		end
	end

	def split_unmapped_bam
		nreads= %x[samtools view #{@unmappedBam} | wc -l]; result=$?.success?
		if result
			chunksize = (nreads.strip.to_f/@batchSize).ceil
			bam = IO.popen("samtools view #{@unmappedBam}")
			splitNo = 0
			until bam.eof?
				splitNo += 1
				out = File.open("#{self.path[:blast]}split#{splitNo}.fa",'w')
				(1..chunksize).each do |chunk|
					next if bam.eof?
					splits = bam.readline.chomp.split("\t")
					out.puts ">#{splits[0]}\n#{splits[9]}"
				end
				out.close
			end
			bam.close
		else
			puts "Cannot count number of reads in #{@unmappedBam}"
		end
	end
	
	# bwa aln: Find the coordinates of the input reads (default options, -b to indicate bam)
	# bwa samse: Generate alignments in sam format
	# samtools view: extract mapped reads (-F, discard unmapped)
	def start_bwa_jobs
		@bwaDbs.each do |bwaDb|
			i = File.basename(bwaDb)
			bwarun = %x[bsub -o #{self.path[:report]}bwa_#{i}.out.txt -e #{self.path[:report]}bwa_#{i}.err.txt -Rrusage[mem=20240] "bwa aln -b #{bwaDb} #{@unmappedBam} | bwa samse #{bwaDb} - #{@unmappedBam} | samtools view -S -F4 - > #{self.path[:bwa]}bwa_#{i}.sam"]
			puts "..Bwa job (#{i}) submitted as: #{bwarun.match(/<([0-9]+)>/)[1]}"; $stdout.flush
		end
	end
	
	def start_repeat_masking
		repeatmasker = %x[bsub -o #{self.path[:report]}repeat_out.txt -e #{self.path[:report]}repeat_err.txt -Rrusage[mem=10240] -J "sort[1-#{@batchSize}]%50" "/lustre/stlab/shanno01/bin/repeatmasker/RepeatMasker/RepeatMasker -pa 4 -species vertebrates -dir #{self.path[:repeat]} #{self.path[:blast]}split\\$LSB_JOBINDEX.fa"]; ; rresult=$?.success?
		puts "..RepeatMasker job submitted as: #{repeatmasker.match(/<([0-9]+)>/)[1]}[1-#{@batchSize}]"; $stdout.flush
	end
	
	# -task megablast
	# -num_threads 4 
	# -evalue 0.0000001 (expectation value (E) threshold for saving hits)
	# -word_size 16 (word size for wordfinder algorithm (length of best perfect match)
	# -max_target_seqs 5 (number of aligned sequences to keep)
	# -outfmt 6 (tabular output format)
	# -dust no (disable filtering of query sequence with DUST)
	def start_blast_jobs
		if @microbialDb
			mblast = %x[bsub -o #{self.path[:report]}blast_m_out.txt -e #{self.path[:report]}blast_m_err.txt -Rrusage[mem=10240] -J "sort[1-#{@batchSize}]%50" "blastn -task megablast -query #{self.path[:blast]}split\\$LSB_JOBINDEX.fa -db #{@microbialDb} -num_threads 4 -evalue 0.0000001 -word_size 16 -max_target_seqs 5 -outfmt 6 -dust no -out #{self.path[:blast]}split\\$LSB_JOBINDEX.mout"]; mresult=$?.success?
			puts "..microbial blast alignment job submitted as: #{mblast.match(/<([0-9]+)>/)[1]}[1-#{@batchSize}]"; $stdout.flush
		end
		if @humanBlastDb
			hblast = %x[bsub -o #{self.path[:report]}blast_h_out.txt -e #{self.path[:report]}blast_h_err.txt -Rrusage[mem=10240] -J "sort[1-#{@batchSize}]%50" "blastn -task megablast -query #{self.path[:blast]}split\\$LSB_JOBINDEX.fa -db #{@humanBlastDb} -num_threads 4 -evalue 0.0000001 -word_size 16 -max_target_seqs 5 -outfmt 6 -dust no -out #{self.path[:blast]}split\\$LSB_JOBINDEX.hout"]; hresult=$?.success?
			puts "..human blast alignment job submitted as: #{hblast.match(/<([0-9]+)>/)[1]}[1-#{@batchSize}]"; $stdout.flush
		end
	end
	
	def rerun_blast(clusterNumber,type)
		puts "Rerunning #{type} blast alignment for split#{clusterNumber}.fa"; $stdout.flush
		if type == 'microbial'
			if @microbialDb
				mblast = %x[bsub -o #{self.path[:report]}blast_m_#{clusterNumber}_out.txt -e #{self.path[:report]}blast_m_#{clusterNumber}_err.txt -Rrusage[mem=20240] "blastn -task megablast -query #{self.path[:blast]}split#{clusterNumber}.fa -db #{@microbialDb} -num_threads 4 -evalue 0.0000001 -word_size 16 -max_target_seqs 5 -outfmt 6 -dust no -out #{self.path[:blast]}split#{clusterNumber}.mout"]; mresult=$?.success?
				puts "..microbial blast alignment job submitted as: #{mblast.match(/<([0-9]+)>/)[1]}"; $stdout.flush
			else
				puts "I don't know a microbial database"; $stdout.flush
			end
		elsif type == 'human'
			if @humanBlastDb
				mblast = %x[bsub -o #{self.path[:report]}blast_m_#{clusterNumber}_out.txt -e #{self.path[:report]}blast_m_#{clusterNumber}_err.txt -Rrusage[mem=20240] "blastn -task megablast -query #{self.path[:blast]}split#{clusterNumber}.fa -db #{@humanBlastDb} -num_threads 4 -evalue 0.0000001 -word_size 16 -max_target_seqs 5 -outfmt 6 -dust no -out #{self.path[:blast]}split#{clusterNumber}.mout"]; mresult=$?.success?
				puts "..human blast alignment job submitted as: #{mblast.match(/<([0-9]+)>/)[1]}"; $stdout.flush
			else
				puts "I don't know a human blast database"; $stdout.flush
			end
		end
	end
	
	def rerun_blast_filter(clusterNumber, type)
		@queries ||= {}
		puts "Rerunning filtering for split#{clusterNumber}"; $stdout.flush
		if type == 'human'
			puts "..Importing human blast results"; $stdout.flush
			File.foreach("#{self.path[:blast]}split#{clusterNumber}.hout").collect do |blastLine|
				query = blastLine.split("\t")[0]
				@queries[query] ||= {:human=>[],:masked=>false,:microbial=>false}
				if @queries[query].include? 'hblast'
				else
					@queries[query][:human] << 'hblast'
				end
			end
		end
		puts "..Stripping blast results"
		out = File.open("#{self.path[:filteredHits]}split#{clusterNumber}.out",'w')
		File.foreach("#{self.path[:blast]}split#{clusterNumber}.mout").collect do |blastLine|
			query = blastLine.split("\t")[0]
			@queries[query] ||= {:human=>[],:masked=>false,:microbial=>false}
			@queries[query][:microbial] = true
			if (@queries[query][:human].length > 0) or (@queries[query][:masked])
			else
				out.puts blastLine
			end
		end
		out.close
	end
	
	def rerun_extract_microbial_reads(clusterNumber)
		puts "Extract microbial reads for split#{clusterNumber}"; $stdout.flush
		microbialQueries = {}
		File.foreach("#{self.path[:filteredHits]}split#{clusterNumber}.out").collect do |blastLine|
			query = blastLine.split("\t")[0]
			microbialQueries[query] = 1
		end
		out = File.open("#{self.path[:filteredReads]}split#{clusterNumber}.fa",'w')
		File.foreach("#{self.path[:blast]}split#{clusterNumber}.fa",sep=">").collect do |fastaLine|
			splits =  fastaLine.split("\n")
			query = splits.shift
			seq = splits.join("")
			if microbialQueries[query]
				out.puts ">" + query
				out.puts seq
			end
		end
		out.close
	end

	def import_bwa_results
		puts "Importing BWA results"; $stdout.flush
		@queries ||= {}
		@bwaDbs.each do |bwaDb|
			i = File.basename(bwaDb)
			bwaMatches = 0
			print ".. reading #{i} .."; $stdout.flush
			File.foreach("#{self.path[:bwa]}bwa_#{i}.sam").collect do |samLine|
				query = samLine.split("\t")[0]
				bwaMatches += 1
				@queries[query] ||= {:human=>[],:masked=>false,:microbial=>false}
				@queries[query][:human] << "bwa_#{i}"
				query = nil
				samLine = nil
			end
			puts "#{bwaMatches} mapped reads"; $stdout.flush
		end
		puts "..found #{@queries.values.select{|q| q[:human].length > 0 ? true : false}.length} human reads in total"; $stdout.flush
	end
	
	def import_repeat_masking
		puts "Importing repeat masking"; $stdout.flush
		@queries ||= {}
		repeatMasked = 0
		(1..@batchSize).each do |batch|
			#print "\r..reading #{batch}/#{@batchSize}"; $stdout.flush # removed for verbosity of LSF output
			if File.exist?("#{self.path[:repeat]}split#{batch}.fa.masked")
				File.foreach("#{self.path[:repeat]}split#{batch}.fa.masked",sep=">").collect do |repeatLine|
					splits =  repeatLine.split("\n")
					query = splits.shift
					#splits.pop # remove ">"
					seq = splits.join("")
					if seq.match(REPEATFILTER)
						repeatMasked += 1
						@queries[query] ||= {:human=>[],:masked=>false,:microbial=>false}
						@queries[query][:masked] = true
					end
					query = nil
					seq = nil
					repeatLine = nil
				end
			end
		end	
		puts "\r..found #{repeatMasked} masked reads"; $stdout.flush
	end
	
	def import_human_blast_results
		puts "Importing human blast results"; $stdout.flush
		@queries ||= {}
		humanBlastMatches = 0
		(1..@batchSize).each do |batch|
			#print "\r..reading #{batch}/#{@batchSize}"; $stdout.flush # removed for verbosity of LSF output
			File.foreach("#{self.path[:blast]}split#{batch}.hout").collect do |blastLine|
				query = blastLine.split("\t")[0]
				humanBlastMatches += 1
				@queries[query] ||= {:human=>[],:masked=>false,:microbial=>false}
				if @queries[query].include? 'hblast'
				else
					@queries[query][:human] << 'hblast'
				end
				query = nil
				blastLine = nil
			end
		end
		humanReads = @queries.values.select{|q| q[:human].include? 'hblast'}.length
		puts "\r..found #{humanBlastMatches} mappings for #{humanReads} reads"; $stdout.flush
	end
	
	def import_microbial_blast_results
		puts "Importing microbial blast results"; $stdout.flush
		queries ||= {}
		processedMatches = 0
		(1..@batchSize).each do |batch|
			File.foreach("#{self.path[:blast]}split#{batch}.mout").collect do |blastLine|
				query = blastLine.split("\t")[0]
				@queries[query] ||= {:human=>[],:masked=>false,:microbial=>false}
				@queries[query][:microbial] = true
				processedMatches += 1
				query = nil
				blastLine = nil
			end
		end
		microbialReads = @queries.values.select{|q| q[:microbial] ? true : false}.length
		puts "\r.. found #{processedMatches} mapping for #{microbialReads} microbial reads"; $stdout.flush
	end
	
	def strip_blast_results
		# unfinished want total input count
		puts "Stripping blast results"; $stdout.flush
		queries ||= {}
		# strip_hash = {} # build this prior to queries? only need to look up each query once
		filteredMatches = 0
		
		excludedReads = @queries.values.select{|q| (q[:human].length > 0) or (q[:masked])}.length
		humanReads = @queries.values.select{|q| q[:human].length > 0 ? true : false}.length
		maskedReads = @queries.values.select{|q| q[:masked] ? true : false}.length
		puts ".. excluding a total of: #{excludedReads} reads (#{humanReads} human, #{maskedReads} masked)"
		
		(1..@batchSize).each do |batch|
			#print "\r..processing #{batch}/#{@batchSize}"; $stdout.flush # removed for verbosity of LSF output
			out = File.open("#{self.path[:filteredHits]}split#{batch}.out",'w')
			File.foreach("#{self.path[:blast]}split#{batch}.mout").collect do |blastLine|
				query = blastLine.split("\t")[0]
				@queries[query] ||= {:human=>[],:masked=>false,:microbial=>false}
				@queries[query][:microbial] = true
				if (@queries[query][:human].length > 0) or (@queries[query][:masked])
				else
					out.puts blastLine
					filteredMatches += 1
				end
				query = nil
				blastLine = nil
				end
			out.close
		end
		microbialReads = @queries.values.select{|q| q[:microbial] ? true : false}.length
		microbialFiltered = @queries.values.select{|q| (q[:microbial]) and !((q[:human].length > 0) or (q[:masked])) }.length
		puts "\r.. total of #{filteredMatches} remaining matches processed from #{microbialFiltered} filtered (of #{microbialReads} microbial reads)"; $stdout.flush
	end
	
	def does_queries_exist
		queryFilename = self.path[:root] + "queries.txt"
		File.exist?(queryFilename)
	end
	
	def read_queries
		puts "Reading queries"; $stdout.flush
		@queries ||= {}
		queryFilename = self.path[:root] + "queries.txt"
		if File.exist?(queryFilename)
			File.foreach(queryFilename).collect do |queryLine|
				query, human, masked, microbial = queryLine.chomp.split("\t")
				@queries ||= {}
				@queries[query] ||= {:human=>[],:masked=>false,:microbial=>false}
				@queries[query][:human] = human.split(",").uniq
				@queries[query][:masked] = masked == "true"
				@queries[query][:microbial] = microbial == "true"
			end
		end
		humanReads = @queries.values.select{|q| q[:human].length > 0 ? true : false}.length
		maskedReads = @queries.values.select{|q| q[:masked] ? true : false}.length
		puts ".. found #{@queries.length} reads, with #{humanReads} human reads, and #{maskedReads} masked reads"
	end
	
	def write_queries
		queryFilename = self.path[:root] + "queries.txt"
		queryFile = File.open(queryFilename,'w')
		@queries.each do |query,match|
			queryFile.puts [query,match[:human].join(","),match[:masked],match[:microbial]].join("\t")
		end
		queryFile.close
	end
	
	def summarise_query_counts
		puts "Summarising query counts"; $stdout.flush
		summary = {}
		@queries.each do |query,q|
			summary_key = "#{q[:human].length > 0 ? q[:human].join(",") : 'non_human'}:#{q[:masked] ? 'masked' : 'non_masked'}:#{q[:microbial]? 'microbial' : 'non_microbial'}"
			summary[summary_key] ||= 0
			summary[summary_key] += 1
		end
		querySummaryFilename = self.path[:output] + "query_summary.txt"
		out = File.open(querySummaryFilename,'w')
		summary.each do |summary_key,count|
			out.puts [summary_key,count].join("\t")
		end
		puts ".. summarised #{@queries.length} queries"
	end
	
	def export_muse
		museFilename = self.path[:root] + "muse.dat"
		museFile = File.open(museFilename,'w')
		Marshal::dump(self, museFile)
	end
	
	def load_muse(museFilename=nil)
		if museFilename == nil
			museFilename = self.path[:root] + "muse.dat"
		end
		museFile = File.open(museFilename)
		rmuse = Marshal::load(museFile)
		@queries = rmuse.queries
		@outputDir = rmuse.outputDir
		@unmappedBam = rmuse.unmappedBam
		@batchSize = rmuse.batchSize
		@bwaDbs = rmuse.bwaDbs
		@microbialDb = rmuse.microbialDb
		@humanBlastDb = rmuse.humanBlastDb
		close museFile
	end
	
	def extract_microbial_reads
		puts "Extract microbial reads"; $stdout.flush
		microbialQueries = {}
		processedMicrobial = 0
		processedDiscarded = 0
		
		puts ".. reading microbial output"; $stdout.flush
		(1..@batchSize).each do |batch|
			#print "\r..processing #{batch}/#{@batchSize}"; $stdout.flush # removed for verbosity of LSF output
			File.foreach("#{self.path[:filteredHits]}split#{batch}.out").collect do |blastLine|
				query = blastLine.split("\t")[0]
				microbialQueries[query] = 1
			end
		end
		puts ".. found #{microbialQueries.length} microbial reads to extract"; $stdout.flush
		
		puts ".. extracting microbial reads"; $stdout.flush
		(1..@batchSize).each do |batch|
			#print "\r..processing #{batch}/#{@batchSize}"; $stdout.flush # removed for verbosity of LSF output
			out = File.open("#{self.path[:filteredReads]}split#{batch}.fa",'w')
			File.foreach("#{self.path[:blast]}split#{batch}.fa",sep=">").collect do |fastaLine|
				splits =  fastaLine.split("\n")
				query = splits.shift
				seq = splits.join("")
				if microbialQueries[query]
					out.puts ">" + query
					out.puts seq
					processedMicrobial += 1
				else
					processedDiscarded += 1
				end
			end
			out.close
		end
		puts "\r.. total of #{processedMicrobial} extracted reads, #{processedDiscarded} discarded reads"; $stdout.flush
	end
	
	def assemble_and_align_contigs
		puts "Building contigs"; $stdout.flush
		puts "..constructing velvet dataset (velveth)"; $stdout.flush
		velveth = %x[velveth #{self.path[:contig]} 19 -fasta -short #{self.path[:filteredReads]}*.fa]; hresult=$?.success?
		abort "Problem constructing velvet dataset (velveth)" unless hresult
		puts "..generating contigs (velvetg)"; $stdout.flush
		velvetg = %x[velvetg #{self.path[:contig]} -min_contig_lgth 75]; gresult=$?.success?
		abort "Problem generating contigs (velvetg)" unless gresult		
		if @humanBlastDb
			puts "Aliging against #{@humanBlastDb}"; $stdout.flush
			hblast = %x[blastn -task megablast -query #{self.path[:contig]}contigs.fa -db #{@humanBlastDb} -num_threads 4 -evalue 0.0000001 -word_size 16 -max_target_seqs 5 -outfmt 6 -dust no -out #{self.path[:contig]}contigs.hout]; hresult=$?.success?
		end
		abort "Problem aligning against #{@humanBlastDb}" unless hresult
		if @microbialDb
			puts "Aligning against #{@microbialDb}"; $stdout.flush
			mblast = %x[blastn -task megablast -query #{self.path[:contig]}contigs.fa -db #{@microbialDb} -num_threads 4 -evalue 0.0000001 -word_size 16 -max_target_seqs 5 -outfmt 6 -dust no -out #{self.path[:contig]}contigs.mout]; mresult=$?.success?
		end
		abort "Problem aligning against #{@microbialDb}" unless mresult
		puts "Finished aligning contigs"; $stdout.flush
	end
	
	def build_contigs # deprecated in favour of running all the contig steps in one
		puts "Building contigs"
		puts "..constructing velvet dataset (velveth)"; $stdout.flush
		velveth = %x[velveth #{self.path[:contig]} 19 -fasta -short #{self.path[:filteredReads]}*.fa]; hresult=$?.success?
		if hresult
			puts "..generating contigs (velvetg)"; $stdout.flush
			velvetg = %x[bsub -o #{self.path[:report]}velvetg_out.txt -e #{self.path[:report]}velvetg_err.txt -Rrusage[mem=10240] "velvetg #{self.path[:contig]} -min_contig_lgth 75"]; gresult=$?.success?
			puts "..contig assembly (velvetg) job submitted as: #{velvetg.match(/<([0-9]+)>/)[1]}"; $stdout.flush
		end
	end
	
	def start_contig_blast_jobs # deprecated in favour of running all the contig steps in one
		if @microbialDb
			mblast = %x[bsub -o #{self.path[:report]}cblast_m_out.txt -e #{self.path[:report]}cblast_m_err.txt -Rrusage[mem=10240] "blastn -task megablast -query #{self.path[:contig]}contigs.fa -db #{@microbialDb} -num_threads 4 -evalue 0.0000001 -word_size 16 -max_target_seqs 5 -outfmt 6 -dust no -out #{self.path[:contig]}contigs.mout"]; mresult=$?.success?
			puts "..microbial blast contig alignment job submitted as: #{mblast.match(/<([0-9]+)>/)[1]}"; $stdout.flush
		end
		if @humanBlastDb
			hblast = %x[bsub -o #{self.path[:report]}cblast_h_out.txt -e #{self.path[:report]}cblast_h_err.txt -Rrusage[mem=10240] "blastn -task megablast -query #{self.path[:contig]}contigs.fa -db #{@humanBlastDb} -num_threads 4 -evalue 0.0000001 -word_size 16 -max_target_seqs 5 -outfmt 6 -dust no -out #{self.path[:contig]}contigs.hout"]; hresult=$?.success?
			puts "..human blast contig alignment job submitted as: #{hblast.match(/<([0-9]+)>/)[1]}"; $stdout.flush
		end
	end
	
	def strip_contig_blast_results
		puts "Stripping contig blast results"; $stdout.flush
		humanBlastMatches = 0
		uniqHumanMatches = {}
		processedMatches = 0
		uniqProcessedContigs = {}
		
		print "..reading contig human blast alignment"; $stdout.flush
		File.foreach("#{self.path[:contig]}contigs.hout").collect do |blastLine|
			query = blastLine.split("\t")[0]
			humanBlastMatches += 1
			uniqHumanMatches[query] = 1
		end
		puts "\r..found #{humanBlastMatches} human mappings for #{uniqHumanMatches.keys.length} contigs"; $stdout.flush
		
		#print "\r..processing contig microbial blast alignment"; $stdout.flush
		out = File.open("#{self.path[:contigHits]}contigs.out",'w')
		File.foreach("#{self.path[:contig]}contigs.mout").collect do |blastLine|
			query = blastLine.split("\t")[0]
			if uniqHumanMatches[query]
			else
				out.puts blastLine
				processedMatches += 1
				uniqProcessedContigs[query] = 1
			end
		end
		out.close
		puts "\r.. total of #{processedMatches} remaining matches processed for #{uniqProcessedContigs.keys.length} contigs"; $stdout.flush
	end
	
	def cleanup_contig_directory
		puts "deleting temporary files from contig directory"; $stdout.flush
		%x[rm `ls #{self.path[:contig]}* | grep -v contigs.fa`]; zresult=$?.success?
	end
	
	def zip_results
		print "compressing output directory ..."; $stdout.flush
		zip = %x[zip -r #{self.path[:root]}output #{self.path[:output]}]; zresult=$?.success?
		puts "complete"; $stdout.flush
	end
	
	def show_query_summary
		querySummaryFilename = self.path[:output] + "query_summary.txt" 
		categories = []
		sumQueries = {}
		File.foreach(querySummaryFilename).collect do |summaryLine|
			summary_key, count = summaryLine.chomp.split("\t")
			sumQueries[summary_key] = count
		end
		humandb = []
		sumQueries.keys.each do |summary_key|
			human, masked, microbial = summary_key.split(":")
			humans = human.split(",")
			humandb << humans
		end
		humandb = humandb.flatten.uniq
		humandb.delete('hblast')
		humandb.delete('non_human')
		categories << humandb
		categories << ['hblast','human','masked','microbial','count']
		puts categories.flatten.join("\t")
		sumQueries.each do |summary_key,count|
			outputs = []
			human, masked, microbial = summary_key.split(":")
			humans = human.split(",")
			humandb.each do |bwa|
				outputs << ((humans.include? bwa)? "Mapped" : "Unmapped")
			end
			outputs << ((humans.include? 'hblast')? "Mapped" : "Unmapped")
			outputs << ((outputs.include? 'Mapped')? "Mapped" : "Unmapped")
			outputs << masked
			outputs << microbial
			outputs << count
			puts outputs.join("\t")
		end
	end
	
	def cleanup
		puts "Carrying out directory cleanup"
		puts ".. deleting bwa alignments directory"; $stdout.flush
		%x[rm -r #{self.path[:bwa]}]
		puts ".. deleting blast alignments directory"; $stdout.flush
		%x[rm -r #{self.path[:blast]}]
		puts ".. deleting repeat masking directory"; $stdout.flush
		%x[rm -r #{self.path[:repeat]}]
	end
end

end
