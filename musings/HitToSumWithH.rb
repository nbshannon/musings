#! /usr/bin/ruby
require 'optparse'
require 'ostruct'
# Blast to SP hits
inputFile = ARGV[0]; #// BLAST results
inputHuman = ARGV[1]; #// Human hits
outputFile = ARGV[2];  #// output file
# ask for hits or match
#ARGV = ['--help','hdilw']

options = {}

OptionParser.new do |opts|
  opts.banner = "Usage: HitToSumWithH.rb [options]"

  #opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
  #  options[:verbose] = v
  #end
  
  opts.on("-h", "--human [HUMANHITS]", "Use the HUMANHITS file to filter hits") do |hitfile|
	options[:human] = hitfile
  end
  
  opts.on("-m", "--microbial [MICROBIALHITS]", "Use the MICROBIALHITS file as a list of microbial hits") do |hitfile|
    options[:microbial] = hitfile
  end
  
  opts.on("-o", "--outfile [OUTFILE]", "Output match count to OUTFILE") do |outfile|
	options[:outfile] = outfile
  end
  
  opts.on("-d", "--db [LISTFILE]", "LISTFILE containing sequence ID to species mapping") do |listfile|
    options[:listfile] = listfile
  end
  
 #opts.on("-h", "--help", "Display help message") do |h|
#	puts "Here's some advice"
	#exit
 #end
end.parse!

### add some checks that we have options and that relevant files exist


# database to lookup contigs
contigs = {}
File.foreach(options[:listfile]).collect do |strLine|
	contig, name, type = strLine.chomp.split("\t")
	contigs[contig] = {:name => name, :type => type}
end

reads = {}
File.foreach(options[:microbial]).collect do |strLine|
	name, read = strLine.chomp.split("\t")
	splits = name.split("|")
	if contigs[splits[3]]
		reads[read] ||= []
		reads[read] << contigs[splits[3]][:name]
	end
end

if options[:human]
## make db from human
mapping = {}
File.foreach(options[:human]).collect do |strLine|
	splits = strLine.chomp.split("\t")
	if splits[2] == "No Hits."
	else
		mapping[splits[0]] = 1
	end
end
## strip human mapping reads
reads.each do |read|
	if mapping[read]
		reads.delete(read)
	end
end
end

spCount = {}
spCountUniq = {}
reads.values.each do |contigs|
	contigs.uniq.each do |sp|
		spCount[sp] ||= 0
		spCount[sp] += 1
		spCountUniq[sp] ||= 0
		if contigs.uniq.length == 1
			spCountUniq[sp] += 1
		end
	end
end

if options[:outfile]
	out = File.open(options[:outfile],'w')
	spCount.each do |species, count|
		out.puts "#{species}\t#{count}\t#{spCountUniq[species]}"
	end
	out.close
else
	spCount.each do |species, count|
		puts "#{species}\t#{count}\t#{spCountUniq[species]}"
	end
end