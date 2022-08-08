#! /usr/bin/ruby
# Blast to SP hits
inputFile = ARGV[0]; #// BLAST results
outputFile = ARGV[1];  #// output file
# ask for hits or match

# database to lookup contigs
DBFILE = "combined_list.txt"
contigs = {}
File.foreach(DBFILE).collect do |strLine|
	contig, name, type = strLine.chomp.split("\t")
	contigs[contig] = {:name => name, :type => type}
end

reads = {}
File.foreach(inputFile).collect do |strLine|
	name, read = strLine.chomp.split("\t")
	splits = name.split("|")
	if contigs[splits[3]]
		reads[read] ||= []
		reads[read] << contigs[splits[3]][:name]
	end
end

spCount = {}
reads.values.each do |contigs|
	contigs.uniq.each do |sp|
		spCount[sp] ||= 0
		spCount[sp] += 1
	end
end

spCount.each do |species, count|
	puts "#{species}\t#{count}"
end