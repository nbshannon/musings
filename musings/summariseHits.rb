#! /usr/bin/ruby
# Summarise hits
inputFile = ARGV[0]; #// BLAST results
outputFile = ARGV[1];  #// output file

DBFILE = "combined_list.txt"
contigs = {}
File.foreach(DBFILE).collect do |strLine|
	contig, name, type = strLine.chomp.split("\t")
	contigs[contig] = {:name => name, :type => type}
end
out = File.open(outputFile,'w')
matchCount = {}
File.foreach(inputFile).collect do |strLine|
	match, query = strLine.chomp.split("\t")
	splits = match.split("|")
	if contigs[splits[3]]
		name = contigs[splits[3]][:name]
		matchCount[name] ||= 0
		matchCount[name] += 1
	end
end
matchCount.each do |name,count|
	out.puts "#{name}\t#{count}"
end