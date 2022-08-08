#! /usr/bin/ruby
require 'optparse'

options = {}

OptionParser.new do |opts|
  opts.banner = "Usage: ExtractUnmapped.rb [options] querymatchfile1 querymatchfile2"

  opts.on("-s", "--seqfile [SEQFILE]", "Extract reads from the SEQFILE (.bam, .sam or .fasta format)") do |seqfile|
	options[:seqfile] = seqfile
  end
  
  opts.on("-o", "--outfile [OUTFILE]", "Output match count to OUTFILE") do |outfile|
	options[:outfile] = outfile
  end
  
 #opts.on("-h", "--help", "Display help message") do |h|
#	puts "Here's some advice"
	#exit
 #end
end.parse!

 
## make db from human
mapping = {}
ARGV.each do|f|
	File.foreach(f).collect do |strLine|
		splits = strLine.chomp.split("\t")
		if splits[2] == "No Hits."
		else
			mapping[splits[0]] = 1
		end
	end
end

out = File.open(options[:outfile],'w')

extension = File.extname(options[:seqfile])
if extension == ".bam"
	IO.popen("samtools view #{options[:seqfile]}").each do |readLine|
		name = readLine.chomp.split("\t")[0]
		if mapping[name]
		else
			out.puts readLine
		end
	end
elsif extension == ".sam"
	File.foreach(options[:seqfile]).collect do |readLine|
		name = readLine.chomp.split("\t")[0]
		if mapping[name]
		else
			out.puts readLine
		end
	end
else # extension == ".fasta" # or assume blank = fasta
	inf = File.open(options[:seqfile])
	until inf.eof?
		name = inf.readline
		seq = inf.readline
		if mapping[name.chomp.delete(">")]
		else
			out.puts name
			out.puts seq
		end
	end
	inf.close
end

out.close