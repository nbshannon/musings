#! /usr/bin/ruby
require 'optparse'

options = {}

OptionParser.new do |opts|
  opts.banner = "Usage: StripRepeasts.rb [options]"

  opts.on("-i", "--infile [INFILE]", "INFILE for removing masked reads (.masked)") do |infile|
	options[:infile] = infile
  end
  
  opts.on("-o", "--outfile [OUTFILE]", "Output sequences to OUTFILE") do |outfile|
	options[:outfile] = outfile
  end
  
  opts.on("-n", "--numrepeats [NUMREPEATS]", OptionParser::DecimalInteger, "NUMREPEATS to exclude sequences on") do |numrepeats|
	options[:numrepeats] = numrepeats
end
  
 #opts.on("-h", "--help", "Display help message") do |h|
#	puts "Here's some advice"
	#exit
 #end
end.parse!

if options[:outfile]
	out = File.open(options[:outfile],'w')
else
	out = $stdout
end 

options[:filter] = "N" * options[:numrepeats]

puts "Excluding sequences with #{options[:filter]}, (#{options[:numrepeats]} repeats)"

inf = File.open(options[:infile])
nextline = inf.readline.chomp
until inf.eof?
	name = nextline
	nextline = inf.readline.chomp
	seq = ""
	until nextline.match(/^>/)
		seq = seq + nextline
		if inf.eof?
			nextline = ">end"
		else
		nextline = inf.readline.chomp
		end
	end
	if seq.match(options[:filter])
	else
		out.puts name
		out.puts seq
	end
end

inf.close
out.close