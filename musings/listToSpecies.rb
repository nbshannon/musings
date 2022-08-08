file = 'c:\scripts\ruby\pathseq\microbial_list.txt'
out = File.open(file+".txt",'w')
File.foreach(file).collect do |line|
	splits = line.chomp.split("|")
	contig = splits[3]
	name = splits[4]
	begin
	nsplit = name.split(" ")
	if nsplit.length > 2
		name = nsplit[0..2].join(" ")
	elsif nsplit.length == 2
		name = nsplit[0..1].join(" ")
	end
	out.puts "#{contig}\t#{name}"
	rescue
	end
end