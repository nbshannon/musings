#! /usr/bin/ruby
#@ARGV = ['match','c:/temp/s2.mout','c:/temp/s2.txt']
type = ARGV[0]
inputFile = ARGV[1]; #// BLAST results
outputFile = ARGV[2];  #// output file
# ask for hits or match


start = Time.now;


totalReadCounter 		= 0;
unmappedReadCounter 	= 0;
mappedReadCounter		= 0;

QUERY  = "Query=";   	   #//	Pattern queryLine= Pattern.compile(QUERY);
NOHIT  = "No hits found";  #//Pattern noHits= Pattern.compile(NOHIT);
HIT    = "  gi\\|";

#br = File.open(inputFile);
out = File.open(outputFile,'w');
curRead = "";
strLine = "";
#while (strLine = br.readline)
File.foreach(inputFile).collect do |strLine|
	if (strLine.match(QUERY))
		if (totalReadCounter == 0)
			column = strLine.split(" ");
			if type== 'query'
				out.print(column[1] + "\t=");
			elsif type== 'match'
				curRead = column[1];
			end
			totalReadCounter+=1 ;
		else
			column = strLine.split(" ");
			if type== 'query'
				out.print("\n" + column[1] + "\t=");
			elsif type== 'match'
				curRead = column[1];
			end
			totalReadCounter+=1 ;
		end
	end
	if (strLine.match(NOHIT))
		if type== 'query'
			out.print( "\tNo Hits." );
		end
		unmappedReadCounter+=1;
	end
	if (strLine.match(HIT))
		column = strLine.split("  ");
		if type== 'match'
			out.print(column[1] + "\t" + curRead + "\n");
		end
	end
end
#br.close;
mappedReadCounter= totalReadCounter - unmappedReadCounter;
elapsedTimeSec = Time.now - start;
$stdout.print( "Process complete. \n\n" +
	"Total reads: " + totalReadCounter.to_s + 
	"\nMapped reads: " + mappedReadCounter.to_s + 
	"\nUnmapped reads: " + unmappedReadCounter.to_s + "\n\n" +
	"CPU time: "+elapsedTimeSec.to_s);
out.flush
out.close
