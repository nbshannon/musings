#process query_sum

QUERYSUMS = []

class QuerySumContainer
attr_accessor :path, :bwas, :queries

def initialize(querySummaryFilename)
	@queries = []
	#= self.path[:output] + "query_summary.txt" 
	@path = querySummaryFilename
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
	@bwas = humandb
	sumQueries.each do |summary_key,count|
		human, masked, microbial = summary_key.split(":")
		humans = human.split(",")
		my_querySum = QuerySum.new
		@bwas.each do |bwa|
			my_querySum.bwa[bwa] = ((humans.include? bwa)? true : false)
		end
		my_querySum.hblast = ((humans.include? 'hblast')? true : false)
		my_querySum.masked = (masked == "masked"? true : false)
		my_querySum.microbial = (microbial == "microbial"? true : false)
		my_querySum.count = count.to_i
		@queries << my_querySum
	end
end

def puts_all
	puts [@bwas,'hblast','human','masked','microbial','count'].flatten.join("\t")
	puts @queries
end

def	microbial
	@queries.select{|q| q.microbial}
end


def puts_bwa_microbial
	bwa_counts = {}
	@bwas.each do |bwa| 
		bwa_counts[bwa] = self.microbial.map{|q| q.bwa[bwa]? q.count : 0}.inject(:+)
	end
	bwa_counts.each do |bwa,count|
		puts [bwa,count].join("\t")
	end
end

def puts_bwa_microbial_sequential
	puts "BWA\tUnmapped Microbial\tRemaining Unmapped Microbial"
	total = self.microbial.map{|q| q.count}.inject(:+)
	puts "Total\t#{total}\t#{total}"
	bwa_counts = {}
	@bwas.each do |bwa| 
		bwa_counts[bwa] = self.microbial.map{|q| q.bwa[bwa]? 0 : q.count}.inject(:+)
	end
	sorted = @bwas.sort_by{|b| bwa_counts[b]}
	used = []
	sorted.each do |sorted_bwa|
		used << sorted_bwa
		count = bwa_counts[sorted_bwa]
		culm = self.microbial.map{|q| q.has_any_bwa(used)? 0 : q.count}.inject(:+)
		puts [sorted_bwa,count,culm].join("\t")
	end
end

end

class QuerySum
attr_accessor :bwa, :hblast, :masked, :microbial, :query_key, :count

def initialize
	@bwa = {}
end

def has_all_bwa(bwas)
	ok = true
	bwas.each do |abwa|
		ok = ok & @bwa[abwa]
	end
	ok
end

def has_any_bwa(bwas)
	ok = false
	bwas.each do |abwa|
		ok = ok || @bwa[abwa]
	end
	ok
end

def human?
	[@bwa.values,@hblast].flatten.include? true
end

def to_s
	[@bwa.values,@hblast,self.human?,@masked,@microbial,@count].flatten.join("\t")
end

end

c = QuerySumContainer.new('c:/temp/muse/7504query_summary.txt')
#out = File.open(c.path + ".sum.txt",'w')
c.puts_bwa_microbial_sequential
puts "#############"
c.puts_all