#!/usr/bin/ruby

if ARGV.size == 2
	kallisto_fn, mapping_fn = ARGV
else
	puts "USAGE: #{$0} kallisto_fn mapping_fn"
	exit -1
end

# store transcript->gene mapping
tx2gene = Hash.new
genes = Hash.new
File.open(mapping_fn).each_line do |line|
	tx, gene = line.chomp.split(/\t/)
	tx2gene[tx] = gene
	genes[gene] = 1
end

# aggregate counts to gene level
gene_level_counts = Hash.new(0)
first = true
File.open(kallisto_fn).each_line do |line|
	if first
		first = false
		next
	end
	
	tx, length, eff_length, est_count, tpm = line.chomp.split(/\t/)
	gene = tx2gene[tx]
	next if gene.nil?
	gene_level_counts[gene] += est_count.to_f
end

genes.keys.sort.each { |gene|
	if !gene.nil?
		puts "#{gene}\t#{gene_level_counts[gene]}"
	end
}

