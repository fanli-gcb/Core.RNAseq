#!/usr/bin/ruby

# parses a GTF file and computes gene lengths

# chr1    unknown exon    11874   12227   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS16107";
fn = ARGV.shift

gene_lengths = Hash.new(0)
File.open(fn).each_line do |line|
	chrom, source, type, start, stop, s1, strand, s2, comments = line.chomp.split(/\t/)
	if comments =~ /gene_id \"(\S+)\";/
		gene = $1
	else
		puts "Failed to find gene_id from #{line.chomp}"
		next
	end
	gene_lengths[gene] += (stop.to_i - start.to_i + 1)
end

gene_lengths.each_key { |gene|
	puts "#{gene}\t#{gene_lengths[gene]}"
}


