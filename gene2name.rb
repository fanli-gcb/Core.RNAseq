#!/usr/bin/ruby

# parses a GTF file and makes a mapping from gene_id to gene_name

# chr1    unknown exon    11874   12227   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS16107";
fn = ARGV.shift

transcripts = Hash.new
File.open(fn).each_line do |line|
	next if line.chomp =~ /^#/
	chrom, source, type, start, stop, s1, strand, s2, comments = line.chomp.split(/\t/)
	if comments =~ /gene_id \"(\S+)\";/
		gene_id = $1
	else
		puts "Failed to find gene_id from #{line.chomp}"
		next
	end
	if comments =~ /gene_name \"(\S+)\";/
		gene_name = $1
	else
		puts "Failed to find gene_name from #{line.chomp}"
		next
	end
	transcripts[gene_id] = gene_name
end

transcripts.each_key { |tx|
	puts "#{tx}\t#{transcripts[tx]}"
}


