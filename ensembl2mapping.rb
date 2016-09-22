#!/usr/bin/ruby

# parses Emsembl FASTA file with detailed header into a mapping from transcript->gene_symbol

# >ENST00000448914.1 cdna:known chromosome:GRCh38:14:22449113:22449125:1 gene:ENSG00000228985.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:TRDD3 description:T cell receptor delta diversity 3 [Source:HGNC Symbol;Acc:HGNC:12256]
fn = ARGV.shift

transcripts = Hash.new
gene = ""
File.open(fn).each_line do |line|
	line.chomp!
	next unless line =~ /^>/
	arr = line.chomp.split(/\s+/)
	tx = arr.shift.gsub(/>/, "")
	arr.each { |el|
		key, value = el.split(":")
		next unless key == "gene_symbol"
		gene = value
	}
	transcripts[tx] = gene
end

transcripts.each_key { |tx|
	puts "#{tx}\t#{transcripts[tx]}"
}


