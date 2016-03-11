#!/usr/bin/ruby

# reverse complement

fn = ARGV.shift

def rc(seq)
	out = seq.reverse.tr("ATGCatgc", "TACGtacg")
	return out
end

File.open(fn).each_line do |line|
	puts rc(line.chomp)
end
