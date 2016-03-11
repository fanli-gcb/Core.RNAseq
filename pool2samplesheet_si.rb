#!/usr/bin/ruby

# input is a data.csv file from the BaseSpace 'Plan Run' tab
# write out a sample sheet
# single-index version

if ARGV.size == 2
	fn = ARGV.shift
	out_fn = ARGV.shift
else
	puts "USAGE: #{$0} fn out_fn"
	exit -1
end
fp = File.new(out_fn, "w")

def prompt(default, *args)
	print(*args)
	result = gets
	result.chomp!
	return result.empty? ? default : result
end

def rc(seq)
	out = seq.reverse.tr("ATGCatgc", "TACGtacg")
	return out
end

fp.puts "[Header]"

out = prompt("4", "IEMFileVersion [default=4]:")
fp.puts "IEMFileVersion,#{out}"
out = prompt("Grace Aldrovandi", "Investigator Name [default=Grace Aldrovandi]:")
fp.puts "Investigator Name,#{out}"
out = prompt("", "Experiment Name [default=none]:")
fp.puts "Experiment Name,#{out}"
out = prompt("", "Date [default=none]:")
fp.puts "Date,#{out}"
out = prompt("GenerateFASTQ", "Workflow [default=GenerateFASTQ]:")
fp.puts "Workflow,#{out}"
out = prompt("NextSeq FASTQ Only", "Application [default=NextSeq FASTQ Only]:")
fp.puts "Application,#{out}"
out = prompt("Nextera XT", "Assay [default=Nextera XT]:")
fp.puts "Assay,#{out}"
out = prompt("", "Description [default=]:")
fp.puts "Description,#{out}"
out = prompt("Amplicon", "Chemistry [default=Amplicon]:")
fp.puts "Chemistry,#{out}"

fp.puts ""
fp.puts "[Reads]"
out = prompt("4", "Reads(1) [default=76]:")
fp.puts "#{out}"
out = prompt("4", "Reads(2) [default=76]:")
fp.puts "#{out}"

fp.puts ""
fp.puts "[Settings]", "ReverseComplement,0"
fp.puts ""
fp.puts "[Data]"
fp.puts "Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_ID,index,Sample_Project,Description"

first = true
File.open(fn).each_line do |line|
	if first
		first = false
		next
	end
	id, project, well, i1_str, i2_str = line.chomp.split(",")
	id = id.gsub("\"", "")
	i7_id, i7_index = i1_str.gsub("\"","").split(" - ")
#	i5_id, i5_index = i2_str.gsub("\"","").split(" - ")
#	i5_index = rc(i5_index)
	fp.puts "#{id},,,,#{i7_id},#{i7_index},,"
end


