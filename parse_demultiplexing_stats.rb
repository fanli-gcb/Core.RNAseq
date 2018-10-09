#!/usr/bin/ruby

require 'nokogiri'

# input is a DemultiplexingStats.xml file from bcl2fastq
# write out a nice table

fn = ARGV.shift
f = File.open(fn)
doc = Nokogiri::XML(f)

puts "SampleID\tBarcode\tLane\tCount"

counts = Hash.new(0)

root = doc.xpath('//Project')[0]
root.xpath('Sample').each do |thing|
	sample = thing.attr("name")
	barcode = thing.at_xpath('Barcode')
	barcode_str = barcode.attr("name")
	barcode.xpath('Lane').each do |lane|
		lane_number = lane.attr("number")
		count = lane.at_xpath('BarcodeCount').content
		counts[sample] += count.to_i
#		puts "#{sample}\t#{barcode_str}\t#{lane_number}\t#{count}"
	end
end

f.close

counts.each_key { |sample|
	puts "#{sample}\t#{counts[sample]}"
}

