#!/usr/bin/ruby

n = 10400
read_len = 25
pool_size = 50

class Array
  def cumulative_sum
    sum = 0
    self.map{|x| sum += x}
  end
end

# sequence pool
pool = Array.new
1.upto(pool_size) do |i|
	pool << (1..read_len).map { ['A','T','G','C'].to_a[rand(4)] }.join
end

# quality probabilities
qarr = [0].concat(Array.new(25, 0.01)).concat(Array.new(15, 0.05))
qarr = qarr.cumulative_sum

1.upto(n) do |i|
	puts "@NS500273:3:H0NGTAGXX:1:11101:18119:#{i} 1:N:0:24"
	puts "#{pool[rand(pool_size)]}"
	puts "+"
	qstr = ""
	1.upto(read_len) do 
		r = rand()
		0.upto(qarr.length-1) do |j|
			if r < qarr[j]
				qstr = "#{qstr}#{(j+33).chr}"
				break
			end
		end
	end
	puts qstr
end

#<A.AAFFFFFF7FFFFFFFAFFFFFFFFFF<FF<AFFAFF<FFFFFAFFF.<FAF7F7FAFFFFFFFAFFF.FFFF
#ASCII 33 to 126
#CTCACACGCGAGTCGCGAGCCTGT
#0?&@)F?3)I$C1@>215IB%BB,C
#TTCTCTAGAACGCATGGGAAGCCTT
