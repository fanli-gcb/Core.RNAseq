#!/usr/bin/ruby

# generate the basic/advanced single cell RNA-seq workflow as a shell script
# requires a sample sheet file (see included example SAMPLE_SHEET.txt) and a comparison sheet (COMPARISON.txt)

flag_transcriptomeOnly = false
flag_comparison = false
flag_kallisto = false
flag_cuffmerge = true
flag_qualitytrim = true # changed to true by default
aligner = "hisat2" # use HISAT2 aligner by default, tophat2 as an option
de_pipeline = "kallisto" # "tuxedo" for Cufflinks/Cuffdiff/cummeRbund, "featureCounts" for featureCounts/scde, "kallisto" for kallisto+monocle
read_mismatches = 2
read_gap_length = 2
read_edit_dist = 2
max_multihits = 10
num_jobs = 1 # used to split tophat/QC/cuffetc steps into multiple batches
num_alignment_threads = 6	# used for Tuxedo suite multithreading

full_command = ARGV.join(" ")

if (ARGV.size >= 6)
	while (arg = ARGV.shift)
		if arg == "-s"
			sample_sheet_fn = ARGV.shift
		elsif arg == "-c"
			comparison_sheet_fn = ARGV.shift
			flag_comparison = true
		elsif arg == "-o"
			output_dir = ARGV.shift
			output_dir = output_dir.gsub(/\/$/, "")
		elsif arg == "--output-tophat-dir"
			output_tophat_dir = ARGV.shift
			output_tophat_dir = output_tophat_dir.gsub(/\/$/, "")
		elsif arg == "--output-hisat2-dir"
			output_hisat2_dir = ARGV.shift
			output_hisat2_dir = output_hisat2_dir.gsub(/\/$/, "")
		elsif arg == "--output-QC-dir"
			output_QC_dir = ARGV.shift
			output_QC_dir = output_QC_dir.gsub(/\/$/, "")
		elsif arg == "--output-cufflinks-dir"
			output_cufflinks_dir = ARGV.shift
			output_cufflinks_dir = output_cufflinks_dir.gsub(/\/$/, "")
		elsif arg == "--output-cuffquant-dir"
			output_cuffquant_dir = ARGV.shift
			output_cuffquant_dir = output_cuffquant_dir.gsub(/\/$/, "")
		elsif arg == "--output-featureCounts-dir"
			output_featureCounts_dir = ARGV.shift
			output_featureCounts_dir = output_featureCounts_dir.gsub(/\/$/, "")
		elsif arg == "--output-deseq-dir"
			output_deseq_dir = ARGV.shift
			output_deseq_dir = output_deseq_dir.gsub(/\/$/, "")
		elsif arg == "--output-code-dir"
			output_code_dir = ARGV.shift
			output_code_dir = output_code_dir.gsub(/\/$/, "")
		elsif arg == "--output_cuffmerge_dir"
			output_cuffmerge_dir = ARGV.shift
			output_cuffmerge_dir = output_cuffmerge_dir.gsub(/\/$/, "")
		elsif arg == "--output_kallisto_dir"
			output_kallisto_dir = ARGV.shift
			output_kallisto_dir = output_kallisto_dir.gsub(/\/$/, "")
		elsif arg == "--output-saturation-dir"
			output_saturation_dir = ARGV.shift
			output_saturation_dir = output_saturation_dir.gsub(/\/$/, "")
		elsif arg == "--fastq-dir"
			fastq_dir = ARGV.shift
			fastq_dir = fastq_dir.gsub(/\/$/, "")
		elsif arg == "-g"
			genome=  ARGV.shift
		elsif arg == "-a"
			analysis_type = ARGV.shift
		elsif arg == "-p"
			project = ARGV.shift
		elsif arg == "--transcriptome-only"
			flag_transcriptomeOnly = true
		elsif arg == "--kallisto"
			flag_kallisto = true
		elsif arg == "--skip-cuffmerge"
			flag_cuffmerge = false
		elsif arg == "--num-jobs"
			num_jobs = ARGV.shift.to_i
		elsif arg == "--num-alignment-threads"
			num_alignment_threads = ARGV.shift.to_i
		elsif arg == "--read-mismatches"
			read_mismatches = ARGV.shift.to_i
		elsif arg == "--read-gap-length"
			read_gap_length = ARGV.shift.to_i
		elsif arg == "--read-edit-dist"
			read_edit_dist = ARGV.shift.to_i
		elsif arg == "--max-multihits"
			max_multihits = ARGV.shift.to_i
		elsif arg == "--de-pipeline"
			de_pipeline = ARGV.shift
		elsif arg == "--genes-gtf"
			genes_gtf_file = ARGV.shift
			genes_gtf_index = genes_gtf_file.gsub("\.gtf","")
		elsif arg == "--kallisto-index"
			kallisto_index_file = ARGV.shift
		elsif arg == "--transcript-to-gene-file"
			transcript_to_gene_file = ARGV.shift
		elsif arg == "--dont-quality-trim"
			flag_qualitytrim = false
		elsif arg == "--aligner"
			aligner = ARGV.shift
		else
			puts "ERROR: Unrecognized option #{arg}"
			exit -1
		end
	end	
else
	puts "USAGE: #{$0} -s sample_sheet_file -c comparison_file -g GRCh38|hg19|mm10 --output-code-dir output_code_dir -o output_dir -a analysis_type -p project [--transcriptome-only] [--num-jobs num_jobs] [--num-alignment-threads num_alignment_threads] [--de-pipeline de_pipeline] [--aligner hisat2|tophat] [--genes-gtf genes_gtf_file]"
	exit -1
end

# error checking
if sample_sheet_fn.nil?
	puts "ERROR: sample_sheet_file must be specified"
	exit -1
end
##if comparison_sheet_fn.nil?
##	puts "ERROR: comparison_file must be specified"
##	exit -1
##end
if output_dir.nil?
	puts "ERROR: output_dir must be specified"
	exit -1
end
# default locations are:
#		#{output_dir}/tophat
#		#{output_dir}/cufflinks
#		#{output_dir}/QC
#		#{output_dir}/cuffquant
#		#{output_dir}/cuffmerge
#		#{output_dir}/featureCounts
#		#{output_dir}/DESeq2
if output_tophat_dir.nil?
	output_tophat_dir = "#{output_dir}/tophat"
end
if output_hisat2_dir.nil?
	output_hisat2_dir = "#{output_dir}/hisat2"
end
if output_QC_dir.nil?
	output_QC_dir = "#{output_dir}/QC"
end
if output_cufflinks_dir.nil?
	output_cufflinks_dir = "#{output_dir}/cufflinks"
end
if output_cuffquant_dir.nil?
	output_cuffquant_dir = "#{output_dir}/cuffquant"
end
if output_cuffmerge_dir.nil?
	output_cuffmerge_dir = "#{output_dir}/cuffmerge"
end
if output_featureCounts_dir.nil?
	output_featureCounts_dir = "#{output_dir}/featureCounts"
end
if output_deseq_dir.nil?
	output_deseq_dir = "#{output_dir}/DESeq2"
end
if output_kallisto_dir.nil?
	output_kallisto_dir = "#{output_dir}/kallisto"
end
if output_saturation_dir.nil?
	output_saturation_dir = "#{output_dir}/saturation"
end
if fastq_dir.nil?
	fastq_dir = output_dir
end
if output_code_dir.nil?
	puts "ERROR: output_code_dir must be specified"
	exit -1
end
if analysis_type.nil? || (!(analysis_type == "basic" || analysis_type == "advanced"))
	puts "analysis type must be \"basic\" or \"advanced\""
	exit -1
end
if project.nil?
	puts "ERROR: project must be specified"
	exit -1
end
if genome.nil?
	puts "ERROR: genome must be specified"
	exit -1
else
	if !(genome == "hg19" || genome == "mm10" || genome == "GRCh38")
		puts "ERROR: genome must be 'hg19', 'GRCh38', or 'mm10'"
		exit -1
	end
end
if !(de_pipeline == "tuxedo" || de_pipeline == "featureCounts" || de_pipeline == "kallisto")
	puts "ERROR: de-pipeline must be 'featureCounts', 'kallisto', or 'tuxedo'"
	exit -1
end
if (aligner == "tophat")
	output_alignment_dir = output_tophat_dir
elsif (aligner == "hisat2") 
	output_alignment_dir = output_hisat2_dir
else
	puts "ERROR: aligner must be 'hisat2' or 'tophat'"
	exit -1
end

if !File.exists?(output_code_dir)
	system("mkdir #{output_code_dir}")
end
system("cp /Lab_Share/fanli/code/Core.RNAseq/* #{output_code_dir}")

out_fn = "#{output_code_dir}/workflow.#{project}.sh"
out_fp = File.new(out_fn, "w")
str_transcriptomeOnly = (flag_transcriptomeOnly) ? "--transcriptome-only" : ""

if genes_gtf_file.nil?
	genes_gtf_file = "/Lab_Share/iGenomes/#{genome}/Annotation/Genes/genes.gtf"
	genes_gtf_index = "/Lab_Share/iGenomes/#{genome}/Annotation/Genes/genes"
end

if (de_pipeline == "kallisto")
	flag_kallisto = true
end
if ((kallisto_index_file.nil? || transcript_to_gene_file.nil?) && flag_kallisto)
	puts "ERROR: --kallisto-index and --transcript-to-gene-file must be specified if --kallisto is flagged!"
	exit -1
end
if (!flag_cuffmerge) && output_cuffmerge_dir.nil?
	puts "ERROR: --output_cuffmerge_dir must be specified if --skip-cuffmerge is selected"
	exit -1
end

# read in sample sheet
samples = Array.new
sample_header = Array.new
first = true
File.open(sample_sheet_fn).each_line do |line|
	if first
		sample_header = line.chomp.split(/\t/)
		first = false
	else
		this = Hash.new
		arr = line.chomp.split(/\t/)
		0.upto(sample_header.length-1) do |i|
			this[sample_header[i]] = arr[i]
		end
		samples << this
	end
end
# read in comparison sheet
if flag_comparison
	comparisons = Array.new
	comparison_header = Array.new
	first = true
	File.open(comparison_sheet_fn).each_line do |line|
		if first
			comparison_header = line.chomp.split(/\t/)
			first = false
		else
			this = Hash.new
			arr = line.chomp.split(/\t/)
			0.upto(comparison_header.length-1) do |i|
				this[comparison_header[i]] = arr[i]
			end
			comparisons << this
		end
	end
end

#### Start writing workflow script
sub_fps = Array.new

out_fp.puts "### #{analysis_type} RNA-seq analysis for #{project} ###"
out_fp.puts "###"
out_fp.puts "#\truby generate_workflow.rb #{full_command}"
out_fp.puts "export PICARDPATH=/usr/local/picard-tools-1.118"
out_fp.puts "export GSEAPATH=/usr/local/GSEA"
out_fp.puts "################################################"
if flag_qualitytrim
	out_fp.puts "###\t0. Quality trimming"
	out_fp.puts "#\tINPUT:"
	samples.each { |sample|
		out = sample["fastq"].split(",")
		out.each { |fastq|
			out_fp.puts "#\t\t#{fastq}"
		}
	}
	out_fp.puts "#\tEXECUTION:"
	out_fp.puts "mkdir -p #{output_QC_dir}"
	sub_fps.clear
	1.upto(num_jobs) do |i|
		sub_fps << File.new("#{output_code_dir}/workflow.#{project}.trim.#{i}.sh", "w")
		out_fp.puts "bash workflow.#{project}.trim.#{i}.sh &> workflow.#{project}.trim.#{i}.log &"
	end
	out_fp.puts "wait"
	i = 0
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
		fastq_str = sample["fastq"].split(",").collect { |x| "#{fastq_dir}/#{x}" }.join(" ")
		out_files = sample["fastq"].gsub(".fastq.gz", "").gsub(".fastq", "") # only the basename
		out_str = out_files.split(",").collect { |x| "#{fastq_dir}/#{x}.trimmed.fastq.gz #{fastq_dir}/#{x}.unpaired.fastq.gz" }.join(" ")
#		cmd = "java -jar /mnt/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads #{num_alignment_threads} -phred33 #{fastq_str} #{out_str} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>> #{fastq_dir}/trim.log"
		cmd = "/mnt/trim_galore/trim_galore --fastqc_args \"--outdir=#{fastq_dir}\" --stringency 1 --output_dir #{fastq_dir} --paired #{fastq_str}"
		sub_fps[(i % num_jobs)].puts(cmd)
		i += 1
		j = 1
		sample["fastq"].split(",").each { |fastq_str|
			sample_name = fastq_str.gsub(".fastq.gz", "").gsub(".fastq","")
			if (File.extname(fastq_str)==".gz")
				extstr = "fq.gz"
				extstr2 = "fastq.gz"
			else
				extstr = "fq"
				extstr2 = "fastq"
			end
			out_fp.puts "mv #{fastq_dir}/#{sample_name}_val_#{j}.#{extstr} #{fastq_dir}/#{sample_name}.trimmed.#{extstr2}"
			out_fp.puts "mv #{fastq_dir}/#{sample_name}.#{extstr2}_trimming_report.txt #{output_QC_dir}/#{sample_name}.#{extstr2}_trimming_report.txt"
			out_fp.puts "mv #{fastq_dir}/#{sample_name}_val_#{j}_fastqc.zip #{output_QC_dir}/#{sample_name}_val_#{j}_fastqc.zip"
			out_fp.puts "mv #{fastq_dir}/#{sample_name}_val_#{j}_fastqc.html #{output_QC_dir}/#{sample_name}_val_#{j}_fastqc.html"
			j += 1
		}
		sample["fastq"] = sample["fastq"].gsub(".fastq", ".trimmed.fastq")
	}
	sub_fps.each { |fp|
		fp.close
	}
end

out_fp.puts "################################################"
out_fp.puts "###\t1. Alignment and QC metrics"
out_fp.puts "###\ta. Run #{aligner}"
out_fp.puts "#\tINPUT:"
samples.each { |sample|
	out = sample["fastq"].split(",")
	out.each { |fastq|
		out_fp.puts "#\t\t#{fastq}"
	}
}
out_fp.puts "#\tREQUIRED DATA FILES"
out_fp.puts "#\t\tHISAT2/tophat index files"
out_fp.puts "#\t\t#{genes_gtf_file}"
out_fp.puts "#\tEXECUTION:"
sub_fps.clear
out_fp.puts "mkdir -p #{output_alignment_dir}"
1.upto(num_jobs) do |i|
	sub_fps << File.new("#{output_code_dir}/workflow.#{project}.align.#{i}.sh", "w")
	out_fp.puts "bash workflow.#{project}.align.#{i}.sh &> workflow.#{project}.align.#{i}.log &"
end
out_fp.puts "wait"
i = 0
samples.each { |sample|
	sample_id = sample["sample_id"]
	libtype = sample["library-type"]
	if aligner == "tophat"
		fastq_str = sample["fastq"].split(",").collect { |x| "#{fastq_dir}/#{x}" }.join(" ")
		cmd = "tophat --read-mismatches #{read_mismatches} --read-gap-length #{read_gap_length} --read-edit-dist #{read_edit_dist} --max-multihits #{max_multihits} --library-type #{libtype} --GTF #{genes_gtf_file} --transcriptome-index #{genes_gtf_index} #{str_transcriptomeOnly} --num-threads #{num_alignment_threads} --output-dir #{output_alignment_dir}/#{sample_id} /Lab_Share/iGenomes/#{genome}/Sequence/Bowtie2Index/genome #{fastq_str}"
	elsif aligner == "hisat2"
		sub_fps[(i % num_jobs)].puts("mkdir -p #{output_alignment_dir}/#{sample_id}")
		arr = sample["fastq"].split(",")
		if arr.length == 1
			fastq_str = "-U #{fastq_dir}/#{arr}"
		elsif arr.length == 2
			fastq_str = "-1 #{fastq_dir}/#{arr[0]} -2 #{fastq_dir}/#{arr[1]}"
		else
			puts "ERROR: invalid fastq string!"
			exit -1
		end
		cmd = "hisat2 --threads #{num_alignment_threads} --no-unal -x /Lab_Share/iGenomes/#{genome}/Sequence/HISAT2Index/genome_tran #{fastq_str} 2> #{output_alignment_dir}/#{sample_id}/align_summary.txt | samtools view -bS -o #{output_alignment_dir}/#{sample_id}/accepted_hits.bam -"
	end
	sub_fps[(i % num_jobs)].puts(cmd)
	i += 1
}
sub_fps.each { |fp|
	fp.close
}

out_fp.puts "#\tOUTPUT:"
samples.each { |sample|
  sample_id = sample["sample_id"]
	out_fp.puts "#\t\t#{output_alignment_dir}/#{sample_id}"
}
out_fp.puts "#\tCHECKPOINT:"
samples.each { |sample|
  sample_id = sample["sample_id"]
	out_fp.puts "[ -f \"#{output_alignment_dir}/#{sample_id}/accepted_hits.bam\" ] || (echo \"ERROR: #{output_alignment_dir}/#{sample_id}/accepted_hits.bam not found!\" && exit)"
}

out_fp.puts "", "###\tb. Library-level QC metrics"
out_fp.puts "#\tINPUT:"
samples.each { |sample|
  out = sample["fastq"].split(",")
  out.each { |fastq|
    out_fp.puts "#\t\t#{fastq}"
  }
}
out_fp.puts "#\tEXECUTION:"
out_fp.puts "mkdir -p #{output_QC_dir}"
sub_fps.clear
1.upto(num_jobs) do |i|
	sub_fps << File.new("#{output_code_dir}/workflow.#{project}.QC.#{i}.sh", "w")
	out_fp.puts "bash workflow.#{project}.QC.#{i}.sh &> workflow.#{project}.QC.#{i}.log &"
end
out_fp.puts "wait"
i = 0
samples.each { |sample|
	fastqs = sample["fastq"].split(",").collect { |x| "#{fastq_dir}/#{x}" }
	sample_id = sample["sample_id"]
	out = "#{output_QC_dir}/#{sample_id}.unaligned.bam"
	if (fastqs.length==1)
		cmd = "java -Xmx8g -Djava.io.tmpdir=/mnt/tmp -jar \$PICARDPATH/FastqToSam.jar FASTQ=#{fastqs[0]} OUTPUT=#{out} SAMPLE_NAME=\"#{sample_id}\""
	else
		cmd = "java -Xmx8g -Djava.io.tmpdir=/mnt/tmp -jar \$PICARDPATH/FastqToSam.jar FASTQ=#{fastqs[0]} FASTQ2=#{fastqs[1]} OUTPUT=#{out} SAMPLE_NAME=\"#{sample_id}\""
	end
	sub_fps[(i % num_jobs)].puts(cmd)
	sub_fps[(i % num_jobs)].puts "java -Xmx8g -Djava.io.tmpdir=/mnt/tmp -jar \$PICARDPATH/QualityScoreDistribution.jar INPUT=#{out} OUTPUT=#{output_QC_dir}/#{sample_id}.QualityScoreDistribution.txt CHART_OUTPUT=#{output_QC_dir}/#{sample_id}.QualityScoreDistribution.pdf"
	sub_fps[(i % num_jobs)].puts "java -Xmx8g -Djava.io.tmpdir=/mnt/tmp -jar \$PICARDPATH/MeanQualityByCycle.jar INPUT=#{out} OUTPUT=#{output_QC_dir}/#{sample_id}.MeanQualityByCycle.txt CHART_OUTPUT=#{output_QC_dir}/#{sample_id}.MeanQualityByCycle.pdf"
	i += 1
}

out_fp.puts "#\tOUTPUT:"
samples.each { |sample|
	sample_id = sample["sample_id"]
	out_fp.puts "#\t\t#{sample_id}.unaligned.bam"
	out_fp.puts "#\t\t#{sample_id}.QualityScoreDistribution.txt"
	out_fp.puts "#\t\t#{sample_id}.QualityScoreDistribution.pdf"
	out_fp.puts "#\t\t#{sample_id}.MeanQualityByCycle.txt"
	out_fp.puts "#\t\t#{sample_id}.MeanQualityByCycle.pdf"
}

i = 0
out_fp.puts "", "###\tc. Insert size metrics"
out_fp.puts "#\tINPUT:"
samples.each { |sample|
	sample_id = sample["sample_id"]
	out_fp.puts "#\t\t#{sample_id}/accepted_hits.bam"
}
out_fp.puts "#\tEXECUTION:"
samples.each { |sample|
	sample_id = sample["sample_id"]
	sub_fps[(i % num_jobs)].puts "java -Xmx8g -Djava.io.tmpdir=/mnt/tmp -jar \$PICARDPATH/CollectInsertSizeMetrics.jar INPUT=#{output_alignment_dir}/#{sample_id}/accepted_hits.bam OUTPUT=#{output_QC_dir}/#{sample_id}.InsertSizeMetrics.txt HISTOGRAM_FILE=#{output_QC_dir}/#{sample_id}.InsertSizeMetrics.pdf"
	sub_fps[(i % num_jobs)].puts "java -Xmx8g -Djava.io.tmpdir=/mnt/tmp -jar \$PICARDPATH/EstimateLibraryComplexity.jar INPUT=#{output_alignment_dir}/#{sample_id}/accepted_hits.bam OUTPUT=#{output_QC_dir}/#{sample_id}.EstimatedLibraryComplexity.txt"
	i += 1
}
out_fp.puts "#\tOUTPUT:"
samples.each { |sample|
	sample_id = sample["sample_id"]
	out_fp.puts "#\t\t#{sample_id}.InsertSizeMetrics.txt"
	out_fp.puts "#\t\t#{sample_id}.InsertSizeMetrics.pdf"
	out_fp.puts "#\t\t#{sample_id}.EstimatedLibraryComplexity.txt"
}

i = 0
out_fp.puts "", "###\td. RNA-seq metrics"
out_fp.puts "#\tINPUT:"
samples.each { |sample|
	sample_id = sample["sample_id"]
	out_fp.puts "#\t\t#{sample_id}/accepted_hits.bam"
}
out_fp.puts "#\tEXECUTION:"
samples.each { |sample|
	sample_id = sample["sample_id"]
	libtype = sample["library-type"]
	if libtype == "fr-unstranded"
		strand_specificity="NONE"
	elsif libtype == "fr-firststrand"
		strand_specificity="FIRST_READ_TRANSCRIPTION_STRAND"
	elsif libtype == "fr-secondstrand"
		strand_specificity="SECOND_READ_TRANSCRIPTION_STRAND"
	else
		out_fp.puts "ERROR: Invalid library-type #{libtype}"
		exit -1
	end
	sub_fps[(i % num_jobs)].puts "java -Xmx8g -Djava.io.tmpdir=/mnt/tmp -jar \$PICARDPATH/CollectRnaSeqMetrics.jar INPUT=#{output_alignment_dir}/#{sample_id}/accepted_hits.bam OUTPUT=#{output_QC_dir}/#{sample_id}.RnaSeqMetrics.txt CHART_OUTPUT=#{output_QC_dir}/#{sample_id}.RnaSeqMetrics.pdf REF_FLAT=/Lab_Share/iGenomes/#{genome}/Annotation/Genes/refFlat.txt STRAND_SPECIFICITY=#{strand_specificity} RIBOSOMAL_INTERVALS=/Lab_Share/iGenomes/#{genome}/Annotation/Genes/rmsk_rRNA.interval_list"
	sub_fps[(i % num_jobs)].puts "./plot_RnaSeqMetrics.R #{output_QC_dir}/#{sample_id}.RnaSeqMetrics.txt #{output_QC_dir}/#{sample_id}.RnaSeqMetrics.class.pdf"
	i += 1
}
sub_fps.each { |fp|
	fp.close
}
out_fp.puts "#\tOUTPUT:"
samples.each { |sample|
	sample_id = sample["sample_id"]
	out_fp.puts "#\t\t#{sample_id}.RnaSeqMetrics.txt"
	out_fp.puts "#\t\t#{sample_id}.RnaSeqMetrics.pdf"
	out_fp.puts "#\t\t#{sample_id}.RnaSeqMetrics.class.pdf"
}

out_fp.puts "", "###\te. Summarize various mapping stats"
out_fp.puts "#\tINPUT:"
out_fp.puts "#\t\talignment output directories", "#\t\t*.RnaSeqMetrics.txt"
out_fp.puts "#\tEXECUTION:"
out_fp.puts "./summ_mapping_stats.R #{output_alignment_dir} #{output_QC_dir}/mapping_stats.#{project}.txt #{output_QC_dir}/mapping_stats.#{project}.pdf"
out_fp.puts "./summ_rnaseq_stats.R #{output_QC_dir} #{output_QC_dir}/rnaseq_stats.#{project}.txt #{output_QC_dir}/rnaseq_stats.#{project}.pdf"
out_fp.puts "./summ_duplication_stats.R #{output_QC_dir} #{output_QC_dir}/duplication_stats.#{project}.txt #{output_QC_dir}/duplication_stats.#{project}.pdf"
out_fp.puts "#\tOUTPUT:"
out_fp.puts "#\t\tmapping_stats.#{project}.txt"
out_fp.puts "#\t\tmapping_stats.#{project}.pdf"
out_fp.puts "#\t\trnaseq_stats.#{project}.txt"
out_fp.puts "#\t\trnaseq_stats.#{project}.pdf"

if flag_kallisto
	out_fp.puts "", "###\tf. Run kallisto"
	out_fp.puts "#\tINPUT:"
	samples.each { |sample|
		out = sample["fastq"].split(",")
		out.each { |fastq|
			out_fp.puts "#\t\t#{fastq}"
		}
	}
	out_fp.puts "#\tREQUIRED DATA FILES"
	out_fp.puts "#\t\t#{kallisto_index_file}"
	out_fp.puts "#\tEXECUTION:"
	sub_fps.clear
	out_fp.puts "mkdir -p #{output_kallisto_dir}"
	out_fp.puts "mkdir -p #{output_deseq_dir}"
	1.upto(num_jobs) do |i|
		sub_fps << File.new("#{output_code_dir}/workflow.#{project}.kallisto.#{i}.sh", "w")
		out_fp.puts "bash workflow.#{project}.kallisto.#{i}.sh &> workflow.#{project}.kallisto.#{i}.log &"
	end
	out_fp.puts "wait"
	i = 0
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
		fastq_str = sample["fastq"].split(",").collect { |x| "#{fastq_dir}/#{x}" }.join(" ")
		cmd = "kallisto quant -i #{kallisto_index_file} -o #{output_kallisto_dir}/#{sample_id} -b 100 -t #{num_alignment_threads} #{fastq_str}"
		sub_fps[(i % num_jobs)].puts(cmd)
		i += 1
	}
	sub_fps.each { |fp|
		fp.close
	}
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "ruby kallisto_to_gene.rb #{output_kallisto_dir}/#{sample_id}/abundance.tsv #{transcript_to_gene_file} > #{output_kallisto_dir}/#{sample_id}/abundance.gene.tsv"
	}
	out_fp.puts "ls #{output_kallisto_dir} -1 | awk 'BEGIN {ORS=\"\\t\";} {print;}' | sed \"s/\\t$/\\n/\" > #{output_deseq_dir}/read_counts.txt"
	out_fp.puts "paste #{output_kallisto_dir}/*/abundance.gene.tsv | awk 'BEGIN {OFS=\"\\t\"} {printf \"%s\\t\",$1; for(i=2;i<=NF;i+=2) {printf \"%s\\t\",$i}; print \"\"}' | sed \"s/\\t\\t$/\\n/\" >> #{output_deseq_dir}/read_counts.txt"
	out_fp.puts "#\tOUTPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{output_kallisto_dir}/#{sample_id}"
	}
	out_fp.puts "#\tCHECKPOINT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "[ -f \"#{output_kallisto_dir}/#{sample_id}/abundance.tsv\" ] || (echo \"ERROR: #{output_kallisto_dir}/#{sample_id}/abundance.tsv not found!\" && exit)"
	}
end


if de_pipeline == "tuxedo"
	out_fp.puts "","","################################################"
	out_fp.puts "###\t2. Differential expression analysis using Tuxedo suite"
	out_fp.puts "###\ta. Cufflinks"
	out_fp.puts "#\tINPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{sample_id}/accepted_hits.bam"
	}
	out_fp.puts "#\tREQUIRED DATA FILES:"
	out_fp.puts "#\t\t#{genes_gtf_file}"
	out_fp.puts "#\t\trmsk_rRNA.gtf"
	out_fp.puts "#\tEXECUTION:"
	sub_fps.clear
	out_fp.puts "mkdir -p #{output_cufflinks_dir}"
	1.upto(num_jobs) do |i|
		sub_fps << File.new("#{output_code_dir}/workflow.#{project}.cufflinks.#{i}.sh", "w")
		out_fp.puts "bash workflow.#{project}.cufflinks.#{i}.sh &> workflow.#{project}.cufflinks.#{i}.log &"
	end
	out_fp.puts "wait"
	i = 0
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
		cmd = "cufflinks --output-dir #{output_cufflinks_dir}/#{sample_id} --num-threads #{num_alignment_threads} --GTF #{genes_gtf_file} --mask-file /Lab_Share/iGenomes/#{genome}/Annotation/Genes/rmsk_rRNA.gtf --multi-read-correct --library-type #{libtype} --upper-quartile-norm --compatible-hits-norm --quiet --no-update-check #{output_alignment_dir}/#{sample_id}/accepted_hits.bam"
		sub_fps[(i % num_jobs)].puts cmd
		i += 1
	}
	sub_fps.each { |fp|
		fp.close
	}
	out_fp.puts "#\tOUTPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{output_cufflinks_dir}/#{sample_id}"
	}
	out_fp.puts "#\tCHECKPOINT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "[ -f \"#{output_cufflinks_dir}/#{sample_id}/genes.fpkm_tracking\" ] || (echo \"ERROR: #{output_cufflinks_dir}/#{sample_id}/genes.fpkm_tracking not found!\" && exit)"
	}
	
	if flag_cuffmerge
		out_fp.puts "", "###\tb. Cuffmerge"
		out_fp.puts "#\tINPUT:"
		samples.each { |sample|
			sample_id = sample["sample_id"]
			out_fp.puts "#\t\t#{sample_id}/transcripts.gtf"
		}
		out_fp.puts "#\tREQUIRED DATA FILES:"
		out_fp.puts "#\t\t#{genes_gtf_file}"
		out_fp.puts "#\tEXECUTION:"
		out_fp.puts "mkdir -p #{output_cufflinks_dir}"
		samples.each { |sample|
			sample_id = sample["sample_id"]
			out_fp.puts "echo \"#{output_cufflinks_dir}/#{sample_id}/transcripts.gtf\" >> #{output_cuffmerge_dir}/assembly_list.txt"
		}
		out_fp.puts "cuffmerge -o #{output_cuffmerge_dir} --ref-gtf #{genes_gtf_file} -p #{num_alignment_threads} --ref-sequence /Lab_Share/iGenomes/#{genome}/Sequence/Chromosomes #{output_cuffmerge_dir}/assembly_list.txt"
		out_fp.puts "#\tOUTPUT:"
		out_fp.puts "#\t\t#{output_cuffmerge_dir}/merged.gtf"
		out_fp.puts "#\tCHECKPOINT:"
		out_fp.puts "[ -f \"#{output_cuffmerge_dir}/merged.gtf\" ] || (echo \"ERROR: #{output_cuffmerge_dir}/merged.gtf not found!\" && exit)"
	else
		out_fp.puts "", "###\tb. Cuffmerge"
		out_fp.puts "### SKIPPED DUE TO --skip-cuffmerge FLAG ###"
	end
	
	out_fp.puts "", "###\tc. Cuffquant"
	out_fp.puts "#\tINPUT:"
	out_fp.puts "#\t\t#{output_cuffmerge_dir}/merged.gtf"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{sample_id}/accepted_hits.bam"
	}
	out_fp.puts "#\tEXECUTION:"
	sub_fps.clear
	out_fp.puts "mkdir -p #{output_cuffquant_dir}"
	1.upto(num_jobs) do |i|
		sub_fps << File.new("#{output_code_dir}/workflow.#{project}.cuffquant.#{i}.sh", "w")
		out_fp.puts "bash workflow.#{project}.cuffquant.#{i}.sh &> workflow.#{project}.cuffquant.#{i}.log &"
	end
	out_fp.puts "wait"
	i = 0
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
		sub_fps[(i % num_jobs)].puts "cuffquant --output-dir #{output_cuffquant_dir}/#{sample_id} -p #{num_alignment_threads} --mask-file /Lab_Share/iGenomes/#{genome}/Annotation/Genes/rmsk_rRNA.gtf --multi-read-correct --library-type #{libtype} --quiet --no-update-check #{output_cuffmerge_dir}/merged.gtf #{output_alignment_dir}/#{sample_id}/accepted_hits.bam"
		i += 1
	}
	sub_fps.each { |fp|
		fp.close
	}
	out_fp.puts "#\tOUTPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{sample_id}/abundances.cxb"
	}
	out_fp.puts "#\tCHECKPOINT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "[ -f \"#{output_cuffquant_dir}/#{sample_id}/abundances.cxb\" ] || (echo \"ERROR: #{output_cuffquant_dir}/#{sample_id}/abundances.cxb not found!\" && exit)"
	}

	if flag_comparison
		out_fp.puts "", "###\td. Cuffdiff"
		out_fp.puts "#\tINPUT:"
		out_fp.puts "#\t\t#{output_cuffmerge_dir}/merged.gtf"
		samples.each { |sample|
			sample_id = sample["sample_id"]
			out_fp.puts "#\t\t#{sample_id}/abundances.cxb"
		}
		out_fp.puts "#\tEXECUTION:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			cxb_condA = Array.new
			cxb_condB = Array.new
			samples.each { |sample|
				sample_id = sample["sample_id"]
				if (sample[characteristic] == condA)
					cxb_condA << "#{output_cuffquant_dir}/#{sample_id}/abundances.cxb"
				elsif (sample[characteristic] == condB)
					cxb_condB << "#{output_cuffquant_dir}/#{sample_id}/abundances.cxb"
				end
			}
			outA = cxb_condA.join(",")
			outB = cxb_condB.join(",")
			out_fp.puts "cuffdiff --output-dir #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffdiff_output --labels #{condA},#{condB} -p #{num_alignment_threads} --compatible-hits-norm --multi-read-correct --library-norm-method geometric --dispersion-method pooled --quiet --no-update-check #{output_cuffmerge_dir}/merged.gtf #{outA} #{outB}"
		}
		out_fp.puts "#\tOUTPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.cuffdiff_output/"
		}
		out_fp.puts "#\tCHECKPOINT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "[ -f \"#{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffdiff_output/gene_exp.diff\" ] || (echo \"ERROR: #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffdiff_output/gene_exp.diff not found!\" && exit)"
		}

		out_fp.puts "", "###\te. cummeRbund"
		out_fp.puts "#\tINPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.cuffdiff_output/"
		}
		out_fp.puts "#\tREQUIRED DATA FILES:"
		out_fp.puts "#\t\t#{output_cuffmerge_dir}/merged.gtf"
		out_fp.puts "#\tEXECUTION:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "./run_cummeRbund.R #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffdiff_output #{condA},#{condB} #{output_cuffmerge_dir}/merged.gtf #{genome} #{output_dir}/#{characteristic}-#{condA}-#{condB}.cummeRbund.pdf"
		}
		out_fp.puts "#\tOUTPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.cummeRbund.pdf"
		}
		out_fp.puts "#\tCHECKPOINT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "[ -s \"#{output_dir}/#{characteristic}-#{condA}-#{condB}.cummeRbund.pdf\" ] || (echo \"ERROR: #{output_dir}/#{characteristic}-#{condA}-#{condB}.cummeRbund.pdf is empty!\" && exit)"
		}

		out_fp.puts "", "###\tf. GSEA"
		out_fp.puts "#\tINPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.cuffdiff_output/gene_exp.diff"
		}
		out_fp.puts "#\tEXECUTION:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "mkdir #{output_dir}/#{characteristic}-#{condA}-#{condB}.GSEA_output/"
			out_fp.puts "if [ -d \"#{output_dir}/tmp.GSEA_output/\" ]; then", "\trmdir #{output_dir}/tmp.GSEA_output/" ,"fi"
			out_fp.puts "mkdir #{output_dir}/tmp.GSEA_output/"
			out_fp.puts "awk 'BEGIN {OFS=\"\\t\";FS=\"\\t\";} {if ($14==\"yes\") {print $3,log(($8+0.1)/($9+0.1))/log(2)}}' #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffdiff_output/gene_exp.diff > #{output_dir}/tmp.GSEA_output/cuffdiff_hits.rnk"
			out_fp.puts "java -cp $GSEAPATH/gsea2-2.1.0.jar -Xmx512m xtools.gsea.GseaPreranked -rnk #{output_dir}/tmp.GSEA_output/cuffdiff_hits.rnk -gmx $GSEAPATH/msigdb.v4.0.symbols.gmt -scoring_scheme classic -collapse false -set_min 15 -out #{output_dir}/tmp.GSEA_output -rpt_label \"tmp\" -zip_report true"
			out_fp.puts "mv #{output_dir}/tmp.GSEA_output/* #{output_dir}/#{characteristic}-#{condA}-#{condB}.GSEA_output/"
			out_fp.puts "mv #{output_dir}/#{characteristic}-#{condA}-#{condB}.GSEA_output/tmp.GseaPreranked.* #{output_dir}/#{characteristic}-#{condA}-#{condB}.GSEA_output/#{characteristic}-#{condA}-#{condB}.GseaPreranked/"
			out_fp.puts "mv #{output_dir}/#{characteristic}-#{condA}-#{condB}.GSEA_output/#{characteristic}-#{condA}-#{condB}.GseaPreranked/*.zip #{output_dir}/#{characteristic}-#{condA}-#{condB}.GSEA_output/report.zip"
			out_fp.puts "rmdir #{output_dir}/tmp.GSEA_output/"
		}
		out_fp.puts "#\tOUTPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.GSEA_output/report.zip"
		}

		# Monocle
		if analysis_type == "advanced" && 0
			out_fp.puts "","","################################################"
			out_fp.puts "###\t3. (Advanced) Single-cell analysis using Monocle"
			out_fp.puts "###\ta. Cuffnorm to get expression matrix"
			out_fp.puts "#\tINPUT:"
			out_fp.puts "#\t\t#{output_cuffmerge_dir}/merged.gtf"
			samples.each { |sample|
				sample_id = sample["sample_id"]
				out_fp.puts "#\t\t#{sample_id}/abundances.cxb"
			}
			out_fp.puts "#\tEXECUTION:"
			comparisons.each { |comparison|
				characteristic = comparison["characteristic"]
				condA = comparison["condition_A"]
				condB = comparison["condition_B"]
				cxb_condA = Array.new
				cxb_condB = Array.new
				samples.each { |sample|
					sample_id = sample["sample_id"]
					if (sample[characteristic] == condA)
						cxb_condA << "#{output_cuffquant_dir}/#{sample_id}/abundances.cxb"
					elsif (sample[characteristic] == condB)
						cxb_condB << "#{output_cuffquant_dir}/#{sample_id}/abundances.cxb"
					end
				}
				outA = cxb_condA.join(",")
				outB = cxb_condB.join(",")
				out_fp.puts "cuffnorm --output-dir #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffnorm_output --labels #{condA},#{condB} -p #{num_alignment_threads} --compatible-hits-norm --library-norm-method geometric --quiet --no-update-check #{output_cuffmerge_dir}/merged.gtf #{outA} #{outB}"
			}
			out_fp.puts "#\tOUTPUT:"
			comparisons.each { |comparison|
				characteristic = comparison["characteristic"]
				condA = comparison["condition_A"]
				condB = comparison["condition_B"]
				out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/"
			}
			out_fp.puts "#\tCHECKPOINT:"
			comparisons.each { |comparison|
				characteristic = comparison["characteristic"]
				condA = comparison["condition_A"]
				condB = comparison["condition_B"]
				out_fp.puts "[ -f \"#{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffnorm_output/genes.fpkm_table\" ] || (echo \"ERROR: #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffnorm_output/genes.fpkm_table not found!\" && exit)"
			}
	
			out_fp.puts "", "###\tb. Generate sample sheets for Monocle"
			out_fp.puts "#\tINPUT:"
			out_fp.puts "#\t\tCOMPARISON.txt"
			out_fp.puts "#\tEXECUTION:"
			comparisons.each { |comparison|
				characteristic = comparison["characteristic"]
				condA = comparison["condition_A"]
				condB = comparison["condition_B"]
				out_fp.puts "echo \"sample_name\tgroup\" > #{output_dir}/#{characteristic}-#{condA}-#{condB}.sample_sheet.txt"
				repA = 0
				repB = 0
				samples.each { |sample|
					sample_id = sample["sample_id"]
					if (sample[characteristic] == condA)
						out_fp.puts "echo \"#{condA}_#{repA}\t#{condA}\" >> #{output_dir}/#{characteristic}-#{condA}-#{condB}.sample_sheet.txt"
						repA += 1
					elsif (sample[characteristic] == condB)
						out_fp.puts "echo \"#{condB}_#{repB}\t#{condB}\" >> #{output_dir}/#{characteristic}-#{condA}-#{condB}.sample_sheet.txt"
						repB += 1
					end
				}
				out_fp.puts "./reorder_cuffnorm_table.R #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffnorm_output/genes.fpkm_table #{output_dir}/#{characteristic}-#{condA}-#{condB}.sample_sheet.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffnorm_output/genes.fpkm_table.reordered"
			}
			out_fp.puts "#\tOUTPUT:"
			comparisons.each { |comparison|
				characteristic = comparison["characteristic"]
				condA = comparison["condition_A"]
				condB = comparison["condition_B"]
				out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.sample_sheet.txt"
				out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.cuffnorm_output/genes.fpkm_table.reordered"
			}
	
			out_fp.puts "", "###\tc. Run Monocle"
			out_fp.puts "#\tINPUT:"
			comparisons.each { |comparison|
				characteristic = comparison["characteristic"]
				condA = comparison["condition_A"]
				condB = comparison["condition_B"]
				out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.cuffnorm_output/genes.fpkm_table.reordered"
				out_fp.puts "#\t\t#{output_dir}/#{characteristic}-#{condA}-#{condB}.sample_sheet.txt"
				out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.cuffnorm_output/genes.attr_table"
			}
			out_fp.puts "#\tEXECUTION:"
			comparisons.each { |comparison|
				characteristic = comparison["characteristic"]
				condA = comparison["condition_A"]
				condB = comparison["condition_B"]
				out_fp.puts "./run_monocle.R #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffnorm_output/genes.fpkm_table.reordered #{output_dir}/#{characteristic}-#{condA}-#{condB}.sample_sheet.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffnorm_output/genes.attr_table #{output_dir}/#{characteristic}-#{condA}-#{condB}.monocle.pdf"
			}
			out_fp.puts "#\tOUTPUT:"
			comparisons.each { |comparison|
				characteristic = comparison["characteristic"]
				condA = comparison["condition_A"]
				condB = comparison["condition_B"]
				out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.monocle.pdf"
			}
			comparisons.each { |comparison|
				characteristic = comparison["characteristic"]
				condA = comparison["condition_A"]
				condB = comparison["condition_B"]
				out_fp.puts "[ -s \"#{output_dir}/#{characteristic}-#{condA}-#{condB}.monocle.pdf\" ] || (echo \"ERROR: #{output_dir}/#{characteristic}-#{condA}-#{condB}.monocle.pdf is empty!\" && exit)"
			}
		end
	end
elsif de_pipeline == "featureCounts"
	out_fp.puts "","","################################################"
	out_fp.puts "###\t2. Differential expression analysis using featureCounts/scde pipeline"
	out_fp.puts "###\ta. Count reads using featureCounts"
	out_fp.puts "#\tINPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{sample_id}/accepted_hits.bam"
	}
	out_fp.puts "#\tREQUIRED DATA FILES:"
	out_fp.puts "#\t\t#{genes_gtf_file}"
	out_fp.puts "#\tEXECUTION:"
	sub_fps.clear
	out_fp.puts "mkdir -p #{output_featureCounts_dir}"
	1.upto(num_jobs) do |i|
		sub_fps << File.new("#{output_code_dir}/workflow.#{project}.featureCounts_count.#{i}.sh", "w")
		out_fp.puts "bash workflow.#{project}.featureCounts_count.#{i}.sh &> workflow.#{project}.featureCounts_count.#{i}.log &"
	end
	out_fp.puts "wait"
	i = 0
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
		if libtype == "fr-unstranded"
			libtype_str = "0"
		elsif libtype == "fr-firststrand"
			libtype_str = "2"
		elsif libtype == "fr-secondstrand"
			libtype_str = "1"
		else
			puts "ERROR: invalid library-type #{libtype}"
			exit -1
		end
		cmd = "samtools sort -n #{output_alignment_dir}/#{sample_id}/accepted_hits.bam -o #{output_alignment_dir}/#{sample_id}/accepted_hits.namesorted.bam"
		sub_fps[(i % num_jobs)].puts cmd
		cmd = "featureCounts -s #{libtype_str} -p -t exon -g gene_id -a #{genes_gtf_file} -o #{output_featureCounts_dir}/#{sample_id}.featureCounts.counts.txt #{output_alignment_dir}/#{sample_id}/accepted_hits.namesorted.bam"
		sub_fps[(i % num_jobs)].puts cmd
		i += 1
	}
	sub_fps.each { |fp|
		fp.close
	}
	out_fp.puts "#\tOUTPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{sample_id}.featureCounts.counts.txt"
	}
	out_fp.puts "#\tCHECKPOINT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "[ -f \"#{output_featureCounts_dir}/#{sample_id}.featureCounts.counts.txt\" ] || (echo \"ERROR: #{output_featureCounts_dir}/#{sample_id}.featureCounts.counts.txt not found!\" && exit)"
	}
	
	if flag_comparison
		out_fp.puts "", "###\tb. Create read count matrix for each comparison"
		out_fp.puts "#\tINPUT:"
		out_fp.puts "#\t\t#{comparison_sheet_fn}"
		out_fp.puts "#\t\t#{sample_sheet_fn}"
		samples.each { |sample|
			sample_id = sample["sample_id"]
			out_fp.puts "#\t\t#{sample_id}.featureCounts.counts.txt"
		}
		out_fp.puts "#\tEXECUTION:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			condA_files = Array.new
			condB_files = Array.new
			condA_samples = Array.new
			condB_samples = Array.new
			factors = Array.new
			samples.each { |sample|
				sample_id = sample["sample_id"]
				if (sample[characteristic] == condA)
					condA_files << "#{output_featureCounts_dir}/#{sample_id}.featureCounts.counts.txt"
					condA_samples << sample_id
					factors << condA
				elsif (sample[characteristic] == condB)
					condB_files << "#{output_featureCounts_dir}/#{sample_id}.featureCounts.counts.txt"
					condB_samples << sample_id
					factors << condB
				end
			}
			header = [condA_samples, condB_samples].join("\t")
			out_fp.puts "mkdir #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/"
			out_fp.puts "echo \"#{header}\" > #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/read_counts.txt"
			header = [condA_files, condB_files].join(" ")
			farr = Array.new
			farr << 1
			1.upto(condA_files.length+condB_files.length) do |i|
				farr << i*2
			end
			fstr = farr.join(",")
			out_fp.puts "paste #{header} | cut -f#{fstr} >> #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/read_counts.txt"
			header = [Array.new(condA_files.length, condA), Array.new(condB_files.length, condB)].join(",")
			out_fp.puts "echo \"#{header}\" > #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/factors.txt"
		}
		out_fp.puts "#\tOUTPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.scde_output/read_counts.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.scde_output/factors.txt"
		}
	
		out_fp.puts "", "###\tc. Run scde"
		out_fp.puts "#\tINPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.scde_output/read_counts.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.scde_output/factors.txt"
		}
		out_fp.puts "#\tEXECUTION:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "./run_scde.R #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/read_counts.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/factors.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/results.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/results.pdf"
		}
		out_fp.puts "#\tOUTPUT:"
			comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.scde_output/results.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.scde_output/results.pdf"
		}
	
		out_fp.puts "", "###\td. Run DESeq2"
		out_fp.puts "#\tINPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.scde_output/read_counts.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.scde_output/factors.txt"
		}
		out_fp.puts "#\tEXECUTION:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "cp #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/read_counts.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.DESeq2_output/read_counts.txt"
			out_fp.puts "cp #{output_dir}/#{characteristic}-#{condA}-#{condB}.scde_output/factors.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.DESeq2_output/factors.txt"
			out_fp.puts "./run_DESeq2.R #{output_dir}/#{characteristic}-#{condA}-#{condB}.DESeq2_output/read_counts.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.DESeq2_output/factors.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.DESeq2_output/results.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.DESeq2_output/results.UP.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.DESeq2_output/results.DOWN.txt #{output_dir}/#{characteristic}-#{condA}-#{condB}.DESeq2_output/results.pdf"
		}
		out_fp.puts "#\tOUTPUT:"
			comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.DESeq2_output/results.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.DESeq2_output/results.pdf"
		}
	end
end

if analysis_type == "advanced"
	saturation_subsets = ["100000", "200000", "300000", "400000", "500000", "750000", "1000000", "1250000", "1500000", "2000000"]
	out_fp.puts "","","################################################"
	out_fp.puts "###\t4. Saturation analysis"
	out_fp.puts "###\ta. Subsample FASTQ files"
	out_fp.puts "#\tINPUT:"
	samples.each { |sample|
		out = sample["fastq"].split(",")
		out.each { |fastq|
			out_fp.puts "#\t\t#{fastq}"
		}
	}
	out_fp.puts "#\tEXECUTION:"
	out_fp.puts "mkdir -p #{output_saturation_dir}/fastq"
	subsets_str = saturation_subsets.join(" ")
	out_fp.puts "for n in #{subsets_str}", "do"
	samples.each { |sample|
		rseed = rand(1000000)
		sample["fastq"].split(",").each { |fastq|
			out_fp.puts "\tseqtk sample -s#{rseed} #{fastq_dir}/#{fastq} \$n > #{output_saturation_dir}/fastq/#{fastq}.sub_\$n"
		}
	}
	out_fp.puts "done"
	out_fp.puts "#\tOUTPUT:"
	out_fp.puts "#\t\t*.fastq.sub_* files"
#	samples.each { |sample|
#		out = sample["fastq"].split(",")
#		out.each { |fastq|
#			saturation_subsets.each { |n|
#				out_fp.puts "#\t\t#{fastq}.sub_#{n}"
#			}
#		}
#	}
#	out_fp.puts "#\tCHECKPOINT:"
#	samples.each { |sample|
#		out = sample["fastq"].split(",")
#		out.each { |fastq|
#			saturation_subsets.each { |n|
#				out_fp.puts "[ -f \"#{fastq_dir}/#{fastq}.sub_#{n}\" ] || (echo \"ERROR: #{fastq}.sub_#{n} not found!\" && exit)"
#			}
#		}
#	}
	
	out_fp.puts "", "###\tb. Realign subsampled reads"
	out_fp.puts "#\tINPUT:"
	out_fp.puts "#\t\t*.fastq.sub_* files"
	out_fp.puts "#\tREQUIRED DATA FILES"
	out_fp.puts "#\t\t*.bt2 files"
	out_fp.puts "#\t\t#{genes_gtf_file}"
	out_fp.puts "#\tEXECUTION:"
	sub_fps.clear
	out_fp.puts "mkdir -p #{output_saturation_dir}/#{aligner}"
	1.upto(num_jobs) do |i|
		sub_fps << File.new("#{output_code_dir}/workflow.#{project}.saturation.#{i}.sh", "w")
		out_fp.puts "bash workflow.#{project}.saturation.#{i}.sh &> workflow.#{project}.saturation.#{i}.log &"
	end
	out_fp.puts "wait"
	i = 0
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
		saturation_subsets.each { |saturation_subset|
			if aligner == "tophat"
				fastq_str = sample["fastq"].split(",").collect { |x| "#{output_saturation_dir}/fastq/#{x}.sub_#{saturation_subset}" }.join(" ")
				cmd = "tophat --read-mismatches #{read_mismatches} --read-gap-length #{read_gap_length} --read-edit-dist #{read_edit_dist} --max-multihits #{max_multihits} --library-type #{libtype} --GTF #{genes_gtf_file} --transcriptome-index #{genes_gtf_index} #{str_transcriptomeOnly} --num-threads #{num_alignment_threads} --output-dir #{output_saturation_dir}/#{aligner}/#{sample_id}.sub_#{saturation_subset} /Lab_Share/iGenomes/#{genome}/Sequence/Bowtie2Index/genome #{fastq_str}"
			elsif aligner == "hisat2"
				puts "mkdir -p #{output_saturation_dir}/#{aligner}/#{sample_id}.sub_#{saturation_subset}"
				arr = sample["fastq"].split(",")
				if arr.length == 1
					fastq_str = "-U #{fastq_dir}/#{arr}"
				elsif arr.length == 2
					fastq_str = "-1 #{fastq_dir}/#{arr[0]} -2 #{fastq_dir}/#{arr[1]}"
				else
					puts "ERROR: invalid fastq string!"
					exit -1
				end
				cmd = "hisat2 --num-threads #{num_alignment_threads} -x /Lab_Share/iGenomes/#{genome}/Sequence/HISAT2Index/genome #{fastq_str} 2> #{output_saturation_dir}/#{aligner}/#{sample_id}.sub_#{saturation_subset}/align_summary.txt | samtools view -bS -o #{output_saturation_dir}/#{aligner}/#{sample_id}.sub_#{saturation_subset}/accepted_hits.bam -"
			end
			sub_fps[(i % num_jobs)].puts(cmd)
		}
		i += 1
	}
	sub_fps.each { |fp|
		fp.close
	}
	out_fp.puts "#\tOUTPUT:"
	out_fp.puts "#\t\t#{output_saturation_dir}/#{aligner}/*"
		
	out_fp.puts "", "###\tc. Cuffquant"
	out_fp.puts "#\tINPUT:"
	out_fp.puts "#\t\t#{output_cuffmerge_dir}/merged.gtf"
	out_fp.puts "#\t\t#{output_saturation_dir}/#{aligner}/*"
	out_fp.puts "#\tEXECUTION:"
	out_fp.puts "mkdir -p #{output_saturation_dir}/cuffquant"
	subset_str = saturation_subsets.join(" ")
	out_fp.puts "for pct in #{subset_str}"
	out_fp.puts "do"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
#		out_fp.puts "\tcufflinks --output-dir #{output_dir}/#{sample_id}.sub_\$pct.cufflinks_output --num-threads #{num_alignment_threads} --GTF #{genes_gtf_file} --mask-file /Lab_Share/iGenomes/#{genome}/Annotation/Genes/rmsk_rRNA.gtf --multi-read-correct --library-type #{libtype} --upper-quartile-norm --compatible-hits-norm --quiet --no-update-check #{output_dir}/#{sample_id}.sub_\$pct.tophat_output/accepted_hits.bam"
		out_fp.puts "\tcuffquant --output-dir #{output_saturation_dir}/cuffquant/#{sample_id}.sub_\$pct -p #{num_alignment_threads} --mask-file /Lab_Share/iGenomes/#{genome}/Annotation/Genes/rmsk_rRNA.gtf --multi-read-correct --library-type #{libtype} --quiet --no-update-check #{output_cuffmerge_dir}/merged.gtf #{output_saturation_dir}/#{aligner}/#{sample_id}.sub_\$pct/accepted_hits.bam"
	}
	out_fp.puts "done"
	
	out_fp.puts "", "###\td. Cuffnorm"
	out_fp.puts "### PLACEHOLDER - insert manual Cuffnorm code here ###"
	out_fp.puts "### Next step requires cuffnorm's genes.fpkm_table.reordered file for each saturation data point ###"
	
	out_fp.puts "", "###\te. Compute saturation numbers and draw figures"
	out_fp.puts "#\tINPUT:"
#	samples.each { |sample|
#		sample_id = sample["sample_id"]
#		out_fp.puts "#\t\t#{sample_id}.cufflinks_output"
#	}
#	samples.each { |sample|
#		sample_id = sample["sample_id"]
#		saturation_subsets.each { |subset|
#			out_fp.puts "#\t\t#{sample_id}.sub_#{subset}.cufflinks_output/"
#		}
#	}
	out_fp.puts "#\tEXECUTION:"
#	prefix_arr = Array.new
#	full_file_arr = Array.new
#	samples.each { |sample|
#		sample_id = sample["sample_id"]
#		prefix_arr << "#{output_dir}/#{sample_id}"
#		full_file_arr << "#{output_dir}/#{sample_id}.cufflinks_output/genes.fpkm_tracking"
#	}
#	prefix_str = prefix_arr.join(",")
#	full_file_str = full_file_arr.join(",")
	out_fp.puts "./saturation_analysis.R [genes.fpkm_table.reordered FILES] #{output_saturation_dir}/saturation_analysis.pdf"
#	out_fp.puts "./saturation_analysis.R #{prefix_str} #{full_file_str} #{output_dir}/saturation_analysis.pdf"
	out_fp.puts "#\tOUTPUT:"
	out_fp.puts "#\t\tsaturation_analysis.pdf"
end




