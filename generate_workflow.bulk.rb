#!/usr/bin/ruby

# generate the basic/advanced traditional RNA-seq workflow as a shell script
# requires a sample sheet file (see included example SAMPLE_SHEET.txt) and a comparison sheet (COMPARISON.txt)

flag_transcriptomeOnly = false
flag_cuffmerge = true
flag_kallisto = false
flag_comparison = false
aligner = "hisat2"
de_pipeline = "kallisto" # "tuxedo" for Cufflinks/Cuffdiff/cummeRbund, "htseq" for HTseq/DESeq2, "kallisto" for kallisto/DESeq2
read_mismatches = 2
read_gap_length = 2
read_edit_dist = 2
max_multihits = 10
num_jobs = 1 # used to split tophat/QC/cuffetc steps into multiple batches
num_alignment_threads = 6	# used for Tuxedo suite multithreading
preprocess_fastq = false

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
		elsif arg == "--output-htseq-dir"
			output_htseq_dir = ARGV.shift
			output_htseq_dir = output_htseq_dir.gsub(/\/$/, "")
		elsif arg == "--output-code-dir"
			output_code_dir = ARGV.shift
			output_code_dir = output_code_dir.gsub(/\/$/, "")
		elsif arg == "--output_cuffmerge_dir"
			output_cuffmerge_dir = ARGV.shift
			output_cuffmerge_dir = output_cuffmerge_dir.gsub(/\/$/, "")
		elsif arg == "--output_cuffdiff_dir"
			output_cuffdiff_dir = ARGV.shift
			output_cuffdiff_dir = output_cuffdiff_dir.gsub(/\/$/, "")
		elsif arg == "--output_deseq2_dir"
			output_deseq2_dir = ARGV.shift
			output_deseq2_dir = output_deseq2_dir.gsub(/\/$/, "")
		elsif arg == "--output_edger_dir"
			output_edger_dir = ARGV.shift
			output_edger_dir = output_edger_dir.gsub(/\/$/, "")
		elsif arg == "--output_gsea_dir"
			output_gsea_dir = ARGV.shift
			output_gsea_dir = output_gsea_dir.gsub(/\/$/, "")
		elsif arg == "--output_kallisto_dir"
			output_kallisto_dir = ARGV.shift
			output_kallisto_dir = output_kallisto_dir.gsub(/\/$/, "")
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
		elsif arg == "--kallisto"
			flag_kallisto = true
		elsif arg == "--skip-cuffmerge"
			flag_cuffmerge = false
		elsif arg == "--genes-gtf"
			genes_gtf_file = ARGV.shift
			genes_gtf_index = genes_gtf_file.gsub("\.gtf","")
		elsif arg == "--preprocess-input-fastq"
			input_fastq_dir = ARGV.shift
			preprocess_fastq = true
		elsif arg == "--preprocess-output-fastq"
			output_fastq_dir = ARGV.shift
		elsif arg == "--kallisto-index"
			kallisto_index_file = ARGV.shift
		elsif arg == "--transcript-to-gene-file"
			transcript_to_gene_file = ARGV.shift
		elsif arg == "--aligner"
			aligner = ARGV.shift
		else
			puts "ERROR: Unrecognized option #{arg}"
			exit -1
		end
	end	
else
	puts "USAGE: #{$0} -s sample_sheet_file -c comparison_file -g hg19|mm10 --output-code-dir output_code_dir -o output_dir -a analysis_type -p project [--transcriptome-only] [--num-jobs num_jobs] [--num-tophat-threads num_tophat_threads] [--de-pipeline de_pipeline] [--genes-gtf genes_gtf_file]"
	exit -1
end

# error checking
if sample_sheet_fn.nil?
	puts "ERROR: sample_sheet_file must be specified"
	exit -1
end
if flag_comparison && comparison_sheet_fn.nil?
	puts "ERROR: comparison_file must be specified"
	exit -1
end
if output_dir.nil?
	puts "ERROR: output_dir must be specified"
	exit -1
end
# default locations are:
#		#{output_dir}/tophat
#		#{output_dir}/hisat2
#		#{output_dir}/cufflinks
#		#{output_dir}/QC
#		#{output_dir}/cuffquant
#		#{output_dir}/cuffmerge
#		#{output_dir}/cuffdiff
#		#{output_dir}/GSEA
#		#{output_dir}/htseq
#		#{output_dir}/DESeq2
#		#{output_dir}/edgeR
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
if output_cuffdiff_dir.nil?
	output_cuffdiff_dir = "#{output_dir}/cuffdiff"
end
if output_gsea_dir.nil?
	output_gsea_dir = "#{output_dir}/GSEA"
end
if output_htseq_dir.nil?
	output_htseq_dir = "#{output_dir}/htseq"
end
if output_deseq2_dir.nil?
	output_deseq2_dir = "#{output_dir}/DESeq2"
end
if output_edger_dir.nil?
	output_edger_dir = "#{output_dir}/edgeR"
end
if output_kallisto_dir.nil?
	output_kallisto_dir = "#{output_dir}/kallisto"
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
	if !(genome == "hg19" || genome == "GRCh38" || genome == "mm10")
		puts "ERROR: genome must be 'hg19', 'GRCh38', or 'mm10'"
		exit -1
	end
end
if !(de_pipeline == "tuxedo" || de_pipeline == "htseq" || de_pipeline == "kallisto")
	puts "ERROR: de-pipeline must be 'htseq', 'tuxedo', or 'kallisto'"
	exit -1
end
if (de_pipeline == "kallisto")
	flag_kallisto = true
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
out_fp.puts "#\truby generate_workflow.bulk.rb #{full_command}"
out_fp.puts "export PICARDPATH=/usr/local/picard-tools-1.118"
out_fp.puts "export GSEAPATH=/usr/local/GSEA"

if preprocess_fastq
	out_fp.puts "################################################"
	out_fp.puts "###\t0. Preprocess FASTQ files"
	out_fp.puts "for d in `find #{input_fastq_dir}/ -mindepth 1 -type d`"
	out_fp.puts "do"
	out_fp.puts "outfn=`find $d -name *R1* -printf \"%f\n\" | cut -d\".\" -f1 | cut -d\"_\" -f1,2,4,5 | uniq`"
  out_fp.puts "gunzip -c $d/*R1*.fastq.gz > #{output_fastq_dir}/$outfn.fastq"
	out_fp.puts "outfn=`find $d -name *R2* -printf \"%f\n\" | cut -d\".\" -f1 | cut -d\"_\" -f1,2,4,5 | uniq`"
  out_fp.puts "gunzip -c $d/*R2*.fastq.gz > #{output_fastq_dir}/$outfn.fastq"
	out_fp.puts "done"
	out_fp.puts "", ""
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
out_fp.puts "#\t\tHISAT2/tophat2 index files"
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
		cmd = "tophat --read-mismatches #{read_mismatches} --read-gap-length #{read_gap_length} --read-edit-dist #{read_edit_dist} --max-multihits #{max_multihits} --library-type #{libtype} --GTF #{genes_gtf_file} --transcriptome-index #{genes_gtf_index} #{str_transcriptomeOnly} --num-threads #{num_alignment_threads} --output-dir #{output_tophat_dir}/#{sample_id} /Lab_Share/iGenomes/#{genome}/Sequence/Bowtie2Index/genome #{fastq_str}"
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
		cmd = "hisat2 --num-threads #{num_alignment_threads} -x /Lab_Share/iGenomes/#{genome}/Sequence/HISAT2Index/genome #{fastq_str} | samtools view -bS -o #{output_alignment_dir}/#{sample_id}/accepted_hits.bam -"
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
sub_fps.clear
out_fp.puts "mkdir -p #{output_QC_dir}"
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

if flag_kallisto
	out_fp.puts "", "###\te. Run kallisto"
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
		samples.each { |sample|
			sample_id = sample["sample_id"]
			out_fp.puts "echo \"#{output_cufflinks_dir}/#{sample_id}/transcripts.gtf\" >> #{output_cuffmerge_dir}/assembly_list.txt"
		}
		out_fp.puts "cuffmerge -o #{output_cuffmerge_dir} --ref-gtf #{genes_gtf_file} -p #{num_alignment_threads} --ref-sequence /Lab_Share/iGenomes/#{genome}/Sequence/Chromosomes #{output_cuffmerge_dir}/assembly_list.txt"
		out_fp.puts "#\tOUTPUT:"
		out_fp.puts "#\t\tcuffmerge_output/merged.gtf"
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
		sub_fps[(i % num_jobs)].puts "cuffquant --output-dir #{output_cuffquant_dir}/#{sample_id}.cuffquant_output -p #{num_alignment_threads} --mask-file /Lab_Share/iGenomes/#{genome}/Annotation/Genes/rmsk_rRNA.gtf --multi-read-correct --library-type #{libtype} --quiet --no-update-check #{output_cuffmerge_dir}/merged.gtf #{output_alignment_dir}/#{sample_id}/accepted_hits.bam"
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
		out_fp.puts "mkdir -p #{output_cuffdiff_dir}"
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
			out_fp.puts "cuffdiff --output-dir #{output_cuffdiff_dir}/#{characteristic}-#{condA}-#{condB} --labels #{condA},#{condB} -p #{num_alignment_threads} --compatible-hits-norm --multi-read-correct --library-norm-method geometric --dispersion-method pooled --quiet --no-update-check #{output_cuffmerge_dir}/merged.gtf #{outA} #{outB}"
		}
		out_fp.puts "#\tOUTPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{output_cuffdiff_dir}/#{characteristic}-#{condA}-#{condB}/"
		}
		out_fp.puts "#\tCHECKPOINT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "[ -f \"#{output_cuffdiff_dir}/#{characteristic}-#{condA}-#{condB}/gene_exp.diff\" ] || (echo \"ERROR: #{output_cuffdiff_dir}/#{characteristic}-#{condA}-#{condB}/gene_exp.diff not found!\" && exit)"
		}

		out_fp.puts "", "###\te. cummeRbund"
		out_fp.puts "#\tINPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/"
		}
		out_fp.puts "#\tREQUIRED DATA FILES:"
		out_fp.puts "#\t\t#{output_cuffmerge_dir}/merged.gtf"
		out_fp.puts "#\tEXECUTION:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "./run_cummeRbund.R #{output_cuffdiff_dir}/#{characteristic}-#{condA}-#{condB} #{condA},#{condB} #{output_cuffmerge_dir}/merged.gtf #{genome} #{output_cuffdiff_dir}/#{characteristic}-#{condA}-#{condB}.cummeRbund.pdf"
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
			out_fp.puts "[ -s \"#{output_cuffdiff_dir}/#{characteristic}-#{condA}-#{condB}.cummeRbund.pdf\" ] || (echo \"ERROR: #{output_cuffdiff_dir}/#{characteristic}-#{condA}-#{condB}.cummeRbund.pdf is empty!\" && exit)"
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
		out_fp.puts "mkdir -p #{output_gsea_dir}"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "mkdir #{output_gsea_dir}/#{characteristic}-#{condA}-#{condB}/"
			out_fp.puts "if [ -d \"#{output_gsea_dir}/tmp.GSEA_output/\" ]; then", "\trmdir #{output_dir}/tmp.GSEA_output/" ,"fi"
			out_fp.puts "mkdir #{output_gsea_dir}/tmp.GSEA_output/"
			out_fp.puts "awk 'BEGIN {OFS=\"\\t\";FS=\"\\t\";} {if ($14==\"yes\") {print $3,log(($8+0.1)/($9+0.1))/log(2)}}' #{output_cuffdiff_dir}/#{characteristic}-#{condA}-#{condB}/gene_exp.diff > #{output_gsea_dir}/tmp.GSEA_output/cuffdiff_hits.rnk"
			out_fp.puts "java -cp $GSEAPATH/gsea2-2.1.0.jar -Xmx512m xtools.gsea.GseaPreranked -rnk #{output_gsea_dir}/tmp.GSEA_output/cuffdiff_hits.rnk -gmx $GSEAPATH/msigdb.v4.0.symbols.gmt -scoring_scheme classic -collapse false -set_min 15 -out #{output_gsea_dir}/tmp.GSEA_output -rpt_label \"tmp\" -zip_report true"
			out_fp.puts "mv #{output_gsea_dir}/tmp.GSEA_output/* #{output_gsea_dir}/#{characteristic}-#{condA}-#{condB}/"
			out_fp.puts "mv #{output_gsea_dir}/#{characteristic}-#{condA}-#{condB}/tmp.GseaPreranked.* #{output_gsea_dir}/#{characteristic}-#{condA}-#{condB}.GSEA_output/#{characteristic}-#{condA}-#{condB}.GseaPreranked/"
			out_fp.puts "mv #{output_gsea_dir}/#{characteristic}-#{condA}-#{condB}/#{characteristic}-#{condA}-#{condB}.GseaPreranked/*.zip #{output_gsea_dir}/#{characteristic}-#{condA}-#{condB}/report.zip"
			out_fp.puts "rmdir #{output_gsea_dir}/tmp.GSEA_output/"
		}
		out_fp.puts "#\tOUTPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/report.zip"
		}
	end
elsif de_pipeline == "htseq"
	out_fp.puts "","","################################################"
	out_fp.puts "###\t2. Differential expression analysis using Tophat/HISAT2+HTSeq+DESeq2/edgeR pipeline"
	out_fp.puts "###\ta. Count reads using HTSeq"
	out_fp.puts "#\tINPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{output_alignment_dir}/#{sample_id}/accepted_hits.bam"
	}
	out_fp.puts "#\tREQUIRED DATA FILES:"
	out_fp.puts "#\t\t#{genes_gtf_file}"
	out_fp.puts "#\tEXECUTION:"
	out_fp.puts "mkdir -p #{output_htseq_dir}"
	sub_fps.clear
	1.upto(num_jobs) do |i|
		sub_fps << File.new("#{output_code_dir}/workflow.#{project}.htseq_count.#{i}.sh", "w")
		out_fp.puts "bash workflow.#{project}.htseq_count.#{i}.sh &> workflow.#{project}.htseq_count.#{i}.log &"
	end
	out_fp.puts "wait"
	i = 0
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
		if libtype == "fr-unstranded"
			libtype_str = "no"
		elsif libtype == "fr-firststrand"
			libtype_str = "reverse"
		elsif libtype == "fr-secondstrand"
			libtype_str = "yes"
		else
			puts "ERROR: invalid library-type #{libtype}"
			exit -1
		end
		cmd = "samtools sort -n #{output_alignment_dir}/#{sample_id}/accepted_hits.bam #{output_alignment_dir}/#{sample_id}/accepted_hits.namesorted"
		sub_fps[(i % num_jobs)].puts cmd
		cmd = "htseq-count -f bam -r name -s #{libtype_str} -m union #{output_alignment_dir}/#{sample_id}/accepted_hits.namesorted.bam #{genes_gtf_file} > #{output_htseq_dir}/#{sample_id}.HTSeq.counts.txt"
		sub_fps[(i % num_jobs)].puts cmd
		i += 1
	}
	sub_fps.each { |fp|
		fp.close
	}
	out_fp.puts "#\tOUTPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{sample_id}.HTSeq.counts.txt"
	}
	out_fp.puts "#\tCHECKPOINT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "[ -f \"#{output_htseq_dir}/#{sample_id}.HTSeq.counts.txt\" ] || (echo \"ERROR: #{output_htseq_dir}/#{sample_id}.HTSeq.counts.txt not found!\" && exit)"
	}
	
	if flag_comparison
		out_fp.puts "", "###\tb. Create read count matrix for each comparison"
		out_fp.puts "#\tINPUT:"
		out_fp.puts "#\t\tCOMPARISON.txt"
		out_fp.puts "#\t\tSAMPLE_SHEET.txt"
		samples.each { |sample|
			sample_id = sample["sample_id"]
			out_fp.puts "#\t\t#{sample_id}.HTSeq.counts.txt"
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
					condA_files << "#{output_htseq_dir}/#{sample_id}.HTSeq.counts.txt"
					condA_samples << sample_id
					factors << condA
				elsif (sample[characteristic] == condB)
					condB_files << "#{output_htseq_dir}/#{sample_id}.HTSeq.counts.txt"
					condB_samples << sample_id
					factors << condB
				end
			}
			header = [condA_samples, condB_samples].join("\t")
			out_fp.puts "mkdir -p #{output_deseq2_dir}"
			out_fp.puts "mkdir -p #{output_edger_dir}"
			out_fp.puts "mkdir -p #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/"
			out_fp.puts "mkdir -p #{output_edger_dir}/#{characteristic}-#{condA}-#{condB}/"
			out_fp.puts "echo \"#{header}\" > #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/read_counts.txt"
			header = [condA_files, condB_files].join(" ")
			farr = Array.new
			farr << 1
			1.upto(condA_files.length+condB_files.length) do |i|
				farr << i*2
			end
			fstr = farr.join(",")
			out_fp.puts "paste #{header} | cut -f#{fstr} >> #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/read_counts.txt"
			header = [Array.new(condA_files.length, condA), Array.new(condB_files.length, condB)].join(",")
			out_fp.puts "echo \"#{header}\" > #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/factors.txt"
		
			out_fp.puts "cp #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/read_counts.txt #{output_edger_dir}/#{characteristic}-#{condA}-#{condB}/read_counts.txt"
			out_fp.puts "cp #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/factors.txt #{output_edger_dir}/#{characteristic}-#{condA}-#{condB}/factors.txt"
		}
		out_fp.puts "#\tOUTPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/read_counts.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/factors.txt"
		}
	
		out_fp.puts "", "###\tc. Run DESeq2/edgeR"
		out_fp.puts "#\tINPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/read_counts.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/factors.txt"
		}
		out_fp.puts "#\tEXECUTION:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "./run_DESeq2.R #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/read_counts.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/factors.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/results.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/results.UP.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/results.DOWN.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/results.pdf"
			out_fp.puts "./run_edgeR.R #{output_edger_dir}/#{characteristic}-#{condA}-#{condB}/read_counts.txt #{output_edger_dir}/#{characteristic}-#{condA}-#{condB}/factors.txt #{output_edger_dir}/#{characteristic}-#{condA}-#{condB}/results.txt #{output_edger_dir}/#{characteristic}-#{condA}-#{condB}/results.UP.txt #{output_edger_dir}/#{characteristic}-#{condA}-#{condB}/results.DOWN.txt #{output_edger_dir}/#{characteristic}-#{condA}-#{condB}/results.pdf"
		}
		out_fp.puts "#\tOUTPUT:"
			comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/results.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/results.pdf"
		}
	end
elsif de_pipeline == "kallisto"
	if flag_comparison
		out_fp.puts "","","################################################"
		out_fp.puts "###\ta. Create read count matrix for each comparison"
		out_fp.puts "#\tINPUT:"
		out_fp.puts "#\t\tCOMPARISON.txt"
		out_fp.puts "#\t\tSAMPLE_SHEET.txt"
		samples.each { |sample|
			sample_id = sample["sample_id"]
			out_fp.puts "#\t\t#{sample_id}/abundance.gene.tsv"
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
					condA_files << "#{output_kallisto_dir}/#{sample_id}/abundance.gene.tsv"
					condA_samples << sample_id
					factors << condA
				elsif (sample[characteristic] == condB)
					condB_files << "#{output_kallisto_dir}/#{sample_id}/abundance.gene.tsv"
					condB_samples << sample_id
					factors << condB
				end
			}
			header = [condA_samples, condB_samples].join("\t")
			out_fp.puts "mkdir -p #{output_deseq2_dir}"
			out_fp.puts "mkdir -p #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/"
			out_fp.puts "echo \"#{header}\" > #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/read_counts.txt"
			header = [condA_files, condB_files].join(" ")
			farr = Array.new
			farr << 1
			1.upto(condA_files.length+condB_files.length) do |i|
				farr << i*2
			end
			fstr = farr.join(",")
			out_fp.puts "paste #{header} | cut -f#{fstr} >> #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/read_counts.txt"
			header = [Array.new(condA_files.length, condA), Array.new(condB_files.length, condB)].join(",")
			out_fp.puts "echo \"#{header}\" > #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/factors.txt"
		}
		out_fp.puts "#\tOUTPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/read_counts.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/factors.txt"
		}
	
		out_fp.puts "", "###\tb. Run DESeq2"
		out_fp.puts "#\tINPUT:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/read_counts.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/factors.txt"
		}
		out_fp.puts "#\tEXECUTION:"
		comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "./run_DESeq2.R #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/read_counts.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/factors.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/results.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/results.UP.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/results.DOWN.txt #{output_deseq2_dir}/#{characteristic}-#{condA}-#{condB}/results.pdf"
		}
		out_fp.puts "#\tOUTPUT:"
			comparisons.each { |comparison|
			characteristic = comparison["characteristic"]
			condA = comparison["condition_A"]
			condB = comparison["condition_B"]
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/results.txt"
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}/results.pdf"
		}
	end
end

if analysis_type == "advanced" && 0 == 1
	saturation_subsets = ["10", "20", "30", "40", "50", "60", "70", "80", "90"]
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
	samples.each { |sample|
		fastq_str = sample["fastq"].split(",").collect { |x| "#{output_dir}/#{x}" }.join(" ")
		out_fp.puts "bash subsample_fastq.sh #{fastq_str}"
	}
	out_fp.puts "#\tOUTPUT:"
	samples.each { |sample|
		out = sample["fastq"].split(",")
		out.each { |fastq|
			saturation_subsets.each { |subset|
				out_fp.puts "#\t\t#{fastq}.sub_#{subset}"
			}
		}
	}
	out_fp.puts "#\tCHECKPOINT:"
	samples.each { |sample|
		out = sample["fastq"].split(",")
		out.each { |fastq|
			saturation_subsets.each { |subset|
				out_fp.puts "[ -f \"#{output_dir}/#{fastq}.sub_#{subset}\" ] || (echo \"ERROR: #{fastq}.sub_#{subset} not found!\" && exit)"
			}
		}
	}
	
	out_fp.puts "", "###\tb. Realign subsampled reads"
	out_fp.puts "#\tINPUT:"
	samples.each { |sample|
		out = sample["fastq"].split(",")
		out.each { |fastq|
			saturation_subsets.each { |subset|
				out_fp.puts "#\t\t#{fastq}.sub_#{subset}"
			}
		}
	}
	out_fp.puts "#\tREQUIRED DATA FILES"
	out_fp.puts "#\t\t*.bt2 files"
	out_fp.puts "#\t\t#{genes_gtf_file}"
	out_fp.puts "#\tEXECUTION:"
	subset_str = saturation_subsets.join(" ")
	out_fp.puts "for pct in #{subset_str}"
	out_fp.puts "do"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
		fastqs = sample["fastq"].split(",")
		if aligner == "tophat"
			if fastqs.length == 1
				out_fp.puts "\ttophat --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 --max-multihits 10 --num-threads #{num_tophat_threads} --output-dir #{output_dir}/#{sample_id}.sub_\$pct.tophat_output /Lab_Share/iGenomes/#{genome}/Sequence/Bowtie2Index/genome #{output_dir}/#{fastqs[0]}.sub_\$pct"
			elsif fastqs.length == 2
				out_fp.puts "\ttophat --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 --max-multihits 10 --library-type #{libtype} --num-threads 6 --output-dir #{output_dir}/#{sample_id}.sub_\$pct.tophat_output /Lab_Share/iGenomes/#{genome}/Sequence/Bowtie2Index/genome #{output_dir}/#{fastqs[0]}.sub_\$pct #{output_dir}/#{fastqs[1]}.sub_\$pct"
			else
				out_fp.puts "ERROR: weird number of FASTQ files for sample #{sample_id} #{fastqs}"
				exit -1
			end
		elsif aligner == "hisat2"
			if fastqs.length == 1
				fastq_str = "-U #{fastq_dir}/#{arr}"
			elsif fastqs.length == 2
				fastq_str = "-1 #{fastq_dir}/#{arr[0]} -2 #{fastq_dir}/#{arr[1]}"
			else
				puts "ERROR: invalid fastq string!"
				exit -1
			end
			puts "mkdir -p #{output_alignment_dir}/#{sample_id}.sub_\$pct.hisat2_output"
			cmd = "hisat2 --num-threads #{num_alignment_threads} -x /Lab_Share/iGenomes/#{genome}/Sequence/HISAT2Index/genome #{fastq_str} | samtools view -bS -o #{output_alignment_dir}/#{sample_id}.sub_\$pct/accepted_hits.bam -"
		end
	}
	out_fp.puts "done"
	out_fp.puts "#\tOUTPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		saturation_subsets.each { |subset|
			out_fp.puts "#\t\t#{output_dir}/#{sample_id}.sub_#{subset}.#{aligner}_output"
		}
	}
	out_fp.puts "#\tCHECKPOINT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		saturation_subsets.each { |subset|
			out_fp.puts "[ -f \"#{output_dir}/#{sample_id}.sub_#{subset}.#{aligner}_output/accepted_hits.bam\" ] || (echo \"ERROR: #{output_dir}/#{sample_id}.sub_#{subset}.#{aligner}_output/accepted_hits.bam not found!\" && exit)"
		}
	}
	
	out_fp.puts "", "###\tc. Cufflinks/Cuffquant/Cuffdiff"
	out_fp.puts "#\tINPUT:"
	out_fp.puts "#\t\t#{output_dir}/cuffmerge_output/merged.gtf"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		saturation_subsets.each { |subset|
			out_fp.puts "#\t\t#{output_dir}/#{sample_id}.sub_#{subset}.#{aligner}_output"
		}
	}
	out_fp.puts "#\tEXECUTION:"
	subset_str = saturation_subsets.join(" ")
	out_fp.puts "for pct in #{subset_str}"
	out_fp.puts "do"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		libtype = sample["library-type"]
		out_fp.puts "\tcufflinks --output-dir #{output_dir}/#{sample_id}.sub_\$pct.cufflinks_output --num-threads #{num_alignment_threads} --GTF #{genes_gtf_file} --mask-file /Lab_Share/iGenomes/#{genome}/Annotation/Genes/rmsk_rRNA.gtf --multi-read-correct --library-type #{libtype} --upper-quartile-norm --compatible-hits-norm --quiet --no-update-check #{output_dir}/#{sample_id}.sub_\$pct.#{aligner}_output/accepted_hits.bam"
		out_fp.puts "\tcuffquant --output-dir #{output_dir}/#{sample_id}.sub_\$pct.cuffquant_output -p #{num_alignment_threads} --mask-file /Lab_Share/iGenomes/#{genome}/Annotation/Genes/rmsk_rRNA.gtf --multi-read-correct --library-type #{libtype} --quiet --no-update-check #{output_dir}/cuffmerge_output/merged.gtf #{output_dir}/#{sample_id}.sub_\$pct.#{aligner}_output/accepted_hits.bam"
	}
	comparisons.each { |comparison|
		characteristic = comparison["characteristic"]
		condA = comparison["condition_A"]
		condB = comparison["condition_B"]
		cxb_condA = Array.new
		cxb_condB = Array.new
		samples.each { |sample|
			sample_id = sample["sample_id"]
			if (sample[characteristic] == condA)
				cxb_condA << "#{output_dir}/#{sample_id}.sub_\$pct.cuffquant_output/abundances.cxb"
			elsif (sample[characteristic] == condB)
				cxb_condB << "#{output_dir}/#{sample_id}.sub_\$pct.cuffquant_output/abundances.cxb"
			end
		}
		outA = cxb_condA.join(",")
		outB = cxb_condB.join(",")
		out_fp.puts "cuffdiff --output-dir #{output_dir}/#{characteristic}-#{condA}-#{condB}.sub_\$pct.cuffdiff_output --labels #{condA},#{condB} -p #{num_alignment_threads} --compatible-hits-norm --multi-read-correct --library-norm-method geometric --dispersion-method pooled --quiet --no-update-check #{output_dir}/cuffmerge_output/merged.gtf #{outA} #{outB}"
	}
	out_fp.puts "done"
	out_fp.puts "#\tOUTPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		saturation_subsets.each { |subset|
			out_fp.puts "#\t\t#{sample_id}.sub_#{subset}.cufflinks_output/"
			out_fp.puts "#\t\t#{sample_id}.sub_#{subset}.cuffquant_output/"
		}
	}
	comparisons.each { |comparison|
		characteristic = comparison["characteristic"]
		condA = comparison["condition_A"]
		condB = comparison["condition_B"]
		saturation_subsets.each { |subset|
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.sub_#{subset}.cuffdiff_output/"
		}
	}
	out_fp.puts "#\tCHECKPOINT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		saturation_subsets.each { |subset|
			out_fp.puts "[ -f \"#{output_dir}/#{sample_id}.sub_#{subset}.cufflinks_output/genes.fpkm_tracking\" ] || (echo \"ERROR: #{output_dir}/#{sample_id}.sub_#{subset}.cufflinks_output/genes.fpkm_tracking not found!\" && exit)"
			out_fp.puts "[ -f \"#{output_dir}/#{sample_id}.sub_#{subset}.cuffquant_output/abundances.cxb\" ] || (echo \"ERROR: #{output_dir}/#{sample_id}.sub_#{subset}.cuffquant_output/abundances.cxb not found!\" && exit)"
		}
	}
	comparisons.each { |comparison|
		characteristic = comparison["characteristic"]
		condA = comparison["condition_A"]
		condB = comparison["condition_B"]
		saturation_subsets.each { |subset|
			out_fp.puts "[ -f \"#{output_dir}/#{characteristic}-#{condA}-#{condB}.sub_#{subset}.cuffdiff_output/gene_exp.diff\" ] || (echo \"ERROR: #{output_dir}/#{characteristic}-#{condA}-#{condB}.sub_#{subset}.cuffdiff_output/gene_exp.diff not found!\" && exit)"
		}
	}
	
	out_fp.puts "", "###\td. Compute saturation numbers and draw figures"
	out_fp.puts "#\tINPUT:"
	samples.each { |sample|
		sample_id = sample["sample_id"]
		out_fp.puts "#\t\t#{sample_id}.cufflinks_output"
	}
	comparisons.each { |comparison|
		characteristic = comparison["characteristic"]
		condA = comparison["condition_A"]
		condB = comparison["condition_B"]
		out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.cuffdiff_output/"
	}
	samples.each { |sample|
		sample_id = sample["sample_id"]
		saturation_subsets.each { |subset|
			out_fp.puts "#\t\t#{sample_id}.sub_#{subset}.cufflinks_output/"
		}
	}
	comparisons.each { |comparison|
		characteristic = comparison["characteristic"]
		condA = comparison["condition_A"]
		condB = comparison["condition_B"]
		saturation_subsets.each { |subset|
			out_fp.puts "#\t\t#{characteristic}-#{condA}-#{condB}.sub_#{subset}.cuffdiff_output/"
		}
	}
	out_fp.puts "#\tEXECUTION:"
	prefix_arr = Array.new
	full_file_arr = Array.new
	samples.each { |sample|
		sample_id = sample["sample_id"]
		prefix_arr << "#{output_dir}/#{sample_id}"
		full_file_arr << "#{output_dir}/#{sample_id}.cufflinks_output/genes.fpkm_tracking"
	}
	prefix_str = prefix_arr.join(",")
	full_file_str = full_file_arr.join(",")
	out_fp.puts "./saturation_analysis.R #{prefix_str} #{full_file_str} #{output_dir}/saturation_analysis.pdf"
	prefix_arr = Array.new
	full_file_arr = Array.new
	comparisons.each { |comparison|
		characteristic = comparison["characteristic"]
		condA = comparison["condition_A"]
		condB = comparison["condition_B"]
		prefix_arr << "#{output_dir}/#{characteristic}-#{condA}-#{condB}"
		full_file_arr << "#{output_dir}/#{characteristic}-#{condA}-#{condB}.cuffdiff_output/gene_exp.diff"
	}
	prefix_str = prefix_arr.join(",")
	full_file_str = full_file_arr.join(",")
	out_fp.puts "./saturation_analysis_DE.R #{prefix_str} #{full_file_str} #{output_dir}/saturation_analysis_DE.pdf"
	out_fp.puts "#\tOUTPUT:"
	out_fp.puts "#\t\tsaturation_analysis.pdf"
	out_fp.puts "#\t\tsaturation_analysis_DE.pdf"
end




