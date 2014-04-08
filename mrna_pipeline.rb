#!/usr/bin/env ruby

### # # # # # # # # # # # # # # #
### sorghum mrna pipeline
### 
### author: Chris Boursnell (cmb211@cam.ac.uk)
### created: 2014-03-12 
###
### # # # # # # # # # # # # # # #

require 'rubygems'
require 'trollop'

opts = Trollop::options do
  version "v0.1"
  banner <<-EOS
  mrna pipeline

  author: Chris Boursnell (cmb211@cam.ac.uk)
  EOS
  opt :input, "Input file describing fastq files", :required => true, :type => String
  opt :reference, "Reference fasta file", :required => true, :type => String
  opt :genome, "Reference is a genome"
  opt :transcriptome, "Reference is a transcriptome"
  opt :annotation, "gff file", :type => String
  opt :threads, "Number of threads", :default => 22, :type => :int
  opt :morehelp, "More help"
  opt :verbose, "Be verbose"
  opt :test, "Don't actually run anything"
end

Trollop::die :input, "must exist" if !File.exist?(opts[:input]) if (opts[:input] or !opts[:morehelp])
Trollop::die "Please specify either --genome or --transcriptome" if (!opts.genome and !opts.transcriptome)

def parse_fastqc_data(path) #rds
  data = {}
  insection = false
  section = nil
  headers = nil
  File.open("#{path}").each_line do |line|
    line.chomp!
    if line.start_with? '>>'
      # sections begin and end with '>>'-prefixed lines
      if insection && line == '>>END_MODULE'
        # section ended
        insection = false
        headers = nil
      else
        # new section
        section = line[2..-1].split(/\t/).first
        if section == 'Basic Statistics'
          data[section] = {}
        else
          data[section] = []
        end
        insection = true
        headers = []
      end
    elsif insection
      if line.start_with? '#'
        # headers, create arrays
        if line.start_with? '#T'
          # one line doesn't follow the same spec
          # as the rest of the file. nice going,
          # fastqc devs ;)
          # data[section]['Total Duplication Percentage'] = line.split(/\t/).last
          next
        end
        line[1..-1].split(/\t/).each do |header|
          # data[section][header] = []
          headers << header
        end
      else
        # data, populate arrays
        if section == 'Basic Statistics'
          (name, value) = line.split(/\t/)
          data[section][name] = value
        else
          h = {}
          line.split(/\t/).each_with_index do |datum, index|
            h[headers[index]] = datum
          end
          data[section] << h
        end
      end
    else
      next
    end
  end
  data
end

if opts.morehelp
  puts "    The input file format should contain 1 line per input fastq file"
  puts "    and contain 3 comma separated columns."
  puts "    1) The fastq file location (absolute is better than relative)"
  puts "    2) Replicate number"
  puts "    3) Cell type or data type"
  puts "    4) Pair, either 1 or 2"
  puts "    For example:"
  puts "    /home/chris/fastq/B_1_R1.fq,1,B,1"
  puts "    /home/chris/fastq/B_1_R2.fq,1,B,2"
  puts "    /home/chris/fastq/M_1_R1.fq,1,M,1"
  puts "    /home/chris/fastq/M_1_R2.fq,1,M,2"
  exit
end
fastqc_path = "/applications/fastqc_v0.10.1/FastQC/fastqc"

### load data description from sorghum/raw_mrna_data

input = []
File.open("#{opts.input}").each_line do |line|
  cols = line.chomp.split(",")
  input << {:file => cols[0], :rep => cols[1], :cell => cols[2], :pair => cols[3]}
end

### run fastqc on the files

files = ""
input.each do |hash|
  files << " #{hash[:file]} "
end
output_dir = "mrna_fastqc_output"
if !Dir.exists?("#{output_dir}")
  fastqc = "#{fastqc_path} --kmers 7 --threads #{input.length} --outdir #{output_dir} #{files}"
  puts fastqc
  puts "Running fastqc" if opts.verbose
  `mkdir #{output_dir}`
  `#{fastqc}`  if !opts.test
end

### parse the fastqc output and get information about
###   number of reads
###   read length distribution
###   overrepresented sequences

input.each_with_index do |hash,index|
  if Dir.exists?("#{output_dir}/#{File.basename(hash[:file])}_fastqc") 
    data = parse_fastqc_data("#{output_dir}/#{File.basename(hash[:file])}_fastqc/fastqc_data.txt")
    name = data["Basic Statistics"]["Filename"]
    input[index][:count] = data["Basic Statistics"]["Total Sequences"].to_i
    input[index][:encoding] = data["Basic Statistics"]["Encoding"]
    input[index][:len] = data["Basic Statistics"]["Sequence length"].to_i
    input[index][:over] = []
    data["Overrepresented sequences"].each do |hash2|
      if hash2["Percentage"].to_f > 1 and hash2["Sequence"].length > 12
        input[index][:over] << hash2["Sequence"]
      end
    end
    input[index][:raw_quals] = []
    data["Per base sequence quality"].each do |hash2|
      input[index][:raw_quals] << hash2
    end
  end
end

### export number of reads per file as a table (R)

name = "mrna_raw_read_counts.txt"
if !File.exists?("#{name}") and !File.zero?("#{name}")
  File.open("#{name}", "w") do |io|
    input.each do |hash|
      io.write "#{hash[:count]}\t#{hash[:cell]}\t#{hash[:rep]}\t#{hash[:pair]}\n" if !opts.test
    end
  end
end

### export quality of reads to a table (R)

input.each_with_index do |hash, index|
  name = "raw_read_qualities_#{hash[:cell]}-#{hash[:rep]}-#{hash[:pair]}.txt"
  text=""
  hash[:raw_quals].each do |hash2|
    text << "#{hash2["Base"].split("-").first.to_i}\t#{hash2["Mean"].to_f}\t#{hash2["Lower Quartile"]}\t#{hash2["Upper Quartile"]}\t#{hash2["10th Percentile"]}\t#{hash2["90th Percentile"]}\n"
  end
  File.open("#{name}", "w") {|io| io.write(text)} if !File.exists?("#{name}") if !opts.test
end


### run minion on the reads to find if there are any adapters on the 3' end?

#adapter_list = []
#if !File.exists?("adapter_list.txt")
#  puts "Running minion" if opts.verbose
#  input.each do |hash|
#    minion_outfile = "adapters_#{hash[:cell]}-#{hash[:rep]}-#{hash[:pair]}.txt"
 #   if !File.exists?("#{minion_outfile}")
##     minion = "minion search-adapter -do #{hash[:count]-1} -i #{hash[:file]} " 
#      minion = "minion search-adapter -do 1000000 -i #{hash[:file]} " 
#      minion += " > #{minion_outfile}"
#      #puts minion
#      #`#{minion}`
#    end
#    adapter_list << minion_outfile
#  end
#  File.open("adapter_list.txt", "w") {|io| io.write(adapter_list.join("\n"))}
#end

### extract adapters from minion output

#count=0
#adapter_file = "mrna_minion_adapters.fasta"
#if !File.exists?("#{adapter_file}") and !File.zero?("#{adapter_file}")
#  File.open(adapter_file,"w") do |adapter_out|
#    puts "Creating adapters from minion output" if opts.verbose
#    adapter_list.each do |file|
#      File.open("#{file}").each_line do |line|
#        if line.match(/sequence=([ACGT]+)/)
#          adapter_out.write ">minion_adapter#{count}\n"
#          adapter_out.write "#{$1}\n"
#          count+=1
#        end
#      end
#    end
#    input.each_with_index do |hash,index|
#      hash[:over].each do |seq|
#        adapter_out.write ">adapter#{count}\n"
#        adapter_out.write "#{seq}\n"
#        count+=1
#      end
#    end
#  end
#end

# # # # # # # # # #
# run trimmomatic

trimmed=[]
run_trimmomatic = true
if run_trimmomatic
  trim_jar = "/home/cmb211/apps/Trimmomatic-0.32/trimmomatic-0.32.jar"
  threads = 12

  phred=" -phred33 "
  if input[0][:encoding].match(/Illumina 1.5/)
    phred=" -phred64 "
  elsif input[0][:encoding].match(/Illumina 1.8/)
    phred=" -phred33 "
  end
  minlen=40 # 
  windowsize=4
  quality=15
  trailing=15
  leading=15
  seed_mismatches=2 
  palindrome_clip_threshold=40 
  simple_clip=11
  #puts "Running Trimmomatic" if opts.verbose
  input.each_with_index.each_slice(2) do |(a,i), (b,j)|
    infile_left = a[:file]
    infile_right = b[:file]
    #outfile_left = "t.#{File.basename(infile_left)}"
    outfile_left = "#{a[:cell]}_#{a[:rep]}-#{a[:pair]}.t.fq"
    outfile_right = "#{a[:cell]}_#{a[:rep]}-#{b[:pair]}.t.fq"
    outfileU_left = "#{a[:cell]}_#{a[:rep]}-#{a[:pair]}.tU.fq"
    outfileU_right = "#{a[:cell]}_#{a[:rep]}-#{b[:pair]}.tU.fq"
    trim_cmd = "java -jar #{trim_jar} PE #{phred} "                     # phred
    # trim_cmd += " -trimlog trim_log_#{hash[:cell]}#{hash[:rep]}.txt "   # log output
    trim_cmd += " -threads #{threads} "
    trim_cmd += " #{infile_left} #{infile_right} "
    trim_cmd += " #{outfile_left} #{outfileU_left} #{outfile_right} #{outfileU_right}"
    #trim_cmd += " ILLUMINACLIP:#{adapter_file}:#{seed_mismatches}:#{palindrome_clip_threshold}:#{simple_clip} "
    trim_cmd += " LEADING:#{leading} TRAILING:#{trailing} "
    trim_cmd += " SLIDINGWINDOW:#{windowsize}:#{quality} MINLEN:#{minlen} "
    trimmed << outfile_left
    trimmed << outfile_right
    input[i][:trimmed] = outfile_left
    input[j][:trimmed] = outfile_right
    input[i][:unpaired] = outfileU_left
    input[j][:unpaired] = outfileU_right
    if !File.exists?("#{outfile_left}")
      puts trim_cmd if opts.verbose
      `#{trim_cmd}` if !opts.test
    else
      puts "trimmomatic already run"
    end
  end
end


## # # # # # # # # # # # # # # # # # #
## Run fastqc on all trimmed files
##

files = ""
trimmed.each do |tfile|
  files << " #{tfile} "
end

output_dir = "trimmed_mrna_fastqc_output"
if !Dir.exists?("#{output_dir}")
  fastqc = "#{fastqc_path} --kmers 7 --threads #{input.length} --outdir #{output_dir} #{files}"
  puts fastqc
  puts "Running fastqc" if opts.verbose
  `mkdir #{output_dir}` if !opts.test
  `#{fastqc}` if !opts.test
else
  puts "fastqc already run on trimmed reads"
end

### parse the fastqc output again and get information about
###   number of reads
###   read length distribution
###   overrepresented sequences

t_counts=""
input.each_with_index do |hash,index|
  dir = "#{output_dir}/#{hash[:cell]}_#{hash[:rep]}-#{hash[:pair]}.t.fq_fastqc"
  if Dir.exists?("#{dir}")
    data = parse_fastqc_data("#{dir}/fastqc_data.txt")
    count = data["Basic Statistics"]["Total Sequences"].to_i
    t_counts <<  "#{count}\t#{hash[:cell]}\t#{hash[:rep]}\t#{hash[:pair]}\n"
    input[index][:trim] = [] 
    data["Sequence Length Distribution"].each do |hash2|
      input[index][:trim] << hash2
    end
    input[index][:t_quals] = []
    data["Per base sequence quality"].each do |hash2|
      input[index][:t_quals] << hash2
    end
  end
end
File.open("mrna_trimmed_read_counts.txt", "w") {|io| io.write(t_counts)} if !File.exists?("mrna_trimmed_read_counts.txt") if !opts.test


# output trimmed read lengths per file
input.each_with_index do |hash, index|
  name = "t_read_length_#{hash[:cell]}-#{hash[:rep]}-#{hash[:pair]}.txt"
  if !File.exists?("#{name}")
    File.open("#{name}", "w") do |out|
      hash[:trim].each do |hash2|
        out.write "#{hash2["Length"]}\t#{hash2["Count"]}\n" if !opts.test
      end
    end
  else
    puts "#{name} already exists"
  end
end

### export quality of reads to a table (R)

input.each_with_index do |hash, index|
  name = "t_read_qualities_#{hash[:cell]}-#{hash[:rep]}-#{hash[:pair]}.txt"
  text=""
  hash[:t_quals].each do |hash2|
    text << "#{hash2["Base"].split("-").first.to_i}\t#{hash2["Mean"].to_f}\t#{hash2["Lower Quartile"]}"
    text << "\t#{hash2["Upper Quartile"]}\t#{hash2["10th Percentile"]}\t#{hash2["90th Percentile"]}\n"
  end
  File.open("#{name}", "w") {|io| io.write(text)} if !File.exists?("#{name}") if !opts.test
end

reference = opts.reference
if opts.genome
  # if the reference is a genome then build a transcriptome
  # using tophat and cufflinks
  left = []
  right = []
  input.each do |hash|
    if hash[:pair]=="1"
      left << hash[:trimmed]
    elsif hash[:pair]=="2"
      right << hash[:trimmed]
    end
  end
  cmd = "ruby /home/cmb211/scripts/sorghum/tophat.rb --left #{left.join(",")} --right #{right.join(",")} --genome #{opts.reference} "
  cmd += " --outputdir transcriptome --threads #{opts.threads} --verbose "
  puts "Building transcriptome from mRNAseq reads and genome" if opts.verbose
  puts " using Tophat and Cufflinks" if opts.verbose
  reference = "transcriptome/#{File.basename(opts.reference).split(".").first}-transcripts.fa"
  puts cmd if opts.verbose
  if !File.exists?("#{reference}")
    `#{cmd}` if !opts.test
  end
  puts "reference = #{reference}"
end

### align the trimmed mrna reads to the sorghum transcriptome
###   --met-file <path>  send metrics to file at <path> (off)
###   does the metrics output file contain information on number of reads aligned etc?
###   if it does put this into a plot in R

# bowtie2-build index
index = "#{File.basename(reference).split(".")[0..-2].join(".")}" 
# index = reference

if !File.exists?("#{index}.1.bt2")
  build = "bowtie2-build #{reference} #{index}"
  puts build if opts.verbose
  `#{build}`  if !opts.test
else
  puts "bowtie2 index already exists" if opts.verbose
end

#bowtie2 align reads to transcriptome
count=0
# input.each_slice(2) do |slice|
input.each_with_index.each_slice(2) do |(a,i), (b,j)|
  sam = "#{a[:cell]}_#{a[:rep]}.sam"
  bowtie = "bowtie2 -t -p #{opts.threads} --met-file met_file_#{count} --very-sensitive "
  # bowtie += " -a " # if this makes the memory go too large switch to -k 500
  bowtie += " -k 500 "
  bowtie += " -x #{index} "
  bowtie += " -1 #{a[:trimmed]} "
  bowtie += " -2 #{b[:trimmed]} "
  bowtie += " -U #{a[:unpaired]},#{b[:unpaired]}"
  bowtie += " -S transcriptome/#{sam}"
  puts bowtie if opts.verbose
  count+=1
  if !File.exists?("transcriptome/#{sam}")
    `#{bowtie}` if !opts.test
  else
    puts "#{sam} already exists" if opts.verbose
  end
  input[i][:sam] = sam
  input[j][:sam] = sam
end

### run expression quantification on the output sam files

input.each_with_index.each_slice(2) do |(a,i),(b,j)|
  out_dir = "express_#{a[:cell]}-#{a[:rep]}"
  sam = "#{a[:sam]}"
  if !File.exists?("#{out_dir}/results.xprs") # if the eXpress output doesn't exist run eXpress
    express_cmd = "express --output-align-prob "
    express_cmd << " -o #{out_dir} "
    express_cmd << " --no-update-check "   
    express_cmd << " -B 2 " 
    express_cmd << " #{reference} " # fasta file to align reads to
    express_cmd << " transcriptome/#{sam}"
    puts express_cmd if opts.verbose
    if !Dir.exists?("#{out_dir}")
      `#{express_cmd}` if !opts.test
    end
  end
end

# maybe do it just by parsing the sam file in ruby
# and just load the genes that we're interested in
# and then export a file for each gene/transcript
# with a list of the number of 

### construct a list of the important genes

key = Hash.new
File.open("key_sorghum.txt").each_line do |line|
  (gene, desc) = line.chomp.split("\t")
  key[gene] = desc
end

# p key

# a cigar string of 21M1D79M means add a gap to the read, or delete from the reference
# so 21M1I77M would mean add a gap to the reference.
# the total coverage of the read over the reference is the sum of M - I + D

### adapt this code to parse the eXpress output sam files and save the output with decimal coverage based on probabilities in XP:f:[0-9]+
input.each_with_index.each_slice(2) do |(a,i),(b,j)|
  count=0
  # sam = "transcriptome/#{a[:sam]}"
  sam = "express_#{a[:cell]}-#{a[:rep]}/hits.1.prob.sam"
  cell = "#{a[:cell]}_#{a[:rep]}"
  out = "#{cell}.mpileup"
  if !File.exists?("#{out}")
    puts "scanning #{sam}..."
    coverage = Hash.new
    File.open("#{sam}").each_line do |line|
      if line[0] != '@'
        count+=1
        print "." if count%10_000==0
        cols = line.chomp.split("\t")
        pos = cols[3].to_i
        gene = cols[2]
        align=0
        if line.match(/XP:f:(\S+)/)
          prob = $1.to_f
          cols[5].scan(/([0-9]+)([DIM])/).each do |a|
            count = a[0].to_i
            letter = a[1]
            align += count if letter=="M" or letter=="D"
            align -= count if letter=="I"
          end
          coverage[gene] = [] if !coverage.has_key?(gene)
          (pos..(pos+align)).each do |i|
            coverage[gene][i] = 0 if !coverage[gene][i]
            coverage[gene][i] += prob
          end
        end
      end
    end
    ### save output
    puts "saving output from #{sam}" if opts.verbose
    File.open("#{out}", "w") do |out|
      coverage.each_key do |gene|
        coverage[gene].each_with_index do |cov, index|
          cov=0 if !cov
          out.write "#{gene}\t#{index}\t#{cov}\n" if index>0
        end
      end
    end
  else
    puts "#{out} already exists" if opts.verbose
  end
  input[i][:coverage] = out
  input[j][:coverage] = out
  puts "...Done" if opts.verbose
end

### add introns into depth file

input.each_with_index.each_slice(2) do |(a,i),(b,j)|
  depth = a[:coverage]
  out = "#{a[:cell]}_#{a[:rep]}.coverage"
  introns = "ruby ~/scripts/depth/lib/depth.rb --depth #{depth} --annotation #{opts.annotation} --output #{out}"
  if !File.exists?("#{out}")
    puts introns if opts.verbose
    `#{introns}` if !opts.test
  else
    puts "#{out} already exists"
  end
end