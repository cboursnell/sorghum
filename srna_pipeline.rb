#!/usr/bin/env ruby

#
# srna pipeline
#
# author: Chris Boursnell (cmb211@cam.ac.uk)
# created: 2014-03-06
#
# performing quality checking on raw srna data
#

require 'rubygems'
require 'trollop'
require 'bio'

opts = Trollop::options do
  version "v0.1"
  banner <<-EOS
  srna pipeline

  author: Chris Boursnell (cmb211@cam.ac.uk)
  EOS
  opt :input, "Input file describing fastq files", :required => true, :type => String
  opt :reference, "Reference fasta file", :required => :true, :type => String
  opt :annotation, "Reference annotation gff file", :required => :true, :type => String
  opt :morehelp, "More help"
  opt :verbose, "Be verbose"
  opt :test, "Don't actually do anything"
end

Trollop::die :input, "must exist" if !File.exists?(opts[:input]) if (opts[:input] or !opts[:morehelp])
Trollop::die :reference, "must exist" if !File.exists?(opts[:reference]) if (opts[:reference] or !opts[:morehelp])
Trollop::die :annotation, "must exist" if !File.exists?(opts[:annotation]) if (opts[:annotation] or !opts[:morehelp])

class Feature
  attr_accessor :name, :type, :features, :start, :stop

  def initialize(name, type,start,stop)
    @name = name
    @type = type
    @features = []
    @start = start
    @stop= stop
  end

  def to_s
    "#{name}\t#{type}\t#{start}\t#{stop}"
  end
end

def parse_fastqc_data(path) #rds
  data = {}
  insection = false
  section = nil
  headers = nil
  # zipfile.read("#{dir}/fastqc_data.txt").split(/\n/).each do |line|
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
            # data[section][headers[index]] << datum
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
  puts "    The first is the fastq file location (absolute is better than relative)"
  puts "    The second is the replicate"
  puts "    The third and last is the cell type or data type"
  puts "    For example:"
  puts "    /home/chris/fastq/B_1.fq,1,1"
  puts "    /home/chris/fastq/B_2.fq,2,1"
  puts "    /home/chris/fastq/B_3.fq,3,1"
  puts "    /home/chris/fastq/M_1.fq,1,2"
  puts "    /home/chris/fastq/M_2.fq,2,2"
  puts "    /home/chris/fastq/M_3.fq,3,2"
  exit
end
fastqc_path = "/applications/fastqc_v0.10.1/FastQC/fastqc"

## # # # # # # # # # # # # # # # # # #
## Open input file and store contents
puts "opening input file" if opts.verbose
input = []
File.open("#{opts.input}").each_line do |line|
  cols = line.chomp.split(",")
  input << {:file => cols[0], :rep => cols[1], :cell => cols[2]}
end

## # # # # # # # # # # # # # # # # # #
## Run minion on the files individually
## and record the adapters

adapter_list = []
if !File.exists?("adapter_list.txt")
  puts "Running minion" if opts.verbose
  input.each do |hash|
    #minion search-adapter -do 3000000 -i ~/fastq_data/sorghum/Mroll-19_3_12.fq > ~/sorghum/minion-output-m-19.txt
    minion = "minion search-adapter -do #{hash[:len]-1} -i #{hash[:file]} > adapters_#{hash[:cell]}-#{hash[:rep]}.txt"
    puts minion
    adapter_list << "adapters_#{hash[:cell]}-#{hash[:rep]}.txt"
    `#{minion}` if !opts.test
  end
  File.open("adapter_list.txt", "w") {|io| io.write(adapter_list.join("\n"))}
else
  puts "Minion already run"
end

## # # # # # # # # # # # # # # # # # #
## Run fastqc on all files
## 

files = ""
input.each do |hash|
  files << " #{hash[:file]} " 
end

if !Dir.exists?("fastqc_output")
  fastqc = "#{fastqc_path} --kmers 7 --threads #{input.length} --outdir fastqc_output #{files}"
  puts fastqc
  puts "Running fastqc" if opts.verbose
  `mkdir fastqc_output`
  `#{fastqc}` if !opts.test
else
  puts "fastqc already run on raw reads"
end

input.each_with_index do |hash,index|
  if Dir.exists?("fastqc_output/#{File.basename(hash[:file])}_fastqc")
    puts "parsing data for #{hash[:file]}"
    data = parse_fastqc_data("fastqc_output/#{File.basename(hash[:file])}_fastqc/fastqc_data.txt")
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
  else
    puts "can't find fastqc directory for #{File.basename(hash[:file])}_fastqc"
  end
end
# {"Basic Statistics"=>{"Filename"=>"BS-1b-1.fq", "File type"=>"Conventional base calls", 
# "Encoding"=>"Illumina 1.5", "Total Sequences"=>"22932596", "Filtered Sequences"=>"0", "Sequence length"=>"49", "%GC"=>"52"}

## # # # # # # # # # # # # # # # # # #
## Find number of reads in each file
## and export as a table

name = "raw_read_counts.txt"
if !File.exists?("#{name}") and !File.zero?("#{name}")
  File.open("#{name}", "w") do |io|
    input.each do |hash|
      io.write "#{hash[:count]}\t#{hash[:cell]}\t#{hash[:rep]}\n"  if !opts.test
    end
  end
else
  puts "#{name} already exists"
end

## # # # # # # # # # # # # # # # # # #
## Extract adapters from minion output
##
count=0
adapter_file = "minion_adapters.fasta"
if !File.exists?("#{adapter_file}") and !File.zero?("#{adapter_file}")
  File.open(adapter_file,"w") do |adapter_out|
    puts "Creating adapters from minion output" if opts.verbose
    input.each_with_index do |hash,index|
      File.open("adapters_#{hash[:cell]}-#{hash[:rep]}.txt").each_line do |line|
        if line.match(/sequence=([ACGT]+)/)
          adapter_out.write ">adapter#{count}\n"  if !opts.test
          adapter_out.write "#{$1}\n"  if !opts.test
          count+=1
        end
      end
      hash[:over].each do |seq|
        adapter_out.write ">adapter#{count}\n"  if !opts.test
        adapter_out.write "#{seq}\n"  if !opts.test
        count+=1
      end
    end
  end
else
  puts "#{adapter_file} already created"
end

if File.zero?("#{adapter_file}")
  abort "Something went wrong with making the fasta file of adapters"
end

## # # # # # # # # # # # # # # # # # #
## Trim reads using trimmomatic
## and adapters from minion

trimmed=[]
run_trimmomatic = true
if run_trimmomatic 
  trim_jar = "/home/cmb211/apps/Trimmomatic-0.32/trimmomatic-0.32.jar"
  threads = 4

  phred=" -phred33 "
  if input[0][:encoding].match(/Illumina 1.5/)
    phred=" -phred64 "
  elsif input[0][:encoding].match(/Illumina 1.8/)
    phred=" -phred33 "
  end
  minlen=17
  windowsize=4
  quality=15
  trailing=15
  leading=15
  seed_mismatches=2 #  specifies the maximum mismatch count which will still allow a full match to be performed
  palindrome_clip_threshold=40 # pecifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment
  simple_clip=11 # specifies how accurate the match between any adapter etc. sequence must be against a read # i guessed 18 bases and multiplied it by 0.6 see below. /shrug

  # These values are (very roughly speaking) log-10 probabilities of getting a match at random. Each perfectly matching base 
  # scores just over 0.6, so 15 requires a perfect 25 base match. Each mismatching base reduces the score by the Q/10 value of
  # that base. So it takes 5 or even 6 additional matching bases to overcome one high quality mismatch, but maybe only 1 or 2 
  # additional bases if the mismatching base is low quality. [http://seqanswers.com/forums/showthread.php?t=11186]

  puts "Running Trimmomatic" if opts.verbose
  input.each_with_index do |hash,index|
    infile = hash[:file]
    outfile = "t.#{File.basename(infile)}"
    trim_cmd = "java -jar #{trim_jar} SE #{phred} "                     # phred
    # trim_cmd += " -trimlog trim_log_#{hash[:cell]}#{hash[:rep]}.txt "   # log output
    trim_cmd += " -threads #{threads} #{infile} #{outfile} "            # threads
    trim_cmd += " ILLUMINACLIP:#{adapter_file}:#{seed_mismatches}:#{palindrome_clip_threshold}:#{simple_clip} " 
    trim_cmd += " LEADING:#{leading} TRAILING:#{trailing} SLIDINGWINDOW:#{windowsize}:#{quality} MINLEN:#{minlen} "
    trimmed << outfile
    input[index][:trimmed_t] = outfile
    if !File.exists?("#{outfile}")
      puts trim_cmd if opts.verbose
      `#{trim_cmd}`  if !opts.test
    else
      puts "#{outfile} already created with trimmomatic"
    end
  end
end

# # FASTQ MCF # # # # # # # # # #

puts "Starting fastq-mcf" if opts.verbose
input.each_with_index do |hash,index|
  mcf = "fastq-mcf "
  # options
  mcf += " -o mcf.#{File.basename(hash[:file])} "
  mcf += " #{adapter_file}"
  mcf += " #{hash[:file]} "
  outfile = "mcf.#{File.basename(hash[:file])}"
  trimmed << outfile
  input[index][:trimmed_m] = outfile
  if !File.exists?("#{outfile}")
    puts mcf if opts.verbose
    `#{mcf}`  if !opts.test
  else
    puts "#{outfile} already created with fastq-mcf"
  end
end

## # # # # # # # # # # # # # # # # # #
## Run fastqc on all trimmed files
## 

files = ""
trimmed.each do |tfile|
  files << " #{tfile} " 
end

if !Dir.exists?("trimmed_fastqc_output")
  fastqc = "#{fastqc_path} --kmers 7 --threads #{input.length} --outdir trimmed_fastqc_output #{files}"
  puts fastqc
  puts "Running fastqc" if opts.verbose
  `mkdir trimmed_fastqc_output`  if !opts.test
  `#{fastqc}`  if !opts.test
else
  puts "fastqc already run on trimmed reads"
end

t_counts=""
mcf_counts=""

input.each_with_index do |hash,index|
  if Dir.exists?("trimmed_fastqc_output/t.#{File.basename(input[index][:file])}_fastqc")
    data = parse_fastqc_data("trimmed_fastqc_output/t.#{File.basename(input[index][:file])}_fastqc/fastqc_data.txt")
    count = data["Basic Statistics"]["Total Sequences"].to_i
    t_counts <<  "#{count}\t#{hash[:cell]}\t#{hash[:rep]}\n"
    input[index][:t_trim] = []
    data["Sequence Length Distribution"].each do |hash2|
      input[index][:t_trim] << hash2
    end
  end
  if Dir.exists?("trimmed_fastqc_output/mcf.#{File.basename(input[index][:file])}_fastqc")
    data = parse_fastqc_data("trimmed_fastqc_output/mcf.#{File.basename(input[index][:file])}_fastqc/fastqc_data.txt")
    count = data["Basic Statistics"]["Total Sequences"].to_i
    mcf_counts <<  "#{count}\t#{hash[:cell]}\t#{hash[:rep]}\n"
    input[index][:mcf_trim] = []
    data["Sequence Length Distribution"].each do |hash2|
      input[index][:mcf_trim] << hash2
    end
  end
end
File.open("t_read_counts.txt", "w") {|io| io.write(t_counts)} if !File.exists?("t_read_counts.txt")  if !opts.test
File.open("mcf_read_counts.txt", "w") {|io| io.write(mcf_counts)}  if !File.exists?("mcf_read_counts.txt")  if !opts.test

## # # # # # # # # # # # # # # # # # #
## Read lengths after trimming
##


input.each_with_index do |hash, index|
  # puts "#{hash[:cell]}\t#{hash[:rep]}"
  name_t = "t_read_length_#{hash[:cell]}-#{hash[:rep]}.txt"
  if !File.exists?("#{name_t}")
    File.open("#{name_t}", "w") do |out|
      hash[:t_trim].each do |hash2|
        out.write "#{hash2["Length"]}\t#{hash2["Count"]}\n"
      end
    end
  end

  name_m = "mcf_read_length_#{hash[:cell]}-#{hash[:rep]}.txt"
  if !File.exists?("#{name_m}")
    File.open("#{name_m}", "w") do |out|
      hash[:mcf_trim].each do |hash2|
        out.write "#{hash2["Length"]}\t#{hash2["Count"]}\n"
      end
    end
  end
end

## # # # # # # # # # # # # # # # # # #
## Run bowtie to align the reads to the genome
##

index = File.basename(opts.reference).split(".")[0..-2].join(".")
build = "bowtie-build #{opts.reference} #{index}"
if !File.exists?("#{index}.1.ebwt")
  puts build if opts.verbose
  `#{build}` if !opts.test
else
  puts "Index #{index} already exists"
end

threads = 22
input.each_with_index do |hash, i|
  sam = "#{hash[:cell]}-#{hash[:rep]}-k.sam"
  bowtie_cmd = "bowtie "
  bowtie_cmd += " --phred64-quals " # illumina 1.5 == phred64
  bowtie_cmd += " -v 0 " # exact matches only
  bowtie_cmd += " --best " # eliminates strand bias
  bowtie_cmd += " -k 50 " 
  bowtie_cmd += " -t -p #{threads} "
  bowtie_cmd += " #{index} "
  bowtie_cmd += " #{hash[:trimmed_t]}"
  bowtie_cmd += " #{sam}"
  input[i][:sam] = sam
  puts bowtie_cmd if opts.verbose
  if !File.exists?("#{sam}") 
    `#{bowtie_cmd}` if !opts.test
  else
    puts "bowtie already run for #{sam}"
  end
end

## # # # # # # # # # # # # # # # # # #
## get size of sam files
## how many reads aligned to the genome
##

puts "calculating sizes of sam files" if opts.verbose
if !File.exists?("sam_sizes-k.txt")
  File.open("sam_sizes-k.txt", "w") do |io|
    input.each do |hash|
      sam = hash[:sam]
      linecount = "wc -l #{sam}"
      puts linecount if opts.verbose
      if !opts.test
        output = `#{linecount}` 
        puts "output: #{output}" if opts.verbose
        l = output.split.first.to_i
        io.write "#{hash[:cell]}\t#{hash[:rep]}\t#{l}\n" 
      end
    end
  end
end

## # # # # # # # # # # # # # # # # # #
## sort the sam files
##
input.each_with_index do |hash, i|
  sorted_sam = "#{hash[:sam].split(".")[0..-2].join(".")}.sorted.sam"
  sort = "sort -k3,3 -k4,4n #{hash[:sam]} > #{sorted_sam}"
  if !File.exists?("#{sorted_sam}")
    puts sort if opts.verbose
    `#{sort}` if !opts.test
  else
    puts "#{hash[:sam]} already sorted"
  end
  input[i][:sorted_sam] = sorted_sam
end

## # # # # # # # # # # # # # # # # # #
## Find loci
##
sam_files = []
input.each do |hash|
  sam_files << hash[:sorted_sam]
end
find_loci_output = "srna_expression_strand.txt"
loci = "ruby /home/cmb211/scripts/sorghum/find_loci_strand.rb --input #{sam_files.join(",")} --output #{find_loci_output} "
loci << " --verbose " if opts.verbose
if !File.exists?("#{find_loci_output}")
  puts loci if opts.verbose
  `#{loci}` if !opts.test
else
  puts "find_loci_strand already run" if opts.verbose
end

abort "stop here to debug"

if !File.exists?("srna_locations.txt")
  cmd = "ruby loci.rb --srna srna_expression_strand.txt --annotation #{opts.annotation} --output srna_locations.txt"
  puts cmd if opts.verbose
  `#{cmd}` if !opts.test
else
  puts "loci.rb already run" if opts.verbose
end

## run bowtie again but with more numbers

# threads = 22
# input.each_with_index do |hash, i|
#   sam = "#{hash[:cell]}-#{hash[:rep]}-k.sam"
#   bowtie_cmd = "bowtie "
#   bowtie_cmd += " --phred64-quals " # illumina 1.5 == phred64
#   bowtie_cmd += " -v 0 " # exact matches only
#   bowtie_cmd += " --best " # eliminates strand bias
#   bowtie_cmd += " -k 50 " 
#   bowtie_cmd += " -t -p #{threads} "
#   bowtie_cmd += " #{index} "
#   bowtie_cmd += " #{hash[:trimmed_t]}"
#   bowtie_cmd += " #{sam}"
#   input[i][:sam] = sam
#   puts bowtie_cmd if opts.verbose
#   if !File.exists?("#{sam}") 
#     `#{bowtie_cmd}` if !opts.test
#   else
#     puts "bowtie already run for #{sam}"
#   end
# end
