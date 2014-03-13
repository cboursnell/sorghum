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
  opt :input, "Input file describing fastq files", :type => String
  opt :reference, "Reference fasta file"
  opt :morehelp, "More help"
  opt :verbose, "Be verbose"
end

Trollop::die :input, "must exist" if !File.exist?(opts[:input]) if (opts[:input] or !opts[:morehelp])

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
    `#{minion}`
  end
  File.open("adapter_list.txt", "w") {|io| io.write(adapter_list.join("\n"))}
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
  `#{fastqc}`
end

input.each_with_index do |hash,index|
  if Dir.exists?("fastqc_output/#{File.basename(hash[:file])}_fastqc")
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
      io.write "#{hash[:count]}\t#{hash[:cell]}\t#{hash[:rep]}\n"
    end
  end
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
          adapter_out.write ">adapter#{count}\n"
          adapter_out.write "#{$1}\n"
          count+=1
        end
      end
      hash[:over].each do |seq|
        adapter_out.write ">adapter#{count}\n"
        adapter_out.write "#{seq}\n"
        count+=1
      end
    end
  end
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
  trim_jar = "/applications/trimmomatic/Trimmomatic-0.30/trimmomatic-0.30.jar"
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
  input.each do |hash|
    infile = hash[:file]
    outfile = "t.#{File.basename(infile)}"
    trim_cmd = "java -jar #{trim_jar} SE #{phred} "                     # phred
    # trim_cmd += " -trimlog trim_log_#{hash[:cell]}#{hash[:rep]}.txt "   # log output
    trim_cmd += " -threads #{threads} #{infile} #{outfile} "            # threads
    trim_cmd += " ILLUMINACLIP:#{adapter_file}:#{seed_mismatches}:#{palindrome_clip_threshold}:#{simple_clip} " 
    trim_cmd += " LEADING:#{leading} TRAILING:#{trailing} SLIDINGWINDOW:#{windowsize}:#{quality} MINLEN:#{minlen} "
    trimmed << outfile
    if !File.exists?("#{outfile}")
      puts trim_cmd
      `#{trim_cmd}`
    end
  end
end

# # FASTQ MCF # # # # # # # # # #

puts "Starting fastq-mcf" if opts.verbose
input.each do |hash|
  mcf = "fastq-mcf "
  # options
  mcf += " -o mcf.#{File.basename(hash[:file])} "
  mcf += " #{adapter_file}"
  mcf += " #{hash[:file]} "
  trimmed << "mcf.#{File.basename(hash[:file])}"
  if !File.exists?("mcf.#{File.basename(hash[:file])}")
    puts mcf
    `#{mcf}` 
  end
end
# usage: fastq-mcf [options] <adapters.fa> <reads.fq> [mates1.fq ...] 

# Detects levels of adapter presence, computes likelihoods and
# locations (start, end) of the adapters.   Removes the adapter
# sequences from the fastq file(s).

# Stats go to stderr, unless -o is specified.

# Specify -0 to turn off all default settings

# If you specify multiple 'paired-end' inputs, then a -o option is
# required for each.  IE: -o read1.clip.q -o read2.clip.fq

# Options:
#     -h       This help
#     -o FIL   Output file (stats to stdout)
#     -s N.N   Log scale for adapter minimum-length-match (2.2)
#     -t N     % occurance threshold before adapter clipping (0.25)
#     -m N     Minimum clip length, overrides scaled auto (1)
#     -p N     Maximum adapter difference percentage (10)
#     -l N     Minimum remaining sequence length (19)
#     -L N     Maximum remaining sequence length (none)
#     -D N     Remove duplicate reads : Read_1 has an identical N bases (0)
#     -k N     sKew percentage-less-than causing cycle removal (2)
#     -x N     'N' (Bad read) percentage causing cycle removal (20)
#     -q N     quality threshold causing base removal (10)
#     -w N     window-size for quality trimming (1)
#     -0       Set all default parameters to zero/do nothing
#     -U|u     Force disable/enable Illumina PF filtering (auto)
#     -P N     Phred-scale (auto)
#     -R       Don't remove N's from the fronts/ends of reads
#     -n       Don't clip, just output what would be done
#     -C N     Number of reads to use for subsampling (300k)
#     -S       Save all discarded reads to '.skip' files
#     -d       Output lots of random debugging stuff

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
  `mkdir trimmed_fastqc_output`
  `#{fastqc}`
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
File.open("t_read_counts.txt", "w") {|io| io.write(t_counts)}
File.open("mcf_read_counts.txt", "w") {|io| io.write(mcf_counts)}

## # # # # # # # # # # # # # # # # # #
## Read lengths after trimming
##


input.each_with_index do |hash, index|
  # puts "#{hash[:cell]}\t#{hash[:rep]}"
  File.open("t_read_length_#{hash[:cell]}-#{hash[:rep]}.txt", "w") do |out|
    hash[:t_trim].each do |hash2|
      out.write "#{hash2["Length"]}\t#{hash2["Count"]}\n"
    end
  end
  File.open("mcf_read_length_#{hash[:cell]}-#{hash[:rep]}.txt", "w") do |out|
    hash[:mcf_trim].each do |hash2|
      out.write "#{hash2["Length"]}\t#{hash2["Count"]}\n"
    end
  end
end

## # # # # # # # # # # # # # # # # # #
## Run bowtie2 to align the reads to the genome
##