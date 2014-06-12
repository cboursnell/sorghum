#!/usr/bin/env ruby

#
# given a genome and some mRNA seq reads make a transcriptome
# and gtf file
#

require 'rubygems'
require 'trollop'

opts = Trollop::options do
  version "v0.1"
  banner <<-EOS
  mrna pipeline

  author: Chris Boursnell (cmb211@cam.ac.uk)
  EOS
  opt :input, "Input file of fastq files, alternating left then right", :type => String
  opt :left, "First fastq files (comma separated list)", :type => String
  opt :right, "Second fastq files (comma separated list)", :type => String
  opt :genome, "Reference fasta file", :required => true, :type => String
  opt :existing_gff, "Validate an existing gff file", :type => String
  opt :outputdir, "Output directory", :required => true, :type => String
  opt :threads, "Number of threads", :default => 1, :type => :int
  opt :verbose, "Be verbose"
  opt :test, "Don't actually do anything"
end

Trollop::die :genome, "must exist" if !File.exist?(opts[:genome]) if opts[:genome] 

# # # # # # # # # # # # # # #
# run tophat

# tophat_path = "/home/cmb211/apps/tophat/tophat2"
tophat_path = `which tophat2`.chomp
abort "Can't find tophat2. Please make sure it is in your PATH" if tophat_path==""

left=[]
right=[]
if opts.input
  count=0
  File.open("#{opts.input}").each_line do |line|
    line.chomp!
    if count%2==0
      left << line
    else
      right << line
    end
    count += 1
  end
else
  left = opts.left.split(",")
  right = opts.right.split(",")
end
if left.length != right.length
  abort "Please give equal numbers of left and right reads"
end

(left+right).each do |file|
  abort "#{file} must exist" if !File.exists?(file)
end

if !Dir.exists?("#{opts.outputdir}")
  mkdir = "mkdir #{opts.outputdir}"
  `#{mkdir}` if !opts.test
end

## # # # # # # # # # # # # # # # # # # # #
## build bowtie index
##
index = "#{File.basename(opts.genome).split(".")[0..-2].join(".")}" 
if !File.exists?("#{index}.1.bt2")
  build = "/home/cmb211/apps/bowtie2/bowtie2-build #{opts.genome} #{index}"
  puts build if opts.verbose
  `#{build}` if !opts.test
end

## # # # # # # # # # # # # # # # # # # # #
## construct tophat cmd to align reads
## and find splice junctions
##
tophat_cmd = "#{tophat_path}"
tophat_cmd += " -o #{opts.outputdir} " # options
tophat_cmd += " -p #{opts.threads} " # options
# tophat_cmd += " --mate-inner-dist 50 " # default:50, distance between end of 
tophat_cmd += " --phred64-quals "  #
tophat_cmd += " --tmp-dir /tmp/cmb211 " if `hostname`.split(".").first == "node8"
tophat_cmd += " --tmp-dir /disk2/tmp/cmb211 " if `hostname`.split(".").first == "node9"
tophat_cmd += " --no-convert-bam "  # Output is <output_dir>/accepted_hit.sam)
tophat_cmd += " --b2-very-sensitive "
tophat_cmd += " -G #{opts.existing_gff} "
tophat_cmd += " #{index} "
tophat_cmd += " #{left.join(",")} "
tophat_cmd += " #{right.join(",")} "

if !File.exists?("#{opts.outputdir}/accepted_hits.sam")  # accepted_hits.sam
  puts tophat_cmd if opts.verbose
  `#{tophat_cmd}` if !opts.test
else
  puts "Tophat2 already run"
end

## # # # # # # # # # # # # # # # # # # # #
## run cufflinks to take sam file
## and produce new output gtf file
##

cufflinks_path = `which cufflinks`.chomp
abort "Can't find cufflinks. Please make sure it is in your PATH" if cufflinks_path==""

cufflinks_cmd = "#{cufflinks_path}"
cufflinks_cmd += " -o #{opts.outputdir} " 
cufflinks_cmd += " -p #{opts.threads} "
cufflinks_cmd += " #{opts.outputdir}/accepted_hits.sam "

puts cufflinks_cmd if opts.verbose
if !File.exists?("#{opts.outputdir}/transcripts.gtf")
  `#{cufflinks_cmd}` if !opts.test
else
  puts "Cufflinks already run"
end

# Extracting transcript sequences
# The gffread utility can be used to generate a FASTA file with the DNA 
# sequences for all transcripts in a GFF file. For this operation a fasta
#  file with the genomic sequences have to be provided as well. For example,
#   one might want to extract the sequence of all transfrags assembled from 
#   a Cufflinks assembly session. This can be accomplished with a command
#    line like this:

gffread_path = `which gffread`.chomp
abort "Can't find gffread. Please make sure it is in your PATH" if gffread_path==""

gffcmd = "#{gffread_path} -w #{opts.outputdir}/#{index}-transcripts.fa -g #{opts.genome} #{opts.outputdir}/transcripts.gtf"

puts gffcmd if opts.verbose
`#{gffcmd}` if !opts.test
