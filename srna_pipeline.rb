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
  opt :morehelp, "More help"
end

Trollop::die :input, "must exist" if !File.exist?(opts[:input]) if (opts[:input] or !opts[:morehelp])

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

## # # # # # # # # # # # # # # # # # #
## Open input file and store contents

input = []
File.open("#{opts.input}").each_line do |line|
  cols = line.chomp.split(",")
  input << {:file => cols[0], :rep => cols[1], :cell => cols[2]}
end

p input

## # # # # # # # # # # # # # # # # # #
## Find number of reads in each file
## and export as a table

if !File.exists?("raw_file_length.txt")
  
end