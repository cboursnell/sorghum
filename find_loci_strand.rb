#!/usr/bin/env ruby

#
# scan through sam files and create a list of srna loci from where
# the reads have mapped. to be run by the sorghum srna pipeline
# script.
#
# author: chris boursnell (cmb211@cam.ac.uk)
# date: 2014-03-18
#

require 'rubygems'
require 'trollop'

opts = Trollop::options do
  version "v0.1"
  banner <<-EOS
  find loci stranded

  author: Chris Boursnell (cmb211@cam.ac.uk)
  EOS
  opt :input, "Comma separated list of sorted sam files", :required => true, :type => String
  opt :output, "Output file", :required => true, :type => String
  opt :verbose, "Be verbose"
end

opts.input.split(",").each do |file|
  abort "#{file} must exist" if !File.exists?(file)
end

def ticker(i,speed)
  n = 10**speed
  if i <= 1
    print " "*9
    print "0"
  end
  if i % n == 0
    print "\b"*10
    string = "#{i}"
    print " "*(10-string.length)
    print string
  end
end

class Locus
  attr_accessor :chromosome, :start, :stop, :pvecount, :nvecount

  def initialize(chromosome, start, stop, size)
    @chromosome = chromosome
    @start = start.to_i  # derived from the mRNA line in the gff file
    @stop = stop.to_i
    @pvecount = Array.new(size).map {|x| x=0}
    @nvecount = Array.new(size).map {|x| x=0}
  end

  def overlaps(other)
    state=0
    if other.start > @start # 1, 3, 5
      if other.stop > @stop #1, 5
        if other.start < @stop
          state=1
        else
          state=5
        end
      else # 3
        state=3
      end
    else # 2,4,6
      if other.stop < @stop # 2,6
        if other.stop > @start # 2
          state=2
        else
          state=6
        end
      else
        state=4
      end
    end
    state
  end

  def contains(start)
    if start.to_i >= @start and start.to_i <= @stop
      return true
    else
      return false
    end
  end

  def to_sym
    "#{@chromosome}:#{@start}".to_sym
  end

  def to_s
    "#{@chromosome}\t#{@start}..#{@stop}"
  end
end

##
# merge the sam files together into one large sorted sam file
#

list = opts.input.split(",")
merged = "merged_stranded.sam"
merge = "sort -m -k3,3 -k4,4n #{list.join(" ")} > #{merged}"
if !File.exists?("#{merged}")
  puts "merging sorted sam files"
  `#{merge}`
else
  puts "#{merged} already exists" if opts.verbose
end


current_locus = Locus.new("-1", 0, 1,1)
loci = Hash.new
line_count=1
print "calculating loci..." if opts.verbose
File.open("#{merged}").each_line do |line|
  ticker(line_count,5)
  line_count += 1
  cols = line.chomp.split("\t")
  chromosome = cols[2]      #  start    stop                        number of files
  read = Locus.new(chromosome, cols[3], cols[3].to_i+cols[4].length, list.length)
  state = current_locus.overlaps(read)
  if state <= 4 # overlaps
    # if the read in the sam file overlaps the current locus then extend the locus
    # and add 1 to the count
    if state==1
      current_locus.stop = read.stop
    elsif state==2
      current_locus.start = read.start
    elsif state==4
      current_locus.start = read.start
      current_locus.stop = read.stop
    end
  else 
    # if the read doesn't overlap the current locus then create a new current locus
    # based on the read and set the count to 1
    current_locus = read
    loci[chromosome] = [] if !loci.has_key?(chromosome)
    loci[chromosome] << current_locus
  end
end
puts "Done" if opts.verbose


### load the srna loci from the srna_expression.txt file that
### has already been created

# loci = Hash.new
# list = opts.input.split(",")
# size = list.length
# count=1
# print "scanning srna expression" if opts.verbose
# File.open("srna_expression_head.txt").each_line do |line|
#   ticker(count,3)
#   count+=1
#   cols = line.chomp.split("\t")
#   chromosome = cols[0]
#   start = cols[1]
#   stop = cols[2]
#   loci[chromosome] = [] if !loci.has_key?(chromosome)
#   loci[chromosome] << Locus.new(chromosome, start, stop, size) #   def initialize(chromosome, start, stop, size)
# end
# puts "Done" if opts.verbose


### write the loci has to a file at this point because
### it's taken so long to generate it

File.open("srna_loci.txt", "w") do |out|
  loci.each_pair do |chrom, list|
    list.each do |locus|
      out.write "#{locus}\n"
    end
  end
end

###

loci_list=[]
list.each_with_index do |sam, sample_number|
  index=0 # current position in the loci list
  print "opening sam file #{sam}..." if opts.verbose
  count = 1
  current_chromosome = ""
  File.open("#{sam}").each_line.with_index do |line, i|
    ticker(count,5)
    count+=1
    cols = line.chomp.split("\t")
    strand = cols[1]
    chrom = cols[2]
    pos = cols[3]
    puts "line = #{i}, strand = #{strand}, chrom = #{chrom}, pos = #{pos}"
    if chrom != current_chromosome
      # load new list
      # puts "previous \"#{current_chromosome}\" new #{chrom}" if opts.verbose
      loci_list = loci[chrom]
      current_chromosome = chrom
      index = 0
      puts "length of loci_list = #{loci_list.length}"
    end
    if index < loci_list.size
      while index < loci_list.size and !loci_list[index].contains(pos)
        index += 1
        # puts "index = #{index}, pos = #{pos}, loci_list[index] = #{loci_list[index].start}..#{loci_list[index].stop}"
      end
      if index >= loci_list.size
        puts "ERROR 1: sample=#{sample_number}\tindex=#{index}\tline_in_sam=#{count}" if opts.verbose
      else
        if strand == "+"
          loci[chrom][index].pvecount[sample_number] +=1
        else
          loci[chrom][index].nvecount[sample_number] +=1
        end
      end
    else
      puts "ERROR 2: sample=#{sample_number}\tindex=#{index}\tline_in_sam=#{count}" if opts.verbose
    end
  end
  puts "Done" if opts.verbose
end

puts "Writing output #{opts.output}"
File.open("#{opts.output}", "w") do |out|
  loci.each_pair do |chromosome, locus_list|
    locus_list.each do |locus|
      out.write "#{locus.chromosome}\t#{locus.start}\t#{locus.stop}"
      locus.pvecount.each do |x|
        out.write "\t#{x}"
      end
      locus.nvecount.each do |x|
        out.write "\t#{x}"
      end
      out.write "\n"
    end
  end
end
puts "Done"
