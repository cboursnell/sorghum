#!/usr/bin/env ruby

#
# having run segmentSeq to get the start and stop of
# loci and bayseq to get the DE likelihoods for all
# the loci we need to find the strandedness of each
# loci.
#
# The output file from segmentSeq/baySeq looks like this
#
# 0            1          2       3    4  5       6        7      8
# seqnames     start     end   width  strand  BS.1 BS.2   BS.3   M.1    /
# chromosome_1 5876608   5876680   73   * 14885   9534    9166   580    /
# chromosome_1 7426440   7426595   156  * 1079962 1162923 976139 289485 /
# chromosome_1 8138794   8138854   61   * 14647   10024   9328   3916   /
# chromosome_1 11365009  11365091  83   * 9721    6618    5955   654    /

# 9       10      11         12      13
# M.2    M.3     Likelihood  DE   FDR.DE
# 819    1174    1           1>2     0
# 317514 296922  1           1>2     0
# 3570   3725    1           1>2     0
# 842    959     1           1>2     0
#
# The input file looks like this
#
# 0             1                         2     3      4     5
# chr           tag                       count start  end   strand
# chromosome_1  ACCCTAAACCCTAAACCCTAAACC  1     139    162   -
# chromosome_1  CGAATCTTTATACGCATGCATAG   2     612    634   -
# chromosome_1  TTACGAGACGAATCTTTTGAGCAT  1     864    887   +
# chromosome_1  GAGACGAATCTTTTGAGCATAGTT  1     868    891   -
# chromosome_1  AATTGTCTAGCAGCAGGAGTTGGA  4     999    1022  +
# chromosome_1  GTTCCATACAAGTGCCACAAGATT  2     1426   1449  -
# chromosome_1  AGTGCCACAAGATTCGATGTGAT   1     1436   1458  -
# chromosome_1  CAAGATTCGATGTGATGGGAA     1     1443   1463  -
#
# The goal is to place these 'tags' inside the loci and record the count
# on each strand
#

class Loci

  attr_accessor :chromosome, :start, :stop, :width, :strand, :counts
  attr_accessor :likelihood, :de, :fdr

  def initialize line
    cols = line.chomp.split("\s")
    @chromosome = cols[0]
    @start = cols[1].to_i
    @stop = cols[2].to_i
    @width = cols[3].to_i
    @strand = cols[4]
    @counts = cols[5..10].map {|x| x.to_i}
    @likelihood = cols[11]
    @de = cols[12]
    @fdr = cols[13]
    @positive = []
    @negative = []
  end

  def contains? tag
    return true if tag.start >= @start and tag.stop <= @stop and
                      tag.chromosome == @chromosome
    return false
  end

  def add_count count, strand, i
    if strand=="+"
      @positive[i] ||= 0
      @positive[i] += count
    elsif strand=="-"
      @negative[i] ||= 0
      @negative[i] += count
    end
  end

  def to_s
    s = "#{@chromosome}\t#{@start}\t#{@stop}\t#{@de}\t"
    @counts.each do |c|
      s << "#{c}\t"
    end
    s << "#{@likelihood}\t#{@fdr}\t"
    positive = @positive.compact.reduce(0) {|sum,x| sum += x if !x.nil?}
    negative = @negative.compact.reduce(0) {|sum,x| sum += x if !x.nil?}
    pp = positive.to_f/(positive+negative)
    np = negative.to_f/(positive+negative)
    s << "#{pp}\t#{np}\n"
  end

  def to_long_s
    s = "#{@chromosome}\t#{@start}\t#{@stop}\t#{@de}\t"
    @counts.each do |c|
      s << "#{c}\t"
    end
    @positive.each do |p|
      s << "#{p}\t"
    end
    @negative.each do |n|
      s << "#{n}\t"
    end
    s << "1\n"
  end
end

class Tag

  attr_accessor :chromosome, :tag, :count, :start, :stop, :strand

  def initialize line
    cols = line.chomp.split("\s")
    @chromosome = cols[0]
    @tag = cols[1]
    @count = cols[2].to_i
    @start = cols[3].to_i
    @stop = cols[4].to_i
    @strand = cols[5]
  end
end

class Finder

  # file_list: array of filenames containing counts for each rep
  # loci_file: filename of bayseq output
  attr_reader :loci

  def initialize file_list, loci_file
    @files = file_list
    @loci_file = loci_file
    @i = 0
    @c = ""
    @loci = []
  end

  def parse_loci
    path = File.expand_path(@loci_file)
    File.open("#{path}").each_line do |line|
      unless line =~ /seqnames/
        @loci << Loci.new(line)
      end
    end
  end

  def parse_tags file, i
    @i = 0
    path = File.expand_path(file)
    File.open("#{path}").each_line do |line|
      unless line =~ /tag/
        tag = Tag.new(line)
        self.add_tag_to_loci tag, i
      end
    end
  end

  def add_tag_to_loci tag, i
    if @i < @loci.length and tag.stop < @loci[@i].start
      # do nothing
    elsif @i < @loci.length and @loci[@i].contains? tag
      # add tag to @loci
      @loci[@i].add_count tag.count, tag.strand, i
    elsif @i < @loci.length and tag.start > @loci[@i].stop
      while @i < @loci.length and tag.start > @loci[@i].stop and tag.chromosome == @loci[@i].chromosome
        @i += 1
      end
      if @i < @loci.length and @loci[@i].contains? tag
        # add tag to @loci
        @loci[@i].add_count tag.count, tag.strand, i
      end
    end
  end

  def output
    s = ""
    @loci.each do |l|
      s << l.to_s
    end
    s
  end

end

list = []
list << "~/sorghum/srna_new/BS-1.sort.txt"
list << "~/sorghum/srna_new/BS-2.sort.txt"
list << "~/sorghum/srna_new/BS-3.sort.txt"
list << "~/sorghum/srna_new/M-1.sort.txt"
list << "~/sorghum/srna_new/M-2.sort.txt"
list << "~/sorghum/srna_new/M-3.sort.txt"
loci_file = "~/sorghum/srna_new/sorghum_DE_srna.sort.txt"
# loci_file = "/home/cmb211/sorghum/srna_new/sorghum_DE_srna.txt"
# loci_file = "/home/cmb211/scripts/sorghum/sorghum_DE_srna.head3.txt"

finder = Finder.new(list, loci_file)

finder.parse_loci

list.each_with_index do |file,i|
  puts file
  finder.parse_tags file, i
end

File.open("/home/cmb211/sorghum/srna_new/sorghum_DE_srna_stranded.txt", "w") do |f|
  f.write finder.output
end