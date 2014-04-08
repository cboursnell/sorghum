#!/usr/bin/env ruby

#
# find where srna loci are 
#

require 'trollop'

class Feature
  attr_accessor :name, :type, :features, :start, :stop

  def initialize(name, type,start,stop)
    @name = name
    @type = type
    @features = []
    @start = start
    @stop= stop
  end

  def contains?(start, stop)
    if start==@start and stop==@stop
      return 0
    elsif start > @start
      if stop > @stop
        if start < @stop
          return 1
        else
          return 5
        end
      else
        return 3
      end
    else # 2, 4, 6
      if stop > @stop
        return 4
      else
        if stop > @start
          return 2
        else
          return 6
        end
      end
    end
  end

  def to_s
    "#{name}\t#{type}\t#{start}\t#{stop}"
  end
end

def ticker count
  
  if count % 1_000_000 == 0
    print "|"
    print "\b"*9
  elsif count % 100_000 == 0
    print "*"
    print "\b"*9
  elsif count % 10_000 == 0
    print "."
  end
end

def loadGff(file)

  annotation = Hash.new
  mrna=nil
  line_count=0
  chromosome=""
  File.open(file).each_line do |line|
    line_count+=1
    ticker line_count #if opts.verbose
    cols = line.chomp.split("\t")
    chromosome = cols[0]
    feature = cols[2] # mrna, exon, cds, utr etc
    start = cols[3].to_i
    stop = cols[4].to_i
    #strand = cols[6] 
    desc = cols[8]
    if feature == "mRNA"
      if mrna # save previous mrna if it exists
        annotation[chromosome] = [] if !annotation.has_key?(chromosome)
        annotation[chromosome] << mrna
      end
      if desc.match(/Name=(.+?);/)       # Name=Sb01g000200.1;
        name = $1
      else
        name = "unknown"
      end
      mrna = Feature.new(name, "mrna", start,stop)
    elsif feature == "CDS"  or feature == "three_prime_UTR" or feature == "five_prime_UTR" # or feature == "exon"
      f = Feature.new(feature,feature, start, stop)
      mrna.features << f
    end
  end
  # puts "Done" #if opts.verbose
  annotation[chromosome] = [] if !annotation.has_key?(chromosome)
  annotation[chromosome] << mrna 
  annotation
end


### binary search function for finding where loci are in genome
def search(hash, chromosome, start, stop)
  # do a binary search for `start`
  list = hash[chromosome]
  found=false
  found2=false
  i = 0
  feature=nil
  mrna=nil
  p = 0
  k = 0
  k2 = 0
  if list
    j = list.length
    k = j
    while p != k
      p = k
      k = (i+j)/2
      o = list[k].contains?(start, stop)
      if o <= 4
        p = k
        found=true
        mrna=list[k].name
      elsif o == 5
        i = k
      elsif o == 6
        j = k
      end
    end
    if found
      list2 = list[k].features
      list2.sort! {|a,b| a.start <=> b.start }
      i=0
      j=list2.length
      while p != k2
        p = k2
        k2 = (i+j)/2
        if list2[k2]
          o = list2[k2].contains?(start, stop)
          if o<=2
            # puts "srna loci overlaps with feature"
            p = k2
            found2 = true
            feature = list2[k2].name
          elsif o==3
            # puts "srna loci is smaller than feature"
            p = k2
            found2 = true
            feature = list2[k2].name
          elsif o==4
            # puts "srna loci is bigger than feature"
            p = k2
            found2 = true
            feature = list2[k2].name
          elsif o==5
            i = k2
          elsif o==6
            j = k2
          end
        end
      end
    else
      # find nearest
      a = list[k] if list[k]
      b = list[k+1] if list[k+1]
      if a && b
        if start - a.stop < b.start - stop
          mrna = a.name
        else
          mrna = b.name
        end
      end
    end
  end
  feature = "na" if !feature
  mrna = "na" if !mrna
  return {:k => k, :found => found, :mrna => mrna, :found2 => found2, :feature => feature, :o => o}
end


if __FILE__ == $0

  opts = Trollop::options do
    version "v0.1"
    banner <<-EOS
    mrna pipeline

    author: Chris Boursnell (cmb211@cam.ac.uk)
    EOS
    opt :srna, "sRNA expression file. This is output by the pipeline", :required => true, :type => String
    opt :annotation, "gff file", :required => true, :type => String
    opt :output, "Output", :required => true, :type => String
    opt :verbose, "Be verbose"
    opt :test, "Don't actually run anything"
  end

  Trollop::die :srna, "must exist" if !File.exist?(opts[:srna]) if opts[:srna]
  Trollop::die :annotation, "must exist" if !File.exist?(opts[:annotation]) if opts[:annotation]

  gff = loadGff("#{opts.annotation}")
  puts "Done" if opts.verbose
  expression_cutoff=12
  count=0
  File.open("#{opts.output}", "w") do |out|
    File.open("#{opts.srna}").each_line do |line|
      count+=1
      ticker count if opts.verbose
      cols = line.chomp.split("\t")
      chromosome = cols[0]
      start = cols[1].to_i
      stop = cols[2].to_i
      length = stop - start + 1
      expression = cols[3..8].map {|x| x.to_i}
      if expression.inject(:+) >= expression_cutoff
        ans =  search(gff, chromosome, start, stop)
        out.write "#{chromosome}\t#{start}\t#{stop}\t#{length}\t#{ans[:found]}\t#{ans[:mrna]}\t#{ans[:found2]}\t#{ans[:feature]}\t#{expression.join("\t")}\n"
      end
    end
  end
  puts "Done" if opts.verbose

end