#!/usr/bin/env ruby
#
# title           : vcftobeacon.rb
# description     : The following script takes VCF files and generates the three files with the columns for the Beacon.
# author          : Toshiaki Katayama, Dietmar Fernández Orth
# date            : 2019-12-18
# version         : 1.0
# usage           : ruby vcftobeacon.rb /path/to/*.vcf
# extended_usage  : BeaconDatasetID=5 BeaconDataSeparator=";" BeaconDataHeader=true ruby vcftobeacon.rb /path/to/*.vcf > vcftobeacon.log
# notes           : Install bcftools and add plugin +fill-tags to the path (try `bcftools +fill-tags -vv` to see if you have the plugin).
# ruby_version    : 2.6
#
# NOTE:
#  * TODO: confirm if variants.matching.sample and samples tables works
#  * Don't need to print headers as I found CSV headers in sample data are dirfferent from the schema (and also from the ones in Dietmar's script) ...
#  * Couldn't make the first autoincrement id column to null in data; if it's possible, we could use the following shortcut, hopefully:
#    `cat file | psql -h localhost -p 5432 -d elixir_beacon_dev -U microaccounts_dev -c "copy beacon_data_table from stdin with (format 'text')"`
#

# To change the dataset ID, set the environmental variable BeaconDatasetID (deafult: 1)
dataset_id = ENV['BeaconDatasetID'] || 1

# To change a separator for columns, set the environmental variable BeaconDataSeparator, e.g., change to ";" for SQL (default: "\t")
separator = ENV['BeaconDataSeparator'] || "\t"

# To enable column headers to be printed, set the environmental variable BeaconDataHeader to any value, e.g., "true" (default: false)
header = ENV['BeaconDataHeader'] || false

# Define primary columns to be extracted from VCF files by the bcftools
# * AC: allele count in genotypes, for each ALT allele, in the same order as listed
# * AN: total number of alleles in called genotypes
# * NS: Number of samples with data
# * AF: allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
columns = %w(CHROM POS ID REF ALT END AC AN NS AF SVLEN SVTYPE TYPE)

require 'date'

class PrintBuffer
  def initialize(io, lines = 10000)
    @io = io
    @buffer = []
    @counter = 0
    @lines = lines
  end

  def push(line)
    @buffer << line
    @counter += 1
    if @counter % @lines == 0
      $stderr.puts "... flush #{@lines} lines ..."
      flush
    end
  end

  def flush
    @io.puts @buffer unless @buffer.empty?
    @buffer = []
  end
end

puts "The following script takes VCF files and generates the three files needed for Beacon."
puts "---------------"

bcfquery = columns.map { |col| "%#{col}"}.join("\t")
colidx_svlen = columns.index("SVLEN")
colidx_svtype = columns.index("SVTYPE")
colidx_type = columns.index("TYPE")

ARGV.each_with_index do |vcf, i|
  count = "#{i+1} / #{ARGV.size} files"
  puts "#{DateTime.now.to_s} START #{count}"

  puts "Normalizing file #{vcf}"
  system("bcftools norm -m -both #{vcf} -o #{vcf}.norm")

  puts "Generating file #{vcf}.variants.data"
  # https://github.com/ga4gh-beacon/beacon-elixir/blob/master/deploy/db/db/data/1_chr21_subset.variants.csv
  # datasetId;chromosome;position;variantId;reference;alternate;end;svType;svLength;variantCount;callCount;sampleCount;frequency;sampleMatchingCount
  # 1;21;9411238;rs559462325;G;A;;SNP;;1;5008;2504;0.000199681;1
  # 1;21;9411244;rs181691356;C;A;;SNP;;4;5008;2504;0.000798722;4
  variants_file = File.open("#{vcf}.variants.data", "w")
  variants_buffer = PrintBuffer.new(variants_file)
  variants_header = %w(datasetId chromosome position variantId reference alternate end svType svLength variantCount callCount sampleCount frequency sampleMatchingCount)
  variants_file.puts variants_header.join(separator) if header
  variants_columns = %w(CHROM POS ID REF ALT END TYPE SVLEN AC AN NS AF)
  variants_slices = variants_columns.map { |col| columns.index(col) }.compact

  puts "Generating file #{vcf}.variants.matching.sample.data"
  # https://github.com/ga4gh-beacon/beacon-elixir/blob/master/deploy/db/db/data/1_chr21_subset.variants.matching.samples.csv
  # datasetId;chromosome;position;variantId;reference;alternate;svType;sampleId
  # 1;21;9411238;rs559462325;G;A;SNP;{HG01029}
  # 1;21;9411244;rs181691356;C;A;SNP;{HG01104,NA11995,NA19922,NA20807}
  var_sam_file = File.open("#{vcf}.variants.matching.samples.data", "w")
  var_sam_buffer = PrintBuffer.new(var_sam_file)
  var_sam_header = %w(datasetId chromosome position variantId reference alternate svType sampleId)
  var_sam_file.puts var_sam_header.join(separator) if header
  var_sam_columns = %w(CHROM POS ID REF ALT TYPE)
  var_sam_slices = var_sam_columns.map { |col| columns.index(col) }.compact

  puts "#{DateTime.now.to_s} start generation of two data files"
  cmd = "bcftools query -f '#{bcfquery}[\t%SAMPLE=%GT]\n' #{vcf}.norm"
  IO.popen(cmd, "r").each do |line|
    ary = line.strip.split("\t")
    # Assign null value for the SVLEN (sv_length) column if bcftools puts '.' there as the column is typed as integer
    if ary[colidx_svlen] == '.'
      ary[colidx_svlen] = ''
    end
    # If there is a value in "SVTYPE" other than ".", replace "TYPE" column value with it
    if ary[colidx_svtype] != '.'
      ary[colidx_type] = ary[colidx_svtype]
    end
    gts = ary[(columns.size)..-1]
    hetero = gts.select { |gt| gt[/=(0\|1|1\|0|0\/1|1\/0)/] }
    homalt = gts.select { |gt| gt[/=(1\|1|1\/1)/] }
    samples = (hetero + homalt).map { |gt| gt.sub(/=.*/, '') }
    num_hetero = hetero.size
    num_homalt = homalt.size
    sum_hethom = samples.size

    variants_buffer.push [ dataset_id, ary.values_at(*variants_slices), sum_hethom ].flatten.join(separator)
    var_sam_buffer.push  [ dataset_id, ary.values_at(*var_sam_slices), "{#{samples.join(",")}}" ].flatten.join(separator)
  end
  variants_buffer.flush
  var_sam_buffer.flush
  puts "#{DateTime.now.to_s} end generation of two data files"

  puts "Generating file #{vcf}.samples.data"
  # https://github.com/ga4gh-beacon/beacon-elixir/blob/master/deploy/db/db/data/1_chr21_subset.samples.csv
  # sampleId;datasetId
  # HG00096;1
  # HG00097;1
  samples_file = File.open("#{vcf}.samples.data", "w")
  sample_header = %w(sampleId datasetId)
  samples_file.puts sample_header.join(separator) if header
  cmd = "bcftools query -l #{vcf}"
  IO.popen(cmd, "r").each do |sample_id|
    samples_file.puts [ sample_id.strip, dataset_id ].join(separator)
  end

  variants_file.close
  var_sam_file.close
  samples_file.close

  puts "#{DateTime.now.to_s} DONE #{count}"
  puts
end
