#!/usr/bin/env ruby

# title           : vcftobeacon.rb
# description     : The following script takes the vcf files and generates a file with the columns for the Beacon.
# author          : Toshiaki Katayama
# date            : 2019-12-17
# version         : 1.0
# usage           : ruby vcftobeacon.rb /path/to/*.vcf
# notes           : Install bcftools and add plugin +fill-tags to the path.
# ruby_version    : 2.6

require 'date'

columns = %w(CHROM POS ID REF ALT END AC AN NS AF SVLEN SVTYPE TYPE)
query = columns.map { |col| "%#{col}"}.join("\t")

puts "The following script takes the vcf files and generates the three files needed for Beacon."
puts "---------------"

ARGV.each_with_index do |vcf, i|
  count = "#{i+1} / #{ARGV.size}"
  puts "#{DateTime.now.to_s} START #{count}"
  puts "Normalizing file #{vcf}"
  `bcftools norm -m -both #{vcf} -o #{vcf}.norm`

  puts "Generating file #{vcf}.variants.tsv"
  # [1]CHROM        [2]POS  [3]ID   [4]REF  [5]ALT  [6]END  TYPE    [7]SVLEN        [8]AC   [9]AN   [10]NS  [11]AF  SUM
  # 21      9411239 rs559462325     G       A       9411239 SNP     .       1       5008    2504    0.000199681     1
  # 21      9411245 rs181691356     C       A       9411245 SNP     .       4       5008    2504    0.000798722     4
  variants_columns = %w(CHROM POS ID REF ALT END TYPE SVLEN AC AN NS AF SUM)
  variants_slices = variants_columns.map { |col| columns.index(col) }.compact
  variants_tsv = File.open("#{vcf}.variants.tsv", "w")
  variants_tsv.puts variants_columns.join("\t")

  puts "Generating file #{vcf}.variants.matching.sample.tsv"
  # CHR POS ID REF ALT      TYPE    sampleId
  # 21 9411239 rs559462325 G A      SNP     HG01029,
  # 21 9411245 rs181691356 C A      SNP     HG01104,NA11995,NA19922,NA20807,
  var_sam_columns = %w(CHROM POS ID REF ALT TYPE sampleId)
  var_sam_slices = var_sam_columns.map { |col| columns.index(col) }.compact
  var_sam_tsv = File.open("#{vcf}.variants.matching.sample.tsv", "w")
  var_sam_tsv.puts var_sam_columns.join("\t")

  puts "#{DateTime.now.to_s} bcftools start"
  `bcftools query -f '#{query}[\t%SAMPLE=%GT]\n' #{vcf}.norm -o #{vcf}.norm.tsv`
  puts "#{DateTime.now.to_s} bcftools end"
  File.open("#{vcf}.norm.tsv", "r").each do |line|
    ary = line.strip.split("\t")
    gts = ary[(columns.size)..-1]
    hetero = gts.select { |gt| gt[/=(0\|1|1\|0|0\/1|1\/0)/] }
    homalt = gts.select { |gt| gt[/=(1\|1|1\/1)/] }
    samples = (hetero + homalt).map { |gt| gt.sub(/=.*/, '') }
    num_hetero = hetero.size
    num_homalt = homalt.size
    sum_hethom = samples.size

    variants_tsv.puts [ ary.values_at(*variants_slices), sum_hethom ].flatten.join("\t")
    var_sam_tsv.puts  [ ary.values_at(*var_sam_slices), samples.join(",") ].flatten.join("\t")
  end

  puts "Generating file #{vcf}.samples.tsv"
  # sampleId
  # HG00096
  # HG00097
  # HG00099
  samples_tsv = File.open("#{vcf}.samples.tsv", "w")
  samples_tsv.puts "sampleId"
  samples = `bcftools query -l #{vcf}`
  samples_tsv.puts samples

  variants_tsv.close
  var_sam_tsv.close
  samples_tsv.close

  puts "#{DateTime.now.to_s} DONE #{count}"
  puts
end
