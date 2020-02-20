#!/usr/bin/env ruby
#
# repository      : https://github.com/ktym/vcftobeacon
# title           : vcftobeacon.rb
# description     : The following script takes VCF files and generates the three files with the columns for the Beacon.
# author          : Toshiaki Katayama, Dietmar Fernandez-Orth
# date            : 2020-02-20
# version         : 1.4
# usage           : ruby vcftobeacon.rb /path/to/*.vcf
# usage_with_opts : ruby vcftobeacon.rb -i 5 -s ';' -c /path/to/conv_file.txt -h /path/to/*.vcf
# usage_with_vars : BeaconDatasetID=5 BeaconDataSeparator=";" BeaconChrsFile=/path/to/conv_file.txt BeaconDataHeader=true ruby vcftobeacon.rb /path/to/*.vcf
#                   (Note that the command line options override the environmental/shell variables)
# notes           : Install bcftools and add plugin +fill-tags to the path (try `bcftools +fill-tags -vv` to see if you have the plugin).
# ruby_version    : Should work with any version of Ruby from 1.8 through 2.7
#

require 'getoptlong'
require 'date'

class VcfToBeacon

  def parse_options
    @opts = {
      :dataset_id => ENV['BeaconDatasetID'] || 1,
      :separator => ENV['BeaconDataSeparator'] || "\t",
      :header => ENV['BeaconDataHeader'] || false,
      :chrs_file => ENV['BeaconChrsFile'] || './chr_name_conv.txt',
    }
    args = GetoptLong.new(
      [ '--datasetid',  '-i',  GetoptLong::REQUIRED_ARGUMENT ],
      [ '--separator',  '-s',  GetoptLong::REQUIRED_ARGUMENT ],
      [ '--header',     '-p',  GetoptLong::NO_ARGUMENT ],
      [ '--chrsfile',   '-c',  GetoptLong::REQUIRED_ARGUMENT ],
    )
    args.each_option do |name, value|
      case name
      when /--datasetid/
        @opts[:dataset_id] = value
      when /--separator/
        @opts[:separator] = value
      when /--header/
        @opts[:header] = true
      when /--chrsfile/
        @opts[:chrs_file] = value
      end
    end
  end

  def initialize(files)
    parse_options
    setup_columns
    convert_vcf(files)
  end

  def setup_columns
    # Define primary columns to be extracted from VCF files by the bcftools
    # * AC: allele count in genotypes, for each ALT allele, in the same order as listed
    # * AN: total number of alleles in called genotypes
    # * NS: Number of samples with data
    # * AF: allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
    columns = %w(CHROM POS ID REF ALT END AC AN NS AF SVLEN SVTYPE TYPE)
    variants_columns = %w(CHROM POS ID REF ALT END TYPE SVLEN AC AN NS AF)
    var_sam_columns = %w(CHROM POS ID REF ALT TYPE)

    @bcf_query = columns.map { |col| "%#{col}"}.join("\t")

    @col_indices = {
      :svlen => columns.index("SVLEN"),
      :svtype => columns.index("SVTYPE"),
      :type => columns.index("TYPE"),
      :gts => columns.size,
      :variants => variants_columns.map { |col| columns.index(col) }.compact,
      :var_sam => var_sam_columns.map { |col| columns.index(col) }.compact
    }
  end

  def get_index(key)
    @col_indices[key]
  end

  def convert_vcf(files)
    files.each_with_index do |vcf, i|
      count = "#{i+1} / #{ARGV.size} files"
      puts "#{DateTime.now.to_s} START vcftobeacon #{count}"

      @vcf_file = vcf

      # create a .norm file
      normalize_vcf_file

      # create .variants.data, .variants.matching.sample.data, and .sampless.data files
      open_data_files
      generate_variants_data
      generate_samples_data
      close_data_files

      puts "#{DateTime.now.to_s} DONE vcftobeacon #{count}"
      puts
    end
  end

  def normalize_vcf_file
    puts "#{DateTime.now.to_s} Filtering/Coding/Normalizing file #{@vcf_file}"
    if File.exists?(@opts[:chrs_file])
      system("bcftools filter -e 'N_ALT == 0' #{@vcf_file} | bcftools annotate --rename-chrs #{@opts[:chrs_file]} | bcftools norm -m -both | bcftools +fill-tags -o #{@vcf_file}.norm")
    else
      puts "Warning: chromosome rename file #{@opts[:chrs_file]} is not found"
      system("bcftools filter -e 'N_ALT == 0' #{@vcf_file} | bcftools norm -m -both | bcftools +fill-tags -o #{@vcf_file}.norm")
    end
  end

  def open_data_files
    prepare_variants_data
    prepare_variants_sample
    prepare_samples_data
  end
  
  def prepare_variants_data
    # https://github.com/ga4gh-beacon/beacon-elixir/blob/master/deploy/db/db/data/1_chr21_subset.variants.csv
    # datasetId;chromosome;position;variantId;reference;alternate;end;svType;svLength;variantCount;callCount;sampleCount;frequency;sampleMatchingCount
    # 1;21;9411238;rs559462325;G;A;;SNP;;1;5008;2504;0.000199681;1
    # 1;21;9411244;rs181691356;C;A;;SNP;;4;5008;2504;0.000798722;4
    variants_header = %w(datasetId chromosome position variantId reference alternate end svType svLength variantCount callCount sampleCount frequency sampleMatchingCount)
    @variants_file = File.open("#{@vcf_file}.variants.data", "w")
    @variants_file.puts variants_header.join(@opts[:separator]) if @opts[:header]
  end

  def prepare_variants_sample
    # https://github.com/ga4gh-beacon/beacon-elixir/blob/master/deploy/db/db/data/1_chr21_subset.variants.matching.samples.csv
    # datasetId;chromosome;position;variantId;reference;alternate;svType;sampleId
    # 1;21;9411238;rs559462325;G;A;SNP;{HG01029}
    # 1;21;9411244;rs181691356;C;A;SNP;{HG01104,NA11995,NA19922,NA20807}
    var_sam_header = %w(datasetId chromosome position variantId reference alternate svType sampleId)
    @var_sam_file = File.open("#{@vcf_file}.variants.matching.samples.data", "w")
    @var_sam_file.puts var_sam_header.join(@opts[:separator]) if @opts[:header]
  end

  def prepare_samples_data
    # https://github.com/ga4gh-beacon/beacon-elixir/blob/master/deploy/db/db/data/1_chr21_subset.samples.csv
    # sampleId;datasetId
    # HG00096;1
    # HG00097;1
    samples_header = %w(sampleId datasetId)
    @samples_file = File.open("#{@vcf_file}.samples.data", "w")
    @samples_file.puts samples_header.join(@opts[:separator]) if @opts[:header]
  end

  def generate_variants_data
    puts "#{DateTime.now.to_s} Genarating variants data files #{@vcf_file}.variants.data, #{@vcf_file}.variants.matching.sample.data"
    cmd = "bcftools query --allow-undef-tags -f '#{@bcf_query}[\t%SAMPLE=%GT]\n' #{@vcf_file}.norm"
    IO.popen(cmd, "r").each do |line|
      ary = line.strip.split("\t")

      # Assign null value for the SVLEN (sv_length) column if bcftools puts '.' there as the column is typed as integer
      if ary[get_index(:svlen)] == '.'
        ary[get_index(:svlen)] = ''
      end
      # If there is a value in "SVTYPE" other than ".", replace "TYPE" column value with it
      if ary[get_index(:svtype)] != '.'
        ary[get_index(:type)] = ary[get_index(:svtype)]
      end

      gts = ary[(get_index(:gts))..-1]
      hetero = gts.select { |gt| gt[/=(0\|1|1\|0|0\/1|1\/0)/] }
      homalt = gts.select { |gt| gt[/=(1\|1|1\/1)/] }
      samples = (hetero + homalt).map { |gt| gt.sub(/=.*/, '') }
      num_hetero = hetero.size
      num_homalt = homalt.size
      sum_hethom = samples.size
      
      @variants_file.puts [ @opts[:dataset_id], ary.values_at(*get_index(:variants)), sum_hethom ].flatten.join(@opts[:separator])
      @var_sam_file.puts  [ @opts[:dataset_id], ary.values_at(*get_index(:var_sam)), "{#{samples.join(",")}}" ].flatten.join(@opts[:separator])
    end
    puts "#{DateTime.now.to_s} Finished generating variants data files"
  end

  def generate_samples_data
    puts "#{DateTime.now.to_s} Generating file #{@vcf_file}.samples.data"
    cmd = "bcftools query -l #{@vcf_file}"
    IO.popen(cmd, "r").each do |sample_id|
      @samples_file.puts [ sample_id.strip, @opts[:dataset_id] ].join(@opts[:separator])
    end
    puts "#{DateTime.now.to_s} Finished generating samples data file"
  end

  def close_data_files
    @variants_file.close
    @var_sam_file.close
    @samples_file.close
  end
      
end

VcfToBeacon.new(ARGV)



