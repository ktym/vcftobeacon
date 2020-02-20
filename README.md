# vcftobeacon

A script to convert VCF files into GA4GH Beacon database files for PostgreSQL with the help of bcftools.

## Prerequisites

* [bcftools](http://samtools.github.io/bcftools/bcftools.html) with a +fill-tags plugin
  * try `bcftools +fill-tags -vv` to see if you have the plugin
* [Ruby](https://www.ruby-lang.org/) 1.8 or later

## Usage

```
% ruby vcftobeacon.rb /path/to/*.vcf > vcftobeacon.log
```

You can also supply options via shell variables.


```
% BeaconDatasetID=5 BeaconDataSeparator=";" BeaconDataHeader=true ruby vcftobeacon.rb /path/to/*.vcf > vcftobeacon.log
```

## Authors

* Toshiaki Katayama <ktym@dbcls.jp> (Database Center for Life Science, Chiba, Japan)
* Dietmar Fernandez-Orth <dietmar.fernandez@crg.eu> (Centre for Genomic Regulation, Barcelona, Spain)

