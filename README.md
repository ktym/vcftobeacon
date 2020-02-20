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

## Options

You can supply options via environmental/shell variables.

```
% BeaconDatasetID=5 BeaconDataSeparator=";" BeaconChrsFile=/path/to/conv_file.txt BeaconDataHeader=true ruby vcftobeacon.rb /path/to/*.vcf > vcftobeacon.log
```

Or you can override default values and environmental/shell variables by the command line options.

```
% ruby vcftobeacon.rb -i 123 -s ";" -c /path/to/conv_file.txt -p /path/to/*.vcf > vcftobeacon.log
```

Each option has a longer form as well.

```
% ruby vcftobeacon.rb --datasetid 123 --separator ";" --chrsfile /path/to/conv_file.txt --header /path/to/*.vcf > vcftobeacon.log
```

Note that the command line options override the environmental/shell variables.

### Dataset ID

To change the dataset ID,

* set the environmental/shell variable `BeaconDatasetID`
* or use the option `--datasetid` / `-i`
* (deafult: 1)

### Column separator

To change the separator (delimiter) of columns, e.g., change to ";" for SQL,

* set the environmental/shell variable `BeaconDataSeparator`
* or use the option `--separator` / `-s`
* (default: "\t")

### Column headers

To print column headers,

* set the environmental/shell variable `BeaconDataHeader` to any value, e.g., "true"
* or use the option `--header` / `-p` with no argument
* (default: false)

### Rename chromosome names

To change the file to rename chromosomes,

* set the environmental/shell variable `BeaconChrsFile`
* or use the option `--chrsfile` / `-c`
* (default: "./chr_name_conv.txt")

## Authors

* Toshiaki Katayama <ktym@dbcls.jp> (Database Center for Life Science, Chiba, Japan)
* Dietmar Fernandez-Orth <dietmar.fernandez@crg.eu> (Centre for Genomic Regulation, Barcelona, Spain)

## License

* MIT

## TODO

Couldn't make the first autoincrement id column to null in data; if it's possible, we could use the following shortcut, hopefully:

```
% cat file | psql -h localhost -p 5432 -d elixir_beacon_dev -U microaccounts_dev -c "copy beacon_data_table from stdin with (format 'text')"
```
