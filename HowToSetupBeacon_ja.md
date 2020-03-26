# Beacon v1 構築手順

## セットアップ

1. ソースコードをクローンする

    ```bash
    $ mkdir beacon
    $ cd beacon
    $ git clone https://gitlab.com/martaferri/beacon-v1
    $ git clone https://github.com/ktym/vcftobeacon.git
    ```

1. データベースを永続化するため、`beacon-v1/deploy/beacon.yml`を以下のように変更する

    ```diff
    $ diff -u beacon.yml.orig beacon.yml
    --- beacon.yml.orig
    +++ beacon.yml
    @@ -18,6 +18,7 @@
         image: postgres:9.6-alpine
         volumes:
           - ./db:/docker-entrypoint-initdb.d
    +      - ./data/postgres:/beacon-db
         networks:
           - priv
    ```

1. 以下のコマンドを実行して、サービスを起動する

    ```bash
    $ cd beacon-v1/deploy
    $ make all-up
    ```

1. 全てのコンテナが正常に起動していることを確認する

    ```bash
    $ docker-compose ps
       Name                  Command               State           Ports
    -----------------------------------------------------------------------------
    beacon-1      python -m beacon_api             Up      0.0.0.0:5050->5050/tcp
    beacon-db-1   docker-entrypoint.sh postgres    Up      5432/tcp
    beacon-ui-1   aiohttp-wsgi-serve beaconu ...   Up      0.0.0.0:8000->8000/tcp
    ```

1. [http://localhost:8000](http://localhost:8000)にアクセスする

    サンプルクエリ: `Y : 2655179 G > A`


## VCFデータの変換

1. BCFtoolsのインストール

    コンバータの内部で[BCFtools](http://samtools.github.io/bcftools)を使用しているため、インストールされていない場合は下記のコマンドでインストールする。

    #### macOSの場合

    ```bash
    $ brew install bcftools
    ```

    #### ソースコードからコンパイルする場合

    ```bash
    $ git clone git://github.com/samtools/htslib.git
    $ git clone git://github.com/samtools/bcftools.git
    $ cd bcftools
    # The following is optional:
    $ autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
    $ make
    $ sudo make install
    ```

    curlヘッダが見つからないエラーが出る場合

    ```bash
    # RedHat系の場合
    $ sudo yum install curl-devel

    # Debian系の場合
    $ sudo apt-get install libcurl-dev
    ```

1. 下記コマンドでVCFを変換する

    ```bash
    $ ruby vcftobeacon.rb path/to/vcf
    ```

    `sample.vcf`というファイル名のVCFを変換した場合、

    - sample.vcf.samples.data
    - sample.vcf.variants.data
    - sample.vcf.variants.matching.samples.data

    の3つのファイルが生成される。

1. ロード対象のデータを`beacon-v1/deploy/db/mydata`にコピーする

    ```bash
    $ mkdir beacon-v1/deploy/db/mydata
    $ mv sample.vcf.*.data beacon-v1/deploy/db/mydata/
    ```

1. ロードスクリプト`beacon-v1/deploy/db/load-data.sh`を編集する

    ```bash
    $ cd beacon-v1/deploy/db
    $ cp load-data.sh load-data.sh.orig
    ```

    `load-data.sh`を下記内容に置き換え、8行目の変数`VCF`に前の手順でコピーしたファイルのプレフィックス部分を設定する。

    ```bash
    #!/usr/bin/env bash

    echo "Loading initial data into the beacon database"

    # Move to the script directory
    pushd $(dirname ${BASH_SOURCE[0]})

    VCF="mydata/sample.vcf"

    VARIANTS="${VCF}.variants.data"
    SAMPLES="${VCF}.samples.data"
    MATCHING="${VCF}.variants.matching.samples.data"

    for f in "$VARIANTS" "$SAMPLES" "$MATCHING"; do
      if [[ ! -f "$f" ]]; then
        echo "$f not found"
        exit 1
      fi
    done

    # Note: we hide the .sql files inside data/ so that the entrypoint does not pick them up
    # and execute them. We are doing it ourselves.
    # But we do want to run the .sh file

    # Using docker_process_sql() from the /docker-entrypoint-initdb.d
    # This file will be sourced because it is not executable
    # See (for postgres 9.6-alpine):
    # https://github.com/docker-library/postgres/blob/34df4665bfdccf28deac2ed2924127b94489a576/9.6/alpine/docker-entrypoint.sh#L139-L145

    # Schemas and Functions
    docker_process_sql < data/schemas.sql
    docker_process_sql < data/query_data_summary_response.sql
    docker_process_sql < data/query_data_response.sql

    docker_process_sql -c "GRANT EXECUTE ON FUNCTION public.query_data_summary_response(text, integer, integer, integer, integer, integer, integer, character varying, text, text, text, text, text) TO ${POSTGRES_USER};"
    docker_process_sql -c "GRANT EXECUTE ON FUNCTION public.query_data_response(text, integer, integer, integer, integer, integer, integer, character varying, text, text, text, text, text) TO ${POSTGRES_USER};"

    # Datasets
    docker_process_sql < data/init.sql

    # Variants
    docker_process_sql -c "copy beacon_data_table (dataset_id,chromosome,start,variant_id,reference,alternate,\"end\",\"type\",sv_length,variant_cnt,call_cnt,sample_cnt,frequency,matching_sample_cnt) from stdin using delimiters E'\t' csv" < "$VARIANTS"

    # # Change datasetId 1 into datasetId 2 and 3 (datasetId is at position 1)
    # awk -F ';' -v OFS=';' '{if ($1 == "1") $1="2"; print}' data/chr21_subset.variants.csv | \
    # docker_process_sql -c \
    #     "copy beacon_data_table (dataset_id,chromosome,start,variant_id,reference,alternate,\"end\",\"type\",sv_length,variant_cnt,call_cnt,sample_cnt,frequency,matching_sample_cnt) from stdin using delimiters ';' csv header"
    # awk -F ';' -v OFS=';' '{if ($1 == "1") $1="3"; print}' data/chr21_subset.variants.csv | \
    # docker_process_sql -c \
    #     "copy beacon_data_table (dataset_id,chromosome,start,variant_id,reference,alternate,\"end\",\"type\",sv_length,variant_cnt,call_cnt,sample_cnt,frequency,matching_sample_cnt) from stdin using delimiters ';' csv header"

    # docker_process_sql -c \
    # "copy beacon_data_table (dataset_id,chromosome,start,variant_id,reference,alternate,\"end\",\"type\",sv_length,variant_cnt,call_cnt,sample_cnt,frequency,matching_sample_cnt) from stdin using delimiters ';' csv header" < data/chrY_subset.variants.csv

    # Sample list
    docker_process_sql -c "copy tmp_sample_table (sample_stable_id,dataset_id) from stdin using delimiters E'\t' csv" < "$SAMPLES"

    # Change datasetId 1 into datasetId 2 and 3 (datasetId is at position 2)
    # awk -F ';' -v OFS=';' '{if ($2 == "1") $2="2"; print}' data/chr21_subset.samples.csv | \
    # docker_process_sql -c \
    # "copy tmp_sample_table (sample_stable_id,dataset_id) from stdin using delimiters ';' csv header"
    # awk -F ';' -v OFS=';' '{if ($2 == "1") $2="3"; print}' data/chr21_subset.samples.csv | \
    # docker_process_sql -c \
    # "copy tmp_sample_table (sample_stable_id,dataset_id) from stdin using delimiters ';' csv header"

    docker_process_sql -c "INSERT INTO beacon_sample_table (stable_id)
             SELECT DISTINCT t.sample_stable_id
             FROM tmp_sample_table t
             LEFT JOIN beacon_sample_table sam ON sam.stable_id=t.sample_stable_id
             WHERE sam.id IS NULL"

    docker_process_sql -c "INSERT INTO beacon_dataset_sample_table (dataset_id, sample_id)
             SELECT DISTINCT dat.id AS dataset_id, sam.id AS sample_id
             FROM tmp_sample_table t
             INNER JOIN beacon_sample_table sam ON sam.stable_id=t.sample_stable_id
             INNER JOIN beacon_dataset_table dat ON dat.id=t.dataset_id
             LEFT JOIN beacon_dataset_sample_table dat_sam ON dat_sam.dataset_id=dat.id AND dat_sam.sample_id=sam.id
             WHERE dat_sam.id IS NULL"

    # Samples where each variant is found
    docker_process_sql -c "copy tmp_data_sample_table (dataset_id,chromosome,start,variant_id,reference,alternate,\"type\",sample_ids) from stdin using delimiters E'\t' csv" < "$MATCHING"

    # docker_process_sql -c \
    # "copy tmp_data_sample_table (dataset_id,chromosome,start,variant_id,reference,alternate,\"type\",sample_ids) from stdin using delimiters ';' csv header" \
    # < data/chr21_subset.variants.matching.samples.csv

    # Change datasetId 1 into datasetId 2 and 3 (datasetId is at position 1)
    # awk -F ';' -v OFS=';' '{if ($1 == "1") $1="2"; print}' data/chr21_subset.variants.matching.samples.csv | \
    # docker_process_sql -c "copy tmp_data_sample_table (dataset_id,chromosome,start,variant_id,reference,alternate,\"type\",sample_ids) from stdin using delimiters ';' csv header"
    # awk -F ';' -v OFS=';' '{if ($1 == "1") $1="3"; print}' data/chr21_subset.variants.matching.samples.csv | \
    # docker_process_sql -c "copy tmp_data_sample_table (dataset_id,chromosome,start,variant_id,reference,alternate,\"type\",sample_ids) from stdin using delimiters ';' csv header"

    # Finally
    docker_process_sql -c "INSERT INTO beacon_data_sample_table (data_id, sample_id)
             SELECT data_sam_unnested.data_id, s.id AS sample_id
             FROM (
                 SELECT dt.id as data_id, unnest(t.sample_ids) AS sample_stable_id
                 FROM tmp_data_sample_table t
                 INNER JOIN beacon_data_table dt ON dt.dataset_id=t.dataset_id
                                                AND dt.chromosome=t.chromosome
                            AND dt.variant_id=t.variant_id
                            AND dt.reference=t.reference
                            AND dt.alternate=t.alternate
                            AND dt.start=t.start
                            AND dt.type=t.type
             )data_sam_unnested
             INNER JOIN beacon_sample_table s ON s.stable_id=data_sam_unnested.sample_stable_id
             LEFT JOIN beacon_data_sample_table ds ON ds.data_id=data_sam_unnested.data_id AND ds.sample_id=s.id
             WHERE ds.data_id IS NULL"

    docker_process_sql -c "TRUNCATE TABLE tmp_sample_table"
    docker_process_sql -c "TRUNCATE TABLE tmp_data_sample_table"

    # Do some updates now that everything is loaded
    #docker_process_sql < data/updates.sql

    popd
    echo "Initial data loaded"
    ```

1. データをロードする

    ```bash
    # コンテナが起動している場合はdownする
    $ docker-compose down

    # PostgreSQLのデータディレクトリが存在する場合、init-dbが実行されないのでディレクトリを削除する
    $ rm -rf data/postgres

    $ docker-compose up -d

    $ docker-compose logs -f db
    ```

    ログ出力に`CREATE TABLE`や`INSERT`の後にPostgreSQLの再起動したメッセージが表示されればロードが完了している。

    ```
    waiting for server to shut down...LOG:  received fast shutdown request
    LOG:  aborting any active transactions
    LOG:  autovacuum launcher shutting down
    LOG:  shutting down
    ...LOG:  database system is shut down
     done
    server stopped

    PostgreSQL init process complete; ready for start up.

    LOG:  database system was shut down at 2020-03-24 16:33:16 UTC
    LOG:  MultiXact member wraparound protections are now enabled
    LOG:  database system is ready to accept connections
    LOG:  autovacuum launcher started
    ```

1. ロードしたデータを確認する

    ```bash
    $ docker-compose exec db psql -h localhost -p 5432 -d beacon -U beacon -c 'select * from beacon_data_table limit 10;'
     id | dataset_id | chromosome |                 variant_id                  | reference | alternate |  start   |   end    | type | sv_length | variant_cnt | call_cnt | sample_cnt | matching_sample_cnt |  frequency
    ----+------------+------------+---------------------------------------------+-----------+-----------+----------+----------+------+-----------+-------------+----------+------------+---------------------+-------------
      1 |          1 | 22         | rs587697622                                 | A         | G         | 16050075 | 16050075 | SNP  |           |           1 |     5008 |       2504 |                   1 | 0.000199681
      2 |          1 | 22         | rs587755077                                 | G         | A         | 16050115 | 16050115 | SNP  |           |          32 |     5008 |       2504 |                  32 |  0.00638978
      3 |          1 | 22         | rs587654921                                 | C         | T         | 16050213 | 16050213 | SNP  |           |          38 |     5008 |       2504 |                  37 |  0.00758786
      4 |          1 | 22         | rs587712275                                 | C         | T         | 16050319 | 16050319 | SNP  |           |           1 |     5008 |       2504 |                   1 | 0.000199681
      5 |          1 | 22         | rs587769434                                 | C         | A         | 16050527 | 16050527 | SNP  |           |           1 |     5008 |       2504 |                   1 | 0.000199681
      6 |          1 | 22         | rs587638893                                 | C         | A         | 16050568 | 16050568 | SNP  |           |           2 |     5008 |       2504 |                   2 | 0.000399361
      7 |          1 | 22         | rs587720402                                 | G         | A         | 16050607 | 16050607 | SNP  |           |           5 |     5008 |       2504 |                   5 | 0.000998403
      8 |          1 | 22         | rs587593704                                 | G         | T         | 16050627 | 16050627 | SNP  |           |           2 |     5008 |       2504 |                   2 | 0.000399361
      9 |          1 | 22         | rs587670191                                 | G         | T         | 16050646 | 16050646 | SNP  |           |           1 |     5008 |       2504 |                   1 | 0.000199681
     10 |          1 | 22         | esv3647175;esv3647176;esv3647177;esv3647178 | A         | <CN0>     | 16050654 | 16063474 | CNV  |           |           9 |     5008 |       2504 |                   9 |  0.00179712
    (10 rows)
    ```

1. [http://localhost:8000](http://localhost:8000)にアクセスする

    サンプルクエリ: `22 : 16050075 A > G`

    #### beaconが起動していない場合

    ```bash
    $ docker-compose ps
       Name                  Command               State            Ports
    ------------------------------------------------------------------------------
    beacon-1      python -m beacon_api             Exit 1
    beacon-db-1   docker-entrypoint.sh postgres    Up       5432/tcp
    beacon-ui-1   aiohttp-wsgi-serve beaconu ...   Up       0.0.0.0:8000->8000/tcp
    ```

    ```bash
    $ docker-compose restart beacon
    ```

    ```bash
    $ docker-compose ps
       Name                  Command               State           Ports
    -----------------------------------------------------------------------------
    beacon-1      python -m beacon_api             Up      0.0.0.0:5050->5050/tcp
    beacon-db-1   docker-entrypoint.sh postgres    Up      5432/tcp
    beacon-ui-1   aiohttp-wsgi-serve beaconu ...   Up      0.0.0.0:8000->8000/tcp
    ```
