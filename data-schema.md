# Data schema

## Nanopore reads

Input data to the Zika pipeline arrives in the `data/` directory. This should be mounted to `data/` in the Docker container.

  - `data`
    - `usvi-library1-2016-12-10` - library
      - `raw_reads` - squiggle graphs in fast5 format
        - `pass` - contains `.fast5` files
      - `basecalled_reads` - basecalled with Metrichor
        - `pass_demultiplex` - demultiplexed basecalled reads
          - `NB01` - contains `.fast5` files for NB01 barcode
          - `NB02` - contains `.fast5` files for NB02 barcode
          - etc...
        - `nonNB_demultiplexed` - demultiplexed basecalled reads
          - `BC01` - contains `.fast5` files for BC01 barcode
          - `BC02` - contains `.fast5` files for BC02 barcode
          - etc...        
        - `fail` - contains `.fast5` files that weren't demultiplexed

## Sample metadata

Sample metadata for the Zika pipeline arrives in the `samples/` directory. This should be mounted to `samples/` in the Docker container.

  - `samples/`
    - `samples.tsv` - line list of sample metadata
    - `runs.tsv` - line list of run metadata

### `samples.tsv`

Must be `tsv` formatted. Keyed off of column headers rather than column order.

sample_id | strain  | collection_date | country | division     | location
--------- | ------- | --------------- | ------- | ---------- | ---------
ZBRD116   | ZBRD116 | 2015-08-28      | brazil  | alagoas    | arapiraca
ZBRC301   | ZBRC301 | 2015-05-13      | brazil  | pernambuco | paulista

### `runs.tsv`

Must be `tsv` formatted. Keyed off of column headers rather than column order.

run_name | barcode_id | sample_id | primer_scheme
-------- | ---------- | --------- | -------------
library1 | NB01       | ZBRD116   | v2_500.amplicons.ver2
