# Data schema

## Nanopore reads

Input data to the Zika pipeline arrives in the `data/` directory. This should be [mounted to `data/` in the Docker container](https://github.com/zibraproject/zika-pipeline/blob/master/docker-notes.md#data-volume).

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
