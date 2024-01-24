# maize_microbiome_transcriptomics

## Description

Pipeline for analyzing maize microbiome metataxonomics and transcriptomics data using Nextflow.

## Researchers

 * Dr. Renato Augusto Corrêa dos Santos
 * Professor Dr. Jason Wallace
 * Professor Dr. Diego M. Riaño-Pachón

## Dependencies

### Metataxonomics pipeline

 * Nextflow version 23.10.0.5889
 * Python 3.8
 * Conda or Miniconda
 * Qiime2 2023.9 ([qiime2-tiny-2023.9-py38-linux-conda.yml](https://data.qiime2.org/distro/tiny/qiime2-tiny-2023.9-py38-linux-conda.yml))
 * 

### Transcriptomics pipeline



## How to run

```bash
nextflow maize_metataxonomics.nf -with-report -with-dag
```

With `-with-report`, Nextflow will generate a HTML report with details about individuals runs, including reasons why some process failed. `-with-dag` enables the generation of a graph with the workflow (great for presentation).

## Institutional Support

 * Computational, Evolutionary, and Systems Biology Laboratory (LabBCES)
 * Nuclear Center of Energy in Agriculture (CENA), University of Sao Paulo, Piracicaba, SP, Brazil
 * University of Georgia, Athens, GA, USA

## Funding

 * São Paulo Research Foundation (FAPESP), São Paulo, SP, Brazil (grant number [2023/11133-3](https://bv.fapesp.br/en/bolsas/212537/integrating-metataxonomics-and-host-transcriptomics-data-in-maize/))

