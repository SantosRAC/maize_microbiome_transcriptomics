# maize_microbiome_transcriptomics

## Description

Repository storing analyses during [Dr. Santos](https://bv.fapesp.br/pt/pesquisador/164909/renato-augusto-correa-dos-santos/) internship at UGA - Wallace Lab.

Results from analysis of RNAseq, 16S and correlation between both will be included here.

Eventually, it will include a pipeline for analyzing maize microbiome metataxonomics and transcriptomics data using Nextflow.

## Researchers

 * Dr. Renato Augusto Corrêa dos Santos
 * Professor Dr. Jason Wallace
 * Professor Dr. Diego M. Riaño-Pachón

## Dependencies

### Virtual env for all analyses (except Nextflow)

```bash
virtualenv --python=python3.10 .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Institutional Support

 * Computational, Evolutionary, and Systems Biology Laboratory (LabBCES)
 * Nuclear Center of Energy in Agriculture (CENA), University of Sao Paulo, Piracicaba, SP, Brazil
 * Wallace Lab
 * University of Georgia (UGA), Athens, GA, USA

## Funding

 * São Paulo Research Foundation (FAPESP), São Paulo, SP, Brazil (grant number [2023/11133-3](https://bv.fapesp.br/en/bolsas/212537/integrating-metataxonomics-and-host-transcriptomics-data-in-maize/))

