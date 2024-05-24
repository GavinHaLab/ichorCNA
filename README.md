# *ichorCNA*
ichorCNA is a tool for estimating the fraction of tumor in cell-free DNA from ultra-low-pass whole genome sequencing (ULP-WGS, 0.1x coverage). This is a version maintained by the laboratory of Gavin Ha.

## How to run ichorCNA

1. WDL pipleine
   - The [WDL pipeline](https://github.com/GavinHaLab/ichorCNA_WDL/tree/main/WDL) uses release [ichorCNA:v0.5.0](https://github.com/GavinHaLab/ichorCNA/releases/tag/v0.5.0)
   - The docker image `gavinhalab/ichorcna:1.0.0` can be found at [DockerHub](https://hub.docker.com/repository/docker/gavinhalab/ichorcna/general)
   - Clone [WDL pipleine](https://github.com/GavinHaLab/ichorCNA_WDL/tree/main/WDL)
   - To run the WDL pipeline, please follow instructions [here](https://github.com/GavinHaLab/ichorCNA_WDL/tree/main/WDL#readme)
2. Snakemake pipeline
   - The [Snakemake pipeline](https://github.com/GavinHaLab/ichorCNA/tree/v0.4.0/scripts/snakemake) is part of release [ichorCNA:v0.4.0](https://github.com/GavinHaLab/ichorCNA/releases/tag/v0.4.0)
   - Clone [ichorCNA:v0.4.0](https://github.com/GavinHaLab/ichorCNA/tree/v0.4.0) or download [v0.4.0.tar.gz](https://github.com/GavinHaLab/ichorCNA/releases/tag/v0.4.0)
   - To run the Snakemake pipeline, please follow instructions [here](https://github.com/broadinstitute/ichorCNA/wiki/SnakeMake-pipeline-for-ichorCNA)

## ichorCNA Wiki Page
**For more details on usage/pipelines, outputs, and FAQs, please visit the [GitHub Wiki page for ichorCNA](https://github.com/broadinstitute/ichorCNA/wiki)**

## Description
ichorCNA uses a probabilistic model, implemented as a hidden Markov model (HMM), to simultaneously segment the genome, predict large-scale copy number alterations, and estimate the tumor fraction of a ultra-low-pass whole genome sequencing sample (ULP-WGS). 

The methodology and probabilistic model are described in:  
Adalsteinsson, Ha, Freeman, et al. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. (2017) Nature Communications Nov 6;8(1):1324. [doi: 10.1038/s41467-017-00965-y](https://doi.org/10.1038/s41467-017-00965-y)

The analysis workflow consists of 2 tasks:  
1. GC-content bias correction (using HMMcopy)  
  a. Computing read coverage from ULP-WGS  
  b. Data correction and normalization  
2. CNA prediction and estimation of tumor fraction of cfDNA

## Contacts
If you have any questions or feedback, please contact us at:  
**Email:** <gha@fredhutch.org> or <pchandra@fredhutch.org>

## Acknowledgements
ichorCNA is maintained by Gavin Ha, Pooja Chandra, and Michael Yang.  

ichorCNA was originally developed in collaboration with  
- **Blood Biopsy Group**, Group Leader **Viktor Adalsteinsson**, Broad Institute of MIT and Harvard
- Laboratory of **Matthew Meyerson**, Medical Oncology, Dana-Farber Cancer Institute
- Laboratory of **J. Christopher Love**, Koch Institute for integrative cancer research at MIT
- Laboratory of **Gad Getz**, Cancer Program, Broad Institute

## Software License
ichorCNA
Copyright (C) 2024 Gavin Ha Lab

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
