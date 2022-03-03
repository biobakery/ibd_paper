Siyuan Ma, Dmitry Shungin, Himel Mallick, Melanie Schirmer, Long Nguyen, Raivo Kolde, Eric Franzosa, Hera Vlamakis, Ramnik Xavier, Curtis Huttenhower *Population Structure Discovery in Meta-Analyzed Microbial Communities and Inflammatory Bowel Disease.*

* This is the analysis repository for the IBD meta-analysis project. Analysis scripts are ordered, by the sequence of tasks performed, in `rmds/`. `functions/` contains helper functions as well as scripts to generate display items of the paper. `assets/` has additional constants, among which Table 1, i.e., additional metadata of the studies, is included.

* The `rmds/` directory contains analysis scripts including combining otu tables and metadata files from multiple meta-analysis datasets. For the ease of maintenance and user access, a cleaned and pre-filtered phyloseq object at the genus level has been saved at `data/physeq/gennera.RData`. If downloading from github,
the user should start with this file and `rm/ds/3-adjust_batch.Rmd`.

* Many of the analyses in this repository use `MMUPHin` for data batch correction and meta-analysis. MMUPHin has had some interface changes since it was used for this paper. If trying to reproduce the analyses here, the user should download and install the correct old version, at https://github.com/biobakery/MMUPHin/releases/tag/v0.9