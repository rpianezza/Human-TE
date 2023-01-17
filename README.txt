Folder organization explained:

- Human-TE-dynamics: contains the code used for processing the genomes and all its accessory files


- Documentation (R markdown): contains all my documentations for human-TEs analysis

-- HGDP: contains all the analysis done using the HGDP data
--- 0. flo: old files made by Florian at the beginning of the project
--- 1. mq10-all-files: analysis made with all the HGDP samples using mapping quality filtered (>10)
--- 2. no-PCR: analysis made excluding the samples producted with PCR protocols, which are detected as problematic in our analysis
--- 3. PRC: analysis made ONLY on the samples produced with PCR, to better understand the bias created by PCR

-- SGDP: contains all the analysis done using the SGDP data (only with no-PCR samples and with mapping quality 10)

-- other-documentation: contains other Rmd documents important to understand my work

-- other-files: other useful files, such as Jupyter Notebook Python scripts used in the analysis or metadata files created for the datasets


Note that that every documented analysis is divided into three:

1) the folder containing the figures, necessary for creating the GitHub document
2) the raw Rmd script, from which the GitHub document is created (.Rmd)
3) the GitHub document (.md)
