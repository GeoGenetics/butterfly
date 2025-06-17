# Butterfly
This project automates the process of downloading, processing, and visualizing data 

## Limitations

* SMDB entries must have "country_ocean", "field_sample_parent_id" and either "Master Depth (cm)" or "Median Master Age" specified.
* ~~SMDB headers have changed, code must be updated~~
* ~~MetaDMG needs to be updated to handle a combination of ss and ds libraries in one file~~
* FilterBAM needs to be reimplemented to run on large files 
* Controls are not analysed, due to upstream issues

### Minor issues

* Some hardcoded variables in metaDMG should be moved to config
* Change from Selenium DB grabber to internal solution
* Must check what values "age", can take e.g. is it CE, BCE, BP? 
* Need to setup the snakemake to run as a cron, taking input from the latest files
* Make the minimum number of reads in the R part an argument 
* ~~Make sure the IDs in the R part refer to library_id or archive_id depending on context ~~
* Does good or bad fit filtering reflect reality?

## Overview
The pipeline consists of three main steps:

1. Metadata gathering 
2. Post-mapping processing
3. Filtering and visualisation


## Components

### 1. Metadata gathering

#### Downloading from SMDB

Every day, a cronjob on dandy03fl runs `scripts/download_db.sh`, which runs `scripts/download.py` which uses Selenium to access [SMDB](http://dandyweb01fl.unicph.domain:5100/search), and download the megatable ("merged table") with all library data. These are stored in `db/` with the current date.

#### Generating lists of finished libraries

The script `scripts/get_results.sh` also runs daily on dandy03fl. It finds libraries for which the aggregated metaDMG results are more than 30 minutes old. These are libraries which have completed the entire mapping pipeline. The basename and dirnames for these files are stored in the `files/` directory.

### 2. Post-mapping processing

The snakemake workflow `scripts/pmp.smk` runs the post-mapping pipeline.

It is run like this.
``` 
snakemake -s scripts/pmp.smk -j10 -pn --configfile config.yaml --config testrun=True smdb=db/2024-11-05.query_result_full.csv finished_libs=files/testfiles.csv

----- configfile -----

nodes:			nodes taxonomy file 
names: 			names taxonomy file 
acc2tax: 		mapping of accession to taxa file 
base: 			base out directory 


----- config arguments -----

testrun:		should we make fake data for testing purposes [False]
smdb: 			smdb file 
finished_libs:			finished lib file 
```

First, it reads SMDB and the finished libraries list and merges these. It finds the path of the final bam file output of the pipeline, and removes rows that do not have the needed metadata (see Limitations). The base directory for output is specified by `base`, and then organised by country `country_ocean` and master core ID `field_sample_parent_id`, e.g. `/projects/caeg/data/pmp/Cambodia/CAM2313M/`.

```
rule save_data: 		Outputs csv files with merged information for each core 
rule reassign: 			Runs filterBAM reassign on each library bam file 
rule filterbam: 		Runs filterBAM filter on the reassigned bams
rule merge_bams: 		Merges each filtered library bam file to sample level (archive_id)
rule sort_bam: 			Sort merged files by QNAME
rule metadmg_lca: 		Run metaDMG-lca on sorted merged files 
rule metadmg_getdmg: 		Run metaDMG-getdmg on sorted merged files 
rule metadmg_dfit: 		Run metaDMG-dfit on sorted merged files 
rule metadmg_aggregate: 	Aggregate all metaDMG results 
rule filterbam_getstats: 	Run filterBAM filter (with no filtering specified) on merged bam files 
rule generate_report: 		Runs the last part of the pipeline, the filtering and visualisation R script 
```

### 3. Filtering and visualisation

This is the final step of the pipeline that does all the filtering and visualisation `final_pipeline.R`. It reads in metadata, metaDMG data, and filterBAM data, creates several plots, constructs a damage model over time, and outputs Excel files of genus identified at each level in the core. It classifies all taxa as either a good or bad fit to the theoretical damage model (smiley plots).

It works on reads assigned to a species level, but aggregates these to the genus level.


#### Damage model

The damage model is constructed by inspecting d-fit (`A_b`) values over time. Currently, it gets the d-fit of all plant taxa with at least 500 reads, and plots this through time (age or depth). For the top 10 and bottom 10 layers, all taxa are used, and for the middle layers, just those with good fit are used. 

We plot this as a line of damage through the core, and smooth with a loess fit. We then get the 5th percentile of this, and use it as an absolute lower point for damage. Taxa with damage below this point (at each level of the core) are removed. This assigns `status = (pass, fail)`.


#### Filtering

All plots refer to Eukaryotes classified to the species level, with ANI of over 95, an expected breadth ratio of over 0.8, and either gini under 0.6 or entropy above 0.7. Strat plots and the final Excel files use just taxa with `status == "pass"` and a minimum of 2 occurences (in different samples) of at least 100 reads.


#### Output

##### Plots
* Plots of number of reads sequenced, collapsed and assigned
* Proportion of reads coming from each superkingdom over time 
* Proprtion of reads classified to different ranks 
* Damage in all eukaryotes (taxa with at least 100 reads)
* Damage in plants with at least 500 reads 
* The damage model and the data used to construct it 
* Examples of "good" or "bad" taxa according to damage fit 
* Plants and animals damage after filtering with the damage model
* Damage, expected breadth, ANI and gini of filtered animals and plants (only top 20 genera)
* Strat plots of top animal and plant genera 

##### Excel
Excel file with four tabs, the columns are genus names, and the rows are dates 
* Number of reads from animal taxa
* Number of reads from plant taxa 
* Proportion of reads from animal taxa 
* Proportion of reads from plant taxa 

##### Files
* `agg.tsv`: aggregated merged statistics from filterBAM and metaDMG
* `saved_taxids.tsv`: mapping of accessions to taxids for this dataset 
* `qdata.tsv`: damage model info

## Flowchart 

![butterfly](https://github.com/user-attachments/assets/06337188-2d75-4bdc-b15c-6865c0eabd68)

