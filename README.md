# Short-read platform comparison pipelines
This folder contains scripts used to assess short-read sequencing performance.

Warning: files in this folder will be publicly shared!
Do not put any sensitive or environment-specific details here.

##Contributing

We are using the following structure:

    parent folder/
	    pipelines/  # this folder
           step1/
                Snakefile
			  step2/
					Snakefile
			  etc.
       results/    # output files
           step1

The intention is to run all pipelines from the parent folder, and that all paths will be relative to the parent folder.
Where possible, conda environments can be used wthin snakemake to include needed software.

Ultimately we will aim to run the pipeline through a top-level configuration file which lists all datasets to be processed.

## Contributors ##

- Gavin Band (WHG)
- David Flores (WHG)

