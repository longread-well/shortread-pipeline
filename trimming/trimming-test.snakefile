include: "trimming-common.snakefile"

rule all:
	input:
		[
			"results/trimming-test/{platform}/{centre}/fqc/{centre}_{platform}_{sample}_{idx}_R1_{trimmethod}_fastqc.html".format(
				centre = elt['centre'],
				platform = elt['platform'],
				sample = elt['sample'],
				idx = elt['index'],
				trimmethod = trimmethod
			)
			for elt in index
			for trimmethod in [ 'untrimmed', 'trimmomatic', 'cutadapt' ]
		],
		"results/trimming-test/multiqc/report.html"

rule make_test_dataset:
	input:
		read1 = "results/fastuniq/results/mergelanes/{centre}_{platform}_{sample}_{idx}_R1.fq.gz",
		read2 = "results/fastuniq/results/mergelanes/{centre}_{platform}_{sample}_{idx}_R2.fq.gz"
	output:
		read1 = "results/trimming-test/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_untrimmed.fq.gz",
		read2 = "results/trimming-test/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_untrimmed.fq.gz"
	params:
		number_of_reads = "1000000" # should be multiple of 4
	shell: """
		# the ending '|| true' parts here are to stop the command exiting with 141
		# cf. https://stackoverflow.com/questions/19120263/why-exit-code-141-with-grep-q
		{tools[seqtk]} sample -s 2143 {input.read1} {params.number_of_reads} | gzip -c > {output.read1}
		{tools[seqtk]} sample -s 2143 {input.read2} {params.number_of_reads} | gzip -c > {output.read2}
	"""

