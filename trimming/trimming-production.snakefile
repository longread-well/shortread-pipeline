include: "common.snakefile"

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

rule make_data_symlinks:
	input:
		read1 = "results/fastuniq/results/mergelanes/{centre}_{platform}_{sample}_{idx}_R1.fq.gz",
		read2 = "results/fastuniq/results/mergelanes/{centre}_{platform}_{sample}_{idx}_R2.fq.gz"
	output:
		read1 = "results/trimming-test/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_untrimmed.fq.gz",
		read2 = "results/trimming-test/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_untrimmed.fq.gz"
	shell: """
		ln -s {input.read1} {output.read1}
		ln -s {input.read2} {output.read2}
	"""

