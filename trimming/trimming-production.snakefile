include: "trimming-common.snakefile"

data = config['data']

localrules: make_data_symlinks

# GB: 15/11/2020:
# Please note this pipeline has been moved over to new configuration file syntax.
# It will likely need some debugging

rule all:
	input:
		[
			"results/trimming/{technology}/{centre}/fqc/{name}{idx}_R1_{trimmethod}_fastqc.html".format(
				name = key,
				centre = data[key]['centre'],
				technology = data[key]['technology'],
				idx = idx,
				trimmethod = trimmethod
			)
			for key in data.keys()
			for idx in data[key]['indices']
			for trimmethod in [ 'trimmomatic' ]
		],
		[ "results/trimming/multiqc/report-{trimmethod}.html".format( trimmethod = t ) for t in [ 'trimmomatic' ] ]

rule make_data_symlinks:
	input:
		read1 = "results/fastuniq/mergelanes/{centre}_{platform}_{sample}_{idx}_R1.fq.gz",
		read2 = "results/fastuniq/mergelanes/{centre}_{platform}_{sample}_{idx}_R2.fq.gz"
	output:
		read1 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_untrimmed.fq.gz",
		read2 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_untrimmed.fq.gz"
	shell: """
		ln -s ../../../../{input.read1} {output.read1}
		ln -s ../../../../{input.read2} {output.read2}
	"""

