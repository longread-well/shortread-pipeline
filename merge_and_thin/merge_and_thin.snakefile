subsets = [ '10X', '20X', '30X', '40X', 'full' ]

data = config['data']

rule all:
	input:
		[
			"results/analysis-ready/{technology}/{centre}/{name}_R{read}_{subset}.fq.gz".format(
				technology = data[name]['technology'],
				centre = data[name]['centre'],
				name = name,
				read = read,
				subset = subset
			)
			for name in config['data'].keys()
			for read in [ "1", "2" ]
			for subset in subsets
		]

rule merge_data:
	input:
		read1 = lambda w: [
			"results/trimming/{technology}/{centre}/{name}{index}_R1_trimmomatic.fq.gz".format(
				technology = '{technology}',
				centre = '{centre}',
				name = '{name}',
				index = index
			)
			for index in data[w.name]['indices']
		],
		read2 = lambda w: [
			"results/trimming/{technology}/{centre}/{name}{index}_R2_trimmomatic.fq.gz".format(
				technology = '{technology}',
				centre = '{centre}',
				name = '{name}',
				index = index
			)
			for index in data[w.name]['indices']
		]
	output:
		read1 = "results/analysis-ready/{technology}/{centre}/{name}_R1_full.fq.gz",
		read2 = "results/analysis-ready/{technology}/{centre}/{name}_R2_full.fq.gz"
	shell: """
		# gzip files can be concatenated
		cat {input.read1} > {output.read1}
		cat {input.read2} > {output.read2}
	"""

rule thin_data:
	input:
		read1 = "results/analysis-ready/{technology}/{centre}/{name}_R1_full.fq.gz",
		read2 = "results/analysis-ready/{technology}/{centre}/{name}_R2_full.fq.gz"
	output:
		read1 = "results/analysis-ready/{technology}/{centre}/{name}_R1_{coverage}X.fq.gz",
		read2 = "results/analysis-ready/{technology}/{centre}/{name}_R2_{coverage}X.fq.gz"
	conda: "envs/seqtk.yaml" 
	params:
		# Seed must be identical between the two reads to keep them in sync.
		seed = "2143",
		# Compute number of reads for (roughly) a given coverage.
		# Formula is: # reads = number of bases needed / read length
		# and number of bases needed = (genome length) * coverage
		# The precise genome length here is taken from GIAB b38 genome.
		number_of_reads = lambda w: int(( int( w.coverage ) * 3209286105.0) / int( data[w.name]['read_length'] ))
	shell: """
		seqtk sample -s {params.seed} {input.read1} {params.number_of_reads} | gzip -c > {output.read1}
		seqtk sample -s {params.seed} {input.read2} {params.number_of_reads} | gzip -c > {output.read1}
	"""

