subsets = [ '5X', '10X', '20X', '30X', 'full' ]

data = config['data']

rule all:
	input:
		[
			"results/analysis-ready/reads/{technology}/{centre}/{name}_R{read}_{subset}.fq.gz".format(
				technology = data[name]['technology'],
				centre = data[name]['centre'],
				name = name,
				read = read,
				subset = subset
			)
			for name in config['data'].keys()
			for read in [ "1", "2" ]
			for subset in subsets
		],
		[
			"results/analysis-ready/counts/{technology}/{centre}/{name}_R{read}_{subset}.fq.gz.count".format(
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
		read1 = "results/analysis-ready/reads/{technology}/{centre}/{name}_R1_full.fq.gz",
		read2 = "results/analysis-ready/reads/{technology}/{centre}/{name}_R2_full.fq.gz"
	shell: """
		# gzip files can be concatenated
		cat {input.read1} > {output.read1}
		cat {input.read2} > {output.read2}
	"""

#rule truncate_data:
#	input:
#		read = "results/analysis-ready/reads/illumina/whg/whg_illumina_NA12878_R{read}_full.fq.gz"
#	output:
#		read = "results/analysis-ready/reads/illumina/whg/whg_illumina_NA12878_100bp_R{read}_full.fq.gz",
#	conda: "envs/seqtk.yaml"
#	shell: """
#		seqtk trimfq -e 51 {input.read} | gzip -c > {output.read}
#	"""

rule thin_data:
	input:
		read = "results/analysis-ready/reads/{technology}/{centre}/{name}_R{read}_full.fq.gz"
	output:
		read = "results/analysis-ready/reads/{technology}/{centre}/{name}_R{read}_{coverage}X.fq.gz",
	conda: "envs/seqtk.yaml" 
	params:
		# Seed must be identical between the two reads to keep them in sync.
		seed = "2143",
		# Compute number of reads for (roughly) a given coverage.
		# Formula is: # reads = number of bases needed / read length
		# and number of bases needed = (genome length) * coverage
		# The precise genome length here is taken from GIAB b38 genome.
		number_of_reads = lambda w: int(( int( w.coverage ) * 3209286105.0) / int( data[w.name]['read_length'] * 2.0 ))
	shell: """
		seqtk sample -2 -s {params.seed} {input.read} {params.number_of_reads} | pigz -c -p 7 > {output.read}
	"""

rule count_file:
	input:
		read = "results/analysis-ready/reads/{technology}/{centre}/{name}_R{read}_{subset}.fq.gz"	
	output:
		countfile = "results/analysis-ready/counts/{technology}/{centre}/{name}_R{read}_{subset}.fq.gz.count"	
	shell: """
		zcat {input.read} | wc -l > {output.countfile}
	"""
