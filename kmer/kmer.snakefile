subsets = [ '30X', 'full' ]
ks = [ '31', '22' ]
bqs = [ 'any', '20' ]

tools = {
	'Rscript': '/apps/well/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0/bin/Rscript'
}

data = config['data']

def getReadsBySample( data ):
	result = {}
	for key in data.keys():
		elt = data[key]
		sample = elt['sample']
		elt['name'] = key
		if sample in result.keys():
			result[ sample ].append( elt )
		else:
			result[ sample ] = [ elt ]
	return result

samples = getReadsBySample( data )
print( samples )

rule all:
	input:
		histogram = [
			"results/kmer/histograms/{technology}/{centre}/{name}_{subset}.k={k}.bq={bq}.histogram.gz".format(
				technology = data[name]['technology'],
				centre = data[name]['centre'],
				name = name,
				subset = subset,
				k = k,
				bq = bq
			)
			for name in config['data'].keys()
			for subset in subsets
			for k in ks
			for bq in bqs
		],
		jf = [
			"results/kmer/histograms/{technology}/{centre}/{name}_full.k=31.bq=any.histogram.gz".format(
				technology = data[name]['technology'],
				centre = data[name]['centre'],
				name = name
			)
			for name in config['data'].keys()
		],
		full_jf = [ "results/kmer/histograms/bysample/{sample}/{sample}.k={k}.bq={bq}.jf".format( sample = sample, k = k, bq = bq ) for sample in samples.keys() for k in ks for bq in bqs ],
		genomescope = [
			"results/kmer/genomescope/{technology}/{centre}/{name}_{subset}.genomescope/{name}_{subset}.k={k}.bq={bq}_summary.txt".format(
            technology = data[name]['technology'],
            centre = data[name]['centre'],
            name = name,
				subset = subset,
				k = k,
				bq = bq
         )
         for name in config['data'].keys()
			for subset in subsets
			for k in ks
			for bq in bqs
		],
		summary = "results/kmer/histograms/summary.tsv"


def makeMinQualOption( bq ):
	# We ignore bases with lower quality than this:
	# NB. Illumina base qualities are binned, the values are 
	# J   F   A   <   7   2   -   (   #
	# 41  37  32  27  22  17  12  7   2
	minBaseQualityChars = { # See https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
		"40": "I",            # quality=40, or P(error) = 0.0001.  (NB. 41 is the max Illumina outputs)
		"30": "?",            # quality=30, or P(error) = 0.001
		"20": "5",            # quality=20, or P(error) = 0.01
		"10": "+"             # quality=10, or P(error) = 0.1
	}
	if bq == "any":
		return ""
	else:
		return '--min-qual-char "%s"' % minBaseQualityChars[bq]

rule count_kmers:
	input:
		files = [
			"results/analysis-ready/reads/{technology}/{centre}/{name}_R{read}_{subset}.fq.gz".format(
				technology = '{technology}',
				centre = '{centre}',
				name = '{name}',
				read = read,
				subset = '{subset}' 
			)
			for read in [ "1", "2" ]
		]
	output:
		jf = temp( "results/kmer/histograms/{technology}/{centre}/{name}_{subset}.k={k}.bq={bq}.jf" )
	params:
		minQualOption = (lambda w: makeMinQualOption( w.bq )),
		threads = 12
	conda: "envs/jellyfish.yaml"
	shell: """
		zcat {input.files} \
		| jellyfish count \
		--threads {params.threads} \
		-m {wildcards.k} \
		-s 10G \
		{params.minQualOption} \
		-C \
		-o {output.jf}.tmp \
		/dev/fd/0
		mv {output.jf}.tmp {output.jf}
	"""

rule compute_histogram:
	input:
		jf = "results/kmer/histograms/{path}/{what}.jf"
	output:
		histogram = "results/kmer/histograms/{path}/{what}.histogram.gz"
	conda: "envs/jellyfish.yaml"
	shell: """
		jellyfish histo {input.jf} | gzip -c > {output.histogram}.tmp
		mv {output.histogram}.tmp {output.histogram}
	"""

rule run_genomescope:
	input:
		histogram = "results/kmer/histograms/{path}/{prefix}.k={k}.bq={bq}.histogram.gz"
	output:
		summary = "results/kmer/genomescope/{path}/{prefix}.genomescope/{prefix}.k={k}.bq={bq}_summary.txt"
	params:
		outputdir = "results/kmer/genomescope/{path}/{prefix}.genomescope",
		script = srcdir( "scripts/genomescope.R" )
	shell: """
      mkdir -p {params.outputdir}
      {tools[Rscript]} --vanilla \
        {params.script} \
        -i {input.histogram} \
        -o {params.outputdir} \
        -p 2 \
        -k {wildcards.k} \
        -n {wildcards.prefix}.k={wildcards.k}.bq={wildcards.bq}
"""

rule count_full_kmers:
	input:
		files = lambda w: [
			"results/analysis-ready/reads/{technology}/{centre}/{name}_R{read}_full.fq.gz".format(
				technology = readSpec['technology'],
				centre = readSpec['centre'],
				name = readSpec['name'],
				read = read
			)
			for readSpec in samples[w.sample]
			for read in [ "1", "2" ]
		]
	output:
		jf = "results/kmer/histograms/bysample/{sample}/{sample}.k={k}.bq={bq}.jf"
	params:
		minQualOption = (lambda w: makeMinQualOption( w.bq )),
		threads = 12
	conda: "envs/jellyfish.yaml"
	shell: """
		zcat {input.files} \
		| jellyfish count \
		--threads {params.threads} \
		-m {wildcards.k} \
		-s 10G \
		{params.minQualOption} \
		-C \
		-o {output.jf}.tmp \
		/dev/fd/0
		mv {output.jf}.tmp {output.jf}
	"""

rule summarise:
	input:
		histograms = [
			"results/kmer/histograms/{technology}/{centre}/{name}_{subset}.k={k}.bq={bq}.histogram.gz".format(
				name = name,
				technology = data[name]['technology'],
				centre = data[name]['centre'],
				subset = subset,
				k = k,
				bq = bq
			)
			for name in data.keys()
         for subset in subsets
         for k in ks
         for bq in bqs
		]
	output:
		summary = "results/kmer/histograms/summary.tsv"
	params:
		path_template = lambda w: "results/kmer/histograms/{technology}/{centre}",
		identifier_template = lambda w: "{name}_{subset}.k={k}.bq={bq}"
	run:
		outputLines = [[ 'identifier', 'name', 'technology', 'centre', 'subset', 'k', 'bq', 'total', 'error_@5', 'error%_@5', 'error_@3', 'error%_@3' ]]
		print( data.keys() )
		for name in data.keys():
			print( name )
			print( data[name] )
			for subset in subsets:
				for k in ks:
					for bq in bqs:
						identifier = params.identifier_template.format(
							name = name,
							technology = data[name]['technology'],
							centre = data[name]['centre'],
							subset = subset,
							k = k,
							bq = bq
						)
						path = params.path_template.format(
							technology = data[name]['technology'],
							centre = data[name]['centre'],
							subset = subset,
							k = k,
							bq = bq
						)
						print( identifier )
						import gzip
						with gzip.open( path + "/" + identifier + ".histogram.gz", 'r' ) as input:
							lines = input.readlines()
							lines = [ elt.strip().split() for elt in lines ]
							counts = [ float(elt[0]) * float(elt[1]) for elt in lines ]
							total = sum( counts )
							error3 = sum( counts[0:3] )
							error5 = sum( counts[0:5] )
							outputLines.append([
								identifier, name, data[name]['technology'], data[name]['centre'], subset, k, bq,
								str(total), str(error5), '%.5f' % (100 * error5/total), str(error3), '%.5f' % ( 100 * error3/total)
							])
		with open( output.summary, 'w' ) as outputFile:
			outputFile.writelines( '%s\n' % '\t'.join( line ) for line in outputLines )


