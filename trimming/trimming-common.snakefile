data = config['data']
adapter_sequences = config['adapter_sequences']

def getAdapterSequences( name ):
	return adapter_sequences[data[name]['adapters']]

def reverseComplement( sequence ):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	return ''.join( [ complement[base] for base in sequence ][::-1] )

wildcard_constraints:
	idx = "(_[0-9A-Za-z]+)!"

rule createAdapterFastas:
	output:
		fasta = "results/{outputdir}/adapters/{name}.fa"
	run:
		name = wildcards.name
		sequences = adapter_sequences[name]
		lines = {}
		# trimmomatic paired-end primers, note naming is important!
		lines[ 'Prefix%s_adapter/1' % name ] = sequences['fwd']
		lines[ 'Prefix%s_adapter/2' % name ] = sequences['rev']
		# trimmomatic single-end primers
		lines[ '%s_adapter1' % name ] = sequences['fwd']
		lines[ '%s_adapter1_rc' % name ] = reverseComplement( sequences['fwd'] )
		lines[ '%s_adapter2' % name ] = sequences['rev']
		lines[ '%s_adapter2_rc' % name ] = reverseComplement( sequences['rev'] )
		o = open( output.fasta, 'w' )
		for key in lines.keys():
			o.writelines([
				'>%s\n' % key,
				"%s\n" % lines[key]
			])
		o.close()

rule run_trimmomatic:
	input:
		read1 = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R1_untrimmed.fq.gz",
		read2 = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R2_untrimmed.fq.gz",
		adapters = lambda w: "results/{outputdir}/%s.fa" % adapter_sequences[data[w.name]['adapters']],
	output:
		fwd_pair = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R1_trimmomatic.fq.gz",	
		rev_pair = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R2_trimmomatic.fq.gz",
		fwd_unpair = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R1_trimmomatic_unpaired.fq.gz",
		rev_unpair = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R2_trimmomatic_unpaired.fq.gz",
		trimlog = "results/{outputdir}/{technology}/{centre}/{name}{idx}_trimmomatic_trimlog.txt.gz"
	params:
		threads = 8,
		trimlog = lambda w, output: output.trimlog.replace( ".gz", "" )
	conda: "envs/trimmomatic.yaml"
	shell: """
		echo running trimmomatic
		trimmomatic PE \
		-threads {params.threads} \
		-trimlog {params.trimlog} \
		{input.read1} {input.read2} \
		{output.fwd_pair} {output.fwd_unpair} \
		{output.rev_pair} {output.rev_unpair} \
		ILLUMINACLIP:{input.adapters}:2:30:10:2:true MINLEN:50
		gzip {params.trimlog}
	"""

rule run_cutadapt:
	input:
		read1 = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R1_untrimmed.fq.gz",
		read2 = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R2_untrimmed.fq.gz"
	output:
		fwd_pair = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R1_cutadapt.fq.gz",	
		rev_pair = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R2_cutadapt.fq.gz",
		log = "results/{outputdir}/{technology}/{centre}/{name}{idx}_cutadapt.log"
		#info = "results/{outputdir}/{technology}/{centre}/{name}{idx}_cutadapt_info.txt.gz"
	params:
		adapters_fwd = lambda w: getAdapterSequences( w.name )['fwd'],
		adapters_rev = lambda w: getAdapterSequences( w.name )['rev'],
		threads = 8
		#info = lambda w, output: output.info.replace( ".gz", "" )
	conda: "envs/cutadapt.yaml"
	shell: """
		cutadapt \
		-j {params.threads} \
		-a {params.adapters_fwd} \
		-A {params.adapters_rev} \
		--nextseq-trim 5 \
		--minimum-length 50 \
		-o {output.fwd_pair} \
		-p {output.rev_pair} \
		{input.read1} \
		{input.read2} 2>&1 > {output.log}
	"""

rule run_fastqc:
	input:
		fq1 = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R1_{trimmethod}.fq.gz",
		fq2 = "results/{outputdir}/{technology}/{centre}/{name}{idx}_R2_{trimmethod}.fq.gz",
		adapters = "reference/adapter_sequences/adapters.tsv"
	output:
		fq1 = "results/{outputdir}/{technology}/{centre}/fqc/{name}{idx}_R1_{trimmethod}_fastqc.html",
		fq2 = "results/{outputdir}/{technology}/{centre}/fqc/{name}{idx}_R2_{trimmethod}_fastqc.html",
		zip1 = "results/{outputdir}/{technology}/{centre}/fqc/{name}{idx}_R1_{trimmethod}_fastqc.zip",
		zip2 = "results/{outputdir}/{technology}/{centre}/fqc/{name}{idx}_R2_{trimmethod}_fastqc.zip"
	params:
		outputdir = "results/{outputdir}/{technology}/{centre}/fqc"
	conda: "envs/fastqc.yaml"
	shell: """
		fastqc -q -a {input.adapters} -o {params.outputdir} {input.fq1} {input.fq2}
	"""

rule run_multiqc:
	input:
		fastqc = [
			"results/{outputdir}/{technology}/{centre}/fqc/{name}{idx}_R{read}_{trimmethod}_fastqc.zip".format(
				name = key,
				centre = data[key]['centre'],
				technology = data[key]['technology'],
				idx = idx,
				trimmethod = '{trimmethod}',
				outputdir = '{outputdir}',
				read = read
			)
			for key in data.keys()
			for idx in data[key]['indices']
			for read in [ "1", "2" ]
		]
	output:
		report = "results/{outputdir}/multiqc/report-{trimmethod}.html"
	params:
		outputdir = "results/{outputdir}/multiqc",
		filename = "report-{trimmethod}.html"
	conda: "envs/multiqc.yaml"
	shell: """
		export LC_ALL=en_GB.utf8
		export LANG=en_GB.utf8
		multiqc -o {params.outputdir} --filename {params.filename} {input.fastqc}
	"""     

