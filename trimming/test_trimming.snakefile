
tools = {
	"fastqc": "~/Projects/Software/3rd_party/FastQC/fastqc",
	"seqtk": "~/Projects/Software/3rd_party/seqtk/seqtk",
	"cutadapt": "cutadapt",
	"multiqc": "multiqc",
	"trimmomatic": "/apps/well/trimmomatic/0.36/trimmomatic-0.36.jar"
}

files = [
	"nygc_illumina_NA12878_100bp_R{read}",
	"nygc_illumina_NA12878_150bp_R{read}",
	"whg_coolmps_HV31_IDX1_R{read}",
	"whg_coolmps_HV31_IDX2_R{read}",
	"whg_coolmps_NA12878_IDX3_R{read}",
	"whg_coolmps_NA12878_IDX4_R{read}",
	"whg_illumina_NA12878_NE725538_R{read}",
	"whg_illumina_NA12878_NE755566_R{read}",
	"whg_illumina_NA12878_NE760511_R{read}",
	"whg_illumina_NA12878_NE796507_R{read}"
]

adapters = {
	"whg_coolmps": "MGI",
	"whg_illumina": "TruSeq_Y",
	"nygc_illumina": "TruSeq_Y"
}

adapter_fastas = {
	"TruSeq_Y": "reference/adapter_sequences/TruSeq_Y.fa",
	"MGI":  "reference/adapter_sequences/MGI.fa"
}

adapter_sequences = {
	"TruSeq_Y": { "fwd": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "rev": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" },
	"Nextera": { "fwd": "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", "rev": "CTGTCTCTTATACACATCTGACGCTGCCGACGA" },
	"10X_linked_reads": { "fwd": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "rev": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATC" },
	"MGI_stLFR": { "fwd": "CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA", "rev": "TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG" },
	"MGI": { "fwd": "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA", "rev": "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG" }
}

def getCentrePlatformIndex( files ):
	return [
		{
			"centre": a[0],
			"platform": a[1],
			"sample": a[2],
			"index": a[3]
		}
		for a in [ elt.split( "_" ) for elt in files ]
	]

index = getCentrePlatformIndex( files )

print( index )

def getAdapters( centre, platform ):
	return adapter_fastas[adapters[ '%s_%s' % ( centre, platform ) ]]

def getAdapterSequences( centre, platform ):
	return adapter_sequences[adapters[ '%s_%s' % ( centre, platform ) ]]

rule all:
	input:
		[
			"results/trimming/{platform}/{centre}/fqc/{centre}_{platform}_{sample}_{idx}_R1_{trimmethod}_fastqc.html".format(
				centre = elt['centre'],
				platform = elt['platform'],
				sample = elt['sample'],
				idx = elt['index'],
				trimmethod = trimmethod
			)
			for elt in index
			for trimmethod in [ 'untrimmed', 'trimmomatic', 'cutadapt' ]
		],
		"results/trimming/multiqc/report.html"

rule make_test_dataset:
	input:
		read1 = "results/fastuniq/results/mergelanes/{centre}_{platform}_{sample}_{idx}_R1.fq.gz",
		read2 = "results/fastuniq/results/mergelanes/{centre}_{platform}_{sample}_{idx}_R2.fq.gz"
	output:
		read1 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_untrimmed.fq.gz",
		read2 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_untrimmed.fq.gz"
	params:
		count = "1000000" # should be multiple of 4
	shell: """
		# the ending '|| true' parts here are to stop the command exiting with 141
		# cf. https://stackoverflow.com/questions/19120263/why-exit-code-141-with-grep-q
		{tools[seqtk]} sample -s 2143 {input.read1} {params.count} | gzip -c > {output.read1}
		{tools[seqtk]} sample -s 2143 {input.read2} {params.count} | gzip -c > {output.read2}
	"""

rule run_trimmomatic:
	input:
		read1 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_untrimmed.fq.gz",
		read2 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_untrimmed.fq.gz"
	output:
		fwd_pair = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_trimmomatic.fq.gz",	
		rev_pair = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_trimmomatic.fq.gz",
		fwd_unpair = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_trimmomatic_unpaired.fq.gz",
		rev_unpair = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_trimmomatic_unpaired.fq.gz",
		trimlog = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_trimmomatic_trimlog.txt.gz"
	params:
		trimmomatic = tools['trimmomatic'],
		adapters = lambda w: getAdapters( w.centre, w.platform ),
		threads = 8,
		trimlog = lambda w, output: output.trimlog.replace( ".gz", "" )
	shell: """
		java -jar {params.trimmomatic} PE \
		-threads {params.threads} \
		-trimlog {params.trimlog} \
		{input.read1} {input.read2} \
		{output.fwd_pair} {output.fwd_unpair} \
		{output.rev_pair} {output.rev_unpair} \
		ILLUMINACLIP:{params.adapters}:2:30:10:2:true MINLEN:50
		gzip {params.trimlog}
	"""

rule run_cutadapt:
	input:
		read1 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_untrimmed.fq.gz",
		read2 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_untrimmed.fq.gz"
	output:
		fwd_pair = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_cutadapt.fq.gz",	
		rev_pair = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_cutadapt.fq.gz",
		log = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_cutadapt.log"
		#info = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_cutadapt_info.txt.gz"
	params:
		adapters_fwd = lambda w: getAdapterSequences( w.centre, w.platform )['fwd'],
		adapters_rev = lambda w: getAdapterSequences( w.centre, w.platform )['rev'],
		threads = 8
		#info = lambda w, output: output.info.replace( ".gz", "" )
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
		fq1 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R1_{trimmethod}.fq.gz",
		fq2 = "results/trimming/{platform}/{centre}/{centre}_{platform}_{sample}_{idx}_R2_{trimmethod}.fq.gz",
		adapters = "reference/adapter_sequences/adapters.tsv"
	output:
		fq1 = "results/trimming/{platform}/{centre}/fqc/{centre}_{platform}_{sample}_{idx}_R1_{trimmethod}_fastqc.html",
		fq2 = "results/trimming/{platform}/{centre}/fqc/{centre}_{platform}_{sample}_{idx}_R2_{trimmethod}_fastqc.html"
	params:
		outputdir = "results/trimming/{platform}/{centre}/fqc"
	shell: """
		{tools[fastqc]} -q -a {input.adapters} -o {params.outputdir} {input.fq1} {input.fq2}
	"""

rule run_multiqc:
	input:
		fastqc = [
			"results/trimming/{platform}/{centre}/fqc/{centre}_{platform}_{sample}_{idx}_R1_{trimmethod}_fastqc.zip".format(
				centre = elt['centre'],
				platform = elt['platform'],
				sample = elt['sample'],
				idx = elt['index'],
				trimmethod = trimmethod
			)
			for elt in index
			for trimmethod in [ 'untrimmed', 'trimmomatic', 'cutadapt' ]
		]
	output:
		report = "results/trimming/multiqc/report.html"
	shell: """
		export LC_ALL=en_GB.utf8
		export LANG=en_GB.utf8
		{tools[multiqc]} -o "results/trimming/multiqc" --filename "report.html" {input.fastqc}
	"""     
