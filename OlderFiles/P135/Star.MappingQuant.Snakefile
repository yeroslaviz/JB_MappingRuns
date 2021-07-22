# the mpping file for either paired-end or songle-end mapping step using the STAR aligner.
############
configfile:"config.yaml"

gz_command="--readFilesCommand zcat" if config["gzipped"] else ""

# tsting for reading in a list pf files from a specific folder:
IDS, = glob_wildcards("rawData/{sample}.conc.R1.fastq.gz")

print(IDS)

rule all:
    input:
#        expand("rawData/fastqc/{sample}_fastqc.zip", sample=IDS ),
#        "samples/fastqc/multiqc.html",
        expand('Analysis/{organism}/star/bamFiles/{sample}.bam', sample = IDS, organism = config['organism']),
        expand('Analysis/{organism}/star/bamFiles/{sample}.bam.bai', sample = IDS, organism = config['organism']),
#        expand('Analysis/star/bamFiles/{sample}.junctions.bed', sample= config["samples"]),
#        'Analysis/diagnostic_plot/diagnostic.pdf',
        expand("Analysis/{organism}/star/bwig/{sample}.bw", sample= IDS, organism = config['organism'])


rule map_star:
    input:
        R1='rawData/{IDS}.conc.R1.fastq.gz',
        R2='rawData/{IDS}.conc.R2.fastq.gz', 
        index=expand("genome/{organisms}/starIndex/", organisms=config['organism'])
    output:
        bam='Analysis/{organism}/star/bamFiles/{IDS}.bam',
        counts="Analysis/{organism}/star/bamFiles/{IDS}.counts.tab",
        SJ = "Analysis/{organism}/star/bamFiles/{IDS}SJ.out.tab",
    benchmark:
        "Analysis/{organism}/benchmarks/{IDS}.run_rRNA_STAR.txt"
    params:
        prefix ="Analysis/{organism}/star/bamFiles/{IDS}", 
        gz_support=gz_command
    threads: 16
    shell:
        r'''
	mkdir -p {params.prefix} 

####### for single-end reads samples data set
#        STAR --runThreadN {threads} --genomeDir {input.index} --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {params.gz_support} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM {config[RAM]} --quantMode GeneCounts --outReadsUnmapped Fastx 
####### for paired-end reads samples data set
        STAR --runThreadN {threads} --genomeDir {input.index} --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {input.R2} {params.gz_support} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM {config[RAM]} --quantMode GeneCounts --outReadsUnmapped Fastx 
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}


        mv {params.prefix}ReadsPerGene.out.tab {output.counts} 
        rmdir {params.prefix}
        '''
#        scripts/STAR_SJtab2JunctionsBed.py -f Analysis/star/bamFiles/{wildcards.sample}SJ.out.tab > \
#        Analysis/star/bamFiles/{wildcards.sample}.junctions.bed # convert the SJ.tab.out  to bed file format

rule index:
    input:
        'Analysis/{organism}/star/bamFiles/{IDS}.bam'
    output:
        'Analysis/{organism}/star/bamFiles/{IDS}.bam.bai'
    shell:
        'samtools index {input}'

rule align_SJ2Bed:
    input:
        "Analysis/star/bamFiles/{sample}SJ.out.tab"
    output:
        "Analysis/star/bamFiles/{sample}.junctions.bed"
    message:
        "Convert STAR's SJ.out.tab to (tophat) junctions.bed BED12 format"
    shell:
        "scripts/STAR_SJtab2JunctionsBed.py -f {input} > {output}"


rule quant:
    input:
        gtf = "genome/genome.gtf",
        bam = expand("Analysis/star/bamFiles/{sample}.bam", sample=IDS)
    output:
        Gene = "Analysis/quantification/GeneCounts.txt",
        mini = "Analysis/quantification/GeneCounts.mini.txt",
    message:
        "running featureCounts"
    threads: config["FeatureCounts"]["threads"]
    log:
        log  = "Analysis/quantification/featureCounts.stat"
    benchmark:
        "Analysis/benchmarks/featureCounts.txt"
    shell:
        """
        featureCounts -T {threads} -a {input.gtf} -t exon -g gene_id -M \
        -o {output.Gene} {input.bam} >>  {log.log} 2>&1
        cut -f 1,7- {output.Gene} | awk 'NR > 1' | awk '{{gsub("Analysis/star/bamFiles/", "", $0); print}}' > {output.mini}
        """

rule diagnostic_plot:
    input:
        input = "Analysis/quantification/GeneCounts.mini.txt"
    output:
        "Analysis/diagnostic_plot/diagnostic.pdf"
    shell:
        "Rscript scripts/diagnositcPlots.R {input.input} {output}"

rule chrom_size:
    input:
        fastA = "genome/{organism}.fa"
    output:
        fai = "genome/{organism}.fa.fai",
        chromSize = "genome/{organism}.chromSize"
    shell:
        """
        samtools faidx {input.fastA} &&  \
        cut -f 1,2  {output.fai} > {output.chromSize}
        """

rule create_bigwig:
    input:
        chromSize = "genome/{organism}.chromSize",
        bam = "Analysis/{organism}/star/bamFiles/{sample}.bam",
        bai = "Analysis/{organism}/star/bamFiles/{sample}.bam.bai"
    output:
        bw = "Analysis/{organism}/star/bwig/{sample}.bw"
    params:
        dir ="Analysis/{organism}/star/bwig", 
        prefix ="Analysis/{organism}/star/bwig/{sample}"
    shell:
        """
        mkdir -p {params.dir} && \
        bam2wig.py  -i {input.bam} -s {input.chromSize} -o {params.prefix} &> {params.prefix}.log
        rm -f {wildcards.sample}.wig
        """
