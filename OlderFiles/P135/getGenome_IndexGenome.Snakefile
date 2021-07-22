
configfile: "config.yaml"

rule all:
    input:
        expand("genome/{organisms}/starIndex/", organisms=config['organism']),
        expand("genome/{organisms}/bwaIndex/", organisms=config['organism']),
        expand("genome/{organisms}/bowtie2Index/", organisms=config['organism'])

#### get reference genomic data for the mapping (fasta and gtf files). Links are added in the config file
rule get_genome:
   output:
         fastA="genome/{organism}.fa",
         gtf="genome/{organism}.gtf"
   message:
          "getting {config[organism]} data - fastA and gtf files from Ensembl"
   shell:
         """
         wget -nc -O - {config[fastA]} | gunzip -c - > {output.fastA}
         wget -nc -O - {config[gtf]}   | gunzip -c - > {output.gtf}
         """

### Indexing the reference genome

rule star_index:
   input:
          fasta="genome/{organism}.fa",
          gtf="genome/{organism}.gtf"
   output:
          directory("genome/{organism}/starIndex/")
   threads: 16
   params:
         prefix = {config["organism"]}
   shell:
          "mkdir -p {output} && "
          "STAR --runThreadN {threads} "
          "--outFileNamePrefix {output}{params.prefix} "
          "--runMode genomeGenerate "
          "--genomeDir {output} "
          "--limitGenomeGenerateRAM {config[RAM]} "
          "--genomeSAindexNbases {config[SAindex]} "
          "--genomeFastaFiles {input.fasta} "
          "--sjdbGTFfile {input.gtf} "
          "--sjdbOverhang 100"

rule bwa_index:
    input:
          fasta = "genome/{organism}.fa"
    output:
          directory("genome/{organism}/bwaIndex/")
    params:
        prefix = {config["organism"]}

    message:
        "Indexing {params.prefix} genome for bwa Mapping"
#    log:
#        "log/bwa_index.log"
    shell:
        "bwa index -b {config[blockSize]} -p {output}{params.prefix} {input.fasta} "


# rule bowtie2_index
rule bowtie2_index:
    input:
          fasta = "genome/{organism}.fa"
    output:
          directory("genome/{organism}/bowtie2Index/")
    params:
        prefix = {config["organism"]}
    message:
        "Indexing {params.prefix} genome for bowtie2 Mapping"
#    log:
#        "log/bwa_index.log"
    shell:
        "bowtie2-build {input.fasta} {output}{params.prefix}"


	

