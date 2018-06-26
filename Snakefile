configfile:  "config.yaml"

targets = []
for sample in config['references']['samples']:
    targets.extend(
        expand(
            'data/counts/{sample}.htseq.counts.txt',
            sample=sample)
        )

rule all:
    input:
        targets

rule trimgalore:
    input:
        r1 = 'data/raw_reads/{sample}_R1.fastq.gz',
        r2 = 'data/raw_reads/{sample}_R2.fastq.gz'
    output:
        r1 = 'data/trimmed_reads/{sample}_R1_trimmed.fq',
        r2 = 'data/trimmed_reads/{sample}_R2_trimmed.fq'
    params:
        output_dir = 'data/trimmed_reads/'
    shell:
        'trim_galore --paired {input.r1} {input.r2} '
        '--output_dir {output_dir}'

rule star_index:
    input:
        ref_folder = 'data/reference/hg38',
        ref_fasta = 'data/reference/hg38/hg38.fa',
        ref_gtf = 'data/reference/hg38/hg38.gtf'
    output:
        'data/reference/hg38/SAindex'
    threads: 8
    shell:
        'STAR --runThreadN {threads} --runMode genomeGenerate '
        '--genomeDir {input.ref_folder} --genomeFastaFiles {input.ref_fasta} '
        '--sjdbOverhang 100 --sjdbGTFfile {input.ref_gtf}'

rule starAlign_first:
    input:
        genomeDir = 'data/reference/hg38',
        SA_index = 'data/reference/hg38/SAindex',
        r1 = 'data/trimmed_reads/{sample}_R1_trimmed.fq',
        r2 = 'data/trimmed_reads/{sample}_R2_trimmed.fq'
    output:
        'data/sam_files/hg38/{sample
