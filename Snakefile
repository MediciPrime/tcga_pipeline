configfile:  "config.yaml"

targets = []
for reference in config['references']:
    for sample in config['references'][reference]:
        targets.extend(
            expand(
                'data/sam_files/{reference}/{sample}/Aligned.out.sam',
                sample=sample,
                reference=reference)
        )

rule all:
    input:
        targets

rule trimgalore:
    input:
        r1 = 'data/raw_reads/{sample}/{sample}_R1.fastq',
        r2 = 'data/raw_reads/{sample}/{sample}_R2.fastq'
    output:
        r1 = 'data/trimmed_reads/{sample}/{sample}_R1_val_1.fq',
        r2 = 'data/trimmed_reads/{sample}/{sample}_R2_val_2.fq'
    params:
        output_dir = 'data/trimmed_reads/{sample}/'
    shell:
        'trim_galore --paired {input.r1} {input.r2} '
        '--output_dir {params.output_dir}'

rule star_index:
    input:
        ref_folder = 'data/references/{reference}',
        ref_fasta = 'data/references/{reference}/{reference}.fa',
        ref_gtf = 'data/references/{reference}/{reference}.gtf'
    output:
        'data/references/{reference}/SAindex'
    threads: 8
    shell:
        'STAR --runThreadN {threads} --runMode genomeGenerate '
        '--genomeDir {input.ref_folder} --genomeFastaFiles {input.ref_fasta} '
        '--sjdbOverhang 100 --sjdbGTFfile {input.ref_gtf}'

rule starAlign_first:
    input:
        genomeDir = 'data/references/{reference}',
        SA_index = 'data/references/{reference}/SAindex',
        r1 = 'data/trimmed_reads/{sample}/{sample}_R1_val_1.fq',
        r2 = 'data/trimmed_reads/{sample}/{sample}_R2_val_2.fq'
    output:
        'data/sam_files/{reference}/{sample}/Aligned.out.sam'
    params:
        sam_prefix = 'data/sam_files/{reference}/{sample}/'
    threads: 4
    shell:
        'STAR --runThreadN {threads} --genomeDir {input.genomeDir} '
        '--readFilesIn {input.r1} {input.r2} --alignIntronMax 500000 '
        '--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 '
        '--outFilterMismatchNmax 10 --alignMatesGapMax 1000000 '
        '--sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory '
        '--sjdbOverhang 100 '
        '--outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 '
        '--outSAMstrandField intronMotif --outSAMtype None --outSAMmode None'
