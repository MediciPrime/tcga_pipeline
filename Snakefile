configfile:  "config.yaml"

targets = []
for reference in config['references']:
    targets.extend(
        expand(
            'data/sam_files/{reference}/Aligned.out.sam',
            reference=reference)
    )

rule all:
    input:
        targets

def get_raw(wc):
    raw_reads = []
    for v in config['references'][wc.reference]:
        raw_reads.extend(v)
    return expand('data/raw_reads/{sample}/{sample}_R{num}.fastq',
                  sample=samples, num=[1,2])

def get_trimmed(wc):
    trimmed_reads = []
    for v in config['references'][wc.reference]:
        trimmed_reads.extend(v)
        return expand('data/trimmed_reads/{sample}/{sample}_R{num}_val_{num}.fastq',
                      sample=samples, num=[1,2])

        
rule trimgalore:
    input:
        raw_reads = get_raw
    output:
        trimmed_reads = get_trimmed
    params:
        output_dir = 'data/trimmed_reads/{sample}/'
    shell:
        'trim_galore --paired {input.raw_reads} '
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

def get_samples(wc):
    samples = []
    for v in config['references'][wc.reference]:
        samples.extend(v)
    return expand('data/trimmed_reads/{sample}/{sample}_R{num}_val_{num}.fq',
                  sample=samples, num=[1,2])
        
rule starAlign_first:
    input:
        genomeDir = 'data/references/{reference}',
        SA_index = 'data/references/{reference}/SAindex',
        trimmed_reads = get_samples
    output:
        'data/sam_files/{reference}/SJ.out.tab'
    params:
        sam_prefix = 'data/sam_files/{reference}/'
    threads: 4
    shell:
        'STAR --runThreadN {threads} --genomeDir {input.genomeDir} '
        '--readFilesIn {input.r1} {input.r2} --alignIntronMax 500000 '
        '--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 '
        '--outFilterMismatchNmax 10 --alignMatesGapMax 1000000 '
        '--sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory '
        '--sjdbOverhang 100 --outFileNamePrefix {params.sam_prefix} '
        '--outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 '
        '--outSAMstrandField intronMotif --outSAMtype None --outSAMmode None'

rule star_index_2:
    input:
        ref_fasta = 'data/references/{reference}/{reference}.fa',
        sj_outtab = 'data/sam_files/{reference}/SJ.out.tab'
    output:
        genomeOut = 'data/sam_files/{reference}/SAindex'
    params:
        genomeDirOut = 'data/sam_files/{reference}/'
    shell:
        'STAR --runThreadN {threads} --runMode genomeGenerate '
        '--genomeDir {params.genomeDirOut} --genomeFastaFiles {input.ref_fasta} '
        '--sjdbOverhang 100 --sjdbFileChrStartEnd {input.sj_outtab}'
        
rule starAlign_second:
    input:
        trimmed_reads = get_samples,
        SA_index = 'data/sam_files/{reference}/SAindex'
    output:
        'data/sam_files/{reference}/Aligned.out.sam'
    params:
        genomeDir = 'data/sam_files/{reference}/'
    threads: 4
    shell:
        'STAR --runThreadN {threads} --genomeDir {params.genomeDir} '
        '--readFilesIn {input.trimmed_reads} --alignIntronMax 500000 '
        '--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 '
        '--outFilterMismatchNmax 10 --alignMatesGapMax 1000000 --sjdbScore 2 '
        '--alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory '
        '--limitBAMsortRAM 0 --outFilterMatchNminOverLread 0.33 '
        '--outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 '
        '--outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS '
        '--outSAMunmapped Within --outSAMtype BAM SortedByCoordinate '
        '--outSAMheaderHD @HD VN:1.4'
