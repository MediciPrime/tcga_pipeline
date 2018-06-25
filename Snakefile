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
    
