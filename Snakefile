configfile : 'config.json' # list of SRA runs by accession number
data = '/data/' # put data here

sets = config['accs']
reference = data + 'GCA_009858895.3.fasta' # place SARS-CoV-2 reference genome here

time = '/usr/bin/time'

#----------------------------------------------------------------------

rule master :
    input :
        # alignment and rates
        data + 'alignment_rates.csv',
        expand(data + '{s}.bwa.fasta', s = sets),

        # fasta versions of fastq files
        expand(data + '{s}.dist.fasta.gz', s = sets)

#----------------------------------------------------------------------

# generate consensus sequence from vcf (against reference)
rule consensus :
    input :
        ref = reference,
        vcf = '{path}.vcf.gz',
        tbi = '{path}.vcf.gz.tbi' # just a check

    output : '{path}.fasta'
    log : '{path}.fasta.log'

    shell : '''

  bcftools consensus -f {input.ref} {input.vcf} \
    -o {output} > {log} 2>&1 '''

# pile the reads upon the reference genome, call variants, normalize
# and index (https://www.biostars.org/p/367960/)
rule call_variants :
    input :
        ref = reference,
        bam = '{path}.bam'

    output :
        vcf = '{path}.vcf.gz',
        tbi = '{path}.vcf.gz.tbi'

    log : '{path}.calls.log'

    shell : '''

  bcftools mpileup -f {input.ref} {input.bam} -Ou 2> {log} \
    | bcftools call -mv -Ou 2>> {log} \
    | bcftools norm -f {input.ref} -Oz -o {output.vcf} >> {log} 2>&1

  touch {output.vcf} && tabix {output.vcf} >> {log} 2>&1
  touch {output.tbi} '''

# align sequences with bwa and store direclty as bam
rule bwa_mem :
    input :
        ref = reference,
        r1 = '{path}_1.fastq.gz',
        r2 = '{path}_2.fastq.gz'

    output : '{path}.bwa.bam'

    log :
        log = '{path}.bwa.log',
        time = '{path}.bwa.time'

    shell : '''

  bwa index {input.ref}
  {time} -vo {log.time} \
    bwa mem {input.ref} {input.r1} {input.r2} 2> {log.log} \
      | samtools sort -o {output} - >> {log.log} 2>&1 '''

# obtain sequences from NCBI (Entrez)
#----------------------------------------------------------------------

# gzip a fastq file
rule gzip :
    input : '{path}.fastq'
    output : '{path}.fastq.gz'
    shell : 'gzip {input} && touch {output}'

# fetch a fastq file
rule fetch :
    params :
        lambda wildcards : wildcards.a.rsplit('/',1)[1]

    output :
        r1 = temp('{a}_1.fastq'),
        r2 = temp('{a}_2.fastq')

    log : '{a}.fetch.log'

    shell : '''

  fasterq-dump {params} --outdir {data} > {log} 2>&1
  touch {output} '''

# convert to fasta (for k-mers computation, etc.)
#----------------------------------------------------------------------

# convert (zipped) fastq to fasta (https://www.biostars.org/p/85929/)
rule convert :
    input :
        r1 = '{path}_1.fastq.gz',
        r2 = '{path}_2.fastq.gz'

    output : '{path}.dist.fasta.gz'
    log : '{path}.dist.log'

    shell : '''

  zcat {input.r1} \
    | awk '{{ if(NR%4 == 1) printf(">%s left\\n", substr($0,2)); \
             else if(NR%4 == 2) print; }}' 2> {log} \
    | gzip > {output} 2>> {log}

  zcat {input.r2} \
    | awk '{{ if(NR%4 == 1) printf(">%s right\\n", substr($0,2)); \
             else if(NR%4 == 2) print; }}' 2>> {log} \
    | gzip >> {output} 2>> {log} '''

# collect some statistics
#----------------------------------------------------------------------

# gather alignment rates
rule alignment_rates :
    input :
        expand(data + '{s}.bwa.bam.arate.csv', s = sets)

    params : data
    output : data + 'alignment_rates.csv'
    log : data + 'alignment_rates.log'
             
    shell : '''

  echo "Sample,AlignmentRate(%)" > {output} 2> {log}
  cat {params}/*.bwa.bam.arate.csv >> {output} 2>> {log}
  echo "---,---" >> {output} 2>> {log}
  awk -F, '{{s += $2; c += 1}} END {{printf "AverageRate,%lf\\n", s/c}}' \
    {output} >> {output} 2>> {log} '''

# alignment rate
rule alignment_rate :
    input : '{path}.bam'
    output : '{path}.bam.arate.csv'
    log : '{path}.bam.arate.log'

    shell : '''

  samtools stats {input} 2> {log} \
    | egrep "^#|^SN" 2>> {log} \
    | python3 scripts/alignment-rate.py > {output} 2>> {log} '''
