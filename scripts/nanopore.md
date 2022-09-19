## this is the script summary for the nanopore data processing

### basecalling

### filt the low quality base

#### nanoplot

    NanoPlot --fastq combined.fastq -t 4 --plots hex,dot,kde -o nanoplot # --maxlength 40000
    NanoPlot --summary sequencing_summary.txt --loglength -o summary
    NanoPlot --fastq clean.nanofilt.fastq -t 4 --plots {dot,kde} -o nanoplot
    # conda install -c bioconda nanoplot -y

#### nanofilt

    NanoFilt -q 7 -l 300 --logfile nanofilt.log combined.fastq > nanofilt.fastq

### convert from fastq to fasta
    
    seqtk seq -a nanofilt.fastq > nanofilt.fasta

### rename nanofilt fasta to oganized id fasta
    
    python rename.py nanofilt.fasta nanofilt.rename.fasta

## blast to identify ARGs in sequences

### make a blast database from your sample
    
    mkdir blast
    makeblastdb -in nanofilt.rename.fasta -dbtype 'nucl' -parse_seqids -out blast/TEM -logfile blast/makeblastdb.txt 

### blast your forward and reverse primer to your database (your sample)

    blastn -db blast/TEM -query primer.fasta -task blastn-short  -outfmt "6 std qlen slen qcovs" -max_target_seqs 10000000 -perc_identity 80 -num_threads 10 -out blast/primer_vs_txt 

### get the fasta containing all the ARGs
#### option1
    
    seqtk subseq input.fasta id.txt > extracted.id.fasta

#### option2
    
    filterbyname.sh in=input.fa out=output.fa names=ID.txt include=t fixjunk overwrite=t

### pipe into emu environment

    conda activate meta
    emu abundance emu_results/r1tem.f3tem.1000.fasta --output-dir emu_results/ --threads 40

