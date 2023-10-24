# Metagenomics pipeline & useful scripts


## Contents
- [sources](#sources)
- [Basic awk & sed](#basic-awk--sed)


## Basic awk & sed
[[back to top](#contents)]
*an awk program works on a line by line basis and for each line of a file it attempts the folloing:*
	
	aw 'CONDITION {ACTIONS}'


Extract columns 2, 4, and 5 from file.txt:

	awk '{print $2, $4, $5}' input.txt

Print each line where the 5th column is equal to 'abc123':

	awk '$5 == "abc123"' file.txt

Print each line where the 5th column is *not* equal to 'abc123':

	awk '$5 != "abc123"' file.txt


Print each line where the 5th column is less than 60:

	awk '$5 < 60' file.txt

sum column 1 of file.txt:

	awk '{sum+=$1}END{print sum}' file.txt

compute the mean of column 2:

	awk '{x+=$2}END{print x/NR}' file.txt

special awk variables:

	$0: original line
	NF: number of fields in the current line (number of columns that awk recognizes)
	NR: number of records, number of lines processed (line number)
	OFS: output fields separator, the character placed between item when printed
	Field separator: -F '\t'


## bbduk

trim adapter and quality control

	bbduk.sh -Xmx40g in1=R1.fastq.gz in2=R2.fsatq.gz ziplevel=9 out1=trimmed/R1.fastq.gz out2=trimmed/R2.fastq.gz ref=database minlen=50 qtrim=rl trimq=20 ktrim=r k=23 mink=11  hdist=1 tpe tbo trimpolya=10 stats=x.log.txt

trim phix reads

	bbduk.sh -Xmx20g in1=R1.fastq.gz in2=R2.fsatq.gz ziplevel=9 out1=NoPhi/R1.fastq.gz out2=Nophi/R2.fastq.gz minlen=50 k=31 ref=phix174_ill.ref.fa.gz hdist=1 stats=log.txt

trim human contamination

	

#### fastqc to check quality
	mkdir fastqc
	fastqc *.fasta.gz -o fastqc -t 40
	multiqc .

## fastax-toolkit


## [seqtk](https://github.com/lh3/seqtk)
*Seqtk is a fast and lighweight tool for processing sequences in the FASTA or FASTQ format. It seamlessly parses both FASTA and FASTQ files which can also be optionally compressed by gzip (most cases lazy to do)*

fastq to fasta convertion -> apply to zipped file

	seqtk seq -a input.fastq.gz > output.fasta

alternatively, using fastax-toolkit

	fastq_to_fasta -i input.fastq -o output.fasta

alternatively, using sed

	sed -n '1~4s/^@/>/p;2~4p' file.fq > file.fa

Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing):

1. seqtk randomly subsample paired fastq or fasta.

		seqtk sample -s 123 read1.fq 100000 > sub_read1.fq
		seqtk sample -s 123 read2.fq 100000 > sub_read2.fq

2. seqtk randomly subsample fastq or fasta
	
		seqtk sample sample.fa 100000 > sub_sample.fa


Interleave PE reads

	seqtk mergepe R1.fasta.gz R2.fasta.gz > interleaved.fasta

Untangle an interleaved PE FASTQ file. 

	seqtk seq -l0 -1 interleaved.fq > deinterleaved_1.fq
	seqtk seq -l0 -2 interleaved.fq > deinterleaved_2.fq



### [bioawk](https://bioinformaticsworkbook.org/Appendix/Unix/bioawk-basics.html#gsc.tab=0)
extract sequences based on headerID  
		
	seqtk subseq input.fasta id.txt > output.fasta

	filterbyname.sh in=input.fa out=output.fa names=ID.txt include=t fixjunk overwrite=t
	> include: t (include name) f (exclude name)

	bioawk -c fastx 'BEGIN{while((getline <"ID.txt")>0)i[k]=1}{if(i[$name])print ">"$name"\n"$seq}' input.fasta > ID_sequence.fasta

extract sequences from virsorter ids

	cat virsorterid.txt | xargs -n 1 samtools faidx all.faa
get sequences length
		
	bioawk -c fastx '{print $name,length($seq)}' input.fa


linearize the fasta file, [benchmark](https://www.biostars.org/p/363676/)

	seqkit seq -w 0 input.fa 
	seqtk seq input.fa
	bioawk -c fastx '{print ">"$name; print $seq}'   

	
### [remove too short reads](http://www.metagenomics.wiki/tools/short-read/filter-remove-too-short-reads-length)

filter sequences length based on certain cut-off
		
	bioawk -c fastx '{ if(length($seq) > 500) { print ">"$name; print $seq }}' input.fasta > outputøø.500bp.fasta 


# Assembly
### Spade  
isolate scripts

	spades.py  --careful -1 R1.fasta.gz -2 R2.fasta.gz -o metaspade_assembly/

metagenomics scripts 

	spades.py  --meta -1 R1.fasta.gz -2 R2.fasta.gz -o metaspade_assembly/

cross-assembly, generate yaml file format

	python metaspades.py --dataset your_yaml_file.yaml -t threads -m MAX_MEMORY_USAGE -k 21,33,55,77,99,127 -o assembled_viral_contigs



filter spades folder - fitler minimum contig length

	for x in *.fasta/;
	do bioawk -c fastx '{ if(length($seq) > 500) { print ">"$name; print $seq }}' ${x}contigs.fasta > ${x}contigs.500bp.fasta
	done	


### rnaspades

### MEGAHIT

filter minimum contig length < 1000bp

	megahit -1 R1.fasta.gz -2 R2.fasta.gz -m 0.75 --presets meta-sensitive -t 50 -o megahit_assembly/ --min-contig-len 1000
	bioawk -c fastx '{ if(length($seq) > 500) { print ">"$name; print $seq }}' contigs.fasta > contigs.500bp.fasta  
	bioawk -c fastx '{ if(length($seq) > 200) { print ">"$name; print $seq }}' masterscaffolds.rnavirome.c95.fasta > masterscaffolds.rnavirome.c95.200bp.fasta


visualize assemblies results

	megahit_toolkit contig2fastg 99 contigs.fa > contigs.fastg
	it can support gz/bz2 extension file*    

	WGS could be 500bp  
	-t/--num-cpu-threads   default auto detect to use up all free CPUs  
	-m/--memory max memory in byte to be used in SdBG construction; default 0.9(if set between 0-1, fraction of the machine's total memory)  
	Presets parameters:<https://wiki.gacrc.uga.edu/wiki/Megahit>

cross assembly

	1. can cat *.R1.fasta.gz > combined.R1.fasta.gz, cat *.R2.fasta.gz > combined.R2.fasta.gz  
	2.`megahit -t 25 -m 0.8  --min-contig-len 1000 -1 pe1_R1.fq,pe2_R1.fq,pe3_R1.fq -2 pe1_R2.fq,pe2_R2.fq,pe3_R2.fq -o outdir` [comma-separated list, option -> input -> output directory, no need add quate]



### Velvet

shuffle Illumina paired reads  

	shuffleSequences_fasta.pl R1.fq R2.fq shuffled.fa  
Find the best assembly for Illumina paired end reads just for k=31

	VelvetOptimiser.pl -s 31 -e 31 -f '-shortPaired -fasta interleaved.fasta' -t 4 --optFuncKmer 'n50'`
	[Velvet assembly](https://github.com/tseemann/VelvetOptimiser/wiki/VelvetOptimiser-Manual)

build the hash index

	velveth assembly 35 -shortPaired -fastq sample.fastq  (under assembly folder contains 3 files, Log, Roadmaps and Sequences).  
	velvetg assembly -cov_cutoff auto -exp_cov auto

### evaluate the assembly quality

quast to check assembly quality

	quast *.fasta -o contig_report/ -m  --min-contig

check contig length

	seqlen.py *.fasta out.txt
	bioawk -c fastx '{print $name,length($seq)}' input.fa

add prefix in sequence header name

	rename.sh in=x.fasta out=y.fasta prefix=XYZ


### CD-HIT-EST

	cd-hit-est -i input.fasta -o output.fasta -c 0.95 -n 9 -T 0[use all threads] -M 58000[RAM<64GB]
	cd-hit-est -i input.fasta -o output.fasta -c 0.95 -T 0

# Taxonomy Annotation
## BLAST

BLAST database: Sequence databases for use with the stand-alone BLAST programs. The files in this directory are pre-formatted databases that are ready to use with BLAST  
 
FASTA BLAST database: Sequence databases in FASTA format for use with the stand-alone BLAST programs. These databases must be formatted using formatdb before they can be used with BLAST.  

Download resources
1. [BLAST website](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) 
2. [FASTA website](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA)

### Download database
NCBInt blast database

	wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz'
	ls nt.*.gz |xargs -n1 tar -xzvf
	rm nt.*.gz

NCBInr blast database

	wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz'
	ls nr.*.gz |xargs -n1 tar -xzvf
	rm nr.*.gz

	for i in nt.*.tar.gz; do parallel -j 75 tar -xzvf ${i}; done

	parallel -j 同时运行数量 tar -xzvf 需要解压文件的路径/{} ::: `ls 需要解压文件的路径`

Diamond fasta database

	wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/fasta/nt.gz'
	gunzip nt.gz
	gzip -d nt.gz (either is ok)


make blast database

	makeblastdb -in OTUs.fasta -dbtype 'nucl' -parse_seqids -out OTUs -logfile logfile.txt

blastn 

	blastn -db BlastformattedDB -query file.fasta  -out file_VS_database.txt -outfmt "6 std slen qlen qcovs staxids" -evalue 1e-5 -perc_identity 50 -num_threads 30 -max_target_seqs N

	blastn -db ucleotide_fasta_protein_homolog_model -query combined_entero.fasta -out entero_vs_card.txt -outfmt "6 std scovs" -evalue 1e-5 -perc_identity 90 -num_threads 40 -max_target_seqs 1

	office workstation: /data1/database/NCBInt/nt 
megablast

	megablast -d \~/tools/NCBIdatabase/blastnt/nt -i scaffolds.500bp.fasta -e 1e-5 -m 8 -a 60 -o scaffolds.500bp_vs_NCBInt.megahit.m8 

blastn very short primer

	makeblastdb -in combined_entero.fasta  -dbtype 'nucl' -parse_seqids -out combined_entero -logfile blast.log.txt ## -task blastn-short

Wordsize:

1. Megablastn 28
2. Blastn 11
3. Blastn-short 7

https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a

blastn header -m8 file

	blastn.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

blast output

	-outfmt <String>
   	alignment view options:
     0 = Pairwise,
     1 = Query-anchored showing identities,
     2 = Query-anchored no identities,
     3 = Flat query-anchored showing identities,
     4 = Flat query-anchored no identities,
     5 = BLAST XML,
     6 = Tabular,
     7 = Tabular with comment lines,
     8 = Seqalign (Text ASN.1),
     9 = Seqalign (Binary ASN.1),
    10 = Comma-separated values,
    11 = BLAST archive (ASN.1),
    12 = Seqalign (JSON),
    13 = Multiple-file BLAST JSON,
    14 = Multiple-file BLAST XML2,
    15 = Single-file BLAST JSON,
    16 = Single-file BLAST XML2,
    17 = Sequence Alignment/Map (SAM),
    18 = Organism Report

entrez direct

	esearch -db assembly -query "GCA_000005845" | elink -target taxonomy | efetch -format native -mode xml | grep ScientificName | awk -F ">|<" 'BEGIN{ORS=", ";}{print $3;}'

	esearch -db nuccore -query "SAMN02441064" | elink -target taxonomy | esummary | xtract -pattern DocumentSummary -element TaxId

	esummary -db nuccore -id NM_002826 | xtract -pattern DocumentSummary -element Caption,Title, TaxId

	cat <filename> | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,Title,TaxId

filter out 

format protein output

	less protein_vs_test.v2.txt|sort -k9 -n|awk '($14>=70)'|column -t


Extract open reading fram Prodigal

	prodigal -p meta -i file -d file-prodigal.fna -a file-prodigal.faa -o file-prodigal.gbk

### Diamond

	diamond makedb --in nr.faa -d nr
	diamond blastp -d \~/tools/NCBIdatabase/Diamondnr/nr -q input.prodigal.faa -o output.prodigal.m8 -e 1e-5
	# cut -f1 input.prodigal.m8|rev|cut -f2- -d"\_"|rev|paste - input.prodigal.m8|cut -f1,3- > input.prodigal.contigs.m8
	cut -f1 input.prodigal.m8|cut -f1-2 -d"_"|paste - input.prodigal.m8|cut -f1,3- > input.prodigal.contigs.m8
	diamond blastp -d /scratch/users/xiaoqiong/database/nr.dmnd -q input.prodigal.faa -o output.prodigal.m8 -e 1e-5

	
	cut -f1 input.prodigal.m8|rev|cut -f2- -d"\_"|rev|paste - input.prodigal.m8|cut -f1,3- > input.prodigal.contigs.m8	
	diamond blastp -d /scratch/database/nr -q meta.scaffolds.faa -o meta.scaffolds.prodigal.m8 -e 1e-5 -p 30
	cut -f1 meta.scaffolds.prodigal.m8|rev|cut -f2- -d"\_"|rev|paste - meta.scaffolds.prodigal.m8|cut -f1,3- > meta.scaffolds.prodigal.contigs.m8

Megan

	import from blast, load accession mapping file, specify mode, e-value
	control+A, Tree-uncollapse all, export to csv, readName_to_taxonPath, assigned, tab

### kraken2

	kraken2 --db ${kraken_db} ${infasta} --output output.kraken --report output.kreport --use-names --threads 40

### [centrifuge](https://gensoft.pasteur.fr/docs/centrifuge/1.0.4-beta/MANUAL)

	centrifuge -x ${centrifuge_db} -q/f ${infastq/a} -S read_classifications.tsv --report-file read_report.summary.tsv -p/--threads 40


# Binning
### metabat2
	metabat2 -i assembly.fa.gz -a depth.txt -o bins_dir/bin -v  
_remember to run bowtie2 to get depth.txt_  this bowtie2 script is for binning quantification
	bowtie2-build -f database.fasta database -p 50
	bowtie2 -x database -1 R1.fastq -2 R2.fastq -S x.sam -p 50 --very-sensitive
	samtools view -bS spades.sam -@ 40 > spades.bam
	samtools sort spades.bam -o spades.sorted.bam -@ 40
	jgi_summarize_bam_contig_depths --outputDepth depth.txt *.sorted.bam
	
### samtools
convert a SAM file to a BAM file

	samtools view -bS SAMPLE.sam > SAMPLE.bam

convert a BAM file to a SAM file

	samtools view -h SAMPLE.bam > SAMPLE.sam #### -h means header

sort a BAM file

	samtools sort SAMPLE.bam -o SAMPLE_sorted.bam

index a sorted BAM file

	samtools index -b SAMPLE_sorted.bam -@ 40 SAMPLE_sorted.bam.bai

bamtobed a sorted BAM file

	bedtools bamtobed -i *.bam > *.bed
	bedtools bamtobed -i *.bam -cigar > *.bed 

sort by readName

	samtools sort -n SAMPLE.bam -o SAMPLE_sorted.bam
	samtools view -f 4 file.bam > unmapped.sam
	samtools view -b -f 4 file.bam > unmapped.bam -@ 40 # get the unmapped reads using parameter f, to get the output in bam using parameter -b
	samtools view -b -F 4 file.bam > mapped.bam -@ 40 # get the mapped reads using parameter F, which works like -v of grep, to get the output in bam using parameter -b

### [virsorter2](https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-bwm5pc86)

Run VirSorter2

	virsorter run --keep-original-seq -i 5seq.fa -w vs2-pass1 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-length 5000 --min-score 0.5 -j 40 all

Run checkV

	checkv end_to_end vs2-pass1/final-viral-combined.fa checkv -t 40 -d /data1/database/checkv-db-v1.0/

Run VirSorter2 again

	cat checkv/proviruses.fna checkv/viruses.fna > checkv/combined.fna 
	virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i checkv/combined.fna -w vs2-pass2 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-length 5000 --min-score 0.5 -j 40 all

Run DRAMv

	1. step 1 annotate
	
	DRAM-v.py annotate -i vs2-pass2/for-dramv/final-viral-combined-for-dramv.fa -v vs2-pass2/for-dramv/viral-affi-contigs-for-dramv.tab -o dramv-annotate --skip_trnascan --threads 28 --min_contig_size 1000
	
	2. step 2 summarize anntotations

	DRAM-v.py distill -i dramv-annotate/annotations.tsv -o dramv-distill

### VirFinder
install

	install.packages("glmnet", dependencies=TRUE)
	install.packages("Rcpp", dependencies=TRUE)
	source("https://bioconductor.org/biocLite.R")
	biocLite("qvalue")
	install.packages("~/tools/VirFinder/linux/VirFinder_1.1.tar.gz", repos = NULL, type="source")  
	library(VirFinder)

To quick start, one can predict the viral contigs using the command

	predResult <- VF.pred("contigs.fa")

sort sequences by p-value in ascending order

	predResult[order(predResult$pvalue),]

estimate q-values (false discovery rates) based on p-values

	predResult$qvalue <- VF.qvalue(predResult$pvalue)

sort sequences by q-value in ascending order

	predResult[order(predResult$qvalue),]


### concoct
	velveth velveth_k71 71 -fasta -shortPaired -separate All_R1.fa All_R2.fa  
	velvetg velveth_k71 -ins_length 400 -exp_cov auto -cov_cutoff auto  

### picard
	
	to run jar file
	java -jar MarkDuplicates.jar

### checkm
	checkm lineage_wf -f CheckM.txt -t 48 -x fa bins_dir/ bins/CheckM

## detect plasmid from WGS

PlasmidSeeker  

	git clone https://github.com/bioinfo-ut/PlasmidSeeker/`
	cd PlasmidSeeker/
	bash plasmidseeker_ecoli_test.sh
	#help download db_w20 database
	cat EC_1_results.txt

PlasmidSpades

	plasmidspades.py -1 R1.fastq.gz -2 R2.fastq.gz -o assembly/ -t 50

Plasflow

	source activate plasflow
	filter_sequences_by_length.pl -input input_dataset.fasta -output filtered_output.fasta -thresh sequence_length_threshold
	PlasFlow.py --input Citrobacter_freundii_strain_CAV1321_scaffolds.fasta --output test.plasflow_predictions.tsv --threshold 0.7

# Quantifiction
## Bowtie2 - short read alignment and samtools, v1.11
	bowtie2-build -f database.fa database -p 40
	bowtie2 -x database -1 R1.fasta.gz -2 R2.fasta.gz -S x.sam -p 40 --very-sensitive-local
	samtools view -bS input.sam > output.bam -@ 40
	samtools sort input.bam -o output.sorted.bam -@ 40
	samtools index -b input_sorted.bam -@ 40 input_sorted.bam.bai
	samtools idxstats x.input_sorted.bam > x.idxstats.txt
	python get_count_table.py *.idxstats.txt > count.txt # sum up all the output.idxstats.txt to the summary table
	python absolutereads2RPKM.py count.txt rpkm.txt # pipeline folder
	
	samtools depth input_sorted.bam > output.depth.txt  # to see the SNP coverage evolution etc. optional




# Ecology analysis
### Primer 7 or R
	PCO, PERMANOVA, Anoism, cluster, heatmap (Clusvis - online tool)
# How to write scripts
	#!/bin/bash  
	for x in *_R1_001.fastq.gz;do echo Repair.sh -Xmx20g in1=$x in2=${x%_R1_001.*}_R2_001.fastq.gz ziplevel=9 out1=repaired\/$x out2=repaired\/${x%_R1_001.*}_R2_001.fastq.gz >> Repair.sh;done;sh Repair.sh
# unassembled viral contigs
Viral contigs with unassembled overlaps or from the same scaffold were merged using the SeqMan program implemented in the Lasergene software package v7.1 (DNAstar)  
Integrated genomics viewer


# getting the mapping coverage
## [bedtools](https://bedtools.readthedocs.io/en/stable/content/tools/genomecov.html)
	bedtools genomecov -ibam input.bam -g output.genome
	samtools depth	-d 0 aln.sorted.bam # set the depth to the maximum, otherwise maximum coverage depth is 8000
	bedtools bamtobed -i *.bam > *.bed # work for caobo's datasets
	samtools bedcov aln.sorted.bam #did not work

## bam2readcount

1. install conda => https://docs.conda.io/en/latest/
2. set up an environment and install bam-readcount using conda in this environment => https://anaconda.org/bioconda/bam-readcount
3. download brc-parser.py, git clone https://github.com/sridhar0605/brc-parser => https://github.com/sridhar0605/brc-parser

example

	bam-readcount -f ref.fa some.bam chrID > some.txt
	-d  --max-count
	-w 1 maximum number of warnings of each type to emit. -1 gives an unlimited number.

add chrID

	RlmH4_FKDL220004689-1a_1_trimmed_cutadapt_bowtie2_unique_header_filter.bam Escherichia_coli_BW25113_K-12_substr_BW25113_tRNA-Leu-TAG-1-1  > RlmH4.txt

	python ../tools/brc-parser/brc-parser.py RlmH4.txt 

to make a for loop in the current folder which contains the bam file:

	for x in *.bam
	do 
		bam-readcount -w 1 -d 10000000 -f E.coli_BW25113_tRNA_reference.renamed.fasta ${x} Escherichia_coli_BW25113_K-12_substr_BW25113_tRNA-Leu-TAG-1-1 
		python ../tools/brc-parser/brc-parser.py ${x%.bam}.txt
	done


for x in *.bam
	do 
		bam-readcount -w 1 -d 10000000 -f /data2/xiaoqiong/augmentin_meta/meta_augmentin_analysis/ecoli_snp/genome_assemblies_genome_fasta/GCF_000008865.2_ASM886v2_genomic.fna ${x} 
		NC_002695.2

		python ../tools/brc-parser/brc-parser.py ${x%.bam}.txt
	done

>NC_002695.2 Escherichia coli O157:H7 str. Sakai DNA, complete genome
>NC_002127.1 Escherichia coli O157:H7 str. Sakai plasmid pOSAK1, complete sequence
>NC_002128.1 Escherichia coli O157:H7 str. Sakai plasmid pO157, complete sequence


## vCONTACT2
### [installation](https://bitbucket.org/MAVERICLab/vcontact2/src/master/)
### installed on c3ddb node, conda activate vContact2, run on 20/May/2019
	vcontact --raw-proteins P1_1.faa --rel-mode ‘Diamond’ --proteins-fp P1.gene2genome.mappingfile.csv --db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /scratch/users/xiaoqiong/anaconda3/envs/vContact2/bin --output-dir VirSorted_Outputs
### getting the mapping quality adapted from bolduc script: Read2ReferenceMapper.py
	bamm filter -b {} --percentage_id {} --percentage_aln {} -o {}
### in the GOV2.0 paper

	BamM (https://github.com/ecogenomics/BamM) was used to remove reads that mapped at < 95% nucleotide identity to the contigs, bedtools genomecov (Quinlan and Hall, 2010) was used to determine how many positions across each genome were covered by reads, and custom Perl scripts were used to further filter out contigs without enough coverage across the length of the contig.
		usage: bamm filter -b BAMFILE [-o OUT_FOLDER]
		                   [--mapping_quality MAPPING_QUALITY]
		                   [--max_distance MAX_DISTANCE] [--length LENGTH]
		                   [--percentage_id PERCENTAGE_ID]
		                   [--percentage_aln PERCENTAGE_ALN] [--use_secondary]
		                   [--use_supplementary] [-v] [-h]
		Apply stringency filter to Bam file reads
		required arguments:
		  -b BAMFILE, --bamfile BAMFILE
		                        bam file to filter
		optional arguments:
		  -o OUT_FOLDER, --out_folder OUT_FOLDER
		                        write to this folder (output file has '_filtered.bam' suffix) (default: .)
		  --mapping_quality MAPPING_QUALITY
		                        mapping quality threshold
		  --max_distance MAX_DISTANCE
		                        maximum allowable edit distance from query to reference (default: 1000)
		  --length LENGTH       minimum allowable query length
		  --percentage_id PERCENTAGE_ID
		                        minimum base identity of mapped region (between 0 and 1)
		  --percentage_aln PERCENTAGE_ALN
		                        minimum fraction of read mapped (between 0 and 1)
		  --use_secondary       use reads marked with the secondary flag
		  --use_supplementary   use reads marked with the supplementary flag
		  -v, --invert_match    select unmapped reads
		  -h, --help            show this help message and exit

	/scratch/users/xiaoqiong/tools/MAVERICLab-vcontact2-0eca9dac02f7/vcontact/data/
	/scratch/users/xiaoqiong/anaconda3/envs/vContact2/lib/python3.7/site-packages/vcontact/data/ViralRefSeq-prokaryotes-v85.faa.gz
	'data/ViralRefSeq-prokaryotes-v88.protein2contig.csv'
	'data/ViralRefSeq-prokaryotes-v88.Merged-reference.csv'
	'data/ViralRefSeq-prokaryotes-v85.faa.gz',
	'data/ViralRefSeq-prokaryotes-v85.protein2contig.csv',
	'data/ViralRefSeq-prokaryotes-v85.ICTV-reference.csv',
	'data/ViralRefSeq-prokaryotes-v85.Merged-reference.csv',
	'data/ViralRefSeq-archaea-v85.faa.gz',
	'data/ViralRefSeq-archaea-v85.protein2contig.csv',
	'data/ViralRefSeq-archaea-v85.Merged-reference.csv'

virsorter 20190523 in c3ddb, https://github.com/simroux/VirSorter

	source activate virsorter
	wrapper_phage_contigs_sorter_iPlant.pl -f contigs.500bp.fasta --db 1 --wdir virsorter_output --ncpu 4 --data-dir scratch/users/xiaoqiong/tools/VirSorter/virsorter-data
	wrapper_phage_contigs_sorter_iPlant.pl -f scaffolds.fasta --db 1 --wdir virsorter_output --ncpu 20 --data-dir /scratch/tools/virsorter-data
	wrapper_phage_contigs_sorter_iPlant.pl -f scaffolds.fasta --db 1 --wdir virsorter_output --ncpu 20 --data-dir /scratch/users/xiaoqiong/tools/VirSorter/virsorter-data
virsorter 20190618 in c3ddb

	wrapper_phage_contigs_sorter_iPlant.pl -f combined.assembly.virsorter.fa --db 1 --wdir virsorter_output  --data-dir scratch/users/xiaoqiong/tools/VirSorter/virsorter-data
working directory

	cd /scratch/users/xiaoqiong/shingiek/kraken2_phages/filter_spades
	#!/bin/bash
	for x in *.fasta/
	do wrapper_phage_contigs_sorter_iPlant.pl -f ${x}contigs.500bp.fasta --db 1 --wdir ${x}virsorter_output --ncpu 40 --data-dir /scratch/users/xiaoqiong/tools/VirSorter/virsorter-data
	done
	can I use reference-based assembly?
	which tools to use? for phages?
	pseudomnas phage - 45.1 kb
### rnavirome 20190523

	samtools depth sorted.Bam > file.coverage
	bedtools genomecov
	how to visualize the samtools depth file
	need to check the RVDB database - if the virus genome is complete? or not?
### [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)

	By default, bedtools genomecov will compute a histogram of coverage for the genome file provided. The default output format is as follows:
	chromosome (or entire genome)
	depth of coverage from features in input file
	number of bases on chromosome (or genome) with depth equal to column 2.
	size of chromosome (or entire genome) in base pairs
	fraction of bases on chromosome (or entire genome) with depth equal to column 2.
	download the pseudomnas phage TL and map all the reads to the reference
	!!!be remember that all the bowtie2 input files should be in the exact correct reachable location!!!
### [How to plot coverage and depth statistics of a bam file](https://www.biostars.org/p/104063/)

To select the coverage for a particular chromosome (Chr#1 in my case)

	awk '$1 == 1 {print $0}' deduped_MA605.coverage > chr1_MA605.coverage

To select coverage from chr #2

	awk '$1 == 2 {print $0}' deduped_MA605.coverage > chr2_MA605.coverage

To plot the data in R this coverage file will need to be imported and the headers need to be added

	awk '$1 == "acc|GENBANK|LS992247.1|Human" {print $0}' 10.txt > 10_top1.txt

### 20190701 crass, http://ctskennerton.github.io/crass/Tutorial.html, https://github.com/ctSkennerton/crisprtools
	anaconda install crass
	crass -o crass_out -l 4 file1 file2    #this is the script to be used.
	crisprtools merge [-hs] [-o FILE] input.crispr{1,n}
	crisprtools filter [-h] [-o FILE] [-s INT] [-f INT] [-d INT] input.crispr
	Unknown option: -h
crass (1.0.1)
crass is a set of smal utilities for manipulating .crispr files
The .crispr file specification is a standard xml based format for describing CRISPRs
Type crass <subcommand> -h for help on each utility
Usage:	crass <subcommand> [options]
subcommand:  merge       combine multiple files
             help        display this message and exit
             extract     extract sequences in fasta
             filter      make new files based on parameters
             sanitise    change the IDs of elements
             stat        show statistics on some or all CRISPRs
             rm          remove a group from a .crispr file
### sequencing data synthetics
	Grinder: a versatile amplicon and shotgun sequence simulator
	This is the first tool to simulate amplicon datasets (e.g. 16S rRNA) widely used by microbial ecologists. Grinder can create sequence libraries with a specific community structure, α and β diversities and experimental biases (e.g. chimeras, gene copy number variation) for commonly used sequencing platforms.

### crispr assembled_viral_contigs
	crass-assembler --velvet -x crass.crispr -g NUM -s s1,s2,s3,s4,s5 -i DIR

### multiple alingnment
tools

	ClustalW
	MUSCLE -> developed by Robert Edgar
	mafft ->  Multiple alignment of a large number of short and highly similar sequences. Typical data size is up to ∼200,000 sequences × ∼5,000 sites (including gaps), but depends on similarity

muscle/mafft alignment and export as *fasta* format

	muscle -in seqs.fa -out seqs.afa
	mafft --thread -40 --auto in > out
	mafft --6merpair --thread -40 --keeplength --addfragments in.seq.fasta SARSCOV2.NC_045512.fasta > in.align.mafft.fasta
	mafft --parttree --thread -40 in > out ## align up to 1.4 million sequences, check out this paper: Large multiple sequence alignments with a root-to-leaf regressive method

muscle/mafft alignment and export to fasta as *CLUSTALW* format (more readable than FASTA)

	muscle -in seqs.fa -clw -out seqs.aln
	mafft --thread -40 --clustalout seq.fasta > align.fa  #https://mafft.cbrc.jp/alignment/software/
	mafft is doing in the interactive mode #https://towardsdatascience.com/
	how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6

count the number of sequences

	zcat my.fastq.gz | echo $((`wc -l`/4))

# STAR usage
## Creating a genome index

	STAR --runThreadN 40 \
	--runMode genomeGenerate \
	--genomeDir /data2/ZenNa/mice_data/STAR/mm39_index \
	--genomeFastaFiles /data2/ZenNa/mice_data/ref/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna \
	--sjdbGTFfile /data2/ZenNa/mice_data/ref/GCF_000001635.27/genomic.gtf


## Aligning reads

	mkdir STAR
	STAR --genomeDir /data2/ZenNa/mice_data/STAR/mm39_index/ \
	--runThreadN 40 \
	--readFilesCommand zcat \
	--readFilesIn /data2/ZenNa/mice_data/trimmed/SRR16633917_1.fastq.gz /data2/ZenNa/mice_data/trimmed/SRR16633917_2.fastq.gz \
	--outFileNamePrefix  SRR16633917 \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard  \
	--quantMode TranscriptomeSAM

## htseq-count 
	
	htseq-count -f bam -s no -t feature_type input.bam annotation.gtf > output.counts
	htseq-count -f bam -s no input.bam annotation.gtf > counts.txt





### taxonkit manual, https://wemp.app/posts/91cab928-7c98-4ea3-b54a-defaff169a08
# convert the accession id to the taxid
	blastdbcmd -db /data/database/NCBI_nt/nt -entry all -outfmt "%a %T"|pigz -c > nt.acc2taxid.txt.gz
### convert the kmer covereage to base coverage
	http://seqanswers.com/forums/showthread.php?t=1529
	http://seqanswers.com/forums/showthread.php?t=6887
	Cx=Ck*L/(L-k+1)
	where k is your kmer setting and L is your read length
	so if I used a kmer of 37 and an average read length of 50
	Cx=202*50/(50-37+1)=721X

## NCBI project download
[download sratools](https://github.com/ncbi/sra-tools)
[sratools wiki](https://github.com/ncbi/sra-tools/wiki)

	Limit download to 1000 reads from file SRR519926.
	fastq-dump -X 10000 --split-files SRR519926

[faster-dump wiki](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)

	faster-dump 

## construct phylogenetic tree

	https://docs.qiime2.org/2019.10/tutorials/phylogeny/


## fragment recruitment plot

	https://github.com/jianshu93/RecruitmentPlot_blast

## Others

compare if two files are identical or different

	cmp --silent $old $new && echo 'files are identical!' || echo "files are different!"

	&&: The AND Operator (&&). The second command will only execute if the first command has executed successfully i.e, its exit status is zero. 
	||: The OR Operator (||).  much like an 'else' statement in programming. The above operator allow you to execute second command only if the execution of first command fails, i.e., the exit status of first command is '1'.

## bcftools

output

	-O, --output-type	b|u|z|v
	Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion.

	
bcftools mpileup to do genotpye likelyhood, [format](https://en.wikipedia.org/wiki/Variant_Call_Format)

	bcftools mpileup --threads 40 -a FMT/AD,FMT/ADF,FMT/ADR,FMT/DP -Ou -d 8000 -A -f $ref $sorted.bam| bcftools call -Ov -A -M -c --threads 40 > snp/$f.raw.vcf

	bcftools mpileup -Ou -d 8000 -f $ref $sorted.bam | bcftools call --threads 40 --ploidy 1 -mv -Ov -o $out_dir/$root.$aligner.$caller.vcf\n"; 
	# the "v" in "-mv" specifies that the output file only contain variant sites; omit this to output all sites
	# mA -> output all sites

	AD: Read depth for each allele -> AD may not always sum to DP
	ADF: Read depth for each allele on the forward strand
	ADR: Read depth for each allele on the reverse strand
	AC: Total alternate allele count for all possible genotypes
	AN: Total number of alleles in all possible genotypes
	AF: allele frequency = AC/AN, If AF < 0.5, then AF is equal to MAF; rare variants generally has AF or MAF < 5 % (0.05)
	DP: Read depth
	v: uncompressed VCF
	u: uncompressed BCF
	z: compressed VCF 
	b: compressed BCF
	-A, --keep-alts
	-M, --keep-masked-ref, output sites where REF allele is N
	-c, --consensus-caller the original samtools/bcftools calling method (conflicts with -m)

merge different samples

	bcftools merge -Ov *.ecoli.raw.vcf.gz > merged.ecoli.raw.vcf.gz

reformat bcftools output to clean up the formats

	bcftools view -i 'QUAL>=20' merged.ecoli.raw.vcf.gz|bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %AD]\t%DP\t%QUAL\t%TYPE\n' -H -i 'TYPE="SNP"' > merged.ecoli.flt.snp.vcf


AD refers to the allele depth. AD reports the informative reads supporting each allele. .
What does informative reads mean?

DP
Generally, markers are retained with DP > 10 or DP > 5 to get high-quality genotypes. This value can be changed based on research applications. For example, in clinical research higher DP is desirable for filtering.

resources

	https://www.reneshbedre.com/blog/vcf-fields.html#an-total-allele-count


### instrain pipeline

	inStrain profile X.sam ref.fa -o sample.fa-vs-ref.IS -p 6 -g X.genes.fna -s X.maxbin2.stb

	inStrain compare -i N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.IS/ N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G2.IS/ -s .N5_271_010G1.maxbin2.stb -p 6 -o N5_271_010G1_scaffold_min1000.fa.IS.COMPARE


	parse_stb.py --reverse -f UHGG_reps/* -o UHGG.stb
	for genome in $(ls *.fna); do echo prodigal -i $genome -o ../UHGG_genes/$genome.genes -a ../UHGG_genes/$genome.gene.faa -d ../UHGG_genes/$genome.gene.fna -m -p single; done | parallel -j 6

	cat UHGG_genes/*.gene.fna > UHGG_reps.genes.fna

	cat UHGG_genes/*.gene.faa > UHGG_reps.genes.faa


	cat genomes-all_metadata.tsv | awk -F "\t" '{if ($17 == $1) print "curl ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_catalogue/" substr($18,0,13) "/" $18 "/genome/" $18 ".fna -o UHGG_reps/" $1 ".fna"}' | bash


## GNU parallel - master scripts
	cat ids.txt | parallel echo spades.py --meta -1 {}_R1.fq.gz -2 {}_R2.fq.gz -t 45 -o {}


## pangenome analysis
## roary plots

	https://github.com/sanger-pathogens/Roary/blob/master/contrib/roary_plots/roary_plots.ipynb

## IGV batch process

	https://github.com/igvteam/igv/wiki/Batch-commands
	https://software.broadinstitute.org/software/igv/batch
	https://github.com/stevekm/IGV-snapshot-automator
	


## Sources
* <http://gettinggeneticsdone.blogspot.com/2013/10/useful-linux-oneliners-for-bioinformatics.html#comments>
* <http://sed.sourceforge.net/sed1line.txt>
* <https://github.com/lh3/seqtk>
* <http://lh3lh3.users.sourceforge.net/biounix.shtml>
* <http://genomespot.blogspot.com/2013/08/a-selection-of-useful-bash-one-liners.html>
* <http://biowize.wordpress.com/2012/06/15/command-line-magic-for-your-gene-annotations/>
* <http://genomics-array.blogspot.com/2010/11/some-unixperl-oneliners-for.html>
* <http://bioexpressblog.wordpress.com/2013/04/05/split-multi-fasta-sequence-file/>
* <http://www.commandlinefu.com/>
* [biostarthandbook-Sep-2021]



