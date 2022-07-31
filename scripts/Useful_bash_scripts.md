# Useful bash scripts
[bioinformatics bash one-liners](https://github.com/stephenturner/oneliners)
### talk between two systems
	nc -l 3333 workstation1
	nc $IP1 3333 workstation2
### how to kill a running nohup
	ps -ef | grep (the command name)
	kill PID
	kill -9 PID (force to kill)
### to check the number of cpus
	cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l
	grep -c ^processor /proc/cpuinfo
	nproc --all
	lscpu -p [linux]
	sysctl [MacOS]
### to count the number of columns in files with different column number
	awk '{print NF}' file.txt| sort -nu | tail -n 1 #NF is the built-in variable
### count the number of reads in fastq files
	echo $(cat yourfile.fastq|wc -l)/4|bc  
 	echo $(zcat yourfile.fastq.gz|wc -l)/4|bc
 	echo $(( $(wc -l < $filename) / 4 )) 

### translate \n with comma in files
	cat file.txt |tr '\n' ',' > file.txt
	cat file.txt |tr '\n' $'\t' > file.txt
	cat negativeid.txt |tr '\n' ' ' > 1.txt
### nohup usage
	nohup command > out.file 2>&1 &
	2: st
### how to tar zip files in multiple threads in linux
	sudo apt install pigz
	tar -zcvf -test.txt|pigz > test.tar.gz
### how to tar zip files
	nohup tar -zcfv trimmed.adapter.tar.gz adapter&
### how to unzip tar.gz files
	tar -zxfv bigfile.tar.gz -C /folder/subfolder/
	z zipped file
	x extract
	c create
	f use the following tar archive for the operation
	v verbose
	-C to extract the files into a specified directory
	
### show the size of hard disk usage of all the directory
	df -h
### show all the sub-directory hard disk size and also the total size
	du -h --max-depth=1 your_user_folder/
### show the total size of current folder sizes (disk usagetar)
	du -sh folder_path

### nano usage
	cut the entire line, ctrl+k

### mkdir the list of directories according to the file list
	xargs mkdir <list.txt
>xargs simply "flattens" your text file by replacing newlines with spaces, thereby invoking mkdir with a long list of arguments containing all your directory names at once instead of one at a time.

### Download NCBInr/nt database
	wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz'
	cat nr.*.tar.gz | tar -zxvi -f - -C .
	
	wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz'

### How to create a sym link
	mv ../ncbi .
	cd ..
	ln -s ~/proc/ncbi/ ncbi
	ln -s source_file myfile

### create the alias for simplify the command line
	alias script='cd /home/xiaoqiong/script'
	alias meta='module add c3ddb/miniconda;source activate meta'

### difference between ctrl+Z vs ctrl+ C
	ctrl+ Z: hang down
	ctrl+ C: kill the programme

### Use filename as the column name
    for file in ./file_*.txt
    do
    sed -i "1i\
    rows $file" "$file"
    done
### for what purpose?
	for x in *.txt ; do sed -i "1i\contigs\t$x" "$x";done
### use file name to fill in a whole column
	nawk '{a=FILENAME;}{print a"\t"$0}' yourFilename > yourFilename.bk && mv yourFilename.bk yourFilename
### change comma delimited to tab
	for file in *.txt; do cat $file | tr '[,]' '[\t]' > new_*.txt; done
### extract lines based on line numbers
	grep -n "keyword" contigs.faa [to show the line number which containing keyword]
	sed -n x,yp contigs.faa > contig1.fa [extract line x to line y]
	
	grep -H "keyword" *.txt #can show the list of file containing keyword as also show the file names
	grep -w "0$" test.txt #can extract the lines that only containing 0 in the end


### Merge multiple files by common field
https://stackoverflow.com/questions/13710876/merge-multiple-files-by-common-field-unix
##
	cut -f 1 file1 > delim   ## use first column as delimiter
	i=0
	for file in file*
	do
        i=$(($i+1))  ## for adding count to distinguish files from original ones
        cut -f 2 $file > ${file}__${i}.temp
	done
	paste -d\\t delim *.temp > output

### paste two files
	paste -d "\0" file1 file2 > file2

### Extract fasta sequences based on sequence ID
	!xargs samtools faidx final.contigs.fa < Id.txt > Id.contigs.fa

### sort usage
[sort by numerical values](https://stackoverflow.com/questions/4856030/how-to-sort-a-file-based-on-its-numerical-values-for-a-field)
	
	sort -u  # increasing order
	sort -ur # decreasing order
	sort -V, --version-sort natural sort of (version) numbers within text


### find and replace the text in files
	sed -i 's/old-text/new-text/g' input.txt
	in MacOS `sed -i '' 's/original/new/g' test.txt`
### find and replace the text in specific line in files
	sed -i '5s/old-text/new-text/' input.txt #do it in place, also specify it is in which line [now is 5 line]

# sed usage

delete lines

	sed '7,9d' in > out 	# from lines 7 to 9
	sed '2d' in > out 		# delete line 2
	sed '/oops/d' in > out  # remove the lines containing “oops”
	sed '/^$/d' in > out	# remove all empty lines
	sed '$d' in > out		# remove the last line

delete the first character each line in a file

	sed 's/^.//' in > out 

add prefix/suffix each line in a file

	sed -i -e 's/^/prefix/' file
	sed -i -e 's/$/suffix/' file
	sed allows in place editing via the -i parameter, or
	sed -e 's/$/suffix/' in > out 
	^ means beginning of the line, $ means end of the line
	sed 's/^X$/d' file #delete all the x

add file name to every sequence header
	
	for f in *.fna; do sed -i "" "s/^>/>${f%.*}_/" "$f"; done
	## note in the mac system, need to add ""-> https://blog.csdn.net/stpeace/article/details/106110704

when trying to copy the data path, be careful to use '\' in front of '/' in order to let '/' have the meaning

### separate taxon file
	echo '#Kingdom#Phylum#Class#Order#Family#Genus#Species' | tr '#' '\t' > taxonomy-table.tsv
	cut -f 1-2 path-to-your-taxonomy-file.tsv | tr ';' '\t' | tail -n +2 >> taxonomy-table.tsv

### tail usage
select all rows except the first row

	tail -n+2 file.txt

### extract fastq sequences based on sequence ID (fastq)
	seqtk subseq in.fq name.lst > out.fq


### awk usage

parse each field of document into each line

	awk '{ print }' file.txt
	awk '{ print $1 }' file.txt #print the first columns
	awk '{ /[a-z]/ print }' file.txt #show every line containing characters
	awk '{ /[0-9]/ print }' file.txt #show every line containing number
	awk '{ /^[0-9]/ print }' file.txt #show every line starting with the number
	awk '{ /[0-9]$/ print }' file.txt #show every line ending with the number
	awk '{ if ($1 ~ /123/) print }' file.txt #show if the first column equal to 123
	awk '{ if ($2 ~/[0-9]/) print }' file.txt #show if the second column containing the numbers

	awk '{ if ($3 ~ /100/) print }'
	awk ' NR==11 && NR<= 17{sum +=$3; if(NR == 17){exit}} END {print sum/7}' $FILE


sum certain numbers in a column using awk. I would like to sum just column 3 of the "smiths" to get a total of 212. 
	
	awk -F '|' '$1 ~ /smiths/ {sum += $3} END {print sum}' inputfilename
	awk '{sum+=$3} END { print "Average = ",sum/NR}'

	
	 -F":" '{print $2 }' file.txt #print the second column based on ":" delimiter
	awk -F":" 'NR==1,NR==10 {print $1}' file.txt #print the first column, between row 1 to row 10 based on ":" delimiter, NR is the built in variables

	awk '$4 > 9' #grep the fourth column which more than 9
	awk '$1==$7 && $6==$12' file #two conditions
	awk 'NR%2==0' file > outfile #grep even lines from a txt file
	
	FS: Input field separator variable
	OFS: Output Field Separator Variable
	
	awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
	
	awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget ",ftpdir,file}' ftpdirpaths > ftpfilepaths.download
	
	awk 'BEGIN{FS=OFS=" ";}{print "wget",$0}' ftpfilepaths > ftpfilepaths.download
	
	awk '{sum+=$1} END {print sum}'


### ls alias - to add color after ls -lh
	alias ls='ls --color=auto'

### open a new R session
	open -n /Applications/RStudio.app

### new bash scripts, chmod a+x script.sh
	#!/bin/bash
	NAME=$1
	echo "your name is $NAME"

$0 is the name of the script, $1 is the first argument
## extract seqs based on seqids
	WPB073-TGCTACAT_S1_L006.kraken.id.filter
	seqtk subseq in.fq name.lst > out.fq
	Note: Use 'samtools faidx' if only a few regions are intended   #this is not good to use
	
	filterbyname.sh in=reads.fq out=filtered.fq names=names.txt #from bbmap
	filterbyname.sh in=WPB073-TGCTACAT_S1_L006_R1_001.fastq.gz in2=WPB073-TGCTACAT_S1_L006_R2_001.fastq.gz out=WPB073-TGCTACAT_S1_R1.filtered.fastq out2=WPB073-TGCTACAT_S1_R2.filtered.fastq names=WPB073-TGCTACAT_S1_L006.kraken.id.filter


# Linearizing the complete fasta file
	while read line;do if [ "${line:0:1}" == ">" ]; then echo -e "\n"$line; else echo $line | tr -d '\n' ; fi; done < input.fasta > output.fasta #too slow
	
	sed -e 's/\(^>.*$\)/#\1#/' scaffolds.fasta | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d'

	bioawk -c fastx '{print ">"$name; print $seq}'   #useful link https://www.biostars.org/p/363676/

# split usage
	split -d -b 200M -l 500 file.txt log
	-b file sizes
	-l line numbers
	-d rename the splited file with log1, log2, log3, ...., logn
	split -l 1000

## sum up specific column numbers
	awk '{sum+=$1;}END{print $1}' infile > outfile
	
	awk '{for (i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}' file
	{for (i=1;i<=NF;i++)         Set field to 1 and increment through
	$i=(a[i]+=$i)                Set the field to the sum + the value in field
	END{print}                   Print the last line which now contains the sums
	
	tail -n+2 | awk '{for (i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}' kraken_table.txt|less
	
	numsum -c FILENAME [function]
	-c   ---    Print out the sum of each column.


1. [How to List all files in a directory and export result to a text file? ](https://stackoverflow.com/questions/14314947/how-to-list-all-files-in-a-directory-and-export-result-to-a-text-file)

2. [Find a file matching with certain pattern and giving that file name as value to a variable in shell script?](https://unix.stackexchange.com/questions/361655/find-a-file-matching-with-certain-pattern-and-giving-that-file-name-as-value-to)

	`find /data/meta_augmentin/meta_raw_data/combined -type f -iname "*_1.fq.gz" -printf '%p\n'|sort > test`



### the difference between bash and python
	bash shell for small stuffs vs python
	you can either run bash inside python or run python inside bash 

### alias usage
1. run 'git status' without typing so many letters
2. alias rm='rm -i'

### ls vs find
	folder_name
	find . list all the things and recursive things in the folder

### xargs
	cat file.txt|xargs mkdir
	cat file.txt|xargs touch
	find | xargs echo

### grep
	grep "" --color -n test.txt
	grep "\\$" test.txt #regular experssion  
	ls -lah | grep -v "Jul 28"

### man usage
	man awk 
	man sed 

### sed
	cat file.txt | sed 's/matcher/replacement/'
	cat file.txt | sed 's/$/.txt/'

### write a function
	$1: means the first argument
	$@: means all the argument
	Use mcd to create a directory and cd to it simultaneously:
	function mcd(){
		mkdir "$1"
		cd "$1"
	}

	function printme() {echo hi}

	function printme(){
	file=$1
	echo hi the length of $file is:
	wc -l $file
	}

	function rm(){
		/bin/rm $@
	}

	unset printme


### chmod
	chmod u+x file.txt

### echo
	echo -n -> 不会产生新的一行