#!/bin/bash
# install 2 June 2020, on the singene workstation
# install 1 Sep 2021, on the old macbook laptop
# install 21 Sep 2023, on 

# Building identical conda environments

	conda list --explicit > spec-file.txt
	conda create --name myenv --file spec-file.txt

# installing in silent mode

	bash ~/anaconda.sh -b -p $HOME/anaconda
	eval "$(/home/kasm-user/anaconda/bin/conda shell.YOUR_SHELL_NAME hook)"
	conda init



# the order of channels should be in file $HOME/.condarc

	conda create --name meta python=3.10 -y
	conda install -n meta -c conda-forge awscli -y
	conda install -n meta -c bioconda bbmap -y
	conda install -n meta -c bioconda bioawk -y
	conda install -n meta -c conda-forge -c bioconda quast -y   #python version 3.6, if not pip install quast
	conda install -n meta -c bioconda fastqc -y
	conda install -n meta -c bioconda -c conda-forge multiqc -y
	conda install -n meta -c bioconda seqtk -y
	conda install -n meta -c bioconda diamond -y
	conda install -n meta -c bioconda bowtie2 -y
	conda install -n meta -c bioconda samtools -y
	conda install -n meta -c bioconda prodigal -y
	conda install -n meta -c bioconda megahit -y
	conda install -n meta -c bioconda spades -y
	conda install -n meta -c bioconda novoalign -y
	conda install -n meta -c bioconda stamp -y
	conda install -n meta -c anaconda biopython -y
	conda install -n meta -c bioconda blast -y
	conda install -n meta -c bioconda trnascan-se -y
	conda install -n meta -c bioconda lefse -y
	conda install -n meta -c bioconda cd-hit -y
	conda install -n meta -c bioconda bedtools -y
	conda install -n meta -c ursky metabat2 -y
	conda install -n meta -c bioconda hmmer -y # no
	conda install -n meta -c bioconda checkm-gen-yome -y
	conda install -n meta -c bioconda muscle -y
	conda install -n meta -c bioconda eggnog-mapper -y
	conda install -n meta -c bioconda stamp -y
	conda install -n meta -c bioconda igv -y
	conda install -n meta -c bioconda pkg -y
	conda install -n meta -c bioconda bioperl -y
	conda install -n meta -c r rpy2 -y
	conda install -n meta -c bioconda mafft -y
	conda install -n meta -c bioconda fasttree -y
	conda install -n meta -c bioconda taxonkit -y
	conda install -n meta -c bioconda kraken2 -y
	conda install -n meta -c bioconda kraken-biom -y
	conda install -n meta -c bioconda star -y
	conda install -c bioconda nanoplot -y

### conda install numpy, pandas etc.

	conda install -c anaconda jupyter -y
	conda install -c anaconda numpy -y
	conda install -c anaconda pandas -y
	conda install -c conda-forge biopython -y
	conda install -c conda-forge xlrd -y
	conda install -c conda-forge scipy -y
	conda install -c anaconda seaborn -y
	conda install -c conda-forge matplotlib -y


# if samtools has conflict problem
# conda install -c bioconda samtools=1.9=h8ee4bcc_1
# conda install -c bioconda samtools openssl=1.0 #downgrade the packages

	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda install samtools==1.11 #samtools version 1.11 is stable

### GNU parallel installation

	sudo apt install parallel 

### phylophlan2 install
	
	conda env create fasnicar/ppa2
	sudo apt install mercurial
	hg clone https://bitbucket.org/nsegata/phylophlan

	conda activate ppa2
	conda install -c anaconda biopython -y
	conda install -c bioconda dendropy -y
	conda install -c bioconda fasttree -y
	conda install -c bioconda muscle -y
	conda install -c bioconda trimal -y
	conda install -c bioconda diamond -y
	conda install -c bioconda blast -y

### phylophlan3 install
	
	conda create -n "phylophlan" -c bioconda phylophlan=3.0
	conda install -c bioconda raxml -y
	conda install -c bioconda iqtree -y

### plasflow install, python=3.5
  
  conda create --name plasflow python=3.5
  conda install plasflow -c bioconda
  source activate plasflow
  source deactivate

## error message libreadline.so.6 issue in ubuntu 18.04 when running plasflow
	
	https://askubuntu.com/questions/1168787/libreadline-so-6-issue-in-ubuntu-18-04
	sudo apt-get install libreadline-dev 
  cd /lib/x86_64-linux-gnu/
	sudo ln -s libreadline.so.7.0 libreadline.so.6 # created  symlink so can continue work with libreadline.so.6.

### anvio install: http://merenlab.org/2016/06/26/installation-v2/
  
  conda create -y --name anvio-6.2 python=3.6
  conda activate anvio-6.2
  conda install -y -c conda-forge -c bioconda anvio==6.2

### crass install: requirment for xerces-c version is 3.1.4
	
	conda install -c usgs-astrogeology xerces-c
	conda install -c bioconda crass

###  crass install
	conda create --name crass python=3.5
	conda install -c bioconda crass -y
	conda activate crass

=============install crisprtools==================
#### dependency, instal libcrispr, https://github.com/ctSkennerton/libcrispr
	
	git clone https://github.com/ctSkennerton/libcrispr.git
	sudo apt-get install libxerces-c3-dev #install libexrces version 3 instead of version 2

	./autogen.sh
	./configure
	[sudo] make install

	export LD_LIBRARY_PATH="/usr/local/lib" #add them permanantly in the ~/.bashrc profile, https://github.com/ctSkennerton/crass/issues/64
	export LD_LIBRARY_PATH="/usr/local/lib:/opt/lib"

#### install crisprtools, http://ctskennerton.github.io/crisprtools/
	[./autogen.sh]
	./configure
	make
### now crisprtools is succesfully installed in the /usr/local/bin
=============install crisprtools==================

## Mac System

============install homebrew on the macOS Catalina==================
# https://docs.brew.sh/Installation
# https://brew.sh/

	/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

## zsh set up, https://juejin.im/post/6844903939121348616
	
	brew install autojump
	git clone git://github.com/zsh-users/zsh-autosuggestions $ZSH_CUSTOM/plugins/zsh-autosuggestions
	git clone git://github.com/zsh-users/zsh-syntax-highlighting $ZSH_CUSTOM/plugins/zsh-syntax-highlighting
	vim ~/.zshrc
	source ~/.zshrc

============install Nextstrain on the macOS Catalina==================
# https://docs.nextstrain.org/en/latest/guides/install/cli-install.html
# install Docker


## Linux System
============install zsh on the linux ==================
	https://novnan.github.io/Linux/install-zsh-shell-ubuntu-18-04/
	git clone https://github.com/zsh-users/zsh-autosuggestions ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-autosuggestions
	git clone https://github.com/zsh-users/zsh-syntax-highlighting ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting

## install ranger on linux/mac system


###install sublime on linux system, https://blog.csdn.net/u012707739/article/details/78148976
#安装GPG
	
	wget -qO - https://download.sublimetext.com/sublimehq-pub.gpg | sudo apt-key add -

#确保apt被设置为https源
	
	sudo apt-get install apt-transport-https

#选择稳定版本
	
	echo "deb https://download.sublimetext.com/ apt/stable/" | sudo tee /etc/apt/sources.list.d/sublime-text.list

#安装sublime-text
	
	sudo apt-get update
	sudo apt-get install sublime-text


# To fix the CONCOCT endless warning messages in metaWRAP=1.2+, run
 
 	conda install -y blas=2.5=mkl
 	conda create --name metawrap-1.3.2 --channel ursky --channel bioconda --channel conda-forge metawrap-mg=1.3.2
 	conda activate metawrap-1.3.2


# conda tutorial and snakemake
	
	https://www.biostars.org/p/335903/

# change the ownership of the harddisk
	
	https://askubuntu.com/questions/527304/how-to-change-permission-of-a-drive-in-an-external-hard-disk
	sudo chown -R $USER:$USER /media/$USER/drivename
	sudo su -   https://askubuntu.com/questions/376199/sudo-su-vs-sudo-i-vs-sudo-bin-bash-when-does-it-matter-which-is-used
	fdisk -l
	nvidia-smi
	watch cat /proc/mdstat

# install R
	
	sudo apt install r-base r-base-dev -y

# copy the identical conda environment
	
	conda list --explicit > conda-spec-file.txt
	conda create --name amr --file conda-spec-file.txt
	conda install --name amr --file conda-spec-file.txt


	

