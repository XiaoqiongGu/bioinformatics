
import os
from metavirus import base
from metavirus import trimming
from metavirus import annotation
from metavirus import assembly
from metavirus import quantification

class pipeline_module(base.Moduler):
    
    def __init__(self, 
                 project_name, 
                 input_dir, 
                 project_dir=None,
                 
                 adapter_init=False,
                 adapter_dir=None,
                 adapter_ref='/data1/tools/bbmap/resources/adapters.fa',
                 adapter_quality=10,
                 
                 phix_init=False,
                 phix_dir=None,
                 phix_ref='/data1/tools/bbmap/resources/phix174_ill.ref.fa.gz',
                 phix_quality=10,
                 
                 primer_init=False,
                 primer_dir=None,
                 primer_ref='GTTTCCCAGTCACGATA',
                 primer_quality=10,

                 human_init=False,
                 human_dir=None,
                 human_ref='/data1/tools/bbmap/resources/human',
                 
                 fastq_init=False,
                 fastq_dir=None,
                 
                 megahit_init=False,
                 megahit_dir=None,
                 megahit_contig_length=1000,
                 
                 metaspade_init=False,
                 metaspade_dir=None,
                 metaspade_contig_length=1000,
                 
                 diamond_init=False,
                 diamond_dir=None,
                 diamond_database='/data1/database/Diamondnr/nr',
                 e_value=0.001,
                 blast2rma_ref='/data1/tools/megan/tools/blast2rma',
                 acc2tax='/data/database/megan/old_mapping_file/prot_acc2tax-Jun2018X1.abin',
                 acc2eggnog='/data/database/megan/old_mapping_file/acc2eggnog-Jul2019X.abin',
                 acc2interpro='/data/database/megan/old_mapping_file/acc2interpro-Jul2019X.abin',
                 acc2kegg='/data/database/megan/old_mapping_file/acc2kegg-Jul2019X.abin',
                 acc2seed='/data/database/megan/old_mapping_file/acc2seed-May2015XX.abin',
                 collapse_rank='Genus',
                 select_rank='Genus',
                 tax_format='readName_to_taxonName',
                 egg_format='readName_to_eggnogPath',
                 seed_format='readName_to_seedPath',
                 kegg_format='readName_to_keggPath',
                 
                 bowtie2_init=False,
                 bowtie2_dir=None
                 ):
        
        self.project_name = project_name
        self.input_dir = input_dir
        if project_dir == None:
            self.project_dir = os.path.join(self.input_dir.split(self.input_dir.split('/')[-1])[0], self.project_name)
        else:
            self.project_dir = project_dir
        self.create_output(self.project_dir)
        
        self.readin_module = base.readin_module(self.project_dir)
        
        if adapter_init:
            if adapter_dir == None:
                self.adapter_dir = os.path.join(self.project_dir, 'adapter')
            else:
                self.adapter_dir = adapter_dir
            self.adapter_module = trimming.adapter_module(self.adapter_dir, adapter_ref, adapter_quality)

        if phix_init:
            if phix_dir == None:
                self.phix_dir = os.path.join(self.project_dir, 'phix')
            else:
                self.phix_dir = phix_dir
            self.phix_module = trimming.phix_module(self.phix_dir, phix_ref, phix_quality)
        
        if primer_init:
            if primer_dir == None:
                self.primer_dir = os.path.join(self.project_dir, 'primer')
            else:
                self.primer_dir = primer_dir
            self.primer_module = trimming.primer_module(self.primer_dir, primer_ref, primer_quality)

        if human_init:
            if human_dir == None:
                self.human_dir = os.path.join(self.project_dir, 'human')
            else:
                self.human_dir = human_dir
            self.human_module = trimming.human_module(self.human_dir, human_ref)
        
        if fastq_init:
            if fastq_dir == None:
                self.fastq_dir = os.path.join(self.project_dir, 'fastq')
            else:
                self.fastq_dir = fastq_dir
            self.fastq_module = trimming.fastq_module(self.fastq_dir)
        
        if megahit_init:
            if megahit_dir == None:
                self.megahit_dir = os.path.join(self.project_dir, 'megahit')
            else:
                self.megahit_dir = megahit_dir
            self.megahit_module = assembly.megahit_module(self.megahit_dir, megahit_contig_length)
            
        if metaspade_init:
            if metaspade_dir == None:
                self.metaspade_dir = os.path.join(self.project_dir, 'metaspade')
            else:
                self.metaspade_dir = metaspade_dir
            self.metaspade_module = assembly.metaspade_module(self.metaspade_dir, metaspade_contig_length)
            
        if diamond_init:
            if diamond_dir == None:
                self.diamond_dir = os.path.join(self.project_dir, 'diamond')
            else:
                self.diamond_dir = diamond_dir
            self.diamond_module = annotation.diamond_module(self.diamond_dir, diamond_database, e_value, blast2rma_ref,
                                                            acc2tax, acc2eggnog, acc2interpro, acc2kegg, acc2seed,
                                                            collapse_rank, select_rank, tax_format, egg_format, seed_format,
                                                            kegg_format)
            
        if bowtie2_init:
            if bowtie2_dir == None:
                self.bowtie2_dir = os.path.join(self.project_dir, 'bowtie2')
            else:
                self.bowtie2_dir = bowtie2_dir
            self.bowtie2_module = quantification.bowtie2_module(self.bowtie2_dir)
        



    
class virus_NO1(pipeline_module):
    
    def __init__(self, 
                 project_name,
                 input_dir,
                 adapter_init=True,
                 adapter_quality=20,
                 phix_init=True,
                 phix_quality=20,
                 primer_init=True,
                 primer_quality=20,
                 human_init=True,
                 fastq_init=True,
                 megahit_init=True,
                 megahit_contig_length=1000,
                 metaspade_init=True,
                 diamond_init=True,
                 bowtie2_init=True
                 ):
        
        pipeline_module.__init__(self, project_name, input_dir, 
                                 adapter_init=adapter_init,
                                 adapter_quality=adapter_quality,
                                 phix_init=phix_init,
                                 phix_quality=phix_quality,
                                 primer_init=primer_init,
                                 primer_quality=primer_quality,
                                 human_init=human_init,
                                 fastq_init=fastq_init,
                                 megahit_init=megahit_init,
                                 megahit_contig_length=megahit_contig_length,
                                 metaspade_init=metaspade_init,
                                 diamond_init=diamond_init,
                                 bowtie2_init=bowtie2_init)
        
        
    
    def forward(self):
        
        input_files = self.readin_module.run(self.input_dir)
        
        '''trimming'''
        print('Start adapter')
        adapter_files = self.adapter_module.run(input_files)
        
        print('Start phix')
        phix_files = self.phix_module.run(adapter_files)
        
        '''assembly'''
        print('start megahit cross assembly')
        megahit_cross = self.megahit_module.run_cross(phix_files)
        megahit_indi = self.megahit_module.run_indiviudal(phix_files)
        
        '''annotation'''
        self.diamond_module.run(megahit_indi)
        
        '''quantification'''
        cross_databse = self.bowtie2_module.build_database(megahit_cross)
        self.bowtie2_module.run(phix_files, cross_databse)
        
        
    
class virus_NO2(pipeline_module):
    
    def __init__(self, 
                 project_name,
                 input_dir,
                 project_dir,
                 adapter_init=True,
                 adapter_quality=20,
                 phix_init=True,
                 phix_quality=20,
                 primer_init=True,
                 primer_quality=20,
                 human_init=True,
                 fastq_init=True,
                 megahit_init=True,
                 metaspade_init=True,
                 diamond_init=True,
                 bowtie2_init=True
                 ):
        
        pipeline_module.__init__(self, project_name, input_dir,
                                 project_dir=project_dir,
                                 adapter_init=adapter_init,
                                 adapter_quality=adapter_quality,
                                 phix_init=phix_init,
                                 phix_quality=phix_quality,
                                 primer_init=primer_init,
                                 primer_quality=primer_quality,
                                 human_init=human_init,
                                 fastq_init=fastq_init,
                                 megahit_init=megahit_init,
                                 metaspade_init=metaspade_init,
                                 diamond_init=diamond_init,
                                 bowtie2_init=bowtie2_init)
        
        
    
    def forward(self):
        
        input_files = self.readin_module.run(self.input_dir)
        
        '''trimming'''
        print('Start adapter')
        adapter_files = self.adapter_module.run(input_files)
        
        print('Start phix')
        phix_files = self.phix_module.run(adapter_files)
        
        '''assembly'''
        print('start megahit cross assembly')
        megahit_cross = self.megahit_module.run_cross(phix_files)
        megahit_indi = self.megahit_module.run_indiviudal(phix_files)
        
        '''annotation'''
        self.diamond_module.run(megahit_indi)
        
        '''quantification'''
        cross_databse = self.bowtie2_module.build_database(megahit_cross)
        self.bowtie2_module.run(phix_files, cross_databse)  
        
        
        
