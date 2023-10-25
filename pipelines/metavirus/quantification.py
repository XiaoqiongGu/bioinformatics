
import os
from metavirus.base import Moduler


class bowtie2_module(Moduler):
    
    def __init__(self, output_dir):
        Moduler.__init__(self, output_dir)
        self.bowtie_result_dir = os.path.join(self.output_dir, 'bowtie_result')
        self.create_output(self.bowtie_result_dir)
        self.bowtie_database_dir = os.path.join(self.output_dir, 'bowtie_database')
        self.create_output(self.bowtie_database_dir)
        
        
    def build_database(self, database_file_list, database_target_dir=None):
        
        database_file_num = len(database_file_list)
        
        if database_file_num == 1:
            assembly_file = database_file_list[0]
            if database_target_dir == None:
                database_dir = os.path.join(self.bowtie_database_dir, 'cross_database')
            else:
                database_dir = database_target_dir
            
            assembly_size = os.path.getsize(assembly_file)/1024**3
        
            if assembly_size > 3:
                build_databse_op = 'bowtie2-build -f ' + assembly_file + ' ' + database_dir + ' -p 40 --large-index'
            else:
                build_databse_op = 'bowtie2-build -f ' + assembly_file + ' ' + database_dir + ' -p 40'
            os.system(build_databse_op)
        else:
            print('Still need working on')
            
        return database_dir
    
    
    def run(self, input_files, database_dir, bowtie_target_dir=None):
        
        if bowtie_target_dir == None:
            bowtie_output_dir = os.path.join(self.bowtie_result_dir, database_dir.split('/')[-1])
        else:
            bowtie_output_dir = bowtie_target_dir
        self.create_output(bowtie_output_dir)
        
        for name in input_files:
            R1_input_dir = name[0]
            R2_input_dir = name[1]
            
            R1_file = R1_input_dir.split('/')[-1]
            subject_name = R1_file.split('R1')[0]
            
            sam_file = os.path.join(bowtie_output_dir, subject_name+'.sam')
            
            bowtie_op = 'bowtie2 -x ' + database_dir + ' -1 ' + R1_input_dir + ' -2 ' + R2_input_dir + ' -S ' + sam_file + ' -p 40'
            os.system(bowtie_op)
            
            bam_file = os.path.join(bowtie_output_dir, subject_name+'.bam')
            os.system('samtools view -bS ' + sam_file + ' > ' + bam_file)
    
            sort_file = os.path.join(bowtie_output_dir, subject_name+'.sorted.bam')
            os.system('samtools sort ' + bam_file + ' -o ' + sort_file)
    
            os.system('samtools index ' + sort_file)
    
            os.system('samtools depth ' + sort_file)
    
            idx_file = os.path.join(bowtie_output_dir, subject_name+'.idxstats.txt')
            os.system('samtools idxstats ' + sort_file + ' > ' + idx_file)
    
            genome_file = os.path.join(bowtie_output_dir, subject_name+'.genome')
            os.system('bedtools genomecov -ibam ' + sort_file + ' -g ' + genome_file)            
        
        idx_star = os.path.join(self.output_dir, '*.idxstats.txt')
        absolute_count_op = 'get_count_table.py ' + idx_star + ' > ' + bowtie_output_dir + 'absolute_count.txt'
        os.system(absolute_count_op)
        
        abs2RPKM_op = 'absolutereads2RPKM.py ' + bowtie_output_dir + 'absolute_count.txt ' + bowtie_output_dir + 'RPKM.txt'
        os.system(abs2RPKM_op)
