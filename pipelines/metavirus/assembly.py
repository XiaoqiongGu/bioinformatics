
import os
from metavirus.base import Moduler


class megahit_module(Moduler):
    
    def __init__(self, output_dir, contig_length):
        Moduler.__init__(self, output_dir)
        self.contig_length = contig_length
        self.create_output()
        
    def run_cross(self, input_files):

        output_cross_dir = os.path.join(self.output_dir, 'cross')
        
        R1_combine = ''
        R2_combine = ''
        for name in input_files:
            
            R1_input_dir = name[0]
            R2_input_dir = name[1]
            
            R1_combine = R1_combine+','+R1_input_dir
            R2_combine = R2_combine+','+R2_input_dir
              
        '''megahit operation'''
        R1_combine = R1_combine[1:]
        R2_combine = R2_combine[1:]

        megahit_op = 'megahit -o '+output_cross_dir+' --min-contig-len '+str(self.contig_length)+\
        ' -1 '+R1_combine+' -2 '+R2_combine
        #print(megahit_op)
        os.system(megahit_op)
        
        '''output file change'''
        original_name = os.path.join(output_cross_dir,'final.contigs.fa')
        changed_name = os.path.join(output_cross_dir,'cross.contigs.fasta')
        self.rename(original_name, changed_name, 'cross_')
        
        output_report_dir = os.path.join(output_cross_dir, 'cross_contig_report')
        self.create_output(output_report_dir)
        quast_op = 'quast.py '+changed_name+' -o '+output_report_dir
        os.system(quast_op)
        
        return [changed_name]
       
        
        
    def run_indiviudal(self, input_files):
        
        output_file_list = []
        
        for name in input_files:
            
            R1_input_dir = name[0]
            R2_input_dir = name[1]
            
            R1_file = R1_input_dir.split('/')[-1]
            subject_name = R1_file.split('R1')[0]
                
            output_sample_dir = os.path.join(self.output_dir, subject_name+'contigs')

            megahit_op = 'megahit -o '+output_sample_dir+' --min-contig-len '+str(self.contig_length)+\
            ' -1 '+R1_input_dir+' -2 '+R2_input_dir
            os.system(megahit_op)
            
            original_name = os.path.join(output_sample_dir,'final.contigs.fa')
            changed_name = os.path.join(output_sample_dir, subject_name+'contigs.fasta')
            self.rename(original_name, changed_name, subject_name)
            
            output_report_dir = os.path.join(output_sample_dir, subject_name+'contig_report')
            self.create_output(output_report_dir)
            quast_op = 'quast.py '+changed_name+' -o '+output_report_dir
            os.system(quast_op)
            
            output_file_list.append(changed_name)
            
        return output_file_list
            

         
class metaspade_module(Moduler):
    
    def __init__(self, output_dir, contig_length):
        Moduler.__init__(self, output_dir)
        self.contig_length = contig_length
        
        
    def run_cross(self, input_files):
        
        output_cross_dir = os.path.join(self.output_dir, 'cross')
        self.create_output(output_cross_dir)
        
        R1_combine_op = 'cat '
        R2_combine_op = 'cat '
        
        R1_combine_file = os.path.join(output_cross_dir, 'R1_combined_file.fastq.gz')
        R2_combine_file = os.path.join(output_cross_dir, 'R2_combined_file.fastq.gz')
        
        size_sum = 0
        
        for name in input_files:
            
            R1_input_dir = name[0]
            R2_input_dir = name[1]
            
            R1_combine_op = R1_combine_op+R1_input_dir+' '
            R2_combine_op = R2_combine_op+R2_input_dir+' '
            
            size_sum = size_sum+os.path.getsize(R1_input_dir)+os.path.getsize(R2_input_dir)
        
        
        R1_combine_op = R1_combine_op+'> '+R1_combine_file
        R2_combine_op = R2_combine_op+'> '+R2_combine_file
        
        if size_sum/1024**3 < 20:
            os.system(R1_combine_op)
            os.system(R2_combine_op)
        else:
            raise ValueError('Combined file too large')
        
        metaspade_op = 'spades.py  --meta -t 70 -m 400'+' -o '+output_cross_dir+' -1 '+R1_combine_file+' -2 '+R2_combine_file
        os.system(metaspade_op)
        
        metaspade_output = os.path.join(output_cross_dir,'scaffolds.fasta')
        bioawk_output = os.path.join(output_cross_dir,'cross_filter.fasta')
        bioawk_op = 'bioawk -c fastx \'{ if(length($seq) > '+str(self.contig_length)+') { print \">\"$name; print $seq }}\' '+\
                        metaspade_output+' > '+bioawk_output
        os.system(bioawk_op)
        
        changed_name = os.path.join(output_cross_dir,'cross_contigs.fasta')
        self.rename(bioawk_output, changed_name, 'cross_')
        
        output_report_dir = os.path.join(output_cross_dir, 'cross_contig_report')
        self.create_output(output_report_dir)
        quast_op = 'quast.py '+changed_name+' -o '+output_report_dir
        os.system(quast_op)
        
        return [changed_name]
        
        
    def run_indiviudal(self, input_files):
        
        output_file_list = []
        
        for name in input_files:
            R1_input_dir = name[0]
            R2_input_dir = name[1]
            
            R1_file = R1_input_dir.split('/')[-1]
            subject_name = R1_file.split('R1')[0]
            
            output_sample_dir = os.path.join(self.output_dir, subject_name+'contigs')
            self.create_output(output_sample_dir)
            
            metaspade_op = 'spades.py  --meta -t 70 -m 400'+' -o '+output_sample_dir+' -1 '+R1_input_dir+' -2 '+R2_input_dir
            os.system(metaspade_op)
            
            metaspade_output = os.path.join(output_sample_dir,'scaffolds.fasta')
            bioawk_output = os.path.join(output_sample_dir, subject_name+'filter.fasta')
            bioawk_op = 'bioawk -c fastx \'{ if(length($seq) > '+str(self.contig_length)+') { print \">\"$name; print $seq }}\' '+\
                        metaspade_output+' > '+bioawk_output
            os.system(bioawk_op)
            
            changed_name = os.path.join(output_sample_dir,subject_name+'contigs.fasta')
            self.rename(bioawk_output, changed_name, subject_name)
            
            output_report_dir = os.path.join(output_sample_dir, subject_name+'contig_report')
            self.create_output(output_report_dir)
            quast_op = 'quast.py '+changed_name+' -o '+output_report_dir
            os.system(quast_op)
            
            output_file_list.append(changed_name)
            
        return output_file_list