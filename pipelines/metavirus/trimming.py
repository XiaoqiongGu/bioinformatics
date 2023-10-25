
import os
from metavirus.base import Moduler
import multiprocessing as mp


class adapter_module(Moduler):
    
    def __init__(self, output_dir, adapter_ref, adapter_quality):
        Moduler.__init__(self, output_dir)
        self.adapter_ref = adapter_ref
        self.adapter_quality = adapter_quality
        self.create_output()
    
    
    def run_module(self, input_files):

        log_dir = os.path.join(self.output_dir,'logs')
        self.create_output(log_dir)
        output_file_list = []
        
        pool = mp.Pool(processes=8)
        
        for name in input_files:
            R1_input_dir = name[0]
            R2_input_dir = name[1]
            
            R1_file = R1_input_dir.split('/')[-1]
            R2_file = R2_input_dir.split('/')[-1]
            
            R1_output_dir = os.path.join(self.output_dir, R1_file)
            R2_output_dir = os.path.join(self.output_dir, R2_file)
            
            log_file = os.path.join(log_dir, R1_file.split('R1')[0]+'log.txt')
            
            '''trim adapter and quality control'''
            adapter_op = 'bbduk.sh -Xmx50g t=40 in1='+R1_input_dir+' in2='+R2_input_dir+' ziplevel=9'+\
            ' out1='+R1_output_dir+' out2='+R2_output_dir+' ref='+self.adapter_ref+\
            ' minlen=50 qtrim=rl trimq='+str(self.adapter_quality)+' ktrim=r k=25 mink=11'+\
            ' hdist=1 stats='+log_file
            
            output_file_list.append([R1_output_dir,R2_output_dir])
            
            pool.apply_async(self.run_operation, (adapter_op,))
            
        pool.close()
        pool.join()
        
        return output_file_list
            
                  
    def run(self, input_files):     
        output_file_list = self.run_module(input_files)
        for file_pair in output_file_list:
            if not os.path.isfile(file_pair[0]):
                print('Missing adapter output: '+file_pair[0])
            if not os.path.isfile(file_pair[1]):
                print('Missing adapter output: '+file_pair[1])
                
        return output_file_list



class phix_module(Moduler):
    
    def __init__(self, output_dir, phix_ref, phix_quality):
        Moduler.__init__(self, output_dir)
        self.phix_ref = phix_ref
        self.phix_quality = phix_quality
        self.create_output()
        
        
    def run_module(self, input_files):
        
        log_dir = os.path.join(self.output_dir,'logs')
        self.create_output(log_dir)
        output_file_list = [];
        
        pool = mp.Pool(processes=8)
        
        for name in input_files:
            R1_input_dir = name[0]
            R2_input_dir = name[1]
            
            R1_file = R1_input_dir.split('/')[-1]
            R2_file = R2_input_dir.split('/')[-1]
            
            R1_output_dir = os.path.join(self.output_dir, R1_file)
            R2_output_dir = os.path.join(self.output_dir, R2_file)
            
            log_file = os.path.join(log_dir, R1_file.split('R1')[0]+'log.txt')
            
            '''phix read'''
            phix_op = 'bbduk.sh -Xmx50g t=40 in1='+R1_input_dir+' in2='+R2_input_dir+' ziplevel=9'+\
            ' out1='+R1_output_dir+' out2='+R2_output_dir+' ref='+self.phix_ref\
            +' trimq='+str(self.phix_quality)+' k=31 hdist=1 stats='+log_file
            
            output_file_list.append([R1_output_dir,R2_output_dir])
            
            pool.apply_async(self.run_operation, (phix_op,))
            
        pool.close()
        pool.join()
        
        return output_file_list
    
    
    def run(self, input_files):
        output_file_list = self.run_module(input_files)
        for file_pair in output_file_list:
            if not os.path.isfile(file_pair[0]):
                print('Missing phix output: '+file_pair[0])
            if not os.path.isfile(file_pair[1]):
                print('Missing phix output: '+file_pair[1])
                
        return output_file_list
            
            
             
class primer_module(Moduler):
    
    def __init__(self, output_dir, primer_ref, primer_quality):
        Moduler.__init__(self, output_dir)
        self.primer_ref = primer_ref
        self.primer_quality = primer_quality
        self.create_output()
        
    def run_module(self, input_files):
        
        log_dir = os.path.join(self.output_dir,'logs')
        self.create_output(log_dir)
        output_file_list = []
        
        pool = mp.Pool(processes=8)
        
        for name in input_files:
            
            R1_input_dir = name[0]
            R2_input_dir = name[1]
            
            R1_file = R1_input_dir.split('/')[-1]
            R2_file = R2_input_dir.split('/')[-1]
            
            R1_output_dir = os.path.join(self.output_dir, R1_file)
            R2_output_dir = os.path.join(self.output_dir, R2_file)
            
            log_file = os.path.join(log_dir, R1_file.split('R1')[0]+'log.txt')
            
            '''primer'''
            primer_op = 'bbduk.sh -Xmx50g t=40 in1='+R1_input_dir+' in2='+R2_input_dir+' ziplevel=9'+\
            ' out1='+R1_output_dir+' out2='+R2_output_dir+' literal='+self.primer_ref\
            +' minlen=50 trimq='+str(self.primer_quality)+' ktrim=l k=17 stats='+log_file
            #print(primer_op)
            
            output_file_list.append([R1_output_dir,R2_output_dir])
            
            pool.apply_async(self.run_operation, (primer_op,))
        
        pool.close()
        pool.join()
        
        return output_file_list
        
        
    def run(self, input_files):
        output_file_list = self.run_module(input_files)
        for file_pair in output_file_list:
            if not os.path.isfile(file_pair[0]):
                print('Missing primer output: '+file_pair[0])
            if not os.path.isfile(file_pair[1]):
                print('Missing primer output: '+file_pair[1])
                
        return output_file_list
    
    
class human_module(Moduler):

    def __init__(self, output_dir, human_ref):
        Moduler.__init__(self, output_dir)
        self.human_ref = human_ref
        self.create_output()


    def run_module(self, input_files):
        output_file_list = []
        
        
        for name in input_files:
            
            R1_input_dir = name[0]
            R2_input_dir = name[1]
            
            R1_file = R1_input_dir.split('/')[-1]
            R2_file = R2_input_dir.split('/')[-1]
            
            R1_output_dir = os.path.join(self.output_dir, R1_file)
            R2_output_dir = os.path.join(self.output_dir, R2_file)
            
            human_file = os.path.join(self.output_dir, R1_file.split('R1')[0]+'human.fastq.gz')
            
            '''human'''
            human_op = 'bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2'+\
            ' path='+self.human_ref+' qtrim=rl trimq=10 untrim -Xmx23g in1='+R1_input_dir+\
            ' in2='+R2_input_dir+' outu1='+R1_output_dir+' outu2='+R2_output_dir+' outm='+human_file
            #print(human_op)
            os.system(human_op)
                    
            output_file_list.append([R1_output_dir,R2_output_dir])
            
        
        return output_file_list


    def run(self, input_files):
        output_file_list = self.run_module(input_files)
        for file_pair in output_file_list:
            if not os.path.isfile(file_pair[0]):
                print('Missing primer output: '+file_pair[0])
            if not os.path.isfile(file_pair[1]):
                print('Missing primer output: '+file_pair[1])
                
        return output_file_list

        
        
        
class fastq_module(Moduler):
    
    def __init__(self, output_dir):
        Moduler.__init__(self, output_dir)
        self.create_output()
        
        
    def run(self, input_files):    
        
        for name in input_files:
            fastqc_op1 = 'fastqc -o '+self.output_dir+' '+name[0]
            os.system(fastqc_op1)
            fastqc_op2 = 'fastqc -o '+self.output_dir+' '+name[1]
            os.system(fastqc_op2)
            
            
            