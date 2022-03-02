
import os

from Bio import SeqIO


class Moduler(object):

    def __init__(self, output_dir):
        self.output_dir = output_dir


    def read_input(self, input_dir):
        name_list = []
        R1_list = []
        R1_count = 0
        R2_count = 0

        directory_list = os.listdir(input_dir)

        for filename in directory_list:
            if filename.find('R1') != -1:
                R1_list.append(filename)
                R1_count = R1_count+1
            elif filename.find('R2') != -1:
                R2_count = R2_count+1

        if len(R1_list) == 0:
            raise ValueError('The input folder is empty or the file names are not follow the correct format.')

        if R1_count != R2_count:
            raise ValueError('File missing: R1 R2 files not match.')

        for R1_file in R1_list:
            file_split = R1_file.split('R1')
            R2_file = file_split[0]+'R2'+file_split[1]
            new_dir = os.path.join(input_dir, R2_file)
            if not os.path.isfile(new_dir):
                raise ValueError('Missing file: '+R2_file)
            name_list.append([R1_file,R2_file])

        return name_list


    def create_output(self, other_output_dir=None):
        if other_output_dir == None:
            if not os.path.isdir(self.output_dir):
                os.makedirs(self.output_dir)
        else:
            if not os.path.isdir(other_output_dir):
                os.makedirs(other_output_dir)


    def rename(self, input_file, output_file, prefix):
        input_file = SeqIO.parse(input_file,'fasta')
        new_seq = []

        count = 0
        for records in input_file:
            records.name = prefix+str(count)
            records.id = prefix+str(count)
            records.description = records.id+'_len='+str(len(records.seq))
            count += 1
            print(records)
            new_seq.append(records)
        SeqIO.write(new_seq, output_file, 'fasta')
        
        
    def run_operation(self, operation):
        os.system(operation)



class readin_module(Moduler):

    def __init__(self, output_dir):
        Moduler.__init__(self, output_dir)


    def run(self, input_dir):
        file_list = self.read_input(input_dir)
        output_file_list = []
        for name in file_list:
            R1_dir = os.path.join(input_dir, name[0])
            R2_dir = os.path.join(input_dir, name[1])
            output_file_list.append([R1_dir, R2_dir])
        return output_file_list



