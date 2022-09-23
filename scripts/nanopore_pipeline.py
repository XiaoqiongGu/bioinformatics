import pandas as pd
import os
import argparse
from Bio import SeqIO

def rename(original_file, rename_file, keyword):
    with open(original_file) as f1, open(rename_file, 'w') as f2:
        sequences = SeqIO.parse(f1, 'fasta')
        n = 1
        for record in sequences:
            length = len(record.seq)
            fname = keyword + str(n) + '_' + str(length)
            n += 1
            record.id = fname
            record.description = fname
            SeqIO.write(record, f2, 'fasta')

def blastn_id_merge(blast_output, merge_id):
    blast = pd.read_csv(blast_output, sep='\t', header=None)
    blast.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',\
         'sstart', 'send','evalue','bitscore','qlen','slen','qcovs']
    r1tem = blast[blast['qseqid'] == 'R1tem']['sseqid'].to_list()
    f3tem = blast[blast['qseqid'] == 'F3tem']['sseqid'].to_list()

    r1tem_f3tem = list(set(r1tem)&set(f3tem))

    with open(merge_id,'w') as f:
        for i in r1tem_f3tem:
            f.write(i + '\n')

def nano_filter(step1_input, step1_nanofilt_output, step1_nanoplot_output, quality, length, thread_num):
    print('******Start Step 1 NanoFilt and NanoPlot********')
    nano_filter_cmd = 'NanoFilt -q '+str(quality)+' -l '+str(length)+\
        ' --logfile nanofilt.log '+step1_input+' > '+step1_nanofilt_output
    os.system(nano_filter_cmd)
    nano_plot_cmd = 'NanoPlot --fastq '+step1_nanofilt_output+' -t '+str(thread_num) +\
        ' --plots {dot,kde} -o '+step1_nanoplot_output
    print('******Finish Step 1 NanoFilt and NanoPlot*******')

def preblast_format(step2_input, step2_original_output, step2_rename_output, keyword):
    print('*******Start Step 2 preBlast formatting*******')
    preblast_cmd = 'seqtk seq -a '+step2_input+' > '+step2_original_output
    os.system(preblast_cmd)
    rename(step2_original_output, step2_rename_output, keyword)
    print('*******Finish Step 2 preBlast formatting*******')

def blast(step3_input, blastdb, step3_blast_logfile, primer_file, step3_blast_output, thread_num):
    print('*******Start Step 3 Blast formatting*******')
    makeblastdb_cmd = 'makeblastdb -in '+step3_input+' -dbtype \'nucl\' -parse_seqids -out '+blastdb+\
        ' -logfile '+step3_blast_logfile
    os.system(makeblastdb_cmd)
    runblast_cmd = 'blastn -db '+blastdb+' -query '+primer_file+\
        ' -task blastn-short -outfmt \"6 std qlen slen qcovs\" -max_target_seqs 10000000'+\
        ' -perc_identity 80 -num_threads '+str(thread_num)+' -out '+step3_blast_output
    os.system(runblast_cmd)
    print('*******Finish Step 3 Blast formatting*******')


def pullingout_ARG_seqs(step3_blast_output, ARG_id_file, step2_rename_output, step4_ARG_seqs):
    print('*******Start Step 4 pulling out ARG seqs*******')
    blastn_id_merge(step3_blast_output, ARG_id_file)
    pullingout_cmd = 'seqtk subseq '+step2_rename_output+' '+ARG_id_file+' > '+step4_ARG_seqs
    os.system(pullingout_cmd)
    print('*******Finish Step 4 pulling out ARG seqs*******')

def emu(step4_ARG_seqs, step5_emu_folder, thread_num):
    print('*******Start Step 5 emu analysis*******')
    emu_cmd = 'emu abundance '+step4_ARG_seqs+' --output-dir '+step5_emu_folder+' --threads '+str(thread_num)+\
        ' --keep-counts'
    os.system(emu_cmd)
    print('*******Finish Step 5 emu analysis*******')

def pipeline():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', type=str, required=True)
    parser.add_argument('--input', type=str, required=True, help='nanopore base calling fastq')
    parser.add_argument('--quality', default=7, required=False, help='quality control value')
    parser.add_argument('--length', default=300, required=False, help='length control value')
    parser.add_argument('--keyword', type=str, required=True, help='rename fasta header keyword, for example: TEM')
    parser.add_argument('--primer_file', type=str, required=True, help='epicpcr primers for blast')
    parser.add_argument('--thread_num', type=int, required=True)
    flags, _ = parser.parse_known_args()

    #input_name = flags.input.spilt('/')[-1].split('.fastq')[0]
    input_name = flags.input.split('.fastq')[0]
    step1_nanofilt_dir = os.path.join(flags.dir, 'nanofilt')
    if not os.path.isdir(step1_nanofilt_dir):
        os.makedirs(step1_nanofilt_dir)
    step1_nanofilt_output = os.path.join(step1_nanofilt_dir, input_name+'.nanofilt.fastq')
    
    step1_nanoplot_dir = os.path.join(step1_nanofilt_dir, 'nanoplot')
    if not os.path.isdir(step1_nanoplot_dir):
        os.makedirs(step1_nanoplot_dir)
    
    nano_filter(flags.input, step1_nanofilt_output,
                step1_nanoplot_dir, flags.quality, flags.length, flags.thread_num)

    step2_original_output = os.path.join(step1_nanofilt_dir, input_name+'.nanofilt.fasta')
    step2_rename_output = os.path.join(step1_nanofilt_dir, input_name+'.nanofilt.rename.fasta')
    preblast_format(step1_nanofilt_output, step2_original_output,
                    step2_rename_output, flags.keyword)
    
    blastdb_dir = os.path.join(flags.dir, 'blast', 'blastdb')
    if not os.path.isdir(blastdb_dir):
        os.makedirs(blastdb_dir)
    step3_blast_logfile = os.path.join(blastdb_dir, 'blastdb.log')
    step3_blast_output = os.path.join(blastdb_dir, input_name+'.blast.txt')
    blast(step2_rename_output, blastdb_dir, step3_blast_logfile,
          flags.primer_file, step3_blast_output, flags.thread_num)

    blast_dir = os.path.join(flags.dir, 'blast')
    ARG_id_file = os.path.join(blast_dir, input_name+'.arg.id.txt')
    step4_ARG_seqs = os.path.join(blast_dir, input_name+'.arg.seqs.fasta')
    pullingout_ARG_seqs(step3_blast_output, ARG_id_file,
                        step2_rename_output, step4_ARG_seqs)
    
    step5_emu_folder = os.path.join(flags.dir, 'emu')
    if not os.path.isdir(step5_emu_folder):
        os.makedirs(step5_emu_folder)
    emu(step4_ARG_seqs, step5_emu_folder, flags.thread_num)

if __name__ == '__main__':
    pipeline()
