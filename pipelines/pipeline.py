
import os
import sys
import numpy as np
from Bio import SeqIO


def rename(input_file, output_file, prefix):
    input_file = SeqIO.parse(input_file,'fasta')
    new_seq = []
    
    count = 0
    for records in input_file:
        records.name = prefix+str(count)
        records.id = prefix+str(count)
        records.description = records.id+' len='+str(len(records.seq))
        count += 1
        print(records)
        new_seq.append(records)
        
    SeqIO.write(new_seq, output_file, 'fasta')


def write_text(input_dir,sample_name,collapse_rank,select_rank,tax_format,egg_format,seed_format,kegg_format):
    exportcode_dir = input_dir+'code/'
    if not os.path.isdir(exportcode_dir):
        os.mkdir(exportcode_dir)
            
    exporttax_dir = input_dir+'Taxonomy/'
    if not os.path.isdir(exporttax_dir):
        os.mkdir(exporttax_dir)
        
    exportegg_dir = input_dir+'Eggnog/'
    if not os.path.isdir(exportegg_dir):
        os.mkdir(exportegg_dir)
        
    exportseed_dir = input_dir+'Seed/'
    if not os.path.isdir(exportseed_dir):
        os.mkdir(exportseed_dir)
        
    exportkegg_dir = input_dir+'Kegg/'
    if not os.path.isdir(exportkegg_dir):
        os.mkdir(exportkegg_dir)
      
    sample_split = sample_name.split('prodigal')
    
    code_file = exportcode_dir+sample_split[0]+'code.txt'
    tax_file = exporttax_dir+sample_split[0]+'taxon.txt'
    egg_file = exportegg_dir+sample_split[0]+'egg.txt'
    seed_file = exportseed_dir+sample_split[0]+'seed.txt'
    kegg_file = exportkegg_dir+sample_split[0]+'kegg.txt'
    
    codeF = open(code_file,'w')
    codeF.write("open file = '"+input_dir+sample_name+"';\n")
    codeF.write('show window=mainViewer;set context=Taxonomy;update;\n')
    codeF.write("collapse rank='"+collapse_rank+"';select rank='"+select_rank+"';\n")
    codeF.write('select nodes = all;\n')
    codeF.write('export what=CSV format='+tax_format+" separator=comma counts=assigned file='"+tax_file+"';\n")
    
    codeF.write('open viewer=EGGNOG;\n')
    codeF.write('set context=EGGNOG;\n')
    codeF.write('update;\n')
    codeF.write('uncollapse nodes=all;\n')
    codeF.write('select nodes=all;\n')
    codeF.write('export what=CSV format='+egg_format+" separator=comma counts=assigned file='"+egg_file+"';\n")
    codeF.write('close what=current;\n')

    codeF.write('open viewer=SEED;\n')
    codeF.write('set context=SEED;\n')
    codeF.write('update;\n')
    codeF.write('uncollapse nodes=all;\n')
    codeF.write('select nodes=all;\n')
    codeF.write('export what=CSV format='+seed_format+" separator=comma counts=assigned file='"+seed_file+"';\n")
    codeF.write('close what=current;\n')
    
    codeF.write('open viewer=KEGG;\n')
    codeF.write('set context=KEGG;\n')
    codeF.write('update;\n')
    codeF.write('uncollapse nodes=all;\n')
    codeF.write('select nodes=all;\n')
    codeF.write('export what=CSV format='+kegg_format+" separator=comma counts=assigned file='"+kegg_file+"';\n")
    codeF.write('close what=current;\n')
    codeF.write('quit;\n')
    codeF.close()
    
    return code_file



def step0_name_extraction(input_dir):
    input_dir = input_dir+'/input'
    for filename in os.listdir(input_dir):
        if filename[-5:len(filename)] == 'fastq':
            operation1 = 'gzip '+filename
            os.system(operation1)
            operation2 = 'rm '+ filename
            os.system(operation2)
        elif filename[-8:len(filename)] != 'fastq.gz':
            raise ValueError('The folder contains unwanted files.')
    
    name_list = []   
    R1_count = 0
    R2_count = 0     
    for filename in os.listdir(input_dir):
        if filename.find('R1') != -1:
            name_list.append(filename)
            R1_count = R1_count+1
        elif filename.find('R2') != -1:
            R2_count = R2_count+1
    
    if R1_count != R2_count:
        raise ValueError('File missing: R1 R2 files not match.')
        
    for filename in name_list:
        file_split = filename.split('R1')
        file_new = file_split[0]+'R2'+file_split[1]
        new_dir = input_dir+'/'+file_new
        if not os.path.isfile(new_dir):
            raise ValueError('Missing file:'+file_new)
        
    return name_list


def step1_trimming_ada(input_dir, name_list, adapter_name, quality_con):
    outtrim_dir = input_dir+'/adapter/'
    if not os.path.isdir(outtrim_dir):
        os.mkdir(input_dir+'/adapter')
    
    for i in range(len(name_list)):
        R1_file = name_list[i]
        R2_split = R1_file.split('R1')
        R2_file = R2_split[0]+'R2'+R2_split[1]
        
        R1_fileM = input_dir+'/input/'+R1_file
        R2_fileM = input_dir+'/input/'+R2_file
        
        '''trim adapter and quality control'''
        bbduk_op = 'bbduk.sh -Xmx20g in1='+R1_fileM+' in2='+R2_fileM+' ziplevel=9'+\
        ' out1='+outtrim_dir+R1_file+' out2='+outtrim_dir+R2_file+' ref='+adapter_name\
        +' minlen=50 qtrim=rl trimq='+str(quality_con)+' ktrim=r k=25 mink=11'+\
        ' hdist=1 stats='+outtrim_dir+R2_split[0]+'log.txt'
        os.system(bbduk_op)
    


def step1_trimming_phix(input_dir, name_list, phix_name, quality_con):
    outphix_dir = input_dir+'/phix/'
    if not os.path.isdir(outphix_dir):
        os.mkdir(outphix_dir)
        
    for i in range(len(name_list)):
        R1_file = name_list[i]
        R2_split = R1_file.split('R1')
        R2_file = R2_split[0]+'R2'+R2_split[1]
        
        R1_fileM = input_dir+'/adapter/'+R1_file
        R2_fileM = input_dir+'/adapter/'+R2_file
        
        '''phix read'''
        phix_op = 'bbduk.sh -Xmx20g in1='+R1_fileM+' in2='+R2_fileM+' ziplevel=9'+\
        ' out1='+outphix_dir+R1_file+' out2='+outphix_dir+R2_file+' ref='+phix_name\
        +' minlen=50 qtrim=rl trimq='+str(quality_con)+' ktrim=r k=31 mink=11'+\
        ' hdist=1 stats='+outphix_dir+R2_split[0]+'.txt'
        os.system(phix_op)
        
           
       

def step1_trimming_primer(input_dir, name_list, primer_name, quality_con):
    outprimer_dir = input_dir+'/primer/'
    if not os.path.isdir(outprimer_dir):
        os.mkdir(outprimer_dir)
        
        
    for i in range(len(name_list)):
        R1_file = name_list[i]
        R2_split = R1_file.split('R1')
        R2_file = R2_split[0]+'R2'+R2_split[1]
        
        R1_fileM = input_dir+'/phix/'+R1_file
        R2_fileM = input_dir+'/phix/'+R2_file
        
        
        '''primer'''
        primer_op = 'bbduk.sh -Xmx20g in1='+R1_fileM+' in2='+R2_fileM+' ziplevel=9'+\
        ' out1='+outprimer_dir+R1_file+' out2='+outprimer_dir+R2_file+' ref='+primer_name\
        +' minlen=50 qtrim=rl trimq='+str(quality_con)+' ktrim=r k=17 mink=11'+\
        ' hdist=1 stats='+outprimer_dir+R2_split[0]+'.txt'
        os.system(primer_op)
   
     
def step1_trimming_fastq(input_dir, mode_sel):
    if mode_sel == 1:
        file_dir = input_dir+'/adapter/'
    elif mode_sel == 2:
        file_dir = input_dir+'/phix/'
    elif mode_sel == 3:
        file_dir = input_dir+'/primer/'
      
    outfastqc_dir = input_dir+'/fastqc/'
    if not os.path.isdir(outfastqc_dir):
        os.mkdir(outfastqc_dir)
    fastqc_op = 'fastqc -o '+outfastqc_dir+' '+file_dir+'*.fastq.gz'
    os.system(fastqc_op)
        


      
def step1_trimming_warper(input_dir, name_list, quality_con, mode_sel, *database_name):
        
    database_env = os.environ['BBMAP_RESOURCE']+'/'
    if mode_sel == 1:
        adapter_env = database_env+database_name[0]
        if not os.path.isfile(adapter_env):
            raise ValueError('Adapter database file does not exsit.')
        step1_trimming_ada(input_dir, name_list, adapter_env, quality_con)
    
    elif mode_sel == 2:
        adapter_env = database_env+database_name[0]
        if not os.path.isfile(adapter_env):
            raise ValueError('Adapter database file does not exsit.')
        phix_env = database_env+database_name[1]
        if not os.path.isfile(phix_env):
            raise ValueError('Phix database file does not exsit')
        step1_trimming_ada(input_dir, name_list, adapter_env, quality_con)
        step1_trimming_phix(input_dir, name_list, phix_env, quality_con)
        
    elif mode_sel == 3:
        adapter_env = database_env+database_name[0]
        if not os.path.isfile(adapter_env):
            raise ValueError('Adapter database file does not exsit.')
        phix_env = database_env+database_name[1]
        if not os.path.isfile(phix_env):
            raise ValueError('Phix database file does not exsit')
        primer_env = database_env+database_name[2]
        if not os.path.isfile(primer_env):
            raise ValueError('Primer databse file does not exsit')
        step1_trimming_ada(input_dir, name_list, adapter_env, quality_con)
        step1_trimming_phix(input_dir, name_list, phix_env, quality_con)
        step1_trimming_primer(input_dir, name_list, primer_env, quality_con)
        

def step2_assembly_megahit(input_dir, name_list, mode_sel, megahit_ass_mode, cong_length):
    if mode_sel == 1:
        get_file = '/adapter/'
    elif mode_sel == 2:
        get_file = '/phix/'
    elif mode_sel == 3:
        get_file = '/primer/'
    else:
        raise ValueError('Invalide mode selection indicator.')
    
    outmegahit_dir = input_dir+'/megahit/'
    if not os.path.isdir(outmegahit_dir):
        os.mkdir(outmegahit_dir)
    
    quality_dir = outmegahit_dir+'contig_quality/'
    if not os.path.isdir(quality_dir):
        os.mkdir(quality_dir)
    
    if megahit_ass_mode == 'cross':
        outmegahit_dir = outmegahit_dir+'cross'
        if os.path.isdir(outmegahit_dir):
            os.rmdir(outmegahit_dir)
        
        R1_file = name_list[0]
        R2_split = R1_file.split('R1')
        R2_file = R2_split[0]+'R2'+R2_split[1]
        R1_fileM = input_dir+get_file+R1_file
        R2_fileM = input_dir+get_file+R2_file
        for i in range(1, len(name_list)):
            R1_file = name_list[i]
            R2_split = R1_file.split('R1')
            R2_file = R2_split[0]+'R2'+R2_split[1]
            
            R1_fileM = R1_fileM+','+input_dir+get_file+R1_file
            R2_fileM = R2_fileM+','+input_dir+get_file+R2_file
            
        megahit_op = 'megahit -o '+outmegahit_dir+' --min-contig-len '+str(cong_length)+\
        ' -1 '+R1_fileM+' -2 '+R2_fileM
        os.system(megahit_op)
        
        '''name change'''
        original_name = outmegahit_dir+'/final.contigs.fa'
        changed_name = outmegahit_dir+'/cross.contigs.fa'
        rename(original_name, changed_name, 'cross_')
        
        '''get report'''
        report_dir = outmegahit_dir+'/cross_contig_report/'
        if not os.path.isdir(report_dir):
            os.mkdir(report_dir)
            
        quast_op = 'quast '+changed_name+' -o '+report_dir
        os.system(quast_op)
        
        report_copy_op = 'cp -R '+report_dir+' '+quality_dir
        os.system(report_copy_op)
            
            
    elif megahit_ass_mode == 'indi':
        
        for i in range(len(name_list)):
            R1_file = name_list[i]
            R2_split = R1_file.split('R1')
            R2_file = R2_split[0]+'R2'+R2_split[1]
            
            R1_fileM = input_dir+get_file+R1_file
            R2_fileM = input_dir+get_file+R2_file
            
            outmegahit_file = outmegahit_dir+R2_split[0]+'contigs'
            if os.path.isdir(outmegahit_file):
                os.rmdir(outmegahit_file)
            
            '''run megahit'''
            megahit_op = 'megahit -o '+outmegahit_file+' --min-contig-len '+str(cong_length)+\
            ' -1 '+R1_fileM+' -2 '+R2_fileM
            os.system(megahit_op)
            
            '''name change'''
            original_name = outmegahit_file+'/final.contigs.fa'
            changed_name = outmegahit_file+'/'+R2_split[0]+'contigs.fa'
            rename(original_name, changed_name, R2_split[0])
            
            '''get report'''
            report_dir = outmegahit_file+'/'+R2_split[0]+'contig_report/'
            if not os.path.isdir(report_dir):
                os.mkdir(report_dir)
                
            quast_op = 'quast '+changed_name+' -o '+report_dir
            os.system(quast_op)
            
            report_copy_op = 'cp -R '+report_dir+' '+quality_dir
            os.system(report_copy_op)
            
    else:
        raise ValueError('Wrong MegaHit assembly mode indicator.')
    

def step2_assembly_metaspade(input_dir, name_list, mode_sel, metaspade_ass_mode, cong_length):
    if mode_sel == 1:
        get_file = '/adapter/'
    elif mode_sel == 2:
        get_file = '/phix/'
    elif mode_sel == 3:
        get_file = '/primer/'
    else:
        raise ValueError('Invalide mode selection indicator.')
    
    outmetaspade_dir = input_dir+'/metaspade/'
    if not os.path.isdir(outmetaspade_dir):
        os.mkdir(outmetaspade_dir)
    
    quality_dir = outmetaspade_dir+'contig_quality/'
    if not os.path.isdir(quality_dir):
        os.mkdir(quality_dir)
    
    
    if metaspade_ass_mode == 'cross':
        print('I am empty')
        
    elif metaspade_ass_mode == 'indi':
        
        for i in range(len(name_list)):
            R1_file = name_list[i]
            R2_split = R1_file.split('R1')
            R2_file = R2_split[0]+'R2'+R2_split[1]
            
            R1_fileM = input_dir+get_file+R1_file
            R2_fileM = input_dir+get_file+R2_file
            
            outmetaspade_dir = outmetaspade_dir+R2_split[0]+'contigs'
            if os.path.isdir(outmetaspade_dir):
                os.rmdir(outmetaspade_dir)
            
            '''run metaspade'''
            metaspade_op = 'spades.py  -- meta'+' -o '+outmetaspade_dir+' -1 '+\
            R1_fileM+' -2 '+R2_fileM
            os.system(metaspade_op)
            
            '''name change'''
            original_name = outmetaspade_dir+'/scaffolds.fasta'
            changed_name = outmetaspade_dir+'/'+R2_split[0]+'scaffolds.fasta'
            os.rename(original_name,changed_name)
            
            
            outmetabioawk_dir = outmetaspade_dir+'/filter/'
            if not os.path.isdir(outmetabioawk_dir):
                os.mkdir(outmetabioawk_dir)
            
            '''run bioawk'''
            filter_file = outmetabioawk_dir+R2_split[0]+'scaffolds.fasta'
            bioawk_op = 'bioawk -c fastx \'{ if(length($seq) > '+str(cong_length)\
                        +') { print \">\"$name; print $seq }}\' '+\
                        changed_name+' > '+filter_file
            os.system(bioawk_op)
            
            report_dir = outmetaspade_dir+'/'+R2_split[0]+'contig_report/'
            if not os.path.isdir(report_dir):
                os.mkdir(report_dir)
                
            quast_op = 'quast '+filter_file+' -o '+report_dir
            os.system(quast_op)
            
            report_copy_op = 'cp -R '+report_dir+' '+quality_dir
            os.system(report_copy_op)
            
    else:
        raise ValueError('Wrong MetaSpade assembly mode indicator.')


def step2_assembly_warper(input_dir, name_list, mode_sel, assembly_sel, assembly_mode, cong_length):
    if assembly_sel == 'megahit':
        if assembly_mode == 'indi':
            step2_assembly_megahit(input_dir, name_list, mode_sel, assembly_mode, cong_length)
        elif assembly_mode == 'cross':
            step2_assembly_megahit(input_dir, name_list, mode_sel, assembly_mode, cong_length)
        elif assembly_mode == 'hybrid':
            step2_assembly_megahit(input_dir, name_list, mode_sel, 'indi', cong_length)
            step2_assembly_megahit(input_dir, name_list, mode_sel, 'cross', cong_length)
        else:
            raise ValueError('Invalid assembly mode.')
            
    elif assembly_sel == 'metaspade':
        if assembly_mode == 'indi':
            step2_assembly_metaspade(input_dir, name_list, mode_sel, assembly_mode, cong_length)
        else:
            raise ValueError('Invalid assembly mode.')
            
    elif assembly_sel == 'mega_meta':
        if assembly_mode == 'indi':
            step2_assembly_megahit(input_dir, name_list, mode_sel, assembly_mode, cong_length)
            step2_assembly_metaspade(input_dir, name_list, mode_sel, assembly_mode, cong_length)
        else:
            raise ValueError('Invalid assembly mode.')




def step3_alignment_diamond_blastp(input_dir, assembly_sel, diamond_database, e_value, blast2rma_ref, \
                                   acc2tax2018,acc2eggnog2016, acc2interpro2018, acc2kegg2017, acc2seed2015,\
                                   collapse_rank, select_rank, tax_format, egg_format, seed_format, kegg_format):
    '''extract open reading frame'''
    prodigal_dir = input_dir+'/prodigal/'
    if not os.path.isdir(prodigal_dir):
        os.mkdir(prodigal_dir)
        
    diamond_dir = input_dir+'/diamond/'
    if not os.path.isdir(diamond_dir):
        os.mkdir(diamond_dir)
        
    blast2rma_dir = input_dir+'/blast2brma/'
    if not os.path.isdir(blast2rma_dir):
        os.mkdir(blast2rma_dir)
    
    if assembly_sel == 'megahit':
        prod_sub_dir = prodigal_dir+'megahit/'
        if not os.path.isdir(prod_sub_dir):
            os.mkdir(prod_sub_dir)
            
        diam_sub_dir = diamond_dir+'megahit/'
        if not os.path.isdir(diam_sub_dir):
            os.mkdir(diam_sub_dir)
            
        blast_sub_dir = blast2rma_dir+'megahit/'
        if not os.path.isdir(blast_sub_dir):
            os.mkdir(blast_sub_dir)
            
        assembly_dir = input_dir+'/'+assembly_sel+'/'
        assemnly_folder = os.listdir(assembly_dir)
        
        for folder_name in assemnly_folder:
            if folder_name != 'contig_quality':
                
                '''prodigal part'''
                megahit_file = assembly_dir+folder_name+'/'+folder_name+'.fa'
                fna_file = prod_sub_dir+folder_name+'_prodigal.fna'
                faa_file = prod_sub_dir+folder_name+'_prodigal.faa'
                pro_mega_op = 'prodigal -p meta -i '+megahit_file\
                +' -d '+fna_file+' -a '+faa_file
                os.system(pro_mega_op)
                print('finish prodigal')
                
                '''diamond part'''
                diamond_raw = diam_sub_dir+folder_name+'_prodigal.m8'
                diamond_op = 'diamond blastp -d '+diamond_database+' -q '+faa_file\
                +' -o '+diamond_raw+' -e '+e_value
                os.system(diamond_op)
                print('finish diamond')
                
                '''cut part'''
                diamond_cut = diam_sub_dir+folder_name+'_prodigal.contigs.m8'
                cut_op = 'cut -f1 '+diamond_raw+'|rev|cut -f2- -d"_"|rev|paste - '+\
                diamond_raw+'|cut -f1,3- > '+diamond_cut
                os.system(cut_op)
                print('finish cut')
                
                '''blast2rma'''
                blast_out_dir = blast_sub_dir+folder_name+'_prodigal.contigs.rma'
                blast2rma_op = blast2rma_ref+' -i '+diamond_cut+' -f BlastTab -bm BlastP -o '+\
                blast_out_dir+' -me '+e_value+' -a2t '+acc2tax2018+' -a2eggnog '+acc2eggnog2016+\
                ' -a2interpro2go '+acc2interpro2018+' -a2kegg '+acc2kegg2017+' -a2seed '+acc2seed2015
                os.system(blast2rma_op)
                print('finish blast2rma')
                
                
                '''export'''
                sample_name = folder_name+'_prodigal.contigs.rma'
                code_file = write_text(blast_sub_dir,sample_name,collapse_rank,select_rank,\
                                       tax_format,egg_format,seed_format,kegg_format)
                export_op = 'xvfb-run --auto-servernum --server-num=1 ~/megan/MEGAN -g -c '+code_file
                os.system(export_op)
                
            
    elif assembly_sel == 'metaspade':
        prod_sub_dir = prodigal_dir+'metaspade/'
        if not os.path.isdir(prod_sub_dir):
            os.mkdir(prod_sub_dir)
               
    else:
        raise ValueError('Invalid assembly mode.')
        
        
        
def step4_quantification_bowtie2_mapping(input_dir, mode_sel):
    
    if mode_sel == 1:
        get_file = '/adapter/'
    elif mode_sel == 2:
        get_file = '/phix/'
    elif mode_sel == 3:
        get_file = '/primer/'
    else:
        raise ValueError('Invalide mode selection indicator.')
        
    cross_ass_dir = input_dir+get_file
    database_dir = input_dir+'/megahit/cross/cross.contigs.fa'
    output_dir = input_dir+'/bowtie2/'
    
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    input_files = os.listdir(cross_ass_dir)
    
    R1_count = 0
    R2_count = 0
    for input_file in input_files:
        if '1.' in input_file:
            R1_count += 1
        if '2.' in input_file:
            R2_count += 1
    
    if R1_count != R2_count:
        raise ValueError('R1,R1 files are not matched.')
    
    #grep patterns catering *1.fastq.gz  
    name_list = []
    for input_file in input_files:
        if '1.' in input_file:
            R1_name = input_file.split('1.')
            R2_name = R1_name[0]+'2.'+R1_name[1]
            if not os.path.isfile(cross_ass_dir+'/'+R2_name):
                raise ValueError('file '+R2_name+'is missing')
            name_list.append(R1_name)
    
    database_size = os.path.getsize(database_dir)
    database_size = database_size/1024/1024/1024
    
    database_file = database_dir.split('/')[-1]
    database_name = database_file.split('.')[0]
    
    if database_size > 3:
        bt1_op = 'bowtie2-build -f '+database_dir+' '+output_dir+' -p 40 --large-index'
    else:
        bt1_op = 'bowtie2-build -f '+database_dir+' '+output_dir+'/'+database_name+' -p 40'
    
    os.system(bt1_op)
    
    for sample_name in name_list:
        R1_name = os.path.join(cross_ass_dir, sample_name[0]+'1.'+sample_name[1])
        R2_name = os.path.join(cross_ass_dir, sample_name[0]+'2.'+sample_name[1])
        sam_file = os.path.join(output_dir,sample_name[0]+'OUT.sam')
        database_path = os.path.join(output_dir, database_name)
        bt2_op = 'bowtie2 -x '+database_path+' -1 '+R1_name+' -2 '+R2_name+' -S '+sam_file+' -p 40'
        os.system(bt2_op)
    
        bam_file = os.path.join(output_dir, sample_name[0]+'OUT.bam')
        os.system('samtools view -bS '+sam_file+' > '+bam_file)
    
        sort_file = os.path.join(output_dir, sample_name[0]+'OUT.sorted.bam')
        os.system('samtools sort '+bam_file+' -o '+sort_file)
    
        os.system('samtools index '+sort_file)
    
        os.system('samtools depth '+sort_file)
    
        idx_file = os.path.join(output_dir, sample_name[0]+'OUT.idxstats.txt')
        os.system('samtools idxstats '+sort_file+' > '+idx_file)
    
        genome_file = os.path.join(output_dir, sample_name[0]+'OUT.genome')
        os.system('bedtools genomecov -ibam '+sort_file+' -g '+genome_file)

    absolute_count_op = 'get_count_table.py '+output_dir+'*.idxstats.txt > '+output_dir+'absolute_count.txt'
    os.system(absolute_count_op)
    
    abs2RPKM_op = 'absolutereads2RPKM.py '+output_dir+'absolute_count.txt '+output_dir+'RPKM.txt'
    os.system(abs2RPKM_op)


input_dir = sys.argv[1]
mode_sel = int(sys.argv[2])
adapter_ref = sys.argv[3]
phix_ref = sys.argv[4]
primer_ref = sys.argv[5]
assembly_sel = sys.argv[6]
assembly_mode = sys.argv[7]
cong_length = int(sys.argv[8])


diamond_database = sys.argv[9]
e_value = sys.argv[10]
blast2rma_ref = sys.argv[11]
acc2tax2018 = sys.argv[12]
acc2eggnog2016 = sys.argv[13]
acc2interpro2018 = sys.argv[14]
acc2kegg2017 = sys.argv[15]
acc2seed2015 = sys.argv[16]
collapse_rank = sys.argv[17]
select_rank = sys.argv[18]
tax_format = sys.argv[19]
egg_format = sys.argv[20]
seed_format = sys.argv[21]
kegg_format = sys.argv[22]



'''Trimming'''
print('**********************************************')
print('start trimming')
name_list = step0_name_extraction(input_dir)
step1_trimming_warper(input_dir, name_list, 30, mode_sel, adapter_ref, phix_ref, primer_ref)
step1_trimming_fastq(input_dir, mode_sel)
print('trimming complete')

'''Assembly'''
print('**********************************************')
print('start assembly')
step2_assembly_warper(input_dir, name_list, mode_sel, assembly_sel, assembly_mode, cong_length)
print('assembly complete')


'''Aligment'''
print('**********************************************')
print('start alignment')
step3_alignment_diamond_blastp(input_dir, assembly_sel, diamond_database, e_value, blast2rma_ref, acc2tax2018,\
                                   acc2eggnog2016, acc2interpro2018, acc2kegg2017, acc2seed2015, collapse_rank,\
                                   select_rank, tax_format, egg_format, seed_format, kegg_format)

'''Quantification'''
print('**********************************************')
print('start quantification')
step4_quantification_bowtie2_mapping(input_dir, mode_sel)




