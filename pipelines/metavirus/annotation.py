
import os
import time
from metavirus.base import Moduler


class diamond_module(Moduler):
    
    def __init__(self, output_dir, diamond_database, e_value, blast2rma_ref,
                 acc2tax, acc2eggnog, acc2interpro, acc2kegg, acc2seed, collapse_rank, 
                 select_rank, tax_format, egg_format, seed_format, kegg_format):
        Moduler.__init__(self, output_dir)
        self.diamond_database = diamond_database
        self.e_value = e_value
        self.blast2rma_ref = blast2rma_ref
        self.acc2tax = acc2tax
        self.acc2eggnog = acc2eggnog
        self.acc2interpro = acc2interpro
        self.acc2kegg = acc2kegg
        self.acc2seed = acc2seed
        self.collapse_rank = collapse_rank
        self.select_rank = select_rank
        self.tax_format = tax_format
        self.egg_format = egg_format
        self.seed_format = seed_format
        self.kegg_format = kegg_format
        self.create_output()
        
        self.prodigal_dir = os.path.join(self.output_dir, 'prodigal')
        self.create_output(self.prodigal_dir)
        
        self.diamond_dir = os.path.join(self.output_dir, 'diamond')
        self.create_output(self.diamond_dir)
        
        self.blast2rma_dir = os.path.join(self.output_dir, 'blast2brma')
        self.create_output(self.blast2rma_dir)
        

    def write_text(self, sample_name):
        exportcode_dir = os.path.join(self.blast2rma_dir, 'Code')
        self.create_output(exportcode_dir)

        exporttax_dir = os.path.join(self.blast2rma_dir, 'Taxonomy')
        self.create_output(exporttax_dir)

        exportegg_dir = os.path.join(self.blast2rma_dir, 'Eggnog')  
        self.create_output(exportegg_dir)

        exportseed_dir = os.path.join(self.blast2rma_dir, 'Seed')
        self.create_output(exportseed_dir)

        exportkegg_dir = os.path.join(self.blast2rma_dir, 'Kegg')
        self.create_output(exportkegg_dir)          
          
        subject_name = sample_name.split('prodigal')[0]
        
        code_file = os.path.join(exportcode_dir, subject_name+'code.txt')
        tax_file = os.path.join(exporttax_dir, subject_name+'taxon.txt')
        egg_file = os.path.join(exportegg_dir, subject_name+'egg.txt')
        seed_file = os.path.join(exportseed_dir, subject_name+'seed.txt')
        kegg_file = os.path.join(exportkegg_dir, subject_name+'kegg.txt')
        
        codeF = open(code_file,'w')
        
        codeF.write("open file = '"+self.blast2rma_dir+subject_name+"';\n")
        codeF.write('show window=mainViewer;set context=Taxonomy;update;\n')
        codeF.write("collapse rank='"+self.collapse_rank+"';select rank='"+self.select_rank+"';\n")
        codeF.write('select nodes = all;\n')
        codeF.write('export what=CSV format='+self.tax_format+" separator=comma counts=assigned file='"+tax_file+"';\n")
        
        codeF.write('open viewer=EGGNOG;\n')
        codeF.write('set context=EGGNOG;\n')
        codeF.write('update;\n')
        codeF.write('uncollapse nodes=all;\n')
        codeF.write('select nodes=all;\n')
        codeF.write('export what=CSV format='+self.egg_format+" separator=comma counts=assigned file='"+egg_file+"';\n")
        codeF.write('close what=current;\n')
    
        codeF.write('open viewer=SEED;\n')
        codeF.write('set context=SEED;\n')
        codeF.write('update;\n')
        codeF.write('uncollapse nodes=all;\n')
        codeF.write('select nodes=all;\n')
        codeF.write('export what=CSV format='+self.seed_format+" separator=comma counts=assigned file='"+seed_file+"';\n")
        codeF.write('close what=current;\n')
        
        codeF.write('open viewer=KEGG;\n')
        codeF.write('set context=KEGG;\n')
        codeF.write('update;\n')
        codeF.write('uncollapse nodes=all;\n')
        codeF.write('select nodes=all;\n')
        codeF.write('export what=CSV format='+self.kegg_format+" separator=comma counts=assigned file='"+kegg_file+"';\n")
        codeF.write('close what=current;\n')
        codeF.write('quit;\n')
        codeF.close()
        
        return code_file

    
    def run(self, input_files):
        print('***** Start Annotation *****')
        for name in input_files:
            file_name = name.split('/')[-1]
            subject_name = file_name.split('contigs')[0]
            print('******** Start Subject '+subject_name+' ********')

            '''prodigal'''
            print('***** Start Prodigal *****')
            prodigal_start = time.time()
            fna_file = os.path.join(self.prodigal_dir, subject_name+'prodigal.fna')
            faa_file = os.path.join(self.prodigal_dir, subject_name+'prodigal.faa')
            if not os.path.isfile(fna_file) and not os.path.isfile(faa_file):
                pro_mega_op = 'prodigal -p meta -i '+name+' -d '+fna_file+' -a '+faa_file
                os.system(pro_mega_op)
            prodigal_time = time.time()-prodigal_start
            print('Finish prodigal, time usage: '+str(prodigal_time))
            print('***** End Prodigal *****')
            
            '''diamond'''
            print('***** Start Diamond *****')
            diamond_start = time.time()
            diamond_file = os.path.join(self.diamond_dir, subject_name+'prodigal.m8')
            if not os.path.isfile(diamond_file):
                diamond_op = 'diamond blastp -d '+self.diamond_database+' -q '+faa_file+' -o '+diamond_file+' -e '+str(self.e_value)
                os.system(diamond_op)
            diamond_time = time.time()-diamond_start
            print('Finish diamond, time usage: '+str(diamond_time))
            print('***** End Diamond *****')
                     
            '''cut'''
            cut_file = os.path.join(self.diamond_dir, subject_name+'prodigal.contigs.m8')
            cut_op = 'cut -f1 '+diamond_file+'|rev|cut -f2- -d"_"|rev|paste - '+diamond_file+'|cut -f1,3- > '+cut_file
            os.system(cut_op)
            
            '''blast2rma'''
            print('***** Start Blast2RMA *****')
            blast2rma_file = os.path.join(self.blast2rma_dir, subject_name+'prodigal.contigs.rma')
            blast2rma_op = self.blast2rma_ref+' -i '+cut_file+' -f BlastTab -bm BlastP -o '+\
                    blast2rma_file+' -me '+str(self.e_value)+' -a2t '+self.acc2tax+' -a2eggnog '+self.acc2eggnog+\
                    ' -a2interpro2go '+self.acc2interpro+' -a2kegg '+self.acc2kegg+' -a2seed '+self.acc2seed
            os.system(blast2rma_op)
            print('***** End Blast2RMA *****')
            
            '''export'''
            blast2rma_name = subject_name+'prodigal.contigs.rma'
            code_file = self.write_text(blast2rma_name)
            export_op = 'xvfb-run --auto-servernum --server-num=1 ~/megan/MEGAN -g -c '+code_file
            os.system(export_op)

            print('******** Start Subject '+subject_name+' ********')
