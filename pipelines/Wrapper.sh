

# input_dir:        directory for the input folder
# mode_sel:         1 -- only adapter and quality control
#                   2 -- adapter and quality control & phix
#                   3 -- adapter and quality control & phix & primer
# adapter_ref:      reference file for adapter and quality control
# phix_ref:         reference file for phix
# primer_ref:       reference file for primer
# assembly_sel:     megahit -- use software MegaHit only
#                   metaspade -- use software MetaSpade only
#                   mega_meta -- use both MegaHit and MetaSpade
# assembly_mode:    indi -- individual assembly (work for megahit, metaspade, mega_meta)
#                   cross -- cross assembly (work for megahit)
#                   hybrid -- both individual and cross assembly (work for megahit)
# cong_length:      contig length
# diamond_database: the absolute path for NCBInr database
# e_value:          e value set-off for blast alignment
# blast2rma_ref:    the absolute path for blast2rma algorithm
# acc2tax2018:      the absolute path for Megan6 Protein accession to NCBI-taxonomy mapping file
# acc2eggnog2016:   the absolute path for Megan6 Protein accession to eggNOG mapping file
# acc2interpro2018: the absolute path for Megan6 Protein accession to InterPro mapping file
# acc2kegg2017:     the absolute path for Megan6 Protein accession to KEGG mapping file
# acc2seed2015:     the absolute path for Megan6 Protein accession to SEED mapping file
# collapse_rank:    Domain, Phylum, Class, Order, Family, Genus, Species
# select_rank:      Domain, Phylum, Class, Order, Family, Genus, Species [keep the same with collapse rank]
# tax_format:       export format, default='readName_to_taxonName'
# egg_format:       export format, default='readName_to_eggnogPath'
# seed_format:      export format, default='readName_to_seedPath'
# kegg_format:      export format, default='readName_to_keggPath'



############## Parameters ################
input_dir=/media/create/DATA3/SCELSE3_viral_HMwastewatersamples/Test
mode_sel=2
adapter_ref=truseq.fa.gz
phix_ref=phix174_ill.ref.fa.gz
primer_ref=A
assembly_sel=megahit
assembly_mode=indi
cong_length=1000
diamond_database=~/tools/DATABASE/Diamondnr/nr
e_value=1e-5
blast2rma_ref=~/megan/tools/blast2rma
acc2tax2018=~/tools/megan/database/prot_acc2tax-Jun2018X1.abin
acc2eggnog2016=~/tools/megan/database/acc2eggnog-Oct2016X.abin
acc2interpro2018=~/tools/megan/database/acc2interpro-June2018X.abin
acc2kegg2017=~/tools/megan/database/acc2kegg-Dec2017X1-ue.abin
acc2seed2015=~/tools/megan/database/acc2seed-May2015XX.abin
collapse_rank=species
select_rank=species
tax_format=readName_to_taxonName
egg_format=readName_to_eggnogPath
seed_format=readName_to_seedPath
kegg_format=readName_to_keggPath
###########################################

python pipeline.py $input_dir $mode_sel $adapter_ref $phix_ref $primer_ref $assembly_sel $assembly_mode $cong_length $diamond_database $e_value $blast2rma $acc2tax2018 $acc2eggnog2016 $acc2interpro2018 $acc2kegg2017 $acc2seed2015 $collapse_rank $select_rank $tax_format $egg_format $seed_format $kegg_format > $input_dir/Test.txt
